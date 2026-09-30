"""Potential-energy curves for the I₂ X¹Σg⁺ and B³Π(0u⁺) states."""

from __future__ import annotations

import json
from dataclasses import dataclass, replace
from importlib import resources

import numpy as np

from .constants import DEFAULT_PARAMETERS
from numpy.polynomial import chebyshev as C
from numpy.polynomial import polynomial as P


@dataclass(frozen=True)
class XRepPotential:
    """Hannover "X-representation" potential with its extensions and BO corrections.

    For R_I <= R <= R_O the potential is a power series in X = (R - Rm)/(R + b Rm).
    Below R_I it continues as A_I exp(-B_I (R - R_I)), above R_O as
    De - sum_n C_n/R^n - A_O exp(-B_O (R - R_O)). The Born-Oppenheimer corrections
    follow Salumbides et al., EPJD 47, 171 (2008), eqs. (5)-(9); like there, they are
    evaluated without extensions over the whole grid.
    """

    Rm: float
    b: float
    a: tuple[float, ...]
    RI: float
    AI: float
    BI: float
    RO: float
    AO: float
    BO: float
    De: float
    C: dict[int, float]
    alpha: tuple[float, ...] = ()
    vad: tuple[float, ...] = ()
    vad_power: int = 1

    def x(self, R):
        return (R - self.Rm) / (R + self.b * self.Rm)

    def series(self, R):
        return P.polyval(self.x(R), self.a)

    def series_derivative(self, R):
        return P.polyval(self.x(R), P.polyder(self.a)) * self.Rm * (1 + self.b) / (R + self.b * self.Rm) ** 2

    def dispersion(self, R):
        return self.De - sum(c / R**n for n, c in self.C.items())

    def __call__(self, R):
        R = np.asarray(R, dtype=float)
        V = np.empty_like(R)
        inner, outer = R < self.RI, R > self.RO
        mid = ~(inner | outer)
        V[inner] = self.AI * np.exp(-self.BI * (R[inner] - self.RI))
        V[mid] = self.series(R[mid])
        V[outer] = self.dispersion(R[outer]) - self.AO * np.exp(-self.BO * (R[outer] - self.RO))
        return V

    def nonadiabatic(self, R, mass_ratio=1.0):
        """alpha(R) in the centrifugal term, for an isotopologue with mu_ref/mu = mass_ratio."""
        if not self.alpha:
            return np.zeros_like(np.asarray(R, dtype=float))
        return mass_ratio * 2 * self.Rm / (R + self.Rm) * P.polyval(self.x(R), self.alpha)

    def adiabatic(self, R, mass_ratio=1.0):
        """V_corr(R) = (1 - mu_ref/mu) V_ad(R); zero for the reference isotopologue."""
        if not self.vad:
            return np.zeros_like(np.asarray(R, dtype=float))
        return (1 - mass_ratio) * (2 * self.Rm / (R + self.Rm)) ** self.vad_power * P.polyval(self.x(R), self.vad)

    def with_continuous_extensions(self) -> XRepPotential:
        """Copy whose A_I, B_I, A_O, B_O make V and dV/dR continuous at R_I and R_O."""
        AI = self.series(self.RI)
        BI = -self.series_derivative(self.RI) / AI
        ddisp = sum(n * c / self.RO ** (n + 1) for n, c in self.C.items())
        AO = self.dispersion(self.RO) - self.series(self.RO)
        BO = (self.series_derivative(self.RO) - ddisp) / AO
        return replace(self, AI=float(AI), BI=float(BI), AO=float(AO), BO=float(BO))


@dataclass(frozen=True)
class MLRPotential:
    """Morse/Long-Range potential (Le Roy & Henderson 2007; Le Roy, Dattani et al. 2011).

        V(R) = De [1 - u(R)/u(Re) exp(-beta(R) y_p^eq(R))]^2,   y_n^r(R) = (R^n - r^n)/(R^n + r^n)

    with u(R) = sum_n C_n/R^n the long-range function, beta(R) = y_p^ref beta_inf +
    (1 - y_p^ref) sum_i beta_i (y_q^ref)^i, and beta_inf = ln(2 De / u(Re)).

    Why this form, in one line: V(Re) = 0 and V(R -> inf) = De - u(R) + u(R)^2/4De hold *by
    construction*, whatever the beta_i do, so the long-range tail is theory, not extrapolated fit.
    The Hannover X-representation instead splices a power series onto a tail at R_O, and outside the
    fitted data that seam is free to drift: it is +7.76 cm-1 at v'' = 48 and -21.1 cm-1 at v'' = 54
    (docs/design/fitting.md). The exponent parameters p and q set how fast the switch happens; p must
    exceed the spread of the u(R) powers for the tail to be reached correctly.

    The Born-Oppenheimer corrections are not part of the MLR: ``bo`` supplies them, so they stay the
    published functions of whichever potential it is taken from. No damping functions yet; they matter
    at small R, well inside the region any of our data reaches.
    """

    De: float
    Re: float
    C: dict[int, float]
    beta: tuple[float, ...]
    p: int = 5
    q: int = 3
    Rref: float | None = None
    bo: XRepPotential | None = None
    #: "power" for the published convention, sum beta_i y^i; "chebyshev" for the same functions in a
    #: basis that a least-squares fit can actually resolve.
    basis: str = "power"
    #: Energy of the potential minimum. 0 for the X state, whose minimum defines the energy origin;
    #: for B it is the electronic origin, so the asymptote sits at Te + De.
    Te: float = 0.0
    #: An added constant q_far in the centrifugal term beyond q_far_R (a tanh switch of width q_far_w).
    #: The published alpha(R) is an effective function fitted inside the well; near dissociation the
    #: rotational energy of the last B levels needs its own value (docs/research/orsay-atlas-19700-20035.md).
    q_far: float = 0.0
    q_far_R: float = 5.0
    q_far_w: float = 0.5

    def far_switch(self, R):
        return 0.5 * (1.0 + np.tanh((np.asarray(R, dtype=float) - self.q_far_R) / self.q_far_w))

    @property
    def reference(self) -> float:
        return self.Re if self.Rref is None else self.Rref

    def u(self, R):
        """Long-range function sum_n C_n / R^n, cm^-1."""
        R = np.asarray(R, dtype=float)
        return sum(c / R ** n for n, c in self.C.items())

    def beta_inf(self) -> float:
        return float(np.log(2 * self.De / self.u(self.Re)))

    def beta_of(self, R):
        yp = _radial_variable(R, self.p, self.reference)
        yq = _radial_variable(R, self.q, self.reference)
        return yp * self.beta_inf() + (1 - yp) * self._exponent(yq)

    def _exponent(self, y):
        """Sum beta_i b_i(y), in whichever polynomial basis this potential carries.

        y runs over [-1, 1] by construction, where powers of y are famously close to collinear: fitting
        thirty of them leaves a Jacobian whose condition number reaches 1e17, and an optimiser crawling
        in a valley it cannot see the floor of. Chebyshev polynomials span exactly the same functions
        and are nearly orthogonal there, which is the whole difference (docs/design/global-fit.md).
        """
        if self.basis == "chebyshev":
            return C.chebval(y, self.beta)
        return P.polyval(y, self.beta)

    def __call__(self, R):
        R = np.asarray(R, dtype=float)
        y = _radial_variable(R, self.p, self.Re)
        shape = self.u(R) / self.u(self.Re) * np.exp(-self.beta_of(R) * y)
        return self.Te + self.De * (1.0 - shape) ** 2

    def nonadiabatic(self, R, mass_ratio=1.0):
        far = mass_ratio * self.q_far * self.far_switch(R) if self.q_far else 0.0
        if self.bo is None:
            return np.zeros_like(np.asarray(R, dtype=float)) + far
        return self.bo.nonadiabatic(R, mass_ratio) + far

    def adiabatic(self, R, mass_ratio=1.0):
        if self.bo is None:
            return np.zeros_like(np.asarray(R, dtype=float))
        return self.bo.adiabatic(R, mass_ratio)


def _radial_variable(R, n, ref):
    """y_n^ref(R) = (R^n - ref^n)/(R^n + ref^n): -1 at R = 0, 0 at ref, +1 as R -> infinity."""
    ref = float(ref)
    if not np.isfinite(ref) or abs(ref) > 1e3:      # an optimiser exploring nonsense, not a molecule
        raise OverflowError(f"reference radius {ref} is not a bond length")
    Rn, rn = np.asarray(R, dtype=float) ** n, ref ** n
    return (Rn - rn) / (Rn + rn)


def parameter_set(name: str = DEFAULT_PARAMETERS) -> dict:
    """The parameter set ``i2spec/data/<name>.json`` as a dict (potentials, provenance, corrections)."""
    return json.loads(resources.files("i2spec").joinpath("data").joinpath(f"{name}.json").read_text())


def load_mlr(name: str, bo: XRepPotential) -> MLRPotential:
    """An MLR potential from ``i2spec/data/<name>.json`` (the data/potentials schema), with the
    Born-Oppenheimer functions of ``bo``."""
    d = json.loads(resources.files("i2spec").joinpath("data").joinpath(f"{name}.json").read_text())
    rref = d.get("Rref") if "Rref" in d else d["Rref_factor"] * d["Re"]
    return MLRPotential(De=d["De"], Re=d["Re"], C={int(k): v for k, v in d["C"].items()}, beta=tuple(d["beta"]),
                        p=d["p"], q=d["q"], Rref=rref, Te=d.get("Te", 0.0), bo=bo,
                        **{k: d[k] for k in ("q_far", "q_far_R", "q_far_w") if k in d})


def load_extended(name: str = DEFAULT_PARAMETERS):
    """The parameter set's extended-range potentials, if it names any: ({state: MLRPotential}, {state: v_from}).

    A set may carry a second pair of potentials for the levels its published curves cannot reach
    (``"extended": {"X": file, "B": file, "from_v": {"X": 18, "B": 44}}``); RovibronicModel takes the
    levels at and above ``from_v`` from them. Returns (None, None) when the set has none.
    """
    d = parameter_set(name)
    if "extended" not in d:
        return None, None
    base = load_potentials(name)
    ext = {st: load_mlr(d["extended"][st], base[st]) for st in ("X", "B")}
    return ext, {st: int(v) for st, v in d["extended"]["from_v"].items()}


def load_potentials(name: str = DEFAULT_PARAMETERS) -> dict[str, XRepPotential]:
    """Load a published parameter set shipped in ``i2spec/data/<name>.json``."""
    states = parameter_set(name)["states"]
    potentials = {}
    for state, p in states.items():
        boc = p.get("boc", {})
        potentials[state] = XRepPotential(
            Rm=p["Rm"], b=p["b"], a=tuple(p["a"]),
            RI=p["inner"]["R"], AI=p["inner"]["A"], BI=p["inner"]["B"],
            RO=p["outer"]["R"], AO=p["outer"]["A"], BO=p["outer"]["B"],
            De=p["outer"]["De"], C={int(n): c for n, c in p["outer"]["C"].items()},
            alpha=tuple(boc.get("alpha", ())), vad=tuple(boc.get("vad", ())),
            vad_power=boc.get("vad_power", 1),
        )
    return potentials
