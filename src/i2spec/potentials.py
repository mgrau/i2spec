"""Published potential-energy curves for the I₂ X¹Σg⁺ and B³Π(0u⁺) states (the released model's own curves are in
i2spec.research and i2spec.research_expanded)."""

from __future__ import annotations

import json
from dataclasses import dataclass, replace
from importlib import resources

import numpy as np

from numpy.polynomial import polynomial as P

from .constants import PUBLISHED_PARAMETERS


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


def parameter_set(name: str = PUBLISHED_PARAMETERS) -> dict:
    """The parameter set ``i2spec/data/<name>.json`` as a dict (potentials and provenance)."""
    return json.loads(resources.files("i2spec").joinpath("data").joinpath(f"{name}.json").read_text())


def load_potentials(name: str = PUBLISHED_PARAMETERS) -> dict[str, XRepPotential]:
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
