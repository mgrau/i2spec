"""Local Dunham models of the I₂ B-X near-infrared bands (v' = 0, v'' = 12-17).

These are *local* models: closed-form polynomials in (v''+1/2) and J(J+1), fitted to one corner of
the quantum-number space, where they beat the global potential model. ``RovibronicModel`` stays the
general tool; this module covers 755-815 nm, where the Hannover group published dedicated fits.

Two parameter sets ship with the package:

* ``knockel2004`` -- Knoeckel, Bodermann & Tiemann, EPJD 28, 199 (2004), Table 1. Fit sigma 40 kHz;
  stated 1 sigma prediction uncertainty below 200 kHz over "roughly 775 to 815 nm".
* ``liao2010`` -- Liao et al., JOSA B 27, 1208 (2010), Table 3: the same functional form and
  parameter count refitted with 27 new comb-referenced (0-12) and (0-13) frequencies, claimed to
  predict 755-815 nm to better than 0.2 MHz.

**The 2010 set as printed does not reproduce its own fit.** Table 3 carries only 3-9 significant
digits for strongly correlated parameters, and the predictions it yields miss the paper's own
measurements by up to 51 MHz -- inside the very window it claims. ``load_local_nir`` therefore
refuses ``liao2010`` unless the caller passes ``trust_printed_digits=True``; the numbers are shipped
for the record, with the paper's domain and uncertainty metadata. The 2004 set has no such problem:
inside its own 775-815 nm window it reproduces the Bodermann et al. 2000 anchors to 93 kHz rms. It
does fail *outside* that window -- by up to 45 MHz against Liao's (0-12) and (0-13) measurements at
755-770 nm -- which is exactly what ``Domain.covers_wavelength`` is for, and what the 2010 refit was
meant to fix. See ``docs/research/nir-model-2010.md`` for the analysis.

The model is
    nu(v'=0, v'', J'', branch) = sum_k S_0k [J'(J'+1)]^k - sum_lk S_lk (v''+1/2)^l [J''(J''+1)]^k
with J' = J'' +/- 1 (Knoeckel 2004 eqs. 1-2). The S are not Dunham Y coefficients: they are valid
only inside the fitted window, which is why both papers rename them.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from importlib import resources

from .constants import MHZ_PER_CM

#: Parameter sets shipped in ``i2spec/data/<name>_nir.json``.
SETS = ("knockel2004", "liao2010")


@dataclass(frozen=True)
class Domain:
    """Where a local model may be used: quantum numbers, wavelengths, and the published uncertainty."""

    v_upper: int
    v_lower: tuple[int, int]
    J_lower_max: int
    wavelength_nm: tuple[float, float]
    uncertainty_MHz: float
    uncertainty_basis: str = ""
    J_lower_recommended: int | None = None

    def covers(self, v_upper: int, v_lower: int, J_lower: int) -> bool:
        """Is this line inside the fitted window? Wavelength is implied by the quantum numbers."""
        return (v_upper == self.v_upper
                and self.v_lower[0] <= v_lower <= self.v_lower[1]
                and 0 <= J_lower <= self.J_lower_max)

    def covers_wavelength(self, nm: float) -> bool:
        """Is this vacuum wavelength (nm) inside the model's published range?"""
        lo, hi = self.wavelength_nm
        return lo <= nm <= hi

    def covers_wavenumber(self, nu: float) -> bool:
        """Is this vacuum wavenumber (cm⁻¹) inside the model's published wavelength range?"""
        return self.covers_wavelength(1e7 / nu)


@dataclass(frozen=True)
class LocalNIRModel:
    """One published local Dunham model, with its domain and how well its printed digits work."""

    name: str
    source: str
    doi: str
    upper: tuple[float, ...]  # S_0k, k = 0, 1, ... (v' = 0 folded into S_00)
    lower: dict[tuple[int, int], float]  # S_lk of the X state
    domain: Domain
    notes: str = ""
    validation: dict | None = None

    def transition(self, v_upper: int, v_lower: int, J_lower: int, branch: str) -> float:
        """Wavenumber (cm⁻¹) of the P or R line (v_upper-v_lower) with lower-state J.

        The polynomial is evaluated wherever it is asked for; use ``domain.covers`` to find out
        whether the answer means anything. Only v' = 0 exists in these models.
        """
        if v_upper != self.domain.v_upper:
            raise ValueError(f"{self.name} covers only v' = {self.domain.v_upper}, not {v_upper}")
        if branch == "R":
            J_upper = J_lower + 1
        elif branch == "P":
            J_upper = J_lower - 1
        else:
            raise ValueError("branch must be 'P' or 'R' (Q lines are forbidden for 0u+ - 0g+)")
        if J_lower < 0 or J_upper < 0:
            raise ValueError(f"{branch}({J_lower}) has no upper level")
        xu, xl, u = J_upper * (J_upper + 1), J_lower * (J_lower + 1), v_lower + 0.5
        upper = sum(c * xu**k for k, c in enumerate(self.upper))
        lower = sum(c * u**l * xl**k for (l, k), c in self.lower.items())
        return upper - lower

    def transition_MHz(self, v_upper: int, v_lower: int, J_lower: int, branch: str) -> float:
        return self.transition(v_upper, v_lower, J_lower, branch) * MHZ_PER_CM

    def covers(self, v_upper: int, v_lower: int, J_lower: int) -> bool:
        return self.domain.covers(v_upper, v_lower, J_lower)

    @property
    def uncertainty_MHz(self) -> float:
        """The 1σ prediction uncertainty the source claims inside the domain."""
        return self.domain.uncertainty_MHz

    @property
    def reproduction_MHz(self) -> float | None:
        """RMS (MHz) by which the printed parameters miss measured frequencies inside the model's window."""
        return None if self.validation is None else self.validation.get("in_window_rms_MHz")


def load_local_nir(name: str = "knockel2004", trust_printed_digits: bool = False) -> LocalNIRModel:
    """Load a published local NIR model shipped in ``i2spec/data/<name>_nir.json``.

    ``liao2010`` needs ``trust_printed_digits=True``: its printed parameters do not reproduce its
    own fit (module docstring). ``knockel2004`` is the usable set.
    """
    if name not in SETS:
        raise ValueError(f"unknown local NIR model {name!r}; choose from {SETS}")
    text = resources.files("i2spec").joinpath("data").joinpath(f"{name}_nir.json").read_text()
    p = json.loads(text)
    d = p["domain"]
    model = LocalNIRModel(
        name=p["name"], source=p["source"], doi=p["doi"],
        upper=tuple(p["upper"][f"0{k}"] for k in range(len(p["upper"]))),
        lower={(int(k[0]), int(k[1])): c for k, c in p["lower"].items()},
        domain=Domain(v_upper=d["v_upper"], v_lower=tuple(d["v_lower"]), J_lower_max=d["J_lower_max"],
                      wavelength_nm=tuple(d["wavelength_nm"]), uncertainty_MHz=d["uncertainty_MHz"],
                      uncertainty_basis=d.get("uncertainty_basis", ""),
                      J_lower_recommended=d.get("J_lower_recommended")),
        notes=p.get("notes", ""), validation=p.get("validation"))
    usable = (model.validation or {}).get("printed_digits_reproduce_the_fit", True)
    if not usable and not trust_printed_digits:
        raise ValueError(
            f"{name}: the parameters as printed miss measured frequencies inside their own window by "
            f"{model.reproduction_MHz} MHz rms, against the {model.uncertainty_MHz} MHz the source claims, "
            f"because the table is rounded too hard to reconstruct a correlated fit. Pass "
            f"trust_printed_digits=True to load them anyway; see docs/research/nir-model-2010.md")
    return model
