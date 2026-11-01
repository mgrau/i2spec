"""Line lookup: the query layer behind the command line, the TUI and the explorer's export.

Line lists (strengths, lower-level energies) are built on the model's potentials and cached (``i2spec.intensity``).
Positions, hyperfine components and uncertainties come from the model's level tables (``i2spec.preferred``), with
the uncertainty of ``ResearchModel.uncertainty``. A line whose levels the tables do not hold is flagged outside the
fitted domain: its position comes from the same potentials, its components from the published hyperfine formulae.
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass
from functools import lru_cache

import numpy as np

from .constants import ISOTOPOLOGUES, MHZ_PER_CM
from .intensity import C2, NEAR_FROM, intensity_model, master_line_list, with_dissociation_lines
from .model import RovibronicModel
from .spectrum import air_to_vacuum, doppler_fwhm, molecular_mass, vacuum_to_air

#: Units accepted for a search range: vacuum nm, standard-air nm, wavenumber, frequency.
UNITS = ("nm", "nm-air", "cm-1", "MHz", "GHz", "THz")
_LABEL = re.compile(r"^([PR])\(?(\d+)\)?\s+(\d+)\s*-\s*(\d+)$", re.IGNORECASE)
_QUANTITY = re.compile(r"^\s*([-+0-9.eE]+)\s*([A-Za-z0-9/^-]*)\s*$")
#: One end of a range: a number and, optionally, its unit. A unit starts with a letter, so the "-" between
#: two numbers is not read as part of one; units such as "cm-1" and "nm-air" keep theirs.
_RANGE_END = re.compile(r"(\d+\.?\d*(?:[eE][-+]?\d+)?|\.\d+(?:[eE][-+]?\d+)?)\s*([A-Za-z]+(?:\^?-1|-air)?)?")


def format_mhz(value):
    """An uncertainty in MHz to two figures or so: 0.42, 4.2, 42, 1,203."""
    return f"{value:.2f}" if value < 1 else f"{value:.1f}" if value < 10 else f"{value:,.0f}"


def to_wavenumber(value, unit="nm"):
    """Vacuum wavenumber (cm⁻¹) from a value in any of UNITS, which must be positive."""
    if unit not in UNITS:
        raise ValueError(f"unknown unit {unit!r}; choose from {UNITS}")
    if not float(value) > 0:  # also catches NaN
        raise ValueError(f"a wavelength or frequency must be positive, not {value} {unit}")
    if unit == "nm":
        return 1e7 / float(value)
    if unit == "nm-air":
        return 1e7 / float(air_to_vacuum(value))
    if unit == "cm-1":
        return float(value)
    if unit == "MHz":
        return float(value) / MHZ_PER_CM
    if unit == "GHz":
        return float(value) * 1e3 / MHZ_PER_CM
    return float(value) * 1e6 / MHZ_PER_CM   # THz


def from_wavenumber(nu, unit="nm"):
    """A vacuum wavenumber (cm⁻¹) expressed in any of UNITS."""
    if unit == "nm":
        return 1e7 / nu
    if unit == "nm-air":
        return float(vacuum_to_air(1e7 / nu))
    if unit == "cm-1":
        return nu
    if unit == "MHz":
        return nu * MHZ_PER_CM
    if unit == "GHz":
        return nu * MHZ_PER_CM / 1e3
    if unit == "THz":
        return nu * MHZ_PER_CM / 1e6
    raise ValueError(f"unknown unit {unit!r}; choose from {UNITS}")


def parse_quantity(text, default_unit="nm"):
    """("532.245 nm") -> (532.245, "nm"). The unit may be attached or omitted."""
    m = _QUANTITY.match(str(text))
    if not m:
        raise ValueError(f"cannot read {text!r} as a number with a unit")
    unit = m[2] or default_unit
    aliases = {"nm": "nm", "nmair": "nm-air", "nm-air": "nm-air", "air": "nm-air", "cm-1": "cm-1", "cm": "cm-1",
               "wn": "cm-1", "mhz": "MHz", "thz": "THz", "ghz": "GHz"}
    unit = aliases.get(unit.lower().replace("^", "").replace("/", ""), unit)
    value = float(m[1])
    to_wavenumber(value, unit)  # rejects an unknown unit and a value that is not positive
    return value, unit


def parse_range(text, default_unit="nm"):
    """("532.2-532.3 nm") -> (532.2, 532.3, "nm"). A unit written once applies to both ends; the ends
    may be separated by "-", "to", ".." or a space."""
    text = re.sub(r"\bto\b|\.\.", " ", str(text))
    ends = list(_RANGE_END.finditer(text))
    if len(ends) != 2 or re.sub(_RANGE_END, "", text).strip(" -"):
        raise ValueError(f"cannot read {text.strip()!r} as a range; write it like '532.2-532.3 nm'")
    units = {m[2] for m in ends if m[2]}
    if len(units) > 1:
        raise ValueError(f"give both ends of the range in the same unit, not {' and '.join(sorted(units))}")
    unit = units.pop() if units else default_unit
    (low, unit_out), (high, _) = (parse_quantity(m[1] + unit) for m in ends)
    return low, high, unit_out


def parse_label(text):
    """("R(56) 32-0") -> ("R", 56, 32, 0)."""
    m = _LABEL.match(str(text).strip())
    if not m:
        raise ValueError(f"cannot read {text!r} as a line; expected e.g. 'R(56) 32-0'")
    return m[1].upper(), int(m[2]), int(m[3]), int(m[4])


@dataclass(frozen=True)
class Line:
    """One B-X line: its assignment, position, strength at a temperature, and how well we know it."""

    isotopologue: str
    branch: str
    J_lower: int
    v_upper: int
    v_lower: int
    nu: float  # vacuum wavenumber, cm⁻¹
    strength: float  # integrated cross section at ``temperature``, cm
    E_lower: float  # cm⁻¹ above X(0, 0)
    temperature: float

    @property
    def label(self):
        return f"{self.branch}({self.J_lower}) {self.v_upper}-{self.v_lower}"

    @property
    def J_upper(self):
        return self.J_lower + (1 if self.branch == "R" else -1)

    def position(self, unit="nm"):
        return from_wavenumber(self.nu, unit)

    @property
    def doppler_fwhm_MHz(self):
        return float(doppler_fwhm(self.nu, self.temperature, molecular_mass(self.isotopologue))) * MHZ_PER_CM

    @property
    def uncertainty_MHz(self):
        return uncertainty(self)[0]

    @property
    def source(self):
        """How the apps name the model behind this line's position."""
        from .preferred import BEYOND, OUTSIDE, PUBLIC
        p = _held(self)
        return PUBLIC if p is not None and p[1] < BEYOND else f"{PUBLIC}, {OUTSIDE}"

    @property
    def flags(self):
        p = _held(self)
        if p is None:
            from .preferred import OUTSIDE
            return (OUTSIDE, QUASI_BOUND) if quasi_bound(self) else (OUTSIDE,)
        from .preferred import flags
        return flags(p[1])+((QUASI_BOUND,) if quasi_bound(self) else ())


def _held(line):
    """(sigma_center MHz, status code, class code) where the model's level tables hold this line, else None. A status
    code of preferred.BEYOND or more: a level lies beyond the fitted domain, its sigma from the release's edge rules."""
    from .preferred import preferred
    sigma, status, cls = preferred(line.isotopologue).uncertainties([line.branch], [line.J_lower], [line.v_upper],
                                                                    [line.v_lower])
    return None if status[0] < 0 else (float(sigma[0]), int(status[0]), int(cls[0]))


#: The 1-sigma quoted for a line the model's tables do not hold (MHz): a rough bound, not an estimate
OUTSIDE_SIGMA_MHZ = 1000.0
#: The flag of a line whose B level lies above the dissociation limit, held only by the centrifugal barrier
QUASI_BOUND = "quasi-bound"
QUASI_NOTE = ("; the upper level lies above the rotationless B dissociation limit, held by the centrifugal barrier "
              "(quasi-bound), and the model treats it as bound, with no tunneling lifetime or broadening")


@lru_cache(maxsize=None)
def b_limit(isotopologue):
    """The B dissociation limit above X(0, 0) (cm-1) of the released model."""
    from .preferred import model
    levels = model(isotopologue)._levels_x
    return float(levels.asymptote["B"]-levels.energy("X", 0, 0))


def quasi_bound(line):
    """Whether the line's upper level lies above the B dissociation limit (it is then behind the barrier)."""
    return line.E_lower+line.nu > b_limit(line.isotopologue)


def uncertainty(line: Line):
    """1σ uncertainty (MHz) of the hyperfine-free position, and the evidence for it (``ResearchModel.uncertainty``)."""
    p = _held(line)
    if p is not None:
        from .preferred import basis
        return round(p[0], 3), basis(p[1], p[2], line.isotopologue)+(QUASI_NOTE if quasi_bound(line) else "")
    if quasi_bound(line):
        return OUTSIDE_SIGMA_MHZ, ("Outside the fitted domain: the upper level lies above the rotationless B dissociation "
                                   "limit, held by the centrifugal barrier (quasi-bound); the model treats it as bound, "
                                   "with no tunneling lifetime or broadening, and the uncertainty is a rough bound")
    return OUTSIDE_SIGMA_MHZ, ("Outside the fitted domain: the model's level tables do not hold a level of this line "
                               "(beyond its centrifugal-barrier rule, B v′ > 86, X v″ > 113, or J″ > 445); the position comes from "
                               "the same potentials, and the uncertainty is a rough bound")


@dataclass(frozen=True)
class Component:
    """One hyperfine component of a line."""

    label: str | None  # a1, a2, ... for the main (ΔF = ΔJ) components
    offset_MHz: float  # from the hyperfine-free line position
    strength: float  # fraction of the line's strength
    F_upper: int
    F_lower: int
    frequency_MHz: float
    I_upper: int | None = None  # total nuclear spin of each level (approximately good at high J)
    I_lower: int | None = None
    upper_level: int | None = None  # hyperfine eigenstate indices (hyperfine.Component)
    lower_level: int | None = None
    uncertainty_MHz: float | None = None  # 1σ of frequency_MHz, where the model's tables hold the line

    @property
    def name(self):
        return self.label or f"F {self.F_lower}→{self.F_upper}"


@dataclass
class Result:
    """What a search found: the lines, how many matched, and whether the list was cut short."""

    lines: list[Line]
    total: int
    temperature: float

    @property
    def truncated(self):
        return len(self.lines) < self.total

    def __len__(self):
        return len(self.lines)

    def __iter__(self):
        return iter(self.lines)


class Catalog:
    """Line lists and hyperfine patterns for the isotopologues, built once and cached.

    The first search for an isotopologue builds its master line list, which takes a couple of
    minutes and is then cached on disk (``~/.cache/i2spec``). ``masters`` and ``models`` let a
    caller inject prepared objects, which the tests use.
    """

    def __init__(self, nu_range=(11000.0, 20100.0), s_min=1e-27, temperature=293.15, masters=None, models=None):
        self.nu_range = (float(nu_range[0]), float(nu_range[1]))
        self.s_min = float(s_min)
        self.temperature = float(temperature)
        self._masters = dict(masters or {})
        self._models = dict(models or {})
        #: isotopologues whose list reaches the B limit; an injected list counts as complete
        self._extended = set(self._masters)
        #: isotopologue -> (master list, the same with the tables' positions), so a search does not redo the lookup
        self._positioned = {}

    def master(self, isotopologue="127I2", near_dissociation=False):
        """The master line list; with ``near_dissociation``, extended to the B limit (a further build, once)."""
        if isotopologue not in self._masters:
            self._masters[isotopologue] = master_line_list(intensity_model(isotopologue), *self.nu_range,
                                                           T_range=(200.0, 600.0), S_min=self.s_min)
        if near_dissociation and self.nu_range[1] > NEAR_FROM and isotopologue not in self._extended:
            self._masters[isotopologue] = with_dissociation_lines(self._masters[isotopologue], S_min=self.s_min,
                                                                  nu_max=self.nu_range[1])
            self._extended.add(isotopologue)
        master = self._masters[isotopologue]
        if self._positioned.get(isotopologue, (None,))[0] is not master:
            self._positioned[isotopologue] = master, preferred_positions(master)
        return self._positioned[isotopologue][1]

    def model(self, isotopologue="127I2"):
        """The B-spline model of the lines the level tables do not hold: their positions, lower-level energies and
        hyperfine components (published formulae)."""
        if isotopologue not in self._models:
            self._models[isotopologue] = RovibronicModel(isotopologue)
        return self._models[isotopologue]

    def search(self, low, high, unit="nm", isotopologue="127I2", temperature=None, min_strength=0.0, limit=200,
               sort="nu") -> Result:
        """Lines between ``low`` and ``high`` (in ``unit``), strongest first or by position."""
        if isotopologue not in ISOTOPOLOGUES:
            raise ValueError(f"unknown isotopologue {isotopologue!r}; choose from {sorted(ISOTOPOLOGUES)}")
        T = float(temperature if temperature is not None else self.temperature)
        nu_lo, nu_hi = sorted(to_wavenumber(v, unit) for v in (low, high))
        lines = self.master(isotopologue, near_dissociation=nu_hi > NEAR_FROM).at(T, nu_lo, nu_hi, S_min=min_strength)
        order = np.argsort(lines.nu if sort == "nu" else -lines.S)
        keep = order[:limit] if limit else order
        found = [Line(isotopologue, "R" if lines.branch[k] > 0 else "P", int(lines.J_lower[k]), int(lines.v_upper[k]),
                      int(lines.v_lower[k]), float(lines.nu[k]), float(lines.S[k]), float(lines.E_lower[k]), T)
                 for k in keep]
        return Result(found, len(lines), T)

    def line(self, label, isotopologue="127I2", temperature=None) -> Line:
        """The line named "R(56) 32-0", found in the list by its assignment."""
        branch, J, v_upper, v_lower = parse_label(label)
        T = float(temperature if temperature is not None else self.temperature)
        master = self.master(isotopologue, near_dissociation=v_upper > 60)
        sel = ((master.v_upper == v_upper) & (master.v_lower == v_lower) & (master.J_lower == J)
               & (master.branch == (1 if branch == "R" else -1)))
        if not sel.any():
            nu = self.transition(isotopologue, branch, J, v_upper, v_lower)
            raise LookupError(f"{label} of {isotopologue} is not in the line list "
                              f"(it would lie at {1e7 / nu:.4f} nm)")
        k = int(np.flatnonzero(sel)[0])
        strength = master.strength0[k] * math.exp(-C2 * master.E_lower[k] / T) / master.partition_function(T)
        return Line(isotopologue, branch, J, v_upper, v_lower, float(master.nu[k]), float(strength),
                    float(master.E_lower[k]), T)

    def transition(self, isotopologue, branch, J, v_upper, v_lower):
        """Hyperfine-free wavenumber (cm⁻¹): from the level tables, or the potentials where the tables stop."""
        line = Line(isotopologue, branch, J, v_upper, v_lower, 0.0, 0.0, 0.0, self.temperature)
        if _held(line) is not None:
            from .preferred import preferred
            return float(preferred(isotopologue).positions([branch], [J], [v_upper], [v_lower])[0])
        return self.model(isotopologue).transition(v_upper, v_lower, J, branch)

    def components(self, line: Line, dJ=2, main_only=False) -> list[Component]:
        """Hyperfine components of a line, with absolute frequencies. Where the level tables hold the line, the
        model's, always with the J ± 2 couplings and each component's 1σ; elsewhere the published formulae's, with
        ``dJ``."""
        p = _held(line)
        from .preferred import BEYOND
        if p is not None and p[1] < BEYOND:
            from .preferred import components, sigma_component
            nu0, comps = components(line.branch, line.J_lower, line.v_upper, line.v_lower, line.isotopologue)
            sigma = lambda c: sigma_component(p[0], line.branch, line.J_lower, line.v_upper, line.v_lower,  # noqa: E731
                                              c.offset, line.isotopologue)
        else:
            nu0, comps = self.model(line.isotopologue).hyperfine_components(line.v_upper, line.v_lower, line.J_lower,
                                                                            line.branch, dJ=dJ)
            sigma = lambda c: None  # noqa: E731
        out = [Component(c.label, c.offset, c.strength, c.F_upper, c.F_lower, nu0 + c.offset, c.I_upper, c.I_lower,
                         c.upper_level, c.lower_level, sigma(c))
               for c in comps if c.label or not main_only]
        return sorted(out, key=lambda c: c.offset_MHz)


    def center_only(self, line: Line, comps=None):
        """Why the model gives this line's center only, or None. The model does not predict the hyperfine of a few
        B levels (where a second-order energy denominator nearly vanishes), nor of a line whose J ± 2 coupling
        partner is dissociated; their lines have a center and no components, and the model says why. ``comps``:
        the line's components, when already computed."""
        if comps is None:
            try:
                comps = self.components(line, main_only=True)
            except ValueError as error:               # e.g. a coupling partner beyond the dissociation limit
                return str(error)
        if any(c.label for c in comps):
            return None
        from .preferred import model
        try:
            model(line.isotopologue).predict(line.label, "a1")
        except ValueError as error:                   # the model's own reason
            return str(error)
        except Exception:                             # no reason to give: the plain statement below
            pass
        return "the model gives no hyperfine components for this line"


def preferred_positions(master):
    """A master line list with the level tables' positions where they hold the line."""
    from dataclasses import replace
    from .preferred import preferred
    if getattr(master, "_preferred", False):
        return master
    nu = preferred(master.isotopologue).positions(np.where(master.branch > 0, "R", "P"), master.J_lower,
                                                  master.v_upper, master.v_lower)
    out = replace(master, nu=np.where(np.isfinite(nu), nu, master.nu))
    out._preferred = True
    return out
