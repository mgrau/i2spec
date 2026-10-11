"""Line lookup: the query layer behind the command line and the TUI.

Positions and strengths come from the cached master line list (``i2spec.intensity``) and the
hyperfine patterns from ``RovibronicModel``. This module adds units, filtering, and the position
uncertainties that the Phase A model supports. The uncertainties are estimates from the validation
in ``docs/research/spectra-validation.md`` and ``docs/research/bipm-hyperfine-tables.md``, not the
output of a fit; the Phase B refit replaces them with values from its covariance.
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass
from functools import lru_cache

import numpy as np

from .constants import ATOMIC_MASS, ISOTOPOLOGUES, MHZ_PER_CM
from .intensity import C2, NEAR_FROM, intensity_model, master_line_list, with_dissociation_lines
from .model import RovibronicModel
from .spectrum import air_to_vacuum, doppler_fwhm, vacuum_to_air

#: Units accepted for a search range: vacuum nm, standard-air nm, wavenumber, frequency.
UNITS = ("nm", "nm-air", "cm-1", "MHz", "THz")
#: Highest v' covered by the data behind the published potentials (docs/research/spectra-validation.md).
V_UPPER_FITTED = 43
#: Highest X level the published potentials were fitted to. Above it they extrapolate, and the
#: measured errors are enormous: +7.76 cm-1 at v'' = 48 (matyugin2012) and -14.64 and -21.12 cm-1 at
#: v'' = 53 and 54 (nesterenko2019), against residuals under 4 MHz for every v'' <= 17.
V_LOWER_FITTED = 17
#: From this v' up, in the BIPM range, the B levels drift by more than the 3 MHz quoted below it: the
#: BIPM 532 nm intervals put v' = 33-37 at -6 to +6 MHz against v' = 32, and Nishiyama 2024 puts v' = 39
#: at 5.5 MHz (docs/design/uncertainty.md).
V_UPPER_DRIFT = 31
#: Wavenumber above which the BIPM tables constrain the model (500-667 nm).
NU_PRECISION = 15000.0
#: 755 nm: the short-wavelength end of the near-infrared anchor data (Liao et al. 2010) and of the
#: band set (v' = 0, v'' = 12-17) that the published local models describe.
NU_NIR = 13250.0
#: 815 nm: the long-wavelength end of those anchors and of both local models. Below it nothing has
#: ever been measured absolutely (docs/research/nir-model-2010.md).
NU_NIR_ANCHORS = 12270.0
_LABEL = re.compile(r"^([PR])\(?(\d+)\)?\s+(\d+)\s*-\s*(\d+)$", re.IGNORECASE)
_QUANTITY = re.compile(r"^\s*([-+0-9.eE]+)\s*([A-Za-z0-9/^-]*)\s*$")
#: One end of a range: a number and, optionally, its unit. A unit starts with a letter, so the "-" between
#: two numbers is not read as part of one; units such as "cm-1" and "nm-air" keep theirs.
_RANGE_END = re.compile(r"(\d+\.?\d*(?:[eE][-+]?\d+)?|\.\d+(?:[eE][-+]?\d+)?)\s*([A-Za-z]+(?:\^?-1|-air)?)?")


def to_wavenumber(value, unit="nm"):
    """Vacuum wavenumber (cm⁻¹) from a value in any of UNITS."""
    if unit == "nm":
        return 1e7 / float(value)
    if unit == "nm-air":
        return 1e7 / float(air_to_vacuum(value))
    if unit == "cm-1":
        return float(value)
    if unit == "MHz":
        return float(value) / MHZ_PER_CM
    if unit == "THz":
        return float(value) * 1e6 / MHZ_PER_CM
    raise ValueError(f"unknown unit {unit!r}; choose from {UNITS}")


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
    if unit == "GHz":
        return float(m[1]) * 1e3, "MHz"
    if unit not in UNITS:
        raise ValueError(f"unknown unit {unit!r}; choose from {UNITS}")
    return float(m[1]), unit


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
        a, b = ISOTOPOLOGUES[self.isotopologue]
        return float(doppler_fwhm(self.nu, self.temperature, ATOMIC_MASS[a] + ATOMIC_MASS[b])) * MHZ_PER_CM

    @property
    def uncertainty_MHz(self):
        return uncertainty(self)[0]

    @property
    def flags(self):
        out = []
        if self.v_lower > V_LOWER_FITTED:
            if 48 <= self.v_lower <= 54:
                out.append("v″ ≥ 48: extended model")
            elif self.v_lower <= 25 and corrected(self) is not None:
                out.append("v″ 18-25: atlas-measured")
            else:
                out.append("v″ > 17: unmeasured" if self.v_lower <= 25 or self.v_lower > 89 else "v″ 26-89: Martin 1986")
        if self.v_upper > V_UPPER_FITTED:
            if self.v_upper <= 50:
                out.append("v′ > 43: extended model")
            else:
                out.append("v′ 51-79: atlas-measured" if corrected(self) is not None else "v′ > 50: extrapolated")
        if self.isotopologue != "127I2":
            out.append("isotope shift")
        if corrected(self) is not None:
            out.append("corrected levels")
        elif self.nu < NU_NIR_ANCHORS:
            out.append("beyond 815 nm")
        elif self.nu < NU_NIR:
            out.append("local NIR model" if local_nir_covers(self) else "few data")
        elif self.nu < NU_PRECISION:
            out.append("few data")
        return tuple(out)


@lru_cache(maxsize=1)
def _local_nir():
    """The Knöckel 2004 local Dunham model, loaded once (i2spec.local_nir)."""
    from .local_nir import load_local_nir
    return load_local_nir("knockel2004")


def local_nir_covers(line: Line) -> bool:
    """Does a published local NIR model describe this line better than the potentials do?

    True inside the Knöckel 2004 domain (¹²⁷I₂, v' = 0, v'' = 12-17, J'' <= 242) and its validated
    775-815 nm window, where it reproduces the measured anchors to 0.23 MHz rms against the
    potentials' 2.31 MHz (docs/research/nir-model-2010.md §6).
    """
    m = _local_nir()
    return (line.isotopologue == "127I2" and m.covers(line.v_upper, line.v_lower, line.J_lower)
            and m.domain.covers_wavenumber(line.nu))


@lru_cache(maxsize=1)
def _corrections():
    """The default parameter set's level corrections, if it names any."""
    from .constants import DEFAULT_PARAMETERS
    from .level_corrections import corrections_for
    return corrections_for(DEFAULT_PARAMETERS)


def corrected(line: Line):
    """Is this ¹²⁷I₂ line inside the measured range of the parameter set's level corrections?

    Returns the estimated 1σ (MHz) of its corrected position, or None. X levels below v″ = 11 and
    B v′ = 0 are the reference the corrections are defined against, so they count as covered.
    """
    c = _corrections()
    if c is None or line.isotopologue != "127I2":
        return None
    if c.band_covers(line.v_upper, line.v_lower, line.J_lower):
        return round(c.band["held_out_MHz"], 2)      # the whole prediction, held out a line at a time
    J_up = line.J_lower + (1 if line.branch == "R" else -1)
    x_ok = line.v_lower <= 10 or c.covers("X", line.v_lower, line.J_lower)
    b_ok = line.v_upper == 0 or c.covers("B", line.v_upper, J_up)
    if not (x_ok and b_ok):
        return None
    parts = [0.3]                                   # the floor the corrections were fitted with
    for state, v, J in (("X", line.v_lower, line.J_lower), ("B", line.v_upper, J_up)):
        if (state == "X" and v <= 10) or (state == "B" and v == 0):
            continue
        # the covariance propagated to this J where the set carries it (it grows away from the data),
        # else the level's held-out rms
        u = c.uncertainty_at(state, v, J)
        if u is None:
            u = c.uncertainty(state, v)
        parts.append(1.0 if u is None else u)       # a level fixed by one line: exact there, 1 MHz assumed nearby
    return round(math.sqrt(sum(p * p for p in parts)), 2)


def gp_corrected(line: Line):
    """1σ (MHz) of a line whose B level the Gaussian-process correction fills in, to X v″ <= 5 (reference levels
    the comb data pin down); None otherwise, or where the region value is smaller. Lines to other
    X levels keep the X level's rule: their position carries the correction, their uncertainty the X level's."""
    c = _corrections()
    if c is None or line.isotopologue != "127I2":
        return None
    J_up = line.J_lower + (1 if line.branch == "R" else -1)
    if not c.gp_covers("B", line.v_upper):
        return None
    if line.v_lower > 5:           # beyond v'' = 5 the reference X levels themselves are uncertain by 1-6 MHz
        return None                # (the bake-off's X surface at v'' = 5-8, J up to 150), so the region rule stays
    _, sd = c.gp(line.v_upper, J_up)
    u = math.hypot(sd, 0.5)
    return round(u, 2) if u < 3.0 else None


def uncertainty(line: Line):
    """Estimated 1σ uncertainty (MHz) of the hyperfine-free position, and the evidence for it.

    This is the uncertainty of the position i2spec reports, which always comes from the potentials.
    Where ``local_nir_covers`` is true, ``i2spec.local_nir`` predicts the same line an order of
    magnitude better; the basis string says so.
    """
    c = _corrections()
    if (line.isotopologue == "127I2" and c is not None and c.band_covers(line.v_upper, line.v_lower, line.J_lower)):
        return corrected(line), (
            "NIR band correction: the v′ = 0 → v″ = 12-17 lines on top of the level corrections, "
            "fitted to 108 comb- and beat-referenced rows of Liao 2010, Bodermann 1998/2000 and Bodermann's thesis "
            "(liao2010a 0.027 MHz in-sample with its 114 kHz pressure shift, 0.053 held out); this line is within "
            "15 in J″ of the data at its v″, and the value is the held-out rms of all 108 rows")
    if (u := corrected(line)) is not None and 18 <= line.v_lower <= 25:
        return u, ("X v″ = 18-25 measured: 3 700 lines of the Orsay atlas part I (Gerstenkorn, Vergès & "
                   "Chevillard 1982, cell at 790 °C) assigned and fitted with the atlas's own scale calibration; this line "
                   "is within 15 in J of the atlas lines that fixed its lower level, and the value is that level's "
                   "held-out uncertainty with the atlas's per-line scatter removed")
    if u is not None and line.v_upper > 50:
        return u, ("B v′ = 51-79 measured: 2 009 unblended lines of the Orsay atlas Partie IV (Gerstenkorn "
                   "& Luc 1983), on the scale of the three atlas lines also measured against a comb, fitted level by level "
                   "on the refitted B potential; this line is within 15 in J′ of the atlas lines of its upper level, and "
                   "the value is that level's held-out uncertainty with the atlas's per-line scatter removed")
    if u is not None:
        return u, ("measured level corrections: every comb-referenced line fitted as polynomials in "
                   "J(J+1) per level, 0.02-0.3 MHz in-sample; this line is within 15 in J of the lines that fixed its "
                   "levels, and the value is their held-out rms (1 MHz where a single line fixed a level)")
    if (g := gp_corrected(line)) is not None:
        return g, ("Gaussian-process B correction: this B level has no correction of its own, so its value "
                   "is interpolated from the corrected levels around it in v′ and J′ (docs/research/model-bakeoff.md), and the lower level is X v″ <= 5; "
                   "held out a level at a time, B v′ = 3-35 is predicted to 1.1 MHz rms, against 1.9 MHz uncorrected, "
                   "82 % of lines within the quoted 1σ (the posterior σ, with a 0.5 MHz floor)")
    if 48 <= line.v_lower <= 54:
        return 20.0, ("v″ = 48-54 comes from the extended-range MLR X potential: the emission lines of "
                      "matyugin2012 and nesterenko2019 are reproduced to 7 and 19 MHz, and that residual is the X-state "
                      "hyperfine model there, not the levels")
    if line.v_lower > V_LOWER_FITTED:
        if line.v_lower <= 25:
            return 300.0, ("v″ = 18-25 outside the J the Orsay atlas part I measured: the extended X potential "
                           "alone, which the atlas levels hold to ~15 MHz within their J")
        if line.v_lower <= 75:
            return 150.0, ("the extended X potential fitted to the level constants of "
                           "Martin et al. 1986: 60 MHz rms at v″ = 26-47, 55 at 49-60, 73 at 61-75 for J ≤ 120, "
                           "twice that quoted for the constants' own accuracy and J extrapolation")
        if line.v_lower <= 89:
            return 400.0, ("the extended X potential fitted to Martin et al. 1986 at v″ = 76-89, "
                           "174 MHz rms; near the X limit, where the levels crowd")
        return 2.1e6, ("v″ > 89: not physical. The extended X potential is fitted to v″ = 89 and dips up to "
                       "170 cm⁻¹ below its dispersion limit at 7-9 Å; its levels above sit 70 cm⁻¹ rms from Martin "
                       "et al. 1986, and fits that follow those cost the measured levels below "
                       "(docs/research/x-levels-martin1986.md)")
    if line.v_upper > 50:
        return 1000.0, ("v′ > 50 outside the J the Orsay atlas Partie IV measured, or v′ = 80-86, which it did not reach: "
                        "the refitted B potential alone sits 150-350 MHz from the atlas levels, and the "
                        "level corrections carry that over only in their constant and J(J+1) terms")
    if line.v_upper > V_UPPER_FITTED:
        return 15.0, ("v′ = 44-50 comes from the extended-range MLR B potential: 12.5 MHz rms on the six "
                      "matsunaga2024a lines and 3.3 on yoshiki2023a before their corrections")
    if line.nu >= NU_PRECISION and line.v_upper >= V_UPPER_DRIFT:
        value, why = 5.0, ("the B levels drift across v′ = 31-43: the BIPM 532 nm intervals put v′ = 33-37 at −6 to "
                           "+6 MHz against v′ = 32, and the v′ = 39 lines of Nishiyama 2024 sit 5.5 MHz off")
    elif line.nu >= NU_PRECISION:
        value, why = 3.0, ("BIPM absolute frequencies (±2 MHz) at 500-667 nm; 672 comb-referenced lines at v′ ≤ 17 "
                           "are reproduced to 1.2 MHz rms")
    elif line.nu >= NU_NIR:
        value, why = 5.0, ("the 671 nm anchor of Huang et al. 2013: the model sits +2.3 MHz above it, and within "
                           "0.33 MHz of IodineSpec 5")
    elif line.nu >= NU_NIR_ANCHORS:
        if line.v_upper == 0:
            value, why = 3.0, ("56 comb-referenced components of the v′ = 0 bands at 755-815 nm (Liao et al. 2010, "
                               "Bodermann et al. 2000): the potentials are 2.4 MHz rms, 6.2 MHz at worst. Above "
                               "775 nm i2spec.local_nir predicts these lines to 0.23 MHz")
        else:
            value, why = 8.0, ("v′ > 0 is the worst-known part of the near infrared: the 1-14 band sits 7.0 MHz off "
                               "in Bodermann et al. 2000, and no published local model covers v′ > 0")
    else:
        value, why = 50.0, ("beyond 815 nm nothing has been measured absolutely, both published local NIR models (Knöckel 2004, Liao 2010) stop, "
                            "and the potentials are unchecked")
    if line.isotopologue != "127I2":
        return math.hypot(value, 3.0), why + ("; with the model's adiabatic correction the ¹²⁹I₂ lines are 3.4 MHz rms "
                                              "and ¹²⁷I¹²⁹I 0.3 MHz against BIPM (isotope shifts)")
    return value, why


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

    def master(self, isotopologue="127I2", near_dissociation=False):
        """The master line list; with ``near_dissociation``, extended to the B limit (a further build, once)."""
        if isotopologue not in self._masters:
            model = intensity_model(isotopologue)
            self._masters[isotopologue] = master_line_list(model, *self.nu_range, T_range=(200.0, 600.0),
                                                           S_min=self.s_min)
        if near_dissociation and self.nu_range[1] > NEAR_FROM and isotopologue not in self._extended:
            self._masters[isotopologue] = with_dissociation_lines(self._masters[isotopologue], S_min=self.s_min,
                                                                  nu_max=self.nu_range[1])
            self._extended.add(isotopologue)
        return self._masters[isotopologue]

    def model(self, isotopologue="127I2"):
        """The B-spline model used for hyperfine structure (more accurate than the intensity grid)."""
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
        model = self.model(isotopologue)
        nu = model.transition(v_upper, v_lower, J, branch)
        master = self.master(isotopologue, near_dissociation=v_upper > 60)
        sel = ((master.v_upper == v_upper) & (master.v_lower == v_lower) & (master.J_lower == J)
               & (master.branch == (1 if branch == "R" else -1)))
        if not sel.any():
            raise LookupError(f"{label} of {isotopologue} is not in the line list "
                              f"(it would lie at {1e7 / nu:.4f} nm)")
        k = int(np.flatnonzero(sel)[0])
        strength = master.strength0[k] * math.exp(-C2 * master.E_lower[k] / T) / master.partition_function(T)
        return Line(isotopologue, branch, J, v_upper, v_lower, float(master.nu[k]), float(strength),
                    float(master.E_lower[k]), T)

    def components(self, line: Line, dJ=2, main_only=False) -> list[Component]:
        """Hyperfine components of a line, with absolute frequencies."""
        model = self.model(line.isotopologue)
        nu0, comps = model.hyperfine_components(line.v_upper, line.v_lower, line.J_lower, line.branch, dJ=dJ)
        out = [Component(c.label, c.offset, c.strength, c.F_upper, c.F_lower, nu0 + c.offset, c.I_upper, c.I_lower,
                         c.upper_level, c.lower_level)
               for c in comps if c.label or not main_only]
        return sorted(out, key=lambda c: c.offset_MHz)
