"""Measured line positions and hyperfine intervals: the data sets of the global fit.

A data set is a directory data/observations/<id>/ holding meta.toml (source, provenance, unit,
conditions) and data.csv (one observation per row); docs/design/observations.md describes the
format. Lines are keyed as "127I2 R(56) 32-0". Hyperfine components keep their published labels
("a10", "b12"): the number is the component's rank by frequency (1 = lowest), and the letter only
names the line in the source.
"""

from __future__ import annotations

import csv
import math
import re
import tomllib
from dataclasses import dataclass, replace
from pathlib import Path

import numpy as np

from .constants import ISOTOPOLOGUES, MHZ_PER_CM, DEFAULT_PARAMETERS
from .model import RovibronicModel

KINDS = ("frequency", "interval")
UNITS = {"MHz": 1.0, "cm-1": MHZ_PER_CM}  # MHz per unit
COLUMNS = ("line", "component", "kind", "value", "uncertainty", "ref_line", "ref_component", "group", "note")
META_REQUIRED = ("id", "citation", "source", "retrieved", "unit", "transcription")
DATA_DIR = Path(__file__).resolve().parents[2] / "data" / "observations"
PARAMETERS_DIR = Path(__file__).resolve().parents[2] / "data" / "hyperfine_parameters"
PARAMETER_COLUMNS = ("line", "eqQ", "eqQ_unc", "C", "C_unc", "d", "d_unc", "delta", "delta_unc", "fit_sd_kHz", "note")
PARAMETER_META_REQUIRED = ("id", "citation", "source", "retrieved", "unit", "transcription")
_LINE = re.compile(r"^(\S+) ([PR])\((\d+)\) (\d+)-(\d+)$")
_COMPONENT = re.compile(r"^[a-z](\d+)$")


@dataclass(frozen=True)
class Line:
    isotopologue: str
    branch: str
    J_lower: int
    v_upper: int
    v_lower: int

    @classmethod
    def parse(cls, text):
        m = _LINE.match(text.strip())
        if not m or m[1] not in ISOTOPOLOGUES:
            raise ValueError(f"bad line {text!r}; expected e.g. '127I2 R(56) 32-0'")
        if m[2] == "P" and int(m[3]) == 0:
            raise ValueError(f"{text!r}: there is no P(0) line")
        return cls(m[1], m[2], int(m[3]), int(m[4]), int(m[5]))

    def __str__(self):
        return f"{self.isotopologue} {self.branch}({self.J_lower}) {self.v_upper}-{self.v_lower}"


def component_rank(label):
    """1 for "a1", 12 for "b12", ...: the rank of a hyperfine component by frequency."""
    m = _COMPONENT.match(label)
    if not m:
        raise ValueError(f"bad component label {label!r}; expected a letter and a number, e.g. 'a10'")
    return int(m[1])


@dataclass(frozen=True)
class Observation:
    """kind "frequency": f(line, component). kind "interval": f(line, component) − f(ref_line, ref_component),
    where ref_line defaults to line. A missing component means the hyperfine-free line centre."""

    line: Line
    component: str | None
    kind: str
    value: float
    uncertainty: float
    ref_line: Line | None = None
    ref_component: str | None = None
    group: str | None = None  # observations that share a nuisance parameter, e.g. a calibration offset
    note: str = ""


@dataclass
class Dataset:
    id: str
    meta: dict
    observations: list[Observation]

    @property
    def unit(self):
        return self.meta["unit"]

    def __len__(self):
        return len(self.observations)


def _observation(row):
    kind = row["kind"]
    if kind not in KINDS:
        raise ValueError(f"kind must be one of {KINDS}, not {kind!r}")
    component, ref_component = row["component"] or None, row["ref_component"] or None
    for label in (component, ref_component):
        if label is not None:
            component_rank(label)
    ref_line = Line.parse(row["ref_line"]) if row["ref_line"] else None
    if kind == "interval" and ref_line is None and ref_component is None:
        raise ValueError("an interval needs ref_line or ref_component")
    if kind == "frequency" and (ref_line or ref_component):
        raise ValueError("a frequency has no reference")
    uncertainty = float(row["uncertainty"])
    if not uncertainty > 0:
        raise ValueError("uncertainty must be positive")
    return Observation(Line.parse(row["line"]), component, kind, float(row["value"]), uncertainty, ref_line,
                       ref_component, row["group"] or None, row["note"])


def read_observations(path):
    with open(path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        if tuple(reader.fieldnames or ()) != COLUMNS:
            raise ValueError(f"{path}: columns must be {','.join(COLUMNS)}")
        observations = []
        for n, row in enumerate(reader, start=2):
            try:
                observations.append(_observation(row))
            except ValueError as e:
                raise ValueError(f"{path}, line {n}: {e}") from None
    return observations


def write_observations(observations, path):
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f, lineterminator="\n")
        writer.writerow(COLUMNS)
        for o in observations:
            writer.writerow([o.line, o.component or "", o.kind, repr(o.value), repr(o.uncertainty), o.ref_line or "",
                             o.ref_component or "", o.group or "", o.note])


#: [shift_correction] keys: the cell pressure, pump power and modulation width (peak-to-peak) a set's absolute
#: frequencies were measured at, and the source's own shift coefficients, which take them to zero pressure,
#: power and modulation; allowance_kHz is a
#: type-B term added in quadrature where the source gives no coefficient; groups, if given, limits the
#: correction to the rows of those groups (one set can hold lines measured in different cells).
SHIFT_KEYS = ("pressure_Pa", "pressure_kHz_per_Pa", "pressure_kHz_per_Pa_unc", "power_mW", "power_kHz_per_mW",
              "power_kHz_per_mW_unc", "modulation_MHz", "modulation_kHz_per_MHz", "modulation_kHz_per_MHz_unc",
              "allowance_kHz", "groups", "basis")


def shift_correction(meta):
    """(shift to add, uncertainty to add in quadrature), in kHz, for the absolute frequencies of a set.

    Pressure and power shifts are negative for I2, so a zero-pressure value is the measured one minus
    slope x pressure. Intervals within a line are left alone: the shift is common to its components to
    well under their uncertainty (docs/design/observations.md, "Pressure and power")."""
    c = meta.get("shift_correction")
    if not c:
        return 0.0, 0.0
    unknown = set(c) - set(SHIFT_KEYS)
    if unknown:
        raise ValueError(f"{meta['id']}: unknown shift_correction keys {sorted(unknown)}")
    terms = [(float(c.get(f"{q}_{u}", 0.0)), float(c.get(f"{q}_kHz_per_{u}", 0.0)), float(c.get(f"{q}_kHz_per_{u}_unc", 0.0)))
             for q, u in (("pressure", "Pa"), ("power", "mW"), ("modulation", "MHz"))]
    shift = -sum(x * slope for x, slope, _ in terms)
    unc = math.sqrt(sum((x * du) ** 2 for x, _, du in terms) + float(c.get("allowance_kHz", 0.0)) ** 2)
    return shift, unc


def _corrected(observations, meta):
    shift, unc = shift_correction(meta)
    if not (shift or unc):
        return observations
    k = 1e-3 / UNITS[meta["unit"]]                  # kHz in the data set's unit
    groups = meta["shift_correction"].get("groups")
    return [replace(o, value=o.value + shift * k, uncertainty=math.hypot(o.uncertainty, unc * k))
            if o.kind == "frequency" and (groups is None or o.group in groups) else o for o in observations]


def load_dataset(path, raw=False) -> Dataset:
    """A data set with its absolute frequencies taken to zero pressure and power ([shift_correction]);
    ``raw`` keeps them as printed."""
    path = Path(path)
    meta = tomllib.loads((path / "meta.toml").read_text(encoding="utf-8"))
    missing = [k for k in META_REQUIRED if k not in meta]
    if missing:
        raise ValueError(f"{path / 'meta.toml'}: missing {', '.join(missing)}")
    if meta["id"] != path.name:
        raise ValueError(f"{path / 'meta.toml'}: id {meta['id']!r} differs from the directory name")
    if meta["unit"] not in UNITS:
        raise ValueError(f"{path / 'meta.toml'}: unit must be one of {tuple(UNITS)}")
    observations = read_observations(path / "data.csv")
    return Dataset(meta["id"], meta, observations if raw else _corrected(observations, meta))


def _what(o):
    """What an observation measures, for defer_to: an absolute component, or an interval between two lines
    (or within one), whatever the components."""
    if o.kind == "frequency":
        return ("frequency", o.line, o.component)
    return ("interval", o.line, o.ref_line or o.line)


def load_all(root=DATA_DIR, include_excluded=False, raw=False) -> list[Dataset]:
    """Every data set under ``root``, except those whose meta.toml has an ``exclude`` reason.

    A set whose meta.toml lists ``defer_to`` (a compilation such as the BIPM tables) loses the rows that one
    of those sets also measures: the same absolute component, or intervals within the same line or between
    the same two lines. Each measurement then enters once, from its source."""
    datasets = [load_dataset(p, raw=raw) for p in sorted(Path(root).iterdir()) if (p / "meta.toml").exists()]
    by_id = {d.id: d for d in datasets}
    for d in datasets:
        sources = [by_id[i] for i in d.meta.get("defer_to", []) if i in by_id]
        if sources:
            held = {_what(o) for src in sources for o in src.observations}
            d.observations = [o for o in d.observations if _what(o) not in held]
    return [d for d in datasets if (include_excluded or "exclude" not in d.meta) and d.observations]


@dataclass(frozen=True)
class MeasuredHyperfine:
    """Effective B-state hyperfine parameters fitted by a source to one transition's spectrum, in MHz."""

    line: Line
    eqQ: float
    eqQ_unc: float
    C: float
    C_unc: float
    d: float
    d_unc: float
    delta: float
    delta_unc: float
    fit_sd_kHz: float | None
    note: str

    @property
    def J_upper(self):
        return self.line.J_lower + (1 if self.line.branch == "R" else -1)


@dataclass
class HyperfineParameterSet:
    id: str
    meta: dict
    rows: list[MeasuredHyperfine]


def load_hyperfine_parameters(path) -> HyperfineParameterSet:
    """A data/hyperfine_parameters/<id>/ directory: parameters a source derived, not raw frequencies.

    Kept apart from data/observations because the numbers depend on the source's own model choices
    (which X-state values it held fixed, how many terms it fitted), which meta.toml must record.
    """
    path = Path(path)
    meta = tomllib.loads((path / "meta.toml").read_text(encoding="utf-8"))
    missing = [k for k in PARAMETER_META_REQUIRED if k not in meta]
    if missing:
        raise ValueError(f"{path / 'meta.toml'}: missing {', '.join(missing)}")
    if meta["id"] != path.name:
        raise ValueError(f"{path / 'meta.toml'}: id {meta['id']!r} differs from the directory name")
    if meta["unit"] != "MHz":
        raise ValueError(f"{path / 'meta.toml'}: unit must be MHz")
    rows = []
    with open(path / "data.csv", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        if tuple(reader.fieldnames or ()) != PARAMETER_COLUMNS:
            raise ValueError(f"{path / 'data.csv'}: columns must be {','.join(PARAMETER_COLUMNS)}")
        for n, row in enumerate(reader, start=2):
            try:
                values = {k: float(row[k]) for k in PARAMETER_COLUMNS[1:9]}
                if any(values[k] <= 0 for k in PARAMETER_COLUMNS[2:9:2]):
                    raise ValueError("uncertainties must be positive")
                sd = float(row["fit_sd_kHz"]) if row["fit_sd_kHz"] else None
                rows.append(MeasuredHyperfine(Line.parse(row["line"]), **values, fit_sd_kHz=sd, note=row["note"]))
            except ValueError as e:
                raise ValueError(f"{path / 'data.csv'}, line {n}: {e}") from None
    return HyperfineParameterSet(meta["id"], meta, rows)


def load_all_hyperfine_parameters(root=PARAMETERS_DIR) -> list[HyperfineParameterSet]:
    return [load_hyperfine_parameters(p) for p in sorted(Path(root).iterdir()) if (p / "meta.toml").exists()]


class Predictor:
    """Model values of observations: one RovibronicModel per isotopologue, with cached line results."""

    def __init__(self, parameters=DEFAULT_PARAMETERS, dJ=2, **model_options):
        self.parameters, self.dJ, self.model_options = parameters, dJ, model_options
        self._models, self._centres, self._offsets = {}, {}, {}

    def model(self, isotopologue):
        if isotopologue not in self._models:
            self._models[isotopologue] = RovibronicModel(isotopologue, self.parameters, **self.model_options)
        return self._models[isotopologue]

    def position(self, line, component=None):
        """Frequency (MHz) of a hyperfine component, or of the hyperfine-free line centre."""
        if line not in self._centres:
            m = self.model(line.isotopologue)
            self._centres[line] = m.transition(line.v_upper, line.v_lower, line.J_lower, line.branch) * MHZ_PER_CM
        if component is None:
            return self._centres[line]
        if line not in self._offsets:
            _, comps = self.model(line.isotopologue).hyperfine_components(line.v_upper, line.v_lower, line.J_lower,
                                                                          line.branch, dJ=self.dJ)
            self._offsets[line] = {component_rank(c.label): c.offset for c in comps if c.label}
        return self._centres[line] + self._offsets[line][component_rank(component)]

    def __call__(self, obs, unit="MHz"):
        f = self.position(obs.line, obs.component)
        if obs.kind == "interval":
            f -= self.position(obs.ref_line or obs.line, obs.ref_component)
        return f / UNITS[unit]


def residuals(predictor, dataset):
    """Observed − model for each observation of a data set, in the data set's unit."""
    return np.array([o.value - predictor(o, dataset.unit) for o in dataset.observations])
