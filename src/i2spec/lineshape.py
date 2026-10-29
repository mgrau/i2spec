"""Measured line widths, pressure coefficients and level lifetimes: data/widths/ and data/lifetimes/.

Kept apart from data/observations (frequencies) because the model does not fit them yet; they are the test
data for its line strengths and widths. Each set is a directory <id>/ holding meta.toml (source, provenance,
conditions) and data.csv. Units are fixed per quantity so that sets can be compared without reading their
meta: a source's own convention (HWHM, angular frequency, MHz/Torr, ...) is converted on transcription and
the printed value kept in the row's note.
"""

from __future__ import annotations

import csv
import math
import tomllib
from dataclasses import dataclass
from pathlib import Path

from .constants import ISOTOPOLOGUES
from .observations import Line

ROOT = Path(__file__).resolve().parents[2] / "data"
WIDTHS_DIR = ROOT / "widths"
LIFETIMES_DIR = ROOT / "lifetimes"
META_REQUIRED = ("id", "catalog", "source", "retrieved", "transcription", "uncertainty")

#: quantity -> unit. Widths are full widths at half maximum in linear frequency.
WIDTH_UNITS = {
    "width": "kHz",          # FWHM at pressure_Pa (empty: extrapolated to zero pressure, as the source states)
    "broadening": "kHz/Pa",  # d(FWHM)/dp
    "shift": "kHz/Pa",       # d(frequency)/dp
}
WIDTH_COLUMNS = ("line", "component", "quantity", "value", "uncertainty", "pressure_Pa", "temperature_K",
                 "perturber", "group", "note")

#: quantity -> unit, for one level (or one hyperfine sublevel) of a state.
LIFETIME_UNITS = {
    "lifetime": "ns",                    # collision-free (zero-pressure) lifetime
    "decay_rate": "s-1",                 # total, 1/lifetime
    "radiative_rate": "s-1",
    "predissociation_rate": "s-1",       # natural (spontaneous), field-free
    "hyperfine_predissociation_a2": "s-1",   # |a_v|^2 of the magnetic-hyperfine predissociation (Broyer, Vigue)
    "gyroscopic_predissociation_c2": "s-1",  # |c_v|^2: Gamma_rot = |c_v|^2 J(J+1) (natural, rotation-induced)
    "quenching_cross_section": "A2",     # self-quenching unless the note names another gas
}
LIFETIME_COLUMNS = ("isotopologue", "state", "v", "J", "sublevel", "quantity", "value", "uncertainty", "group",
                    "note")


@dataclass(frozen=True)
class Width:
    line: Line
    component: str | None
    quantity: str
    value: float
    uncertainty: float
    pressure_Pa: float | None
    temperature_K: float | None
    perturber: str
    group: str
    note: str

    @property
    def unit(self):
        return WIDTH_UNITS[self.quantity]


@dataclass(frozen=True)
class LevelDecay:
    isotopologue: str
    state: str
    v: int
    J: int | None
    sublevel: str | None
    quantity: str
    value: float
    uncertainty: float
    group: str
    note: str

    @property
    def unit(self):
        return LIFETIME_UNITS[self.quantity]


@dataclass
class LineShapeSet:
    id: str
    meta: dict
    rows: list


def _optional(text, kind=float):
    return kind(text) if text.strip() else None


def _meta(path):
    meta = tomllib.loads((path / "meta.toml").read_text(encoding="utf-8"))
    missing = [k for k in META_REQUIRED if k not in meta]
    if missing:
        raise ValueError(f"{path / 'meta.toml'}: missing {', '.join(missing)}")
    if meta["id"] != path.name:
        raise ValueError(f"{path / 'meta.toml'}: id {meta['id']!r} differs from the directory name")
    return meta


def _rows(path, columns, parse):
    rows = []
    with open(path / "data.csv", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        if tuple(reader.fieldnames or ()) != columns:
            raise ValueError(f"{path / 'data.csv'}: columns must be {','.join(columns)}")
        for n, row in enumerate(reader, start=2):
            try:
                r = parse(row)
                if not (math.isfinite(r.value) and r.uncertainty > 0):
                    raise ValueError("value must be finite and uncertainty positive")
                rows.append(r)
            except (ValueError, KeyError) as e:
                raise ValueError(f"{path / 'data.csv'}, line {n}: {e}") from None
    return rows


def _width(row):
    if row["quantity"] not in WIDTH_UNITS:
        raise ValueError(f"quantity must be one of {tuple(WIDTH_UNITS)}, not {row['quantity']!r}")
    return Width(Line.parse(row["line"]), row["component"] or None, row["quantity"], float(row["value"]),
                 float(row["uncertainty"]), _optional(row["pressure_Pa"]), _optional(row["temperature_K"]),
                 row["perturber"] or "self", row["group"], row["note"])


def _decay(row):
    if row["quantity"] not in LIFETIME_UNITS:
        raise ValueError(f"quantity must be one of {tuple(LIFETIME_UNITS)}, not {row['quantity']!r}")
    if row["isotopologue"] not in ISOTOPOLOGUES:
        raise ValueError(f"unknown isotopologue {row['isotopologue']!r}")
    return LevelDecay(row["isotopologue"], row["state"], int(row["v"]), _optional(row["J"], int),
                      row["sublevel"] or None, row["quantity"], float(row["value"]), float(row["uncertainty"]),
                      row["group"], row["note"])


def load_widths(path) -> LineShapeSet:
    """A data/widths/<id>/ directory."""
    path = Path(path)
    meta = _meta(path)
    return LineShapeSet(meta["id"], meta, _rows(path, WIDTH_COLUMNS, _width))


def load_lifetimes(path) -> LineShapeSet:
    """A data/lifetimes/<id>/ directory."""
    path = Path(path)
    meta = _meta(path)
    return LineShapeSet(meta["id"], meta, _rows(path, LIFETIME_COLUMNS, _decay))


def _load_all(root, load, include_excluded):
    root = Path(root)
    if not root.exists():
        return []
    sets = [load(p) for p in sorted(root.iterdir()) if (p / "meta.toml").exists()]
    return [s for s in sets if include_excluded or "exclude" not in s.meta]


def load_all_widths(root=WIDTHS_DIR, include_excluded=False) -> list[LineShapeSet]:
    """Every set under ``root``, except those whose meta.toml has an ``exclude`` reason."""
    return _load_all(root, load_widths, include_excluded)


def load_all_lifetimes(root=LIFETIMES_DIR, include_excluded=False) -> list[LineShapeSet]:
    """Every set under ``root``, except those whose meta.toml has an ``exclude`` reason."""
    return _load_all(root, load_lifetimes, include_excluded)
