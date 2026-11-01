"""The training table: every measured B-X line center, interval between lines and X level energy, in MHz.

Every row is a linear functional of term values, value = sum_k sign_k E(isotopologue_k, state_k, v_k, J_k); the
functional always differences the origin away. Three kinds:

  abs       hyperfine-free line center: +E(B, v', J') - E(X, v'', J'')
  interval  difference of two line centers in one data set: (line 1) - (line 2)
  level     a measured level energy: +E(X, v, J) - E(X, 0, 0)

Categories, each with its own weight floor in the fit:

  comb      precision sets whose median sigma is below 0.3 MHz (comb- or standard-referenced)
  doppler   precision sets with median sigma >= 0.3 MHz (etalon or wavemeter, or large unmeasured shifts)
  emission  comb-referenced B -> X emission lines (X v'' = 48-54)
  atlas_sr, atlas_apo   the two refitted Fourier-transform atlases, on scales calibrated against the precision lines
  orsay1, orsay4        the Orsay atlas, part I and Partie IV, assigned and unblended lines
  martin, martin_dissoc X level constants of Martin et al., evaluated at J = 0, 20, ...; v'' >= 90 is martin_dissoc

`new` marks the 29 precision sets transcribed last, held out of the core fit as its benchmark. `nu` is an
approximate wavenumber, a label only: no fit uses it.
"""
from __future__ import annotations

import csv
import gzip
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
TRAINING = ROOT / "data" / "training"
RELEASE = "20261021"

#: the precision sets transcribed last: the core fit's benchmark, held out of every fit
NEW_SETS = ("huang2013a", "hong2009a", "ye1999a", "holzwarth2001a", "jones2002a", "goncharov2004a", "simonsen2000a",
            "huang2018a", "manzoor2024a", "fan2014a", "kobayashi2015a", "hauden2024a", "grieser1994a", "badr2006a",
            "hong2001b", "sakagami2020a", "yoshii2019a", "tanabe2022a", "arie1993a", "arie1994a", "cheng2001a",
            "zhang2001a", "hong2001a", "hong2002a", "hong2000a", "hong2004a", "hong2004b", "hong2002b", "fang2006a")
#: the core: X v'' <= 17, B v' <= 43
CORE_V = dict(X=17, B=43)
#: the line left out of every fit and headline score (a diagnostic target only)
EXCLUSION = "127I2 P(13) 43-0"
PRECISION = ("comb", "doppler")
COLUMNS = ("id", "source", "kind", "category", "isotopologue", "branch", "J_lower", "v_upper", "v_lower",
           "isotopologue2", "branch2", "J_lower2", "v_upper2", "v_lower2", "state", "v", "J",
           "value", "sigma", "sigma_stated", "n_components", "nu", "new")
_INT = {"J_lower", "v_upper", "v_lower", "J_lower2", "v_upper2", "v_lower2", "v", "J", "n_components"}
_FLOAT = {"value", "sigma", "sigma_stated", "nu"}


class Rows:
    """Columns as numpy arrays, plus `levels`: the functional of each row as a tuple of ((isotopologue, state, v,
    J), sign). rows[mask or indices] is a sub-table."""

    def __init__(self, cols: dict, levels: list):
        self.cols = cols
        self.levels = levels

    def __len__(self):
        return len(self.levels)

    def __getattr__(self, name):
        cols = self.__dict__.get("cols")
        if cols is not None and name in cols:
            return cols[name]
        raise AttributeError(name)

    def __getitem__(self, idx):
        if isinstance(idx, str):
            return self.cols[idx]
        idx = np.arange(len(self))[idx] if not (isinstance(idx, np.ndarray) and idx.dtype.kind == "i") else idx
        return Rows({k: v[idx] for k, v in self.cols.items()}, [self.levels[i] for i in idx])

    def copy(self):
        return Rows({k: v.copy() for k, v in self.cols.items()}, list(self.levels))

    def select(self, ids):
        """The rows with these ids, in that order."""
        pos = {k: i for i, k in enumerate(self.cols["id"])}
        return self[np.array([pos[i] for i in ids], dtype=int)]

    def record(self, i):
        return {k: v[i].item() if hasattr(v[i], "item") else v[i] for k, v in self.cols.items()} | {"levels": self.levels[i]}

    def level_set(self):
        """Every (isotopologue, state, v, J) the rows touch."""
        return {k for lv in self.levels for k, _ in lv}


def merge(levels):
    """One signed term per level, sorted: the canonical form of a functional."""
    acc = {}
    for k, s in levels:
        acc[k] = acc.get(k, 0.) + s
    return tuple(sorted(((k, s) for k, s in acc.items() if s != 0), key=str))


def line_levels(isotopologue, branch, J, v_upper, v_lower, sign=1.):
    J_upper = J + (1 if branch == "R" else -1)
    return [((isotopologue, "B", v_upper, J_upper), sign), ((isotopologue, "X", v_lower, J), -sign)]


def _levels_text(levels):
    return ";".join(f"{iso} {state} {v} {J} {sign:g}" for (iso, state, v, J), sign in levels)


def _parse_levels(text):
    out = []
    for term in text.split(";"):
        iso, state, v, J, sign = term.split()
        out.append(((iso, state, int(v), int(J)), float(sign)))
    return tuple(out)


def write(rows: Rows, path: Path, extra: dict | None = None):
    """rows as gzipped CSV, floats to full precision; `extra` adds named columns (one value per row)."""
    extra = extra or {}
    with gzip.open(path, "wt", newline="", encoding="utf-8") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(COLUMNS + ("levels",) + tuple(extra))
        for i in range(len(rows)):
            out = []
            for k in COLUMNS:
                z = rows.cols[k][i]
                out.append(repr(float(z)) if k in _FLOAT else int(bool(z)) if k == "new" else
                           int(z) if k in _INT else z)
            w.writerow(out + [_levels_text(rows.levels[i])] + [extra[k][i] for k in extra])


def read(path: Path) -> tuple[Rows, dict]:
    """(rows, extra columns) of a table written by `write`."""
    with gzip.open(path, "rt", newline="", encoding="utf-8") as f:
        records = list(csv.DictReader(f))
    cols = {}
    for k in COLUMNS:
        z = [r[k] for r in records]
        cols[k] = (np.array(z, float) if k in _FLOAT else np.array([bool(int(a)) for a in z]) if k == "new" else
                   np.array(z, int) if k in _INT else np.array(z, dtype=str))
    extra = {k: np.array([r[k] for r in records], dtype=str)
             for k in (records[0] if records else {}) if k not in COLUMNS and k != "levels"}
    return Rows(cols, [_parse_levels(r["levels"]) for r in records]), extra


def release(stage: str, role: str = "train", name: str = RELEASE) -> Rows:
    """The rows a stage of a release was fitted to ('train') or scored on ('benchmark'), in the fit's order.

    stage: 'core' (X v''<=17, B v'<=43; the 29 newest precision sets are its benchmark) or 'continuation' (every
    127I2 row up to X v''=108, B v'=80, with the same benchmark held out)."""
    rows, extra = read(TRAINING / name / "rows.csv.gz")
    order = extra[stage]
    keep = np.flatnonzero(np.char.startswith(order, role + ":"))
    position = np.array([int(z.split(":")[1]) for z in order[keep]])
    return rows[keep[np.argsort(position)]]
