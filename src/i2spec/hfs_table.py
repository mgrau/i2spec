"""Measured B-state hyperfine parameters, as a correction to the published formulae (lever 3).

Stage 2 (docs/design/hyperfine-fit.md) found that the BKT02/S06 formulae miss the measured hyperfine
parameters in a way that is smooth in J' within one v' but steps from one v' to the next, so no smooth
global correction could be validated. What does validate is a table: at each v' where lines have been
measured, fit the difference (measured - formula) as a low-order polynomial in y = J'(J'+1)/10^4 and
apply it; at unmeasured v' keep the formula, except above v' = 53, where the formula is frozen and
interpolating the corrections of the neighbouring measured v' beats it twenty-fold (Chen 2004 test).

The table itself is src/i2spec/data/b_state_lines.json (shipped with the package), built by
prototypes/hfs_measured_table.py from Chen 2004 and per-line fits of every precise set we hold.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

import numpy as np

from .hyperfine import HyperfineParameters

PARAMS = ("eqQ", "C", "d", "delta")
TABLE_PATH = Path(__file__).resolve().parent / "data" / "b_state_lines.json"   # shipped with the package
#: Above this v' the published formulae are held at a fixed energy (hfs_params.E_B_MAX_S06), and the
#: corrections of neighbouring measured v' are interpolated instead of falling back to zero.
V_FROZEN = 53
#: Floors on the parameter uncertainties used as weights, so one 1 kHz line does not fix the polynomial.
FLOOR = {"eqQ": 0.001, "C": 0.02, "d": 0.5, "delta": 0.5}
#: Chen 2004 quotes eqQ_B to 10-25 kHz, but its values scatter by 0.1-0.2 MHz from one J' to the next
#: at the same v' (at v' = 43: -105, +99, +336, +368 kHz from the formula at J' = 24, 28, 37, 41), so
#: weighted at face value they would pull the correction away from the kHz-level BIPM lines. Its eqQ
#: rows get this floor instead; its C, d and delta are consistent with their quoted errors.
SOURCE_FLOOR = {"chen2004a": {"eqQ": 0.15}}


def _y(J):
    return J * (J + 1) / 1e4


@dataclass
class VibrationalCorrection:
    """Per-parameter polynomial corrections in y = J'(J'+1)/10^4 for one v'."""

    v: int
    J: tuple[int, ...]
    coefficients: dict[str, np.ndarray]      # highest degree last, as numpy.polynomial expects

    def __call__(self, param, J):
        # Outside the measured J' range the polynomial is held at its edge value: a line fitted to
        # J' = 52-59 says nothing about J' = 150, and the formula's own J dependence is kept there.
        y = float(np.clip(_y(J), _y(min(self.J)), _y(max(self.J))))
        return float(np.polynomial.polynomial.polyval(y, self.coefficients[param]))


class HyperfineTable:
    def __init__(self, rows, exclude=()):
        """``exclude`` names lines ("R(56) 32-0") to leave out, for cross-validation."""
        self.rows = [r for r in rows if r["line"] not in exclude]
        by_v: dict[int, list] = {}
        for r in self.rows:
            by_v.setdefault(r["v"], []).append(r)
        self.corrections = {v: self._fit(v, lines) for v, lines in by_v.items()}
        self.measured_v = sorted(self.corrections)

    @classmethod
    def load(cls, path=TABLE_PATH, exclude=()):
        return cls(json.loads(Path(path).read_text())["rows"], exclude)

    @staticmethod
    def _fit(v, lines):
        J = np.array([r["J"] for r in lines], dtype=float)
        y = _y(J)
        span = y.max() - y.min()
        coefficients = {}
        for k in PARAMS:
            d = np.array([r["measured"][k] - r["formula"][k] for r in lines])
            floors = [SOURCE_FLOOR.get(r["source"], {}).get(k, FLOOR[k]) for r in lines]
            w = 1.0 / np.hypot([r["uncertainty"][k] for r in lines], floors)
            # degree: 0 with one line or no J spread, 1 with two or three, 2 with four or more over a
            # wide spread. Stage 2 found a straight line predicts held-out C_B to 0.03 kHz under BKT02
            # and a quadratic to 0.15 kHz under S06.
            degree = 0 if len(lines) < 2 or span < 0.01 else (1 if len(lines) < 4 or span < 0.2 else 2)
            X = np.vander(y, degree + 1, increasing=True)
            coefficients[k] = np.linalg.lstsq(X * w[:, None], d * w, rcond=None)[0]
        return VibrationalCorrection(v, tuple(int(j) for j in J), coefficients)

    def correction(self, param, v, J):
        """measured - formula for parameter ``param`` at (v', J'); 0 where nothing constrains it."""
        if v in self.corrections:
            return self.corrections[v](param, J)
        if v <= V_FROZEN:
            return 0.0
        below = [m for m in self.measured_v if m < v]
        above = [m for m in self.measured_v if m > v]
        if not below:
            return 0.0
        if not above:
            return self.corrections[below[-1]](param, J)
        v0, v1 = below[-1], above[0]
        f0, f1 = self.corrections[v0](param, J), self.corrections[v1](param, J)
        return f0 + (f1 - f0) * (v - v0) / (v1 - v0)

    def apply(self, p: HyperfineParameters, v, J) -> HyperfineParameters:
        return HyperfineParameters(p.eqQ + self.correction("eqQ", v, J), p.C + self.correction("C", v, J),
                                   p.d + self.correction("d", v, J), p.delta + self.correction("delta", v, J))


@lru_cache(maxsize=1)
def default_table() -> HyperfineTable:
    return HyperfineTable.load()
