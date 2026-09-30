"""A Gaussian-process correction for the visible B levels that have no level correction of their own.

The per-level corrections (level_corrections.py) fit a polynomial in J(J+1) to each measured level and
leave every other level at the published curve, which the lookup then quotes at a region value of
3-5 MHz. A held-out bake-off (docs/research/model-bakeoff.md) found that a Gaussian process over (v', J')
that shares information between levels predicts an unmeasured visible B level to 1.7 MHz rms, against
2.8 for no correction, with calibrated uncertainties. This module evaluates that posterior at a B level;
the file, fitted by prototypes/gp_corrections_fit.py in the gauge of the level corrections (X v'' <= 10 and
B v' = 0 fixed at zero), holds the B training levels, beta = S^T K^-1 r and M = S^T K^-1 S.
"""
from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path

import numpy as np

DATA = Path(__file__).resolve().parent / "data"


class GPCorrection:
    def __init__(self, meta: dict, v, J, beta, M):
        self.meta = meta
        self.v, self.y = np.asarray(v, float), np.asarray(J, float) * (np.asarray(J, float) + 1) / 1e4
        self.beta, self.M = np.asarray(beta, float), np.asarray(M, float)
        self.a1, self.len_v, self.len_y = meta["a1"], meta["len_v"], meta["len_y"]
        self.a2, self.len_y2 = meta["a2"], meta["len_y2"]
        self.v_min, self.v_max = int(meta.get("v_min", 3)), int(meta.get("v_max", 35))

    @classmethod
    def load(cls, name: str) -> "GPCorrection":
        meta = json.loads((DATA / f"{name}.json").read_text())
        with np.load(DATA / f"{name}.npz") as f:
            return cls(meta, f["v"], f["J"], f["beta"], f["M"])

    def _k(self, v, y):
        q = np.sqrt(5.0) * np.abs(v - self.v) / self.len_v
        k = self.a1**2 * (1 + q + q * q / 3) * np.exp(-q) * np.exp(-(y - self.y) ** 2 / (2 * self.len_y**2))
        return k + self.a2**2 * (self.v == v) * np.exp(-(y - self.y) ** 2 / (2 * self.len_y2**2))

    def mean(self, v: int, J: int) -> float:
        """The B-level correction at (v', J'), MHz."""
        return float(self._k(float(v), J * (J + 1) / 1e4) @ self.beta)

    def __call__(self, v: int, J: int) -> tuple[float, float]:
        """(mean, 1 sigma) of the B-level correction at (v', J'), MHz."""
        y = J * (J + 1) / 1e4
        k = self._k(float(v), y)
        var = self.a1**2 + self.a2**2 - k @ self.M @ k        # the level's own prior, less what the data explain
        return float(k @ self.beta), float(np.sqrt(max(var, 0.0)))

    def covers(self, v: int) -> bool:
        return self.v_min <= v <= self.v_max


@lru_cache(maxsize=4)
def load_gp_correction(name: str) -> GPCorrection:
    return GPCorrection.load(name)
