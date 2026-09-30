"""Measured corrections to the term values the potentials give: δE(state, v, J) = Σ c_k y^k, y = J(J+1)/10⁴.

The published potentials miss the comb-referenced near-infrared lines (755-815 nm, bands (0-7)-(11-17))
by 2-6 MHz where the data are 0.02-0.3 MHz, with a shape that is smooth in J within a band. A parameter
set may name a corrections file; the model adds them to ¹²⁷I₂ levels, so every band sharing a corrected
level moves together. Fitted by prototypes/level_corrections_fit.py (docs/research/nir-model-2010.md §9).

A file may also carry a ``band`` section: a correction to the *lines* of the NIR bands v′ = 0 → v″ =
12-17, not to levels, applied only inside the J″ range measured at each v″ (prototypes/nir_band_fit.py).
It exists because what those lines still need after the level corrections - a J′-dependence of B v′ = 0
- cannot be put on the levels without moving every other band that shares one of them, where nothing
has been measured (docs/research/iodinespec5.md, "The NIR loss").
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from functools import lru_cache
from importlib import resources

import numpy as np

#: How far beyond the measured J range a correction is still trusted (uncertainty rule).
J_MARGIN = 15


@dataclass(frozen=True)
class LevelCorrections:
    name: str
    coefficients: dict          # state -> {v: (c0, c1, ...)} in MHz
    coverage: dict              # state -> {v: (J_min, J_max)} of the lines that fixed them
    fit: dict                   # the fit's own statistics, for the record
    held_out: dict              # state -> {v: rms (MHz) of the precise lines predicted with the level's own lines left out}
    band: dict | None = None    # the NIR band correction, if any (see the module docstring)
    clamp_J: bool = False       # beyond coverage +- J_MARGIN, extrapolate only the linear part (level_corrections_2026h on)
    covariance: dict | None = None   # state -> {v: covariance of the coefficients, MHz^2} (level_corrections_2026k on)
    discrepancy: dict | None = None  # state -> {v: model-discrepancy term, MHz}: held-out scatter the covariance leaves

    def band_covers(self, v_upper: int, v_lower: int, J_lower: int) -> bool:
        """Is this line one the band correction applies to: in its bands and within J_MARGIN of their data?"""
        b = self.band
        if b is None or v_upper != b["v_upper"] or not b["v_lower"][0] <= v_lower <= b["v_lower"][1]:
            return False
        r = b["coverage"].get(str(v_lower))
        return r is not None and r[0] - J_MARGIN <= J_lower <= r[1] + J_MARGIN

    def band_shift(self, v_upper: int, v_lower: int, J_lower: int, branch: str) -> float:
        """Correction (MHz) to the line (v'-v'') P/R(J''), zero outside the band correction's reach."""
        if not self.band_covers(v_upper, v_lower, J_lower):
            return 0.0
        b = self.band
        J_upper = J_lower + (1 if branch == "R" else -1)
        yu, yl, dv = J_upper * (J_upper + 1) / 1e4, J_lower * (J_lower + 1) / 1e4, v_lower - b["centre"]
        upper = sum(c * yu ** (k + 1) for k, c in enumerate(b["B"]))
        lower = sum(c * dv**l * yl**k for l, row in enumerate(b["X"]) for k, c in enumerate(row))
        return float(upper - lower)

    def shift(self, state: str, v: int, J: int) -> float:
        """Correction in MHz; zero where none was fitted."""
        c = self.coefficients.get(state, {}).get(v)
        if c is None:
            return 0.0
        y = J * (J + 1) / 1e4
        if self.clamp_J and len(c) > 2:
            r = self.coverage.get(state, {}).get(v)
            if r is not None:
                # Outside the measured J the constant and J(J+1) terms - a term-value and a rotational-constant
                # error, smooth by nature - keep extrapolating; the curvature terms, fitted to the data's own
                # scatter, are held at the edge. Checked against velchev1998a (B v' = 16 beyond its J range).
                Jc = min(max(J, max(r[0] - J_MARGIN, 0)), r[1] + J_MARGIN)
                yc = Jc * (Jc + 1) / 1e4
                return float(c[0] + c[1] * y + sum(ck * yc**k for k, ck in enumerate(c) if k >= 2))
        return float(np.polyval(c[::-1], y))

    def shift_levels(self, state: str, J: int, n: int) -> np.ndarray:
        """Corrections (MHz) for v = 0..n-1 at one J, to add to a solver's level array."""
        out = np.zeros(n)
        for v in self.coefficients.get(state, {}):
            if v < n:
                out[v] = self.shift(state, v, J)
        return out

    def uncertainty(self, state: str, v: int) -> float | None:
        """Held-out rms (MHz) of the correction, or None where a single line fixed the level."""
        return self.held_out.get(state, {}).get(v)

    def _basis(self, state, v, J, n):
        """The polynomial terms of shift() at J, with the same clamping of the curvature terms."""
        y = J * (J + 1) / 1e4
        yc = y
        r = self.coverage.get(state, {}).get(v)
        if self.clamp_J and n > 2 and r is not None:
            Jc = min(max(J, max(r[0] - J_MARGIN, 0)), r[1] + J_MARGIN)
            yc = Jc * (Jc + 1) / 1e4
        return np.array([y ** k if k < 2 else yc ** k for k in range(n)])

    def uncertainty_at(self, state: str, v: int, J: int) -> float | None:
        """1-sigma (MHz) of the corrected level at J: the coefficients' covariance propagated to J, with the
        level's model-discrepancy term in quadrature. None where the set carries no covariance for the level."""
        cov = (self.covariance or {}).get(state, {}).get(v)
        if cov is None:
            return None
        cov = np.asarray(cov, dtype=float)
        phi = self._basis(state, v, J, cov.shape[0])
        md = (self.discrepancy or {}).get(state, {}).get(v) or 0.0
        return float(np.sqrt(max(phi @ cov @ phi, 0.0) + md * md))

    def covers(self, state: str, v: int, J: int) -> bool:
        """Inside (or within J_MARGIN of) the J range the correction was fitted on."""
        r = self.coverage.get(state, {}).get(v)
        return r is not None and r[0] - J_MARGIN <= J <= r[1] + J_MARGIN


@lru_cache(maxsize=None)
def load_level_corrections(name: str) -> LevelCorrections:
    text = resources.files("i2spec").joinpath("data").joinpath(f"{name}.json").read_text()
    d = json.loads(text)
    return LevelCorrections(name=d["id"],
                            coefficients={st: {int(v): tuple(c) for v, c in d[st].items()} for st in ("X", "B")},
                            coverage={st: {int(v): tuple(r) for v, r in d.get("coverage", {}).get(st, {}).items()}
                                      for st in ("X", "B")},
                            fit=d.get("fit", {}),
                            held_out={st: {int(v): float(u) for v, u in d.get("held_out_MHz", {}).get(st, {}).items()}
                                      for st in ("X", "B")},
                            band=d.get("band"), clamp_J=bool(d.get("clamp_J", False)),
                            covariance={st: {int(v): c for v, c in d["covariance_MHz2"].get(st, {}).items()} for st in ("X", "B")}
                            if "covariance_MHz2" in d else None,
                            discrepancy={st: {int(v): float(u) for v, u in d["discrepancy_MHz"].get(st, {}).items()} for st in ("X", "B")}
                            if "discrepancy_MHz" in d else None)
