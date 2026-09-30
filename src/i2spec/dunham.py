"""The Gerstenkorn & Luc 1985 Dunham description of the B-X system (J. Phys. 46, 867): the model of the
Orsay atlases, and the starting point of i2spec's own potential fit.

Term values are E(state, v, J) = T + G(v) + B_v x - D_v x^2 + H_v x^3 + L_v x^4 + M_v x^5, x = J(J+1),
with the X energies from the X potential minimum's v'' = 0, J = 0 level (E_X(0, 0) = 0) and the B
energies above it (T00 = 15 724.5871 cm-1). X: Dunham series for G, B, D, H (valid to v'' = 19).
B: Dunham series for G and B (to v' = 80) and exponential polynomials for D, H, L and the effective M*.
Data: i2spec/data/gerstenkorn1985_dunham.json (its citation and transcription record are inside).
"""

from __future__ import annotations

import json
from functools import lru_cache
from importlib import resources

import numpy as np


@lru_cache(maxsize=1)
def constants() -> dict:
    return json.loads(resources.files("i2spec").joinpath("data/gerstenkorn1985_dunham.json").read_text())


def _series(coefficients, v, first_power):
    u = v + 0.5
    return sum(c * u ** (i + first_power) for i, c in enumerate(coefficients))


def G(state: str, v: int) -> float:
    """Vibrational term G(v) = sum_{i>=1} Y_i0 (v+1/2)^i, cm-1, from the potential minimum."""
    return _series(constants()[state]["Y_i0"], v, 1)


def B(state: str, v: int) -> float:
    return _series(constants()[state]["Y_i1"], v, 0)


def cdc(state: str, v: int) -> tuple[float, float, float, float]:
    """(D_v, H_v, L_v, M_v) with the paper's signs: each term lowers the level as -D x^2 - |H| x^3 ..."""
    c = constants()[state]
    if state == "X":
        return _series(c["Y_i2"], v, 0), _series(c["Y_i3"], v, 0), 0.0, 0.0
    ex = lambda C: float(np.exp(_series(C, v, 0)))          # noqa: E731  exp sum_{i>=1} C_i (v+1/2)^(i-1)
    return ex(c["C_d"]), -ex(c["C_h"]), -ex(c["C_l"]), -ex(c["C_m"])


def energy(state: str, v: int, J: int) -> float:
    """Term value (cm-1) above X(v'' = 0, J = 0)."""
    x = J * (J + 1)
    D, H, L, M = cdc(state, v)
    e = G(state, v) - G(state, 0) + B(state, v) * x - D * x**2 + H * x**3 + L * x**4 + M * x**5
    return e + (constants()["T00"] if state == "B" else 0.0)


def transition(v_upper: int, v_lower: int, J_lower: int, branch: str) -> float:
    """Wavenumber (cm-1) of the P or R line, hyperfine-free."""
    J_upper = J_lower + (1 if branch == "R" else -1)
    return energy("B", v_upper, J_upper) - energy("X", v_lower, J_lower)
