"""Wigner 3j, 6j and 9j symbols in floating point (Racah's formulae with log-factorials).

For the ranks used in i2spec (k <= 2) the sums have few terms, so the results agree with exact
SymPy values to ~1e-10 relative even for J ~ 300 (tests/test_wigner.py).
"""

from functools import lru_cache
from math import exp, lgamma


def _lf(n):
    return lgamma(n + 1.0)


def _is_int(x):
    return abs(x - round(x)) < 1e-9


def _triad(a, b, c):
    return abs(a - b) <= c <= a + b and _is_int(a + b + c)


def _log_delta(a, b, c):
    return 0.5 * (_lf(a + b - c) + _lf(a - b + c) + _lf(b + c - a) - _lf(a + b + c + 1))


def _phase(x):
    return -1.0 if round(x) % 2 else 1.0


@lru_cache(maxsize=None)
def wigner_3j(j1, j2, j3, m1, m2, m3):
    if abs(m1 + m2 + m3) > 1e-9 or not _triad(j1, j2, j3):
        return 0.0
    if abs(m1) > j1 or abs(m2) > j2 or abs(m3) > j3 or not all(_is_int(j - m) for j, m in ((j1, m1), (j2, m2), (j3, m3))):
        return 0.0
    pre = _log_delta(j1, j2, j3) + 0.5 * (_lf(j1 + m1) + _lf(j1 - m1) + _lf(j2 + m2) + _lf(j2 - m2) + _lf(j3 + m3) + _lf(j3 - m3))
    total = 0.0
    for t in range(round(max(0, j2 - j3 - m1, j1 - j3 + m2)), round(min(j1 + j2 - j3, j1 - m1, j2 + m2)) + 1):
        total += _phase(t) * exp(pre - _lf(t) - _lf(j3 - j2 + t + m1) - _lf(j3 - j1 + t - m2)
                                 - _lf(j1 + j2 - j3 - t) - _lf(j1 - t - m1) - _lf(j2 - t + m2))
    return _phase(j1 - j2 - m3) * total


@lru_cache(maxsize=None)
def wigner_6j(j1, j2, j3, j4, j5, j6):
    if not (_triad(j1, j2, j3) and _triad(j1, j5, j6) and _triad(j4, j2, j6) and _triad(j4, j5, j3)):
        return 0.0
    pre = _log_delta(j1, j2, j3) + _log_delta(j1, j5, j6) + _log_delta(j4, j2, j6) + _log_delta(j4, j5, j3)
    lower = (j1 + j2 + j3, j1 + j5 + j6, j4 + j2 + j6, j4 + j5 + j3)
    upper = (j1 + j2 + j4 + j5, j2 + j3 + j5 + j6, j3 + j1 + j6 + j4)
    total = 0.0
    for t in range(round(max(lower)), round(min(upper)) + 1):
        total += _phase(t) * exp(pre + _lf(t + 1) - sum(_lf(t - x) for x in lower) - sum(_lf(x - t) for x in upper))
    return total


@lru_cache(maxsize=None)
def wigner_9j(j1, j2, j3, j4, j5, j6, j7, j8, j9):
    """Sum over three 6j symbols (Edmonds 6.4.3)."""
    x = max(abs(j1 - j9), abs(j4 - j8), abs(j2 - j6))
    total = 0.0
    while x <= min(j1 + j9, j4 + j8, j2 + j6) + 1e-9:
        total += (_phase(2 * x) * (2 * x + 1) * wigner_6j(j1, j4, j7, j8, j9, x) * wigner_6j(j2, j5, j8, j4, x, j6)
                  * wigner_6j(j3, j6, j9, x, j1, j2))
        x += 1
    return total
