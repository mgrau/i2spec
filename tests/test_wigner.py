"""Floating-point Wigner symbols against exact SymPy values."""

import random

import numpy as np
from sympy import Rational
from sympy.physics.wigner import wigner_3j as exact_3j
from sympy.physics.wigner import wigner_6j as exact_6j
from sympy.physics.wigner import wigner_9j as exact_9j

from i2spec.wigner import wigner_3j, wigner_6j, wigner_9j


def _r(x):
    return Rational(round(2 * x), 2)


def _check(fast, exact, cases):
    got = np.array([fast(*c) for c in cases])
    want = np.array([float(exact(*map(_r, c))) for c in cases])
    np.testing.assert_allclose(got, want, rtol=1e-9, atol=1e-14)


def test_3j_rotational_symbols_up_to_J_300():
    cases = [(J + d, k, J, 0, 0, 0) for J in (0, 1, 2, 17, 56, 57, 150, 299) for k in (1, 2) for d in range(-k, k + 1)
             if J + d >= 0]
    _check(wigner_3j, exact_3j, cases)


def test_3j_random_small():
    rng = random.Random(3)
    cases = []
    while len(cases) < 200:
        j1, j2 = rng.randint(0, 12) / 2, rng.randint(0, 12) / 2
        j3 = abs(j1 - j2) + rng.randint(0, round(j1 + j2 - abs(j1 - j2)))
        m1, m2 = rng.randint(-round(2 * j1), round(2 * j1)) / 2, rng.randint(-round(2 * j2), round(2 * j2)) / 2
        cases.append((j1, j2, j3, m1, m2, -m1 - m2))
    _check(wigner_3j, exact_3j, cases)


def test_6j_hyperfine_symbols_up_to_J_300():
    rng = random.Random(1)
    cases = []
    while len(cases) < 300:
        J, k, I = rng.randint(0, 300), rng.choice([1, 2]), rng.randint(0, 7)
        case = (J + rng.randint(-I, I), I + rng.choice([-2, 0, 2]), J + rng.choice([-2, -1, 0, 1, 2]), k, J, I)
        if min(case) >= 0:
            cases.append(case)
    cases += [(2.5, Ip, 2.5, I, 2.5, 2) for I in range(6) for Ip in range(6)]
    cases += [(1, 1, 2, 2.5, 2.5, 2.5), (1, 1, 2, 3.5, 3.5, 3.5)]
    _check(wigner_6j, exact_6j, cases)


def test_9j_spin_symbols():
    cases = [(i, i, 1, i, i, 1, Ip, I, 2) for i in (2.5, 3.5) for I in range(round(2 * i) + 1) for Ip in range(round(2 * i) + 1)]
    cases += [(2.5, 2.5, 1, 3.5, 3.5, 1, Ip, I, 2) for I in range(1, 7) for Ip in range(1, 7)]
    _check(wigner_9j, exact_9j, cases)
