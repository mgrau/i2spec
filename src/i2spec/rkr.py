"""RKR: classical turning points from vibrational and rotational constants.

The first-order Rydberg-Klein-Rees inversion turns G(v) and B(v) into the potential that produced
them, as a pair of turning points per level:

    f(v) = C int_{v_min}^{v} [G(v) - G(w)]^(-1/2) dw          = (R_out - R_in) / 2
    g(v) = (1/C) int_{v_min}^{v} B(w) [G(v) - G(w)]^(-1/2) dw = (1/R_in - 1/R_out) / 2

with C = (hbar / sqrt(2 mu h c)) = 4.1058045 / sqrt(mu[amu]) Angstrom cm^(1/2), so that

    R_out, R_in = sqrt(f/g + f^2) +/- f.

v_min = -1/2 is the bottom of the well, where both integrands are singular as (v - w)^(-1/2); the
substitution w = v - u^2 removes the singularity exactly, which is why the quadrature below is in u.
This is the classical starting point for a direct-potential fit: Gerstenkorn & Luc's constants
(i2spec.dunham) give a curve that owes nothing to the Hannover parameters.
"""

from __future__ import annotations

import numpy as np
from numpy.polynomial.legendre import leggauss

from . import dunham
from .constants import reduced_mass

#: hbar / sqrt(2 mu h c) in Angstrom cm^(1/2), with mu in atomic mass units.
RKR_CONSTANT = 4.1058045

V_MIN = -0.5


def _integrals(state, v, isotopologue="127I2", n=160, dunham_model=dunham):
    """(f, g) of the RKR equations, in Angstrom and 1/Angstrom."""
    C = RKR_CONSTANT / np.sqrt(reduced_mass(isotopologue))
    Gv = dunham_model.G(state, v)
    # w = v - u^2, u from 0 to sqrt(v - V_MIN): dw = -2u du and the singularity cancels
    u_max = np.sqrt(v - V_MIN)
    x, wq = leggauss(n)
    u = 0.5 * u_max * (x + 1.0)
    w = v - u**2
    dG = Gv - np.array([dunham_model.G(state, float(t)) for t in w])
    integrand = 2.0 * u / np.sqrt(np.maximum(dG, 1e-12))
    B = np.array([dunham_model.B(state, float(t)) for t in w])
    scale = 0.5 * u_max * wq
    return C * float(np.sum(scale * integrand)), float(np.sum(scale * integrand * B)) / C


def turning_points(state, v, isotopologue="127I2", **kw):
    """(R_in, R_out) in Angstrom of the classical turning points of level v."""
    f, g = _integrals(state, v, isotopologue, **kw)
    root = np.sqrt(f / g + f * f)
    return root - f, root + f


def curve(state, v_max=None, n_v=60, isotopologue="127I2", **kw):
    """Turning points on a grid of v: (R, V) with V the term value above the potential minimum, cm-1.

    Both branches of every level, sorted in R, with the minimum (R_e, 0) included. R_e is where the
    two branches meet as v -> V_MIN, taken as the limit of (R_in + R_out)/2 at the lowest v used.
    """
    if v_max is None:
        v_max = dunham.constants()[state]["valid_v_max"]
    vs = np.concatenate([np.linspace(V_MIN + 0.01, 0.5, 8), np.linspace(1.0, v_max, n_v)])
    R, V = [], []
    for v in vs:
        r_in, r_out = turning_points(state, float(v), isotopologue, **kw)
        e = dunham.G(state, float(v)) - dunham.G(state, V_MIN)
        R += [r_in, r_out]
        V += [e, e]
    order = np.argsort(R)
    return np.array(R)[order], np.array(V)[order]


def re(state, isotopologue="127I2", **kw):
    """Equilibrium bond length, Angstrom: the midpoint of the turning points as v -> the well bottom."""
    r_in, r_out = turning_points(state, V_MIN + 1e-4, isotopologue, **kw)
    return 0.5 * (r_in + r_out)
