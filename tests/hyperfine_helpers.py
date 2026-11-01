"""The brute-force hyperfine Hamiltonian of tests/test_hyperfine.py and tests/test_isotopologue_hyperfine.py.

The brute force works in the uncoupled basis |J M>|m1>|m2>. It uses explicit spin matrices and
C²_q(n) from quadrature over spherical harmonics, so it shares no angular-momentum algebra
(6j/9j symbols, reduced matrix elements, phases) with i2spec.hyperfine.
"""

import numpy as np
from scipy.special import roots_legendre, sph_harm_y
from sympy.physics.quantum.cg import CG

from i2spec.hyperfine import HyperfineParameters

SPIN = 2.5
P = HyperfineParameters(eqQ=-600.0, C=5000.0, d=-3000.0, delta=2000.0)  # exaggerated magnetic terms


def spin_ops(s):
    m = np.arange(s, -s - 1, -1)
    raise_ = np.zeros((m.size, m.size))
    for k in range(1, m.size):
        raise_[k - 1, k] = np.sqrt(s * (s + 1) - m[k] * (m[k] + 1))
    return [(raise_ + raise_.T) / 2, (raise_ - raise_.T) / 2j, np.diag(m).astype(complex)]


def rotation_ops(J_values):
    labels = [(J, M) for J in J_values for M in range(J, -J - 1, -1)]
    n = len(labels)
    J_ops = [np.zeros((n, n), complex) for _ in range(3)]
    start = 0
    for J in J_values:
        for a, op in enumerate(spin_ops(J)):
            J_ops[a][start:start + 2 * J + 1, start:start + 2 * J + 1] = op
        start += 2 * J + 1
    x, w = roots_legendre(40)
    phi = np.linspace(0, 2 * np.pi, 48, endpoint=False)
    theta, phi = np.meshgrid(np.arccos(x), phi, indexing="ij")
    weights = np.outer(w, np.full(phi.shape[1], 2 * np.pi / phi.shape[1]))
    Y = np.array([sph_harm_y(J, M, theta, phi) for J, M in labels])
    C2 = {q: np.einsum("aij,ij,bij->ab", Y.conj(), np.sqrt(4 * np.pi / 5) * sph_harm_y(2, q, theta, phi) * weights, Y)
          for q in range(-2, 3)}
    return labels, J_ops, C2


def spherical(v):
    return {1: -(v[0] + 1j * v[1]) / np.sqrt(2), 0: v[2], -1: (v[0] - 1j * v[1]) / np.sqrt(2)}


def rank2(a, b):
    a, b = spherical(a), spherical(b)
    return {q: sum(float(CG(1, q1, 1, q - q1, 2, q).doit()) * (a[q1] @ b[q - q1]) for q1 in (-1, 0, 1) if abs(q - q1) <= 1)
            for q in range(-2, 3)}


def dot2(t, u):
    return sum((-1) ** q * t[q] @ u[-q] for q in range(-2, 3))


def product_space(J_values, i1=SPIN, i2=SPIN):
    """The rotation and the two nuclear spins (i1, i2) on |J M>|m1>|m2>: labels, J, I1, I2 and C²_q."""
    labels, J_ops, C2 = rotation_ops(J_values)
    one_r, one_1, one_2 = np.eye(len(labels)), np.eye(round(2 * i1 + 1)), np.eye(round(2 * i2 + 1))
    J = [np.kron(op, np.kron(one_1, one_2)) for op in J_ops]
    I1 = [np.kron(one_r, np.kron(op, one_2)) for op in spin_ops(i1)]
    I2 = [np.kron(one_r, np.kron(one_1, op)) for op in spin_ops(i2)]
    C = {q: np.kron(op, np.kron(one_1, one_2)) for q, op in C2.items()}
    return labels, J, I1, I2, C


def dot(a, b):
    return sum(x @ y for x, y in zip(a, b))
