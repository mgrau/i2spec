"""The coupled-basis hyperfine Hamiltonian against brute-force constructions.

The brute force works in the uncoupled basis |J M>|m1>|m2>. It uses explicit spin matrices and
C²_q(n) from quadrature over spherical harmonics, so it shares no angular-momentum algebra
(6j/9j symbols, reduced matrix elements, phases) with i2spec.hyperfine.
"""

import numpy as np
import pytest
from scipy.special import roots_legendre, sph_harm_y
from sympy.physics.quantum.cg import CG

from i2spec.hyperfine import K_Q, K_T, HyperfineParameters, build_blocks, level_structure, line_components

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


def product_space(J_values):
    labels, J_ops, C2 = rotation_ops(J_values)
    nr, ns = len(labels), round(2 * SPIN + 1)
    one_r, one_s = np.eye(nr), np.eye(ns)
    S = spin_ops(SPIN)
    J = [np.kron(op, np.kron(one_s, one_s)) for op in J_ops]
    I1 = [np.kron(one_r, np.kron(op, one_s)) for op in S]
    I2 = [np.kron(one_r, np.kron(one_s, op)) for op in S]
    C = {q: np.kron(op, np.kron(one_s, one_s)) for q, op in C2.items()}
    return labels, J, I1, I2, C


def dot(a, b):
    return sum(x @ y for x, y in zip(a, b))


def coupled_spectrum(J_values, rot):
    blocks = build_blocks(J_values, lambda J: P, rot, symmetry=None)
    return np.sort(np.concatenate([np.repeat(np.linalg.eigvalsh(H), 2 * F + 1) for F, (_, H) in blocks.items()]))


@pytest.mark.parametrize("J_values", [[1, 3], [2, 4], [0, 2, 4]])
def test_coupled_basis_matches_brute_force_with_dJ_2_mixing(J_values):
    rot = lambda J: 40.0 * J * (J + 1)  # small rotational constant, so ΔJ = ±2 mixing is strong
    labels, J, I1, I2, C = product_space(J_values)
    quad = SPIN * (2 * SPIN - 1)
    H = np.diag(np.kron([rot(Jl) for Jl, _ in labels], np.ones(36))).astype(complex)
    H += 1e-3 * P.C * (dot(J, I1) + dot(J, I2)) + 1e-3 * P.delta * dot(I1, I2)
    H += K_Q * P.eqQ * np.sqrt(6) / (4 * quad) * (dot2(rank2(I1, I1), C) + dot2(rank2(I2, I2), C))
    H += K_T * 1e-3 * P.d * dot2(rank2(I1, I2), C)
    assert np.allclose(H, H.conj().T)
    np.testing.assert_allclose(np.linalg.eigvalsh(H), coupled_spectrum(J_values, rot), atol=1e-8)


@pytest.mark.parametrize("Jv", [2, 3, 6])
def test_dJ0_blocks_match_literature_operator_forms(Jv):
    _, J, I1, I2, _ = product_space([Jv])
    x, i, n = Jv * (Jv + 1), SPIN, 36 * (2 * Jv + 1)
    one = np.eye(n)
    IJ1, IJ2 = dot(I1, J), dot(I2, J)
    eq = sum(3 * a @ a + 1.5 * a - i * (i + 1) * x * one for a in (IJ1, IJ2)) / (2 * i * (2 * i - 1) * (2 * Jv - 1) * (2 * Jv + 3))
    tss = (3 * IJ1 @ IJ2 + 3 * IJ2 @ IJ1 - 2 * x * dot(I1, I2)) / ((2 * Jv - 1) * (2 * Jv + 3))
    H = -P.eqQ * eq + 1e-3 * (P.C * (IJ1 + IJ2) + P.d * tss + P.delta * dot(I1, I2))
    np.testing.assert_allclose(np.linalg.eigvalsh(H), coupled_spectrum([Jv], lambda J: 0.0), atol=1e-8)


@pytest.mark.parametrize("J, symmetry, n_levels", [(56, "g", 15), (57, "g", 21), (57, "u", 15), (56, "u", 21)])
def test_level_counts_follow_nuclear_spin_statistics(J, symmetry, n_levels):
    levels = level_structure(J, lambda J: P, lambda J: 1100.0 * J * (J + 1), symmetry=symmetry)
    assert len(levels) == n_levels


def test_line_strengths_and_main_component_labels():
    lower = level_structure(56, lambda J: P, lambda J: 1120.0 * J * (J + 1), symmetry="g")
    upper = level_structure(57, lambda J: P, lambda J: 700.0 * J * (J + 1), symmetry="u")
    comps = line_components(upper, lower, 57, 56)
    assert sum(c.strength for c in comps) == pytest.approx(1.0)
    main = [c for c in comps if c.label]
    assert [c.label for c in sorted(main, key=lambda c: c.offset)] == [f"a{n}" for n in range(1, 16)]
    assert sum(c.strength for c in main) > 0.99  # ΔF ≠ ΔJ components are O(1/J²)
