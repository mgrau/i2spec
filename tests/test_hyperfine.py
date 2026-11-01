"""The coupled-basis hyperfine Hamiltonian against brute-force constructions (tests/hyperfine_helpers.py), which
share no angular-momentum algebra with i2spec.hyperfine.
"""

import numpy as np
import pytest

from i2spec.hyperfine import K_Q, K_T, build_blocks, level_structure, line_components
from hyperfine_helpers import SPIN, P, dot, dot2, product_space, rank2


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
