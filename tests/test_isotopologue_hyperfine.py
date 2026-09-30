"""¹²⁹I₂ and ¹²⁷I¹²⁹I hyperfine structure (Salumbides et al., Mol. Phys. 104, 2641 (2006), S06).

Targets: peaks of the calculated profiles in S06 Figs. 4 and 5, digitized from the PDF vector paths
to about 0.1 MHz (docs/research/isotopologue-hyperfine.md §4.1b). The profiles come from S06's
per-line fits, not from its Table 2-4 formulae, so differences of a few 0.1 MHz are expected.
"""

import numpy as np
import pytest
from scipy.signal import find_peaks

from i2spec import RovibronicModel, hfs_params
from i2spec.hyperfine import K_Q, K_T, HyperfineParameters, build_blocks, level_structure
from test_hyperfine import dot, dot2, rank2, rotation_ops, spin_ops

P = HyperfineParameters(eqQ=-600.0, C=5000.0, d=-3000.0, delta=2000.0)  # exaggerated magnetic terms
EQQ_RATIO, C_RATIO = 0.701213, 0.6655

# (Lorentzian FWHM, peak intervals from the first peak), MHz.
FIG4_129I2_R34_11_3 = (11.9, [0.00, 136.75, 156.20, 217.41, 238.39, 247.08, 269.61, 292.54, 388.20, 407.88, 448.08,
                              460.07, 484.24, 524.82, 543.62, 582.81])
FIG5_127I129I_R37_11_3 = (12.2, [0.00, 21.99, 44.59, 136.44, 170.54, 196.05, 236.42, 265.81, 291.95, 304.25, 323.94,
                                 344.28, 415.16, 459.13, 504.70, 518.48, 551.13, 562.74, 586.89, 618.42, 630.52, 682.00,
                                 695.49, 706.09, 726.34, 745.87])


@pytest.fixture(scope="module")
def models():
    return {iso: RovibronicModel(iso) for iso in ("127I2", "129I2", "127I129I")}


def product_space(J_values, i1, i2):
    labels, J_ops, C2 = rotation_ops(J_values)
    one_r, one_1, one_2 = np.eye(len(labels)), np.eye(round(2 * i1 + 1)), np.eye(round(2 * i2 + 1))
    J = [np.kron(op, np.kron(one_1, one_2)) for op in J_ops]
    I1 = [np.kron(one_r, np.kron(op, one_2)) for op in spin_ops(i1)]
    I2 = [np.kron(one_r, np.kron(one_1, op)) for op in spin_ops(i2)]
    C = {q: np.kron(op, np.kron(one_1, one_2)) for q, op in C2.items()}
    return labels, J, I1, I2, C


@pytest.mark.parametrize("i1, i2", [(2.5, 3.5), (3.5, 2.5)])
@pytest.mark.parametrize("J_values", [[3], [1, 3], [2, 4]])
def test_heteronuclear_hamiltonian_matches_brute_force(J_values, i1, i2):
    """Per-nucleus eqQ and C, and ΔI = ±1 couplings, against the uncoupled product space."""
    rot = lambda J: 40.0 * J * (J + 1)
    labels, J, I1, I2, C = product_space(J_values, i1, i2)
    n = round(2 * i1 + 1) * round(2 * i2 + 1)
    H = np.diag(np.kron([rot(Jl) for Jl, _ in labels], np.ones(n))).astype(complex)
    H += 1e-3 * P.C * (dot(J, I1) + C_RATIO * dot(J, I2)) + 1e-3 * P.delta * dot(I1, I2)
    H += K_Q * P.eqQ * np.sqrt(6) / 4 * (dot2(rank2(I1, I1), C) / (i1 * (2 * i1 - 1))
                                          + EQQ_RATIO * dot2(rank2(I2, I2), C) / (i2 * (2 * i2 - 1)))
    H += K_T * 1e-3 * P.d * dot2(rank2(I1, I2), C)
    blocks = build_blocks(J_values, lambda J: P, rot, i1, i2, None, EQQ_RATIO, C_RATIO)
    coupled = np.sort(np.concatenate([np.repeat(np.linalg.eigvalsh(h), 2 * F + 1) for F, (_, h) in blocks.items()]))
    np.testing.assert_allclose(np.linalg.eigvalsh(H), coupled, atol=1e-8)


@pytest.mark.parametrize("J", [5, 6])
def test_distinguishable_nuclei_code_reproduces_homonuclear_levels(J):
    """S06 §3.1: with identical nuclei and constants, symmetry=None gives the union of the g and u levels."""
    rot = lambda J: 1100.0 * J * (J + 1)
    homonuclear = [lv.energy for s in ("g", "u") for lv in level_structure(J, lambda J: P, rot, symmetry=s)]
    distinguishable = [lv.energy for lv in level_structure(J, lambda J: P, rot, symmetry=None)]
    np.testing.assert_allclose(sorted(distinguishable), sorted(homonuclear), atol=1e-9)


@pytest.mark.parametrize("iso, J_lower, n_main", [("127I2", 34, 15), ("127I2", 35, 21), ("129I2", 34, 28),
                                                  ("129I2", 35, 36), ("127I129I", 36, 48), ("127I129I", 37, 48)])
def test_main_component_counts(models, iso, J_lower, n_main):
    _, comps = models[iso].hyperfine_components(11, 3, J_lower, "R", dJ=0)
    assert sum(c.label is not None for c in comps) == n_main
    assert sum(c.strength for c in comps) == pytest.approx(1.0)


@pytest.mark.parametrize("v, J, eqQ, C, delta, d", [(32, 57, -544.814, 89.935, -6.75, -41.30),
                                                    (53, 88, -569.446, 490.376, 100.41, -271.57)])
def test_s06_b_state_values(models, v, J, eqQ, C, delta, d):
    """docs/research/isotopologue-hyperfine.md §8(b), Ω = 0 sign flipped, E at J = 0."""
    p = hfs_params.b_state_s06(v, J, models["127I2"]._reference_term("B", v))
    np.testing.assert_allclose([p.eqQ, p.C, p.delta, p.d], [eqQ, C, delta, d], atol=0.01)


def interval_residuals(components, fwhm, intervals):
    """Model minus tabulated peak intervals, for components convolved with a Lorentzian (strength^1)."""
    offset = np.array([c.offset for c in components])
    strength = np.array([c.strength for c in components])
    f = np.arange(offset.min() - 2 * fwhm, offset.max() + 2 * fwhm, 0.01)
    y = sum(s / (1 + (2 * (f - o) / fwhm) ** 2) for o, s in zip(offset, strength))
    peaks, _ = find_peaks(y, prominence=0.03 * np.ptp(y))
    p = f[peaks] - f[peaks[0]]
    return np.array([p[np.argmin(np.abs(p - t))] - t for t in intervals])


def test_s06_fig4_129I2_R34_11_3(models):
    _, comps = models["129I2"].hyperfine_components(11, 3, 34, "R")
    assert np.abs(interval_residuals(comps, *FIG4_129I2_R34_11_3)).max() < 0.5


def test_s06_fig5_127I129I_R37_11_3(models):
    _, comps = models["127I129I"].hyperfine_components(11, 3, 37, "R")
    assert np.abs(interval_residuals(comps, *FIG5_127I129I_R37_11_3)).max() < 0.5


def test_one_spin_rotation_constant_misses_fig5(models, monkeypatch):
    """Negative control: C(¹²⁹I) = C(¹²⁷I) instead of 0.6655 C(¹²⁷I)."""
    monkeypatch.setitem(hfs_params.GAMMA_MU, 129, 1.0)
    _, comps = models["127I129I"].hyperfine_components(11, 3, 37, "R")
    assert np.abs(interval_residuals(comps, *FIG5_127I129I_R37_11_3)).max() > 1.0
