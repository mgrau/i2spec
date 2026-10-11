"""Master line list (temperature scaling, caching, box-state filter) and apparent cross sections."""

import numpy as np
import pytest

from i2spec.constants import MHZ_PER_CM
from i2spec.intensity import C2, bound, intensity_model, line_list, master_line_list, match_levels, shared_offset
from i2spec.model import RovibronicModel
from i2spec.spectrum import apparent_cross_section, ils_kernel

WINDOW = (18780.0, 18800.0)


@pytest.fixture(scope="module")
def model():
    return intensity_model(nlev_x=12, nlev_b=45, grid=dict(step=0.008))  # coarser grid: fast, fine for intensities


@pytest.fixture(scope="module")
def master(model):
    return master_line_list(model, *WINDOW, T_range=(250.0, 450.0), S_min=1e-26, cache=False)


def test_temperature_scaling_is_exact(master):
    a, b = master.at(300.0), master.at(400.0)
    expected = (np.exp(-C2 * a.E_lower * (1 / 400.0 - 1 / 300.0))
                * master.partition_function(300.0) / master.partition_function(400.0))
    np.testing.assert_allclose(b.S / a.S, expected, rtol=1e-12)


def test_master_list_agrees_with_single_temperature_list(model, master):
    direct = line_list(model, 300.0, *WINDOW, S_min=1e-24)
    via = master.at(300.0, *WINDOW, S_min=1e-24)
    assert len(direct) == len(via) > 50
    np.testing.assert_allclose(direct.nu, via.nu)
    np.testing.assert_allclose(direct.S, via.S, rtol=1e-7)


def test_cache_round_trip(model, tmp_path, monkeypatch):
    monkeypatch.setenv("I2SPEC_CACHE", str(tmp_path))
    built = master_line_list(model, 18780.0, 18782.0, T_range=(290.0, 300.0), S_min=1e-25)
    assert len(list(tmp_path.glob("linelist_127I2_*.npz"))) == 1
    loaded = master_line_list(model, 18780.0, 18782.0, T_range=(290.0, 300.0), S_min=1e-25)
    np.testing.assert_array_equal(built.nu, loaded.nu)
    np.testing.assert_array_equal(built.strength0, loaded.strength0)
    assert loaded.isotopologue == "127I2"


def test_positions_come_from_the_b_spline_solver(master):
    """Line positions equal RovibronicModel's transitions, not the sinc-DVR's.

    The DVR was 0.24 MHz off at R(56) 32-0 and a few MHz elsewhere; the list used to carry that error
    into the lookup, the TUI and the web app. Here the DVR grid is deliberately coarse (0.008 Å), and the
    positions must still agree to the Hz.
    """
    m = RovibronicModel("127I2")
    d = np.array([(master.nu[i] - m.transition(int(master.v_upper[i]), int(master.v_lower[i]), int(master.J_lower[i]),
                                                "R" if master.branch[i] > 0 else "P")) * MHZ_PER_CM
                  for i in range(len(master))])
    assert len(d) > 100 and np.abs(d).max() < 1e-3                      # MHz
    i = np.flatnonzero((master.v_upper == 32) & (master.v_lower == 0) & (master.J_lower == 56) & (master.branch == 1))
    assert i.size == 1 and master.nu[i[0]] == pytest.approx(m.transition(32, 0, 56, "R"), abs=1e-7)


def test_levels_are_matched_by_energy_not_index():
    """A box state that one solver puts between two levels must not shift the rest."""
    dvr = np.array([1.0, 2.0005, 3.0, 5.0])
    exact = np.array([1.0001, 1.5, 2.0, 3.0002])                          # 1.5 is a box state the DVR lacks
    matched, ok = match_levels(dvr, exact)
    np.testing.assert_allclose(matched, [1.0001, 2.0, 3.0002, 5.0])
    assert list(ok) == [True, True, True, False]                          # 5.0 has no partner: kept, flagged


def test_x_grid_reaches_the_inner_wall_of_high_v():
    """At 2.33 Å the grid edge was a hard wall for X: v'' = 30 came out 1.6 GHz high."""
    model = intensity_model("127I2")
    X, B = model.states["X"], model.states["B"]
    assert X.R[0] == pytest.approx(2.10, abs=0.005) and B.R[0] == pytest.approx(2.33)
    assert np.allclose(X.R[shared_offset(X, B):], B.R)
    exact = RovibronicModel("127I2").states["X"].levels(0)
    assert abs(X.levels(0)[30] - exact[30]) * MHZ_PER_CM < 5.0
    with pytest.raises(ValueError):
        shared_offset(B, X)


def test_box_states_are_flagged():
    R = np.linspace(2.3, 7.0, 500)
    localized = np.exp(-(((R - 3.0) / 0.1) ** 2))
    spread = np.ones_like(R)
    vectors = np.column_stack([localized / np.linalg.norm(localized), spread / np.linalg.norm(spread)])
    assert list(bound(vectors, R)) == [True, False]


@pytest.mark.parametrize("kind, width", [("gauss", 0.5), ("boxcar", 0.5), ("sinc", 0.2), ("sinc2", 0.2)])
def test_ils_kernels_have_unit_sum(kind, width):
    assert ils_kernel(0.002, kind, width).sum() == pytest.approx(1.0)


def test_apparent_cross_section_thin_limit_and_saturation(master):
    lines = master.at(300.0, 18781.0, 18799.0, S_min=1e-24)
    nu = np.linspace(18785.0, 18795.0, 21)
    thin_a = apparent_cross_section(nu, lines, 300.0, 1e8, ils=("gauss", 0.5))
    thin_b = apparent_cross_section(nu, lines, 300.0, 1e10, ils=("gauss", 0.5))
    np.testing.assert_allclose(thin_a, thin_b, rtol=1e-4)  # optically thin: independent of N
    thick = apparent_cross_section(nu, lines, 300.0, 3e17, ils=("gauss", 0.5))
    assert np.all(thick <= thin_b * (1 + 1e-9))
    assert thick.mean() < 0.9 * thin_b.mean()  # unresolved strong lines saturate


def test_apparent_cross_section_stays_finite_when_saturated(master):
    # a sinc's negative lobes, and round-off in a black core, take the convolved transmission to zero or below
    lines = master.at(300.0, 18781.0, 18799.0, S_min=1e-24)
    nu = np.linspace(18785.0, 18795.0, 2001)
    with np.errstate(invalid="raise", divide="raise"):
        sigma = apparent_cross_section(nu, lines, 300.0, 1e18, ils=("sinc", 0.01))
    assert np.all(np.isfinite(sigma))   # (a sinc also overshoots: negative values are its ringing, not an error)
