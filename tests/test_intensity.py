"""Intensities against Tellinghuisen, J. Chem. Phys. 134, 084301 (2011)."""

import numpy as np
import pytest

from i2spec.intensity import (S_UNIT, intensity_model, line_list, mu_tellinghuisen2011, nuclear_spin_weight,
                              partition_function, shared_offset)
from i2spec.spectrum import cross_section, doppler_fwhm, vapor_pressure


@pytest.fixture(scope="module")
def model():
    return intensity_model("127I2")


def test_strength_constant_matches_tellinghuisen():
    assert S_UNIT == pytest.approx(108.862 * 3.8235e-21, rel=2e-4)  # eq. (3) times ln10·1000/N_A


@pytest.mark.parametrize("iso, even, odd", [("127I2", 15, 21), ("129I2", 28, 36), ("127I129I", 48, 48)])
def test_nuclear_spin_weights(iso, even, odd):
    assert (nuclear_spin_weight(56, iso), nuclear_spin_weight(57, iso)) == (even, odd)


def test_transition_moment_matches_linear_fit_over_its_range():
    R = np.linspace(2.63, 2.88, 11)
    linear = (1.1123 + 0.712 * (R - 2.85)) ** 2  # eq. (9)
    np.testing.assert_allclose(mu_tellinghuisen2011(R) ** 2, linear, rtol=0.02)


@pytest.mark.parametrize("v_upper", [20, 26, 32, 40])
def test_r_centroids_match_tellinghuisen_fit(model, v_upper):
    """R-centroids of v''=0 bands vs Tellinghuisen's fit to λ: tests wavefunctions and overlaps."""
    X, B = model.states["X"], model.states["B"]
    e_x, c_x = X.states(10)
    e_b, c_b = B.states(11)
    c_x = c_x[shared_offset(X, B):]                  # X extends B's grid inward; B vanishes there
    overlap = c_b[:, v_upper] @ c_x[:, 0]
    r_centroid = (c_b[:, v_upper] * B.R) @ c_x[:, 0] / overlap
    lam = 1e7 / (e_b[v_upper] - e_x[0])
    assert r_centroid == pytest.approx(2.6322 + 0.001548 * (lam - 500) - 9.1e-7 * (lam - 500) ** 2, abs=0.005)


def test_partition_function_close_to_rigid_rotor_harmonic_estimate(model):
    T = 295.0
    kT = 0.6950348 * T
    estimate = 18 * kT / 0.03737 / (1 - np.exp(-214.5 / kT))  # mean nuclear-spin weight 18
    assert partition_function(model, T) == pytest.approx(estimate, rel=0.05)


def test_doppler_width_and_vapor_pressure():
    assert doppler_fwhm(18788.0, 300.0) == pytest.approx(4.494e-8 * 18788.0 * np.sqrt(300.0), rel=1e-3)
    assert vapor_pressure(293.15) / 133.322368 == pytest.approx(0.201, rel=0.01)  # Torr at 20 °C


def test_line_list_and_cross_section_near_532_nm(model):
    lines = line_list(model, 295.0, 18787.0, 18790.0, S_min=1e-24)
    strongest = int(np.argmax(lines.S))
    assert lines.label(strongest) in {"R(56) 32-0", "P(53) 32-0", "R(57) 32-0", "P(54) 32-0"} or lines.v_lower[strongest] == 0
    grid = np.arange(18787.0, 18790.0, 0.0005)
    sigma = cross_section(grid, lines, 295.0)
    # integrating the cross section recovers the summed line strengths of lines well inside the window
    inside = (lines.nu > 18787.2) & (lines.nu < 18789.8)
    window = (grid > 18787.1) & (grid < 18789.9)
    assert np.trapezoid(sigma[window], grid[window]) == pytest.approx(lines.S[inside].sum(), rel=0.05)
