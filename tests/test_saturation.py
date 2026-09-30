"""Structure of Doppler-free spectra: Lamb dips, crossovers and the modulated line shapes."""

import numpy as np
import pytest

from i2spec import RovibronicModel
from i2spec.saturation import (doppler_width, lineshape, line_resonances, resonances, shared_level, signal,
                               transit_time_width)


@pytest.fixture(scope="module")
def model():
    return RovibronicModel("127I2")


@pytest.fixture(scope="module")
def line(model):
    """R(56) 32-0, the 532 nm frequency standard: 15 main components at J'' = 56."""
    nu0, components = model.hyperfine_components(32, 0, 56, "R")
    return nu0, components, doppler_width("127I2", nu0, 300.0)


def test_doppler_width(line):
    _, _, width = line
    assert width == pytest.approx(438.0, rel=0.02)  # MHz at 532 nm, 300 K


def test_one_dip_per_component(line):
    _, components, width = line
    dips = [r for r in resonances(components, width, threshold=0.0) if r.kind == "dip"]
    assert len(dips) == len(components)
    for r in dips:
        c = components[r.components[0]]
        assert r.offset == c.offset and r.labels == (c.label, c.label) and r.sharing is None


@pytest.mark.parametrize("exponent", [1.0, 2.0])
def test_dip_amplitudes_follow_the_strength_law(line, exponent):
    _, components, width = line
    dips = {r.components[0]: r.amplitude for r in resonances(components, width, threshold=0.0, exponent=exponent)
            if r.kind == "dip"}
    strong = sorted(dips, key=lambda i: -components[i].strength)[:6]
    ratio = [dips[i] / components[i].strength**exponent for i in strong]
    assert np.allclose(ratio, ratio[0], rtol=1e-12)


def test_crossovers_sit_at_midpoints_of_level_sharing_pairs(line):
    _, components, width = line
    crossovers = [r for r in resonances(components, width, threshold=0.0) if r.kind == "crossover"]
    assert crossovers
    for r in crossovers:
        i, j = r.components
        assert i != j
        assert r.offset == pytest.approx((components[i].offset + components[j].offset) / 2, abs=1e-9)
        assert r.sharing == shared_level(components[i], components[j])
        assert r.sharing in ("lower", "upper")


def test_main_components_never_produce_crossovers(line):
    """Main (ΔF = ΔJ) components take one upper and one lower level each, so no pair shares one."""
    _, components, width = line
    labelled = [r for r in resonances(components, width, threshold=0.0)
                if r.kind == "crossover" and r.labels[0] and r.labels[1]]
    assert labelled == []


def test_crossovers_need_a_populated_velocity_class(line):
    _, components, width = line
    wide = [r for r in resonances(components, width, threshold=0.0) if r.kind == "crossover"]
    narrow = [r for r in resonances(components, width / 50, threshold=0.0) if r.kind == "crossover"]
    assert len(narrow) < len(wide)
    for r, w in ((wide, width), (narrow, width / 50)):
        for x in r:
            i, j = x.components
            assert abs(components[i].offset - components[j].offset) <= 3.0 * w


def test_lambda_weight_scales_only_the_upper_level_crossovers(line):
    _, components, width = line
    both = resonances(components, width, threshold=0.0)
    without = resonances(components, width, threshold=0.0, lambda_weight=0.0)
    assert {r.sharing for r in both} == {None, "lower", "upper"}
    assert {r.sharing for r in without} == {None, "lower"}


def test_amplitudes_are_normalized(line):
    _, components, width = line
    assert sum(r.amplitude for r in resonances(components, width)) == pytest.approx(1.0)


def test_threshold_keeps_the_strongest(line):
    _, components, width = line
    all_res = resonances(components, width, threshold=0.0)
    kept = resonances(components, width, threshold=1e-3)
    assert 0 < len(kept) < len(all_res)
    assert kept[0].offset == all_res[0].offset


def test_lineshape_is_a_lorentzian_without_modulation():
    x = np.array([0.0, 0.5, -0.5, 5.0])
    assert lineshape(x, 1.0) == pytest.approx([1.0, 0.5, 0.5, 1 / 101])


def test_third_harmonic_is_odd_and_crosses_zero_at_the_centre():
    x = np.linspace(-6, 6, 241)
    s = lineshape(x, 1.0, harmonic=3, modulation=1.0)
    assert abs(s[x == 0][0]) < 1e-12 * np.abs(s).max()
    assert s == pytest.approx(-s[::-1], abs=1e-12)


def test_third_harmonic_tends_to_the_third_derivative():
    x = np.linspace(-4, 4, 801)
    h = 1e-3
    lorentz = lambda u: 1 / (1 + (2 * u) ** 2)  # noqa: E731
    third = (lorentz(x + 2 * h) - 2 * lorentz(x + h) + 2 * lorentz(x - h) - lorentz(x - 2 * h)) / (2 * h**3)
    s = lineshape(x, 1.0, harmonic=3, modulation=0.02)
    assert np.corrcoef(s, third)[0, 1] > 0.9999
    assert s @ third > 0  # same sign convention as the derivative


def test_signal_adds_the_resonances(line):
    _, components, width = line
    res = resonances(components, width)
    nu = np.linspace(-800, 500, 1301)
    direct = sum(r.amplitude / (1 + (2 * (nu - r.offset) / 3.0) ** 2) for r in res)
    assert signal(nu, res, 3.0) == pytest.approx(direct)


def test_modulated_signal_matches_a_single_shape(line):
    _, components, width = line
    res = resonances(components, width)[:1]
    nu = np.linspace(-800, 500, 1301)
    expected = res[0].amplitude * lineshape(nu - res[0].offset, 3.0, harmonic=3, modulation=2.0)
    tolerance = 1e-5 * res[0].amplitude  # the kernel is interpolated; see signal()
    assert signal(nu, res, 3.0, harmonic=3, modulation=2.0) == pytest.approx(expected, abs=tolerance)


def test_dJ0_and_dJ2_give_the_same_main_dips(model):
    out = {}
    for dJ in (0, 2):
        nu0, res = line_resonances(model, 32, 0, 56, "R", 300.0, dJ=dJ, threshold=0.0)
        out[dJ] = sorted(r.offset for r in res if r.kind == "dip" and r.labels[0])
    assert len(out[0]) == len(out[2]) == 15
    assert np.allclose(out[0], out[2], atol=2.0)  # dJ = 0 is good to about 1 MHz


def test_transit_time_width_is_tens_of_kilohertz():
    assert transit_time_width(1.0, 300.0) == pytest.approx(0.12, rel=0.1)  # MHz for a 1 mm beam
