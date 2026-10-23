"""The fit's hyperfine beyond the core: the second-order engine's physics, and the release's tables."""
import json

import numpy as np
import pytest

from i2spec.constants import HBAR2_2U, MHZ_PER_CM, reduced_mass
from i2spec.fit import hyperfine_beyond as hb, hyperfine_rows, second_order as so, x_hyperfine as xh
from i2spec.hyperfine import _coefficient_terms, allowed_spins
from i2spec.research import _entry, _read

RESULT = dict(u=42.32383417447296, g=30., gamma=59.59378708522672)      # the release's out-of-core terms


def test_couplings_reproduce_pique_table_ii():
    """Separated-atom couplings against Pique et al. 1986 Table II (printed to the MHz), with the normalization
    factors of TABLE2's comment."""
    rot10 = HBAR2_2U/reduced_mass("127I2")/so.R_TABLE**2*MHZ_PER_CM
    k2 = {0: 4., 1: -4*so.SQRT2, 2: 4*so.SQRT2}
    ours = so.separated_atom()
    for name, (omega, parity, v0, v1, v2) in so.TABLE2.items():
        o = ours[name]
        assert o[2]*rot10 == pytest.approx(v0, abs=3.)
        assert o[3]/so.IA*2*(-so.SQRT2)**omega == pytest.approx(v1, abs=1.)
        assert o[4]*k2[omega] == pytest.approx(v2, abs=1.)


def test_b_quadrupole_element_is_the_eqq_operator():
    J, q = 20, so.COUPLINGS["B"][4]
    worst = 0.
    for F in range(J-5, J+6):
        for I in allowed_spins(J, so.IA, so.IA, "u"):
            for Ip in allowed_spins(J, so.IA, so.IA, "u"):
                for Jp in (J-2, J, J+2):
                    if abs(J-I) <= F <= J+I and abs(Jp-Ip) <= F <= Jp+Ip:
                        ref = _coefficient_terms(Jp, Ip, J, I, F, so.IA, so.IA, 1., 1.)[2]*4*q
                        worst = max(worst, abs(so.element(2, Jp, Ip, J, I, F, 0, "u", 0, q)-ref))
    assert 4*q == pytest.approx(-so.ATOMIC["B32"]/2)
    assert worst < 1e-9


@pytest.mark.parametrize("name", ["0u-", "0g-", "1'g", "1'u"])
def test_projection_equals_the_closed_forms(name):
    """States without a quadrupole coupling: H2 projected onto the four constants equals Chen et al. 2005 eqs 13-16."""
    J, G = 20, -1/50.
    blocks = {F: so.second_order_block(J, F, name, lambda k, kk, Jp: G) for F in range(J-5, J+6)}
    fit, misfit, size = so.project(J, {F: b for F, b in blocks.items() if len(b[0])})
    closed = so.chen_constants(name, lambda k, kk: G)
    assert misfit < 1e-12*size
    for k in ("C", "delta", "d", "eqQ"):
        assert fit[k] == pytest.approx(closed[k], rel=1e-9, abs=1e-6)


def test_frozen_hyperfine_rows_match_the_builders_counts():
    intervals, constants = hyperfine_rows.release()
    assert (len(intervals), sum(r["inside"] for r in intervals), len(constants)) == (1426, 1151, 280)


@pytest.fixture(scope="module")
def table():
    return _read(_entry("balanced_expanded"), "hyperfine.json")


def test_b_constants_at_three_levels_are_the_release_table(table):
    constants = table["B"]["model"]["curve_constants"]
    b = np.array(table["B"]["model"]["background"])
    rows = {(int(r[0]), int(r[1])): r[2:] for r in table["B"]["rows"]}
    levels = hb.Levels()
    picks = [(45, 30), (58, 99), (70, 12)]          # model indices equal v here (no level index shift)
    model = hb.Model(levels, hb.fill_terms(picks, constants, workers=4))
    for v, J in picks:
        assert np.allclose(np.round(model.constants(b, RESULT, v, J), 7), rows[v, J], rtol=0, atol=1e-7)


def test_x_continuation_refit(table):
    lv = xh.XLevels()
    m = xh.fit(lv, xh.KNOT, xh.ORDER)
    release = table["X"]["model"]
    assert m["lam"] == pytest.approx(release["eqQ"]["lam"], rel=1e-5)
    assert m["slope_c"] == pytest.approx(release["C"]["slope"], rel=1e-9)
    rows = {(int(r[0]), int(r[1])): r[2:] for r in table["X"]["rows"]}
    for v, J in ((20, 0), (60, 40), (100, 10)):
        q, c = xh.constants(lv, m, v, J)
        assert q == pytest.approx(rows[v, J][0], abs=1e-3) and c == pytest.approx(rows[v, J][1], abs=1e-6)


@pytest.mark.slow
def test_beyond_core_fit_reproduces_the_release(table):
    rows, params = hb.load_data()
    constants = hb.state_constants(False, hb.FIT_OVERRIDES)
    assert constants == table["B"]["model"]["curve_constants"]
    model = hb.Model(hb.Levels(), hb.fill_terms(hb.needed_levels(rows, params), constants))
    b, theta = hb.fit(model, rows, params, dict(u=30., g=30.), False, keys=list(hb.FIT_KEYS))
    assert theta["u"] == pytest.approx(RESULT["u"], abs=1e-7) and theta["gamma"] == pytest.approx(RESULT["gamma"], abs=1e-6)
    assert np.allclose(b, table["B"]["model"]["background"], atol=1e-6)
    assert model.training["out_cost"] == pytest.approx(74346.0011055406, rel=1e-8)
