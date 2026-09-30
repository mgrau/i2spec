"""RKR inversion of the Gerstenkorn & Luc constants, against the published potentials."""
import numpy as np
import pytest
from scipy.optimize import brentq

from i2spec import dunham, rkr
from i2spec.potentials import load_potentials


@pytest.fixture(scope="module")
def published():
    return load_potentials("hannover2008")


def _turning_points(p, E, r_min, r_max, re):
    return (brentq(lambda r: p(r) - p(re) - E, r_min, re), brentq(lambda r: p(r) - p(re) - E, re, r_max))


def test_rkr_reproduces_the_published_x_potential(published):
    """The X potential and the atlas's Dunham constants describe the same data, so the RKR turning
    points of the one must land on the other: they agree to 0.05 mA at every v'' the constants cover."""
    p = published["X"]
    re = brentq(lambda r: (p(r + 1e-4) - p(r - 1e-4)) / 2e-4, 2.4, 3.2)
    assert rkr.re("X") == pytest.approx(re, abs=1e-4)
    for v in (0, 5, 10, 19):
        ri, ro = rkr.turning_points("X", v)
        E = dunham.G("X", v) - dunham.G("X", -0.5)
        lo, hi = _turning_points(p, E, 2.1, 5.0, re)
        assert ri == pytest.approx(lo, abs=2e-4)
        assert ro == pytest.approx(hi, abs=2e-4)


def test_rkr_b_state_parts_from_the_published_potential_where_it_fails(published):
    """The B state agrees to 0.6 mA up to v' = 30 and parts above v' = 44, where the published curve
    is known to be wrong (docs/design/mlr-x.md)."""
    p = published["B"]
    re = brentq(lambda r: (p(r + 1e-4) - p(r - 1e-4)) / 2e-4, 2.8, 3.5)
    for v, tol in ((0, 6e-4), (10, 6e-4), (30, 8e-4)):
        ri, ro = rkr.turning_points("B", v)
        E = dunham.G("B", v) - dunham.G("B", -0.5)
        lo, hi = _turning_points(p, E, 2.4, 7.5, re)
        assert ri == pytest.approx(lo, abs=tol) and ro == pytest.approx(hi, abs=tol)
    ri, ro = rkr.turning_points("B", 60)
    E = dunham.G("B", 60) - dunham.G("B", -0.5)
    lo, hi = _turning_points(p, E, 2.4, 7.5, re)
    assert abs(ro - hi) > 3e-3                      # 4.5 mA apart at v' = 60


def test_the_curve_is_monotonic_and_brackets_re():
    R, V = rkr.curve("X", v_max=19, n_v=20)
    assert np.all(np.diff(R) > 0)
    assert R.min() < rkr.re("X") < R.max()
    assert V.min() >= 0.0
