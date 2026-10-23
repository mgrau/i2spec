"""The physical coordinates (i2spec.fit.physical): exact Taylor coefficients and an exact round trip."""
import numpy as np
import pytest

from i2spec.fit import core, physical
from i2spec.fit.curves import Config, Curve


@pytest.fixture(scope="module")
def curves():
    spec = core.release_spec()
    config = Config.from_spec(spec["config"])
    return {s: (Curve(s, config, spec["values"][s]), spec["values"][s]) for s in ("X", "B")}


@pytest.mark.parametrize("state", ["X", "B"])
def test_jet_expansion_matches_a_fit_of_the_curve(curves, state):
    curve, values = curves[state]
    re = np.exp(values["logRe"])
    A, q = physical.expansion(curve, values)
    fit = physical._taylor(lambda R: curve.potential_values(values, R)-values.get("Te", 0.), re)
    assert abs(A[1]) < 1e-9*A[2]                       # V'(Re) = 0: Re is the minimum
    assert np.allclose(A[2:5], fit[2:5], rtol=1e-8)


def test_mechanical_constants_are_the_levels_dunham_constants(curves):
    """Within the Dunham corrections and the level fit's truncation (v <= 12)."""
    m = physical.mechanical(*curves["X"])
    assert m["Be"] == pytest.approx(0.03736843, rel=1e-5)
    assert m["omega_e"] == pytest.approx(214.52140, rel=1e-5)
    assert m["omega_e_x_e"] == pytest.approx(0.608449, rel=1e-4)
    assert m["alpha_e"] == pytest.approx(1.13881e-4, rel=1e-4)


@pytest.mark.parametrize("state", ["X", "B"])
def test_round_trip_is_exact(curves, state):
    curve, values = curves[state]
    p = physical.to_physical(curve, values)
    start = dict(values, logRe=values["logRe"]+1e-3, c0=values["c0"]+.01, c1=values["c1"]-.01, c2=values["c2"]+.01)
    back = physical.from_physical(curve, p, start)
    assert max(abs(back[k]-values[k]) for k in physical.NATIVE) < 1e-13
    assert all(back[k] == values[k] for k in values if k not in physical.NATIVE)
