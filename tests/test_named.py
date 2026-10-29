"""Named coordinates (Be, T0, D0): an exact change of variables (i2spec.named)."""
import json
from importlib import resources

import numpy as np
import pytest

from i2spec import named


@pytest.fixture(scope="module")
def tied():
    spec = json.loads(resources.files("i2spec").joinpath("data", "research", "20261018", "balanced", "model.json")
                      .read_text())
    spec["asymptote"] = "tied"
    named.tie(spec)
    return spec


def test_named_coordinates_return_the_coefficients(tied):
    back = named.to_native(named.to_named(tied))
    for state, keys in (("X", ("logRe", "logDe")), ("B", ("logRe", "logDe", "Te"))):
        for k in keys:
            assert back["values"][state][k] == pytest.approx(tied["values"][state][k], abs=1e-9)
    assert back["named"]["B"]["T0"] == pytest.approx(named.ground(tied, "B")-named.ground(tied, "X"), abs=1e-9)


def test_the_tie_puts_the_B_limit_at_D0_plus_the_atomic_splitting(tied):
    from i2spec.constants import ATOMIC_SPLITTING
    v = tied["values"]
    limit_b = v["B"]["Te"]+np.exp(v["B"]["logDe"])
    assert limit_b-np.exp(v["X"]["logDe"]) == pytest.approx(ATOMIC_SPLITTING, abs=1e-9)


def test_a_continued_model_stores_D0_of_the_whole_model(tied):
    core_named = named.named(tied)
    spec = dict(values={"X": {"dDe": .79}, "B": {"dDe": .79}}, structure={"de_unit": .01})
    out = named.continued_to_named(spec, core_named)
    assert out["values"]["X"]["D0"] == pytest.approx(core_named["X"]["D0"]+.0079)
    assert named.continued_to_native(out, core_named)["values"]["B"]["dDe"] == pytest.approx(.79)


def test_named_coordinates_need_the_tie():
    with pytest.raises(ValueError):
        named.to_named(dict(values={}, config={}))
