"""The Gaussian-process B correction for levels without their own (src/i2spec/gp_corrections.py)."""

import pytest

from i2spec.constants import MHZ_PER_CM
from i2spec.level_corrections import corrections_for
from i2spec.lookup import uncertainty
from i2spec.observations import Line
from i2spec.model import RovibronicModel


@pytest.fixture(scope="module")
def c():
    return corrections_for("i2spec2026n")


def test_it_fills_only_uncorrected_visible_levels(c):
    filled = [v for v in range(0, 60) if c.gp_covers("B", v)]
    assert filled == [7, 12, 13, 19, 27, 29]
    assert not c.gp_covers("X", 7) and corrections_for("i2spec2026m").gp is None


def test_positions_carry_it_and_earlier_sets_do_not(c):
    """R(20) 19-0 moves by exactly the GP mean at B(19, 21); i2spec2026m leaves the level uncorrected."""
    new, old = RovibronicModel("127I2", "i2spec2026n"), RovibronicModel("127I2", "i2spec2026m")
    d = (new.transition(19, 0, 20, "R") - old.transition(19, 0, 20, "R")) * MHZ_PER_CM
    assert d == pytest.approx(c.gp.mean(19, 21), abs=1e-3)
    assert 0.01 < abs(d) < 5


def test_the_lookup_quotes_the_posterior(c):
    line = Line("127I2", "R", 20, 19, 0)
    u, why = uncertainty(line)
    assert why.startswith("Gaussian-process B correction")
    assert 0.5 <= u < 3.0
