"""The preferred model (src/i2spec/preferred.py): balanced-expanded's positions, hyperfine and uncertainties in the
line list, the command line and the explorer, i2spec2026o elsewhere. Uses the cached level table (built once)."""
import numpy as np
import pytest

from i2spec import preferred as pf
from i2spec.constants import MHZ_PER_CM
from i2spec.intensity import MasterLineList
from i2spec.lookup import Catalog, Line, preferred_positions, uncertainty


@pytest.fixture(scope="module")
def table():
    return pf.preferred()


def test_the_table_is_the_research_model(table):
    m = pf.model()
    for branch, J, vu, vl in (("R", 56, 32, 0), ("P", 13, 43, 0), ("R", 41, 51, 0), ("R", 86, 33, 48), ("R", 3, 80, 0)):
        nu = table.positions([branch], [J], [vu], [vl])[0]
        assert nu == m.transition(vu, vl, J, branch)
        sigma, status, cls = table.uncertainties([branch], [J], [vu], [vl])
        u = m.uncertainty(f"{branch}({J}) {vu}-{vl}")
        assert sigma[0] == pytest.approx(u.sigma_MHz, rel=1e-6)
        assert pf.STATUS[status[0]] == u.status and ("precise", "atlas", "beyond")[cls[0]] == u.line_class


def test_lines_outside_the_domain_are_not_held(table):
    sigma, status, _ = table.uncertainties(["R", "R", "R", "P"], [201, 60, 20, 0], [32, 80, 81, 32], [0, 0, 0, 0])
    assert list(status) == [-1, -1, -1, -1] and np.isnan(sigma).all()      # J > 200, dissociated, v' > 80, P(0)
    assert np.isnan(table.positions(["R"], [20], [32], [109])[0])


def test_lookup_takes_the_preferred_line_and_says_so(table):
    line = Line("127I2", "R", 56, 32, 0, 0.0, 0.0, 0.0, 300.0)
    value, why = uncertainty(line)
    assert line.source == "balanced-expanded" and line.flags == ("balanced-expanded", "levels measured")
    assert 0.1 < value < 1.0 and why.startswith("balanced-expanded")
    other = Line("129I2", "R", 56, 32, 0, 0.0, 0.0, 0.0, 300.0)
    assert other.source == "i2spec2026o" and "i2spec2026o" in other.flags
    hot = Line("127I2", "R", 230, 10, 2, 0.0, 0.0, 0.0, 300.0)
    assert hot.source == "i2spec2026o" and "i2spec2026o: outside balanced-expanded" in hot.flags


def test_master_lists_take_the_preferred_positions(table):
    nu = np.array([18788.0, 18788.0, 18000.0])
    master = MasterLineList(nu=nu, strength0=np.ones(3), v_upper=np.array([32, 81, 10]), v_lower=np.array([0, 0, 2]),
                            J_lower=np.array([56, 20, 230]), branch=np.array([1, 1, 1]), E_lower=np.zeros(3),
                            x_levels=np.zeros((2, 2)), upper_cut=np.full(2, np.nan), isotopologue="127I2")
    out = preferred_positions(master)
    assert out.nu[0] == table.positions(["R"], [56], [32], [0])[0] and list(out.nu[1:]) == [18788.0, 18000.0]
    assert preferred_positions(out) is out and list(out.nu_2026o) == list(nu)


def test_catalog_components_carry_their_uncertainty(table):
    catalog = Catalog()
    line = Line("127I2", "R", 56, 32, 0, 0.0, 0.0, 0.0, 300.0)
    comps = catalog.components(line, main_only=True)
    center = pf.model().transition(32, 0, 56, "R")*MHZ_PER_CM
    assert len(comps) == 15 and all(c.uncertainty_MHz > uncertainty(line)[0] for c in comps)
    a10 = next(c for c in comps if c.label == "a10")
    assert a10.frequency_MHz == pytest.approx(pf.model().frequency("R(56) 32-0", "a10"), abs=1e-6)
    assert a10.frequency_MHz-a10.offset_MHz == pytest.approx(center, abs=1e-6)
    assert catalog.transition("127I2", "R", 56, 32, 0) == pf.model().transition(32, 0, 56, "R")


def test_i2spec_model_switch(monkeypatch):
    monkeypatch.setenv("I2SPEC_MODEL", "i2spec2026o")
    line = Line("127I2", "R", 56, 32, 0, 18788.3356, 0.0, 0.0, 300.0)
    assert not pf.active() and line.source == "i2spec2026o" and line.flags == ("corrected levels",)
