"""The preferred model (src/i2spec/preferred.py): the i2spec model's positions, hyperfine and uncertainties in the
line list, the command line and the explorer, for every isotopologue; i2spec2026o only for lines the model does not
hold. Uses the cached level tables (built once)."""
import numpy as np
import pytest

from i2spec import preferred as pf
from i2spec.constants import MHZ_PER_CM
from i2spec.intensity import MasterLineList
from i2spec.lookup import Catalog, Line, preferred_positions, uncertainty


@pytest.fixture(scope="module")
def table():
    return pf.preferred()


#: The table is built in worker processes, whose eigensolves can differ from this one's in the last bit
#: (1e-11 cm-1 on Linux), and the statistical sigma comes from finite-difference level derivatives, which
#: carry that to about 1e-6 of sigma: far below anything physical, but not exact.
NU_TOL, SIGMA_REL = 1e-9, 1e-4


def test_the_table_is_the_research_model(table):
    m = pf.model()
    for branch, J, vu, vl in (("R", 56, 32, 0), ("P", 13, 43, 0), ("R", 41, 51, 0), ("R", 86, 33, 48), ("R", 3, 80, 0)):
        nu = table.positions([branch], [J], [vu], [vl])[0]
        assert nu == pytest.approx(m.transition(vu, vl, J, branch), abs=NU_TOL)
        sigma, status, cls = table.uncertainties([branch], [J], [vu], [vl])
        u = m.uncertainty(f"{branch}({J}) {vu}-{vl}")
        assert sigma[0] == pytest.approx(u.sigma_MHz, rel=SIGMA_REL)
        assert pf.STATUS[status[0]] == u.status and ("precise", "atlas", "beyond")[cls[0]] == u.line_class


def test_lines_outside_the_domain_are_not_held(table):
    sigma, status, _ = table.uncertainties(["R", "R", "R", "P"], [446, 60, 20, 0], [32, 80, 81, 32], [0, 0, 0, 0])
    assert list(status) == [-1, -1, -1, -1] and np.isnan(sigma).all()      # J > 445, dissociated, v' > 80, P(0)
    assert np.isnan(table.positions(["R"], [20], [32], [109])[0])


def test_lines_above_the_fit_J_are_held_and_extrapolated(table):
    m = pf.model()
    sigma, status, _ = table.uncertainties(["R", "P"], [230, 330], [10, 15], [2, 4])
    assert list(status) == [2, 2] and (sigma > 1).all()                    # extrapolated, and said to be
    assert table.positions(["R"], [230], [10], [2])[0] == pytest.approx(m.transition(10, 2, 230, "R"), abs=NU_TOL)


def test_lookup_takes_the_preferred_line_and_says_so(table):
    line = Line("127I2", "R", 56, 32, 0, 0.0, 0.0, 0.0, 300.0)
    value, why = uncertainty(line)
    assert line.source == "i2spec" and line.flags == ("levels measured",)
    assert 0.1 < value < 1.0 and why.startswith("The fit's statistical uncertainty")
    # the apps present one model: the internal names appear in no flag, source or reason
    other = Line("129I2", "R", 56, 32, 0, 0.0, 0.0, 0.0, 300.0)       # every isotopologue is the model's
    assert other.source == "i2spec" and other.flags == ("levels extrapolated",)
    hot = Line("127I2", "R", 230, 10, 2, 0.0, 0.0, 0.0, 300.0)         # and every J the line lists reach
    assert hot.source == "i2spec" and hot.flags == ("levels extrapolated",)
    far = Line("127I2", "R", 20, 81, 0, 0.0, 0.0, 0.0, 300.0)          # beyond B v' = 80: not the model's
    assert far.source == "i2spec, outside the fitted domain" and far.flags[-1] == "outside the fitted domain"
    for shown in (line, other, hot, far):
        text = " ".join([shown.source, *shown.flags, uncertainty(shown)[1]]).lower()
        assert "2026o" not in text and "balanced" not in text, text


def test_master_lists_take_the_preferred_positions(table):
    nu = np.array([18788.0, 18788.0, 18000.0])
    master = MasterLineList(nu=nu, strength0=np.ones(3), v_upper=np.array([32, 81, 10]), v_lower=np.array([0, 0, 2]),
                            J_lower=np.array([56, 20, 230]), branch=np.array([1, 1, 1]), E_lower=np.zeros(3),
                            x_levels=np.zeros((2, 2)), upper_cut=np.full(2, np.nan), isotopologue="127I2")
    out = preferred_positions(master)
    assert out.nu[0] == table.positions(["R"], [56], [32], [0])[0] and out.nu[1] == 18788.0    # v' 81: kept
    assert out.nu[2] == table.positions(["R"], [230], [10], [2])[0]
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
    assert catalog.transition("127I2", "R", 56, 32, 0) == pytest.approx(pf.model().transition(32, 0, 56, "R"), abs=NU_TOL)


def test_i2spec_model_switch(monkeypatch):
    monkeypatch.setenv("I2SPEC_MODEL", "i2spec2026o")
    line = Line("127I2", "R", 56, 32, 0, 18788.3356, 0.0, 0.0, 300.0)
    assert not pf.active() and line.source == "i2spec2026o" and line.flags == ("corrected levels",)
