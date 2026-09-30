"""Stage 2 of the refit: corrections to the published hyperfine formulae (src/i2spec/hyperfine_fit.py)."""

import numpy as np
import pytest

from i2spec.hfs_params import NO_CORRECTIONS, Corrections
from i2spec.hyperfine_fit import MODEL_FLOOR, HyperfineFit, rms
from i2spec.model import RovibronicModel
from i2spec.observations import load_dataset


@pytest.fixture(scope="module")
def bipm532():
    return HyperfineFit(datasets=[load_dataset("data/observations/bipm2012a")])


def test_zero_corrections_reproduce_the_published_formulae():
    m = RovibronicModel("127I2")
    _, plain = m.hyperfine_components(32, 0, 56, "R", table=None)
    _, zero = m.hyperfine_components(32, 0, 56, "R", corrections=NO_CORRECTIONS, table=None)
    assert [c.offset for c in plain] == [c.offset for c in zero]
    _, moved = m.hyperfine_components(32, 0, 56, "R", corrections=Corrections(eqQ_b=0.1), table=None)
    assert max(abs(a.offset - b.offset) for a, b in zip(plain, moved)) > 1e-3


def test_rms_is_not_the_standard_deviation():
    assert rms([3.0, 3.0, 3.0]) == 3.0      # a pure offset has zero scatter but is not zero error
    assert rms([]) == 0.0


def test_splittings_only_and_units(bipm532):
    """Only intra-line intervals are used, the line-to-line ones carry the rovibronic centre."""
    assert all(d.dataset == "bipm2012a" for d in bipm532.data)
    assert len(bipm532.data) == 309                      # of 329 rows: 1 frequency, 19 inter-line intervals
    assert len(bipm532.lines) == 20
    d = bipm532.data[0]
    assert d.weight == pytest.approx(1 / np.hypot(d.uncertainty, MODEL_FLOOR))


def test_the_best_line_sets_the_model_floor(bipm532):
    """R(56) 32-0, the CIPM 532 nm line, is the best-reproduced line in the set, at ~25 kHz rms.

    That is the demonstrated accuracy of the formulae, and MODEL_FLOOR is taken from it.
    """
    by_line = bipm532.rms_by_line(bipm532.x0)
    best = min(by_line, key=by_line.get)
    assert by_line[("127I2", "R", 56, 32, 0)] == pytest.approx(24.7, abs=0.5)
    assert by_line[best] < 30.0
    assert MODEL_FLOOR * 1e3 == pytest.approx(by_line[("127I2", "R", 56, 32, 0)], rel=0.05)


def test_residuals_grow_toward_the_dissociation_limit(bipm532):
    """The formulae degrade steeply as the B level approaches dissociation, rotation included.

    R(145) 37-0 lies at 19 680 cm-1, 180 cm-1 above the E_b < 19 500 cm-1 validity limit once rotation is
    counted, and is out by 1.2 MHz rms against 2 kHz measurements. Measurement-weighted, that one line
    was 91% of chi-squared for all 642 splittings, and it is what drove the first Stage 2 fit wrong.
    """
    by_line = bipm532.rms_by_line(bipm532.x0)
    core = by_line[("127I2", "R", 56, 32, 0)]
    corner = by_line[("127I2", "R", 145, 37, 0)]
    assert corner > 40 * core
    m = RovibronicModel("127I2")
    assert m.energy("B", 37, 146) > 19500.0 > m.energy("B", 32, 57)

    raw = (bipm532.deviations(bipm532.x0) / bipm532.sigma) ** 2
    share = sum(c for c, d in zip(raw, bipm532.data) if d.J_lower == 145) / raw.sum()
    assert share > 0.9


def test_per_line_fit_brings_the_worst_line_to_its_noise(bipm532):
    """Freeing the B-state parameters of R(145) 37-0 alone takes it from 1.2 MHz to ~10 kHz: the
    Hamiltonian reproduces the pattern, and the published parameter formulae are what miss it.
    C_B carries most of that by itself."""
    key = ("127I2", "R", 145, 37, 0)
    four = bipm532.per_line_fit(key)
    c_only = bipm532.per_line_fit(key, terms=("C_b",))
    assert four["rms_before"] > 1000 and four["rms_after"] < 15
    assert c_only["rms_after"] < 50
    assert four["x"]["C_b"] == pytest.approx(1.42, abs=0.05) and four["error"]["C_b"] < 0.05


def test_the_linearised_problem_is_exact(bipm532):
    """cross_validate and linear_fit work on the Jacobian; check it against the full calculation.

    Up to half a prior step on three terms at once the two agree exactly. At a full step, C_B moves by
    about 6 kHz at J' = 146 -- four times the fitted correction -- and two components of R(145) 37-0
    swap order, which the rank matching reads as an 840 kHz jump. What matters is agreement at the
    coefficients a fit actually produces, and there it is exact too.
    """
    fit = HyperfineFit(datasets=[load_dataset("data/observations/bipm2012a")], free=("C_b", "C_b_y2", "d_b_y"))
    J = fit.jacobian()
    dev0 = fit.deviations(fit.x0)
    step = 0.5 * fit.prior
    assert rms((fit.deviations(step) - (dev0 + J @ step)) * 1e3) < 0.01          # kHz
    assert rms((J @ step) * 1e3) > 100                       # while the step itself moves things a lot
    coef, linear = fit.linear_fit(list(fit.free), J, dev0)
    full = fit.deviations(np.array([coef[t] for t in fit.free]))
    assert rms((full - linear) * 1e3) < 0.1


def test_cross_validation_holds_out_whole_groups(bipm532):
    fit = HyperfineFit(datasets=[load_dataset("data/observations/bipm2012a")], free=("C_b", "C_b_y"))
    J, dev0 = fit.jacobian(), fit.deviations(fit.x0)
    groups = fit.groups("v_upper")
    assert len(np.unique(groups)) == 6                        # v' = 32 ... 37
    assert np.array_equal(fit.cross_validate([], groups, J, dev0), dev0)
    _, in_sample = fit.linear_fit(["C_b", "C_b_y"], J, dev0)
    held_out = fit.cross_validate(["C_b", "C_b_y"], groups, J, dev0)
    assert rms(in_sample) < rms(dev0)                         # a fit always helps the data it saw
    assert rms(held_out) > rms(in_sample)                     # and predicts unseen groups worse
    with pytest.raises(ValueError):
        HyperfineFit(datasets=[], free=("not_a_term",))
