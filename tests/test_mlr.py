"""The Morse/Long-Range potential form and the Phase C trial fit of the X state."""

import json
from pathlib import Path

import numpy as np
import pytest

from i2spec.bspline import BSplineSolver
from i2spec.constants import MHZ_PER_CM, reduced_mass
from i2spec.model import RovibronicModel
from i2spec.observations import Predictor, load_dataset, residuals
from i2spec.potentials import MLRPotential, load_potentials

ROOT = Path(__file__).resolve().parents[1]
SOLVER = dict(rmin=2.10, rmax=6.0, h=0.01, order=10, nlev=70)


@pytest.fixture(scope="module")
def published():
    return load_potentials()["X"]


def _load(name, published):
    d = json.loads((ROOT / "data/potentials" / name).read_text())
    return MLRPotential(De=d["De"], Re=d["Re"], C={int(k): v for k, v in d["C"].items()},
                        beta=tuple(d["beta"]), p=d["p"], q=d["q"], Rref=d["Rref"], bo=published)


@pytest.fixture(scope="module", params=["mlr_x_2026a.json", "mlr_x_2026b.json"])
def mlr(request, published):
    """Both trial sets: 2026a fitted to the published levels, 2026b refitted to the observations."""
    return _load(request.param, published)


def test_the_form_holds_its_limits_whatever_beta_does(published):
    """V(Re) = 0 and V -> De - u(R) + u(R)^2/4De by construction: that is the point of the form."""
    for beta in ((2.0,), (1.0, -3.0, 5.0), (0.5, 0.0, 0.0, 12.0)):
        m = MLRPotential(De=published.De, Re=2.6664, C=published.C, beta=beta, p=5, q=3)
        assert m(np.array([m.Re]))[0] == pytest.approx(0.0, abs=1e-9)
        R = np.array([40.0, 60.0, 100.0])
        expected = m.De - m.u(R) + m.u(R) ** 2 / (4 * m.De)
        assert np.abs(m(R) - expected).max() < 1e-7
        # the approach is O(R^-p): at 15 A the same difference is still about 1e-4 cm^-1
        near = np.array([15.0])
        assert 1e-6 < np.abs(m(near) - (m.De - m.u(near) + m.u(near) ** 2 / (4 * m.De)))[0] < 1e-2
    assert MLRPotential(De=1.0, Re=1.0, C={6: 1.0}, beta=(0.0,)).u(np.array([2.0]))[0] == pytest.approx(1 / 64)


def test_a_potential_without_joins_works_everywhere(published, mlr):
    """MLR has no R_I or R_O; the B-spline solver and RovibronicModel must not require them."""
    assert not hasattr(mlr, "RI")
    s = BSplineSolver(mlr, reduced_mass("127I2"), **SOLVER)          # joins=None must not raise
    assert 100.0 < s.energy(0, 0) < 115.0             # X zero point above the minimum, omega_e/2 ~ 107
    m = RovibronicModel("127I2", potentials={"X": mlr})
    assert m.potentials["X"] is mlr and m.potentials["B"] is not None
    # a replaced potential carries no level corrections, so compare with the uncorrected set
    assert m.energy("B", 32, 57) == pytest.approx(RovibronicModel("127I2", "i2spec2026a").energy("B", 32, 57), abs=1e-9)


def test_the_two_trial_sets_agree_where_no_data_exist(published):
    """v'' = 18-47 is unmeasured, so the prediction there must not depend on the low-v'' target used."""
    a, b = _load("mlr_x_2026a.json", published), _load("mlr_x_2026b.json", published)
    sa = BSplineSolver(a, reduced_mass("127I2"), **SOLVER, joins=())
    sb = BSplineSolver(b, reduced_mass("127I2"), **SOLVER, joins=())
    shift = lambda s, v: s.energy(v, 0) - s.energy(0, 0)                     # noqa: E731
    d = np.array([shift(sa, v) - shift(sb, v) for v in range(18, 55)])
    assert np.abs(d).max() < 0.002                                           # cm-1
    assert shift(sa, 42) - shift(BSplineSolver(published, reduced_mass("127I2"), **SOLVER), 42) < -20.0


def test_the_trial_fit_keeps_the_published_low_v_levels(published, mlr):
    """v'' <= 17 is where the published potential is anchored by data; the MLR must not move it."""
    a = BSplineSolver(published, reduced_mass("127I2"), **SOLVER)
    b = BSplineSolver(mlr, reduced_mass("127I2"), **SOLVER, joins=())
    d = np.array([(b.energy(v, 0) - b.energy(0, 0)) - (a.energy(v, 0) - a.energy(0, 0)) for v in range(18)])
    assert np.sqrt(np.mean(d ** 2)) * MHZ_PER_CM < 5.0              # MHz
    assert mlr.De == pytest.approx(published.De, abs=0.128)          # De stayed inside its published uncertainty
    for n, c in published.C.items():
        assert mlr.C[n] == pytest.approx(c, rel=0.02)


def test_the_trial_fit_reaches_the_high_v_measurements(mlr):
    """matyugin2012 (v'' = 48) and nesterenko2019 (v'' = 53, 54) are 234 and 545 GHz out with the
    published X potential; the MLR brings both under 50 MHz without touching the B state."""
    pub, new = Predictor("i2spec2026a"), Predictor(potentials={"X": mlr})     # the published X, not the default
    for name, before in (("matyugin2012", 2.3e5), ("nesterenko2019", 5.4e5)):
        ds = load_dataset(ROOT / "data/observations" / name)
        scale = MHZ_PER_CM if ds.meta["unit"] == "cm-1" else 1.0
        rms = lambda p: float(np.sqrt(np.mean((residuals(p, ds) * scale) ** 2)))    # noqa: E731
        assert rms(pub) > before
        assert rms(new) < 50.0


def test_the_b_state_trial_fit(published):
    """The published B potential is out by up to 2.26 cm-1 above v' = 44 against the atlas; the B MLR
    is within 0.037. Here we only check the parameter file's own invariants and its low-v' levels,
    because the atlas band shifts are a generated file (prototypes/out/, not in the repository)."""
    d = json.loads((ROOT / "data/potentials/mlr_b_2026a.json").read_text())
    B0 = load_potentials()["B"]
    mlr = MLRPotential(De=d["De"], Re=d["Re"], C={int(k): v for k, v in d["C"].items()}, beta=tuple(d["beta"]),
                       p=d["p"], q=d["q"], Rref=d["Rref_factor"] * d["Re"], Te=d["Te"], bo=B0)
    assert mlr.Te + mlr.De == pytest.approx(published.De + 7602.9762, abs=1e-6)   # tied to the atomic splitting
    assert mlr(np.array([mlr.Re]))[0] == pytest.approx(mlr.Te, abs=1e-9)          # V(Re) = Te
    for n, c in B0.C.items():
        assert mlr.C[n] == pytest.approx(c, rel=0.01)                             # long-range values held
    grid = dict(rmin=2.35, rmax=8.0, h=0.01, order=10, nlev=70)
    a = BSplineSolver(B0, reduced_mass("127I2"), **grid)
    b = BSplineSolver(mlr, reduced_mass("127I2"), **grid, joins=())
    d44 = np.array([b.energy(v, J) - a.energy(v, J) for J in (0, 80, 160) for v in range(45)])
    assert np.sqrt(np.mean(d44 ** 2)) * MHZ_PER_CM < 400.0            # MHz: not yet MHz-level, see the doc
    assert abs(b.energy(54, 26) - a.energy(54, 26)) > 1.0             # v' = 54 moves by more than 1 cm-1
