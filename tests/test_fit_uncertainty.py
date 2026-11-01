"""The fit's uncertainty build (i2spec.fit.uncertainty) reproduces the release's uncertainty model."""
import numpy as np
import pytest

from i2spec.fit import continuation, levels, uncertainty
from i2spec.research import _entry, _read


@pytest.fixture(scope="module")
def release():
    return _read(_entry("balanced_expanded"), "uncertainty.json")


@pytest.fixture(scope="module")
def covariance():
    return uncertainty.Covariance(uncertainty.Jacobian())


def test_statistical_covariance(release, covariance):
    stat = release["statistical"]
    assert covariance.names == stat["parameters"] and len(covariance.names) == 97
    assert covariance.info["dropped"] == stat["dropped_directions"] == 11
    C, Cr = covariance.L@covariance.L.T, np.array(stat["factor"])@np.array(stat["factor"]).T
    assert np.max(abs(C-Cr)) < 1e-5*np.max(abs(Cr))
    assert np.allclose(np.sqrt(np.diag(C)), np.sqrt(np.diag(Cr)), rtol=1e-4)
    for c, b in covariance.ratios.items():
        assert b["birge"] == pytest.approx(stat["birge"][c]["birge"], rel=1e-5)
    assert {covariance.names[j]: s for j, s in covariance.prior.items()} == stat["priors"]


def test_coverage_tables(release, covariance):
    every, precise = uncertainty.coverage_tables(covariance.jac.rows)
    assert every == release["coverage"]["all"] and precise == release["coverage"]["precise"]


@pytest.mark.slow
def test_model_error_and_hyperfine_sigma(release, monkeypatch):
    # the release this reproduces predates the barrier rule on short grids
    monkeypatch.setattr(levels, "ABOVE_LIMIT_ON_LONG_GRID", False)
    record, _ = uncertainty.build()
    for key, ref in release["model_error"].items():
        for k, z in ref.items():
            if isinstance(z, float):
                assert record["model_error"][key][k] == pytest.approx(z, rel=1e-4, abs=1e-6), (key, k)
    for key in ("core", "B beyond", "X beyond"):
        for k, z in release["hyperfine"][key].items():
            if isinstance(z, float):
                assert record["hyperfine"][key][k] == pytest.approx(z, rel=1e-6), (key, k)


@pytest.mark.slow
def test_a_continuation_fold_refits_to_its_frozen_predictions():
    frozen = uncertainty.cross_validation()["continuation_folds"][0]
    ids, pred, observed, result = continuation.fold(0)
    assert list(ids) == frozen["ids"]
    assert result.cost == pytest.approx(frozen["cost"], rel=1e-9)
    sigma = dict(zip(continuation.training.release("continuation").id,
                     continuation.training.release("continuation").sigma))
    assert np.max(abs(pred-np.array(frozen["pred"]))/np.array([sigma[k] for k in ids])) < 1e-2


@pytest.mark.slow
def test_a_hyperfine_fold_refits_to_its_frozen_result():
    frozen = uncertainty.cross_validation()["hyperfine_folds"]["0.65"]["folds"][0]
    b, theta, training = uncertainty.hyperfine_fold(.65, frozen["holdout"][0])
    assert theta["u"] == pytest.approx(frozen["theta"]["u"], abs=1e-7)
    assert np.allclose(b, frozen["b"], atol=1e-6)
    assert training["out_cost"] == pytest.approx(frozen["training"]["out_cost"], rel=1e-8)
