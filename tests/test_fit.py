"""The fit (i2spec.fit) reproduces the released coefficients from the release's own training rows."""
import numpy as np
import pytest

from i2spec.fit import continuation, core, core_hyperfine, rows
from i2spec.fit.curves import Config, Curve, initial


def test_training_rows_are_frozen_with_both_stages():
    core_rows = rows.release("core")
    cont = rows.release("continuation")
    assert (len(core_rows), len(rows.release("core", "benchmark"))) == (41788, 81)
    assert (len(cont), len(rows.release("continuation", "benchmark"))) == (46964, 81)
    assert not core_rows.new.any() and rows.release("core", "benchmark").new.all()
    assert all(all(k[0] == "127I2" for k, _ in lv) for lv in cont.levels)
    # the core rows enter the continuation with the same values
    value = dict(zip(cont.id, cont.value))
    assert all(value[i] == z for i, z in zip(core_rows.id, core_rows.value))


@pytest.fixture(scope="module")
def core_problem():
    spec = core.release_spec()
    config = Config.from_spec(spec["config"])
    values, names, *_ = initial(config)
    prob = core.Problem(rows.release("core"), config, values, names, core.Settings(), core.shape_nodes())
    prob.base = core.start(names, spec)
    return prob


def test_core_cost_at_the_released_coefficients(core_problem):
    cost = core_problem.cost(core_problem.base[core_problem.active])
    assert cost["cost"] == pytest.approx(5165.943842402221, rel=1e-8)      # the release's own fit
    assert cost["physical"] == pytest.approx(2.266087174307829e-4, rel=1e-6)


def test_core_hyperfine_at_the_released_coefficients():
    spec = core.release_spec()
    config = Config.from_spec(spec["config"])
    curves = {s: Curve(s, config, spec["values"][s]) for s in ("X", "B")}
    records = core_hyperfine.records()
    assert (len(records), sum(not r["future"] for r in records)) == (1151, 559)
    energies, means = core_hyperfine.moments(curves, spec["values"], [r["line"] for r in records])
    prob = core_hyperfine.Problem(energies, means, records)
    theta = core_hyperfine.release_coefficients(spec, prob.names)
    out = core_hyperfine.fit(prob, start=theta)
    assert out.cost == pytest.approx(2686.183859376119, rel=1e-8)
    assert np.max(abs(out.theta-theta)) < 1e-9


def test_continuation_is_converged_at_the_released_coefficients():
    spec = continuation.release_spec()
    curves = continuation.build_curves()
    every = [(s, n) for s, c in curves.items() for n in c.tail_names]
    for s, n in every:
        curves[s].values[n] = spec["values"][s][n]
    prob = continuation.Problem(rows.release("continuation"), curves, every, tie={("B", "dDe"): ("X", "dDe")})
    x = continuation.release_start(prob, spec)
    out = continuation.gauss_newton(prob, x, 4)
    assert out["converged"] and np.array_equal(out["x"], x)
    assert out["cost"] == pytest.approx(7235.011561544405, rel=1e-8)


def _same(frozen, rebuilt, tolerance):
    from i2spec.fit import reduction
    d = reduction.compare(frozen, rebuilt, tolerance=tolerance)
    assert not d["only_a"] and not d["only_b"], (d["only_a"][:5], d["only_b"][:5])
    assert not d["changed"], d["changed"][:5]


def test_core_rows_rebuild_from_the_data_with_the_model_hyperfine():
    """The data files, reduced with the released core hyperfine, give the frozen core rows to 1 kHz (the frozen
    centers came from the previous round of the potential-hyperfine alternation, 50 Hz away)."""
    from i2spec.fit import reduction
    core_train, core_bench, *_ = reduction.release_rows(beyond=False)
    _same(rows.release("core"), core_train, 1e-3)
    _same(rows.release("core", "benchmark"), core_bench, 1e-3)


@pytest.mark.slow
def test_continuation_rows_rebuild_from_the_data_with_the_model_hyperfine():
    from i2spec.fit import reduction
    *_, cont_train, cont_bench, calibration, _ = reduction.release_rows()
    _same(rows.release("continuation"), cont_train, 1e-3)
    _same(rows.release("continuation", "benchmark"), cont_bench, 1e-3)
