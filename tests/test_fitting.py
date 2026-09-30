"""The global fit: residual construction, the analytic Jacobian and conditioning (docs/design/fitting.md)."""

from dataclasses import replace

import numpy as np
import pytest

from i2spec.bspline import BSplineSolver
from i2spec.constants import reduced_mass
from i2spec.fitting import GlobalFit
from i2spec.observations import Predictor, load_dataset, residuals
from i2spec.potentials import load_potentials


@pytest.fixture(scope="module")
def fit():
    return GlobalFit(datasets=[load_dataset("data/observations/velchev1998a"),
                               load_dataset("data/observations/bodermann1998b")])


def test_wavefunctions_are_normalised():
    pots = load_potentials("hannover2008")
    s = BSplineSolver(pots["X"], reduced_mass("127I2"), rmin=2.10, rmax=4.0, h=0.01, order=10, nlev=12)
    energies, psi = s.wavefunctions(30)
    assert np.allclose((s.W[:, None] * psi ** 2).sum(axis=0), 1.0, atol=1e-10)
    assert np.allclose(energies, s.levels(30))


def test_hellmann_feynman_matches_finite_differences():
    """dE/da from <psi|dV/da|psi> must agree with differencing the eigenvalue itself."""
    pots = load_potentials("hannover2008")
    mu = reduced_mass("127I2")
    grid = dict(rmin=2.10, rmax=4.0, h=0.01, order=10, nlev=20)
    base = pots["X"]
    s = BSplineSolver(base, mu, **grid)
    _, psi = s.wavefunctions(40)
    for k in (1, 2, 5):
        # the same step rule as i2spec.fitting.DERIVATIVE_STEP: scaling by |a_k| alone makes the
        # *finite difference* noise-limited for small coefficients, because the eigenvalue then
        # moves only a few times the eigensolver's own round-off
        step = 1e-5 * max(abs(base.a[k]), 1.0)
        def shifted(delta):
            a = list(base.a); a[k] += delta
            return replace(base, a=tuple(a)).with_continuous_extensions()
        fd = (BSplineSolver(shifted(+step), mu, **grid).levels(40)
              - BSplineSolver(shifted(-step), mu, **grid).levels(40)) / (2 * step)
        dV = (shifted(+step)(s.R) - shifted(-step)(s.R)) / (2 * step)
        hf = (s.W[:, None] * dV[:, None] * psi ** 2).sum(axis=0)
        for v in (0, 5, 15):
            if abs(fd[v]) < 1e-7:             # derivative indistinguishable from round-off
                continue
            assert hf[v] == pytest.approx(fd[v], rel=5e-3)


def test_residuals_match_the_observations_predictor(fit):
    """The fit's forward model must agree with the one used to report residuals elsewhere."""
    ds = load_dataset("data/observations/velchev1998a")
    # the fit refits the potentials themselves, so it starts from the set without level corrections
    reference = np.array(residuals(Predictor("i2spec2026a"), ds))          # observed - model, MHz
    mine = fit.residuals(fit.x0) * np.array([d.uncertainty for d in fit.data])
    mine = np.array([m for m, d in zip(mine, fit.data) if d.dataset == "velchev1998a"])
    assert len(mine) == len(reference)
    assert np.allclose(np.sort(-mine), np.sort(reference), atol=2e-3)


def test_intra_line_splittings_are_separated(fit):
    """A splitting within one line cannot constrain the potentials, so it must not enter the fit."""
    assert all(d.terms for d in fit.data)
    assert any(not d.terms for d in fit.hyperfine_only)
    for d in fit.hyperfine_only:
        assert d.terms == ()


def test_energy_origin_is_held(fit):
    """X.a[0] is the energy zero and is degenerate with B.a[0]; it must not be fitted."""
    assert fit.fixed == (("X", 0),)
    assert ("X", "a0") not in fit.knobs
    assert ("X", "a1") in fit.knobs and ("B", "a0") in fit.knobs
    # the parameter vector is [potential knobs, one offset per large data group]
    assert fit.n_parameters == fit.n_knobs + fit.n_offsets
    # whatever the optimiser does, X.a[0] stays at the published zero
    assert fit._potentials(fit.x0 + 1.0)["X"].a[0] == 0.0


def test_long_range_parameters_are_fitted_with_their_own_priors(fit):
    """De and the dispersion coefficients carry the uncertainties of the measurement they come from."""
    from i2spec.fitting import LONG_RANGE_PRIORS
    names = [fit.knobs[i][1] for i in sorted(fit.knob_prior)]
    assert names == ["De", "C6", "C8", "C10"]
    for i, sigma in fit.knob_prior.items():
        assert sigma == LONG_RANGE_PRIORS["X"][fit.knobs[i][1]]
    # nothing at low v'': with only v'' <= 17 data the dispersion terms are invisible
    J = fit.jacobian(fit.x0)[:, :fit.n_knobs]
    i_c6 = next(i for i, (s, n) in enumerate(fit.knobs) if n == "C6")
    sig = np.array([d.uncertainty for d in fit.data])
    assert np.abs(J[:, i_c6] * sig * fit.knob_prior[i_c6]).max() < 1.0    # MHz


def test_long_range_moves_high_v_and_not_low_v():
    """C6 at its own 1 sigma must move v'' = 53 by GHz and v'' <= 17 by nothing.

    This is the degeneracy that blocks a fit of the series coefficients alone: without the
    long-range terms the response of high v'' is 94% collinear with that of low v''.
    """
    g = GlobalFit(datasets=[load_dataset("data/observations/velchev1998a"),
                            load_dataset("data/observations/nesterenko2019")])
    J = g.jacobian(g.x0)[:, :g.n_knobs]
    sig = np.array([d.uncertainty for d in g.data])
    i_c6 = next(i for i, (s, n) in enumerate(g.knobs) if n == "C6")
    effect = np.abs(J[:, i_c6] * sig * g.knob_prior[i_c6])          # MHz per 1 sigma of C6
    high = np.array([d.dataset == "nesterenko2019" for d in g.data])
    assert effect[~high].max() < 1.0            # v'' = 1: untouched
    assert effect[high].min() > 1e3             # v'' = 53, 54: GHz


def test_adaptive_boxes_match_a_single_large_box(fit):
    """Solving each J on the cheapest exact grid must not change any answer."""
    groups, which = fit._boxes()
    assert len(groups) > 1                       # more than one box is actually in use
    assert set(which) == set(fit.levels_needed)
    reference = {}
    for k, vmax in fit.levels_needed.items():
        iso, state, J = k
        grid = dict(rmin=2.10 if state == "X" else 2.35,
                    rmax=6.0 if state == "X" else 8.0, h=0.01, order=10,
                    nlev=70 if state == "X" else 60)
        key = (iso, state, grid["rmax"], grid["nlev"])
        reference.setdefault(key, (grid, []))[1].append(J)
    saved = fit._box_map
    adaptive = np.array(fit.residuals(fit.x0))
    try:
        fit._box_map = (reference, {k: (k[0], k[1], 6.0 if k[1] == "X" else 8.0,
                                        70 if k[1] == "X" else 60) for k in fit.levels_needed})
        single = np.array(fit.residuals(fit.x0))
    finally:
        fit._box_map = saved
    sigma = np.array([d.uncertainty for d in fit.data])
    assert np.abs((adaptive - single) * sigma).max() < 1e-3    # under 1 kHz


def test_analytic_jacobian_matches_numerical(fit):
    """The whole Jacobian, not just one level: predicted residual changes must match actual ones."""
    J = fit.jacobian(fit.x0)
    r0 = fit.residuals(fit.x0)
    rng = np.random.default_rng(1)
    dx = rng.normal(size=fit.n_parameters) * 1e-8 * np.maximum(np.abs(fit.x0), 1.0)
    predicted = J @ dx
    actual = fit.residuals(fit.x0 + dx) - r0
    big = np.abs(actual) > 1e-6
    assert big.sum() > 10
    assert np.allclose(predicted[big], actual[big], rtol=2e-3, atol=1e-6)


def test_prior_makes_the_problem_full_rank(fit):
    """Without the prior the data alone cannot determine every coefficient."""
    fit._set_scale()
    scaled = fit.scaled_jacobian(np.zeros(fit.n_parameters))
    sv = np.linalg.svd(scaled, compute_uv=False)
    assert np.sum(sv > sv[0] * 1e-12) == fit.n_parameters
    # the prior block alone bounds the smallest singular value from below: 1/prior_tau for the free
    # coefficients, scale/sigma for the ones with a measured prior (the long-range constants)
    diag = np.full(fit.n_parameters, 1.0 / fit.prior_tau)
    for i, sigma in fit.knob_prior.items():
        diag[i] = fit.scale[i] / sigma
    assert sv[-1] >= diag.min() * (1 - 1e-9)
    data_only = fit.jacobian(fit.x0) * fit.scale
    sv_data = np.linalg.svd(data_only, compute_uv=False)
    assert np.sum(sv_data > sv_data[0] * 1e-12) < fit.n_parameters
