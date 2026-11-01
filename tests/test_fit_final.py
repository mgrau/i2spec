"""The rebuild of release 20261031 (i2spec.fit.final): its frozen inputs, and with --slow every stage from them
against the installed release."""
import gzip
import json

import numpy as np
import pytest

from i2spec.fit import calibrate, final, hyperfine_beyond as hb, package as pk


@pytest.fixture(scope="module")
def table():
    return final.installed("expanded", "hyperfine.json")


@pytest.fixture(scope="module")
def unc():
    return final.installed("expanded", "uncertainty.json")


def test_frozen_hyperfine_folds_are_the_release_form(table):
    with gzip.open(final.TRAINING/"hyperfine_folds.json.gz", "rt") as f:
        folds = json.load(f)
    model = table["B"]["model"]
    assert set(folds) == {str(s) for s in hb.FORM_SCALES}
    chosen = folds[str(model["c5_scale"])]
    assert chosen["orders"] == model["orders"] and chosen["fit"] == ",".join(final.FORM["keys"])
    assert hb.state_constants(False, hb.parse_overrides(chosen["override"])) == model["curve_constants"]
    assert [f["holdout"] for f in chosen["folds"]] == [[list(b)] for b in hb.BLOCKS]
    # every fold keeps the release's scale (the nested choice: the lowest training cost)
    cv = final.cross_validation()
    from i2spec.fit import uncertainty as U
    assert {U.chosen_scale(cv, i) for i in range(len(hb.BLOCKS))} == {model["c5_scale"]}
    # the rest is release 20261029's cross-validation, unchanged
    assert cv["core_folds"] == U.cross_validation("20261029")["core_folds"]


def test_frozen_edges_are_the_release_rules(unc):
    rules = final.edges()
    released = json.loads(json.dumps(unc["edges"]))
    d0 = released["X above 108"].pop("sigma_D0_MHz")
    assert rules == released
    assert d0 == pytest.approx(unc["named"]["sigma"]["D0(X)"]*final.CM_TO_MHZ, rel=1e-15)


def test_spike_rule_finds_the_release_levels(table):
    rows = table["B"]["rows"]+table["B"]["spikes"]["rows"]
    assert pk.spike_levels(rows, table["resonance"]) == [(d["v"], d["J"]) for d in table["B"]["not_predicted"]]


def test_second_order_terms_of_the_release():
    """The release's 1g(1Pig) terms at its scale are the release continuation's; every other state's, and the
    1g(1Pig)'s at release 20261030's scale, release 20261021's hyperfine potential's; the cache keeps them apart."""
    form = hb.state_constants(False, hb.overrides_at(.85))
    assert not final.s2_frozen("1g(1Pig)", form["1g(1Pig)"])
    assert all(final.s2_frozen(s, c) for s, c in form.items() if s != "1g(1Pig)")
    assert final.s2_frozen("1g(1Pig)", hb.state_constants(False, hb.FIT_OVERRIDES)["1g(1Pig)"])
    core, cont = final.released()
    from i2spec.fit import base
    with base.based_on(core=core, continuation=cont):
        frozen, current = hb.potential_key(True), hb.potential_key(False)
    assert frozen == hb.potential_key() and frozen != current
    # the frozen terms: every level the fit reads, of the 1g(1Pig) at the fitted scale
    fixed = final.s2_fixed("1g(1Pig)", form["1g(1Pig)"])
    rows, params = hb.load_data()
    assert len(fixed) == 390 and set(hb.needed_levels(rows, params)) <= set(fixed)
    assert final.s2_fixed("b'2u", form["b'2u"]) is None


def test_calibrate_rows_are_the_release_count():
    """The calibrate rows of step 4 (2,432 of them line rows with a model sigma: the slow test)."""
    rows, left_out = calibrate.calibrate_rows()
    assert len(rows) == 2711
    assert left_out["printed again in a later set"] == 359 and left_out["inside a training CIPM mean"] == 115


@pytest.fixture(scope="module")
def rebuilt(tmp_path_factory):
    return final.rebuild(tmp_path_factory.mktemp("release"))


@pytest.mark.slow
def test_rebuild_reproduces_the_release(rebuilt):
    """Every release file from the repository: the same text and structure, and the numbers within the refits'
    convergence. The hyperfine fit beyond the core has a flat optimum (its three starts end 1e-3 deg apart in
    theta_u, 3e-5 apart in a cost of 41 000), and the covariance's weakest directions move with the last digits of
    the coefficients (the release stores them in named coordinates), so those two set the tolerances beyond the core
    and of the uncertainty model."""
    files = rebuilt["files"]
    assert set(files) == {f"{r}/{n}" for r, n in (("core", "model.json"), ("core", "reference_predictions.json"),
                                                  ("expanded", "model.json"), ("expanded", "hyperfine.json"),
                                                  ("expanded", "reference_predictions.json"),
                                                  ("expanded", "uncertainty.json"), ("i2spec", "isotopologues.json"),
                                                  ("i2spec", "reference_predictions.json"))}
    report = final.compare(files)
    assert not [(k, r["other"][:3]) for k, r in report.items() if r["other"]]
    # the core refit: a flat optimum (cost 1e-8 below the release's), every coefficient within 1e-6 of itself
    core, ref = files["core/model.json"], final.installed("core", "model.json")
    for s in ref["values"]:
        for n, z in ref["values"][s].items():
            assert core["values"][s][n] == pytest.approx(z, rel=1e-6, abs=1e-9), (s, n)
    for p, c in ref["hyperfine"]["coefficients"]["B"].items():
        assert np.allclose(core["hyperfine"]["coefficients"]["B"][p], c, rtol=0, atol=1e-6), p
    assert report["core/reference_predictions.json"]["max_abs"] < 1e-6                         # MHz
    # the continuation on the released core: the released tails
    exp, eref = files["expanded/model.json"], final.installed("expanded", "model.json")
    for s in eref["values"]:
        for n, z in eref["values"][s].items():
            assert exp["values"][s][n] == pytest.approx(z, rel=1e-9, abs=1e-9), (s, n)
    assert exp["fit"]["cost"] == pytest.approx(eref["fit"]["cost"], rel=1e-7)
    # the hyperfine beyond the core: the scale, the spike levels, the terms within the fit's flat optimum
    hf, href = files["expanded/hyperfine.json"], final.installed("expanded", "hyperfine.json")
    m, mref = hf["B"]["model"], href["B"]["model"]
    assert (m["c5_scale"], m["curve_constants"]) == (mref["c5_scale"], mref["curve_constants"])
    assert hf["B"]["not_predicted"] == href["B"]["not_predicted"]
    for k, z in mref["theta_deg"].items():
        assert m["theta_deg"][k] == pytest.approx(z, rel=1e-4, abs=2e-3), k
    B, Bref = np.array(hf["B"]["rows"]), np.array(href["B"]["rows"])
    assert np.all(abs(B-Bref) <= 2e-3*abs(Bref)+1e-3)                          # max 9e-4 MHz eqQ, 0.09 kHz delta
    assert np.max(abs(np.array(hf["X"]["rows"])-np.array(href["X"]["rows"]))) < 1e-4
    for v, d in href["resonance"]["levels"].items():
        r = hf["resonance"]["levels"][v]
        assert r["J"] == d["J"] and r["model_level"] == d["model_level"]
        assert max(abs(x-y) for J in d["constants"] for x, y in zip(r["constants"][J], d["constants"][J])) < 1e-3
        assert max(abs(r["overlaps"][k]-z) for k, z in d["overlaps"].items()) < 1e-9
    # references: centers and core components to Hz, components beyond the core within the fit's optimum
    refs, rref = files["expanded/reference_predictions.json"], final.installed("expanded", "reference_predictions.json")
    d = np.array([abs(a["predicted_MHz"]-b["predicted_MHz"]) for a, b in zip(refs, rref)])
    beyond = np.array([a["component"] is not None and a["line"][2] > 43 for a in refs])
    assert d[~beyond].max() < 1e-5 and d[beyond].max() < 2e-3                                  # MHz
    assert report["i2spec/reference_predictions.json"]["max_abs"] < 1e-4
    # the uncertainty model
    u, uref = files["expanded/uncertainty.json"], final.installed("expanded", "uncertainty.json")
    assert u["statistical"]["parameters"] == uref["statistical"]["parameters"]
    assert u["statistical"]["dropped_directions"] == uref["statistical"]["dropped_directions"]
    assert u["coverage"] == uref["coverage"] and u["edges"]["B above 80"] == uref["edges"]["B above 80"]
    F, Fref = np.array(u["statistical"]["factor"]), np.array(uref["statistical"]["factor"])
    C, Cref = F@F.T, Fref@Fref.T
    assert np.max(abs(C-Cref)) < 1e-5*np.max(abs(Cref))
    assert np.allclose(np.sqrt(np.diag(C)), np.sqrt(np.diag(Cref)), rtol=1e-3)
    for group in ("model_error", "hyperfine"):
        for key, ref in uref[group].items():
            for k, z in (ref.items() if isinstance(ref, dict) else ()):
                if isinstance(z, float):
                    assert u[group][key][k] == pytest.approx(z, rel=5e-3, abs=1e-3), (group, key, k)
    near = u["model_error"]["X beyond"]["near_dissociation"]
    for k, z in uref["model_error"]["X beyond"]["near_dissociation"].items():
        assert near[k] == pytest.approx(z, rel=1e-5), k
    assert u["named"]["values"] == pytest.approx(uref["named"]["values"], rel=1e-12)
    assert u["named"]["sigma"] == pytest.approx(uref["named"]["sigma"], rel=5e-3)
    assert rebuilt["step4"]["features"]["n"] == 2432 and rebuilt["step4"]["features"]["max_relative"] < 1e-12
    # the isotopologues
    iso, iref = files["i2spec/isotopologues.json"], final.installed("i2spec", "isotopologues.json")
    assert iso["domains"] == iref["domains"] and iso["hyperfine"] == iref["hyperfine"]
    s, sref = iso["isotope_shift"], iref["isotope_shift"]
    assert s["coefficients"] == pytest.approx(sref["coefficients"], rel=1e-6)
    for k, z in sref["conditioned_line_sigma_MHz"].items():
        assert s["conditioned_line_sigma_MHz"][k] == pytest.approx(z, rel=1e-4), k


@pytest.mark.slow
def test_a_hyperfine_fold_refits_to_its_frozen_result(table):
    with gzip.open(final.TRAINING/"hyperfine_folds.json.gz", "rt") as f:
        frozen = json.load(f)["0.85"]["folds"][0]
    core, cont = final.released()
    out = final.hyperfine_fold(core, cont, .85, frozen["holdout"][0], table["B"]["model"]["theta_deg"])
    for k, z in frozen["theta"].items():
        assert out["theta"][k] == pytest.approx(z, rel=1e-4, abs=1e-3), k
    assert np.allclose(out["b"], frozen["b"], rtol=0, atol=1e-3)
    assert out["training"]["out_cost"] == pytest.approx(frozen["training"]["out_cost"], rel=1e-7)
    held = [(a["resid"], b["resid"]) for x, y in zip(out["eqQ"], frozen["eqQ"]) for a, b in zip(x["rows"], y["rows"])]
    assert len(held) and max(abs(a-b) for a, b in held) < 1e-4                                 # MHz
