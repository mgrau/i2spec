"""The model's frozen releases: installed assets, selection, and numerical API contracts."""
from dataclasses import FrozenInstanceError
import json

import numpy as np
import pytest

from i2spec import ResearchModel
from i2spec import cli, hfs_params
from i2spec.constants import DEFAULT_PARAMETERS, MHZ_PER_CM
from i2spec.research import _Potential, _entries, _entry, _info, _read


def test_releases_and_ledgers():
    models = [_info(entry) for entry in _entries()]
    assert [m.role for m in models] == ["balanced", "balanced_expanded", "i2spec"]*4
    assert [(m.inventory.fitted_molecular, m.inventory.inherited_molecular, m.inventory.fitted_calibration)
            for m in models] == ([(61, 18, 4), (105, 73, 5), (107, 75, 5)]+[(63, 18, 4), (106, 73, 5), (106, 75, 5)]*2
                                 + [(63, 18, 4), (110, 73, 5), (110, 75, 5)])
    # release 20261031 (version 1.0.0) is the released model: the role "i2spec" resolves to it; the parts' roles
    # resolve to the first release that has them (release 20261021's, which the fit rebuilds)
    assert [_entry(m.role)["id"] for m in models[:2]] == [m.id for m in models[:2]]
    assert [m.id for m in models[3:]] == [f"research-{r}-{n}" for r in ("20261029", "20261030", "20261031")
                                          for n in ("core", "expanded", "i2spec")]
    assert DEFAULT_PARAMETERS == _entry("i2spec")["id"] == "research-20261031-i2spec"
    assert all(m.inventory.line_specific == m.inventory.level_specific == 0 for m in models)
    with pytest.raises(FrozenInstanceError):
        models[0].inventory.fitted_molecular = 1


def test_all_packaged_source_artifacts_match_frozen_snapshot_hashes():
    for entry in _entries():
        assert entry["files_sha256"] == entry.get("source_snapshot_sha256", entry["files_sha256"])
        for name in entry["files_sha256"]:
            _read(entry, name)


def test_artifact_tampering_is_rejected(tmp_path, monkeypatch):
    import i2spec.research as research

    folder = tmp_path/"role"
    folder.mkdir()
    (folder/"model.json").write_text('{}')
    monkeypatch.setattr(research, "_data", lambda release=None: tmp_path)
    with pytest.raises(ValueError, match="SHA256"):
        _read(dict(id="test", directory="role", files_sha256={"model.json": "bad"}), "model.json")


def test_unknown_model_and_wrong_isotope_do_not_fall_back():
    with pytest.raises(ValueError, match="unknown research model"):
        ResearchModel("unrecognized")
    with pytest.raises(ValueError, match="only 127I2"):
        ResearchModel("balanced", "129I2")
    assert ResearchModel("balanced").info.id == "research-20261018-balanced"


def test_rovibronic_model_rejects_a_core_fit():
    from i2spec.model import RovibronicModel
    with pytest.raises(ValueError, match="core fit.*research-20261031-i2spec"):
        RovibronicModel("127I2", "research-20261030-core")


@pytest.mark.parametrize("line", ["P(0) 32-0", "R(201) 32-0", "R(20) 44-0", "R(20) 32-18",
                                  ("R", 2.5, 32, 0), ("R", True, 32, 0), ("Q", 20, 32, 0)])
def test_outside_domain_is_rejected_before_radial_work(line, monkeypatch):
    model = ResearchModel("balanced")
    monkeypatch.setattr(model, "_moments", lambda *a: pytest.fail("unexpected radial solve"))
    with pytest.raises(ValueError):
        model.frequency(line)


def test_inspecting_coefficients_cannot_change_predictions():
    model = ResearchModel("balanced")
    coefficients = model.coefficients
    original = coefficients["values"]["X"]["c0"]
    coefficients["values"]["X"]["c0"] = 999.
    assert model.coefficients["values"]["X"]["c0"] == original
    assert model.info.inventory.line_specific == model.info.inventory.level_specific == 0


def test_excluded_line_is_marked_as_diagnostic(monkeypatch):
    model = ResearchModel("balanced")
    monkeypatch.setattr(model, "transition", lambda *args: 1.)
    prediction = model.predict("P(13) 43-0")
    assert prediction.excluded_from_fit
    assert prediction.frequency_MHz == MHZ_PER_CM
    assert not model.predict("P(53) 32-0").excluded_from_fit


def test_frozen_fixed_x_formula_and_anchored_ground_values(monkeypatch):
    fixed = ResearchModel("balanced")
    monkeypatch.setattr(fixed, "_energy", lambda state, v, J: float(v*100))
    actual = fixed._couplings("X", 4, 53)
    expected = hfs_params.x_state(4, 53, 400.)
    np.testing.assert_allclose([actual.eqQ, actual.C, actual.d, actual.delta],
                              [expected.eqQ, expected.C, expected.d, expected.delta], atol=1e-12, rtol=0.)


def test_shared_potential_minimum_and_rotation_are_preserved():
    model = ResearchModel("balanced")
    spec = model.coefficients
    for state in ("X", "B"):
        potential = _Potential(state, spec)
        re = np.exp(spec["values"][state]["logRe"])
        assert float(potential(np.array([re]))[0]) == spec["values"][state].get("Te", 0.)
        assert np.all(1+potential.nonadiabatic(np.array([2.2, 3., 6.])) > 0)


@pytest.mark.parametrize("role", ["balanced"])
def test_public_api_matches_frozen_absolute_and_interval_references(role):
    model = ResearchModel(role)
    references = _read(_entry(role), "reference_predictions.json")
    absolute = next(r for r in references if "ref_component" not in r)
    interval = next(r for r in references if "ref_component" in r)
    assert abs(model.frequency(absolute["line"], absolute["component"])-absolute["predicted_MHz"]) < 100e-6
    value = model.frequency(interval["line"], interval["component"])
    value -= model.frequency(interval["line"], interval["ref_component"])
    assert abs(value-interval["predicted_MHz"]) < 100e-6
    branch, J, vu, vl = absolute["line"]
    center, components = model.hyperfine_components(vu, vl, J, branch)
    assert center/model.transition(vu, vl, J, branch) == pytest.approx(MHZ_PER_CM, abs=1e-10)
    component = next(c for c in components if c.label == absolute["component"])
    assert center+component.offset == model.frequency(absolute["line"], absolute["component"])
    component.offset = 999.
    _, fresh = model.hyperfine_components(vu, vl, J, branch)
    assert next(c for c in fresh if c.label == absolute["component"]).offset != 999.


def test_cli_model_info(capsys):
    assert cli.main(["model", "--json"]) == 0
    value = json.loads(capsys.readouterr().out)
    assert value["role"] == "i2spec" and value["domain"]["J_lower_max"] == 445
    assert cli.main(["model"]) == 0
    assert "none specific to a line or level" in capsys.readouterr().out
    assert cli.build_parser().parse_args([]).func is cli.cmd_tui


def test_step4_uncertainty_terms_only_in_release_20261030():
    """Release 20261030 is release 20261029 with the uncertainty refitted: the same frequencies, the two terms beyond
    the core (the B hyperfine floor's growth, B sigma_model's (v - v0)^2) present and >= 0, absent before."""
    old, new = ResearchModel("research-20261029-i2spec"), ResearchModel("research-20261030-i2spec")
    t_old, t_new = old._uncertainty_model().table, new._uncertainty_model().table
    assert "quadratic_per_v2" not in t_old["model_error"]["B beyond"]
    assert "a_growth_per_v" not in t_old["hyperfine"]["B beyond"]
    assert t_new["model_error"]["B beyond"]["quadratic_per_v2"] >= 0
    assert t_new["hyperfine"]["B beyond"]["a_growth_per_v"] >= 0
    for line in ("R(56) 32-0", "P(53) 65-0"):
        assert new.frequency(line) == old.frequency(line)
    # the floor beyond the core grows with v'; the core's is one value
    a = [new.uncertainty(f"R(20) {v}-0", "a2", "a1").hyperfine_MHz for v in (45, 60)]
    assert a[1] > a[0]
    b = new._uncertainty_model().hyperfine_terms(new._spec["core"], 50, 0)[0]
    g = t_new["hyperfine"]["B beyond"]
    assert b == pytest.approx(g["a_MHz"]*np.exp(g["a_growth_per_v"]*6))


def test_release_20261031_barrier_rule_on_every_grid():
    """Release 20261031 solves a level that a short box puts above the limit on the long box (its barrier rule):
    release 20261030 counted standing waves of the 8.6 A box as B v' 46-59 at J' 132-200. At J' = 199 the barrier
    holds B to v' = 45 (133.03 cm-1 above the limit); release 20261030 listed box states as v' = 46-59."""
    old, new = ResearchModel("research-20261030-i2spec"), ResearchModel("research-20261031-i2spec")
    lo, ln = old._levels_x, new._levels_x
    assert ln.energy("B", 45, 199)-ln.asymptote["B"] == pytest.approx(133.03, abs=0.01)
    assert ln.energy("B", 45, 199) == pytest.approx(lo.energy("B", 45, 199), abs=1e-9)
    with pytest.raises(ValueError):
        ln.energy("B", 46, 199)
    assert lo.energy("B", 59, 199) > lo.asymptote["B"]          # the box states the old release held
    for v, J in ((30, 100), (50, 120), (60, 100), (70, 20)):      # bound and long-grid levels: unchanged
        assert ln.energy("B", v, J) == pytest.approx(lo.energy("B", v, J), abs=1e-9)


def test_release_20261031_hyperfine_not_predicted_and_edges():
    """Where a second-order denominator nearly vanishes the release gives a line's center only; the uncertainty
    model carries the X near-dissociation term, the rules beyond the tables and the tunneling width."""
    m = ResearchModel("research-20261031-i2spec")
    missing = {(r["v"], r["J"]) for r in _read(_entry("research-20261031-expanded"), "hyperfine.json")["B"]["not_predicted"]}
    assert (71, 20) in missing
    assert np.isfinite(m.frequency("R(19) 71-0"))
    with pytest.raises(ValueError, match="does not predict"):
        m.frequency("R(19) 71-0", "a1")
    table = m._uncertainty_model().table
    assert table["model_error"]["X beyond"]["near_dissociation"]["v0"] == 85
    assert table["edges"]["B above 80"]["growth_per_v"] == 3.0 and table["edges"]["X above 108"]["v_max"] == 113
    unc = m._uncertainty_model()
    core = m._spec["core"]
    s80, s82 = (unc.model_sigma(core, "B", v, 10, "beyond") for v in (80, 82))
    assert s82 == pytest.approx(9*s80)
    x84, x90 = (unc.model_sigma(core, "X", v, 10, "beyond") for v in (84, 90))
    assert x90 > x84
    assert m._levels_x.width("B", 30, 50) == 0.
    assert m._levels_x.width("B", 61, 125)*MHZ_PER_CM > 10.


def test_release_20261031_core_patterns_unchanged():
    """In the core the release changes no hyperfine pattern, at any J (above J = 200 a J +- 2 partner beyond the
    barrier is left out, as before)."""
    old, new = ResearchModel("research-20261030-i2spec"), ResearchModel("research-20261031-i2spec")
    for vu, vl, J, branch in ((43, 0, 208, "R"), (32, 0, 56, "R"), (20, 3, 210, "P")):
        a = [c.offset for c in old.hyperfine_components(vu, vl, J, branch)[1]]
        b = [c.offset for c in new.hyperfine_components(vu, vl, J, branch)[1]]
        assert b == pytest.approx(a, abs=1e-6)
