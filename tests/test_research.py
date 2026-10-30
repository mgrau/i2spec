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
    assert [m.role for m in models] == ["balanced", "balanced_expanded", "i2spec"]*2
    assert [(m.inventory.fitted_molecular, m.inventory.inherited_molecular, m.inventory.fitted_calibration)
            for m in models] == [(61, 18, 4), (105, 73, 5), (107, 75, 5), (63, 18, 4), (106, 73, 5), (106, 75, 5)]
    # release 20261029 (the final form) is reached by id; each role resolves to the first release that has it
    assert [_entry(m.role)["id"] for m in models[:3]] == [m.id for m in models[:3]]
    assert [m.id for m in models[3:]] == [f"research-20261029-{n}" for n in ("core", "expanded", "i2spec")]
    assert DEFAULT_PARAMETERS == _entry("i2spec")["id"] == "research-20261021-i2spec"
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
