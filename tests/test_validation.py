"""Forecast immutability, measurement identity and prospective-test accounting."""
from datetime import datetime, timedelta, timezone
import json

import pytest

from i2spec.validation import freeze_registry, main, read_registry, score_measurements


def factory(name):
    offset = {"primary": 0., "other": .2}[name]
    def frequency(line, component):
        return 500000000. + offset + int(component[1:])*.1
    return name, {"id": name}, {"offset": offset}, frequency


def plan():
    return dict(models=["primary", "other"], primary_model="primary", outcome_phase="prospective_laboratory",
        thresholds_MHz=[.1, .3, 1.],
        targets=[dict(id="p53-a1", line="127I2 P(53) 32-0", component="a1", stratum="routine"),
                 dict(id="p53-a3", line="127I2 P(53) 32-0", component="a3", stratum="routine")],
        envelopes={"routine": dict(half_width_MHz=.6, nominal_probability=None, scope="Reused benchmark; descriptive")})


@pytest.fixture
def registry(tmp_path):
    folder = tmp_path/"frozen"
    freeze_registry(folder, plan(), factory=factory)
    return folder


def measurement(target="p53-a1", *, identity="new-1", value=500000000.2):
    return dict(id=identity, target_id=target, value_MHz=value, comparable=True,
        source_id="lab-run-1", source_url="urn:lab:run-1", conditions={"pressure_Pa": 1., "corrections": "source supported"},
        frequency_basis="zero_pressure_unperturbed", standard_uncertainty_MHz=.001,
        independent_of_fitting_archive=True, outcome_unseen_before_freeze=True,
        acquired_utc=datetime.now(timezone.utc).isoformat())


def test_frozen_registry_refuses_overwrite_and_detects_changed_bytes(registry):
    protocol, rows, digest = read_registry(registry)
    assert protocol["primary_model"] == "primary" and len(rows) == 4 and len(digest) == 64
    with pytest.raises(FileExistsError):
        freeze_registry(registry, plan(), factory=factory)
    with (registry/"forecasts.json").open("a") as stream:
        stream.write(" ")
    with pytest.raises(ValueError, match="hash mismatch"):
        read_registry(registry)


@pytest.mark.parametrize("change", ["excluded", "duplicate", "probability", "same_alias"])
def test_invalid_design_is_rejected_before_freezing(tmp_path, change):
    protocol = plan()
    predictor = factory
    if change == "excluded":
        protocol["targets"][0]["line"] = "127I2 P(13) 43-0"
    elif change == "duplicate":
        protocol["targets"][1]["component"] = "a1"
    elif change == "probability":
        protocol["envelopes"]["routine"]["nominal_probability"] = .95
    else:
        predictor = lambda name: factory("primary")
    with pytest.raises(ValueError):
        freeze_registry(tmp_path/"bad", protocol, factory=predictor)
    assert not (tmp_path/"bad").exists()


def test_repeats_do_not_create_independent_targets_and_intraline_errors_are_scored(registry):
    rows = [measurement(), measurement(identity="repeat", value=500000000.4),
            measurement("p53-a3", identity="third", value=500000000.6)]
    report = score_measurements(registry, rows)
    primary = report["summaries"]["primary"]
    assert primary["n"] == 2 and primary["n_measurements"] == 3
    assert primary["n_lines"] == primary["n_sources"] == 1
    assert primary["n_line_source_connected_blocks"] == 1
    assert primary["per_stratum"]["routine"]["n"] == 2
    assert primary["bias_MHz"] == pytest.approx(-.25, abs=1e-7)
    assert primary["maximum_repeat_error_MHz"] == pytest.approx(.3, abs=1e-7)
    assert primary["derived_intraline_intervals"]["metrics"]["rms_MHz"] == pytest.approx(.1, abs=1e-7)
    assert report["prospective_measurement_plan_complete"] is True
    # The other model scores better but the primary selection is immutable.
    assert report["summaries"]["other"]["rms_MHz"] < primary["rms_MHz"]
    assert report["primary_model"] == "primary" and report["refitted"] is False
    assert "Requires scientific review" in report["model_adequacy"]


@pytest.mark.parametrize("identity", [{"line": "127I2 R(51) 32-0"}, {"component": "a3"}])
def test_measurement_line_and_component_must_match_preregistered_target(registry, identity):
    row = measurement()
    row.update(identity)
    with pytest.raises(ValueError, match="registered identity"):
        score_measurements(registry, [row])


def test_component_letter_alias_cannot_create_a_new_identity(registry):
    row = measurement()
    row["component"] = "b1"
    assert score_measurements(registry, [row])["scored_unique_targets"] == 1


@pytest.mark.parametrize("change", ["archive", "exposed", "old", "future", "naive", "conditions", "basis"])
def test_prospective_results_need_new_measurements_and_comparable_conditions(registry, change):
    row = measurement()
    if change == "archive":
        row["independent_of_fitting_archive"] = False
    elif change == "exposed":
        row["outcome_unseen_before_freeze"] = False
    elif change == "old":
        row["acquired_utc"] = "2020-01-01T00:00:00+00:00"
    elif change == "future":
        row["acquired_utc"] = (datetime.now(timezone.utc)+timedelta(days=1)).isoformat()
    elif change == "naive":
        row["acquired_utc"] = "2026-10-18T10:00:00"
    elif change == "conditions":
        row["conditions"] = {}
    else:
        row["frequency_basis"] = "unknown"
    with pytest.raises(ValueError):
        score_measurements(registry, [row])


def test_missing_or_noncomparable_rows_do_not_silently_become_validation(registry):
    pending = dict(id="not-yet", target_id="p53-a1", value_MHz=None)
    unknown = measurement("p53-a3", identity="unknown-convention")
    unknown.update(comparable=False, frequency_basis="unknown")
    report = score_measurements(registry, [pending, unknown])
    assert report["scored_unique_targets"] == 0
    assert len(report["pending_measurements"]) == len(report["noncomparable_measurements"]) == 1
    assert report["prospective_measurement_plan_complete"] is False
    with pytest.raises(ValueError, match="unique ids"):
        score_measurements(registry, [pending, pending])


def test_expanded_uncertainty_requires_explicit_coverage_factor(registry):
    row = measurement()
    row.update(standard_uncertainty_MHz=None, expanded_uncertainty_MHz=.004)
    assert score_measurements(registry, [row])["rows"][0]["measurement_standard_uncertainty_MHz"] is None
    row["coverage_factor"] = 2
    assert score_measurements(registry, [row])["rows"][0]["measurement_standard_uncertainty_MHz"] == .002
    primary = next(r for r in score_measurements(registry, [row])["rows"] if r["primary"])
    assert primary["residual_over_measurement_sigma"] == pytest.approx(-50, abs=1e-4)


def test_pending_score_cli_creates_separate_output_once(registry, tmp_path, capsys):
    measurements = tmp_path/"measurements.json"
    measurements.write_text(json.dumps([dict(id="pending", target_id="p53-a1", value_MHz=None)]))
    output = tmp_path/"report.json"
    digest = read_registry(registry)[2]
    assert main(["score", str(registry), str(measurements), "--out", str(output)]) == 0
    assert json.loads(capsys.readouterr().out)["scored_unique_targets"] == 0
    assert read_registry(registry)[2] == digest
    with pytest.raises(SystemExit):
        main(["score", str(registry), str(measurements), "--out", str(output)])
    assert "File exists" in capsys.readouterr().err


def test_installed_prospective_registry_is_complete_and_has_no_measurements():
    from importlib.resources import files
    folder = files("i2spec").joinpath("data/research/20261018/validation/prospective_v1")
    protocol, forecasts, _ = read_registry(folder)
    assert protocol["primary_model"] == "research-20261018-balanced"
    assert len(protocol["targets"]) == 48 and len(forecasts) == 192
    assert len({t["line"] for t in protocol["targets"]}) == 16
    assert len({t["pair_id"] for t in protocol["targets"]}) == 8
    assert protocol["outcome_phase"] == "prospective_laboratory"
    assert all(e["nominal_probability"] is None for e in protocol["envelopes"].values())
    score = score_measurements(folder, [])
    assert score["scored_unique_targets"] == 0 and len(score["missing_target_ids"]) == 48
    assert score["prospective_measurement_plan_complete"] is False
