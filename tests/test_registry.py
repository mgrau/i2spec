"""Blind-test registries (i2spec.fit.registry).

The two frozen registries of release 20261020 are intact; scoring them again from their frozen forecasts and
today's data reproduces their results; freezing the second one's sets again reproduces its targets and, to
numerical noise, its forecasts. Nothing is written into a frozen registry: every new file goes to tmp_path, and
no pending set outside the two registries is read.
"""
import hashlib
import json
import math
import shutil
from pathlib import Path

import pytest

from i2spec.fit import registry
from i2spec.observations import DATA_DIR, Line

VALIDATION = registry.PACKAGE / "data" / "research" / "20261020" / "validation"
V1, V2 = VALIDATION / "pending_blind_v1", VALIDATION / "pending_blind_v2"
#: v1 sets whose meta.toml gained an [overlaps] declaration after freezing; every v1 and v2 meta.toml was since
#: reworded too (source and exclude lines), so the input-hash check flags all of them
OVERLAPS_ADDED = ("jungner1995a", "picard2003a", "zhang1999a")
#: The registries were written on one platform; the arithmetic of another differs in the last digits. Forecasts of an
#: unchanged model differ by a few ulp of the line center (up to 1.1e-6 MHz); the 1-sigma by up to 1.4e-6 of itself,
#: well under the smallest intended change of v2's band intervals (1.3e-5); rescored results by a few ulp.
PREDICTION_TOL_MHZ = 1e-5
SIGMA_RTOL = 5e-6
SCORE_RTOL = 1e-9


def _load(path):
    return json.loads(Path(path).read_text())


def _differences(new, old, path="", rel=0.0):
    """Every leaf where two JSON values differ (floats by more than rel of themselves): (path, new, old)."""
    if isinstance(new, dict) and isinstance(old, dict):
        return [d for k in sorted(set(new) | set(old), key=str)
                for d in (_differences(new[k], old[k], f"{path}/{k}", rel) if k in new and k in old else
                          [(f"{path}/{k}", new.get(k, "<missing>"), old.get(k, "<missing>"))])]
    if isinstance(new, list) and isinstance(old, list):
        head = [(path, f"{len(new)} items", f"{len(old)} items")] if len(new) != len(old) else []
        return head + [d for i, (a, b) in enumerate(zip(new, old)) for d in _differences(a, b, f"{path}[{i}]", rel)]
    if isinstance(new, float) and isinstance(old, float):
        return [] if math.isclose(new, old, rel_tol=rel) else [(path, new, old)]
    return [] if new == old else [(path, new, old)]


@pytest.fixture(scope="module", autouse=True)
def frozen_registries_untouched():
    """Every file of the frozen registries, byte for byte, is the same after this module's tests."""
    def files():
        return {p: p.read_bytes() for d in (V1, V2) for p in sorted(d.iterdir())}
    before = files()
    yield
    assert files() == before


def test_frozen_registries_are_intact():
    for folder, n in ((V1, 580), (V2, 698)):
        protocol, targets, forecasts, snapshot, _ = registry.read_registry(folder)
        assert len(targets) == n and len(forecasts) == len(protocol["models"])*n
        assert {t["set"] for t in targets} == set(snapshot["data"])
    # the second registry holds only sets the first does not, and records the first as earlier
    assert _load(V2 / "snapshot.json")["earlier_sets"] == registry.registered_sets(V1)
    assert not set(registry.registered_sets(V1)) & set(registry.registered_sets(V2))


def test_a_changed_registry_file_is_refused(tmp_path):
    copy = tmp_path / "v2"
    shutil.copytree(V2, copy)
    forecasts = _load(copy / "forecasts.json")
    forecasts[0]["predicted_MHz"] += 1e-3
    (copy / "forecasts.json").write_text(json.dumps(forecasts, indent=1))
    with pytest.raises(ValueError, match="forecasts.json changed after freezing"):
        registry.read_registry(copy)
    with pytest.raises(ValueError, match="changed after freezing"):
        registry.score_registry(copy, tmp_path / "results.json")
    assert not (tmp_path / "results.json").exists()


def test_a_scored_registry_is_not_scored_again():
    for folder in (V1, V2):
        with pytest.raises(FileExistsError, match="scored once"):
            registry.score_registry(folder)


def test_freezing_refuses_existing_folders_and_sets_that_are_not_blind(tmp_path):
    with pytest.raises(FileExistsError):
        registry.freeze_registry(tmp_path, ["acef1993a"], after=[V1])
    active = next(p.name for p in sorted(DATA_DIR.iterdir())
                  if (p / "meta.toml").exists() and "exclude" not in registry._meta(DATA_DIR, p.name))
    for sets in ([active], ["baird1979a"]):         # an active set; a set the earlier registry holds
        with pytest.raises(ValueError, match="only pending sets"):
            registry.freeze_registry(tmp_path / "new", sets, after=[V1])
    assert not (tmp_path / "new").exists()


def _check_scores(new, folder):
    old = _load(folder / "results.json")
    diffs = [d for key in old if key not in ("scored_utc", "meta_changed_since_freeze")
             for d in _differences(new[key], old[key], key, SCORE_RTOL)]
    assert not diffs, diffs[:20]


def test_rescoring_v1_reproduces_its_results_and_reports_the_changed_metas(tmp_path):
    snapshot = _load(V1 / "snapshot.json")
    # v1 hashed each set's meta.toml with its data.csv: the strict check stops at the changed metas ...
    with pytest.raises(ValueError, match="inputs changed after freezing") as error:
        registry.score_registry(V1, tmp_path / "strict.json")
    assert all(f"{s}/meta.toml" in str(error.value) for s in OVERLAPS_ADDED)
    assert not (tmp_path / "strict.json").exists()
    # ... and without it they are listed: every meta.toml, no data.csv
    new = registry.score_registry(V1, tmp_path / "v1.json", strict=False)
    assert new["inputs_changed_since_freeze"] == sorted(f"{s}/meta.toml" for s in snapshot["data"])
    assert _load(tmp_path / "v1.json")["inputs_changed_since_freeze"] == new["inputs_changed_since_freeze"]
    _check_scores(_load(tmp_path / "v1.json"), V1)


def test_rescoring_v2_reproduces_its_results(tmp_path):
    snapshot = _load(V2 / "snapshot.json")
    new = registry.score_registry(V2, tmp_path / "v2.json")       # v2 hashed data.csv only: strict passes
    assert new["inputs_changed_since_freeze"] == []
    # the metas, hashed for information, were all reworded after scoring (the frozen results list none)
    assert _load(V2 / "results.json")["meta_changed_since_freeze"] == []
    assert new["meta_changed_since_freeze"] == sorted(snapshot["meta_informational"])
    _check_scores(_load(tmp_path / "v2.json"), V2)


@pytest.fixture(scope="module")
def refrozen_v2(tmp_path_factory):
    """The v2 sets frozen again, after v1, with the model released then (release 20261021, primary) and the release
    v2 forecast with."""
    from i2spec import research
    folder = tmp_path_factory.mktemp("registry") / "v2"
    with pytest.MonkeyPatch.context() as patch:
        patch.setattr(research, "_RELEASED", "20261021")
        registry.freeze_registry(folder, list(_load(V2 / "snapshot.json")["data"]), after=[V1],
                                 models=("i2spec", "balanced_expanded"))
    return registry.read_registry(folder)


def test_refreezing_v2_reproduces_its_targets_and_snapshot(refrozen_v2):
    protocol, targets, _, snapshot, _ = refrozen_v2
    old = _load(V2 / "snapshot.json")
    assert targets == _load(V2 / "targets.json")        # identity columns, status and training flags alike
    assert snapshot["data"] == old["data"]
    for key in ("earlier_registries", "earlier_sets"):
        assert snapshot[key] == old[key]
    # the model v2 forecast with is the same frozen release: its coefficients hash as v2 recorded once the key
    # the coefficients gained after freezing (extension, None for this release) is left out, and its
    # uncertainty.json is the one v2 recorded
    frozen = _load(V2 / "protocol.json")["models"]["balanced_expanded"]
    from i2spec import ResearchModel
    coefficients = ResearchModel("balanced_expanded").coefficients
    assert coefficients.pop("extension") is None
    assert hashlib.sha256(registry._canonical(coefficients)).hexdigest() == frozen["coefficients"]
    assert protocol["models"]["balanced_expanded"]["id"] == frozen["id"]
    release = protocol["models"]["i2spec"]["files_sha256"]["20261020"]
    assert release["uncertainty.json"] == frozen["uncertainty_json_sha256"]


def test_refreezing_v2_reproduces_its_forecasts(refrozen_v2):
    _, targets, forecasts, _, _ = refrozen_v2
    by_id = {t["id"]: t for t in targets}
    new = {(f["model"], f["target"]): f for f in forecasts}
    old = {f["target"]: f for f in _load(V2 / "forecasts.json") if f["model"] == "balanced_expanded"}
    predicted = [k for k, f in old.items() if f["predicted_MHz"] is not None]
    assert len(predicted) == 502 and all(by_id[k]["isotopologue"] == "127I2" for k in predicted)
    # the released model is the same as v2's for 127I2 (J'' <= 150 here): bit-identical forecasts
    assert all(new[("i2spec", k)]["predicted_MHz"] == new[("balanced_expanded", k)]["predicted_MHz"] for k in predicted)
    # and both match v2's to numerical noise (intervals within a line to 1 mHz)
    for k in predicted:
        assert new[("i2spec", k)]["predicted_MHz"] == pytest.approx(old[k]["predicted_MHz"], abs=PREDICTION_TOL_MHZ)
        if by_id[k]["kind"] == "intraline":
            assert new[("i2spec", k)]["predicted_MHz"] == pytest.approx(old[k]["predicted_MHz"], abs=1e-9)
    # the targets v2 had no forecast of (the other isotopologues) the released model now forecasts
    others = sorted(set(old) - set(predicted))
    assert len(others) == 196 and {by_id[k]["isotopologue"] for k in others} == {"129I2", "127I129I"}
    assert all(new[("i2spec", k)]["predicted_MHz"] is not None for k in others)
    assert all(new[("balanced_expanded", k)]["predicted_MHz"] is None for k in others)


def test_refreezing_v2_reproduces_its_uncertainties_but_for_one_band_intervals(refrozen_v2):
    """The 1-sigma matches v2's to SIGMA_RTOL except for intervals between two lines of one band, whose hyperfine part
    takes the pattern-scale error as common to both lines since v2 was frozen (at most 0.2% of the total)."""
    _, targets, forecasts, _, _ = refrozen_v2
    new = {f["target"]: f for f in forecasts if f["model"] == "i2spec"}
    old = {f["target"]: f for f in _load(V2 / "forecasts.json") if f["model"] == "balanced_expanded"}
    changed = set()
    for t in targets:
        if old[t["id"]]["predicted_MHz"] is None:
            continue
        a, b = new[t["id"]], old[t["id"]]
        assert (a["level_status"], a["hyperfine_status"]) == (b["level_status"], b["hyperfine_status"])
        assert a["sigma_parts_MHz"][:2] == pytest.approx(b["sigma_parts_MHz"][:2], rel=SIGMA_RTOL)
        if a["sigma_model_MHz"] != pytest.approx(b["sigma_model_MHz"], rel=SIGMA_RTOL):
            changed.add(t["id"])
            assert a["sigma_model_MHz"] == pytest.approx(b["sigma_model_MHz"], rel=2e-3)
    def band(text):
        line = Line.parse(text)
        return line.isotopologue, line.v_upper, line.v_lower
    one_band = {t["id"] for t in targets if t["kind"] == "interline" and t["component"] and t["ref_component"]
                and t["isotopologue"] == "127I2" and band(t["line"]) == band(t["ref_line"])}
    assert changed == one_band and len(changed) == 7


def _pending_set(folder, name, files, meta=""):
    path = folder / name
    path.mkdir(parents=True)
    (path / "meta.toml").write_text(f'id = "{name}"\nexclude = "pending: test"\n{meta}')
    for file, text in files.items():
        (path / file).write_text(text)


def test_constants_are_registered_from_identity_columns_and_scored(tmp_path):
    """Hyperfine constants and X intervals: freezing reads no value (the value columns hold text here), and
    scoring reads them, per quantity, with the second sublevel's I from a printed label in the note."""
    data, empty = tmp_path / "data", tmp_path / "observations"
    empty.mkdir()
    hf = "line,eqQ,eqQ_unc,C,C_unc,d,d_unc,delta,delta_unc,fit_sd_kHz,note\n"
    _pending_set(data / "hyperfine_parameters", "hfa", {"data.csv": hf + "127I2 R(56) 32-0,x,x,x,x,x,x,x,x,,\n"})
    _pending_set(data / "hyperfine_parameters", "hfb", {"data.csv": hf + "127I2 R(56) 32-0,x,x,x,x,x,x,x,x,,\n"},
                 'kind = "B-X hyperfine difference constants"\n')
    xi = "v,J,eps,F_from,F_to,value_kHz,uncertainty_kHz,table,note\n"
    _pending_set(data / "x_hyperfine", "xa", {
        "constants.csv": "v,J,quantity,value,uncertainty,unit,table,note\n0,13,eqQ,x,x,MHz,,\n0,13,C1,x,x,kHz,,\n0,13,Bo,x,x,MHz,,\n",
        "data.csv": xi + "0,13,3,14,13,x,x,,\n0,13,3,14,13,x,x,,\"printed 3,14 -> 5,13\"\n"
                    "0,13,1,12,11,x,x,,\"printed 1,12 -> 3,11\"\n"})
    folder = tmp_path / "registry"
    registry.freeze_registry(folder, [], models=("research-20261029-i2spec", "i2spec"), root=empty, constants=True,
                             constants_root=data)
    _, targets, forecasts, snapshot, _ = registry.read_registry(folder)
    assert [(t["kind"], t["quantity"]) for t in targets] == [("hf_constant", q) for q in registry.HF_QUANTITIES]*2 + [
        ("x_constant", "eqQ"), ("x_constant", "C"), ("x_interval", None), ("x_interval", None), ("x_interval", None)]
    assert [t["convention"] for t in targets[:8:4]] == ["B", "B - X"]
    from i2spec import ResearchModel
    model = ResearchModel("research-20261029-i2spec")
    p = {f["target"]: f for f in forecasts if f["model"] == "research-20261029-i2spec"}
    B, X = model.level_hyperfine("B", 32, 57), model.level_hyperfine("X", 0, 56)
    assert p["hyperfine_parameters/hfa:0:C"]["predicted_MHz"] == pytest.approx(B.C*1e-3)
    assert p["hyperfine_parameters/hfb:0:eqQ"]["predicted_MHz"] == pytest.approx(B.eqQ-X.eqQ)
    # the values written after freezing: scoring reads them (strict=False: the data files changed)
    (data / "hyperfine_parameters" / "hfa" / "data.csv").write_text(
        hf + "127I2 R(56) 32-0,-544.7,0.01,0.2,0.0001,-0.1,0.001,0.005,0.001,,\n")
    (data / "hyperfine_parameters" / "hfb" / "data.csv").write_text(hf + "127I2 R(56) 32-0,1908,0.01,,,,,,,,\n")
    (data / "x_hyperfine" / "xa" / "constants.csv").write_text(
        "v,J,quantity,value,uncertainty,unit,table,note\n0,13,eqQ,-2452.6,0.01,MHz,,\n0,13,C1,3.1,0.1,kHz,,\n0,13,Bo,1,1,MHz,,\n")
    (data / "x_hyperfine" / "xa" / "data.csv").write_text(xi + "0,13,3,14,13,1000,1,,\n0,13,3,14,13,2000,1,,\"printed 3,14 -> 5,13\"\n"
                                                           "0,13,1,12,11,3000,1,,\"printed 1,12 -> 3,11\"\n")
    result = registry.score_registry(folder, tmp_path / "results.json", root=empty, strict=False)
    rows = {r["id"]: r for r in result["rows"]}
    assert rows["x_hyperfine/xa:constants.csv:1:C"]["observed_MHz"] == pytest.approx(3.1e-3)
    alt = {f["target"]: f["alternatives_MHz"] for f in forecasts if f["model"] == "research-20261029-i2spec"
           and "alternatives_MHz" in f}
    second = rows["x_hyperfine/xa:data.csv:1"]
    assert second["I_to_printed"] == 5
    assert second["residual_MHz"]["research-20261029-i2spec"] == pytest.approx(2.-alt["x_hyperfine/xa:data.csv:1"]["5"])
    # at J = 13, (1, 11) does not exist (I = 1 has F = 12-14): the same-I pair has no forecast, but the pair the
    # note prints, (1, 12) -> (3, 11), is forecast and scored
    third = {f["model"]: f for f in forecasts if f["target"] == "x_hyperfine/xa:data.csv:2"}
    assert third["research-20261029-i2spec"]["predicted_MHz"] is None
    assert rows["x_hyperfine/xa:data.csv:2"]["residual_MHz"]["research-20261029-i2spec"] == pytest.approx(
        3.-third["research-20261029-i2spec"]["alternatives_MHz"]["3"])
    assert sorted(result["not_printed"]) == [f"hyperfine_parameters/hfb:0:{q}" for q in ("C", "d", "delta")]
    assert "blind | hf_constant eqQ | all" in result["summary"] and "blind | x_interval | all" in result["summary"]
    assert "hyperfine_parameters/hfa | C" in result["by_set"]


def test_named_constant_sets_and_protocol_fields(tmp_path):
    """``constants`` may name the sets to register; ``protocol_extra`` is frozen into the hashed protocol."""
    data, empty = tmp_path / "data", tmp_path / "observations"
    empty.mkdir()
    hf = "line,eqQ,eqQ_unc,C,C_unc,d,d_unc,delta,delta_unc,fit_sd_kHz,note\n"
    for name in ("hfa", "hfb"):
        _pending_set(data / "hyperfine_parameters", name, {"data.csv": hf + "127I2 R(56) 32-0,x,x,x,x,x,x,x,x,,\n"})
    (data / "x_hyperfine").mkdir()
    folder = tmp_path / "registry"
    registry.freeze_registry(folder, [], models=("i2spec",), root=empty, constants=["hfb"], constants_root=data,
                             protocol_extra=dict(criteria="fixed"))
    protocol, targets, *_ = registry.read_registry(folder)
    assert {t["set"] for t in targets} == {"hfb"} and protocol["criteria"] == "fixed"
    with pytest.raises(ValueError, match="pending constant sets"):
        registry.freeze_registry(tmp_path / "other", [], models=("i2spec",), root=empty, constants=["hfc"],
                                 constants_root=data)
