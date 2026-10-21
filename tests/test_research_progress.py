"""Live accounting must include helper results without mistaking failures for acceptance."""
from datetime import datetime, timedelta, timezone
import json
import os
from pathlib import Path

import pytest

from i2spec import research_progress as progress


NOW = datetime(2026, 10, 11, 20, tzinfo=timezone.utc)
POTENTIAL, HYPERFINE = (item[0] for item in progress.CAMPAIGNS)
HELPER = progress.CAMPAIGNS[1][2]


def write(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value), encoding="utf-8")


def job(seed, cell="candidate", fold="physical_level:B:39"):
    return dict(phase="validate", cell=dict(id=cell), fold=fold, seed=seed, suffix=f"validation_s{seed}")


def completion(j, accepted=True, hours_ago=1, prefix=HYPERFINE):
    return dict(j, tag=progress._tag(j, prefix), accepted=accepted,
                started_utc=(NOW-timedelta(hours=30)).isoformat(),
                finished_utc=(NOW-timedelta(hours=hours_ago)).isoformat())


@pytest.fixture
def out(tmp_path):
    for prefix, workers in ((POTENTIAL, 6), (HYPERFINE, 2), (HELPER, 4)):
        write(tmp_path/prefix/"state.json", dict(status="running", workers=workers, completed=0))
    write(tmp_path/HELPER/"launch.json", dict(started_utc=(NOW-timedelta(hours=18)).isoformat()))
    write(tmp_path/POTENTIAL/"jobs_validate.json", [job(0)])
    write(tmp_path/HYPERFINE/"jobs_validate.json", [job(0), job(2), job(4), job(6)])
    return tmp_path


def test_canonical_helper_results_count_once_despite_lagging_state(out):
    j = job(0)
    row = completion(j, prefix=HELPER)
    # The canonical filename is original; the honest record tag names the helper.
    write(out/HYPERFINE/"attempts"/(progress._tag(j, HYPERFINE)+".json"), row)
    write(out/HELPER/"attempts"/(row["tag"]+".json"), row)
    j2 = job(2)
    write(out/HYPERFINE/"attempts"/(progress._tag(j2, HYPERFINE)+".json"), completion(j2, accepted=False))
    # A result outside the frozen queue must not inflate its progress.
    orphan = job(99)
    write(out/HYPERFINE/"attempts"/(progress._tag(orphan, HYPERFINE)+".json"), completion(orphan))
    data = progress.ProgressMonitor(out).snapshot(now=NOW)
    c = data["campaigns"][1]
    assert (c["finished"], c["accepted"], c["failed"], c["remaining"]) == (2, 1, 1, 2)
    assert c["progress_pct"] == 50
    assert data["totals"]["configured_workers"] == 12
    assert len(data["recent_completions"]) == 2


def test_active_helper_fit_and_partial_log_line(out):
    j = job(4);tag = progress._tag(j, HELPER)
    write(out/HELPER/"launches"/(tag+".json"), dict(j, tag=tag, started_utc=(NOW-timedelta(hours=2)).isoformat()))
    path = out/(tag+"_potential.log")
    path.write_text('{"evaluation": 10, "jacobian_evaluations": 3, "cost": 5.0}\n{"evaluation":')
    stamp = (NOW-timedelta(seconds=30)).timestamp();os.utime(path, (stamp, stamp))
    c = progress.ProgressMonitor(out).snapshot(now=NOW)["campaigns"][1]
    assert (c["active"], c["queued"], c["recently_active"]) == (1, 3, 1)
    active = c["active_jobs"][0]
    assert active["stage"] == "Fitting potential"
    assert active["evaluation"]["evaluation"] == 10
    assert active["evaluation"]["jacobian_evaluations"] == 3
    assert active["elapsed_hours"] == 2


def test_completed_launch_is_not_an_active_worker(out):
    j = job(0)
    for prefix in (HYPERFINE, HELPER):
        write(out/prefix/"launches"/(progress._tag(j, prefix)+".json"), dict(j, tag=progress._tag(j, prefix)))
    write(out/HYPERFINE/"attempts"/(progress._tag(j, HYPERFINE)+".json"), completion(j))
    c = progress.ProgressMonitor(out).snapshot(now=NOW)["campaigns"][1]
    assert c["active"] == 0
    assert c["remaining"] == 3


def test_overlap_counts_one_active_job_with_original_precedence(out):
    j = job(0)
    for prefix in (HYPERFINE, HELPER):
        write(out/prefix/"launches"/(progress._tag(j, prefix)+".json"), dict(j, tag=progress._tag(j, prefix)))
    c = progress.ProgressMonitor(out).snapshot(now=NOW)["campaigns"][1]
    assert c["active"] == 1
    assert c["active_jobs"][0]["tag"] == progress._tag(j, HYPERFINE)
    assert c["recently_active"] == 0  # No evidence of live log activity.


def test_rate_windows_exclude_completions_before_worker_change(out):
    for seed, hours_ago in ((0, 20), (2, 10)):
        j = job(seed)
        write(out/HYPERFINE/"attempts"/(progress._tag(j, HYPERFINE)+".json"), completion(j, hours_ago=hours_ago))
    c = progress.ProgressMonitor(out).snapshot(now=NOW)["campaigns"][1]
    rate = c["rates"][-1]
    assert rate["observed_hours"] == 18
    assert rate["count"] == 1
    assert rate["per_day"] == 1.33
    assert c["eta"]["days_low"] <= c["eta"]["days_high"]


def test_new_or_stalled_queue_has_no_invented_eta(out):
    monitor = progress.ProgressMonitor(out)
    assert monitor.snapshot(now=NOW)["campaigns"][1]["eta"] is None
    j = job(0)
    write(out/HYPERFINE/"attempts"/(progress._tag(j, HYPERFINE)+".json"), completion(j, hours_ago=48))
    assert monitor.snapshot(now=NOW)["campaigns"][1]["eta"] is None


def test_missing_queue_is_unavailable_not_finished(tmp_path):
    data = progress.ProgressMonitor(tmp_path).snapshot(now=NOW)
    assert data["available"] is False
    assert data["totals"]["progress_pct"] is None
    assert data["totals"]["eta"] is None
    assert not any(c["available"] for c in data["campaigns"])


def test_partial_queue_data_has_no_overall_finish_estimate(out):
    j = job(0)
    write(out/POTENTIAL/"attempts"/(progress._tag(j, POTENTIAL)+".json"), completion(j, prefix=POTENTIAL))
    (out/HYPERFINE/"jobs_validate.json").unlink()
    data = progress.ProgressMonitor(out).snapshot(now=NOW)
    assert data["available"] is False
    assert data["totals"]["eta"] is None


def test_corrupt_or_unfinished_results_do_not_count_as_completed(out):
    path = out/HYPERFINE/"attempts"/(progress._tag(job(0), HYPERFINE)+".json")
    path.parent.mkdir(parents=True);path.write_text('{"accepted":')
    write(out/HYPERFINE/"attempts"/(progress._tag(job(2), HYPERFINE)+".json"), dict(job(2), accepted=True))
    data = progress.ProgressMonitor(out).snapshot(now=NOW)
    assert data["campaigns"][1]["finished"] == 0
    assert data["warnings"]


def test_read_cache_invalidates_and_monitor_does_not_modify_artifacts(out):
    j = job(0);path=out/HYPERFINE/"attempts"/(progress._tag(j, HYPERFINE)+".json")
    write(path, completion(j))
    before = {p: (p.read_bytes(), p.stat().st_mtime_ns) for p in out.rglob('*.json')}
    monitor = progress.ProgressMonitor(out)
    assert monitor.snapshot(now=NOW)["campaigns"][1]["accepted"] == 1
    assert all((p.read_bytes(), p.stat().st_mtime_ns) == value for p, value in before.items())
    write(path, completion(j, accepted=False))
    data = monitor.snapshot(now=NOW)
    assert data["campaigns"][1]["accepted"] == 0
    assert data["campaigns"][1]["failed"] == 1
    assert data["daily_completions"][-1]["failed"] == 1


def test_leader_validation_uses_only_frozen_jobs(out):
    last = dict(free=61, inherited=11, score=dict(unique_absolute_comb=dict(rms_MHz=.156)),
                hyperfine_held_out=dict(rms_MHz=.108), potential_training_cost=100)
    write(out/HYPERFINE/"comparison.json", dict(selected=[dict(cell=dict(id="candidate"), last=last)]))
    for seed in (0, 99):
        j=job(seed)
        write(out/HYPERFINE/"attempts"/(progress._tag(j, HYPERFINE)+".json"), dict(completion(j), last=last))
    leader = progress.ProgressMonitor(out).snapshot(now=NOW)["leader"]
    assert leader["free"] == 61
    assert leader["calibration"] == 4
    assert leader["validation"]["accepted"] == 1
    assert len(leader["validation"]["folds"]) == 1


def test_finalization_is_separate_and_reads_live_worker_logs(out):
    prefix = "followup"
    root = out/prefix
    j = dict(job(2), phase="study")
    tag = progress._tag(j, prefix)
    write(out/"finalization_registry.json", dict(worker_budget=12, tracks=[
        dict(id="stability", label="Starting values", root=prefix, workers=2, detail="Training-only study")]))
    write(root/"jobs_study.json", [dict(job(0), phase="study"), j])
    write(root/"state.json", dict(status="running", phase="study", active=1, workers=2, active_jobs=[tag]))
    write(root/"jobs"/(tag+".json"), j)
    write(root/"worker_runs"/(tag+".json"), dict(started_utc=(NOW-timedelta(hours=1)).isoformat()))
    write(root/"attempts"/"failed.json", dict(phase="study", accepted=False, finished_utc=NOW.isoformat()))
    path = root/"logs"/(tag+"_potential.log")
    path.parent.mkdir(parents=True)
    path.write_text('{"evaluation": 7, "cost": 10}\n')
    os.utime(path, (NOW.timestamp(), NOW.timestamp()))
    data = progress.ProgressMonitor(out).snapshot(now=NOW)
    assert data["totals"]["finished"] == 0  # Original validation percentage stays historical.
    f = data["finalization"]
    assert (f["finished"], f["total"], f["active"], f["recently_active"]) == (1, 2, 1, 1)
    assert f["tracks"][0]["failed"] == 1
    assert f["tracks"][0]["active_jobs"][0]["evaluation"]["evaluation"] == 7
    assert data["follow_up"][1]["status"] == "running"


def test_completed_joint_pilot_is_not_an_accepted_final_model(out):
    root = out/"joint"
    write(out/"finalization_registry.json", dict(tracks=[
        dict(id="joint", root="joint", label="Joint pilots", kind="pilot")]))
    write(root/"state.json", dict(status="running", total=2, finished=1, active=0, workers=2,
        jobs=[dict(status="finished"), dict(status="queued")]))
    track = progress.ProgressMonitor(out).snapshot(now=NOW)["finalization"]["tracks"][0]
    assert (track["finished"], track["accepted"], track["failed"]) == (1, 0, 0)


def test_finalization_paths_cannot_escape_the_artifact_directory(out):
    write(out/"finalization_registry.json", dict(tracks=[dict(root="../", id="bad", label="Bad")]))
    data = progress.ProgressMonitor(out).snapshot(now=NOW)
    assert data["finalization"]["tracks"] == []
    assert any("outside" in warning for warning in data["warnings"])


def test_joint_recovery_gate_keeps_followup_track_open(out):
    root = out/"joint"
    write(out/"finalization_registry.json", dict(tracks=[
        dict(id="joint", root="joint", label="Joint pilots", kind="pilot")]))
    write(root/"state.json", dict(status="finished", total=6, finished=6, active=0, workers=0))
    write(root/"review_gate.json", dict(ready_for_review=False, additional_track_states=["recovery.json"]))
    write(root/"recovery.json", dict(status="running", total=7, finished=1, active=2, workers=2))
    finalization = progress.ProgressMonitor(out).snapshot(now=NOW)["finalization"]
    assert finalization["completed_tracks"] == 0
    assert (finalization["total"], finalization["finished"], finalization["active"]) == (7, 1, 2)


def test_supplemental_workers_and_terminal_aliases_are_counted_once(out):
    prefix="p13_finalization_stability_v1"
    root=out/prefix
    j0,j1=job(0),job(2)
    tag0,tag1=(progress._tag(j,prefix) for j in (j0,j1))
    legacy=tag0.replace(prefix,"p13_model_exploration_v1",1)
    write(out/"finalization_registry.json",dict(tracks=[dict(id="stability",label="Starts",root=prefix,
        supplemental_states=["supplemental_state.json"])]))
    write(root/"manifest.json",dict(jobs=[j0,j1,job(4)]))
    write(root/"state.json",dict(status="running",workers=2,active_jobs=[legacy]))
    write(root/"supplemental_state.json",dict(status="running",workers=6,active_jobs=[tag0,tag1]))
    write(root/"jobs"/(legacy+".json"),j0)
    write(root/"jobs"/(tag0+".json"),j0)
    write(root/"jobs"/(tag1+".json"),j1)
    terminal=dict(tag="ended",accepted=True,finished_utc=NOW.isoformat())
    write(root/"attempts"/"canonical.json",terminal)
    write(root/"attempts"/"alias.json",terminal)
    track=progress.ProgressMonitor(out).snapshot(now=NOW)["finalization"]["tracks"][0]
    assert (track["total"],track["finished"],track["active"],track["workers"])==(3,1,2,8)


def test_reproduction_success_does_not_imply_numerical_finalization(out):
    write(out/"finalization_registry.json",dict(tracks=[],review_root="review"))
    write(out/"review"/"review.json",dict(status="ready_for_review",model_finalization_ready=False))
    write(out/"review"/"verification_state.json",dict(status="finished",passed=True))
    review=progress.ProgressMonitor(out).snapshot(now=NOW)["finalization"]["review"]
    assert review["snapshot_verified"] and not review["numerical_review_ready"]


def test_exact_leader_waits_for_validation_and_snapshot_replay(out):
    root=out/"exact"
    write(out/"finalization_registry.json",dict(tracks=[dict(id="atlas_exact",label="Exact",root="exact")],review_root="review"))
    folds=["physical_level:B:1","physical_level:B:16","physical_level:B:33","physical_level:B:39",
           "physical_level:X:4","physical_source:kobayashi2016a","physical_source:sakagami2020a"]
    jobs=[dict(cell=dict(id="exact_model"),phase="validate",fold=fold,seed=seed,suffix=f"s{seed}")
          for fold in folds for seed in (0,2)]
    write(root/"jobs_validate.json",jobs)
    write(root/"state.json",dict(status="finished",phase="validate",workers=0))
    for j in jobs:
        tag=progress._tag(j,"exact")
        write(root/"attempts"/(tag+".json"),dict(j,tag=tag,accepted=True,finished_utc=NOW.isoformat(),
            last=dict(score=dict(unique_absolute_comb=dict(rms_MHz=.2)))))
    candidate=dict(model="exact_model",free=61,inherited=11,eligible=True,retained_atlas_exact_J=True,
        absolute_rms_MHz=.15,intraline_rms_MHz=.06,fold_results={fold:{} for fold in folds})
    write(out/"review"/"review.json",dict(status="ready_for_review",model_finalization_ready=True,
        recommendations=dict(absolute_accuracy=candidate)))
    write(out/"review"/"verification_state.json",dict(status="running",passed=False))
    assert progress.ProgressMonitor(out).snapshot(now=NOW)["leader"] is None
    write(out/"review"/"verification_state.json",dict(status="finished",passed=True))
    leader=progress.ProgressMonitor(out).snapshot(now=NOW)["leader"]
    assert leader["model"]=="exact_model" and leader["validation"]["accepted"]==14
    assert len(leader["validation"]["folds"])==14


def completion_gates():
    return dict(schema_version=1, status="computational_release_ready_awaiting_new_measurements",
        computational_release_ready=True, scientific_validation_complete=False, gates=[
            dict(id="numerical", label="Numerical checks and replay", status="complete", detail="Frozen predictions reproduced"),
            dict(id="blind", label="Fresh experimental validation", status="awaiting_measurements", detail="Acquire new frequencies"),
            dict(id="primary", label="Source documentation", status="requires_primary_source", detail="Obtain original source"),
            dict(id="coverage", label="Probability coverage", status="not_established", detail="Do not claim 95% coverage")])


def test_finished_queues_do_not_complete_fresh_experimental_validation(out):
    for prefix, seeds in ((POTENTIAL, (0,)), (HYPERFINE, (0, 2, 4, 6))):
        for seed in seeds:
            row = job(seed)
            write(out/prefix/"attempts"/(progress._tag(row, prefix)+".json"), completion(row, prefix=prefix))
    path = out/"completion_state.json"
    payload = completion_gates()
    write(path, payload)
    before = path.read_bytes(), path.stat().st_mtime_ns
    monitor = progress.ProgressMonitor(out)
    data = monitor.snapshot(now=NOW)
    assert data["totals"]["progress_pct"] == 100
    assert data["completion"]["computational_release_ready"]
    assert not data["completion"]["scientific_validation_complete"]
    assert data["follow_up"] == [dict(title=g["label"], status=g["status"], detail=g["detail"]) for g in payload["gates"]]
    assert (path.read_bytes(), path.stat().st_mtime_ns) == before
    payload["gates"][1]["detail"] = "Changed acquisition plan"
    write(path, payload)
    assert monitor.snapshot(now=NOW)["completion"]["gates"][1]["detail"] == "Changed acquisition plan"


def test_missing_or_invalid_completion_state_cannot_invent_completion(out):
    monitor = progress.ProgressMonitor(out)
    assert monitor.snapshot(now=NOW)["completion"] is None
    payload = completion_gates()
    payload["gates"][1]["status"] = "unknown"
    write(out/"completion_state.json", payload)
    data = monitor.snapshot(now=NOW)
    assert data["completion"] is None
    assert any("completion_state.json" in warning for warning in data["warnings"])


def test_unresolved_gates_override_a_contradictory_scientific_completion_flag(out):
    payload = completion_gates()
    payload["scientific_validation_complete"] = True
    write(out/"completion_state.json", payload)
    data = progress.ProgressMonitor(out).snapshot(now=NOW)
    assert not data["completion"]["scientific_validation_complete"]
    assert any("unresolved gates" in warning for warning in data["warnings"])


