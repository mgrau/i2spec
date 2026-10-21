"""Read-only progress reporting for the local spectroscopy research queues.

Only saved artifacts are inspected. No optimizer imports, process signals,
coordinator writes, or model calculations are needed to update this dashboard.
"""
from __future__ import annotations

from datetime import datetime, timedelta, timezone
import argparse
import json
import math
from pathlib import Path
from threading import Lock
import time

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "prototypes/out/explore/physical_search"
CAMPAIGNS = (
    ("p13_model_exploration_v1", "Potential models", None),
    ("p13_model_exploration_hf_v1", "Hyperfine allocations", "p13_model_exploration_hf_parallel_v1"),
)


def _date(value):
    try:
        result = datetime.fromisoformat(value)
        return result.replace(tzinfo=timezone.utc) if result.tzinfo is None else result
    except (ValueError, TypeError):
        return None


def _finite(value):
    return value if isinstance(value, (int, float)) and math.isfinite(value) else None


def _tag(job, prefix):
    return "_".join((prefix, job["cell"]["id"], job["fold"].replace(":", "_"), job["suffix"]))


def _stage(name):
    if "fshape" in name:
        return "Refining physical shape"
    if any(word in name for word in ("audit", "mesh", "quadrature", "actual_j")):
        return "Checking numerical accuracy"
    if "hyperfine" in name or name.endswith(("_before.log", "_after.log")):
        return "Fitting hyperfine structure"
    if "potential" in name:
        return "Fitting potential"
    return "Fitting and checking model"


def _last_evaluation(path):
    """Read a bounded log tail, even while the optimizer is appending to it."""
    try:
        with path.open("rb") as stream:
            stream.seek(0, 2)
            stream.seek(max(0, stream.tell() - 16384))
            lines = stream.read().decode("utf-8", errors="replace").splitlines()
        for line in reversed(lines):
            try:
                row = json.loads(line)
            except ValueError:
                continue
            if isinstance(row, dict) and "evaluation" in row:
                return {k: _finite(row.get(k)) for k in
                        ("evaluation", "jacobian_evaluations", "cost", "physical_cost")}
    except OSError:
        pass
    return {}


class ProgressMonitor:
    """Cache unchanged JSON files and refresh a compact snapshot at most every 10 s."""

    def __init__(self, out=OUT):
        self.out = Path(out)
        self._files = {}
        self._snapshot = None
        self._refreshed = 0.
        self._lock = Lock()

    def _read(self, path, default=None):
        try:
            stat = path.stat()
            version = (stat.st_mtime_ns, stat.st_size, stat.st_ino)
            cached = self._files.get(path)
            if cached and cached[0] == version:
                return cached[1]
            result = json.loads(path.read_text(encoding="utf-8"))
            self._files[path] = (version, result)
            return result
        except FileNotFoundError:
            return default
        except (OSError, ValueError) as error:
            self.warnings.append(f"Could not read {path.name}: {type(error).__name__}")
            return default

    def snapshot(self, *, now=None, force=False):
        with self._lock:
            if not force and now is None and self._snapshot and time.monotonic() - self._refreshed < 10:
                return self._snapshot
            now = now or datetime.now(timezone.utc)
            self.warnings = []
            # One directory scan is shared by every active fit.
            logs = list(self.out.glob("*.log"))
            helper = self._read(self.out / CAMPAIGNS[1][2] / "launch.json", {})
            changed = _date(helper.get("started_utc"))
            control = self.out / CAMPAIGNS[1][2] / "control.json"
            if control.exists():
                modified = datetime.fromtimestamp(control.stat().st_mtime, timezone.utc)
                changed = max(changed, modified) if changed else modified
            campaigns = [self._campaign(prefix, label, extra, now, changed, logs)
                         for prefix, label, extra in CAMPAIGNS]
            totals = {key: sum(c[key] for c in campaigns) for key in
                      ("total", "finished", "accepted", "failed", "remaining", "active", "queued",
                       "configured_workers", "recently_active")}
            totals["progress_pct"] = 100 * totals["finished"] / totals["total"] if totals["total"] else None
            etas = [c["eta"] for c in campaigns if c["remaining"]]
            totals["eta"] = (dict(days_low=max(e["days_low"] for e in etas),
                                  days_high=max(e["days_high"] for e in etas))
                             if etas and all(e is not None for e in etas) else None)
            if totals["total"] and not totals["remaining"]:
                totals["eta"] = dict(days_low=0, days_high=0)
            if not all(c["available"] for c in campaigns):
                totals["eta"] = None
            completed = [r for c in campaigns for r in c.pop("_completed")]
            selected = [r for c in campaigns for r in c.pop("_selected")]
            leader = self._leader(selected, campaigns)
            daily = []
            for offset in range(6, -1, -1):
                day = (now - timedelta(days=offset)).date()
                rows = [r for r in completed if _date(r["finished_utc"]).date() == day]
                daily.append(dict(date=day.isoformat(), total=len(rows),
                                  accepted=sum(bool(r["accepted"]) for r in rows),
                                  failed=sum(not r["accepted"] for r in rows)))
            activity = sorted(completed, key=lambda r: r["finished_utc"], reverse=True)[:12]
            legacy = self._read(self.out / "p13_grid_validation_comparison.json", {})
            milestones = [dict(label="Potential architecture screening", finished=self._phase(campaigns[0], "screen")["finished"],
                               total=self._phase(campaigns[0], "screen")["total"]),
                          dict(label="Initial potential refits", **self._phase(campaigns[0], "complete")),
                          dict(label="Initial hyperfine allocation refits", **self._phase(campaigns[1], "complete")),
                          dict(label="Earlier 72 / 76 / 80 validation", finished=legacy.get("accepted_attempts", 0), total=42)]
            result = dict(generated_utc=now.isoformat(), refresh_seconds=30,
                available=all(c["available"] for c in campaigns), totals=totals, campaigns=campaigns,
                completion=self._completion(),
                finalization=self._finalization(now, logs),
                leader=leader, daily_completions=daily, recent_completions=activity, milestones=milestones,
                follow_up=[
                    dict(title="Compare final models", status="After validation",
                         detail="Compare fitted coefficient count, absolute errors, hyperfine intervals and each withheld band/source."),
                    dict(title="Investigate failures and unstable starts", status="Needs review",
                         detail="Review important rejected fits and sensitivity to initialization, including B(v=39). Recovery fits are outside these queue totals."),
                    dict(title="Joint fit to raw measurements", status="Not queued",
                         detail="Fit potential and hyperfine coefficients together, with measurement correlations and no double counting. Needs implementation and a timing pilot."),
                    dict(title="Extend promising combinations", status="Not queued",
                         detail="Test promising hyperfine layouts on other potential families, then assess regularization and broader vibrational/isotope coverage."),
                ], warnings=list(dict.fromkeys(self.warnings)),
                notes=["Progress measures the two frozen validation queues; accepted and failed attempts both count as finished.",
                       "Worker activity comes from saved launches and log timestamps, not an operating-system process check.",
                       "ETA uses recent completion rates at the current allocation; it excludes recovery fits and unqueued studies.",
                       "Accuracy uses a reused exploratory benchmark: 60 distinct absolute components and 592 hyperfine intervals. P(13) 43–0 remains excluded.",
                       "Finished fit queues and numerical replay do not establish fresh or blind experimental validation.",
                       "Fitted, inherited and calibration coefficients are counted separately. No independent line or level shifts are fitted."])
            if result["finalization"]:
                tracks = {track["id"]: track for track in result["finalization"]["tracks"]}
                result["follow_up"][0].update(status="Baseline comparison complete",
                    detail="The original campaigns are complete; finalization studies and their acceptance reviews are tracked above.")
                for index, key in ((1, "stability"), (2, "joint"), (3, "combinations")):
                    if key in tracks:
                        result["follow_up"][index].update(status=tracks[key]["status"], detail=tracks[key]["detail"])
                review = result["finalization"].get("review") or {}
                if review.get("leader"):
                    result["leader"] = review["leader"]
                    result["follow_up"][0].update(status="Validated research snapshots ready",
                        detail="The exact-J finalists passed numerical checks, seven-fold validation and frozen prediction replay.")
            if result["completion"]:
                result["follow_up"] = [dict(title=gate["label"], status=gate["status"], detail=gate["detail"])
                                       for gate in result["completion"]["gates"]]
            result["warnings"] = list(dict.fromkeys(self.warnings))
            self._snapshot, self._refreshed = result, time.monotonic()
            return result

    def _completion(self):
        record = self._read(self.out / "completion_state.json")
        if record is None:
            return None
        statuses = {"complete", "awaiting_measurements", "requires_primary_source", "not_established"}
        valid = (isinstance(record, dict) and record.get("schema_version") == 1
                 and isinstance(record.get("status"), str)
                 and isinstance(record.get("computational_release_ready"), bool)
                 and isinstance(record.get("scientific_validation_complete"), bool)
                 and isinstance(record.get("gates"), list)
                 and all(isinstance(gate, dict) and all(isinstance(gate.get(key), str)
                         for key in ("id", "label", "status", "detail"))
                         and gate["status"] in statuses for gate in record.get("gates", [])))
        if not valid:
            self.warnings.append("Ignored completion_state.json: unsupported or incomplete completion gates")
            return None
        completion = dict(record, gates=[dict(gate) for gate in record["gates"]])
        if completion["scientific_validation_complete"] and (
                not completion["gates"] or any(gate["status"] != "complete" for gate in completion["gates"])):
            completion["scientific_validation_complete"] = False
            self.warnings.append("Scientific completion was declared with unresolved gates; fresh validation remains open")
        return completion

    def _finalization(self, now, logs):
        """Read follow-up studies without mixing pilots into accepted-model scores."""
        registry = self._read(self.out / "finalization_registry.json")
        if not registry:
            return None
        tracks = []
        exact_records = []
        for declaration in registry.get("tracks", []):
            root = self.out / declaration["root"]
            if not root.resolve().is_relative_to(self.out.resolve()):
                self.warnings.append("Ignored a finalization path outside the artifact directory")
                continue
            state = self._read(root / "state.json", {})
            pilot = declaration.get("kind") == "pilot"
            gate = self._read(root / "review_gate.json", {}) if pilot else {}
            if pilot and state.get("status") in ("finished", "complete"):
                for filename in gate.get("additional_track_states", []):
                    path = Path(filename)
                    if not path.is_absolute():
                        path = root / path
                    if path.resolve().is_relative_to(self.out.resolve()):
                        additional = self._read(path, {})
                        if additional:
                            state = additional
            manifest = self._read(root / "manifest.json", {})
            phase = state.get("phase")
            phase_jobs = self._read(root / f"jobs_{phase}.json") if phase else None
            jobs = phase_jobs if isinstance(phase_jobs, list) else manifest.get("jobs", [])
            if not jobs and isinstance(state.get("jobs"), list):
                jobs = state["jobs"]
            if not jobs and state.get("status") == "finished":
                jobs = [job for path in root.glob("jobs_*.json") for job in self._read(path, [])]
            terminal = [row for path in (root / "attempts").glob("*.json")
                        if isinstance(row := self._read(path), dict) and _date(row.get("finished_utc"))
                        and (not phase or row.get("phase", phase) == phase)]
            terminal = list({row.get("tag", str(index)): row for index, row in enumerate(terminal)}.values())
            if declaration["id"] == "atlas_exact":
                exact_records = terminal
            total = max(len(jobs), int(state.get("total", state.get("total_jobs", 0))), len(terminal))
            finished = len(terminal)
            accepted = sum(bool(row.get("accepted")) for row in terminal)
            failed = len(terminal)-accepted
            if pilot:
                finished = int(state.get("finished", 0))
                total = int(state.get("total", len(jobs)))
                accepted = 0  # A bounded joint pilot is not a validated final model.
                failed = sum(row.get("status") in ("failed", "error") for row in jobs)
            states = [state]
            for filename in declaration.get("supplemental_states", []):
                path = root / filename
                if path.resolve().is_relative_to(self.out.resolve()):
                    states.append(self._read(path, {}))
            live_states = [s for s in states if s.get("status") in ("running", "auditing", "waiting")]
            running = bool(live_states)
            labels = [item for s in live_states for item in s.get("active_jobs", [])]
            if not labels and running:
                ended_tags = {row.get("tag") for row in terminal}
                labels = [row for path in (root / "launches").glob("*.json")
                          if isinstance(row := self._read(path), dict) and row.get("tag") not in ended_tags]
            activity = []
            observed_tags = set()
            track_logs = logs + list((root / "logs").glob("*.log")) + list(root.glob("*.log"))
            for item in labels if running else []:
                if isinstance(item, str):
                    tag = item
                    job = self._read(root / "jobs" / (tag+".json"), {})
                    launch = self._read(root / "worker_runs" / (tag+".json"), {})
                    if not launch:
                        launch = self._read(root / "worker_runs" / (tag+"_supplemental.json"), {})
                    if root.name == "p13_finalization_stability_v1" and tag.startswith("p13_model_exploration_v1_"):
                        tag = tag.replace("p13_model_exploration_v1_", root.name+"_", 1)
                    launch = dict(launch, tag=tag)
                else:
                    job = item
                    tag = item.get("tag", item.get("id", "unknown"))
                    launch = dict(item, tag=tag)
                if tag in observed_tags:
                    continue
                observed_tags.add(tag)
                job = dict(job, cell=job.get("cell", {"id": tag}), fold=job.get("fold", "future"))
                active = self._active(job, launch, declaration["label"], now, track_logs)
                if pilot:
                    active["stage"] = "Joint fit" if job.get("stage") == "fit" else "Independent numerical audits"
                    active["tag"] = job.get("id", tag)
                activity.append(active)
            active_count = (int(state.get("active", len(activity))) if pilot else len(activity)) if running else 0
            status = state.get("status", "prepared" if manifest else "not started")
            if live_states and status in ("finished", "complete"):
                status = "running"
            if pilot and gate and not gate.get("ready_for_review") and status in ("finished", "complete"):
                status = "awaiting numerical review"
            tracks.append(dict(id=declaration["id"], label=declaration["label"],
                               status=status,
                               phase=phase or ("Joint fits and audits" if pilot else "Refit and audit"),
                               total=total, finished=finished, accepted=accepted, failed=failed,
                               remaining=max(0, total-finished), active=active_count,
                               recently_active=sum(row["activity"] == "recent" for row in activity),
                               workers=sum(int(s.get("workers", 0)) for s in live_states),
                               active_jobs=activity, detail=declaration.get("detail", ""),
                               conditional_validation=len(manifest.get("candidate_validation", [])) if phase in ("combine", "exact_screen") else 0,
                               updated_utc=state.get("updated_utc"), pilot=pilot))
        review = None
        if registry.get("review_root"):
            folder = self.out / registry["review_root"]
            if folder.resolve().is_relative_to(self.out.resolve()):
                report = self._read(folder / "review.json", {})
                verification = self._read(folder / "verification_state.json", {})
                review = dict(status=report.get("status", "waiting"),
                    numerical_review_ready=bool(report.get("model_finalization_ready")),
                    snapshot_verification=verification.get("status", "waiting"),
                    snapshot_verified=bool(verification.get("passed")),
                    report_path=str(folder / "review.md"),
                    detail="Final selection requires numerical repair, seven-fold evidence and prediction replay.")
                candidate = report.get("recommendations", {}).get("absolute_accuracy", {})
                if review["numerical_review_ready"] and review["snapshot_verified"] and candidate.get("eligible") and candidate.get("retained_atlas_exact_J"):
                    rows = [row for row in exact_records if row.get("cell", {}).get("id") == candidate.get("model")
                            and row.get("phase") == "validate" and row.get("seed") in (0, 2)]
                    identities = {(row.get("fold"), row["seed"]) for row in rows}
                    expected = {(fold, seed) for fold in candidate.get("fold_results", {}) for seed in (0, 2)}
                    if len(expected) == 14 and identities == expected:
                        folds = [dict(fold=row["fold"], seed=row["seed"], accepted=bool(row.get("accepted")),
                                 absolute_rms_MHz=_finite(row.get("last", {}).get("score", {}).get("unique_absolute_comb", {}).get("rms_MHz")))
                                 for row in rows]
                        review["leader"] = dict(model=candidate["model"], free=candidate["free"],
                            inherited=candidate["inherited"], calibration=4,
                            absolute_rms_MHz=candidate["absolute_rms_MHz"], intraline_rms_MHz=candidate["intraline_rms_MHz"],
                            validation=dict(total=14, accepted=sum(row["accepted"] for row in folds), folds=folds),
                            source="verified_exact_J_research_snapshot")
        return dict(tracks=tracks, review=review, worker_budget=registry.get("worker_budget", 12),
                    active=sum(track["active"] for track in tracks),
                    recently_active=sum(track["recently_active"] for track in tracks),
                    total=sum(track["total"] for track in tracks), finished=sum(track["finished"] for track in tracks),
                    remaining=sum(track["remaining"] for track in tracks),
                    completed_tracks=sum(track["status"] == "finished" for track in tracks),
                    notes=registry.get("notes", ""))

    @staticmethod
    def _phase(campaign, phase):
        return next((dict(finished=p["finished"], total=p["total"]) for p in campaign["phases"] if p["id"] == phase),
                    dict(finished=0, total=0))

    def _campaign(self, prefix, label, extra, now, changed, logs):
        root = self.out / prefix
        queue_path = root / "jobs_validate.json"
        jobs = self._read(queue_path, [])
        available = queue_path.exists() and bool(jobs)
        records = {}
        phases = {}
        for path in (root / "attempts").glob("*.json"):
            row = self._read(path)
            if not isinstance(row, dict) or not _date(row.get("finished_utc")):
                continue
            records[path.stem] = row
            phase = row.get("phase", "unknown")
            counts = phases.setdefault(phase, dict(id=phase, finished=0, accepted=0, failed=0, total=0))
            counts["finished"] += 1
            counts["accepted" if row.get("accepted") else "failed"] += 1
        for phase, counts in phases.items():
            planned = self._read(root / f"jobs_{phase}.json", [])
            counts["total"] = len(planned) if planned else counts["finished"]
        state = self._read(root / "state.json", {})
        states = [state]
        if extra:
            states.append(self._read(self.out / extra / "state.json", {}))
        workers = sum(int(s.get("workers", 0)) for s in states if s.get("status") == "running")
        completed, accepted, active, pending = [], [], [], []
        for job in jobs:
            tag = _tag(job, prefix)
            record = records.get(tag)
            if record:
                item = dict(tag=tag, model=job["cell"]["id"], fold=job["fold"], seed=job.get("seed"),
                            accepted=bool(record.get("accepted")), finished_utc=record["finished_utc"],
                            campaign=label, error=record.get("error"))
                completed.append(item)
                if record.get("accepted"):
                    accepted.append(record)
                continue
            pending.append(job)
            for owner in (prefix, extra):
                if owner is None:
                    continue
                launch = self._read(self.out / owner / "launches" / (_tag(job, owner) + ".json"))
                if launch:
                    active.append(self._active(job, launch, label, now, logs))
                    break
        rates = []
        validation_records = [records[_tag(j, prefix)] for j in jobs if _tag(j, prefix) in records]
        first = min((_date(r["started_utc"]) for r in validation_records if _date(r.get("started_utc"))), default=now)
        for hours in (6, 12, 24):
            start = max(now - timedelta(hours=hours), first)
            if changed and changed > start:
                start = min(changed, now)
            duration = (now - start).total_seconds() / 3600
            count = sum(_date(r["finished_utc"]) >= start for r in completed)
            rates.append(dict(window_hours=hours, observed_hours=round(duration, 2), count=count,
                              per_day=round(count * 24 / duration, 2) if duration >= 1 else None))
        usable = [r["per_day"] for r in rates if r["observed_hours"] >= 6 and r["per_day"] and r["per_day"] > 0]
        remaining = len(pending)
        # A spread across observed windows plus 20% slack exposes long-fit risk.
        eta = (dict(days_low=round(remaining / max(usable), 2),
                    days_high=round(remaining / min(usable) * 1.2, 2)) if usable and remaining else None)
        if available and not remaining:
            eta = dict(days_low=0, days_high=0)
        comparison = self._read(root / "comparison.json", {})
        validation = {}
        for job in jobs:
            group = validation.setdefault(job["cell"]["id"], dict(total=0, finished=0, accepted=0, folds=[]))
            group["total"] += 1
            record = records.get(_tag(job, prefix))
            if record:
                group["finished"] += 1
                group["accepted"] += bool(record.get("accepted"))
        for row in validation_records:
            group = validation.get(row["cell"]["id"])
            if group is not None:
                last = row.get("last", {})
                score = last.get("score", {}).get("unique_absolute_comb", {})
                group["folds"].append(dict(fold=row["fold"], seed=row.get("seed"), accepted=bool(row.get("accepted")),
                    absolute_rms_MHz=_finite(score.get("rms_MHz")),
                    training_cost=_finite(last.get("potential_training_cost"))))
        return dict(id=prefix, label=label, available=available, status=state.get("status", "unavailable"),
                    total=len(jobs), finished=len(completed), accepted=len(accepted), failed=len(completed)-len(accepted),
                    remaining=remaining, active=len(active), queued=max(0, remaining-len(active)),
                    configured_workers=workers, recently_active=sum(r["activity"] == "recent" for r in active),
                    progress_pct=round(100*len(completed)/len(jobs), 2) if jobs else None,
                    latest_completion_utc=max((r["finished_utc"] for r in completed), default=None),
                    rates=rates, eta=eta, active_jobs=active, phases=list(phases.values()), validation=validation,
                    _completed=completed, _selected=[dict(r, _campaign=prefix) for r in comparison.get("selected", [])])

    def _active(self, job, launch, campaign, now, logs):
        tag = launch["tag"]
        latest, stamp = None, None
        for path in logs:
            if path.name.startswith(tag + "_") or path.name == tag + ".log":
                try:
                    modified = path.stat().st_mtime
                except OSError:
                    continue
                if stamp is None or modified > stamp:
                    latest, stamp = path, modified
        last_activity = datetime.fromtimestamp(stamp, timezone.utc) if stamp else None
        started = _date(launch.get("started_utc"))
        seconds = (now-last_activity).total_seconds() if last_activity else None
        return dict(tag=tag, model=job["cell"]["id"], fold=job["fold"], seed=job.get("seed"), campaign=campaign,
                    started_utc=launch.get("started_utc"), elapsed_hours=max(0, (now-started).total_seconds()/3600) if started else None,
                    last_activity_utc=last_activity.isoformat() if last_activity else None,
                    activity="recent" if seconds is not None and seconds < 900 else "quiet",
                    stage=_stage(latest.name[len(tag):]) if latest else "Waiting for a log update",
                    evaluation=_last_evaluation(latest) if latest else {})

    @staticmethod
    def _leader(selected, campaigns):
        candidates = []
        for row in selected:
            last = row.get("last", {})
            score = last.get("score", {}).get("unique_absolute_comb", {})
            rms = _finite(score.get("rms_MHz"))
            if rms is not None:
                candidates.append((rms, row))
        if not candidates:
            return None
        rms, row = min(candidates, key=lambda item: item[0])
        last, model = row["last"], row["cell"]["id"]
        validation = next((c["validation"][model] for c in campaigns
                           if c["id"] == row["_campaign"] and model in c["validation"]), None)
        return dict(model=model, absolute_rms_MHz=rms,
                    intraline_rms_MHz=_finite(last.get("hyperfine_held_out", {}).get("rms_MHz")),
                    free=last.get("free"), inherited=last.get("inherited"), calibration=4,
                    validation=validation)


def main():
    parser = argparse.ArgumentParser(description="Live terminal dashboard for the model search")
    parser.add_argument("--out", type=Path, default=OUT, help="research artifact directory")
    parser.add_argument("--interval", type=float, default=30, help="refresh interval in seconds (minimum 5)")
    parser.add_argument("--json", action="store_true", help="print a single JSON snapshot and exit")
    args = parser.parse_args()
    if not math.isfinite(args.interval) or args.interval < 5:
        parser.error("Use an interval of at least 5 seconds")
    if args.json:
        print(json.dumps(ProgressMonitor(args.out).snapshot(), allow_nan=False, indent=2))
    else:
        from .progress_tui import run
        run(out=args.out, interval=args.interval)


if __name__ == "__main__":
    main()
