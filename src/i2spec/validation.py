"""Freeze research forecasts and score later measurements without refitting.

Empirical envelopes have no assigned coverage probability. Reused calibration
data, independently published checks and prospective laboratory measurements
are kept distinct. Registries are immutable; results are separate new files.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import statistics

EXCLUSION = "127I2 P(13) 43-0"


def _utc():
    return datetime.now(timezone.utc).isoformat()


def _datetime(value):
    try:
        parsed = datetime.fromisoformat(value.replace("Z", "+00:00"))
        if parsed.tzinfo is None:
            raise ValueError("missing timezone")
        return parsed
    except (AttributeError, TypeError, ValueError) as error:
        raise ValueError("measurement timestamps need ISO8601 with an explicit timezone") from error


def _sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()


def _new_json(path, value):
    path = Path(path)
    with path.open("x", encoding="utf-8") as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.write("\n")


def _number(value, name, *, positive=False):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValueError(f"{name} must be a finite number")
    if positive and value <= 0:
        raise ValueError(f"{name} must be positive")
    return float(value)


def _target(raw):
    from .observations import Line, component_rank
    target = dict(raw)
    if not isinstance(target.get("id"), str) or not target["id"].strip():
        raise ValueError("each target needs a nonempty id")
    line = Line.parse(target["line"])
    if line.isotopologue != "127I2" or str(line) == EXCLUSION:
        raise ValueError("target violates the reference-isotope domain or P(13) exclusion")
    rank = component_rank(target["component"])
    if rank < 1:
        raise ValueError("component rank must be positive")
    target.update(line=str(line), component=f"a{rank}")
    return target


def _envelopes(raw):
    result = {}
    for name, value in raw.items():
        item = dict(value)
        item["half_width_MHz"] = _number(item["half_width_MHz"], "envelope half-width", positive=True)
        if item.get("nominal_probability") is not None:
            raise ValueError("these empirical envelopes do not establish nominal coverage probabilities")
        if not item.get("scope"):
            raise ValueError("an empirical envelope needs a declared scope")
        item["nominal_probability"] = None
        result[name] = item
    return result


def _research_factory(name):
    from .research import ResearchModel
    model = ResearchModel(name)
    return model.info.id, model.info.as_dict(), model.coefficients, model.frequency


def freeze_registry(destination, plan, *, factory=None):
    """Create a new registry; model choices and scoring precede outcomes.

    ``factory`` allows an external predictor to provide the same frozen
    metadata/coefficients/frequency contract; normal use selects ResearchModel.
    A failed preparation leaves an incomplete directory without manifest.json.
    """
    destination = Path(destination)
    if destination.exists():
        raise FileExistsError("a forecast registry cannot be overwritten; choose a new destination")
    if not plan.get("models") or len(set(plan["models"])) != len(plan["models"]):
        raise ValueError("declare unique model names before prediction")
    if plan.get("primary_model") not in plan["models"]:
        raise ValueError("the primary model must be declared before outcomes")
    if plan.get("outcome_phase") not in ("prospective_laboratory", "external_published"):
        raise ValueError("declare prospective_laboratory or external_published outcome phase")
    targets = [_target(raw) for raw in plan["targets"]]
    if not targets or len({t["id"] for t in targets}) != len(targets):
        raise ValueError("declare nonempty unique target ids")
    if len({(t["line"], t["component"]) for t in targets}) != len(targets):
        raise ValueError("repeated physical identities do not create additional targets")
    envelopes = _envelopes(plan.get("envelopes", {}))
    thresholds = [_number(value, "accuracy threshold", positive=True) for value in plan["thresholds_MHz"]]
    if len(set(thresholds)) != len(thresholds) or not thresholds:
        raise ValueError("declare unique, nonempty accuracy thresholds")
    factory = factory or _research_factory
    models, forecasts, model_ids = [], [], []
    primary = None
    for name in plan["models"]:
        identity, info, coefficients, frequency = factory(name)
        if identity in model_ids:
            raise ValueError("different model aliases resolve to the same frozen model")
        model_ids.append(identity)
        if name == plan["primary_model"]:
            primary = identity
        models.append(dict(id=identity, info=info,
            coefficient_sha256=hashlib.sha256(_canonical(coefficients)).hexdigest()))
        for target in targets:
            label = target["line"].split(" ", 1)[1]
            predicted = _number(frequency(label, target["component"]), "predicted frequency", positive=True)
            forecasts.append(dict(model=identity, target_id=target["id"], predicted_MHz=predicted,
                envelopes={key:dict(lower_MHz=predicted-rule["half_width_MHz"],
                    upper_MHz=predicted+rule["half_width_MHz"], nominal_probability=None)
                    for key, rule in envelopes.items()}))
    # Data/coefficients and all point forecasts are computed before creating the
    # destination. The manifest is the commit marker for a complete registry.
    destination.mkdir(parents=True, exist_ok=False)
    protocol = dict(plan, schema_version=1, created_utc=_utc(), primary_model=primary,
        models=model_ids, targets=targets, envelopes=envelopes, thresholds_MHz=thresholds,
        refitting_allowed=False, post_outcome_primary_selection_allowed=False,
        exclusion=EXCLUSION, status="frozen",
        interval_interpretation="Empirical discrepancy envelopes; nominal probability unspecified")
    _new_json(destination/"protocol.json", protocol)
    _new_json(destination/"models.json", models)
    _new_json(destination/"forecasts.json", forecasts)
    files = ("protocol.json", "models.json", "forecasts.json")
    manifest = dict(schema_version=1, created_utc=_utc(), files_sha256={name:_sha(destination/name) for name in files})
    _new_json(destination/"manifest.json", manifest)
    return manifest


def read_registry(folder):
    """Read a frozen registry, including the first dated literature registry."""
    folder = Path(folder)
    if (folder/"manifest.json").exists():
        manifest = json.loads((folder/"manifest.json").read_text())
        expected = manifest["files_sha256"]
        if set(expected) != {"protocol.json", "models.json", "forecasts.json"}:
            raise ValueError("unexpected or incomplete registry manifest")
        for name, digest in expected.items():
            if _sha(folder/name) != digest:
                raise ValueError("frozen registry hash mismatch: "+name)
        protocol = json.loads((folder/"protocol.json").read_text())
        forecasts = json.loads((folder/"forecasts.json").read_text())
        fingerprint = _sha(folder/"manifest.json")
    else:
        # This pre-integration registry was frozen before inspecting the new
        # Kondo2025 outcome. Preserve its exact bytes and exposure chronology.
        expected = json.loads((folder/"frozen_registry_sha256.json").read_text())
        required = {"point_predictions_v1.json", "validation_protocol_v1.json", "interval_endpoints_v1.json"}
        if set(expected) != required:
            raise ValueError("unexpected legacy registry manifest")
        for name, digest in expected.items():
            if _sha(folder/name) != digest:
                raise ValueError("frozen registry hash mismatch: "+name)
        points = json.loads((folder/"point_predictions_v1.json").read_text())
        source = json.loads((folder/"validation_protocol_v1.json").read_text())
        endpoints = json.loads((folder/"interval_endpoints_v1.json").read_text())["rows"]
        protocol = dict(source, primary_model=source["primary_role"], outcome_phase="external_published",
            exclusion=points["exclusion"], refitting_allowed=source["coefficient_refitting_allowed"],
            post_outcome_primary_selection_allowed=source["post_outcome_model_selection_allowed"],
            thresholds_MHz=source["external_accuracy_thresholds_MHz"],
            targets=[dict(t,line=f"127I2 {t['line'][0]}({t['line'][1]}) {t['line'][2]}-{t['line'][3]}") for t in points["targets"]])
        forecasts = [dict(model=r["role"],target_id=r["target_id"],predicted_MHz=r["predicted_MHz"],
            envelopes=r["envelopes"]) for r in endpoints]
        fingerprint = _sha(folder/"frozen_registry_sha256.json")
    if protocol.get("exclusion") != EXCLUSION:
        raise ValueError("registry exclusion changed")
    if protocol.get("refitting_allowed") or protocol.get("post_outcome_primary_selection_allowed"):
        raise ValueError("registry permits an invalid post-outcome model change")
    _envelopes(protocol["envelopes"])
    for target in protocol["targets"]:
        _target(target)
    targets = {target["id"] for target in protocol["targets"]}
    models = {row["model"] for row in forecasts}
    if protocol["primary_model"] not in models:
        raise ValueError("registry has no primary-model forecasts")
    if {(row["model"], row["target_id"]) for row in forecasts} != {(m,t) for m in models for t in targets} or len(forecasts) != len(models)*len(targets):
        raise ValueError("registry forecast identities are incomplete or duplicated")
    for row in forecasts:
        _number(row["predicted_MHz"], "frozen frequency", positive=True)
        if set(row["envelopes"]) != set(protocol["envelopes"]):
            raise ValueError("registry forecast envelopes are incomplete")
        for name, rule in protocol["envelopes"].items():
            endpoints = row["envelopes"][name]
            for side, sign in (("lower_MHz", -1), ("upper_MHz", 1)):
                if abs(endpoints[side]-(row["predicted_MHz"]+sign*rule["half_width_MHz"])) > 1e-6:
                    raise ValueError("registry envelope endpoints differ from the frozen protocol")
    return protocol, forecasts, fingerprint


def _standard_uncertainty(row):
    if row.get("standard_uncertainty_MHz") is not None:
        return _number(row["standard_uncertainty_MHz"], "standard measurement uncertainty", positive=True)
    if row.get("expanded_uncertainty_MHz") is not None:
        expanded = _number(row["expanded_uncertainty_MHz"], "expanded measurement uncertainty", positive=True)
        if row.get("coverage_factor") is not None:
            return expanded/_number(row["coverage_factor"], "coverage factor", positive=True)
    return None


def _metrics(rows, thresholds, widths):
    if not rows:
        return dict(n=0, bias_MHz=None, rms_MHz=None, maximum_abs_error_MHz=None,
            thresholds={}, empirical_envelopes={})
    errors = [row["error_MHz"] for row in rows]
    return dict(n=len(errors), bias_MHz=statistics.mean(errors),
        rms_MHz=math.sqrt(statistics.mean(error*error for error in errors)),
        maximum_abs_error_MHz=max(abs(error) for error in errors),
        thresholds={str(value):dict(within=sum(abs(error)<=value for error in errors), total=len(errors)) for value in thresholds},
        empirical_envelopes={name:dict(within=sum(abs(error)<=width for error in errors), total=len(errors),
            half_width_MHz=width, nominal_probability=None) for name,width in widths.items()})


def _connected_blocks(measurements, targets):
    """Count connected line/source groups; repeats are not independent trials."""
    parent = {}
    def root(key):
        parent.setdefault(key, key)
        while parent[key] != key:
            parent[key] = parent[parent[key]]
            key = parent[key]
        return key
    for row in measurements:
        line = ("line", targets[row["target_id"]]["line"])
        source = ("source", row["source_id"])
        parent[root(line)] = root(source)
    return len({root(key) for key in parent})


def score_measurements(folder, measurements):
    """Compare later outcomes with immutable forecasts; never fit a correction."""
    protocol, forecasts, fingerprint = read_registry(folder)
    targets = {target["id"]:target for target in protocol["targets"]}
    if not isinstance(measurements, list) or any(not isinstance(row, dict) for row in measurements):
        raise ValueError("measurements must be a list of objects")
    if len({row.get("id") for row in measurements}) != len(measurements):
        raise ValueError("measurement rows need unique ids")
    usable, pending, noncomparable = [], [], []
    for original in measurements:
        row = dict(original)
        if not isinstance(row.get("id"), str) or not row["id"]:
            raise ValueError("measurement needs a nonempty id")
        if row.get("target_id") not in targets:
            raise ValueError("measurement identity was not preregistered: "+str(row.get("target_id")))
        target = targets[row["target_id"]]
        checked = _target(dict(id="check", line=row.get("line",target["line"]),
            component=row.get("component",target["component"])))
        if any(checked[key] != target[key] for key in ("line", "component")):
            raise ValueError("measurement line or component differs from its registered identity")
        if row.get("value_MHz") is None:
            pending.append(row)
            continue
        row["value_MHz"] = _number(row["value_MHz"], "measured frequency", positive=True)
        row["standard_uncertainty_MHz"] = _standard_uncertainty(row)
        if row.get("comparable") is not True:
            noncomparable.append(row)
            continue
        if not row.get("source_id") or not row.get("source_url") or not row.get("conditions"):
            raise ValueError("a comparable measurement needs source_id, source_url and measurement conditions")
        if row.get("frequency_basis") != "zero_pressure_unperturbed":
            raise ValueError("comparable outcomes need a source-supported zero_pressure_unperturbed frequency basis")
        if protocol["outcome_phase"] == "prospective_laboratory":
            if row.get("independent_of_fitting_archive") is not True or row.get("outcome_unseen_before_freeze") is not True:
                raise ValueError("prospective outcomes need explicit independence and pre-freeze exposure declarations")
            acquired = _datetime(row.get("acquired_utc"))
            if acquired < _datetime(protocol["created_utc"]) or acquired > datetime.now(timezone.utc):
                raise ValueError("prospective measurement must be acquired after forecast freezing and before scoring")
        usable.append(row)
    by_target = defaultdict(list)
    for row in usable:
        by_target[row["target_id"]].append(row)
    scored, summaries = [], {}
    thresholds = protocol["thresholds_MHz"]
    widths = {name:rule["half_width_MHz"] for name,rule in protocol["envelopes"].items()}
    models = sorted({forecast["model"] for forecast in forecasts})
    for model in models:
        unique = []
        for forecast in forecasts:
            if forecast["model"] != model or forecast["target_id"] not in by_target:
                continue
            target = targets[forecast["target_id"]]
            observed = by_target[forecast["target_id"]]
            predicted = forecast["predicted_MHz"]
            errors = [predicted-row["value_MHz"] for row in observed]
            unique.append(dict(target_id=target["id"],line=target["line"],
                error_MHz=predicted-statistics.median(row["value_MHz"] for row in observed),
                maximum_repeat_error_MHz=max(abs(error) for error in errors)))
            for row, error in zip(observed, errors):
                sigma = row["standard_uncertainty_MHz"]
                scored.append(dict(model=model, primary=model==protocol["primary_model"],measurement_id=row["id"],
                    target_id=target["id"],line=target["line"],component=target["component"],
                    predicted_MHz=predicted, observed_MHz=row["value_MHz"],error_MHz=error,
                    measurement_standard_uncertainty_MHz=sigma,
                    residual_over_measurement_sigma=error/sigma if sigma else None,
                    inside_empirical_envelopes={name:abs(error)<=width for name,width in widths.items()}))
        line_groups, interval_rows = defaultdict(list), []
        for row in unique:
            line_groups[row["line"]].append(row)
        for line, group in line_groups.items():
            ordered = sorted(group, key=lambda r: int(targets[r["target_id"]]["component"][1:]))
            for row in ordered[1:]:
                interval_rows.append(dict(line=line, target_id=row["target_id"], ref_target_id=ordered[0]["target_id"],
                    error_MHz=row["error_MHz"]-ordered[0]["error_MHz"]))
        summary = _metrics(unique, thresholds, widths)
        summary.update(n_measurements=sum(len(by_target[row["target_id"]]) for row in unique),
            n_lines=len(line_groups), n_sources=len({r["source_id"] for r in usable}),
            n_line_source_connected_blocks=_connected_blocks(usable, targets),
            per_stratum={name:_metrics([r for r in unique if targets[r["target_id"]].get("stratum", "unspecified")==name],thresholds,widths)
                for name in sorted({targets[r["target_id"]].get("stratum", "unspecified") for r in unique})},
            line_block_thresholds={str(value):dict(within=sum(all(abs(r["error_MHz"])<=value for r in group) for group in line_groups.values()),
                total=len(line_groups)) for value in thresholds},
            maximum_repeat_error_MHz=max((r["maximum_repeat_error_MHz"] for r in unique),default=None),
            derived_intraline_intervals=dict(rows=interval_rows, metrics=_metrics(interval_rows, thresholds, {}),
                interpretation="Differences of median component observations sharing the lowest measured rank; correlated, with no independence or measurement-uncertainty propagation assumed"))
        summaries[model] = summary
    return dict(schema_version=1, created_utc=_utc(), registry_sha256=fingerprint,
        measurements_sha256=hashlib.sha256(_canonical(measurements)).hexdigest(),
        primary_model=protocol["primary_model"],outcome_phase=protocol["outcome_phase"],
        registered_targets=len(targets),scored_unique_targets=len(by_target),
        missing_target_ids=sorted(set(targets)-set(by_target)),
        rows=scored,summaries=summaries,pending_measurements=pending,noncomparable_measurements=noncomparable,
        refitted=False,model_selected_using_outcomes=False,
        interpretation="Descriptive external errors and empirical-envelope inclusion; no nominal confidence coverage inferred",
        prospective_measurement_plan_complete=protocol["outcome_phase"]=="prospective_laboratory" and len(by_target)==len(targets),
        model_adequacy="Requires scientific review of residuals, strata, conditions and tolerances; collecting all targets alone does not validate the model")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    freeze = commands.add_parser("freeze", help="freeze a plan before outcomes")
    freeze.add_argument("plan",type=Path)
    freeze.add_argument("destination",type=Path)
    verify = commands.add_parser("verify", help="check registry file integrity")
    verify.add_argument("registry",type=Path)
    score = commands.add_parser("score",help="score later measurements without altering forecasts")
    score.add_argument("registry",type=Path)
    score.add_argument("measurements",type=Path)
    score.add_argument("--out",required=True,type=Path)
    args = parser.parse_args(argv)
    try:
        if args.command=="freeze":
            result=freeze_registry(args.destination,json.loads(args.plan.read_text()))
        elif args.command=="verify":
            protocol,forecasts,digest=read_registry(args.registry)
            result=dict(passed=True,registry_sha256=digest,primary_model=protocol["primary_model"],forecasts=len(forecasts))
        else:
            result=score_measurements(args.registry,json.loads(args.measurements.read_text()))
            _new_json(args.out,result)
        print(json.dumps(result,indent=2,allow_nan=False))
    except (OSError,ValueError,KeyError,TypeError) as error:
        parser.exit(2,f"validation: {error}\n")
    return 0


if __name__=="__main__":
    raise SystemExit(main())
