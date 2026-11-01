"""Rebuild release 20261031, the released model, from the repository: its fits in their final form, each stage on the
previous one's output, written as a candidate release (a folder) whose models the later stages read by role, as the
release's own build did.

  1  core          potentials (core.fit) on the core rows and the six CIPM isotopologue line centers
                   (data/training/20261029/isotopologue_centers.json): comb floor 0.3 MHz, the B isotope shift of power
                   5, the B limit tied to D0(X) + the atomic splitting; then the B hyperfine functions with eqQ order 6
  2  continuation  tail terms (B rotation order 2, the release's structure) on that core, comb floor 0.3 MHz
  3  hyperfine     beyond the core in the release's form (hyperfine_beyond.FORM_20261031) at every 1g(1Pig) C5 scale,
                   the lowest training cost kept; release 20261030's fit (FIT_KEYS at the scale 0.65), which the B~1g
                   resonance at B v' 77, 78 keeps; the X eqQ and C (x_hyperfine)
  4  package       levels with the barrier rule on every grid, the hyperfine tables, the resonance block, the levels a
                   near-degenerate second-order term spoils (package.spike_levels: not predicted), the references
  5  uncertainty   uncertainty.build on release 20261029's cross-validation with this form's nested hyperfine folds
                   (data/training/20261031/hyperfine_folds.json.gz), and the named parameters' covariance
  6  isotopologues the isotope shift (isotopologues.fit_shift) and its record
  7  calibrate     the calibration refit (calibrate.refit, release 20261030's forms) on the candidate with the
                   build's table; the X term near dissociation (calibrate.near_dissociation) on release 20261030's
                   table, which it was fitted to; the rules beyond the tables (data/training/20261031/edges.json);
                   then the isotopologue record conditioned on the CIPM lines and the references of every isotopologue

The refits start from the release's own coefficients (a refit returns them: tests/test_fit_final.py), and, as in
release.rebuild, the stages after the core build on the released core: the core refit moves no coefficient by more
than 1e-6 of itself (its cost falls by 1e-8 of itself, a flat optimum), and the continuation refit on the released core
returns the released tails exactly. The text of the release files (ids, titles, descriptions, inventories) is the
installed release's; every number is rebuilt.

The release read its second-order terms (hyperfine_beyond.fill_terms) from a cache whose key did not hold the B
potential, so each came from whichever earlier computation had filled it: every state's (and the 1g(1Pig)'s at
release 20261030's scale) on the potential release 20261021's hyperfine stage was fitted on, the 1g(1Pig)'s at the
fitted scale on the release's continuation (`s2_frozen` states that rule), except 390 levels of the 1g(1Pig) at the
fitted scale (among them every level the fit reads) and one level, B v' 43 J' 12, of six states, which no committed
potential reproduces (they differ by up to 38 MHz in eqQ or kHz in C, d, delta): those are frozen
(data/training/20261031/s2_terms.json.gz, `s2_fixed`). The cache key now holds the potential.
"""
from __future__ import annotations

import copy
import gzip
import hashlib
import json
import tempfile
import time
from dataclasses import replace
from functools import lru_cache
from pathlib import Path

import numpy as np

from .. import named
from . import b1g_resonance as rz, base, calibrate, continuation, core, core_hyperfine, hyperfine_beyond as hb
from . import isotopologues as iso_fit, package as pk, reduction, rows as training, uncertainty as U, x_hyperfine as xh
from .curves import Config, Curve, initial
from .release import hyperfine_rows_release, x_model_record

RELEASE = "20261031"
IDS = {k: f"research-{RELEASE}-{k}" for k in ("core", "expanded", "i2spec")}
INSTALLED = training.ROOT/"src"/"i2spec"/"data"/"research"/RELEASE
TRAINING = training.TRAINING/RELEASE
COMB_FLOOR = .3
ISO_SHIFT = 5
CORE_ORDERS = {"X": {}, "B": dict(core_hyperfine.ORDERS["B"], eqQ=6)}
FORM = hb.FORM_20261031
NEAR_DISSOCIATION_TABLE = training.ROOT/"src"/"i2spec"/"data"/"research"/"20261030"/"expanded"/"uncertainty.json"
SPIKE_REASON = ("a second-order energy denominator nearly vanishes at this level (a perturber level of the inherited "
                "curves lies within a few GHz), so the second-order sum is not valid and the two levels would have to "
                "be diagonalized; unlike B v' 77/78 no measured perturber level exists to place them, so the model "
                "gives this level's center only")
SPIKE_RULE = ("an entry departs from the median of its J +- 1..3 neighbors at the same v by more than 5 x max(10 MAD, "
              "5 % of the median, 5 MHz eqQ / 50 kHz C d delta)")
RESONANCE_NOTE = ("release 20261030's hyperfine fit (its 1g mixing angle 30 deg, 1g(1Pig) C5 scale 0.65 and "
                  "backgrounds): the 1g coupling of the earlier fit, kept at the resonance; the other levels beyond "
                  "the core take this release's form")
B_MODEL_NOTE = ("the release's hyperfine fit (form of release 20261031: theta_g, curvatures of C, d, delta, background "
                "orders d 5, delta 4)")
CALIBRATION_RECIPE = (
    " Release 20261030 (the calibration refit): sigma_model and sigma_HF refitted on the cross-validation residuals "
    "above and the calibrate sets' residuals together (the final form's residuals on rows it never trained on: "
    "2,432 rows, without compilation rows a held set also gives, lines the training holds inside a CIPM mean, record "
    "rows, diagnostic exclusions and values printed again in a later set). In the calibrate rows each level is a "
    "random effect of variance sigma_level^2 and each lab family a random offset of the absolute frequencies "
    "(measurement error, not in predictions); in the intraline intervals the pattern-scale error is one per line "
    "and the floor one per component. Two terms beyond the core, each >= 0: the B hyperfine floor grows as "
    "exp(g (v' - 44)) up to the last validated v' (a_growth_per_v); B sigma_model beyond the core has a (v - v0)^2 "
    "term (quadratic_per_v2: the likelihood is flat in it, so it is chosen by an explicit profile) and does not fall "
    "below the core's per-level share at v' 44. The isotopologue floors are the release's.")
RECIPE_20261031 = (
    " Release 20261031: the same calibration refit on the release's hyperfine beyond the core (theta_g, curvatures, "
    "background orders; its nested cross-validation) with the calibrate residuals recomputed with it. Added: X "
    "v'' >= 85 a per-level term exp(a + c (v - 85)) in quadrature (near_dissociation, fitted by maximum likelihood "
    "to the cross-validation records at v'' >= 85); the rules for the levels beyond the tables (edges: B v' 81-86 "
    "the v' = 80 share times 3 per quantum, a rule calibrated on three levels of one source; X v'' 109-113 the "
    "near-dissociation term continued, capped at sqrt(binding^2 + sigma(D0)^2)); and the tunneling width of a "
    "quasi-bound B level (parabolic barrier) in quadrature with its share.")
ISOTOPOLOGUE_FLOORS_KEPT = "the release's: the calibrate rows are covered with or without them"
ATLAS_CHECK = ("the 1978 atlas volumes (B v' <= 62) on their own scale (offset and slope on their core-region lines, "
               "scatter 35-51 MHz) are within 2 sigma for 94-98 % of the lines at B v' 44-65; they cannot set "
               "sigma_model")
CM_TO_MHZ = 29979.2458


def s2_frozen(state, constants):
    """The B potential of release 20261031's second-order terms (module note): the frozen 20261021 one for every
    state but the 1g(1Pig) at a C5 scale other than release 20261030's."""
    return not (state == rz.STATE and constants.get("C5") != hb.FIT_OVERRIDES[rz.STATE]["C5"])


@lru_cache(maxsize=None)
def _s2_terms():
    with gzip.open(TRAINING/"s2_terms.json.gz", "rt", encoding="utf-8") as f:
        return [dict(state=x["state"], constants=x["constants"],
                     levels={(v, J): {k: np.array(t) for k, t in terms.items()} for v, J, terms in x["levels"]})
                for x in json.load(f)["terms"]]


def s2_fixed(state, constants):
    """The release's own second-order terms of the levels no committed input reproduces (s2_terms.json.gz)."""
    return next((x["levels"] for x in _s2_terms() if x["state"] == state and x["constants"] == constants), None)


# ---------------------------------------------------------------------------------------------------------------
# the installed release and the frozen inputs

def installed(role, filename):
    """A file of the installed release, checked against its manifest's SHA256 (read directly, so a candidate with
    the same ids cannot stand in for it)."""
    manifest = json.loads((INSTALLED/"manifest.json").read_text())
    entry = next(e for e in manifest["models"] if e["id"] == IDS[role])
    raw = (INSTALLED/entry["directory"]/filename).read_bytes()
    if hashlib.sha256(raw).hexdigest() != entry["files_sha256"][filename]:
        raise ValueError(f"{IDS[role]}/{filename} failed its SHA256 check")
    return json.loads(raw)


def released():
    """(core, continuation) of the installed release in native coordinates: the refits' starts and the later stages'
    bases. The core is as its fit left it for the later stages: the isotope-shift coefficients ride along in B's values
    (the isotopologue record holds the shift), so the covariance holds the potential's own coefficients."""
    core_spec = named.to_native(installed("core", "model.json"))
    core_spec["optimizer_parameters"]["potential"] = [p for p in core_spec["optimizer_parameters"]["potential"]
                                                      if p["name"] not in ("iso0", "iso1")]
    return core_spec, named.continued_to_native(installed("expanded", "model.json"), core_spec["named"])


def cross_validation():
    """Release 20261029's cross-validation (continuation and core folds, benchmark) with this form's nested hyperfine
    folds."""
    cv = U.cross_validation("20261029")
    with gzip.open(TRAINING/"hyperfine_folds.json.gz", "rt", encoding="utf-8") as f:
        cv["hyperfine_folds"] = json.load(f)
    cv["note"] = f"Cross-validation of release {RELEASE}, read by fit.uncertainty"
    return cv


def edges():
    """The rules for the levels beyond the tables and the quasi-bound width, with the evidence they were set on."""
    return json.loads((TRAINING/"edges.json").read_text())["edges"]


# ---------------------------------------------------------------------------------------------------------------
# stages

def core_rows():
    """The core training rows and the six CIPM isotopologue line centers."""
    rows = training.release("core")
    data = json.loads((training.TRAINING/"20261029"/"isotopologue_centers.json").read_text())["data"]
    extra = iso_fit.core_rows(data)
    return training.Rows({k: np.concatenate([rows.cols[k], extra.cols[k]]) for k in rows.cols},
                         list(rows.levels)+list(extra.levels))


def stage_core(start):
    """Potentials, then the B hyperfine functions, from the core `start` (native): the core spec with the refitted
    values and coefficients, and the fits' records."""
    with base.based_on(core=start, comb_floor=COMB_FLOOR):
        settings = core.Settings(isotope_shift=ISO_SHIFT, shared_asymptote=True)
        prob, result = core.fit(core_rows(), settings=settings)
    names = initial(replace(Config.from_spec(start["config"]), iso_shift=ISO_SHIFT))[1]
    values = {s: dict(start["values"][s]) for s in ("X", "B")}
    for (s, n), z in zip(names, result.theta):
        values[s][n] = float(z)
    config = Config.from_spec(start["config"])
    curves = {s: Curve(s, config, values[s]) for s in ("X", "B")}
    records = [r for r in hyperfine_rows_release() if r["inside"]]
    energies, means = core_hyperfine.moments(curves, values, [r["line"] for r in records])
    hp = core_hyperfine.Problem(energies, means, records, orders=CORE_ORDERS)
    old = start["hyperfine"]["coefficients"]
    hf = core_hyperfine.fit(hp, start=np.array([old[s][p][k] if k < len(old[s][p]) else 0. for s, p, k in hp.names]))
    spec = copy.deepcopy(start)
    spec["values"] = values
    for (s, p, k), z in zip(hp.names, hf.theta):
        spec["hyperfine"]["coefficients"][s][p][k] = float(z)
    return spec, dict(cost=float(result.cost), nfev=int(result.nfev), hyperfine_cost=float(hf.cost))


def stage_continuation(core_spec, start):
    """Tail terms on `core_spec` from the continuation `start` (native; its structure): the continuation spec with
    the refitted values, offset and fit record."""
    with base.based_on(core=core_spec, continuation=start, comb_floor=COMB_FLOOR):
        prob, result = continuation.fit(core=dict(config=core_spec["config"], values=core_spec["values"]))
        train = training.release("continuation")
        stats = continuation.statistics(train, train.value-result.pred)
    spec = copy.deepcopy(start)
    for k, z in result.parameters.items():
        s, n = k.split(":", 1)
        if s in spec["values"]:
            spec["values"][s][n] = float(z)
    spec["values"]["B"]["dDe"] = spec["values"]["X"]["dDe"]          # one shared shift
    spec["calibration"]["orsay4_offset_MHz"] = result.parameters["calibration:orsay4"]
    spec["fit"] = dict(spec["fit"], cost=float(result.cost), converged=bool(result.converged), nfev=int(result.nfev),
                       training=stats)
    return spec


def _bases(core_spec, cont_spec, candidate=None):
    return base.based_on(core=core_spec, continuation=cont_spec, comb_floor=COMB_FLOOR, s2_frozen=s2_frozen,
                         s2_fixed=s2_fixed, candidate=candidate)


def stage_hyperfine(core_spec, cont_spec, workers=12, cache=None):
    """The hyperfine beyond the core: FORM's fit at every scale (the lowest training cost kept), release 20261030's
    fit for the resonance, and the X model."""
    with _bases(core_spec, cont_spec):
        rows, params = hb.load_data()
        levels = hb.Levels()
        needed = hb.needed_levels(rows, params)
        scales = {}
        for scale in hb.FORM_SCALES:
            constants = hb.state_constants(False, hb.overrides_at(scale))
            model = hb.Model(levels, hb.fill_terms(needed, constants, workers, cache=cache), orders=FORM["orders"])
            b, theta, rec = hb.fit_starts(model, rows, params, FORM["keys"], FORM["starts"])
            scales[scale] = dict(constants=constants, b=np.asarray(b), theta=theta, core_cost=rec["core_cost"],
                                 out_cost=rec["out_cost"], n_out=rec["n_out"], converged=rec["converged"])
        chosen = min(scales, key=lambda s: scales[s]["out_cost"])
        constants = hb.state_constants(False, hb.FIT_OVERRIDES)
        model = hb.Model(levels, hb.fill_terms(needed, constants, workers, cache=cache))
        b, theta = hb.fit(model, rows, params, dict(u=30., g=30.), False, keys=list(hb.FIT_KEYS))
        resonance = dict(constants=constants, b=np.asarray(b), theta=theta, training=dict(model.training))
        xlv = xh.XLevels()
        xmodel = xh.fit(xlv, xh.KNOT, xh.ORDER)
    return dict(scales=scales, chosen=chosen, resonance=resonance, xmodel=xmodel)


def hyperfine_fold(core_spec, cont_spec, scale, block, start, workers=12, cache=None):
    """One fold of FORM's nested cross-validation (frozen in data/training/20261031/hyperfine_folds.json.gz): the fit
    at the 1g(1Pig) C5 `scale` with the v' `block` (lo, hi) held out, from the full fit's terms `start` and FORM's
    first start, and its held-out eqQ scored; the record as frozen."""
    with _bases(core_spec, cont_spec):
        rows, params = hb.load_data()
        levels = hb.Levels()
        constants = hb.state_constants(False, hb.overrides_at(scale))
        model = hb.Model(levels, hb.fill_terms(hb.needed_levels(rows, params), constants, workers, cache=cache),
                         orders=FORM["orders"])
        starts = [{k: start[k] for k in FORM["keys"]}]+list(FORM["starts"][:1])
        b, theta, _ = hb.fit_starts(model, rows, params, FORM["keys"], starts, holdout=(tuple(block),))
        b = np.asarray(b)
        held = [r for r in params if r["parameter"] == "eqQ" and block[0] <= r["v"] <= block[1]]
        scores = hb.score([], held, lambda v, J: model.constants(b, theta, v, J), levels)
    return dict(holdout=[list(block)], b=b.tolist(), theta={k: float(x) for k, x in theta.items()},
                training=dict(model.training),
                eqQ=[dict(key=k, rows=x["rows"]) for k, x in scores.items() if " eqQ " in k])


def _resonance(hf, levels=None, cache=None):
    r = hf["resonance"]
    return rz.Resonance(r["constants"], r["b"], r["theta"], levels=levels, cache=cache)


def stage_package(core_spec, cont_spec, hf, candidate, workers=12, cache=None):
    """Levels, tables, the resonance block, the spike levels and the references (the core's from the candidate);
    the expanded model's structure and values (native) and its hyperfine.json."""
    every = {f"{s}:{n}": z for s in "XB" for n, z in cont_spec["values"][s].items()}
    every["B:dDe"] = every["X:dDe"]
    with _bases(core_spec, cont_spec, candidate):
        from ..research import ResearchModel, _entry, _read
        structure, values = pk.structure_spec(every, cont_spec["structure"])
        structure["above_limit"] = "long grid"
        E = pk.levels(values, workers=workers)
        fit = hf["scales"][hf["chosen"]]
        constants, b, theta = fit["constants"], fit["b"], fit["theta"]
        blv = hb.Levels()
        xlv = xh.XLevels()
        xmodel = xh.fit(xlv, xh.KNOT, xh.ORDER)
        tables = pk.hyperfine_tables(E, b, theta, constants, xmodel, xlv, blv=blv, workers=workers, cache=cache,
                                     orders=FORM["orders"])
        resonance = pk.resonance_table(E, _resonance(hf, blv, cache))
        resonance["form"] = RESONANCE_NOTE
        found = pk.spike_levels(tables["B"], resonance)
        kept = dict(tables, B=[r for r in tables["B"] if (int(r[0]), int(r[1])) not in set(found)])
        own = reduction.OwnHyperfine(constants, b, theta, xmodel, xlv=xlv, levels=blv, orders=FORM["orders"],
                                     cache=cache)
        core_model = ResearchModel("balanced")
        refs = pk.references(E, dict(B=pk.rounded(kept["B"]), X=pk.rounded(kept["X"])), resonance, own,
                             _resonance(hf, cache=cache), core_model, _read(_entry("balanced"),
                                                                            "reference_predictions.json"))
        b_model = dict(result=B_MODEL_NOTE, theta_deg=theta, potential="the release's continuation", grid=hb.so.GRID,
                       curve_constants=constants, background=b.tolist(), orders=dict(FORM["orders"]),
                       c5_scale=hf["chosen"], training=dict(core_cost=fit["core_cost"], out_cost=fit["out_cost"],
                                                            n_out=fit["n_out"]), holdout=[])
        x_model = dict(x_model_record(xmodel), potential="the release's continuation", grid=dict(xh.GRID, nlev=xh.NLEV))
        hyperfine = pk.hyperfine_record(kept, resonance, b_model, x_model,
                                        core_model.coefficients["hyperfine"]["x_spin_spin"])
    hyperfine["B"]["not_predicted"] = [dict(v=v, J=J, reason=SPIKE_REASON) for v, J in found]
    hyperfine["B"]["spikes"] = dict(rule=SPIKE_RULE,
                                    rows=[r for r in pk.rounded(tables["B"]) if (r[0], r[1]) in set(found)])
    return dict(structure=structure, values=values, levels=E, hyperfine=hyperfine, references=refs)


# ---------------------------------------------------------------------------------------------------------------
# the candidate release

class Candidate:
    """The rebuild's release folder: the installed release's manifest and texts with the rebuilt files, rewritten
    as the stages finish."""

    def __init__(self, folder):
        self.folder = Path(folder)
        self.manifest = json.loads((INSTALLED/"manifest.json").read_text())
        self.files = {role: {} for role in IDS}
        self.content = {role: {} for role in IDS}

    def put(self, role, name, value, text=None):
        compact = name in ("hyperfine.json", "uncertainty.json")
        text = text or json.dumps(value, separators=(",", ":") if compact else None, indent=None if compact else 1,
                                  allow_nan=False)+"\n"
        entry = next(e for e in self.manifest["models"] if e["id"] == IDS[role])
        path = self.folder/entry["directory"]/name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text)
        self.files[role][name] = hashlib.sha256(text.encode()).hexdigest()
        self.content[role][name] = value
        manifest = copy.deepcopy(self.manifest)
        manifest["models"] = []
        for e in self.manifest["models"]:
            role_e = e["id"].rsplit("-", 1)[1]
            if not self.files[role_e]:
                continue
            e = dict(e, files_sha256=dict(self.files[role_e]))
            if "extends" in e:
                have = self.files["expanded"]
                e["extends"] = dict(e["extends"], files_sha256={k: have[k] for k in e["extends"]["files_sha256"]
                                                               if k in have})
            manifest["models"].append(e)
        (self.folder/"manifest.json").write_text(json.dumps(manifest, indent=1)+"\n")
        return self.files[role][name]


def core_record(spec):
    """core/model.json: the installed record with the refitted values and coefficients, in named coordinates."""
    out = named.to_named(spec)
    out.pop("named", None)
    template = installed("core", "model.json")
    out = dict(template, values=out["values"], hyperfine=dict(template["hyperfine"],
                                                              coefficients=spec["hyperfine"]["coefficients"]))
    out["optimizer_parameters"] = dict(
        potential=[dict(p, value=spec["values"][p["state"]][p["name"]])
                   for p in template["optimizer_parameters"]["potential"]],
        hyperfine=[dict(p, value=spec["hyperfine"]["coefficients"][p["state"]][p["name"]][p["power"]])
                   for p in template["optimizer_parameters"]["hyperfine"]])
    return out


def core_references(model_id):
    """The core's reference predictions: release 20261018's replay rows with the core's own values."""
    from ..research import ResearchModel, _data
    m = ResearchModel(model_id)
    refs = json.loads(_data("20261018").joinpath("balanced", "reference_predictions.json").read_text())
    for row in refs:
        value = m.frequency(row["line"], row["component"])
        if "ref_component" in row:
            value -= m.frequency(row["line"], row["ref_component"])
        row["predicted_MHz"] = value
    return refs


def expanded_record(core_spec, cont_spec, pkg, core_sha):
    """expanded/model.json: the installed record with the rebuilt structure, values (named), offset, fit and counts."""
    template = installed("expanded", "model.json")
    spec = dict(template, structure=pkg["structure"], values=pkg["values"], coordinates=None)
    out = named.continued_to_named(spec, named.named(core_spec))
    out["base"] = dict(template["base"], model_sha256=core_sha)
    out["calibration"] = dict(template["calibration"], orsay4_offset_MHz=cont_spec["calibration"]["orsay4_offset_MHz"])
    out["fit"] = dict(template["fit"], **{k: cont_spec["fit"][k] for k in ("cost", "converged", "nfev", "training")})
    out["reference_counts"] = pk.reference_counts(pkg["references"])
    detail = copy.deepcopy(template["inventory_detail"])
    detail["hyperfine_derived"]["B_background"] = len(pkg["hyperfine"]["B"]["model"]["background"])
    out["inventory_detail"] = detail
    out["coordinates"] = "named"
    return out


def stage_uncertainty(core_spec, cont_spec, cv, candidate, workers=12, cache=None):
    """uncertainty.build on the candidate, with the named parameters' covariance: (record, parts)."""
    with _bases(core_spec, cont_spec, candidate):
        record, parts = U.build(cv=cv, cache=cache, workers=workers)
        spec = U.core_spec()
        spec["values"]["X"]["dDe"] = cont_spec["values"]["X"]["dDe"]
        record["named"] = named.covariance(spec, record["statistical"]["parameters"], record["statistical"]["factor"],
                                           de_unit=cont_spec["structure"]["de_unit"])
    record["model"] = IDS["expanded"]
    return record, parts


def stage_isotopologues(core_spec, cont_spec, candidate):
    """The isotope shift and the isotopologues' record (not yet conditioned)."""
    with _bases(core_spec, cont_spec, candidate):
        # the 2008 formulae's comparison on the potentials of the default model when the isotope shift was fitted
        fitted = iso_fit.fit_shift(s06_parameters="research-20261021-i2spec")
        return fitted, iso_fit.extension_record(fitted, model_id=IDS["i2spec"], base_id=IDS["expanded"])


def stage_calibrate(core_spec, cont_spec, record, parts, cv, extension, candidate):
    """The calibration refit on the candidate (its uncertainty the build's record) and the X term near
    dissociation: the release's uncertainty.json."""
    table = json.loads(NEAR_DISSOCIATION_TABLE.read_text())
    near = calibrate.near_dissociation(calibrate.near_dissociation_records(parts, table))
    with _bases(core_spec, cont_spec, candidate):
        theta = calibrate.Theta(record["model_error"], record["hyperfine"], 1., extension["hyperfine"]["floor_MHz"])
        step4 = calibrate.refit(IDS["i2spec"], parts, cv, theta)
    fits, check = step4["fits"], step4["check"]
    unc = dict(record)
    unc["recipe"] = record["recipe"]+CALIBRATION_RECIPE+RECIPE_20261031
    unc["model_error"] = copy.deepcopy(fits["model_error"])
    unc["model_error"]["X beyond"]["near_dissociation"] = near
    unc["hyperfine"] = dict(record["hyperfine"], **{k: fits["hyperfine"][k] for k in ("core", "B beyond", "X beyond")})
    unc["checks"] = dict(record["checks"], step4=dict(
        calibrate_coverage=check["calibrate"], heldout_calibration_refit=check["cross_validation"]["refit"],
        isotopologue_floors=dict(fits["isotopologues"], kept=ISOTOPOLOGUE_FLOORS_KEPT), atlas=ATLAS_CHECK))
    out = edges()
    out["X above 108"]["sigma_D0_MHz"] = record["named"]["sigma"]["D0(X)"]*CM_TO_MHZ
    unc["edges"] = out
    return unc, step4


def stage_conditioned(core_spec, cont_spec, extension, candidate):
    """The isotopologue record conditioned on the CIPM lines, and the references of every isotopologue."""
    with _bases(core_spec, cont_spec, candidate):
        from ..research import ResearchModel
        ext = copy.deepcopy(extension)
        ext["isotope_shift"] = {k: v for k, v in ext["isotope_shift"].items() if not k.startswith("conditioned")}
        return iso_fit.condition(ext, {i: ResearchModel(IDS["i2spec"], i) for i in iso_fit.ISOS})


def stage_references(core_spec, cont_spec, candidate):
    with _bases(core_spec, cont_spec, candidate):
        return iso_fit.references(IDS["i2spec"])


def rebuild(folder=None, workers=12, cache=None, log=print):
    """Every stage from the repository; the candidate release is written to `folder` (a new temporary folder by
    default). Returns {stage: record} with "files" {"<role>/<file>": the rebuilt content of each release file (the
    core's model.json from the core refit)} and "folder". About four minutes with a warm second-order cache, ten
    with a cold one."""
    folder = Path(folder or tempfile.mkdtemp(prefix=f"i2spec-{RELEASE}-"))
    cand = Candidate(folder)
    clock = time.time()

    def done(stage):
        log(f"{stage} done: {time.time()-clock:.0f} s")
    core_spec, cont_start = released()
    core_refit, core_fit = stage_core(core_spec)
    done("core")
    # the later stages build on the released core, which the core refit returns (as release.rebuild)
    cont_spec = stage_continuation(core_spec, cont_start)
    done("continuation")
    hf = stage_hyperfine(core_spec, cont_spec, workers, cache)
    done("hyperfine")
    core_sha = cand.put("core", "model.json", None, (INSTALLED/"core"/"model.json").read_text())
    with _bases(core_spec, cont_spec, folder):
        cand.put("core", "reference_predictions.json", core_references(IDS["core"]))
    pkg = stage_package(core_spec, cont_spec, hf, folder, workers, cache)
    done("package")
    cand.put("expanded", "model.json", expanded_record(core_spec, cont_spec, pkg, core_sha))
    cand.put("expanded", "hyperfine.json", pkg["hyperfine"])
    cand.put("expanded", "reference_predictions.json", pkg["references"])
    cv = cross_validation()
    record, parts = stage_uncertainty(core_spec, cont_spec, cv, folder, workers, cache)
    done("uncertainty")
    cand.put("expanded", "uncertainty.json", record)
    fitted, extension = stage_isotopologues(core_spec, cont_spec, folder)
    cand.put("i2spec", "isotopologues.json", extension)
    cand.put("i2spec", "reference_predictions.json", [])
    cand.put("i2spec", "isotopologues.json", stage_conditioned(core_spec, cont_spec, extension, folder))
    done("isotopologues")
    unc, step4 = stage_calibrate(core_spec, cont_spec, record, parts, cv, extension, folder)
    done("calibrate")
    cand.put("expanded", "uncertainty.json", unc)
    cand.put("i2spec", "isotopologues.json", stage_conditioned(core_spec, cont_spec, extension, folder))
    cand.put("i2spec", "reference_predictions.json", stage_references(core_spec, cont_spec, folder))
    done("references")
    files = {f"{role}/{name}": v for role, files in cand.content.items() for name, v in files.items()}
    files["core/model.json"] = core_record(core_refit)
    return dict(folder=folder, core=core_fit, core_spec=core_refit, continuation=cont_spec, hyperfine=hf,
                package=pkg, uncertainty=record, isotope_shift=fitted, step4=step4, files=files)


# ---------------------------------------------------------------------------------------------------------------
# against the release

def _numbers(a, b, path=""):
    """[(path, |a - b|, |b|)] over the numbers of two JSON values of the same shape, and the paths where the shape or
    a non-number differs."""
    if isinstance(a, dict) and isinstance(b, dict):
        if set(a) != set(b):
            return [], [f"{path}: keys {sorted(set(a) ^ set(b))[:6]}"]
        out, bad = [], []
        for k in a:
            o, x = _numbers(a[k], b[k], f"{path}/{k}")
            out += o
            bad += x
        return out, bad
    if isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b):
            return [], [f"{path}: length {len(a)} != {len(b)}"]
        out, bad = [], []
        for i, (x, y) in enumerate(zip(a, b)):
            o, z = _numbers(x, y, f"{path}[{i}]")
            out += o
            bad += z
        return out, bad
    if isinstance(a, bool) or isinstance(b, bool) or not isinstance(a, (int, float)) or not isinstance(b, (int, float)):
        return [], ([] if a == b else [f"{path}: {str(a)[:60]!r} != {str(b)[:60]!r}"])
    return [(path, abs(float(a)-float(b)), abs(float(b)))], []


def compare(files, skip=("created_utc", "/fit/nfev")):
    """Per release file: the largest absolute and relative differences of the rebuilt numbers from the installed
    ones (with their paths), and every path where anything else differs (`skip`: paths not compared, the creation
    time and the refit's evaluation count)."""
    out = {}
    for key, value in files.items():
        role, name = key.split("/")
        ref = installed(role, name)
        diffs, bad = _numbers(json.loads(json.dumps(value)), ref)
        diffs = [d for d in diffs if not any(s in d[0] for s in skip)]
        bad = [x for x in bad if not any(s in x for s in skip)]
        worst = max(diffs, key=lambda d: d[1], default=("", 0., 0.))
        rel = max(diffs, key=lambda d: d[1]/max(d[2], 1e-300), default=("", 0., 1.))
        out[key] = dict(n=len(diffs), max_abs=worst[1], max_abs_at=worst[0], max_rel=rel[1]/max(rel[2], 1e-300),
                        max_rel_at=rel[0], other=bad)
    return out
