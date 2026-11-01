"""The uncertainty model refitted on the cross-validation and the calibrate sets together (the calibration refit,
from release 20261030), and the X term near dissociation (release 20261031).

The calibrate residuals are a release's on rows it never trained on: the line targets of the third blind registry
(src/i2spec/data/research/20261029/validation/pending_blind_v3, scored) of the sets whose role is calibrate
(data/training/roles_20261029), without the record rows, the rows of a compilation that a held set also gives, the
lines the training holds inside a CIPM mean (outcome.json held_out_in_training), the excluded line, and values printed
again in a later set; every line row's residual is recomputed with the release being built (observed - predicted),
a constant's residual stays the registry's. The form of sigma is the release's (sigma_stat + sigma_model + sigma_HF):

  sigma_model  each level a random effect of variance sigma_level^2 (the model's per-level shares, so lines sharing a
               level are correlated as the model assumes), each lab family a random offset of the absolute
               frequencies (measurement error: not part of sigma_model); the cross-validation part as
               uncertainty.fit_model. With FORMS "quadratic": a (v - v0)^2 term in the B model error beyond the core,
               chosen by an explicit profile (the likelihood is flat in it); B beyond does not fall below the core's
               per-level share at v' 44.
  sigma_HF     the pattern-scale error kappa is one per line, shared by its components, so a line counts once however
               many of its components were measured; the floor a is per component. With FORMS "a_growth": the B floor
               beyond the core grows as exp(g (v' - 44)) up to the last validated v'.
  isotopologues  the shift floor (a factor on the release's) and the hyperfine floor, on the calibrate rows of the
               other isotopologues (reported; the release keeps its own).

near_dissociation: for X v'' >= 85 an added per-level term exp(a + c (v - 85)) in quadrature, fitted by maximum
likelihood to the cross-validation records there with a given table's other terms, checked leaving one v'' out.

Only scored registries are read: no pending data file.
"""
from __future__ import annotations

import copy
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

from ..observations import Line
from . import hyperfine_rows, uncertainty as U
from .rows import CORE_V, EXCLUSION

ROOT = Path(__file__).resolve().parents[3]
V3 = ROOT/"src"/"i2spec"/"data"/"research"/"20261029"/"validation"/"pending_blind_v3"/"results.json"
V3_MODEL = "research-20261029-i2spec"           # the registry's primary forecast
ROLES = ROOT/"data"/"training"/"roles_20261029"
FORMS = ("a_growth", "quadratic")               # the release's two terms beyond the core (release 20261030)
NEAR_DISSOCIATION_V0 = 85


def roles(folder=ROLES):
    """(sets by id, final role of each set, lines inside a training CIPM mean per set, record (set, group) rows)."""
    r = json.loads((Path(folder)/"roles.json").read_text())
    o = json.loads((Path(folder)/"outcome.json").read_text())
    sets = {s["set"]: s for s in r["sets"]}
    final = {s: o["final_roles"].get(s, sets[s]["role"]) for s in sets}
    inside = defaultdict(set)
    for d in o["held_out_in_training"]["calibrate"]:
        inside[d["set"]] |= set(d["lines"]) if isinstance(d["lines"], list) else {"all"}
    records = {(d["set"], d["group"]) for d in o["row_records"]}
    return sets, final, inside, records


def calibrate_rows(results=V3, folder=ROLES):
    """The registry rows of calibrate sets that enter (module note), and a count of what is left out and why."""
    sets, final, inside, records = roles(folder)
    kept, why = [], Counter()
    for r in json.loads(Path(results).read_text())["rows"]:
        if final.get(r["set"]) != "calibrate":
            continue
        lines = {r.get("line"), r.get("ref_line")} - {None}
        reason = ("superseded" if r["status"] == "superseded" else
                  f"compilation row ({r['overlap_relation']})" if r["overlap_relation"] else
                  "diagnostic exclusion" if r["diagnostic_exclusion"] else
                  "record rows" if (r["set"], r.get("group")) in records else
                  "inside a training CIPM mean" if "all" in inside[r["set"]] or lines & inside[r["set"]] else
                  "P(13) 43-0" if EXCLUSION in lines else
                  "no residual" if r["residual_MHz"].get(V3_MODEL) is None else None)
        if reason:
            why[reason] += 1
        else:
            kept.append(r)
    # the same quantity with the same value in two sets is one measurement printed twice: the later printing stays
    # (the year in the set id; a tie keeps the first id)
    year = lambda sid: int("".join(ch for ch in sid if ch.isdigit())[:4] or 0)  # noqa: E731
    fields = ("kind", "line", "component", "ref_line", "ref_component", "quantity", "state", "v", "J", "I_from",
              "F_from", "I_to", "F_to")
    groups = defaultdict(list)
    for r in kept:
        groups[tuple(r.get(k) for k in fields)+(round(r["observed_MHz"], 4),)].append(r)
    once, printed = [], Counter()
    for g in groups.values():
        g = sorted(g, key=lambda r: (-year(r["set"]), r["id"]))
        once.append(g[0])
        for r in g[1:]:
            printed[f"{r['set']} (kept: {g[0]['set']})"] += 1
    why["printed again in a later set"] = sum(printed.values())
    return sorted(once, key=lambda r: r["id"]), dict(why, printed_again=dict(printed.most_common()))


# ---------------------------------------------------------------------------------------------------------------
# features: each row's sigma terms, from which sigma is recomputed for any refitted table

def _label(line):
    return str(line).split(" ", 1)[1]


def _level_features(m, terms):
    """[(state, v, J, sign, v_e, J_e, class, growth, B shift floor)] of one line's levels (the model's pieces)."""
    unc = m._uncertainty_model()
    out = []
    for (state, v, J), (sign, _) in terms["levels"].items():
        ve, Je = m._equivalent(state, v, J)
        growth, floor = 1., 0.
        if m._extension is not None:
            reach = unc.reach_J(state, ve)
            if reach is not None and Je > reach:
                growth = Je*(Je+1)/(reach*(reach+1))
            if m._iso is not None and state == "B":
                floor = m._shift_unc().floor
        out.append((state, v, J, sign, ve, Je, terms["cls"], growth, floor))
    return out


def _hf_features(terms):
    _, _, vu, vl = terms["kappa_floor"]
    return dict(vu=vu, vl=vl, offset=terms["offset"], component=terms["component"])


def features(model_id, rows, sets):
    """Each calibrate row with its residual recomputed with the model `model_id` and the pieces of its sigma (each
    level's coefficient and class, each component's offset), and its sigma from the model; (rows, refused)."""
    from ..research import ResearchModel
    models = {}

    def model(iso):
        if iso not in models:
            models[iso] = ResearchModel(model_id, iso)
        return models[iso]
    out, skipped = [], Counter()
    for r in rows:
        f = dict(id=r["id"], set=r["set"], lab="+".join(sets[r["set"]]["lab_families"]), kind=r["kind"],
                 status=r["status"], residual=float(r["residual_MHz"][V3_MODEL]), sigma_meas=float(r["sigma_MHz"]),
                 observed=float(r["observed_MHz"]), line=r.get("line"), ref_line=r.get("ref_line"),
                 quantity=r.get("quantity"))
        if r["kind"] not in ("absolute", "interline", "intraline"):
            f["sigma"] = None
            out.append(f)
            continue
        line = Line.parse(r["line"])
        ref = Line.parse(r["ref_line"]) if r["ref_line"] else None
        f["iso"] = line.isotopologue
        if ref is not None and ref.isotopologue != line.isotopologue:
            # no 1-sigma of an interval between isotopologues in the model: each line's own, composed
            one = model(line.isotopologue).uncertainty(_label(line), r["component"])
            two = model(ref.isotopologue).uncertainty(_label(ref), r["ref_component"])
            f.update(sigma=None, sigma_composed=math.hypot(one.sigma_MHz, two.sigma_MHz), cross_iso=True)
            out.append(f)
            continue
        m = model(line.isotopologue)
        try:
            pred = m.frequency(_label(line), r["component"])
            if r["kind"] == "intraline":
                pred -= m.frequency(_label(line), r["ref_component"])
            elif r["kind"] == "interline":
                pred -= m.frequency(_label(ref), r["ref_component"])
        except (ValueError, KeyError) as e:
            skipped[str(e)[:60]] += 1
            continue
        f["residual_v3"] = f["residual"]
        f["residual"] = float(r["observed_MHz"]-pred)
        args = (_label(line), r["component"]) + ((r["ref_component"],) if r["kind"] != "absolute" else ()) \
            + ((_label(ref),) if r["kind"] == "interline" else ())
        try:
            u = m.uncertainty(*args)
        except (ValueError, KeyError) as e:
            skipped[str(e)[:60]] += 1
            continue
        with m._lock:
            one = m._uncertainty_terms(_label(line), r["component"])
            two = None
            if r["kind"] == "intraline":
                two = m._uncertainty_terms(_label(line), r["ref_component"])
            elif r["kind"] == "interline":
                two = m._uncertainty_terms(_label(ref), r["ref_component"])
            f["levels"] = [_level_features(m, one)] + ([_level_features(m, two)] if r["kind"] == "interline" else [])
            f["hf"] = [_hf_features(one)] + ([_hf_features(two)] if two else [])
        f.update(sigma=u.sigma_MHz, stat=u.statistical_MHz, model=u.model_MHz, hyperfine=u.hyperfine_MHz,
                 level_status=u.status, line_class=u.line_class,
                 region="B beyond" if any(lv[0] == "B" and lv[4] > CORE_V["B"] for lvs in f["levels"] for lv in lvs)
                 else "X beyond" if any(lv[0] == "X" and lv[4] > CORE_V["X"] for lvs in f["levels"] for lv in lvs)
                 else "core")
        out.append(f)
    return out, dict(skipped)


# ---------------------------------------------------------------------------------------------------------------
# sigma of a row for a table (research_expanded's formulas)

def model_sigma(table, state, v, J, cls):
    """research_expanded's model_sigma, with the (v - v0)^2 term beyond the core (quadratic_per_v2 >= 0)."""
    m = table["model_error"]
    if v <= CORE_V[state]:
        return m[f"core {'atlas' if cls == 'beyond' else cls}"]["sigma_MHz"]/math.sqrt(2.)
    f = m[f"{state} beyond"]
    J = min(J, f["J_clip"])
    x = v-f["v0"]
    return float(np.exp(f["log_sigma0"]+f["slope_per_v"]*x+f.get("quadratic_per_v2", 0.)*x*x
                        + f["J2_coefficient"]*(J/100)**2))


def kappa(table, state, v):
    hf = table["hyperfine"]
    if v <= CORE_V[state]:
        return hf["core"]["kappa"]
    if state == "B":
        b = hf["B beyond"]
        return float(np.exp(b["log_kappa0"]+b["power"]*np.log1p((min(v, b["v_range"][1])-b["v_ref"])/b["v_scale"])))
    x = hf["X beyond"]
    return x["kappa_17"]*np.exp(x["slope_per_v"]*(min(v, x["validated_v"][1])-17))


def a_beyond(b, v):
    """The B floor beyond the core, with its growth exp(g (v - 44)) (a_growth_per_v >= 0)."""
    return b["a_MHz"]*math.exp(b.get("a_growth_per_v", 0.)*(min(v, b["v_range"][1])-44))


def hyperfine_terms(table, vu, vl):
    """research_expanded's hyperfine_terms with the floor's growth beyond the core: (a, kappa)."""
    hf = table["hyperfine"]
    kc = hf["core"]["kappa"]
    k2 = kappa(table, "B", vu)**2+kappa(table, "X", vl)**2-kc**2
    a = a_beyond(hf["B beyond"], vu) if vu > CORE_V["B"] else hf["core"]["a_MHz"]
    return a, math.sqrt(max(k2, kc**2))


class Theta:
    """An uncertainty table (model_error and hyperfine as in uncertainty.json) and the isotopologues' two floors (a
    factor on the release's shift floor per B level, the hyperfine floor)."""

    def __init__(self, model_error, hyperfine, shift_factor=1., iso_hf_floor=None):
        self.table = dict(model_error=model_error, hyperfine=hyperfine)
        self.shift_factor, self.iso_hf_floor = shift_factor, iso_hf_floor

    def level(self, state, ve, Je, cls, growth, floor):
        s = model_sigma(self.table, state, ve, Je, cls)*growth
        return math.hypot(s, floor*self.shift_factor) if floor else s

    def hf_terms(self, h, iso):
        a, k = hyperfine_terms(self.table, h["vu"], h["vl"])
        if iso != "127I2":
            a = math.hypot(a, self.iso_hf_floor)
        return a, k

    def parts(self, f):
        """(sigma_model^2 per level {key: (coefficient, sigma)}, sigma_HF) of a row."""
        net = {}
        for lvs, sign in zip(f["levels"], (1., -1.)):
            for state, v, J, s, ve, Je, cls, growth, floor in lvs:
                c, _ = net.get((state, v, J), (0., 0.))
                net[(state, v, J)] = (c+sign*s, self.level(state, ve, Je, cls, growth, floor))
        if f["kind"] == "absolute":
            net = {k: (1., s) for k, (_, s) in net.items()}
        h = f["hf"]
        hf = 0.
        if f["kind"] == "intraline":
            a, k = self.hf_terms(h[0], f["iso"])
            hf = math.sqrt(2*a*a+(k*(h[0]["offset"]-h[1]["offset"]))**2)
            net = {}
        elif f["kind"] == "absolute":
            if h[0]["component"]:
                a, k = self.hf_terms(h[0], f["iso"])
                hf = math.hypot(a, k*h[0]["offset"])
        elif h[0]["component"] and h[1]["component"] and (h[0]["vu"], h[0]["vl"]) == (h[1]["vu"], h[1]["vl"]):
            a, k = self.hf_terms(h[0], f["iso"])          # one band: the scale error is common
            hf = math.sqrt(2*a*a+(k*(h[0]["offset"]-h[1]["offset"]))**2)
        else:
            for x in h:
                if x["component"]:
                    a, k = self.hf_terms(x, f["iso"])
                    hf = math.hypot(hf, math.hypot(a, k*x["offset"]))
        return net, hf

    def sigma(self, f):
        net, hf = self.parts(f)
        return math.sqrt(f["stat"]**2+sum((c*s)**2 for c, s in net.values())+hf**2)


def check_features(rows, theta):
    """The sigma recomputed from each row's pieces at `theta` against the model's own: their largest relative
    difference."""
    d = [abs(theta.sigma(f)-f["sigma"])/f["sigma"] for f in rows if f.get("sigma")]
    return dict(n=len(d), max_relative=float(max(d)), over_1e6=int(sum(x > 1e-6 for x in d)))


# ---------------------------------------------------------------------------------------------------------------
# the fits

def _nll_dense(r, C):
    """Gaussian negative log-likelihood of residuals r with covariance C."""
    try:
        L = np.linalg.cholesky(C)
    except np.linalg.LinAlgError:
        return 1e30
    z = np.linalg.solve(L, r)
    return float(np.sum(np.log(np.diag(L)))+.5*z@z)


class HFGroup:
    """Intraline intervals of one line (or one fold's line): component floors a per component, the pattern-scale
    error kappa * Delta shared by the line's intervals, measurement sigma per row."""

    def __init__(self, r, sm, comp, ref):
        self.r, self.sm2 = np.asarray(r, float), np.asarray(sm, float)**2
        eq = lambda x, y: np.equal.outer(np.asarray(x), np.asarray(y)).astype(float)  # noqa: E731
        self.A = eq(comp, comp)-eq(comp, ref)-eq(ref, comp)+eq(ref, ref)

    def nll(self, a2, kD):
        C = np.diag(self.sm2)+a2*self.A+np.outer(kD, kD)
        return _nll_dense(self.r, C)


def hf_groups_cv(hf_rows, region):
    """The held-out intraline records of a region, one group per (fold, line); components identified by their
    offsets."""
    groups = defaultdict(list)
    for o in hf_rows:
        if o["region"] == region and o["offset"] is not None and o["ref_offset"] is not None:
            groups[(o["fold"], o["line"])].append(o)
    return [(HFGroup([o["predicted"]-o["observed"] for o in g], [o["sigma"] for o in g],
                     [round(o["offset"], 6) for o in g], [round(o["ref_offset"], 6) for o in g]),
             np.array([o["predicted"] for o in g]), g[0]["v"]) for g in groups.values()]


def hf_groups_cal(rows, region, iso="127I2"):
    """The calibrate intraline rows of a region, one group per line (all labs)."""
    groups = defaultdict(list)
    for f in rows:
        if f["kind"] == "intraline" and f.get("sigma") and f["region"] == region and f["iso"] == iso:
            groups[f["line"]].append(f)
    return [(HFGroup([f["residual"] for f in g], [f["sigma_meas"] for f in g], [f["hf"][0]["component"] for f in g],
                     [f["hf"][1]["component"] for f in g]), g) for g in groups.values()]


def kappa_B(t, v):
    return np.exp(t[1]+t[2]*np.log1p((np.asarray(v, float)-42)/2.))


def prepare(parts, cv):
    """The cross-validation records of uncertainty.build (its intermediate products `parts`, on `cv`) as the fits
    read them: the core folds' records grouped as fit_model groups them (residual, noise without the hf floor, is
    component) and the records beyond the core per state."""
    out = dict(held=parts["heldout"], hf_rows=parts["hf_heldout"], fits=parts["fits"],
               ratios=parts["covariance"].ratios, rows=parts["jacobian"].rows, cv=cv,
               chen={(int(r["v"]), int(r["J"])): r["observed"] for r in hyperfine_rows.release()[1]
                     if r["parameter"] == "eqQ"})
    precise_cover = U.fold_cover(out["rows"])
    core, beyond = {"precise": defaultdict(list), "atlas": defaultdict(list)}, {"B": [], "X": []}
    for o in out["held"]:
        lv = [(s, v, J) for s, v, J, _ in o["levels"]]
        sm = o["sigma_meas"]*(out["ratios"][o["category"]]["birge"] if o["category"] in U.NOISE_FROM_BIRGE else 1.)
        var = sm**2+o["sigma_stat"]**2
        far = [(s, v, J) for s, v, J in lv if v > CORE_V[s]]
        if not far:
            cls = U.line_class(precise_cover(o["fold"]), lv)
            core[cls][(o["fold"], tuple(lv))].append((o["residual"], var, float(o["kind"] == "component")))
        else:
            s, v, J = max(far, key=lambda x: (x[0] == "X", x[1]))
            beyond[s].append((o["residual"], var, (o["fold"], s, v, J//20), v, J, len(lv)-1,
                              float(o["kind"] == "component")))
    out["core_groups"] = {c: [np.array(g) for g in d.values()] for c, d in core.items()}
    out["beyond"] = beyond
    return out


def fit_hf(cv, rows, forms=FORMS):
    """sigma_HF: core (a, kappa) on the core folds' and the calibrate core intervals; B beyond (a, kappa0, power and
    the floor's growth) on the nested hyperfine folds, their held-out eqQ, the calibrate intervals beyond and the
    calibrate eqQ beyond (an eqQ' error dq scales the pattern by dq / |eqQ' - eqQ''|); X beyond as the release."""
    rel = cv["fits"]["hyperfine"]
    out = dict(rel)
    gcv, gcal = hf_groups_cv(cv["hf_rows"], "core"), hf_groups_cal(rows, "core")

    def offsets(g):                       # the model's Delta of each calibrate interval (signed, as the residual)
        return np.array([f["hf"][0]["offset"]-f["hf"][1]["offset"] for f in g])
    dcal = [offsets(g) for _, g in gcal]

    def nll_core(t):
        a2, k = np.exp(2*t[0]), np.exp(t[1])
        return (sum(G.nll(a2, k*np.abs(D)*np.sign(D)) for G, D, _ in gcv)
                + sum(G.nll(a2, k*D) for (G, _), D in zip(gcal, dcal)))
    (la, lk), nll = U._minimize(nll_core, [np.log(rel["core"]["a_MHz"]), np.log(max(rel["core"]["kappa"], 1e-8))])
    out["core"] = dict(a_MHz=float(np.exp(la)), kappa=float(np.exp(lk)), n=sum(len(G.r) for G, *_ in gcv),
                       n_calibrate=sum(len(G.r) for G, _ in gcal), groups_cv=len(gcv), groups_calibrate=len(gcal),
                       nll=nll)
    bcv, bcal = hf_groups_cv(cv["hf_rows"], "B beyond"), hf_groups_cal(rows, "B beyond")
    eq = []
    for i in range(len(cv["cv"]["hyperfine_folds"][str(U.HF_SCALES[1])]["folds"])):
        g = U.chosen_scale(cv["cv"], i)
        for score in cv["cv"]["hyperfine_folds"][str(g)]["folds"][i]["eqQ"]:
            eq += [(x["v"], x["resid"], x["sigma"]) for x in score["rows"]]
    eq = np.array(eq)
    chen = cv["chen"]
    dq = np.array([abs(U.EQQ_X-np.mean([q for (v, _), q in chen.items() if v == int(x)])) for x in eq[:, 0]])
    ceq = [f for f in rows if f["kind"] == "hf_constant" and f["quantity"] == "eqQ"
           and Line.parse(f["line"]).v_upper > CORE_V["B"]]
    cal_eq = np.array([[Line.parse(f["line"]).v_upper, f["residual"], f["sigma_meas"]] for f in ceq]) \
        if ceq else np.zeros((0, 3))
    cdq = np.array([abs(U.EQQ_X-f["observed"]) for f in ceq])
    bdcal = [offsets(g) for _, g in bcal]
    bv = [np.array([Line.parse(f["line"]).v_upper for f in g]) for _, g in bcal]
    growth = "a_growth" in forms
    vs = [v for *_, v in bcv]+[int(x) for x in eq[:, 0]]+[int(x) for v in bv for x in v]+[int(x) for x in cal_eq[:, 0]]
    vmax = max(vs)

    def nll_b(t):
        g = t[3] if growth else 0.
        a2 = lambda v: np.exp(2*(t[0]+g*(min(np.max(v), vmax)-44)))  # noqa: E731
        total = 0.
        for G, D, v in bcv:
            total += G.nll(a2(v), kappa_B(t, v)*D)
        for (G, _), D, v in zip(bcal, bdcal, bv):
            total += G.nll(a2(v), kappa_B(t, v)*D)
        total += U._gauss_nll(eq[:, 1], eq[:, 2]**2+(kappa_B(t, eq[:, 0])*dq)**2)
        if len(cal_eq):
            total += U._gauss_nll(cal_eq[:, 1], cal_eq[:, 2]**2+(kappa_B(t, cal_eq[:, 0])*cdq)**2)
        return total
    b = rel["B beyond"]
    t0 = [np.log(b["a_MHz"]), b["log_kappa0"], b["power"]]+([.1] if growth else [])
    bounds = [(-12., 5.), (-30., 0.), (-5., 10.)]+([(0., 2.)] if growth else [])
    t, nll = U._minimize(nll_b, t0, bounds)
    a, k0, pw = t[:3]
    out["B beyond"] = dict(b, a_MHz=float(np.exp(a)), log_kappa0=float(k0), power=float(pw),
                           a_growth_per_v=float(t[3]) if growth else 0.,
                           n_intervals=sum(len(G.r) for G, *_ in bcv), n_eqQ=len(eq),
                           n_calibrate_intervals=sum(len(G.r) for G, _ in bcal), n_calibrate_eqQ=len(cal_eq),
                           v_range=[int(min(vs)), int(max(vs))], nll=nll)
    out["X beyond"] = dict(rel["X beyond"], kappa_17=out["core"]["kappa"],
                           slope_per_v=float(np.log(rel["X beyond"]["kappa_51"]/out["core"]["kappa"])/(51-17)))
    return out


class LevelCov:
    """The calibrate rows' sigma_model as level random effects: A (rows x (level, class) keys) and each key's
    pieces, so that C_model = A diag(sigma_key^2) A^T for any table; lab offsets on the absolute rows."""

    def __init__(self, rows):
        self.rows = rows
        keys, entries = {}, []
        for i, f in enumerate(rows):
            for lvs, sign in zip(f["levels"], (1., -1.)):
                for state, v, J, s, ve, Je, cls, growth, floor in lvs:
                    k = keys.setdefault((state, v, J, cls), (len(keys), (state, ve, Je, cls, growth, floor)))[0]
                    entries.append((i, k, 1. if f["kind"] == "absolute" else sign*s))
        self.keys = sorted(keys.values())
        self.A = np.zeros((len(rows), len(keys)))
        for i, k, c in entries:
            self.A[i, k] += c
        if any(f["kind"] == "absolute" for f in rows):
            self.A = np.abs(self.A)*(np.array([f["kind"] == "absolute" for f in rows])[:, None]) + \
                self.A*(np.array([f["kind"] != "absolute" for f in rows])[:, None])
        labs = sorted({f["lab"] for f in rows})
        self.Z = np.array([[f["kind"] == "absolute" and f["lab"] == lab for lab in labs] for f in rows], float)
        self.r = np.array([f["residual"] for f in rows])
        self.stat2 = np.array([f["stat"]**2+f["sigma_meas"]**2 for f in rows])

    def C(self, th, tau_lab):
        s2 = np.array([th.level(*p)**2 for _, p in self.keys])
        hf2 = np.array([th.parts(f)[1]**2 for f in self.rows])
        return np.diag(self.stat2+hf2)+(self.A*s2)@self.A.T+tau_lab**2*self.Z@self.Z.T

    def nll(self, th, tau_lab):
        return _nll_dense(self.r, self.C(th, tau_lab))


def _cv_core_nll(cv, sig, hf_a):
    """uncertainty.fit_model's core likelihood on the core folds' records, per class, at sigma_class."""
    total = 0.
    for cls, groups in cv["core_groups"].items():
        tau2 = sig[cls]**2
        for g in groups:
            r, s2 = g[:, 0], g[:, 1]+g[:, 2]*hf_a**2
            a, b = np.sum(1/s2), np.sum(r/s2)
            total += .5*(np.sum(np.log(s2))+np.log(1+tau2*a)+np.sum(r*r/s2)-tau2*b*b/(1+tau2*a))
    return total


def fit_model_error(cv, rows, hf, base, forms=FORMS):
    """sigma_model: the core's two classes and a lab offset on the calibrate core rows; then B beyond (a, b, c and
    the quadratic term) with the calibrate rows beyond and their own lab offset; X beyond as the release's form on
    its folds and the calibrate rows. `base`: the release's Theta (its other terms)."""
    rel = cv["fits"]["model"]
    me = copy.deepcopy(rel)
    hft = dict(base.table["hyperfine"], **{k: hf[k] for k in ("core", "B beyond", "X beyond")})
    core_rows = [f for f in rows if f["kind"] in ("absolute", "interline") and f.get("sigma") and f["iso"] == "127I2"
                 and f["region"] == "core"]
    lc = LevelCov(core_rows)
    a_core = hf["core"]["a_MHz"]

    def theta(t):
        m = copy.deepcopy(me)
        m["core precise"]["sigma_MHz"], m["core atlas"]["sigma_MHz"] = float(np.exp(t[0])), float(np.exp(t[1]))
        return Theta(m, hft, 1., base.iso_hf_floor)

    def nll_core(t):
        sig = dict(precise=np.exp(t[0]), atlas=np.exp(t[1]))
        return _cv_core_nll(cv, sig, a_core)+lc.nll(theta(t), np.exp(t[2]))
    t0 = [np.log(rel["core precise"]["sigma_MHz"]), np.log(rel["core atlas"]["sigma_MHz"]), np.log(.05)]
    (lp, la, lt), nll = U._minimize(nll_core, t0)
    me["core precise"].update(sigma_MHz=float(np.exp(lp)), n_calibrate=int(sum(f["line_class"] == "precise"
                                                                                  for f in core_rows)))
    me["core atlas"].update(sigma_MHz=float(np.exp(la)), n_calibrate=int(sum(f["line_class"] != "precise"
                                                                                for f in core_rows)))
    me["lab offset (core absolute, not in predictions)"] = dict(sigma_MHz=float(np.exp(lt)), nll=nll,
                                                                n_labs=int(lc.Z.shape[1]))
    core_level_var = me["core atlas"]["sigma_MHz"]**2/2
    for s in ("B", "X"):
        brows = cv["beyond"][s]
        groups = defaultdict(list)
        for r, var, key, v, J, others, comp in brows:
            groups[key].append((r, var+others*core_level_var+comp*a_core**2, v, J))
        groups = [(v, np.array(g)) for (_, _, v, _), g in groups.items()]
        cal = [f for f in rows if f["kind"] in ("absolute", "interline") and f.get("sigma") and f["iso"] == "127I2"
               and f["region"] == f"{s} beyond"]
        J_top = int(max([J for *_, J, _, _ in brows]+[lv[5] for f in cal for lvs in f["levels"] for lv in lvs
                                                         if lv[0] == s and lv[4] > CORE_V[s]]))
        lcb = LevelCov(cal) if cal else None
        quad = "quadratic" in forms and s == "B"
        nq = 4 if quad else 3

        def theta_b(t, s=s, J_top=J_top, quad=quad):
            m = copy.deepcopy(me)
            m[f"{s} beyond"].update(log_sigma0=float(t[0]), slope_per_v=float(t[1]), J2_coefficient=float(t[2]),
                                    J_clip=J_top, quadratic_per_v2=float(t[3]) if quad else 0.)
            return Theta(m, hft, 1., base.iso_hf_floor)

        def nll_b(t, groups=groups, s=s, J_top=J_top, lcb=lcb, quad=quad, nq=nq):
            a, b, c = t[:3]
            d = t[3] if quad else 0.
            total = 0.
            for v, x in groups:
                r, s2, J = x[:, 0], x[:, 1], x[:, 3]
                tau2 = np.exp(2*(a+b*(v-U.V0[s])+d*(v-U.V0[s])**2+c*(min(np.mean(J), J_top)/100)**2))
                A, B = np.sum(1/s2), np.sum(r/s2)
                total += .5*(np.sum(np.log(s2))+np.log(1+tau2*A)+np.sum(r*r/s2)-tau2*B*B/(1+tau2*A))
            if lcb is not None:
                total += lcb.nll(theta_b(t), np.exp(t[nq]))
            return total
        m0 = rel[f"{s} beyond"]
        t0 = [max(m0["log_sigma0"], -1.5), m0["slope_per_v"], m0["J2_coefficient"]]+([0.] if quad else []) \
            + ([np.log(5.)] if lcb else [])
        # model error may not fall below the core's per-level share at the core boundary (B), nor toward dissociation
        # or with J (X: Martin's noise hides it beyond v''25)
        floor_b = float(np.log(me["core atlas"]["sigma_MHz"]/np.sqrt(2.)))
        bounds = ([(-10., 10.), (0., 1.), (0., 5.)] if s == "X" else [(floor_b, 10.), (-1., 1.), (0., 10.)]) + \
            ([(0., .1)] if quad else [])+([(-10., 10.)] if lcb else [])
        t, nll = U._minimize(nll_b, t0, bounds)
        if quad:
            # the likelihood is nearly flat in d and the optimizer settles anywhere along it: an explicit profile,
            # each d from the release's curve, the best kept
            best = None
            for dfix in np.arange(0., .0105, .001):
                fixed = [b if i != 3 else (dfix, dfix) for i, b in enumerate(bounds)]
                t1, n1 = U._minimize(nll_b, list(t0[:3])+[dfix]+list(t0[4:]), fixed)
                if best is None or n1 < best[1]:
                    best = (t1, n1)
            t, nll = best
        vs = [v for *_, v, J, _, _ in brows]+[lv[4] for f in cal for lvs in f["levels"] for lv in lvs if lv[0] == s
                                               and lv[4] > CORE_V[s]]
        me[f"{s} beyond"] = dict(m0, log_sigma0=float(t[0]), slope_per_v=float(t[1]), J2_coefficient=float(t[2]),
                                 quadratic_per_v2=float(t[3]) if quad else 0., J_clip=J_top, n_calibrate=len(cal),
                                 v_range=[int(min(vs)), int(max(vs))], nll=nll)
        if lcb:
            me[f"lab offset ({s} beyond absolute, not in predictions)"] = dict(sigma_MHz=float(np.exp(t[nq])),
                                                                              n_labs=int(lcb.Z.shape[1]))
    return me


def fit_iso(rows, me, hf, base):
    """The isotopologues' two floors: the shift floor (a factor on the release's, per B level) on the calibrate
    intervals within one isotopologue, and the hyperfine floor on their intraline intervals (one group per line)."""
    hft = dict(base.table["hyperfine"], **{k: hf[k] for k in ("core", "B beyond", "X beyond")})
    pos = [f for f in rows if f["kind"] in ("absolute", "interline") and f.get("sigma") and f["iso"] != "127I2"]
    lc = LevelCov(pos)
    hfg = [(G, g) for iso in ("129I2", "127I129I") for G, g in hf_groups_cal(rows, "core", iso)]

    def nll(t):
        th = Theta(me, hft, float(np.exp(t[0])), float(np.exp(t[1])))
        total = lc.nll(th, np.exp(t[2]))
        for G, g in hfg:
            a, k = th.hf_terms(g[0]["hf"][0], g[0]["iso"])
            D = np.array([f["hf"][0]["offset"]-f["hf"][1]["offset"] for f in g])
            total += G.nll(a*a, k*D)
        return total
    (ls, lf, lt), val = U._minimize(nll, [0., np.log(base.iso_hf_floor), np.log(.05)])
    return dict(shift_floor_factor=float(np.exp(ls)), hyperfine_floor_MHz=float(np.exp(lf)),
                lab_offset_MHz=float(np.exp(lt)), n_intervals=len(pos), n_intraline=sum(len(G.r) for G, _ in hfg),
                release_hyperfine_floor_MHz=base.iso_hf_floor, nll=val)


def check(rows, cv, base, fits):
    """Coverage of the calibrate rows (per isotopologue, kind and region; release table against the refit, with and
    without the lab offset) and of the cross-validation records (uncertainty.calibration) under both tables."""
    table = copy.deepcopy(base.table)
    table["model_error"] = fits["model_error"]
    table["hyperfine"] = dict(table["hyperfine"], **{k: fits["hyperfine"][k] for k in ("core", "B beyond",
                                                                                        "X beyond")})
    iso = fits["isotopologues"]
    new = Theta(table["model_error"], table["hyperfine"], iso["shift_floor_factor"], iso["hyperfine_floor_MHz"])
    lab = {"core": fits["model_error"]["lab offset (core absolute, not in predictions)"]["sigma_MHz"],
           "B beyond": fits["model_error"].get("lab offset (B beyond absolute, not in predictions)",
                                               {}).get("sigma_MHz", 0.)}
    strata = defaultdict(lambda: defaultdict(list))
    for f in rows:
        if not f.get("sigma"):
            continue
        v = max(lv[4] for lvs in f["levels"] for lv in lvs if lv[0] == "B")
        region = f["region"] if f["region"] != "B beyond" else f"B beyond v' {v//10*10}-{v//10*10+9}"
        key = (f["iso"] if f["iso"] != "127I2" else "", f["kind"], region)
        tau = lab.get(f["region"], 0.) if f["kind"] == "absolute" else 0.
        for name, s_pred in (("release", f["sigma"]), ("refit", new.sigma(f))):
            z = f["residual"]/math.hypot(s_pred, f["sigma_meas"])
            strata[key][name].append(z)
            if tau:
                strata[key][name+" +lab"].append(f["residual"]/math.sqrt(s_pred**2+f["sigma_meas"]**2+tau**2))
    report = {}
    for key in sorted(strata):
        d = strata[key]
        stats = {k: (float(np.sqrt(np.mean(np.square(z)))), float(np.mean(np.abs(z) <= 1)),
                     float(np.mean(np.abs(z) <= 2))) for k, z in d.items()}
        report[" | ".join(x for x in key if x)] = dict(n=len(d["release"]), **stats)
    for name, sel in (("all", lambda k: True), ("127I2 core", lambda k: k[0] == "" and k[2] == "core")):
        for which in ("release", "refit"):
            z = np.concatenate([strata[k][which] for k in strata if sel(k)])
            report[f"total {name} | {which}"] = dict(n=len(z), within_1=float(np.mean(abs(z) <= 1)),
                                                    within_2=float(np.mean(abs(z) <= 2)))
    cal = {which: U.calibration(cv["held"], dict(model=me, hyperfine=hf), cv["ratios"], cv["rows"])
           for which, me, hf in (("release", cv["fits"]["model"], cv["fits"]["hyperfine"]),
                                 ("refit", fits["model_error"], fits["hyperfine"]))}
    return dict(calibrate=report, cross_validation=cal)


def refit(model_id, parts, cv, base, results=V3, folder=ROLES, forms=FORMS):
    """The calibration refit on the release `model_id` (its uncertainty table the build's, `base`): the calibrate
    rows and their features, the hyperfine, model-error and isotopologue fits, and the coverage check. `parts`:
    uncertainty.build's intermediate products on the cross-validation `cv`."""
    kept, left_out = calibrate_rows(results, folder)
    rows, skipped = features(model_id, kept, roles(folder)[0])
    data = prepare(parts, cv)
    hf = fit_hf(data, rows, forms)
    me = fit_model_error(data, rows, hf, base, forms)
    iso = fit_iso(rows, me, hf, base)
    fits = dict(hyperfine=hf, model_error=me, isotopologues=iso)
    return dict(fits=fits, check=check(rows, data, base, fits), left_out=left_out, skipped=skipped,
                features=check_features(rows, base), rows=rows)


# ---------------------------------------------------------------------------------------------------------------
# X v'' >= 85: the term near dissociation

def near_dissociation_records(parts, table, v0=NEAR_DISSOCIATION_V0):
    """The cross-validation records whose highest level beyond the core is X v'' >= v0, with the sigma of the table
    `table` (uncertainty.calibration's variance) and the measurement sigma (Birge-scaled where the build scales it)."""
    held, ratios = parts["heldout"], parts["covariance"].ratios
    me, hf = table["model_error"], table["hyperfine"]
    out = []
    for o in held:
        lv = [(s, v, J) for s, v, J, _ in o["levels"]]
        beyond = [(s, v, J) for s, v, J in lv if v > CORE_V[s]]
        if not beyond:
            continue
        s, v, J = max(beyond, key=lambda x: (x[0] == "X", x[1]))
        if s != "X" or v < v0:
            continue
        sm = o["sigma_meas"]*(ratios[o["category"]]["birge"] if o["category"] in U.NOISE_FROM_BIRGE else 1.)
        pred2 = o["sigma_stat"]**2+sum(model_sigma(dict(model_error=me), st, w, K, "atlas")**2 for st, w, K in lv)
        if o["kind"] == "component":
            pred2 += hf["core"]["a_MHz"]**2
        out.append(dict(v=int(v), residual=float(o["residual"]), sigma_meas=float(sm), sigma=math.sqrt(pred2)))
    return out


def _near_nll(t, rows, v0):
    a, c = t
    return sum(.5*(math.log(var)+r["residual"]**2/var) for r in rows
               for var in (r["sigma"]**2+r["sigma_meas"]**2+math.exp(2*(a+c*(r["v"]-v0))),))


def _near_fit(rows, v0):
    from scipy.optimize import minimize
    best = None
    for a0 in (2., 4., 6.):
        res = minimize(_near_nll, [a0, .1], args=(rows, v0), bounds=[(-5., 9.), (0., 1.)])
        if best is None or res.fun < best.fun:
            best = res
    return best.x


def near_dissociation(records, v0=NEAR_DISSOCIATION_V0):
    """The X term near dissociation fitted to `records` (near_dissociation_records), with the coverage of the
    records within 1 and 2 sigma without it and with it, each v'' left out in turn."""
    t = _near_fit(records, v0)

    def z(r, t):
        a, c = t
        return r["residual"]/math.sqrt(r["sigma"]**2+r["sigma_meas"]**2+math.exp(2*(a+c*(r["v"]-v0))))
    zl, z0 = [], []
    for v in sorted({r["v"] for r in records}):
        tv = _near_fit([r for r in records if r["v"] != v], v0)
        for r in records:
            if r["v"] == v:
                zl.append(z(r, tv))
                z0.append(r["residual"]/math.hypot(r["sigma"], r["sigma_meas"]))
    zl, z0 = np.abs(zl), np.abs(z0)
    return dict(v0=v0, log_sigma_add=float(t[0]), slope_per_v=float(t[1]), n_cv=len(records),
                form="added in quadrature to an X level's share for v'' >= v0: exp(log_sigma_add + slope_per_v "
                     "(v - v0))",
                coverage_cv_released=[float(np.mean(z0 <= 1)), float(np.mean(z0 <= 2))],
                coverage_cv_leave_one_v_out=[float(np.mean(zl <= 1)), float(np.mean(zl <= 2))])
