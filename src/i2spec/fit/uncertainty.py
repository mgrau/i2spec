"""The uncertainty model: built from the final fit and its cross-validation only.

  sigma_center^2 = sigma_stat^2 + sigma_model^2          (hyperfine-free line center)
  sigma_component^2 = sigma_center^2 + sigma_HF^2        (absolute hyperfine component)
  sigma_interval = sigma_HF of the interval              (two components of one line share the center)

  sigma_stat   linearized covariance at the final continuation fit, propagated per line: g^T Cov g, g = d(line)/
               d(parameter). The parameters are the free tail terms (B:dDe = X:dDe; 37 from release 20261029), the
               core potential coefficients (46 with the B limit tied; frozen in the fit, refitted jointly in the
               linearization, so their covariance and cross
               terms come from the same rows), the 7 inherited dispersion coefficients (as log factors, prior 0.1),
               the Orsay IV offset and the 4 atlas calibration terms (priors: their own fits). Cov is the sandwich of
               the weighted estimator actually used, each data category scaled by its Birge ratio (at least 1).
  sigma_model  a smooth function of (state, v) and J, fitted by marginal likelihood to held-out residuals (the core
               fit's 7 folds, the continuation's 4 v folds), each residual's measurement sigma and its fold's
               sigma_stat removed in quadrature.
  sigma_HF     sqrt(a^2 + (kappa |offset|)^2) per region (core, B beyond v'43, X beyond v''17), fitted the same way to
               held-out intraline intervals (core folds, nested hyperfine cross-validation).

The cross-validation outputs the release used are frozen in data/training/<release>/cross_validation.json.gz
(`cross_validation`); `continuation_fold` and `hyperfine_fold` refit a fold with the ported fits.
"""
from __future__ import annotations

import datetime
import gzip
import json

import numpy as np
import scipy.sparse as sp
from scipy.optimize import minimize

from ..constants import MHZ_PER_CM
from ..observations import Line, load_all
from . import continuation, hyperfine_beyond as hb, hyperfine_rows, rows as training
from .core import release_spec as core_spec
from .levels import Levels
from .reduction import EMISSION_SETS

DISPERSION_PRIOR = .1        # relative 1-sigma of each inherited dispersion coefficient (no source quotes one)
FD_STEP = 1e-5               # relative central-difference step for the frozen coefficients
RCOND = 1e-8
J_MARGIN = 15                # a level is measured if a row reaches its v within this many J
CALIBRATED = ("atlas_sr", "atlas_apo")
PRECISE = ("comb", "doppler", "emission")
NOISE_FROM_BIRGE = ("atlas_sr", "atlas_apo", "orsay1", "orsay4", "martin", "martin_dissoc")
CORE_V = dict(X=17, B=43)
V0 = dict(X=18, B=44)        # first level beyond the core
EQQ_X = -2452.2916           # MHz, X v''=0 eqQ (the core's fixed X constants): eqQ' - eqQ'' sets the pattern's scale
HF_SCALES = (.6, .65, .7, .75, .8, .85, .9)
HF_STATUS = dict(
    B=[[43, "core"], [62, "validated: held-out intervals"], [75, "validated: held-out constants (Chen 2004, Pique 1986)"],
       [76, "extrapolated"], [78, "resonance with 1g(1Pig) (Pique 1986), not validated"], [80, "extrapolated"]],
    X=[[17, "core"], [47, "interpolated: no hyperfine data"], [54, "validated: emission components"],
       [108, "extrapolated"]])


def cross_validation(name=training.RELEASE):
    with gzip.open(training.TRAINING/name/"cross_validation.json.gz", "rt", encoding="utf-8") as f:
        return json.load(f)


# ---------------------------------------------------------------------------------------------------------------
# the joint Jacobian at the final fit

def core_names(spec):
    """The 47 core potential coefficients the core fit's optimizer fitted, and the 7 inherited dispersion log factors."""
    fitted = [(p["state"], p["name"]) for p in spec["optimizer_parameters"]["potential"]]
    inherited = [(s, f"logC{n}") for s in ("X", "B") for n in sorted(map(int, spec["inherited_potential_coefficients"][s]))]
    return fitted, inherited


def joint(curve, structure):
    """Patch `curve.derivatives` to differentiate the frozen core coefficients too (central differences of the
    whole tailed curve, whose Taylor continuations follow the core); the tails keep their analytic columns."""
    analytic, cache = curve.derivatives, {}

    def step(name):
        return FD_STEP*max(1., abs(curve.frozen.get(name, 0.)))

    def pair(name):
        if name not in cache:
            out = []
            for sign in (1, -1):
                frozen = dict(curve.frozen)
                frozen[name] = frozen.get(name, 0.)+sign*step(name)
                out.append(continuation.TailCurve(curve.state, curve.config, frozen, structure))
            cache[name] = out
        return cache[name]

    def derivatives(values, R, ratio, names):
        dV, df = np.zeros((len(R), len(names))), np.zeros((len(R), len(names)))
        tail = [n for n in names if n in curve.tail_names]
        if tail:
            a, b = analytic(values, R, ratio, tail)
            for j, n in enumerate(tail):
                dV[:, names.index(n)], df[:, names.index(n)] = a[:, j], b[:, j]
        for j, n in enumerate(names):
            if n not in curve.tail_names:
                plus, minus = pair(n)
                vp, fp = plus.curves(values, R, ratio)
                vm, fm = minus.curves(values, R, ratio)
                dV[:, j], df[:, j] = (vp-vm)/(2*step(n)), (fp-fm)/(2*step(n))
        return dV, df

    curve.derivatives = derivatives


class Jacobian:
    """Levels and their joint derivatives at the release's continuation (or `values`), and the training rows'
    functional M, offsets N, weights w and the fit's predictions."""

    def __init__(self, rows=None, extra_keys=(), values=None):
        spec = continuation.release_spec()
        structure = spec["structure"]
        self.rows = rows if rows is not None else training.release("continuation")
        curves = continuation.build_curves(structure)
        values = values or spec["values"]
        every = [(s, n) for s, c in curves.items() for n in c.tail_names]
        for s, n in every:
            curves[s].values[n] = values[s][n]
        prob = continuation.Problem(self.rows, curves, every, tie={("B", "dDe"): ("X", "dDe")})
        self.problem = prob
        x = continuation.release_start(prob, spec)
        self.pred, _ = prob.predict(x)
        for c in curves.values():
            joint(c, structure)
        fitted, inherited = core_names(core_spec())
        self.names = every+fitted+inherited
        self.theta = np.array([curves[s].values.get(n, curves[s].frozen.get(n, 0.)) for s, n in self.names])
        self.tie = {"B:dDe": "X:dDe"}
        keys = sorted(self.rows.level_set() | set(extra_keys))
        engine = Levels(keys, joins={s: c.joins for s, c in curves.items()}, box_padding=.6, refine_x=True)
        if core_spec().get("asymptote") == "tied":
            # De(B) = De(X) + the atomic splitting - Te is no free parameter: its column folds into theirs
            de_b = ("B", "logDe")
            theta = np.r_[self.theta, curves["B"].values.get("logDe", curves["B"].frozen.get("logDe"))]
            self.E, D = engine.evaluate(curves, theta, self.names+[de_b], jac=True)
            dx, db = np.exp(theta[self.names.index(("X", "logDe"))]), np.exp(theta[-1])
            D[:, self.names.index(("X", "logDe"))] += D[:, -1]*dx/db
            D[:, self.names.index(("B", "Te"))] -= D[:, -1]/db
            self.D = D[:, :-1]
        else:
            self.E, self.D = engine.evaluate(curves, self.theta, self.names, jac=True)
        self.keys, self.index = engine.keys, engine.index
        rr, cc, vv = [], [], []
        for i, lv in enumerate(self.rows.levels):
            for key, sign in lv:
                rr.append(i)
                cc.append(engine.index[key])
                vv.append(sign)
        self.M = sp.csr_matrix((vv, (rr, cc)), shape=(len(self.rows), len(keys)))
        self.N, self.offsets, self.w = prob.N, prob.offsets, prob.w
        self.x_offsets = x[len(prob.free):]
        self.category, self.value, self.sigma = self.rows.category, self.rows.value, self.rows.sigma


# ---------------------------------------------------------------------------------------------------------------
# covariance

def parameter_space(jac):
    """Level Jacobian columns on the free parameters (B:dDe folded into X:dDe), and their names."""
    names = [f"{s}:{n}" for s, n in jac.names]
    keep = [n for n in names if n not in jac.tie]
    T = np.zeros((len(names), len(keep)))
    for i, n in enumerate(names):
        T[i, keep.index(jac.tie.get(n, n))] = 1.
    return jac.D@T, keep


def row_jacobian(jac, D):
    """MHz per unit of each free parameter, the Orsay IV offset and the atlas calibrations (offset, slope)."""
    G = np.asarray(jac.M@(D*MHZ_PER_CM))
    extra, names = [], []
    for j, c in enumerate(jac.offsets):
        extra.append(jac.N[:, j])
        names.append(f"calibration:{c}")
    x = (jac.value/MHZ_PER_CM-17000.)/1000.
    for c in CALIBRATED:
        m = (jac.category == c).astype(float)
        extra += [m, m*x]
        names += [f"calibration:{c}:offset", f"calibration:{c}:slope"]
    return np.column_stack([G]+extra), names


def priors(names, spec):
    """1-sigma priors: inherited dispersion log factors and the atlas calibrations (their own fits' sigmas)."""
    out = {}
    for n in names:
        if ":logC" in n:
            out[n] = DISPERSION_PRIOR
        elif n.startswith("calibration:atlas"):
            _, c, part = n.split(":")
            cal = spec["calibration"][c]
            out[n] = cal["offset_unc_MHz"] if part == "offset" else cal["slope_unc"]
    return out


def birge(category, w, r):
    """Per category: n, f_c (the fit's count down-weighting), rms of r/sigma' (sigma' with the fit's floor)."""
    out = {}
    for c in sorted(set(category)):
        m = category == c
        f = min(1., 1000./m.sum())
        zw = r[m]*w[m]/np.sqrt(f)
        out[c] = dict(n=int(m.sum()), f=f, birge=float(np.sqrt(np.mean(zw**2))))
    return out


def sandwich(G, w, category, ratios, prior_sigma, mask=None, rcond=None, constraints=()):
    """Factor L (Cov = L L^T) of the weighted Gauss-Newton estimator's covariance A^-1 B A^-1, A = G'W^2G + P,
    B = sum_c f_c max(1, B_c^2) G_c'W_c^2G_c + P (P the prior precisions, f_c the fit's count down-weighting, B_c the
    category's Birge ratio). From the SVD of the column-scaled weighted Jacobian (priors as extra rows), never the
    normal equations: with A = V S^2 V^T, Cov = V S^-1 K S^-1 V^T, K = U^T diag(beta) U. Directions with singular value
    below rcond * the largest are dropped (counted): low-order core Chebyshev combinations no level depends on.
    constraints: further prior rows, each a gradient over the parameters divided by its sigma (linear combinations)."""
    rcond = RCOND if rcond is None else rcond
    m = np.ones(len(w), bool) if mask is None else mask
    Gw = G[m]*w[m, None]
    beta = np.array([ratios[c]["f"]*max(1., ratios[c]["birge"]**2) for c in category[m]])
    rows = [Gw]
    for j, sig in prior_sigma.items():
        e = np.zeros(G.shape[1])
        e[j] = 1/sig
        rows.append(e[None, :])
    rows += [np.asarray(c, float)[None, :] for c in constraints]
    aug = np.vstack(rows)
    beta = np.r_[beta, np.ones(len(prior_sigma)+len(constraints))]
    scale = np.linalg.norm(aug, axis=0)
    scale[scale == 0] = 1.
    U, S, Vt = np.linalg.svd(aug/scale, full_matrices=False)
    keep = S > rcond*S[0]
    Uk = U[:, keep]
    K = (Uk*beta[:, None]).T@Uk
    C = np.linalg.cholesky((K+K.T)/2)
    L = (Vt[keep].T/S[keep])@C/scale[:, None]
    leverage = (Uk[:len(Gw)]@C)/w[m, None]
    return L, leverage, dict(dropped=int((~keep).sum()), smallest_kept=float(S[keep][-1]/S[0]))


def dissociation_constraints(names, spec, settings=None):
    """The core fit's two dissociation priors as rows over the parameters: De(X) = its center +- dissociation_prior,
    and Te + De(B) - De(X) = the atomic splitting +- 0.0065 cm-1 (core.Settings). The release's covariance left them
    out, so its two limits could move apart although the fit held their difference."""
    from .core import ATOMIC_GAP_SIGMA, Settings
    settings = settings or Settings()
    dx, db = np.exp(spec["values"]["X"]["logDe"]), np.exp(spec["values"]["B"]["logDe"])
    ix, ib, it = (names.index(n) for n in ("X:logDe", "B:logDe", "B:Te"))
    rows = []
    if settings.dissociation_prior:
        r = np.zeros(len(names))
        r[ix] = dx/settings.dissociation_prior
        rows.append(r)
        r = np.zeros(len(names))
        r[ix], r[ib], r[it] = -dx/ATOMIC_GAP_SIGMA, db/ATOMIC_GAP_SIGMA, 1/ATOMIC_GAP_SIGMA
        rows.append(r)
    return rows


class Covariance:
    """The statistical part: G, the parameter names, priors and Birge ratios of the final fit, and its factor L.
    With `ties` the core fit's dissociation priors enter too (dissociation_constraints)."""

    def __init__(self, jac, ties=False):
        spec = core_spec()
        self.jac = jac
        self.D, self.free = parameter_space(jac)
        self.G, extra = row_jacobian(jac, self.D)
        self.names = self.free+extra
        self.ratios = birge(jac.category, jac.w, jac.pred-jac.value)
        self.prior = {self.names.index(n): s for n, s in priors(self.names, spec).items()}
        self.constraints = dissociation_constraints(self.names, spec) if ties else []
        self.L, self.leverage, self.info = sandwich(self.G, jac.w, jac.category, self.ratios, self.prior,
                                                    constraints=self.constraints)

    def masked(self, mask):
        return sandwich(self.G, self.jac.w, self.jac.category, self.ratios, self.prior, mask=mask,
                        constraints=self.constraints)[0]

    def level_sigma(self, L, levels):
        """sigma_stat of a functional [(state, v, J, sign)] of 127I2 levels, with factor L."""
        g = np.zeros(len(self.names))
        for s, v, J, sign in levels:
            g[:len(self.free)] += sign*self.D[self.jac.index[("127I2", s, v, J)]]*MHZ_PER_CM
        return float(np.linalg.norm(g@L))


# ---------------------------------------------------------------------------------------------------------------
# held-out residuals

def _row_levels(rows, i):
    return [(k[1], k[2], k[3], float(s)) for k, s in rows.levels[i]]


def heldout(cov, cv, benchmark_rows):
    """Every held-out center residual with its measurement sigma and its own fold's sigma_stat: the continuation's
    v folds (vmod:4:i), the final fit's benchmark rows (future), and the core fit's folds (component frequencies;
    the precision rows of that level or source removed from the covariance as from the fit)."""
    jac = cov.jac
    rows = jac.rows
    ids = {k: i for i, k in enumerate(rows.id)}
    out = []
    for fold in cv["continuation_folds"]:
        held = np.array([ids[k] for k in fold["ids"]])
        mask = np.ones(len(rows), bool)
        mask[held] = False
        L = cov.masked(mask)
        sig = np.sqrt(np.sum((cov.G[held]@L)**2, axis=1))
        for k, r, s in zip(held, np.array(fold["pred"])-np.array(fold["observed"]), sig):
            out.append(dict(fold=fold["fold"], id=str(rows.id[k]), source=str(rows.source[k]),
                            category=str(rows.category[k]), kind="center", residual=float(r),
                            sigma_meas=float(rows.sigma[k]), sigma_stat=float(s), levels=_row_levels(rows, k)))
    bench = cv["benchmark"]
    pos = {k: i for i, k in enumerate(benchmark_rows.id)}
    # the release took these rows' sigma from the table before its re-reduction (bench["sigma_meas"]); without that
    # record, the rows' own sigma
    stated = bench.get("sigma_meas") or [float(benchmark_rows.sigma[pos[k]]) for k in bench["ids"]]
    for k, r, sm in zip(bench["ids"], np.array(bench["pred"])-np.array(bench["observed"]), stated):
        j = pos[k]
        lv = _row_levels(benchmark_rows, j)
        out.append(dict(fold="future", id=k, source=str(benchmark_rows.source[j]),
                        category=str(benchmark_rows.category[j]), kind="center", residual=float(r),
                        sigma_meas=float(sm), sigma_stat=cov.level_sigma(cov.L, lv), levels=lv))
    sigma = _stated_sigmas()
    precise = np.isin(jac.category, PRECISE)
    touched = [{(k[1], k[2]) for k, _ in lv} for lv in rows.levels]
    for r in cv["core_folds"]:
        if r["kind"] != "absolute" or not r["fold"].startswith("physical"):
            continue
        line = Line.parse(r["line"])
        Ju = line.J_lower+(1 if line.branch == "R" else -1)
        lv = [("B", line.v_upper, Ju, 1.), ("X", line.v_lower, line.J_lower, -1.)]
        what = r["fold"].split(":", 1)[1]
        if r["fold"].startswith("physical_level"):
            st, v = what.split(":")
            drop = precise & np.array([(st, int(v)) in t for t in touched])
        else:
            drop = precise & (jac.rows.source == what)
        out.append(dict(fold=r["fold"], id=r["raw_ids"], source=r["sources"], category="comb", kind="component",
                        residual=float(r["error_MHz"]),
                        sigma_meas=float(np.median([sigma[i] for i in r["raw_ids"].split(";")])),
                        sigma_stat=cov.level_sigma(cov.masked(~drop), lv), levels=lv,
                        line=r["line"].replace("127I2 ", ""), component=json.loads(r["identity"])[1]))
    return out


def _stated_sigmas():
    sigma = {}
    for ds in load_all():
        scale = MHZ_PER_CM if ds.unit == "cm-1" else 1.
        for k, obs in enumerate(ds.observations):
            sigma[f"{ds.id}:{k}"] = obs.uncertainty*scale
    return sigma


def chosen_scale(cv, i):
    """The 1g(1Pig) C5 scale with the lowest training cost in hyperfine fold i (the nested choice)."""
    return min(HF_SCALES, key=lambda g: cv["hyperfine_folds"][str(g)]["folds"][i]["training"]["out_cost"])


def hf_heldout(cv, levels=None, workers=12, cache=None):
    """Per-interval residuals of the nested hyperfine cross-validation (B v' > 43): in each v' block the C5 scale
    with the lowest training cost, its fitted b and angles, every held-out interval recomputed exactly; and the core
    folds' intraline residuals, with each component's offset from the line center under the core fit."""
    from ..research import ResearchModel
    rows, params = hb.load_data()
    levels = levels or hb.Levels()
    out, models = [], {}
    for i in range(len(cv["hyperfine_folds"][str(HF_SCALES[1])]["folds"])):
        g = chosen_scale(cv, i)
        run = cv["hyperfine_folds"][str(g)]
        r = run["folds"][i]
        if g not in models:
            constants = hb.state_constants(False, hb.parse_overrides(run["override"]))
            terms = hb.fill_terms(hb.needed_levels(rows, params), constants, workers, cache=cache)
            models[g] = hb.Model(levels, terms, use_s2=True, asymptotes=hb.ASYMPTOTES, orders=run.get("orders"))
        model, b, theta = models[g], np.array(r["b"]), r["theta"]
        held = [(lo, hi) for lo, hi in r["holdout"]]
        offsets = {}
        for row in rows:
            line = row["line"]
            if row["inside"] or not any(lo <= line.v_upper <= hi for lo, hi in held):
                continue
            if line not in offsets:
                offsets[line] = hb.line_offsets(levels, line, {j: model.constants(b, theta, line.v_upper, j)
                                                               for j in hb.neighbours(hb.upper_J(line))})
            o = offsets[line]
            out.append(dict(region="B beyond", fold=f"hf:{held[0][0]}-{held[0][1]}", source=row["source"],
                            line=str(line), v=line.v_upper, J=hb.upper_J(line), observed=row["observed"],
                            sigma=row["sigma"], predicted=o[row["rank"]]-o[row["ref_rank"]],
                            offset=o[row["rank"]], ref_offset=o[row["ref_rank"]]))
    sigma = _stated_sigmas()
    core = ResearchModel("balanced")
    offsets = {}
    for r in cv["core_folds"]:
        if r["kind"] != "intraline":
            continue
        line = Line.parse(r["line"])
        ranks = json.loads(r["identity"])[1:]
        if line not in offsets:
            _, comps = core.hyperfine_components(line.v_upper, line.v_lower, line.J_lower, line.branch)
            offsets[line] = {int(c.label[1:]): c.offset for c in comps if c.label}
        o = offsets[line]
        lo, hi = ranks           # identity [line, lo, hi]: predicted and observed are a_hi - a_lo
        out.append(dict(region="core", fold=r["fold"], source=r["sources"], line=r["line"].replace("127I2 ", ""),
                        v=line.v_upper, J=line.J_lower+(1 if line.branch == "R" else -1),
                        observed=float(r["observed_median_MHz"]),
                        sigma=float(np.median([sigma[i] for i in r["raw_ids"].split(";")])),
                        predicted=float(r["predicted_MHz"]), offset=o.get(hi), ref_offset=o.get(lo),
                        n=int(r["n_measurements"])))
    return out


# ---------------------------------------------------------------------------------------------------------------
# sigma_model, sigma_HF and the data flags

def _minimize(fun, x0, bounds=None):
    if bounds:
        r = minimize(fun, x0, method="L-BFGS-B", bounds=bounds)
        r = minimize(fun, r.x, method="Powell", bounds=bounds, options=dict(xtol=1e-6, ftol=1e-10, maxfev=40000))
    else:
        r = minimize(fun, x0, method="Nelder-Mead", options=dict(xatol=1e-6, fatol=1e-8, maxiter=20000, maxfev=40000))
    return r.x, float(r.fun)


def _gauss_nll(r, var):
    return float(.5*np.sum(np.log(var)+r*r/var))


def fit_hyperfine(hf_rows, cv):
    """sigma of an interval Delta between two components: sqrt(sigma_meas^2 + 2 a^2 + (kappa Delta)^2); a component's
    offset o from the center carries sqrt(a^2 + (kappa o)^2). Core: one (a, kappa). B beyond: kappa(v') =
    k0 (1 + (v' - 42)/2)^p, fitted to the held-out intervals and to the held-out Chen 2004 / Pique 1986 eqQ (an eqQ'
    error dq scales the pattern by dq / |eqQ' - eqQ''|). X beyond: kappa_X(v'') log-linear from the core value at
    v'' = 17 to the within-line scatter of the emission components at v'' = 48-54."""
    from ..research import ResearchModel
    out = {}
    core = [o for o in hf_rows if o["region"] == "core"]
    r = np.array([o["predicted"]-o["observed"] for o in core])
    sm = np.array([o["sigma"] for o in core])
    D = np.array([abs(o["predicted"]) for o in core])
    (la, lk), nll = _minimize(lambda t: _gauss_nll(r, sm**2+2*np.exp(2*t[0])+(np.exp(t[1])*D)**2),
                              [np.log(.03), np.log(1e-4)])
    out["core"] = dict(a_MHz=float(np.exp(la)), kappa=float(np.exp(lk)), n=len(r), nll=nll)
    beyond = [o for o in hf_rows if o["region"] == "B beyond"]
    rb = np.array([o["predicted"]-o["observed"] for o in beyond])
    smb = np.array([o["sigma"] for o in beyond])
    Db = np.array([abs(o["predicted"]) for o in beyond])
    vb = np.array([o["v"] for o in beyond])
    eq = []
    for i in range(len(cv["hyperfine_folds"][str(HF_SCALES[1])]["folds"])):
        g = chosen_scale(cv, i)
        for score in cv["hyperfine_folds"][str(g)]["folds"][i]["eqQ"]:
            eq += [(x["v"], x["resid"], x["sigma"]) for x in score["rows"]]
    eq = np.array(eq)
    chen = {(int(row["v"]), int(row["J"])): row["observed"] for row in hyperfine_rows.release()[1]
            if row["parameter"] == "eqQ"}
    dq = np.array([abs(EQQ_X-np.mean([q for (v, _), q in chen.items() if v == int(x)])) for x in eq[:, 0]])

    def nll_b(t):
        a, k0, p = t
        kap = k0+p*np.log1p((vb-42)/2.)
        ke = k0+p*np.log1p((eq[:, 0]-42)/2.)
        return (_gauss_nll(rb, smb**2+2*np.exp(2*a)+(np.exp(kap)*Db)**2)
                + _gauss_nll(eq[:, 1], eq[:, 2]**2+(np.exp(ke)*dq)**2))
    (a, k0, pw), nll = _minimize(nll_b, [np.log(.05), np.log(2e-4), 1.])
    out["B beyond"] = dict(a_MHz=float(np.exp(a)), log_kappa0=float(k0), power=float(pw), v_ref=42, v_scale=2.,
                           n_intervals=len(rb), n_eqQ=len(eq), v_range=[int(min(vb.min(), eq[:, 0].min())),
                                                                       int(max(vb.max(), eq[:, 0].max()))], nll=nll)
    model = ResearchModel("balanced_expanded")
    dev, off, sig = [], [], []
    for ds in load_all():
        if ds.id not in EMISSION_SETS:
            continue
        unit = MHZ_PER_CM if ds.unit == "cm-1" else 1.
        lines = {}
        for o in ds.observations:
            if o.kind != "frequency" or not o.component:
                continue
            label = f"{o.line.branch}({o.line.J_lower}) {o.line.v_upper}-{o.line.v_lower}"
            center = model.frequency(label)
            p = model.frequency(label, o.component)
            lines.setdefault(label, []).append((o.value*unit-p, p-center, o.uncertainty*unit))
        for rows in lines.values():
            if len(rows) < 3:
                continue
            x = np.array(rows)
            dev += list(x[:, 0]-x[:, 0].mean())
            off += list(x[:, 1]-x[:, 1].mean())
            sig += list(x[:, 2])
    dev, off, sig = map(np.array, (dev, off, sig))
    (lk,), nll = _minimize(lambda t: _gauss_nll(dev, sig**2+(np.exp(t[0])*off)**2), [np.log(.01)])
    kx, kc = float(np.exp(lk)), out["core"]["kappa"]
    out["X beyond"] = dict(kappa_17=kc, kappa_51=kx, slope_per_v=float(np.log(kx/kc)/(51-17)), n_components=len(dev),
                           validated_v=[48, 54], nll=nll)
    return out


def _training_levels(rows):
    return [[(k[1], k[2], k[3]) for k, _ in lv] for lv in rows.levels]


def _coverage(levels, rows):
    out = {}
    for i in rows:
        for s, v, J in levels[i]:
            out.setdefault((s, v), set()).add(J)
    return {k: np.array(sorted(v)) for k, v in out.items()}


def _near(cov, s, v, J, margin=J_MARGIN):
    Js = cov.get((s, v))
    return Js is not None and bool(np.any(abs(Js-J) <= margin))


def fold_cover(rows):
    """fold -> {(state, v): J of the precision training rows' levels} (the fold's own rows removed)."""
    levels = _training_levels(rows)
    precise_all = [i for i in range(len(levels)) if rows.category[i] in PRECISE]
    cache = {}

    def cover(fold):
        if fold not in cache:
            keep = precise_all
            if fold.startswith("physical_level"):
                s, v = fold.split(":")[1:]
                keep = [i for i in keep if not any(a == s and b == int(v) for a, b, _ in levels[i])]
            elif fold.startswith("physical_source"):
                src = fold.split(":", 1)[1]
                keep = [i for i in keep if rows.source[i] != src]
            cache[fold] = _coverage(levels, keep)
        return cache[fold]
    return cover


def line_class(cover, levels):
    return "precise" if all(_near(cover, s, v, J) for s, v, J in levels) else "atlas"


def fit_model(held, hf, ratios, rows):
    """sigma_model: core, two classes (both levels within 15 in J of a precision training row at their v, or not),
    one value per line; beyond the core, a random effect per held-out level, sigma_m(state, v, J) =
    exp(a + b (v - v0) + c (J/100)^2), fitted by marginal likelihood (rows of one level and J bin share it)."""
    precise_cover = fold_cover(rows)

    def noise(o):
        sm = o["sigma_meas"]
        if o["category"] in NOISE_FROM_BIRGE:
            sm *= ratios[o["category"]]["birge"]
        var = sm**2+o["sigma_stat"]**2
        if o["kind"] == "component":
            var += hf["core"]["a_MHz"]**2
        return var

    core = {"precise": [], "atlas": []}
    beyond = {"B": [], "X": []}
    for o in held:
        lv = [(s, v, J) for s, v, J, _ in o["levels"]]
        out = [(s, v, J) for s, v, J in lv if v > CORE_V[s]]
        if not out:
            cls = line_class(precise_cover(o["fold"]), lv)
            core[cls].append((o["residual"], noise(o), o["fold"], tuple(lv)))
        else:
            s, v, J = max(out, key=lambda x: (x[0] == "X", x[1]))
            beyond[s].append((o["residual"], noise(o), (o["fold"], s, v, J//20), v, J, len(lv)-1))
    result = {}
    for cls, crows in core.items():
        groups = {}
        for r, var, fold, lv in crows:
            groups.setdefault((fold, lv), []).append((r, var))

        def nll(t, groups=groups):
            tau2 = np.exp(2*t[0])
            total = 0.
            for g in groups.values():
                x = np.array(g)
                r, s2 = x[:, 0], x[:, 1]
                a, b = np.sum(1/s2), np.sum(r/s2)
                total += .5*(np.sum(np.log(s2))+np.log(1+tau2*a)+np.sum(r*r/s2)-tau2*b*b/(1+tau2*a))
            return total
        (lt,), val = _minimize(nll, [np.log(.2)])
        result[f"core {cls}"] = dict(sigma_MHz=float(np.exp(lt)), n_rows=len(crows), n_groups=len(groups), nll=val)
    core_level_var = {cls: result[f"core {cls}"]["sigma_MHz"]**2/2 for cls in core}
    for s, brows in beyond.items():
        groups = {}
        for r, var, key, v, J, others in brows:
            groups.setdefault(key, []).append((r, var+others*core_level_var["atlas"], v, J))
        J_top = int(max(J for *_, J, _ in brows))

        def nll(t, groups=groups, s=s, J_top=J_top):
            a, b, c = t
            total = 0.
            for (_, _, v, _), g in groups.items():
                x = np.array(g)
                r, s2, J = x[:, 0], x[:, 1], x[:, 3]
                tau2 = np.exp(2*(a+b*(v-V0[s])+c*(min(np.mean(J), J_top)/100)**2))
                A, B = np.sum(1/s2), np.sum(r/s2)
                total += .5*(np.sum(np.log(s2))+np.log(1+tau2*A)+np.sum(r*r/s2)-tau2*B*B/(1+tau2*A))
            return total
        # model error may not fall toward dissociation or with J (X: Martin's noise hides it beyond v''25)
        bounds = [(-10., 10.), (0., 1.), (0., 5.)] if s == "X" else None
        (a, b, c), val = _minimize(nll, [np.log(5.), .05, .1], bounds)
        vs = [v for *_, v, J, _ in brows]
        result[f"{s} beyond"] = dict(log_sigma0=float(a), slope_per_v=float(b), J2_coefficient=float(c), v0=V0[s],
                                     J_clip=J_top, n_rows=len(brows), n_groups=len(groups),
                                     v_range=[int(min(vs)), int(max(vs))], nll=val,
                                     constraints="slope, J2_coefficient >= 0" if bounds else None)
    return result


def sigma_m(model, state, v, J, line_class="atlas"):
    """sigma_model share of one level (MHz)."""
    if v <= CORE_V[state]:
        return model[f"core {line_class}"]["sigma_MHz"]/np.sqrt(2)
    m = model[f"{state} beyond"]
    J = min(J, m["J_clip"])
    x = v-m["v0"]
    return float(np.exp(m["log_sigma0"]+m["slope_per_v"]*x+m.get("quadratic_per_v2", 0.)*x*x
                        + m["J2_coefficient"]*(J/100)**2))


def calibration(held, fits, ratios, rows):
    """In-sample check of the fitted functions on the held-out rows they were fitted to: rms z and the fraction
    within 1 and 2 sigma, per region, J bin and data class; per held-out level (group means) as well."""
    hf, model = fits["hyperfine"], fits["model"]
    cover = fold_cover(rows)
    out, groups = {}, {}

    def add(key, zval):
        out.setdefault(key, []).append(zval)
    for o in held:
        lv = [(s, v, J) for s, v, J, _ in o["levels"]]
        beyond = [(s, v, J) for s, v, J in lv if v > CORE_V[s]]
        sm = o["sigma_meas"]*(ratios[o["category"]]["birge"] if o["category"] in NOISE_FROM_BIRGE else 1.)
        cls = "atlas" if beyond else line_class(cover(o["fold"]), lv)
        var = sm**2+o["sigma_stat"]**2+sum(sigma_m(model, s, v, J, cls)**2 for s, v, J in lv)
        if o["kind"] == "component":
            var += hf["core"]["a_MHz"]**2
        zval = o["residual"]/np.sqrt(var)
        if beyond:
            s, v, J = max(beyond, key=lambda x: (x[0] == "X", x[1]))
            add((f"{s} {v//10*10}-{v//10*10+9}", "all J"), zval)
            add((f"{s} beyond", f"J {J//50*50}-{J//50*50+49}"), zval)
            add((f"{s} beyond", o["category"]), zval)
            groups.setdefault((o["fold"], s, v, J//20), []).append((o["residual"], var-sigma_m(model, s, v, J)**2,
                                                                     sigma_m(model, s, v, J)**2))
        else:
            add(("core", o["fold"].split(":")[0]), zval)
            add(("core", cls), zval)
    report = {f"{a} | {b}": dict(n=len(x), rms_z=float(np.sqrt(np.mean(np.square(x)))),
                                 within_1=float(np.mean(np.abs(x) <= 1)), within_2=float(np.mean(np.abs(x) <= 2)))
              for (a, b), x in sorted(out.items())}
    gz = {}
    for (_, s, v, _), g in groups.items():
        x = np.array(g)
        w = 1/x[:, 1]
        mean = np.sum(w*x[:, 0])/np.sum(w)
        zg = mean/np.sqrt(1/np.sum(w)+x[0, 2])
        for key in ("all", f"{s} {v//10*10}-{v//10*10+9}"):
            gz.setdefault(key, []).append(zg)
    for key, x in sorted(gz.items()):
        x = np.array(x)
        report[f"per held-out level | {key}"] = dict(n=len(x), rms_z=float(np.sqrt(np.mean(x**2))),
                                                     within_1=float(np.mean(abs(x) <= 1)),
                                                     within_2=float(np.mean(abs(x) <= 2)))
    return report


def _runs(J):
    """Sorted J values as [first, last] runs of consecutive integers."""
    out = []
    for j in J:
        if out and j == out[-1][1]+1:
            out[-1][1] = int(j)
        else:
            out.append([int(j), int(j)])
    return out


def coverage_tables(rows):
    """J of every training level per (state, v), all rows and precision rows (comb, Doppler-free, emission)."""
    levels = _training_levels(rows)
    every = _coverage(levels, range(len(levels)))
    precise = _coverage(levels, [i for i in range(len(levels)) if rows.category[i] in PRECISE])
    return ({f"{s}:{v}": _runs(J) for (s, v), J in sorted(every.items())},
            {f"{s}:{v}": _runs(J) for (s, v), J in sorted(precise.items())})


RECIPE = (
    "sigma_center = sqrt(sigma_stat^2 + sigma_model^2); a hyperfine component adds sigma_HF = sqrt(a^2 + (kappa o)^2) "
    "(o its offset from the center); an intraline interval Delta has sqrt(2 a^2 + (kappa Delta)^2) alone. sigma_stat: "
    "the sandwich covariance of the weighted Gauss-Newton estimator at the final fit over the {tails} free tail terms "
    "(B:dDe = X:dDe), the {core} frozen core potential coefficients, 7 inherited dispersion log factors (prior 0.1) and 5 "
    "calibrations (priors: their own fits), each data category scaled by its Birge ratio (at least 1), from the SVD "
    "of the weighted Jacobian (directions below 1e-8 of the largest singular value dropped); propagated through "
    "Hellmann-Feynman level derivatives. sigma_model: fitted by marginal likelihood to held-out center residuals (the "
    "core fit's folds; the continuation's four v folds), measurement sigmas removed in quadrature (atlases and "
    "constants scaled by their in-sample Birge ratios), each held-out level a random effect; core per line in two "
    "classes (both levels within 15 in J of a precision training row at their v, or not); beyond the core per level, "
    "exp(a + b (v - v0) + c (min(J, J_clip)/100)^2), X with b, c >= 0. sigma_HF: Gaussian likelihood on held-out "
    "intraline intervals (core folds; nested hyperfine cross-validation beyond, with its held-out Chen 2004 / Pique "
    "1986 eqQ converted through |eqQ' - eqQ''|) and, for X beyond, the emission components' within-line scatter; "
    "kappa held at its last validated v beyond. Status of a level: measured (training rows within 15 in J at its v), "
    "interpolated (inside its v's J range, or between neighbors within 2 in v), else extrapolated.")


def recipe(names):
    """RECIPE with the counts of the covariance's parameters: tail terms (dDe among them) and core coefficients."""
    tails = {f"{st}:{n}" for st, c in continuation.build_curves().items() for n in c.tail_names}-{"B:dDe"}
    rest = [n for n in names if not n.startswith("calibration") and ":logC" not in n]
    n_tail = sum(1 for n in rest if n in tails)
    return RECIPE.format(tails=n_tail, core=len(rest)-n_tail)


def build(cv=None, cache=None, workers=12):
    """Every stage: (uncertainty record without the spread check, intermediate products). Minutes: the joint
    Jacobian solves every training level with 93 derivative columns."""
    cv = cv or cross_validation()
    bench = training.release("continuation", "benchmark")
    extra = {k for lv in bench.levels for k, _ in lv}
    for r in cv["core_folds"]:
        if r["kind"] == "absolute":
            line = Line.parse(r["line"])
            Ju = line.J_lower+(1 if line.branch == "R" else -1)
            extra |= {("127I2", "B", line.v_upper, Ju), ("127I2", "X", line.v_lower, line.J_lower)}
    jac = Jacobian(extra_keys=extra)
    cov = Covariance(jac)
    held = heldout(cov, cv, bench)
    hf_rows = hf_heldout(cv, workers=workers, cache=cache)
    hf = fit_hyperfine(hf_rows, cv)
    fits = dict(hyperfine=hf, model=fit_model(held, hf, cov.ratios, jac.rows))
    every, precise = coverage_tables(jac.rows)
    record = dict(
        schema_version=1, model="research-20261020-balanced-expanded",
        created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
        recipe=recipe(cov.names), pending_outcomes_used=False,
        statistical=dict(parameters=cov.names, factor=cov.L.tolist(), fd_step=FD_STEP, rcond=RCOND,
                         dropped_directions=cov.info["dropped"], birge=cov.ratios,
                         priors={cov.names[j]: s for j, s in cov.prior.items()}),
        model_error=fits["model"], hyperfine=dict(fits["hyperfine"], status=HF_STATUS),
        coverage=dict(all=every, precise=precise, margin_J=J_MARGIN, neighbor_v=2,
                      note="J of the levels of the final fit's training rows, as [first, last] runs, per state:v"),
        checks=dict(heldout_calibration=calibration(held, fits, cov.ratios, jac.rows)))
    return record, dict(jacobian=jac, covariance=cov, heldout=held, hf_heldout=hf_rows, fits=fits)


def hyperfine_fold(scale, block, levels=None, workers=12, cache=None):
    """The nested hyperfine cross-validation's fit for one 1g(1Pig) C5 scale with the v' block (lo, hi) held out:
    (b, theta, training costs)."""
    rows, params = hb.load_data()
    saute = hb.so.perturber_constants()["1g(1Pig)"]["C5"]
    overrides = dict(hb.FIT_OVERRIDES, **{"1g(1Pig)": dict(C5=round(scale*saute, 6))})
    constants = hb.state_constants(False, overrides)
    model = hb.Model(levels or hb.Levels(), hb.fill_terms(hb.needed_levels(rows, params), constants, workers, cache=cache))
    b, theta = hb.fit(model, rows, params, dict(u=30., g=30.), False, holdout=(tuple(block),), keys=list(hb.FIT_KEYS))
    return b, theta, model.training
