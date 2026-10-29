"""The continuation: tail terms past the core's edge, fitted with every core coefficient frozen.

The core fit's levels sample its potentials only between their classical turning points (X 2.43-3.06 A, B 2.65-4.66
A). At knots on that edge each frozen series continues as its own Taylor polynomial (X exponent C^2, B exponent C^1,
rotational functions C^1) plus new terms t^(k+1) B_j(t), B_j the Bernstein polynomials of degree m-1 in t, the
normalized distance past the knot in the function's bounded coordinate. Outer exponent terms carry a further
(1-t)^2, so the asymptote stays the dispersion series and De. Also fitted: a change dC6 of the X C6 (switched on past
the knot as t^(k+1)), a shift dDe of both dissociation energies (shared, so the limits keep the atomic splitting;
units of 0.01 cm-1), cubic B-splines on V(X) over 4.2-6.6 A (units of 0.1 cm-1), and the offset of the Orsay IV
atlas. The settings are the release's `structure` record (research_expanded reads the same one).

Objective: the core fit's weights over every 127I2 training row up to X v''=108, B v'=80, the core rows included (they
keep the core where it was), plus the single-well penalty (-ln of the Morse factor increasing, to 40 A) and f > 0.2.
Solved by Levenberg-Marquardt, then undamped Gauss-Newton to convergence.
"""
from __future__ import annotations

import json
import os
from math import comb, factorial
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import scipy.sparse as sp
from numpy.polynomial.chebyshev import chebder, chebval
from scipy.interpolate import BSpline

from ..constants import MHZ_PER_CM
from . import rows as training
from .core import FLOORS, release_spec as core_spec, weights
from .curves import Config, Curve, long_range
from .levels import Levels
from .optimize import gauss_newton, marquardt

ROOT = Path(__file__).resolve().parents[3]
RELEASE = Path(os.environ.get("I2SPEC_FIT_CONTINUATION") or ROOT / "src" / "i2spec" / "data" / "research" / "20261020" / "balanced_expanded")
CORE_V = dict(X=17, B=43)
CORE_J = dict(X=200, B=201)
DATA_V = dict(X=108, B=80)
OFFSETS = ("orsay4",)


def release_spec():
    return json.loads((RELEASE / "model.json").read_text())


def core_level(key):
    _, state, v, J = key
    return v <= CORE_V[state] and J <= CORE_J[state]


class Tail:
    """A Chebyshev series in a mapped coordinate z, continued past z_k by Taylor(k) + sum_j a_j t^(k+1) B_j(t)."""

    def __init__(self, coefficients, shift, scale, zk, side, k, names, end=None, vanish=0):
        self.c, self.shift, self.scale = np.asarray(coefficients, float), shift, scale
        self.zk, self.side, self.k, self.names, self.end, self.vanish = zk, side, k, list(names), end, vanish
        x = (zk-shift)/scale
        self.taylor = [chebval(x, chebder(self.c, i) if i else self.c)/scale**i/factorial(i) for i in range(k+1)]
        self.span = abs((end if end is not None else side)-zk)       # t = 1 at `end` (default: the coordinate's limit)

    def region(self, z):
        return (z-self.zk)*self.side > 0

    def basis(self, z):
        t = np.where(self.region(z), (z-self.zk)*self.side/self.span, 0.)
        if self.end is not None:
            t = np.minimum(t, 1.)
        m = len(self.names)
        if not m:
            return np.zeros(np.shape(z)+(0,))
        return np.stack([t**(self.k+1)*(1-t)**self.vanish*comb(m-1, j)*t**j*(1-t)**(m-1-j) for j in range(m)], axis=-1)

    def continuation(self, z, a):
        dz = z-self.zk
        return sum(d*dz**i for i, d in enumerate(self.taylor))+(self.basis(z)@a if len(a) else 0.)


class TailCurve(Curve):
    """The frozen core curve of one state, with its tails (structure: the release's structure record)."""

    def __init__(self, state, config: Config, values, structure):
        super().__init__(state, config, values)
        st = structure
        knots, orders = st["knots"], st["orders"][state]
        self.frozen = dict(values)
        self.free_de = st["free_de"][state]
        self.de_unit, self.window_unit, self.new_cn_unit = st["de_unit"], st["window_unit"], st["new_cn_unit"]
        self.window = tuple(st["window_knots"][state]) if st["window_knots"][state] else None
        self.window_names = [f"w{j}" for j in range(1, len(self.window)-4+1)] if self.window else []
        self.dispersion = tuple(int(n) for n in st["dispersion"][state])
        re = np.exp(values["logRe"])
        q = config.powers(state)[1]
        self.ref = re*config.ref(state)
        yq = lambda R: (R**q-self.ref**q)/(R**q+self.ref**q)     # noqa: E731
        c = [values[f"c{i}"] for i in range(config.order(state))]
        shift, scale = (.2, .6) if state == "X" and config.mapped_x else (0., 1.)
        k = st["continuity"][state]
        self.exponent_tails = [Tail(c, shift, scale, yq(knots[f"{state}_out"]), 1, k,
                                    [f"to{j}" for j in range(1, orders["outer"]+1)], vanish=st["vanish"][state])]
        if f"{state}_in" in knots:
            self.exponent_tails.append(Tail(c, shift, scale, yq(knots[f"{state}_in"]), -1, k,
                                            [f"ti{j}" for j in range(1, orders["inner"]+1)], end=yq(st["inner_end"])))
        a = [values[f"alpha{i}"] for i in range(config.nalpha(state))]
        ashift, ascale = (.3, .4) if config.mapped_alpha else (0., 1.)
        kr = knots[f"{state}_rot"]
        self.rotation_tail = Tail(a, ashift, ascale, (kr-re)/(kr+re), 1, st["rotation_continuity"],
                                  [f"ro{j}" for j in range(1, orders["rotation"]+1)])
        self.inherited = long_range(state, values)
        self.dispersion_names = {n: f"dC{n}" if n in self.inherited else f"C{n}" for n in self.dispersion}
        self.tail_names = ([n for t in self.tails() for n in t.names]+list(self.dispersion_names.values())
                           + (["dDe"] if self.free_de else [])+self.window_names)
        for name in self.tail_names:
            self.values.setdefault(name, 0.)
        self.joins = tuple(sorted({knots[key] for key in knots if key.startswith(state)}
                                  | ({self.window[0], self.window[-1]} if self.window else set())))

    def tails(self):
        return self.exponent_tails+[self.rotation_tail]

    def window_basis(self, R):
        k = np.asarray(self.window, float)
        return np.nan_to_num(np.stack([BSpline.basis_element(k[i:i+5], extrapolate=False)(np.asarray(R, float))
                                       for i in range(len(k)-4)], axis=-1))

    def _exponent(self, values, yq):
        first = self.exponent_tails[0]
        out = chebval((yq-first.shift)/first.scale, first.c)
        for tail in self.exponent_tails:
            out = np.where(tail.region(yq), tail.continuation(yq, np.array([values[n] for n in tail.names])), out)
        return out

    def _unit(self, n):
        return self.inherited.get(n, self.new_cn_unit)

    def _pieces(self, values, R):
        cfg, v = self.config, self.frozen
        re, de = np.exp(v["logRe"]), np.exp(v["logDe"])
        p, q = cfg.powers(self.state)
        u = sum(c/R**n for n, c in self.inherited.items())
        ue = float(sum(c/re**n for n, c in self.inherited.items()))
        y = (R**p-self.ref**p)/(R**p+self.ref**p)
        yq = (R**q-self.ref**q)/(R**q+self.ref**q)
        yp = (R**p-re**p)/(R**p+re**p)
        outer = self.exponent_tails[0]
        S = None
        if self.dispersion:
            S = np.where(outer.region(yq), ((yq-outer.zk)/outer.span)**(outer.k+1), 0.)
            u = u+S*sum(values[self.dispersion_names[n]]*self._unit(n)/R**n for n in self.dispersion)
        beta = y*np.log(2*de/ue)+(1-y)*self._exponent(values, yq)
        argument = -beta*yp
        shape = u/ue*np.exp(np.clip(argument, -100, 100))
        self._last = (u, S)
        V = v.get("Te", 0.)+de*(1-shape)**2
        if self.free_de:
            self._de_switch = np.where(outer.region(yq), np.minimum((yq-outer.zk)/outer.span, 1.)**(outer.k+1), 0.)
            V = V+values["dDe"]*self.de_unit*self._de_switch
        if self.window:
            V = V+self.window_unit*self.window_basis(R)@np.array([values[n] for n in self.window_names])
        return V, de, shape, y, yq, yp, argument

    def potential_values(self, values, R):
        return self._pieces(values, np.asarray(R, float))[0]

    def _rotation(self, values, R):
        re = np.exp(self.frozen["logRe"])
        ya = (R-re)/(R+re)
        cfg = self.config
        ashift, ascale = (.3, .4) if cfg.mapped_alpha else (0., 1.)
        A = chebval((ya-ashift)/ascale, [self.frozen[f"alpha{i}"] for i in range(cfg.nalpha(self.state))])
        tail = self.rotation_tail
        A = np.where(tail.region(ya), tail.continuation(ya, np.array([values[n] for n in tail.names])), A)
        return 2*re/(R+re), ya, A

    def curves(self, values, R, ratio):
        R = np.asarray(R, float)
        envelope, _, A = self._rotation(values, R)
        return self.potential_values(values, R), 1+ratio*envelope*A

    def derivatives(self, values, R, ratio, names):
        R = np.asarray(R, float)
        V, de, shape, y, yq, yp, argument = self._pieces(values, R)
        factor = 2*de*(1-shape)*shape*yp*(abs(argument) < 100)*(1-y)
        envelope, ya, _ = self._rotation(values, R)
        dV, df = np.zeros((len(R), len(names))), np.zeros((len(R), len(names)))
        for tail in self.exponent_tails:
            basis = tail.basis(yq)
            for j, name in enumerate(tail.names):
                if name in names:
                    dV[:, names.index(name)] = factor*basis[:, j]
        basis = self.rotation_tail.basis(ya)
        for j, name in enumerate(self.rotation_tail.names):
            if name in names:
                df[:, names.index(name)] = ratio*envelope*basis[:, j]
        if self.dispersion:
            u, S = self._last
            for n in self.dispersion:
                name = self.dispersion_names[n]
                if name in names:
                    dV[:, names.index(name)] = -2*de*(1-shape)*shape*S*self._unit(n)/R**n/u*(abs(argument) < 100)
        if self.free_de and "dDe" in names:
            dV[:, names.index("dDe")] = self.de_unit*self._de_switch
        if self.window:
            basis = self.window_basis(R)
            for j, name in enumerate(self.window_names):
                if name in names:
                    dV[:, names.index(name)] = self.window_unit*basis[:, j]
        unknown = set(names)-set(self.tail_names)
        if unknown:
            raise ValueError(f"frozen core coefficients cannot be fitted here: {sorted(unknown)}")
        return dV, df

    def log_shape_values(self, values, R):
        """-ln(Morse factor), increasing in R for a single well. With the spline window the factor is read back
        from V on the outer branch (less any dDe), so the check covers it."""
        V, de, shape, *_ = self._pieces(values, R)
        if self.window:
            outer = R > np.exp(self.frozen["logRe"])
            if self.free_de:
                V = V-values["dDe"]*self.de_unit*self._de_switch
            shape = np.where(outer, 1-np.sqrt(np.maximum(V-self.frozen.get("Te", 0.), 0.)/de), shape)
        return -np.log(np.maximum(shape, 1e-300))


class SolverPotential:
    """BSplineSolver adapter for a TailCurve at fixed parameter values (127I2 only)."""

    def __init__(self, curve, values):
        self.curve, self.values = curve, values

    def __call__(self, R):
        return self.curve.potential_values(self.values, R)

    def adiabatic(self, R, mass_ratio=1.):
        return np.zeros_like(R)

    def nonadiabatic(self, R, mass_ratio=1.):
        return self.curve.curves(self.values, R, 1.)[1]-1


def build_curves(structure=None, core=None):
    """{state: TailCurve} on the frozen core, tail values from the structure's release (zeros if absent)."""
    core = core or core_spec()
    structure = structure or release_spec()["structure"]
    config = Config.from_spec(core["config"])
    return {s: TailCurve(s, config, core["values"][s], structure) for s in ("X", "B")}


class Problem:
    """Weighted residuals of rows against the tailed curves, both states solved."""

    def __init__(self, rows, curves, names, offsets=OFFSETS, tie=None, physical=True, box_padding=.6, refine_x=True):
        self.rows, self.curves, self.names, self.offsets = rows, curves, list(names), list(offsets)
        # tie: (state, name) -> the (state, name) whose value it takes; T maps the free terms (and the offsets)
        # onto every named term, so Jacobian columns of tied terms add up
        self.tie = dict(tie or {})
        self.free = [k for k in self.names if k not in self.tie]
        cols = [self.free.index(self.tie.get(k, k)) for k in self.names]
        cols += list(range(len(self.free), len(self.free)+len(self.offsets)))
        self.T = np.zeros((len(cols), len(self.free)+len(self.offsets)))
        self.T[np.arange(len(cols)), cols] = 1.
        self.levels = Levels(rows.level_set(), joins={s: c.joins for s, c in curves.items()}, box_padding=box_padding,
                             refine_x=refine_x)
        rr, cc, vv = [], [], []
        for i, lv in enumerate(rows.levels):
            for k, sign in lv:
                rr.append(i)
                cc.append(self.levels.index[k])
                vv.append(sign)
        self.M = sp.csr_matrix((vv, (rr, cc)), shape=(len(rows), len(self.levels.keys)))
        self.N = np.array([(rows.category == c).astype(float) for c in self.offsets]).T.reshape(len(rows), -1)
        self.w = weights(rows, FLOORS)
        self.physical_on = physical
        self.grid = {s: np.unique(np.r_[np.linspace(2.1 if s == "X" else 2.3, 6., 1500), np.geomspace(6.01, 40., 400)])
                     for s in curves}
        self.calls, self.cache = 0, None

    @property
    def parameter_names(self):
        return [f"{s}:{n}" for s, n in self.free]+[f"calibration:{c}" for c in self.offsets]

    def compress(self, x_full):
        """The free-parameter vector of a vector over every named term (and the offsets)."""
        return np.array([x_full[self.names.index(k)] for k in self.free]+list(x_full[len(self.names):]))

    def values(self, x):
        values = {s: dict(c.values) for s, c in self.curves.items()}
        for (s, n), z in zip(self.names, self.T@x):
            values[s][n] = float(z)
        return values

    def physical(self, x):
        if not self.physical_on:
            return np.zeros(0)
        values, out = self.values(x), []
        for s, curve in self.curves.items():
            R = self.grid[s]
            out.append(np.minimum(np.diff(curve.log_shape_values(values[s], R))/np.diff(R), 0.)/1e-5)
            out.append(np.minimum(curve.curves(values[s], R, 1.)[1]-.2, 0.)*1e4)
        return np.concatenate(out)

    def predict(self, x, jac=False):
        n = len(self.names)
        z = self.T@x
        E, D = self.levels.evaluate(self.curves, z[:n], self.names, jac=jac)
        pred = self.M@(E*MHZ_PER_CM)+self.N@z[n:]
        return pred, (np.column_stack((self.M@(D*MHZ_PER_CM), self.N))@self.T if jac else None)

    def evaluate(self, x):
        if self.cache is not None and np.array_equal(self.cache[0], x):
            return self.cache[1:]
        pred, J = self.predict(x, jac=True)
        residual = (pred-self.rows.value)*self.w
        missing = ~np.isfinite(residual)
        residual[missing], J[missing] = 1e6, 0.
        J = J*self.w[:, None]
        phys = self.physical(x)
        if len(phys):
            h = 1e-6
            P = np.column_stack([(self.physical(x+h*e)-self.physical(x-h*e))/(2*h) for e in np.eye(len(x))])
            residual, J = np.r_[residual, phys], np.vstack((J, P))
        self.calls += 1
        self.cache = (x.copy(), residual, J)
        return residual, J


def release_start(problem, spec=None):
    """The released tail values and offsets as a free-parameter vector of `problem`."""
    spec = spec or release_spec()
    full = ([spec["values"].get(s, {}).get(n, 0.) for s, n in problem.names]      # a term the release lacks: 0
            + [spec["calibration"][f"{c}_offset_MHz"] for c in problem.offsets])
    return problem.compress(np.array(full))


def fit(rows=None, structure=None, start=None, max_nfev=1, finish=40, fixed=None, initial=None, core=None):
    """Fit the tail terms. Defaults: the release's continuation rows and structure, started from (and, with
    max_nfev=1 and finish=40, the settings of) the release's own final run. `fixed`: {(state, name): value} held;
    `initial`: {state: {name: value}} starting values in place of the release's (e.g. `project`'s); `core`: a core
    model record ({"config", "values"}) to build on in place of the released core."""
    spec = release_spec()
    structure = structure or spec["structure"]
    curves = build_curves(structure, core)
    rows = rows if rows is not None else training.release("continuation")
    every = [(s, n) for s, c in curves.items() for n in c.tail_names]
    values = initial or spec["values"]
    for (s, n) in every:
        curves[s].values[n] = values.get(s, {}).get(n, 0.)
    fixed = dict(fixed or {})
    if structure.get("shared_de"):           # one shift for both states: fixing either fixes both
        for a, b in ((("X", "dDe"), ("B", "dDe")), (("B", "dDe"), ("X", "dDe"))):
            if a in fixed:
                fixed.setdefault(b, fixed[a])
    for (s, n), z in fixed.items():
        curves[s].values[n] = z
    names = [k for k in every if k not in fixed]
    tie = {("B", "dDe"): ("X", "dDe")} if structure.get("shared_de") and ("B", "dDe") in names else {}
    prob = Problem(rows, curves, names, tie=tie)
    if start is not None:
        x0 = np.asarray(start, float)
    elif initial is not None:
        x0 = prob.compress(np.array([curves[s].values[n] for s, n in prob.names]
                                    + [spec["calibration"][f"{c}_offset_MHz"] for c in prob.offsets]))
    else:
        x0 = release_start(prob, spec)
    scale = np.array([.01 if n.startswith("ro") else .1 if n.startswith("dC") else 1. for _, n in prob.free]+[50.]*len(prob.offsets))
    result = marquardt(prob, x0, max_nfev, scale)
    if finish:
        done = gauss_newton(prob, result["x"], finish)
        result = dict(done, nfev=result["nfev"]+done["nfev"])
    x = result["x"]
    pred, _ = prob.predict(x)
    return prob, SimpleNamespace(x=x, cost=result["cost"], nfev=result["nfev"], converged=result["converged"],
                                 parameters=dict(zip(prob.parameter_names, x.tolist())), pred=pred,
                                 physical_cost=float(np.sum(prob.physical(x)**2)))


def groups(rows):
    """Masks by where a row's levels lie: core, or the band of its highest out-of-core level."""
    def top(lv, state):
        return max((k[2] for k, _ in lv if k[1] == state and not core_level(k)), default=-1)
    vx = np.array([top(lv, "X") for lv in rows.levels])
    vb = np.array([top(lv, "B") for lv in rows.levels])
    core = np.array([all(core_level(k) for k, _ in lv) for lv in rows.levels])
    return [("core", core), ("X v18-25", (vx >= 18) & (vx <= 25)), ("X v26-47", (vx >= 26) & (vx <= 47)),
            ("X v48-60", (vx >= 48) & (vx <= 60)), ("X v61-89", (vx >= 61) & (vx <= 89)), ("X v90+", vx >= 90),
            ("B v44-50", (vb >= 44) & (vb <= 50)), ("B v51-62", (vb >= 51) & (vb <= 62)),
            ("B v63-72", (vb >= 63) & (vb <= 72)), ("B v73-77", (vb >= 73) & (vb <= 77)), ("B v78+", vb >= 78),
            ("core v, J>200", ~core & (vx < 0) & (vb < 0))]


def statistics(rows, residual):
    out = {}
    for name, g in groups(rows):
        for category in sorted(set(rows.category[g])):
            r = residual[g & (rows.category == category)]
            r = r[np.isfinite(r)]
            if len(r):
                out[f"{name} | {category}"] = dict(n=len(r), median_abs_MHz=float(np.median(abs(r))),
                                                   mean_MHz=float(r.mean()), rms_MHz=float(np.sqrt(np.mean(r*r))))
    return out


def holdout_mask(rows, K, i):
    """The rows of fold i of K: every row that touches a level beyond the core whose v % K == i."""
    return np.array([any(k[2] > CORE_V[k[1]] and k[2] % K == i for k, _ in lv) for lv in rows.levels])


def fold(i, K=4, rows=None, finish=40, **settings):
    """Refit without fold i of K (from the released values, or `settings` as for `fit`) and predict the held-out
    rows: (ids, pred, observed, result)."""
    rows = rows if rows is not None else training.release("continuation")
    held = holdout_mask(rows, K, i)
    prob, result = fit(rows=rows[~held], finish=finish, **settings)
    test = rows[held]
    curves = prob.curves
    values = prob.values(result.x)
    for s in curves:
        curves[s].values.update(values[s])
    pred, _ = Problem(test, curves, prob.names, prob.offsets, tie=prob.tie, physical=False).predict(result.x)
    return test.id, pred, test.value, result


def project(structure, values=None, old=None, reach=dict(X=10., B=20.)):
    """Starting values on another structure: each exponent tail of the old fit (default: the release) projected by
    weighted least squares onto the new structure's terms, over the radii the levels reach (to `reach` A outside,
    the inner end inside), each point weighted by how much the exponent moves V there; rotational tails and the spline
    correction projected over their span; other terms carried over, new ones 0. Returns {state: {name: value}}."""
    spec = release_spec()
    old = old or spec["structure"]
    values = values or spec["values"]
    a, b = build_curves(old), build_curves(structure)
    out = {}
    for s in ("X", "B"):
        va = dict(a[s].values, **values[s])
        out[s] = {n: va.get(n, 0.) for n in b[s].tail_names if n in va}
        re = float(np.exp(a[s].frozen["logRe"]))
        q = a[s].config.powers(s)[1]
        ref = a[s].ref
        for ta, tb in zip(a[s].exponent_tails, b[s].exponent_tails):
            if not tb.names:
                continue
            knot = (ref**q*(1+tb.zk)/(1-tb.zk))**(1/q)
            R = np.linspace(knot, reach[s], 800) if tb.side > 0 else np.linspace(old["inner_end"], knot, 200)
            z = (R**q-ref**q)/(R**q+ref**q)
            *_, factor = _exponent_weight(a[s], va, R)
            target = ta.basis(z)@np.array([va.get(n, 0.) for n in ta.names]) if ta.names else np.zeros(len(z))
            coef = np.linalg.lstsq(tb.basis(z)*factor[:, None], target*factor, rcond=None)[0]
            out[s].update(zip(tb.names, coef))
        ra, rb = a[s].rotation_tail, b[s].rotation_tail
        if rb.names:
            R = np.linspace(re*(1+rb.zk)/(1-rb.zk), reach[s], 400)
            ya = (R-re)/(R+re)
            target = ra.basis(ya)@np.array([va.get(n, 0.) for n in ra.names]) if ra.names else np.zeros(len(R))
            out[s].update(zip(rb.names, np.linalg.lstsq(rb.basis(ya), target, rcond=None)[0]))
        if b[s].window:
            R = np.linspace(b[s].window[0], b[s].window[-1], 600)
            target = (a[s].window_unit*a[s].window_basis(R)@np.array([va[n] for n in a[s].window_names])
                      if a[s].window else np.zeros(len(R)))
            coef = np.linalg.lstsq(b[s].window_unit*b[s].window_basis(R), target, rcond=None)[0]
            out[s].update(zip(b[s].window_names, coef))
    return out


def _exponent_weight(curve, values, R):
    """|dV/d(exponent)| on R: how much a change of the MLR exponent moves the potential there."""
    V, de, shape, y, yq, yp, argument = curve._pieces(values, np.asarray(R, float))
    return V, np.abs(2*de*(1-shape)*shape*yp*(1-y))
