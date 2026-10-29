"""The B hyperfine constants beyond the core: a smooth background plus explicit second-order sums.

B constants p (eqQ in MHz; C, d, delta in kHz) at (v, J):
    p(v, J) = <b_p(R)>_vJ + S2_p(v, J)
b_p, the background (first-order couplings and the far states of the other dissociation limits): the core fit's
radial form, a polynomial in z = (R^2 - re^2)/(R^2 + re^2)/scale (orders eqQ 5, C 5, d 3, delta 1), inside the
knot and its first-order Taylor continuation in z outside it; eqQ continues to its separated-atom limit with a
fitted curvature gamma. S2: second_order's constants of the nine states of 2P3/2 + 2P1/2 (couplings of Chen, de
Jong & Ye 2005, curves of de Jong et al. 1997 and Saute & Aubert-Frecon 1982), the two Omega = 1 pairs mixed by
angles theta_u, theta_g. A mixed state cos(t)|1'> + sin(t)|1''> on a fixed curve has S2 = c^2 A + s^2 B + 2cs X
exactly, so the angles cost nothing once A, B, X are cached per curve.

Data (hyperfine_rows): 127I2 intra-line intervals (core training; core benchmark, scored only; out of the core,
B v' 44-62), Chen et al. 2004's effective constants (v' 42-70) and Pique et al. 1986's (v' 71, 75). Each line's
intervals are linearized in its constants (four at each of J', J' +- 2) about the current model, which makes
every residual exactly linear in the background coefficients and in the trigonometric weights of the angles; the
linearization is renewed until the offsets stop moving.

Roles of the data. The background is anchored to the core fit's own predictions of the core intervals (so the core
keeps its hyperfine, with S2 inside). Out-of-core rows set theta_u and gamma, with model-error floors (OUT_FLOOR,
RELATIVE_FLOOR, EQQ_FLOOR). The curves are inputs; the release scales the 1g(1Pig) C5 (by 0.65, chosen by the
training cost) and sets the (3)0u- C5 to zero (FIT_OVERRIDES).

S2 terms are computed in parallel processes and cached by state and curve constants (CACHE).
"""
from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor
import hashlib
import json
import os
from multiprocessing import get_context
from pathlib import Path
import pickle
import time

import numpy as np
from scipy.optimize import least_squares

from ..bspline import BSplineSolver
from ..constants import MHZ_PER_CM, reduced_mass
from ..hyperfine import HyperfineParameters, level_structure, line_components
from ..observations import component_rank
from ..research import ResearchModel
from . import hyperfine_rows, rows as training
from . import second_order as so
from .continuation import SolverPotential, build_curves
from .core import release_spec as core_spec

CACHE = Path.home()/'.cache'/'i2spec'/'second_order'
POTENTIAL = 'fjc_dec6'
PARAMS = ('eqQ', 'C', 'd', 'delta')
ORDERS = {p: len(c) for p, c in core_spec()['hyperfine']['coefficients']['B'].items()}   # the core fit's B radial functions
KNOT = 4.45                                        # A, the B potential's outer knot (core edge)
# Separated-atom limits of the first-order couplings (R -> infinity): eqQ_B = B(2P3/2)/2 = -573.2 MHz (Bacis et
# al. 1980 Table I; their extrapolation -574 +- 3). Outside the knot a coupling with a limit continues as
#   b_k + b'_k D + gamma (D^2 - D^3/S) + (L - b_k - b'_k S) D^3/S^3,   D = z - z_k, S = 1/scale - z_k:
# the polynomial's value and slope at the knot, the limit L at R -> infinity, and a curvature gamma at the knot
# (gamma = (L - b_k - b'_k S)/S^2 is the plain quadratic). The core hardly sees gamma: it is the shape of the
# LCAO -> separated-atom transition (Bacis et al.: 3.5-4.5 A), fitted with the angles (--fit-shape).
ASYMPTOTES = dict(eqQ=4*so.COUPLINGS['B'][4])
STEPS = np.array([.5, .5, .5, .5])                 # finite-difference steps of the constants (MHz, kHz, kHz, kHz)
MIXED = {n: s['basis'][0] for n, s in so.PERTURBERS.items() if s['basis'][1] is not None}
PAIRS = dict(u=("1u(3Sigmau+)", "C1u"), g=("1g(1Pig)", "1g(3Sigmag-)"))


# ---------------------------------------------------------------------------------------------------------------
# S2 cache: per state, per curve (its constants), per (v, J): the projected constants of each coupling variant

def state_constants(experimental=False, overrides=None):
    """{state: dict of the curve constants actually used} (perturber_curves' rules)."""
    out = {}
    for name, c in so.perturber_constants().items():
        c = {k: v for k, v in c.items() if k != 'exp'} | ({k: v for k, v in c['exp'].items() if v is not None}
                                                           if experimental else {})
        c.update((overrides or {}).get(name, {}))
        out[name] = c
    return out


def _key(state, constants):
    text = json.dumps(dict(potential=POTENTIAL, state=state, constants=constants, grid=so.GRID, join=7.,
                           width=.5), sort_keys=True)
    return hashlib.sha1(text.encode()).hexdigest()[:16]


def _variants(state):
    """Coupling tuples whose projected constants give S2 of any mixture: a = |1'>, b = |1''>, m = (a + b)/sqrt2."""
    if state not in MIXED:
        return {'s': so.molecular_couplings()[state]}
    a, b = so.COUPLINGS["1'"+MIXED[state]], so.COUPLINGS["1''"+MIXED[state]]
    return dict(a=a, b=b, m=(a[0], a[1], *((x+y)/so.SQRT2 for x, y in zip(a[2:], b[2:]))))


def level_terms(radial, state, v, J):
    """{variant: constants (eqQ MHz, C, d, delta kHz)} of one state at B (v, J); the mixed variant m is turned
    into the cross term x = m - (a + b)/2."""
    G = {Jp: radial.green(state, v, J, Jp) for Jp in range(max(J-2, so.OMEGA[state]), J+3)}
    out = {}
    for key, c in _variants(state).items():
        blocks = {}
        for F in range(max(0, J-int(2*so.IA)), J+int(2*so.IA)+1):
            b = so.second_order_block(J, F, state, lambda k, kk, Jp: G[Jp][k, kk], couplings={state: c})
            if len(b[0]):
                blocks[F] = b
        fit, _, _ = so.project(J, blocks)
        out[key] = np.array([fit[p] for p in PARAMS])
    if 'm' in out:
        out['x'] = out.pop('m')-(out['a']+out['b'])/2
    return out


_WORKER = {}


def _init(states_constants):
    curves, values = potential()
    b = SolverPotential(curves['B'], values['B'])
    limit = float(b(np.array([300.]))[0])
    overrides = {s: c for s, c in states_constants.items()}
    _WORKER['radial'] = so.Radial(b, curves['B'].joins, perturber_curves(limit, overrides))


def perturber_curves(limit, constants):
    """so.perturber_curves with every state's constants given explicitly."""
    return so.perturber_curves(limit, overrides={s: {k: c[k] for k in ('C5', 'C6', 'De', 're', 'we')}
                                                for s, c in constants.items()})


def _task(args):
    state, J, vs = args
    radial = _WORKER['radial']
    return {(v, J): level_terms(radial, state, v, J) for v in vs}


def fill_terms(levels, constants, workers=16, cache=None):
    """Cache S2 terms for every (v, J) in `levels` and every state with the given curve constants."""
    cache = Path(cache or CACHE)
    cache.mkdir(parents=True, exist_ok=True)
    todo = []
    for state, c in constants.items():
        path = cache/f'{_key(state, c)}.pkl'
        have = pickle.loads(path.read_bytes()) if path.exists() else {}
        missing = sorted({lv for lv in levels if lv not in have})
        byJ = {}
        for v, J in missing:
            byJ.setdefault(J, []).append(v)
        todo += [(state, J, vs) for J, vs in byJ.items()]
    if todo:
        t = time.time()
        with ProcessPoolExecutor(workers, mp_context=get_context('spawn'), initializer=_init,
                                 initargs=(constants,)) as pool:
            results = list(pool.map(_task, todo, chunksize=4))
        new = {}
        for (state, _, _), r in zip(todo, results):
            new.setdefault(state, {}).update(r)
        for state, r in new.items():
            path = cache/f'{_key(state, constants[state])}.pkl'
            have = pickle.loads(path.read_bytes()) if path.exists() else {}
            have.update(r)
            path.write_bytes(pickle.dumps(have))
        print(f'S2 terms: {len(todo)} tasks in {time.time()-t:.0f} s', flush=True)
    out = {}
    for state, c in constants.items():
        have = pickle.loads((cache/f'{_key(state, c)}.pkl').read_bytes())
        out[state] = {lv: have[lv] for lv in levels}
    return out


# ---------------------------------------------------------------------------------------------------------------
# B levels and background moments

class Levels:
    """B energies (cm-1) and background moments <g_n(R)> per (v, J) on second_order_hyperfine's grid."""

    def __init__(self):
        curves, values = potential()
        self.solver = BSplineSolver(SolverPotential(curves['B'], values['B']), reduced_mass('127I2'),
                                    joins=curves['B'].joins, nlev=90, **so.GRID)
        spec = core_spec()
        re = np.exp(spec['values']['B']['logRe'])
        scale = spec['hyperfine']['radial_scale']['B']
        self.core = {p: np.array(c) for p, c in spec['hyperfine']['coefficients']['B'].items()}
        R = self.solver.R
        z = (R**2-re**2)/(R**2+re**2)/scale
        zk = (KNOT**2-re**2)/(KNOT**2+re**2)/scale
        n = max(ORDERS.values())
        inside = R <= KNOT
        dz, zinf = z-zk, 1/scale
        span = zinf-zk
        linear = [zk**k+k*zk**max(k-1, 0)*dz for k in range(n)]
        to_limit = [linear[k]-(zk**k+k*zk**max(k-1, 0)*span)*(dz/span)**3 for k in range(n)]
        self.g = np.array([np.where(inside, z**k, linear[k]) for k in range(n)]
                          + [z**k for k in range(n)]                    # the core fit's own form
                          + [np.where(inside, z**k, to_limit[k]) for k in range(n)]
                          + [np.where(inside, 0., (dz/span)**3)]       # times the limit
                          + [np.where(inside, 0., dz**2-dz**3/span)])  # times gamma
        self.n = n
        self.cache = {}

    def level(self, J):
        if J not in self.cache:
            E, psi = self.solver.wavefunctions(J)
            density = self.solver.W[:, None]*psi**2
            self.cache[J] = (E, self.g@density)            # (nlev,), (n, nlev)
        return self.cache[J]

    def energy(self, v, J):
        return float(self.level(J)[0][v])

    def moments(self, v, J, limit=False):
        """<g_k>, or with `limit` (moments of the continuation to an asymptote, the limit's weight, gamma's)."""
        m = self.level(J)[1]
        if limit:
            return m[2*self.n:3*self.n, v], float(m[3*self.n, v]), float(m[3*self.n+1, v])
        return m[:self.n, v]

    def plain_moments(self, v, J):
        return self.level(J)[1][self.n:2*self.n, v]


# ---------------------------------------------------------------------------------------------------------------
# the model: constants = M b + S0 + T w(theta) at each (v, J)

TRIG = ('u cc', 'u ss', 'u 2cs', 'g cc', 'g ss', 'g 2cs', 'eqQ gamma')


def trig(theta):
    """w: (cos^2, sin^2, 2 sin cos) of theta['u'], then of theta['g'] (degrees), then theta.get('gamma', 0)."""
    out = []
    for p in 'ug':
        c, s = np.cos(np.radians(theta[p])), np.sin(np.radians(theta[p]))
        out += [c*c, s*s, 2*c*s]
    return np.array(out+[theta.get('gamma', 0.)])


class Model:
    """B constants (eqQ MHz; C, d, delta kHz) = M b + S0 + T w(theta): M the background moments, S0 the unmixed
    states' S2, T the mixed pairs' S2 per trigonometric weight (a state cos|1'> + sin|1''> and its partner
    -sin|1'> + cos|1''> give c^2 (A1 + B2) + s^2 (B1 + A2) + 2cs (X1 - X2))."""

    def __init__(self, levels, terms, use_s2=True, asymptotes=ASYMPTOTES):
        self.levels, self.terms, self.use_s2, self.asymptotes = levels, terms, use_s2, dict(asymptotes)
        self.names = [f'{p}:{k}' for p in PARAMS for k in range(ORDERS[p])]
        self.nb = len(self.names)
        self._parts = {}

    def start(self):
        return np.concatenate([self.levels.core[p][:ORDERS[p]] for p in PARAMS])

    def parts(self, v, J):
        if (v, J) not in self._parts:
            m = self.levels.moments(v, J)
            ml, wl, wg = self.levels.moments(v, J, limit=True)
            M = np.zeros((4, self.nb))
            S0, T = np.zeros(4), np.zeros((4, len(TRIG)))
            i = 0
            for k, p in enumerate(PARAMS):
                if p in self.asymptotes:
                    M[k, i:i+ORDERS[p]] = ml[:ORDERS[p]]
                    S0[k] += self.asymptotes[p]*wl
                    T[k, 6] = wg
                else:
                    M[k, i:i+ORDERS[p]] = m[:ORDERS[p]]
                i += ORDERS[p]
            if self.use_s2:
                for state, t in self.terms.items():
                    if 's' in t[(v, J)]:
                        S0 = S0+t[(v, J)]['s']
                for i, (s1, s2) in enumerate(PAIRS.values()):
                    a, b = self.terms[s1][(v, J)], self.terms[s2][(v, J)]
                    T[:, 3*i] = a['a']+b['b']
                    T[:, 3*i+1] = a['b']+b['a']
                    T[:, 3*i+2] = a['x']-b['x']
            self._parts[(v, J)] = (M, S0, T)
        return self._parts[(v, J)]

    def constants(self, b, theta, v, J):
        M, S0, T = self.parts(v, J)
        return M@b+S0+T@trig(theta)


# ---------------------------------------------------------------------------------------------------------------
# lines

_X_CACHE = {}


def x_structure(v, J):
    if (v, J) not in _X_CACHE:
        if 'core' not in _X_CACHE:
            _X_CACHE['core'] = ResearchModel('balanced')
        _X_CACHE[(v, J)] = _X_CACHE['core']._structure('X', v, J)
    return _X_CACHE[(v, J)]


def upper_J(line):
    return line.J_lower+(1 if line.branch == 'R' else -1)


def neighbours(Jup):
    return [j for j in (Jup-2, Jup, Jup+2) if j >= 0]


def line_offsets(levels, line, constants_by_J):
    """{rank: offset MHz} of a line with B constants {J: array(4)}."""
    Jup, v = upper_J(line), line.v_upper
    E0 = levels.energy(v, Jup)

    def params(j):
        c = constants_by_J[j]
        return HyperfineParameters(eqQ=c[0], C=c[1], d=c[2], delta=c[3])

    s = level_structure(Jup, params, lambda j: (levels.energy(v, j)-E0)*MHZ_PER_CM, symmetry='u', dJ=2)
    comps = line_components(s, x_structure(line.v_lower, line.J_lower), Jup, line.J_lower)
    return {component_rank(c.label): c.offset for c in comps if c.label}


def linearize(levels, line, constants_fn):
    """(Js, c0, ranks, base, D): the line's offsets base + D (c - c0) in its constants c (4 per J of Js)."""
    Js = neighbours(upper_J(line))
    c0 = {j: constants_fn(line.v_upper, j) for j in Js}
    base = line_offsets(levels, line, c0)
    ranks = sorted(base)
    D = np.zeros((len(ranks), 4*len(Js)))
    for a, j in enumerate(Js):
        for k in range(4):
            c = {jj: cc.copy() for jj, cc in c0.items()}
            c[j][k] += STEPS[k]
            o = line_offsets(levels, line, c)
            D[:, 4*a+k] = [(o[r]-base[r])/STEPS[k] for r in ranks]
    return Js, np.concatenate([c0[j] for j in Js]), ranks, np.array([base[r] for r in ranks]), D


class System:
    """Weighted residuals of a set of rows, exactly linear in the background b and the weights w(theta):
    r = (y0 + K b + Q w) * scale (intervals through their lines' linearizations)."""

    def __init__(self, model, linear, rows, params, interval_floor, relative_floor=0., eqq_floor=0.):
        y0, K, Q, scale = [], [], [], []
        lines = {}
        for r in rows:
            line = r['line']
            if line not in lines:
                Js, c0, ranks, base, D = linear[line]
                parts = [model.parts(line.v_upper, j) for j in Js]
                Ms = np.vstack([p[0] for p in parts])
                S0 = np.concatenate([p[1] for p in parts])
                Ts = np.vstack([p[2] for p in parts])
                lines[line] = ({rk: i for i, rk in enumerate(ranks)}, base-D@c0+D@S0, D@Ms, D@Ts)
            index, off0, offK, offQ = lines[line]
            a, b = index[r['rank']], index[r['ref_rank']]
            y0.append(off0[a]-off0[b]-r['observed'])
            K.append(offK[a]-offK[b])
            Q.append(offQ[a]-offQ[b])
            scale.append(1/np.hypot(r['sigma'], interval_floor))
        for r in params:
            M, S0, T = model.parts(r['v'], r['J'])
            k = r['k']
            y0.append(S0[k]-r['observed'])
            K.append(M[k])
            Q.append(T[k])
            floor = eqq_floor if k == 0 else relative_floor*abs(r['observed'])
            scale.append(1/np.hypot(r['sigma'], floor))
        n = len(scale)
        self.y0 = np.array(y0) if n else np.zeros(0)
        self.K = np.array(K) if n else np.zeros((0, model.nb))
        self.Q = np.array(Q) if n else np.zeros((0, len(TRIG)))
        self.scale = np.array(scale)

    def residuals(self, b, w):
        return (self.y0+self.K@b+self.Q@w)*self.scale


# ---------------------------------------------------------------------------------------------------------------
# data

def load_data(frozen=True):
    rows, params = hyperfine_rows.release() if frozen else (hyperfine_rows.intervals(),
                                                            hyperfine_rows.chen_rows()+hyperfine_rows.pique_rows())
    for r in params:
        r['k'] = PARAMS.index(r['parameter'])
    return rows, params


def needed_levels(rows, params):
    out = set()
    for r in rows:
        out |= {(r['line'].v_upper, j) for j in neighbours(upper_J(r['line']))}
    out |= {(r['v'], r['J']) for r in params}
    return sorted(out)


BLOCKS = ((42, 44), (45, 46), (47, 48), (49, 51), (52, 53), (55, 58), (59, 61), (62, 63), (65, 70), (71, 75))


def block(v):
    for lo, hi in BLOCKS:
        if lo <= v <= hi:
            return f'{lo}-{hi}'
    return 'core'


def score(rows, params, constants_fn, levels):
    """Exact per-group scores: intervals rms and median |residual| (kHz) per v' (core training/future pooled);
    constants rms and median |r|/sigma per parameter and v' block."""
    table, cache = {}, {}
    for r in rows:
        line = r['line']
        if line not in cache:
            cache[line] = line_offsets(levels, line, {j: constants_fn(line.v_upper, j)
                                                      for j in neighbours(upper_J(line))})
        o = cache[line]
        key = (f"core {'future' if r['future'] else 'training'}" if r['inside'] else f"intervals v'{r['line'].v_upper}")
        table.setdefault(key, []).append(o[r['rank']]-o[r['ref_rank']]-r['observed'])
    out = {k: dict(n=len(v), rms_kHz=float(np.sqrt(np.mean(np.square(v)))*1e3),
                   median_kHz=float(np.median(np.abs(v))*1e3)) for k, v in sorted(table.items())}
    for r in params:
        key = f"{r['source']} {r['parameter']} {block(r['v'])}"
        value = constants_fn(r['v'], r['J'])[r['k']]
        out.setdefault(key, []).append(dict(v=r['v'], J=r['J'], resid=float(value-r['observed']),
                                            sigma=float(r['sigma'])))
    for k, v in out.items():
        if isinstance(v, list):
            res = np.array([x['resid'] for x in v])
            sig = np.array([x['sigma'] for x in v])
            out[k] = dict(n=len(v), rms=float(np.sqrt(np.mean(res**2))), median_z=float(np.median(abs(res)/sig)),
                          rows=v)
    return out


# ---------------------------------------------------------------------------------------------------------------
# fitting: background from the core (robust), angles from the out-of-core training rows

CORE_FLOOR = .005            # MHz
CORE_F_SCALE = 3.            # soft-L1 knee, in units of the weighted residual
OUT_FLOOR = .02              # MHz: model-error floor of out-of-core intervals
RELATIVE_FLOOR = .03         # of |value|, constants C, d, delta
EQQ_FLOOR = .3               # MHz


def solve_background(core, w, b0):
    """Robust (soft-L1) on measured core intervals; plain linear least squares on the anchored ones."""
    if not getattr(core, 'robust', True):
        A = core.K*core.scale[:, None]
        y = -(core.y0+core.Q@w)*core.scale
        return np.linalg.lstsq(A, y, rcond=None)[0]
    return least_squares(lambda b: core.residuals(b, w), b0, loss='soft_l1', f_scale=CORE_F_SCALE,
                         method='trf', x_scale='jac').x


def anchored(rows, levels):
    """Core rows with 'observed' replaced by the core fit's own prediction (its B radial functions and X constants):
    fitting the background to them keeps the core's hyperfine as it is, with S2 inside."""
    fn = core_constants(levels, continued=False)
    cache, out = {}, []
    for r in rows:
        line = r['line']
        if line not in cache:
            cache[line] = line_offsets(levels, line, {j: fn(line.v_upper, j) for j in neighbours(upper_J(line))})
        o = cache[line]
        out.append(dict(r, observed=o[r['rank']]-o[r['ref_rank']], sigma=.001))
    return out


def fit(model, rows, params, theta, fit_theta, holdout=(), max_outer=5, anchor=True, fit_shape=False, keys=None):
    """Background b (robust, core training intervals) and, if fit_theta, theta_u, theta_g (out-of-core training
    rows not in a held-out v' block). Renews the linearization until the offsets stop moving."""
    def held(v):
        return any(lo <= v <= hi for lo, hi in holdout)

    core_rows = [r for r in rows if r['inside'] and not r['future']]
    if anchor:
        core_rows = anchored(core_rows, model.levels)
    out_rows = [r for r in rows if not r['inside'] and not held(r['line'].v_upper)]
    out_params = [r for r in params if not held(r['v'])]
    lines = sorted({r['line'] for r in core_rows+out_rows}, key=str)
    b, theta = model.start(), dict(theta, gamma=theta.get('gamma', 0.))
    linear, last = {}, None
    for outer in range(max_outer):
        fn = lambda v, J: model.constants(b, theta, v, J)        # noqa: E731
        linear = {line: linearize(model.levels, line, fn) for line in lines}
        core = System(model, linear, core_rows, [], 0. if anchor else CORE_FLOOR)
        out = System(model, linear, out_rows, out_params, OUT_FLOOR, RELATIVE_FLOOR, EQQ_FLOOR)
        core.robust = not anchor
        b = solve_background(core, trig(theta), b)
        if keys is None:
            keys = (['u', 'g'] if fit_theta and model.use_s2 else [])+(['gamma'] if fit_shape else [])
        if keys:
            state = dict(b=b)

            def outer_res(t):
                w = trig(dict(theta, **dict(zip(keys, t))))
                state['b'] = solve_background(core, w, state['b'])
                return out.residuals(state['b'], w)
            # the angles are periodic (180 deg) with several local minima: start each from a 5-degree grid
            t0 = [theta[k] for k in keys]
            for i, k in enumerate(keys):
                if k in ('u', 'g'):
                    grid = np.arange(-90., 90., 5.)
                    costs = [np.sum(outer_res([*t0[:i], th, *t0[i+1:]])**2) for th in grid]
                    t0[i] = float(grid[int(np.argmin(costs))])
            t = least_squares(outer_res, t0, diff_step=1e-3).x
            theta = dict(theta, **{k: float(x) for k, x in zip(keys, t)})
            b = solve_background(core, trig(theta), state['b'])
        offsets = np.concatenate([linearize_offsets(model, linear, line, b, theta) for line in lines])
        change = np.inf if last is None else np.max(np.abs(offsets-last))
        print(f'outer {outer}: theta {theta}, core cost {0.5*np.sum(core.residuals(b, trig(theta))**2):.4g}, '
              f'out cost {0.5*np.sum(out.residuals(b, trig(theta))**2):.4g}, max offset change {change*1e3:.3g} kHz',
              flush=True)
        model.training = dict(core_cost=float(0.5*np.sum(core.residuals(b, trig(theta))**2)),
                              out_cost=float(0.5*np.sum(out.residuals(b, trig(theta))**2)), n_out=len(out.scale))
        if change < 1e-4:
            break
        last = offsets
    return b, theta


def linearize_offsets(model, linear, line, b, theta):
    """Exact offsets of a line at (b, theta), for the convergence test of the linearization."""
    Js = linear[line][0]
    o = line_offsets(model.levels, line, {j: model.constants(b, theta, line.v_upper, j) for j in Js})
    return np.array([o[r] for r in linear[line][2]])


def core_constants(levels, continued=True):
    """(v, J) -> the core fit's B radial functions (no S2): continued past KNOT as here, or as the core
    evaluates them (the polynomial everywhere)."""
    def fn(v, J):
        m = levels.moments(v, J) if continued else levels.plain_moments(v, J)
        return np.array([levels.core[p][:ORDERS[p]]@m[:ORDERS[p]] for p in PARAMS])
    return fn


def parse_overrides(items):
    out = {}
    for item in items or ():
        state, rest = item.split(':', 1)
        key, value = rest.split('=')
        out.setdefault(state, {})[key] = float(value)
    return out


#: the release's fitted out-of-core terms and curve overrides (second-order fit "scan2_g0.65_u0")
FIT_KEYS = ("u", "gamma")
FIT_OVERRIDES = {"1g(1Pig)": dict(C5=-155350.0), "(3)0u-": dict(C5=0.)}


def potential():
    """(curves, values) of the B continuation the hyperfine stage was fitted on: frozen with the release, or the
    current continuation when a rebuild names one (I2SPEC_FIT_CONTINUATION)."""
    if os.environ.get("I2SPEC_FIT_CONTINUATION"):
        from .continuation import release_spec
        frozen = {f"{s}:{n}": z for s, v in release_spec()["values"].items() for n, z in v.items()}
    else:
        frozen = json.loads((training.TRAINING/training.RELEASE/"hyperfine_potential.json").read_text())["parameters"]
    curves = build_curves()
    values = {s: dict(c.values) for s, c in curves.items()}
    for name, z in frozen.items():
        s, n = name.split(":")
        if s in values:
            values[s][n] = z
    return curves, values
