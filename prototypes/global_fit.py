"""i2spec's own direct-potential fit: MLR potentials for X and B, from Gerstenkorn & Luc to our data.

Nothing here starts from the Hannover parameters. The start is the RKR curve of the Gerstenkorn & Luc
1985 Dunham constants (prototypes/mlr_from_rkr.py), the long-range coefficients are the published
measurements of Bacis 1986 (X) and Gerstenkorn, Luc & Amiot 1985 (B), and the data are ours:

  atlas     67 729 line centres from two FTS atlases (data/atlas_lines), 10-25 MHz, v' 2-58, v'' 0-8,
            each atlas with a fitted calibration offset
  precision  the comb-referenced and interferometric sets (data/observations), 5 kHz - 2 MHz
  dunham     term values of the Gerstenkorn & Luc constants themselves, as a *stage*, to bring the
            potentials into the right basin before the real data are shown to them

Parameters: X (Re, beta_0..), B (Te, Re, beta_0..), one calibration offset per atlas. Derivatives are
Hellmann-Feynman, <psi|dV/dtheta|psi> on the solver's own quadrature grid, so an iteration costs one
eigen-solve per J per state and no extra solves for the 30-50 parameters.

An iteration costs one eigen-solve per J per state. Three things make that affordable: LAPACK's simple
and divide-and-conquer drivers instead of the subset one (2.0x and 1.7x), interpolating the X levels
and their derivatives in J(J+1) from 16 node solves instead of 198 (exact to 1e-5 MHz below v'' = 43,
with the few J above it solved directly), and a pool of spawned workers that build their own matrices
and send back contracted derivative tables. Together: 33.9 s an iteration down to 6.6 s, with results
identical to the serial exact computation.

usage: global_fit.py [--stage=dunham|data|both] [--nx=10] [--nb=14] [--jmax=213] [--max-nfev=40]
                     [--tag=] [--start=prototypes/out/mlr_from_rkr.json] [--loss=] [--floors=]
       I2SPEC_FIT_WORKERS=1 turns the pool off.
"""
from __future__ import annotations

import csv
import json
import multiprocessing as mp
import os
import sys
import time
from pathlib import Path

import numpy as np
from scipy.interpolate import make_interp_spline
from scipy.optimize import least_squares

sys.path.insert(0, str(Path(__file__).parent))
from mlr_from_rkr import (B_LONG_RANGE, TAIL_WEIGHT, X_LR_KEYS, X_LR_PRIOR, X_LONG_RANGE,  # noqa: E402
                          mlr, tail_points)

from i2spec import dunham                                                     # noqa: E402
from i2spec.bspline import BSplineSolver                                      # noqa: E402
from i2spec.constants import MHZ_PER_CM, reduced_mass                         # noqa: E402
from i2spec.observations import Line, Predictor, load_all                     # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
GRIDS = {"X": dict(rmin=2.10, rmax=6.0, h=0.01, order=10, nlev=70),
         "B": dict(rmin=2.35, rmax=8.0, h=0.01, order=10, nlev=85)}
ATLASES = ("salami_ross_2005", "apo_nist_2009")
#: J values at which the X levels are solved before being splined in J(J+1) (see Fit.levels).
X_NODES = 16
#: Worker processes for the eigen-solves. Threads do not help -- LAPACK holds the GIL, 1.02x -- and
#: forking a process that has already done heavy BLAS deadlocks its children on macOS, so the pool is
#: spawned and each worker builds its own matrices from the parameter vector.
WORKERS = int(os.environ.get("I2SPEC_FIT_WORKERS", min(8, os.cpu_count() or 1)))


def _potential(state, x, nx, nb, n_lr):
    """One state's MLR from a parameter vector."""
    lr = x[3 + nx + nb:3 + nx + nb + n_lr] if n_lr else None
    return mlr("X", x[:1 + nx], nx, lr=lr) if state == "X" else mlr("B", x[1 + nx:3 + nx + nb], nb, lr=lr)


_solver_cache: dict = {}


def _build(state, x, nx, nb, n_lr):
    """The solver for one state from a parameter vector, in the parent or in a worker.

    Cached on the parameter vector: the optimiser asks for the residuals and then the Jacobian at the
    same point, and building the basis and its matrices costs about as much as thirty eigen-solves.
    """
    key = (state, np.asarray(x, dtype=float).tobytes())
    hit = _solver_cache.get(key)
    if hit is None:
        hit = BSplineSolver(_potential(state, x, nx, nb, n_lr), reduced_mass("127I2"),
                            **GRIDS[state], joins=())
        _solver_cache.clear()                      # one parameter vector at a time, in every process
        _solver_cache[key] = hit
    return hit


def _chunk(task):
    """One worker's share of the J values: levels, and if ``idx`` is given the Hellmann-Feynman
    derivatives for those parameters, contracted to (n_par, nlev) so that no wavefunction is sent back."""
    state, x, nx, nb, n_lr, js, idx, step = task
    x = np.asarray(x, dtype=float)
    s = _build(state, x, nx, nb, n_lr)
    if idx is None:
        return {J: s.levels(J) for J in js}
    wdv = []
    for i in idx:
        h = step * max(abs(x[i]), 1.0)
        xp, xm = x.copy(), x.copy()
        xp[i] += h
        xm[i] -= h
        up = _potential(state, xp, nx, nb, n_lr)
        um = _potential(state, xm, nx, nb, n_lr)
        wdv.append(s.W * ((up(s.R) - um(s.R)) / (2 * h)))        # dV/dtheta on the quadrature grid
    return {J: _derivatives(s, J, wdv) for J in js}


def _derivatives(s, J, wdv):
    e, v = s.wavefunctions(J)
    return e, np.array([(w[:, None] * v ** 2).sum(axis=0) for w in wdv])


def solve_many(state, x, nx, nb, n_lr, js, idx=None, step=1e-4, workers=WORKERS, pool=None):
    """{J: levels} or {J: (levels, dE/dtheta)}, spread over the pool when there is one."""
    task = (state, list(map(float, x)), nx, nb, n_lr, list(js), idx, step)
    if pool is None or workers <= 1 or len(js) < 4 * workers:   # small sets are not worth the dispatch
        return _chunk(task)
    tasks = [(state, task[1], nx, nb, n_lr, list(js)[i::workers], idx, step) for i in range(workers)]
    out = {}
    for part in pool.map(_chunk, tasks):
        out.update(part)
    return out


def make_pool(workers=WORKERS):
    """A spawn pool, created before any solving so no worker inherits a busy BLAS.

    Spawned workers import the model modules afresh, so anything that changes the model must reach
    them through the environment or the task, never through a module global set in the parent.
    """
    return None if workers <= 1 else mp.get_context("spawn").Pool(workers)
#: model-error floor added in quadrature to every observed uncertainty, cm-1 (0.5 MHz): the fit is not
#: asked to chase differences smaller than the hyperfine offsets we hold fixed.
FLOOR = 0.5 / MHZ_PER_CM


def atlas_rows(max_sigma=1e9, min_depth=0.0):
    """Atlas line centres, optionally cut on quality. The cut matters: lines seen by one atlas only,
    with a fit sigma above 22 MHz, sit 265 MHz from any model, while those at sigma <= 14 MHz sit at
    40-67 MHz. sigma is the discriminator, and it is already in the file."""
    rows = []
    for k, name in enumerate(ATLASES):
        for r in csv.DictReader(open(ROOT / "data/atlas_lines" / f"{name}.csv")):
            sigma, depth = float(r["uncertainty"]), float(r["note"].split()[1].rstrip(";"))
            if sigma * MHZ_PER_CM > max_sigma or depth < min_depth:
                continue
            rows.append((Line.parse(r["line"]), float(r["value"]), sigma, k))
    return rows


def precision_rows(predictor, v_lower_max=54):
    """Absolute positions of the precision sets, hyperfine-free, with the hyperfine offset removed
    from the model we are replacing (the hyperfine Hamiltonian is not part of this fit).

    v'' = 48, 53 and 54 (matyugin2012, nesterenko2019) are the only measurements above v'' = 17 and the
    only thing that holds the X potential together through the gap: left out, our X is 12 THz wrong there.
    """
    rows = []
    for ds in load_all():
        scale = MHZ_PER_CM if ds.unit == "cm-1" else 1.0
        for o in ds.observations:
            if o.line.isotopologue != "127I2" or o.kind != "frequency" or o.line.v_lower > v_lower_max:
                continue
            offset = 0.0
            if o.component is not None:
                offset = predictor.position(o.line, o.component) - predictor.position(o.line)
            nu = (o.value * scale - offset) / MHZ_PER_CM
            rows.append((o.line, nu, o.uncertainty * scale / MHZ_PER_CM, None))
    return rows


def dunham_rows(j_max, v_max_x=19, v_max_b=70, step_j=8):
    """Term values of the Gerstenkorn & Luc model as synthetic transitions, 0.002 cm-1 each (its own
    stated accuracy). Every (v', v'') pair is redundant: this is a fit to two sets of levels."""
    rows = []
    for vl in range(0, v_max_x + 1):
        for vu in range(0, v_max_b + 1, 2):
            for J in range(0, j_max + 1, step_j):
                for br in ("R", "P"):
                    if br == "P" and J == 0:
                        continue
                    nu = dunham.transition(vu, vl, J, br)
                    if 9000.0 < nu < 20050.0:
                        rows.append((Line("127I2", br, J, vu, vl), nu, 0.002, None))
    return rows


class StateFit:
    """One state's MLR against its own term values: E(v, J) - E(0, 0) of the Gerstenkorn & Luc model.

    A transition fit ties the two states together through every line and is badly conditioned before
    either potential is right; a state fit is not, because the X levels depend on the X potential alone.
    This is how the pair is brought into the basin, after which the real data see them together.
    """

    def __init__(self, state, n_beta, x0, v_max, j_max, step_j=10, sigma=0.002, t00=None, zpe_x=None,
                 absolute=()):
        self.state, self.n_beta = state, n_beta
        #: For B, one extra row ties the placement to T00 = Te + ZPE_B - ZPE_X. Without it the fit is
        #: blind to Te (every target is a difference from v = 0) while De, which is tied to Te through
        #: the atomic asymptote, drifts: the well ends up 44 cm-1 too deep and the pair 1.3 GHz out.
        self.t00, self.zpe_x = t00, zpe_x
        self.mu = reduced_mass("127I2")
        self.x0 = np.array(x0, dtype=float)
        self.js = list(range(0, j_max + 1, step_j))
        self.targets = [(v, J, dunham.energy(state, v, J) - dunham.energy(state, 0, 0))
                        for v in range(0, v_max + 1) for J in self.js]
        #: Levels the measurements imply directly, above this state's own potential minimum: the only
        #: thing that holds X together above v'' = 19, where the atlas constants stop. An MLR fitted to
        #: v'' <= 19 alone reproduces them to 0.001 cm-1 and is still 408 cm-1 low at v'' = 48.
        self.absolute = list(absolute)
        self.n_rel = len(self.targets)
        self.v = np.array([t[0] for t in self.targets] + [a[0] for a in self.absolute])
        self.j = np.array([t[1] for t in self.targets] + [a[1] for a in self.absolute])
        self.value = np.array([t[2] for t in self.targets] + [a[2] for a in self.absolute])
        self.sigma = np.array([sigma] * len(self.targets) + [a[3] for a in self.absolute])
        for J in {a[1] for a in self.absolute}:
            if J not in self.js:
                self.js.append(J)
        self.js.sort()
        if t00 is not None:
            self.value = np.append(self.value, t00)
            self.sigma = np.append(self.sigma, sigma / 4)
        rt, vt = tail_points(state, mlr(state, x0, n_beta).De)
        self.n_tail = len(rt)
        self.value = np.append(self.value, vt)
        self.sigma = np.append(self.sigma, np.full(self.n_tail, TAIL_WEIGHT))
        self.scale = np.array(([0.005] if state == "X" else [0.05, 0.005]) + [0.2] * n_beta)

    def potential(self, x):
        return mlr(self.state, x, self.n_beta)

    def solve(self, x, want_psi=False):
        s = BSplineSolver(self.potential(x), self.mu, **GRIDS[self.state], joins=())
        levels = {}
        psi = {}
        for J in self.js:
            if want_psi:
                e, v = s.wavefunctions(J)
                psi[J] = v
            else:
                e = s.levels(J)
            levels[J] = e
        return (levels, s, psi) if want_psi else levels

    def model(self, x, levels=None):
        lev = levels if levels is not None else self.solve(x)
        e = np.array([lev[J][v] for v, J in zip(self.v, self.j)])
        out = e - np.where(np.arange(len(e)) < self.n_rel, lev[0][0], 0.0)   # relative rows only
        if self.t00 is not None:
            out = np.append(out, lev[0][0] - self.zpe_x)           # T00 = E_B(0,0) - E_X(0,0)
        return np.append(out, self.tail_model(x))

    def tail_model(self, x):
        """The potential itself, where the long-range form is the potential. Without these rows an
        exponent centred on the well (R_ref = R_e) is free to do anything beyond the data and does."""
        p = self.potential(x)
        rt, _ = tail_points(self.state, p.De)
        with np.errstate(over="ignore", invalid="ignore"):
            v = p(rt) - (0.0 if self.state == "X" else x[0])
        return np.where(np.isfinite(v), v, 1e6)

    def residuals(self, x):
        try:
            px = self.potential(x)
            if px.De <= 0 or not np.isfinite(px(np.arange(2.2, 7.0, 0.05))).all():
                return np.full(len(self.value), 1e4)
            return (self.model(x) - self.value) / self.sigma
        except (ValueError, np.linalg.LinAlgError):
            return np.full(len(self.value), 1e4)

    def jacobian(self, x, step=1e-4):
        lev, s, psi = self.solve(x, want_psi=True)
        Jm = np.zeros((len(self.value), len(x)))
        for i in range(len(x)):
            h = step * max(abs(x[i]), 1.0)
            xp, xm = x.copy(), x.copy()
            xp[i] += h; xm[i] -= h
            dV = (self.potential(xp)(s.R) - self.potential(xm)(s.R)) / (2 * h)
            wdv = s.W * dV
            d = {J: (wdv[:, None] * psi[J] ** 2).sum(axis=0) for J in self.js}
            col = np.array([d[J][v] for v, J in zip(self.v, self.j)])
            col = col - np.where(np.arange(len(col)) < self.n_rel, d[0][0], 0.0)
            if self.t00 is not None:
                col = np.append(col, d[0][0])
            xp, xm = x.copy(), x.copy()
            xp[i] += h; xm[i] -= h
            col = np.append(col, (self.tail_model(xp) - self.tail_model(xm)) / (2 * h))
            Jm[:, i] = col / self.sigma
        return Jm

    def report(self, x, label):
        rt = (self.model(x) - self.value)[-self.n_tail:]
        r = (self.model(x) - self.value)[:len(self.v)]
        if self.absolute:
            ra = r[self.n_rel:]
            print(f"  {self.state} {label}: implied levels n={len(ra)} rms {np.sqrt(np.mean(ra**2)) * MHZ_PER_CM:.4g} MHz")
        r = r[:self.n_rel]
        v_rel = self.v[:self.n_rel]
        if np.abs(rt).max() > 2 * TAIL_WEIGHT:
            print(f"    (long-range tail off by up to {np.abs(rt).max():.2f} cm-1)")
        print(f"  {self.state} {label}: rms {r.std():.5f} cm-1 ({MHZ_PER_CM * np.sqrt(np.mean(r**2)):.1f} MHz), "
              f"max |{np.abs(r).max():.4f}| cm-1")
        for lo, hi in ((0, 5), (5, 20), (20, 45), (45, 71)):
            m = (v_rel >= lo) & (v_rel < hi)
            if m.any():
                print(f"      v {lo:2d}-{hi:2d}: n={m.sum():4d} rms {np.sqrt(np.mean(r[m]**2)) * MHZ_PER_CM:9.1f} MHz")
        return r


class Fit:
    def __init__(self, rows, nx, nb, start, fit_offsets=True, high_v_floor=0.0, fit_long_range=False,
                 pool=None):
        self.rows = rows
        self.nx, self.nb = nx, nb
        self.mu = reduced_mass("127I2")
        self.fit_offsets = fit_offsets
        #: per atlas: an additive offset (cm-1) and a multiplicative scale error, which is what an FTS
        #: wavenumber scale actually gets wrong -- the two atlases differ by +20 MHz at 15 000 cm-1 and
        #: +36 MHz at 18 500, a slope no constant can absorb.
        n_cal = 2 * len(ATLASES) if fit_offsets else 0
        self.n_cal, self.n_lr = n_cal, (len(X_LR_KEYS) if fit_long_range else 0)
        # a result file carries its calibration too: restarting from one should not throw it away
        cal = (list(start.get("offsets", [])) + list(start.get("scales", [])))[:n_cal] if n_cal else []
        self.x0 = np.array(start["X"]["x"][:1 + nx] + start["B"]["x"][:2 + nb]
                           + list(start.get("lr", [0.0] * len(X_LR_KEYS)))[:self.n_lr]
                           + cal + [0.0] * (n_cal - len(cal)))
        self.scale = np.array([0.005] + [0.2] * nx + [0.05, 0.005] + [0.2] * nb
                              + [X_LR_PRIOR[k] for k in X_LR_KEYS][:self.n_lr]
                              + ([1e-4] * len(ATLASES) + [1e-8] * len(ATLASES) if fit_offsets else []))
        self.needed = {"X": sorted({r[0].J_lower for r in rows}),
                       "B": sorted({r[0].J_lower + (1 if r[0].branch == "R" else -1) for r in rows})}
        self.value = np.array([r[1] for r in rows])
        self.group = np.array([-1 if r[3] is None else r[3] for r in rows])
        #: v'' > 17 is measured to kHz but starts 93 cm-1 from any potential fitted below it: at its own
        #: uncertainty those 36 rows outweigh the other 40 000 by 10^10 and the fit goes nowhere (the same
        #: wall as docs/design/fitting.md). The floor comes down in steps as the potential catches up.
        self.high_v = np.array([r[0].v_lower > 17 for r in rows])
        self.set_high_v_floor(high_v_floor, np.array([r[2] for r in rows]))
        self.vx = np.array([r[0].v_lower for r in rows])
        self.vb = np.array([r[0].v_upper for r in rows])
        self.jx = np.array([r[0].J_lower for r in rows])
        self.jb = np.array([r[0].J_lower + (1 if r[0].branch == "R" else -1) for r in rows])
        #: X may be interpolated in J(J+1); B may not (its levels reorder above v' = 30). calibrate()
        #: fills v_safe and exact_js, the levels and J the spline cannot be trusted with.
        self.nodes, self.v_safe, self.exact_js = {}, {}, {}
        if len(self.needed["X"]) > 3 * X_NODES:
            lo, hi = self.needed["X"][0], self.needed["X"][-1]
            self.nodes["X"] = sorted({int(round(v)) for v in np.linspace(lo, hi, X_NODES)})
            self.v_safe["X"], self.exact_js["X"] = int(max(self.vx)), []
        self.pool = pool
        self._cache = {}

    def bounds(self):
        """Physical bounds on every parameter, for the optimiser to respect.

        A flat penalty where the parameters stop being a molecule is no use: it has no gradient, so an
        optimiser that steps into it cannot tell which way is back, and with Jacobian scaling one did
        walk out to R_e = 1e98. Bounds keep it in the region where the model means something.
        """
        lo = [2.5] + [-1e4] * self.nx + [15600.0, 2.8] + [-1e4] * self.nb
        hi = [2.9] + [+1e4] * self.nx + [15900.0, 3.3] + [+1e4] * self.nb
        lo += [-5 * X_LR_PRIOR[k] for k in X_LR_KEYS][:self.n_lr]
        hi += [+5 * X_LR_PRIOR[k] for k in X_LR_KEYS][:self.n_lr]
        if self.fit_offsets:
            n_a = len(ATLASES)
            lo += [-0.02] * n_a + [-3e-6] * n_a          # 600 MHz of offset, 3 ppm of scale
            hi += [+0.02] * n_a + [+3e-6] * n_a
        return np.array(lo), np.array(hi)

    def set_high_v_floor(self, floor_mhz, raw=None):
        if raw is not None:
            self._raw_sigma = raw
        self.sigma = np.hypot(self._raw_sigma, FLOOR)
        if floor_mhz:
            self.sigma = np.where(self.high_v, np.hypot(self.sigma, floor_mhz / MHZ_PER_CM), self.sigma)

    def potentials(self, x):
        lr = x[3 + self.nx + self.nb:3 + self.nx + self.nb + self.n_lr] if self.n_lr else None
        px = mlr("X", x[:1 + self.nx], self.nx, lr=lr)
        pb = mlr("B", x[1 + self.nx:3 + self.nx + self.nb], self.nb, lr=lr)
        return px, pb

    def calibrate(self, x, tol=0.05 / MHZ_PER_CM):
        """Where interpolation in J(J+1) stops being exact, and which J must therefore be solved.

        The levels are smooth in J(J+1) until they approach the top of the well, where they reorder as
        they cross the centrifugal barrier. This solves at a few J away from the nodes, finds the first
        level that the spline misses by more than ``tol``, and marks every J carrying a row from above
        it for an exact solve. For the potentials seen here that is v'' = 43-46 and six J out of 198.
        """
        if not self.valid(x):
            raise ValueError("cannot calibrate at parameters that are not a pair of potentials")
        for state, nodes in self.nodes.items():
            s = _build(state, x, self.nx, self.nb, self.n_lr)
            table = np.array([s.levels(J) for J in nodes])
            spline = make_interp_spline(np.array(nodes, float) * (np.array(nodes) + 1.0), table, k=5, axis=0)
            probe = [J for J in self.needed[state] if J not in nodes]
            probe = probe[::max(1, len(probe) // 5)][:5]
            worst = np.zeros(table.shape[1])
            for J in probe:
                worst = np.maximum(worst, np.abs(spline(J * (J + 1.0)) - s.levels(J)))
            bad = np.nonzero(worst > tol)[0]
            self.v_safe[state] = int(bad[0]) - 1 if len(bad) else table.shape[1] - 1
            v_rows, j_rows = (self.vx, self.jx) if state == "X" else (self.vb, self.jb)
            self.exact_js[state] = sorted(set(j_rows[v_rows > self.v_safe[state]].tolist()))
            print(f"  {state}: {len(nodes)} node solves cover v <= {self.v_safe[state]}; "
                  f"{len(self.exact_js[state])} of {len(self.needed[state])} J solved exactly above it",
                  flush=True)
        self._cache = {}

    def levels(self, x, want_psi=False):
        """{state: {J: energies}}, and with want_psi the Hellmann-Feynman derivative of every level
        with respect to every parameter, contracted to (n_par, nlev) per J.

        X is solved at X_NODES values of J and splined in J(J+1) -- exact to 10^-5 MHz below v'' = 43
        -- with the few J that carry a higher level solved directly. B is solved at every J it needs:
        above v' = 30 its levels reorder with J, and 195 of its 199 J carry such a line. The result is
        cached on the parameter vector, because the optimiser asks for the residuals and the Jacobian
        at the same point.
        """
        key = (x.tobytes(), want_psi)
        if key in self._cache:
            return self._cache[key]
        out, deriv = {}, {}
        lr_slice = list(range(3 + self.nx + self.nb, 3 + self.nx + self.nb + self.n_lr))
        for state, off, n in (("X", 0, 1 + self.nx), ("B", 1 + self.nx, 2 + self.nb)):
            idx = (list(range(off, off + n)) + lr_slice) if want_psi else None
            js, nodes = self.needed[state], self.nodes.get(state)
            solve_js = list(js) if nodes is None else \
                sorted(set(nodes) | set(self.exact_js.get(state, ())))
            got = solve_many(state, x, self.nx, self.nb, self.n_lr, solve_js, idx, pool=self.pool)
            take = (lambda J: got[J][0]) if want_psi else (lambda J: got[J])
            if nodes is None:
                out[state] = {J: take(J) for J in js}
                if want_psi:
                    deriv[state] = {J: got[J][1] for J in js}
                continue
            xn = np.array(nodes, float) * (np.array(nodes) + 1.0)
            spline = make_interp_spline(xn, np.array([take(J) for J in nodes]), k=5, axis=0)
            out[state] = {J: spline(J * (J + 1.0)) for J in js}
            if want_psi:
                dspline = make_interp_spline(xn, np.array([got[J][1] for J in nodes]), k=5, axis=0)
                deriv[state] = {J: dspline(J * (J + 1.0)) for J in js}
            for J in self.exact_js.get(state, ()):
                out[state][J] = take(J)
                if want_psi:
                    deriv[state][J] = got[J][1]
        result = (out, deriv) if want_psi else out
        # both keys: the optimiser asks for the residuals and the Jacobian at the same point
        self._cache = {key: result, (x.tobytes(), False): out}
        return result

    def align(self, x):
        """Place the B minimum so that T00 = E_B(0,0) - E_X(0,0) matches Gerstenkorn & Luc.

        In the MLR the B state's Te is tied to its well depth through the atomic asymptote, so the fit
        of the RKR *shape* leaves the absolute placement free: without this the start is 10 cm-1 out and
        no trust region can recover it.
        """
        x = np.array(x, dtype=float)
        i_te = 1 + self.nx
        px, pb = self.potentials(x)
        sx = BSplineSolver(px, self.mu, **GRIDS["X"], joins=())
        sb = BSplineSolver(pb, self.mu, **GRIDS["B"], joins=())
        t00 = sb.energy(0, 0) - sx.energy(0, 0)
        x[i_te] += t00 - dunham.constants()["T00"]
        return x

    def model(self, x, levels=None):
        if levels is None and not self.valid(x):
            raise ValueError("these parameters are not a pair of potentials")
        lev = levels if levels is not None else self.levels(x)
        e_b = np.array([lev["B"][J][v] for v, J in zip(self.vb, self.jb)])
        e_x = np.array([lev["X"][J][v] for v, J in zip(self.vx, self.jx)])
        nu = e_b - e_x
        if self.fit_offsets:
            n_a = len(ATLASES)
            for k in range(n_a):
                m = self.group == k
                nu = nu + m * (x[-2 * n_a + k] + x[-n_a + k] * self.value)
        return nu

    def valid(self, x):
        """Is this parameter vector a pair of real potentials? The optimiser will try vectors that are
        not: a B state above its own asymptote, or an exponent that overflows."""
        if not np.isfinite(x).all():
            return False
        try:
            px, pb = self.potentials(x)
        except (ValueError, ZeroDivisionError, OverflowError, FloatingPointError):
            return False
        if px.De <= 0 or pb.De <= 0 or not (2.0 < px.Re < 4.0) or not (2.4 < pb.Re < 5.0):
            return False
        with np.errstate(over="ignore", invalid="ignore"):
            for p, grid in ((px, np.arange(2.1, 6.0, 0.02)), (pb, np.arange(2.35, 8.0, 0.02))):
                v = p(grid)
                if not np.isfinite(v).all() or np.abs(v).max() > 1e7:
                    return False
        return True

    def prior(self, x):
        """Gaussian priors on the long-range corrections, in units of the papers' own uncertainties."""
        if not self.n_lr:
            return np.zeros(0)
        lr = x[3 + self.nx + self.nb:3 + self.nx + self.nb + self.n_lr]
        return np.array([v / X_LR_PRIOR[k] for v, k in zip(lr, X_LR_KEYS)])

    def residuals(self, x):
        n = len(self.rows) + self.n_lr
        if not self.valid(x):
            return np.full(n, 1e4)                           # a wall the optimiser can back away from
        try:
            return np.concatenate([(self.model(x) - self.value) / self.sigma, self.prior(x)])
        except (ValueError, OverflowError, FloatingPointError, np.linalg.LinAlgError):
            return np.full(n, 1e4)          # a wall, not a crash: the optimiser backs away from it

    def jacobian(self, x, step=1e-4):
        """Hellmann-Feynman: dE/dtheta = <psi|dV/dtheta|psi> on the quadrature grid."""
        if not self.valid(x):
            return np.zeros((len(self.rows) + self.n_lr, len(x)))
        try:
            _, deriv = self.levels(x, want_psi=True)
        except (ValueError, OverflowError, FloatingPointError, np.linalg.LinAlgError):
            return np.zeros((len(self.rows) + self.n_lr, len(x)))
        n_rows = len(self.rows)
        Jm = np.zeros((n_rows + self.n_lr, len(x)))
        lr_slice = list(range(3 + self.nx + self.nb, 3 + self.nx + self.nb + self.n_lr))
        for state, off, n, sign in (("X", 0, 1 + self.nx, -1.0), ("B", 1 + self.nx, 2 + self.nb, +1.0)):
            idx = list(range(off, off + n)) + lr_slice
            v_, j_ = (self.vb, self.jb) if state == "B" else (self.vx, self.jx)
            d = deriv[state]
            block = np.array([d[J][:, v] for v, J in zip(v_, j_)])          # (rows, len(idx))
            for k, i in enumerate(idx):
                Jm[:n_rows, i] += sign * block[:, k] / self.sigma
        for k, i in enumerate(lr_slice):
            Jm[n_rows + k, i] = 1.0 / X_LR_PRIOR[X_LR_KEYS[k]]
        if self.fit_offsets:
            n_a = len(ATLASES)
            for k in range(n_a):
                m = (self.group == k).astype(float)
                n_rows = len(self.rows)
                Jm[:n_rows, -2 * n_a + k] = m / self.sigma
                Jm[:n_rows, -n_a + k] = m * self.value / self.sigma
        return Jm

    def report(self, x, label=""):
        r = (self.model(x) - self.value) * MHZ_PER_CM
        if self.n_lr:
            lr = x[3 + self.nx + self.nb:3 + self.nx + self.nb + self.n_lr]
            print("     long range: " + ", ".join(f"{k} {100 * v:+.2f}% ({v / X_LR_PRIOR[k]:+.1f} sigma)"
                                                  for k, v in zip(X_LR_KEYS, lr)))
        w = r / (self.sigma * MHZ_PER_CM)
        print(f"  {label}: rms {np.sqrt(np.mean(r**2)):9.2f} MHz, normalised {np.sqrt(np.mean(w**2)):7.2f}")
        for k, name in enumerate(ATLASES):
            m = self.group == k
            if m.any():
                extra = (f"  offset {x[-2 * len(ATLASES) + k] * MHZ_PER_CM:+8.2f} MHz"
                         f"  scale {x[-len(ATLASES) + k] * 1e9:+7.2f} ppb" if self.fit_offsets else "")
                print(f"     {name:18s} n={m.sum():6d}  rms {np.sqrt(np.mean(r[m]**2)):8.2f} MHz" + extra)
        m = self.group == -1
        if m.any():
            print(f"     {'other':18s} n={m.sum():6d}  rms {np.sqrt(np.mean(r[m]**2)):8.2f} MHz")
        return r


def save(x, f, nx, nb, tag):
    path = ROOT / f"prototypes/out/global_fit{tag}.json"
    path.write_text(json.dumps(
        {"X": dict(x=list(x[:1 + nx])), "B": dict(x=list(x[1 + nx:3 + nx + nb])),
         "offsets": list(x[-2 * len(ATLASES):-len(ATLASES)]), "scales": list(x[-len(ATLASES):]),
         "lr": list(x[3 + nx + nb:3 + nx + nb + f.n_lr]), "nx": nx, "nb": nb}, indent=1))
    return path


def main(argv):
    opts = {a.split("=")[0]: a.split("=", 1)[1] for a in argv if "=" in a}
    stage = opts.get("--stage", "both")
    nx, nb = int(opts.get("--nx", 10)), int(opts.get("--nb", 14))
    j_max = int(opts.get("--jmax", 213))
    max_nfev = int(opts.get("--max-nfev", 40))
    tag = opts.get("--tag", "")
    start = json.loads(Path(opts.get("--start", ROOT / "prototypes/out/mlr_from_rkr.json")).read_text())
    for state, n in (("X", nx), ("B", nb)):
        have = len(start[state]["x"]) - (1 if state == "X" else 2)
        if have < n:
            start[state]["x"] = start[state]["x"] + [0.0] * (n - have)

    if stage in ("dunham", "both"):
        print("stage 1: each state's MLR against its own Gerstenkorn-Luc term values")
        out, zpe_x = {}, None
        for state, n, v_max in (("X", nx, 19), ("B", nb, 60)):   # v' <= 60: R_out = 5.7 A, inside the grid, and past our data
            x0 = start[state]["x"][:(1 if state == "X" else 2) + n]
            kw = dict(t00=dunham.constants()["T00"], zpe_x=zpe_x) if state == "B" else {}
            f = StateFit(state, n, x0, v_max, j_max, **kw)
            t = time.time()
            f.report(f.x0, "start")
            res = least_squares(f.residuals, f.x0, jac=f.jacobian, x_scale=f.scale, max_nfev=max_nfev)
            f.report(res.x, f"fitted ({res.nfev} evaluations, {time.time() - t:.0f} s)")
            out[state] = dict(x=list(map(float, res.x)))
            if state == "X":
                zpe_x = float(f.solve(res.x)[0][0])
                print(f"      Re = {f.potential(res.x).Re:.6f} A, zero point {zpe_x:.4f} cm-1")
                x_first = res.x
            else:
                pb = f.potential(res.x)
                t00 = float(f.solve(res.x)[0][0]) - zpe_x
                print(f"      Re = {pb.Re:.6f} A, Te = {res.x[0]:.4f}, De = {pb.De:.3f}, "
                      f"T00 = {t00:.4f} (target {dunham.constants()['T00']:.4f})")
        # the X potential again, now also against the levels the v'' > 17 measurements imply through
        # the B potential just fitted
        pb = mlr("B", out["B"]["x"], nb)
        sb = BSplineSolver(pb, reduced_mass("127I2"), **GRIDS["B"], joins=())
        implied = []
        for line, nu, sig, _ in precision_rows(Predictor()):
            if line.v_lower <= 17:
                continue
            Ju = line.J_lower + (1 if line.branch == "R" else -1)
            implied.append((line.v_lower, line.J_lower, sb.energy(line.v_upper, Ju) - nu,
                            float(np.hypot(sig, 0.05))))          # 50 MHz: the B potential's own error
        print(f"\n  X again, with {len(implied)} levels implied by the v'' > 17 measurements")
        f = StateFit("X", nx, out["X"]["x"], 19, j_max, absolute=implied)
        f.report(f.x0, "start")
        res = least_squares(f.residuals, f.x0, jac=f.jacobian, x_scale=f.scale, max_nfev=max_nfev)
        f.report(res.x, f"fitted ({res.nfev} evaluations)")
        out["X"] = dict(x=list(map(float, res.x)))
        start = out
        (ROOT / f"prototypes/out/global_fit_dunham{tag}.json").write_text(json.dumps(out, indent=1))

    if stage in ("data", "both"):
        pred = Predictor()
        rows = atlas_rows(float(opts.get("--max-sigma", 20.0)), float(opts.get("--min-depth", 0.15))) \
            + precision_rows(pred)
        rows = [r for r in rows if r[0].J_lower <= j_max]
        print(f"\nstage 2, our data: {len(rows)} lines "
              f"({sum(1 for r in rows if r[3] is not None)} atlas, {sum(1 for r in rows if r[3] is None)} precision)")
        pool = make_pool()
        f = Fit(rows, nx, nb, start, fit_offsets=True, fit_long_range="--no-long-range" not in argv,
                pool=pool)
        f.calibrate(f.x0)
        t = time.time()
        if stage == "data":
            f.x0 = f.align(f.x0)
        f.report(f.x0, "start")
        loss = opts.get("--loss", "linear")
        floors = [float(v) for v in opts.get("--floors", "3e6,3e5,3e4,3e3,300,0").split(",")]
        x = f.x0
        # A joint transition fit only ever sees E_B - E_X, so raising both states together costs
        # nothing and the optimiser wanders along that direction. Fitting one state at a time removes
        # it: the X parameters are then constrained by the X levels alone, as in a level fit.
        rounds = int(opts.get("--alternate", "0"))
        if rounds:
            n_all = len(f.x0)
            blocks = {"X": list(range(0, 1 + nx)) + list(range(3 + nx + nb, n_all)),
                      "B": list(range(1 + nx, 3 + nx + nb)) + list(range(3 + nx + nb, n_all))}
            f.set_high_v_floor(floors[-1])
            for r in range(rounds):
                for state, keep in blocks.items():
                    mask = np.zeros(n_all, bool)
                    mask[keep] = True

                    def part(u, mask=mask, base=x):
                        full = base.copy(); full[mask] = u
                        return f.residuals(full)

                    def part_jac(u, mask=mask, base=x):
                        full = base.copy(); full[mask] = u
                        return f.jacobian(full)[:, mask]

                    res = least_squares(part, x[mask], jac=part_jac, x_scale=f.scale[mask],
                                        max_nfev=max_nfev // 2, loss=loss)
                    x = x.copy(); x[mask] = res.x
                    f.report(x, f"round {r + 1}, {state} block ({res.nfev} evaluations)")
            save(x, f, nx, nb, f"{tag}_alt")
            floors = []
        for floor in floors:
            f.set_high_v_floor(floor)
            lo, hi = f.bounds()
            res = least_squares(f.residuals, np.clip(x, lo + 1e-9, hi - 1e-9), jac=f.jacobian,
                                bounds=(lo, hi), x_scale="jac", max_nfev=max_nfev, loss=loss,
                                f_scale=float(opts.get("--f-scale", 5.0)))
            x = res.x
            f.report(x, f"floor {floor:8.1f} MHz on v''>17 ({res.nfev} evaluations)")
            hv = ((f.model(x) - f.value) * MHZ_PER_CM)[f.high_v]
            print(f"       v''>17 rows: rms {np.sqrt(np.mean(hv**2)):.3g} MHz")
            # every stage is saved: the last one is not always the best, because dropping the floor
            # to nothing lets 36 rows the model cannot reach spend 40 000 that it can
            save(x, f, nx, nb, f"{tag}_floor{floor:g}")
        res.x = x
        print(f"  {time.time() - t:.0f} s")
        if pool is not None:
            pool.close()
        save(res.x, f, nx, nb, tag)


if __name__ == "__main__":
    main(sys.argv[1:])
