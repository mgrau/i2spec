"""The core fit: both potentials and rotational functions, to the core rows (X v''<=17, B v'<=43, J<=200).

Objective: sum of (w_i (prediction_i - value_i))^2 over the rows, plus physical penalties and two priors.

  weights   w_i = sqrt(min(1, 1000/n_c)) / sqrt(sigma_i^2 + f_c^2), n_c the rows of the category, f_c its floor
            (0.5 MHz for comb and Doppler-free rows, 20 MHz for emission, 0 for the atlases and level constants)
  shape     on a dense grid (X 1.8-6 A, B 2.0-6 A, then to 40 A, plus the extra nodes of the shape audit): the
            slope of -ln(Morse factor) must be positive (single well; scale 1e-5), V >= Te, f(R) >= 0.2
  priors    De(X) = 12547.340 +- `dissociation_prior` cm-1, and Te + De(B) - De(X) = the atomic splitting
            +- 0.0065 cm-1. With `shared_asymptote` De(B) is instead computed from that identity (one fewer
            parameter); with `fixed_de_x` De(X) is held.

Solved by damped Gauss-Newton (optimize.damped_gauss_newton) with Hellmann-Feynman Jacobians.
"""
from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from pathlib import Path

import numpy as np
import scipy.sparse as sp

from ..constants import ATOMIC_SPLITTING, MHZ_PER_CM
from . import rows as training
from .curves import Config, Curve, initial
from .levels import Levels
from .optimize import damped_gauss_newton

ROOT = Path(__file__).resolve().parents[3]
CORE_RELEASE = ROOT / "src" / "i2spec" / "data" / "research" / "20261018" / "balanced"
FLOORS = dict(comb=.5, doppler=.5, emission=20., atlas_sr=0., atlas_apo=0., orsay1=0., orsay4=0., martin=0.,
              martin_dissoc=0.)
DE_X_CENTER = 12547.340         # cm-1, the published De(X) the prior is centered on
ATOMIC_GAP_SIGMA = .0065        # cm-1, the atomic splitting's uncertainty
EXCLUSION = "127I2 P(13) 43-0"


@dataclass
class Settings:
    """The release's core settings (its request); change one to refit a variant."""
    dissociation_prior: float = .05
    de_x_center: float = DE_X_CENTER
    shared_asymptote: bool = False
    fixed_de_x: float | None = None
    slope_scale: float = 1e-5
    category_cap: float = 1000.
    floors: dict = field(default_factory=lambda: dict(FLOORS))
    box_padding: float = .6
    max_nfev: int = 400
    ftol: float = 1e-10
    gtol: float = 1e-6


def weights(rows, floors=FLOORS, category_cap=1000.):
    counts = {c: int(np.sum(rows.category == c)) for c in set(rows.category)}
    return np.array([np.sqrt(min(1., category_cap/counts[c]))/np.hypot(s, floors[c])
                     for c, s in zip(rows.category, rows.sigma)])


def release_spec():
    return json.loads((CORE_RELEASE / "model.json").read_text())


def shape_nodes():
    """The extra radial nodes of the core's shape constraints (from the shape audit of its final fit)."""
    return json.loads((training.TRAINING / training.RELEASE / "core_shape_nodes.json").read_text())


class Problem:
    """Weighted residuals of `rows` against the two curves, with the physical penalties and priors."""

    def __init__(self, rows, config: Config, values: dict, names: list, settings: Settings, nodes: dict):
        self.rows, self.config, self.settings, self.nodes = rows, config, settings, nodes
        self.curves = {s: Curve(s, config, values[s]) for s in ("X", "B")}
        self.base_names = list(names)
        s = settings
        self.active = np.array([i for i, k in enumerate(names)
                                if not (s.shared_asymptote and k == ("B", "logDe"))
                                and not (s.fixed_de_x is not None and k == ("X", "logDe"))])
        self.names = [names[i] for i in self.active]
        self.base = np.array([values[st][n] for st, n in names])
        self.levels = Levels(rows.level_set(), box_padding=s.box_padding, cache_spectrum=True)
        rr, cc, vv = [], [], []
        for i, lv in enumerate(rows.levels):
            for k, sign in lv:
                rr.append(i)
                cc.append(self.levels.index[k])
                vv.append(sign)
        self.M = sp.csr_matrix((vv, (rr, cc)), shape=(len(rows), len(self.levels.keys)))
        self.w = weights(rows, s.floors, s.category_cap)
        self._cached_x, self._cached = None, None
        self.calls = 0

    def expand(self, x):
        """Every named parameter from the free ones, and d(theta)/dx."""
        theta = self.base.copy()
        theta[self.active] = x
        n = len(self.names)
        if self.settings.fixed_de_x is not None:
            theta[self.base_names.index(("X", "logDe"))] = np.log(self.settings.fixed_de_x)
        chain = np.zeros((len(theta), n))
        chain[self.active, np.arange(n)] = 1.
        if self.settings.shared_asymptote:
            ix, ib, it = (self.base_names.index(k) for k in (("X", "logDe"), ("B", "logDe"), ("B", "Te")))
            dx = np.exp(theta[ix])
            db = dx+ATOMIC_SPLITTING-theta[it]
            if db <= 0:
                raise ValueError("Negative B dissociation energy")
            theta[ib] = np.log(db)
            if ix in self.active:
                chain[ib, np.flatnonzero(self.active == ix)[0]] = dx/db
            chain[ib, np.flatnonzero(self.active == it)[0]] = -1/db
        return theta, chain

    def values(self, theta):
        values = {s: dict(c.values) for s, c in self.curves.items()}
        for (s, n), z in zip(self.base_names, theta):
            values[s][n] = float(z)
        return values

    def grid(self, state, re):
        R = np.unique(np.r_[np.linspace(1.8 if state == "X" else 2., 6., 700), np.geomspace(6.01, 40., 180), re])
        if self.nodes.get(state):
            R = np.unique(np.r_[R, self.nodes[state]])
        return R

    def physical(self, theta, jac=False):
        values, s = self.values(theta), self.settings
        result, blocks, active_shape = [], [], False
        for state in ("X", "B"):
            curve, val = self.curves[state], values[state]
            re = np.exp(val["logRe"])
            R = self.grid(state, re)
            V, f = curve.curves(val, R, 1.)
            signed = curve.log_slope(val, R)
            te = val.get("Te", 0.)
            result.extend(np.clip(np.minimum(signed, 0.), -1e12, 0.)/s.slope_scale)
            result.extend(np.minimum(V-te, 0.))
            result.extend(np.minimum(f-.2, 0.)*1e4)
            if jac:
                columns = [i for i, (st, _) in enumerate(self.base_names) if st == state]
                active_slope = (signed < 0) & (signed > -1e12)
                active_v, active_f = V < te, f < .2
                block = np.zeros((3*len(R), len(theta)))
                if np.any(active_slope) | np.any(active_v) | np.any(active_f):
                    active_shape = True
                    ns = [self.base_names[i][1] for i in columns]
                    dv, df = curve.derivatives(val, R, 1., ns)
                    _, ds = curve.log_slope(val, R, ns)
                    block[:len(R), columns] = ds*active_slope[:, None]/s.slope_scale
                    if "Te" in ns:
                        dv[:, ns.index("Te")] -= 1.
                    block[len(R):2*len(R), columns] = dv*active_v[:, None]
                    block[2*len(R):, columns] = df*active_f[:, None]*1e4
                blocks.append(block)
        if s.dissociation_prior:
            dx, db = np.exp(values["X"]["logDe"]), np.exp(values["B"]["logDe"])
            result.append((dx-s.de_x_center)/s.dissociation_prior)
            if jac:
                row = np.zeros((1, len(theta)))
                row[0, self.base_names.index(("X", "logDe"))] = dx/s.dissociation_prior
                blocks.append(row)
            if not s.shared_asymptote:
                result.append((values["B"]["Te"]+db-dx-ATOMIC_SPLITTING)/ATOMIC_GAP_SIGMA)
                if jac:
                    row = np.zeros((1, len(theta)))
                    for name, z in [(("X", "logDe"), -dx), (("B", "logDe"), db), (("B", "Te"), 1.)]:
                        row[0, self.base_names.index(name)] = z/ATOMIC_GAP_SIGMA
                    blocks.append(row)
        result = np.nan_to_num(np.asarray(result), nan=1e12, posinf=1e12, neginf=-1e12)
        if not jac:
            return result
        derivative = np.vstack(blocks) if blocks else np.zeros((0, len(theta)))
        if active_shape:
            # the grid holds Re and moves with it: full finite differences for the two logRe columns
            for i, (_, name) in enumerate(self.base_names):
                if name == "logRe":
                    h = max(abs(theta[i])*1e-7, 1e-6)
                    plus, minus = theta.copy(), theta.copy()
                    plus[i] += h
                    minus[i] -= h
                    derivative[:, i] = (self.physical(plus)-self.physical(minus))/(2*h)
        return result, derivative

    def evaluate(self, x, jac=True):
        if self._cached_x is not None and np.array_equal(x, self._cached_x) and (not jac or self._cached[1] is not None):
            return self._cached
        theta, chain = self.expand(x)
        E, D = self.levels.evaluate(self.curves, theta, self.base_names, jac=jac)
        pred = self.M@(E*MHZ_PER_CM)
        residual = (pred-self.rows.value)*self.w
        missing = ~np.isfinite(residual)
        residual[missing] = 1e8
        if jac:
            J = self.M@(D@chain*MHZ_PER_CM)
            J[missing] = 0.
            J *= self.w[:, None]
            physical, base_jac = self.physical(theta, jac=True)
            J = np.vstack((J, base_jac@chain))
        else:
            J = None
            physical = self.physical(theta)
        residual = np.r_[residual, physical]
        self.calls += 1
        self._cached_x, self._cached = x.copy(), (residual, J, pred)
        return self._cached

    def fun(self, x):
        return self.evaluate(x, jac=False)[0]

    def cost(self, x):
        r = self.fun(x)
        n = len(self.rows)
        return dict(cost=float(r@r), data=float(r[:n]@r[:n]), physical=float(r[n:]@r[n:]))


def start(names, spec=None):
    """The released core coefficients as a parameter vector over `names`."""
    spec = spec or release_spec()
    return np.array([spec["values"][s][n] for s, n in names])


def fit(rows=None, settings: Settings | None = None, theta=None, verbose=False):
    """Fit the core. Defaults: the release's core training rows, settings and coefficients (as the start).
    Returns the problem and the solver's result (result.theta: every named coefficient)."""
    spec = release_spec()
    config = Config.from_spec(spec["config"])
    settings = settings or Settings()
    rows = rows if rows is not None else training.release("core")
    values, names, lower, upper, _ = initial(config)
    prob = Problem(rows, config, values, names, settings, shape_nodes())
    prob.base = start(names, spec) if theta is None else np.asarray(theta, float)
    x = np.clip(prob.base[prob.active], lower[prob.active]+1e-9, upper[prob.active]-1e-9)
    result = damped_gauss_newton(prob, x, lower[prob.active], upper[prob.active], max_nfev=settings.max_nfev,
                                 ftol=settings.ftol, gtol=settings.gtol, lazy_jacobian=True)
    result.theta, _ = prob.expand(result.x)
    result.names = names
    result.settings = asdict(settings)
    return prob, result
