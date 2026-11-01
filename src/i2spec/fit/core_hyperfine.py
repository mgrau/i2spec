"""The core's B hyperfine functions, fitted to the intervals between components of one line.

Each B coupling (eqQ, C, d, delta) is the expectation of one polynomial in z = ((R^2-Re^2)/(R^2+Re^2))/0.35 over
the level's own wavefunction from the core potentials; the X couplings are the inherited formulae of Bodermann,
Knoeckel and Tiemann (hfs_params.x_state). The angular Hamiltonian, with J +- 2 mixing, is the production one
(hyperfine._block_coefficients); eigenvalue derivatives are analytic. Fitted with the potentials held fixed, to
the training intervals (weights 1/sqrt(sigma^2 + 0.025^2 MHz^2)). Nothing is indexed by line or level.
"""
from __future__ import annotations

from types import SimpleNamespace

import numpy as np
from numpy.polynomial.polynomial import polyvander
from scipy.linalg import eigh
from scipy.optimize import least_squares

from .. import hfs_params as hp
from ..constants import ISOTOPOLOGUES, MHZ_PER_CM, NUCLEAR_SPIN
from ..hyperfine import HyperfineLevel, HyperfineLevels, _block_coefficients, line_components
from ..observations import component_rank
from . import hyperfine_rows
from .levels import Levels
from .rows import EXCLUSION

PARAMETERS = ("C", "delta", "eqQ", "d")          # the order of the angular coefficient matrices
RADIAL_SCALE = {"X": .15, "B": .35}
ORDERS = {"X": {}, "B": dict(eqQ=5, C=5, d=3, delta=1)}
FLOOR = .025                                     # MHz
DEGREE = 8                                       # radial moments computed per level


def records(exclusion=(EXCLUSION,), datasets=None):
    """The core's intra-line intervals between two labelled components: one dict per measurement (the core's rows
    of hyperfine_rows.intervals, without `inside`)."""
    return [{k: z for k, z in r.items() if k != "inside"} for r in hyperfine_rows.intervals(datasets, exclusion)
            if r["inside"]]


def line_levels(line):
    return ((line.isotopologue, "X", line.v_lower, line.J_lower),
            (line.isotopologue, "B", line.v_upper, line.J_lower+(1 if line.branch == "R" else -1)))


def moments(curves, values, lines, box_padding=.6):
    """Energies (cm-1) and radial moments <z^k>, k < DEGREE, of every level the lines' hyperfine structure
    needs (their J and J +- 2, and J = 0 of each v), from the potentials' own wavefunctions."""
    keys = {("127I2", "X", 0, 0)}
    for line in lines:
        for iso, state, v, J in line_levels(line):
            keys.add(("127I2", state, v, 0))
            keys.update((iso, state, v, j) for j in range(max(J-2, J % 2), J+3, 2))
    engine = Levels(keys, box_padding=box_padding, refine_x=True)
    energies, means = {}, {}
    for group, jobs in engine.jobs.items():
        iso, state, _ = group
        bs = engine.bases[group]
        V, f = curves[state].curves(values[state], bs.R, bs.ratio)
        vm, rm = bs.project(V), bs.project(bs.c_rot*f)
        re = np.exp(values[state]["logRe"])
        z = ((bs.R**2-re**2)/(bs.R**2+re**2))/RADIAL_SCALE[state]
        functions = polyvander(z, DEGREE-1)
        for J, entries in jobs:
            ev, vectors = eigh(bs.T+vm+J*(J+1)*rm, bs.S, driver="gvd", check_finite=False)
            vmax = max(v for v, _ in entries)
            density = bs.W[:, None]*(bs.B@vectors[:, :vmax+1])**2
            mean = density.T@functions
            for v, _ in entries:
                energies[iso, state, v, J], means[iso, state, v, J] = float(ev[v]), mean[v]
    worst = max(abs(m[0]-1) for m in means.values())
    if worst > 1e-10:
        raise ValueError(f"radial wavefunction normalization error {worst}")
    return energies, means


class Problem:
    """Component offsets of the lines and the intervals of the records, as functions of the B coefficients."""

    def __init__(self, energies, means, records, extra_lines=(), orders=ORDERS):
        self.records, self.orders = records, orders
        self.energies, self.moments = energies, means
        self.lines = sorted({r["line"] for r in records}.union(extra_lines), key=str)
        self.names = [(s, p, k) for s in ("X", "B") for p in PARAMETERS for k in range(orders[s].get(p, 0))]
        self.index = {k: i for i, k in enumerate(self.names)}
        self.blocks = {key: self._blocks(key) for key in sorted({k for line in self.lines for k in line_levels(line)})}
        self.x0 = self._formula_start()
        self._last = None

    def _formula_start(self):
        """The published B formulae projected onto the radial polynomials, over the training levels."""
        x0 = np.zeros(len(self.names))
        train = {k for r in self.records if r.get("train", not r["future"]) for k in line_levels(r["line"])}
        e00 = self.energies["127I2", "X", 0, 0]
        for state in ("X", "B"):
            chosen = [k for k in sorted({("127I2", s, v, J) for _, s, v, J in train if s == state}) if k in self.moments]
            for param in PARAMETERS:
                indices = [k for k in range(self.orders[state].get(param, 0))]
                if not indices:
                    continue
                design = np.array([self.moments[k][indices] for k in chosen])
                target = []
                for _, s, v, J in chosen:
                    energy = self.energies["127I2", s, v, 0]-e00
                    p = hp.x_state(v, J, energy) if s == "X" else hp.b_state(v, J, energy)
                    target.append(getattr(p, param))
                for k, z in zip(indices, np.linalg.lstsq(design, target, rcond=1e-10)[0]):
                    x0[self.index[state, param, k]] = z
        return x0

    def _blocks(self, key):
        iso, state, v, J = key
        a, b = ISOTOPOLOGUES[iso]
        q, m, spin = hp._scales(iso)
        multipliers = dict(C=m, delta=spin, eqQ=q, d=spin)
        qr, cr = hp.nucleus_ratios(iso)
        symmetry = ("g" if state == "X" else "u") if a == b else None
        js = tuple(j for j in range(J-2, J+3, 2) if j >= 0)
        out = []
        for F, (basis, Js, coef) in _block_coefficients(js, NUCLEAR_SPIN[a], NUCLEAR_SPIN[b], symmetry, qr, cr).items():
            matrix = np.diag([(self.energies[iso, state, v, j]-self.energies[key])*MHZ_PER_CM for j in Js])
            deriv = np.zeros((len(self.names), len(basis), len(basis)))
            for ip, param in enumerate(PARAMETERS):
                if state == "X" and not self.orders["X"]:
                    ex = self.energies["127I2", "X", v, 0]-self.energies["127I2", "X", 0, 0]
                    fixed = np.array([getattr(hp.x_state(v, j, ex), param) for j in Js])
                    matrix += coef[ip]*(fixed[:, None]+fixed[None, :])*.5*multipliers[param]
                elif state == "X" and param in ("d", "delta"):
                    matrix += coef[ip]*multipliers[param]*(hp.D_X if param == "d" else hp.DELTA_X)
                for k in range(self.orders[state].get(param, 0)):
                    means = np.array([self.moments[iso, state, v, j][k] for j in Js])
                    deriv[self.index[state, param, k]] = coef[ip]*(means[:, None]+means[None, :])*.5*multipliers[param]
            out.append((F, basis, Js == J, matrix, deriv))
        return out

    def evaluate(self, theta):
        """(interval predictions, their Jacobian, {line: {rank: component offset MHz}})."""
        if self._last is not None and np.array_equal(theta, self._last[0]):
            return self._last[1:]
        structures = {}
        for key, blocks in self.blocks.items():
            selected = []
            for F, basis, central, base, deriv in blocks:
                ev, vec = np.linalg.eigh(base+np.einsum("k,knm->nm", theta, deriv))
                weight = np.where(central[:, None], vec**2, 0.)
                for column in np.flatnonzero(weight.sum(axis=0) > .5):
                    vector = vec[:, column]
                    level = HyperfineLevel(F, basis[int(np.argmax(weight[:, column]))][1], float(ev[column]), basis, vector)
                    selected.append((level, np.einsum("n,knm,m->k", vector, deriv, vector)))
            selected.sort(key=lambda pair: pair[0].energy)
            structures[key] = (HyperfineLevels([z[0] for z in selected]), np.array([z[1] for z in selected]))
        offsets, gradients = {}, {}
        for line in self.lines:
            lower, upper = line_levels(line)
            low, dl = structures[lower]
            up, du = structures[upper]
            components = line_components(up, low, upper[-1], lower[-1])
            offsets[line] = {component_rank(c.label): c.offset for c in components if c.label}
            gradients[line] = {component_rank(c.label): du[c.upper_level]-dl[c.lower_level] for c in components if c.label}
        pred = np.array([offsets[r["line"]][r["rank"]]-offsets[r["line"]][r["ref_rank"]] for r in self.records])
        jac = np.array([gradients[r["line"]][r["rank"]]-gradients[r["line"]][r["ref_rank"]] for r in self.records])
        self._last = (theta.copy(), pred, jac, offsets)
        return self._last[1:]

    def center(self, line):
        lower, upper = line_levels(line)
        return (self.energies[upper]-self.energies[lower])*MHZ_PER_CM


def release_coefficients(spec, names):
    """The released B coefficients as a vector over `names`."""
    coefficients = spec["hyperfine"]["coefficients"]
    return np.array([coefficients[s][p][k] for s, p, k in names])


def fit(problem, start=None, floor=FLOOR, max_nfev=600):
    """Weighted least squares on the training records (record["train"], else not record["future"])."""
    target = np.array([r["observed"] for r in problem.records])
    weight = 1/np.hypot(np.array([r["sigma"] for r in problem.records]), floor)
    train = np.array([r.get("train", not r["future"]) for r in problem.records])
    start = problem.x0.copy() if start is None else np.asarray(start, float)
    _, jac0, _ = problem.evaluate(start)
    scale = 1/np.maximum(np.linalg.norm(jac0[train]*weight[train, None], axis=0), 1e-8)

    def residual(q):
        return (problem.evaluate(start+q*scale)[0]-target)[train]*weight[train]

    def jacobian(q):
        return problem.evaluate(start+q*scale)[1][train]*weight[train, None]*scale

    result = least_squares(residual, np.zeros(len(start)), jac=jacobian, max_nfev=max_nfev, ftol=1e-10, xtol=1e-10,
                           gtol=1e-7)
    theta = start+scale*result.x
    pred, _, _ = problem.evaluate(theta)
    return SimpleNamespace(theta=theta, names=problem.names, cost=float(2*result.cost), nfev=int(result.nfev),
                           success=bool(result.success), message=str(result.message), train=train,
                           residual=pred-target)
