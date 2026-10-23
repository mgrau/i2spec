"""Level energies E(isotopologue, state, v, J) and their parameter derivatives, for every level a fit needs at once.

One B-spline basis per (isotopologue, state, grid tier) is built once. For a parameter vector the potential and
the centrifugal term are projected onto it, and each J costs one generalized eigensolve. Derivatives are
Hellmann-Feynman, <psi| dV/dp + J(J+1) d(rot)/dp |psi>, with the curves' own dV/dp and df/dp on the quadrature
grid. A J goes to the smallest tier that holds its highest v; the last tiers are long boxes with graded meshes
and keep quasi-bound levels that sit inside their centrifugal barrier (90 % of the density inside).
"""
from __future__ import annotations

import copy

import numpy as np
import scipy.sparse as sp
from scipy.linalg import eigh

from ..bspline import BSplineSolver
from ..constants import HBAR2_2U, REFERENCE_ISOTOPOLOGUE, reduced_mass

TIERS = {
    "X": [(17, dict(rmin=2.10, rmax=3.6, h=0.01, order=10, nlev=20)),
          (24, dict(rmin=2.10, rmax=4.0, h=0.01, order=10, nlev=30)),
          (33, dict(rmin=2.10, rmax=4.4, h=0.01, order=10, nlev=40)),
          (54, dict(rmin=2.10, rmax=5.0, h=0.01, order=10, nlev=60)),
          (69, dict(rmin=2.10, rmax=6.0, h=0.01, order=10, nlev=70)),
          (115, dict(rmin=2.10, rmax=40.0, h=0.01, order=10, nlev=118, mesh=((5.0, 0.02), (8.0, 0.05), (15.0, 0.1))))],
    "B": [(17, dict(rmin=2.35, rmax=4.0, h=0.01, order=10, nlev=26)),
          (33, dict(rmin=2.35, rmax=4.5, h=0.01, order=10, nlev=40)),
          (43, dict(rmin=2.35, rmax=5.5, h=0.01, order=10, nlev=55)),
          (59, dict(rmin=2.35, rmax=8.0, h=0.01, order=10, nlev=60)),
          (90, dict(rmin=2.35, rmax=40.0, h=0.01, order=10, nlev=260, mesh=((5.0, 0.02), (8.0, 0.05), (15.0, 0.1))))],
}


class _Flat:
    def __call__(self, R):
        return np.zeros_like(np.asarray(R, float))

    def adiabatic(self, R, mass_ratio=1.):
        return np.zeros_like(np.asarray(R, float))

    def nonadiabatic(self, R, mass_ratio=1.):
        return np.zeros_like(np.asarray(R, float))


class Basis:
    def __init__(self, isotopologue, state, grid):
        g = dict(grid)
        self.nlev = g.pop("nlev")
        self.mesh = "mesh" in g
        self.mu = reduced_mass(isotopologue)
        self.ratio = reduced_mass(REFERENCE_ISOTOPOLOGUE)/self.mu
        s = BSplineSolver(_Flat(), self.mu, nlev=self.nlev, joins=g.pop("joins", ()), **g)
        self.R, self.W = s.R, s.W
        self.B = sp.csr_matrix(np.where(np.abs(s._B) > 0, s._B, 0.0))
        self.S, self.T = s.S, s._T
        self.c_rot = HBAR2_2U/self.mu/self.R**2
        self.BT = self.B.T.tocsr()

    def project(self, f):
        return (self.BT@sp.diags(self.W*f)@self.B).toarray()


def tiers(joins=None, box_padding=0., refine_x=False):
    """The grid tiers with the outer box extended by `box_padding` (A), knots at each state's `joins`, and with
    `refine_x` the X tiers above v''=33 at half spacing."""
    out = copy.deepcopy(TIERS)
    for grids in out.values():
        for _, grid in grids:
            grid["rmax"] += box_padding
    if refine_x:
        for limit, grid in out["X"]:
            if limit >= 33:
                grid["h"] *= .5
                if "mesh" in grid:
                    grid["mesh"] = tuple((r, h*.5) for r, h in grid["mesh"])
    for state, radii in (joins or {}).items():
        for _, grid in out[state]:
            grid["joins"] = radii
    return out


class Levels:
    """E and dE/dp of a fixed set of levels. `cache_spectrum` keeps each J's spectrum until the parameters change,
    so that a residual-only evaluation and the Jacobian at the same point share their eigensolves."""

    def __init__(self, keys, joins=None, box_padding=0., refine_x=False, cache_spectrum=False):
        self.keys = sorted(set(keys))
        self.index = {k: i for i, k in enumerate(self.keys)}
        self.tiers = tiers(joins, box_padding, refine_x)
        need = {}
        for iso, st, v, J in self.keys:
            need.setdefault((iso, st, J), []).append(v)
        self.jobs = {}          # (iso, state, tier) -> [(J, [(v, index)])]
        for (iso, st, J), vs in need.items():
            tier = next(i for i, (limit, _) in enumerate(self.tiers[st]) if max(vs) <= limit)
            self.jobs.setdefault((iso, st, tier), []).append((J, [(v, self.index[(iso, st, v, J)]) for v in sorted(vs)]))
        self.bases = {k: Basis(k[0], k[1], self.tiers[k[1]][k[2]][1]) for k in self.jobs}
        self.cache_spectrum = cache_spectrum
        self._signature, self._spectra, self._assemblies = None, {}, {}

    def evaluate(self, curves, theta, names, jac=True):
        """curves: {state: curve}; names: [(state, name)] of the free parameters, theta their values. Returns E
        (cm-1, nan where a level is missing) and dE/dtheta (levels x parameters, cm-1)."""
        values = {s: dict(c.values) for s, c in curves.items()}
        for (s, n), x in zip(names, theta):
            values[s][n] = float(x)
        if self.cache_spectrum:
            signature = tuple((s, id(curves[s]), tuple(sorted(v.items()))) for s, v in sorted(values.items()))
            if signature != self._signature:
                self._signature, self._spectra, self._assemblies = signature, {}, {}
        E = np.full(len(self.keys), np.nan)
        D = np.zeros((len(self.keys), len(names))) if jac else None
        for key, jobs in self.jobs.items():
            _, state, _ = key
            bs, curve = self.bases[key], curves[state]
            if self.cache_spectrum and key in self._assemblies:
                V, f, Vm, Rm, asym = self._assemblies[key]
            else:
                V, f = curve.curves(values[state], bs.R, bs.ratio)
                Vm, Rm = bs.project(V), bs.project(bs.c_rot*f)
                asym = curve.asymptote(values[state])
                if self.cache_spectrum:
                    self._assemblies[key] = (V, f, Vm, Rm, asym)
            columns = [i for i, (s, _) in enumerate(names) if s == state]
            if jac:
                dV, df = curve.derivatives(values[state], bs.R, bs.ratio, [names[i][1] for i in columns])
                df = df*bs.c_rot[:, None]
            for J, entries in jobs:
                y = J*(J+1)
                if self.cache_spectrum and (key, J) in self._spectra:
                    ev, C = self._spectra[key, J]
                else:
                    ev, C = eigh(bs.T+Vm+y*Rm, bs.S, driver="gvd", check_finite=False)
                    n = min(len(ev), bs.nlev)
                    ev, C = ev[:n], C[:, :n].copy()
                    if self.cache_spectrum:
                        self._spectra[key, J] = (ev, C)
                if jac or bs.mesh:
                    density = bs.W[:, None]*(bs.B@C)**2
                if bs.mesh:
                    effective = V+y*bs.c_rot*f
                    outer = bs.R > 4.5
                    if effective[outer].max() > asym:
                        barrier = bs.R[outer][np.argmax(effective[outer])]
                        inside = density[bs.R < barrier].sum(axis=0)/density.sum(axis=0)
                        keep = (ev < asym) | (inside >= .9)
                    else:
                        keep = ev < asym
                    ev, density = ev[keep], density[:, keep]
                for v, index in entries:
                    if v < len(ev):
                        E[index] = ev[v]
                        if jac and columns:
                            D[index, columns] = density[:, v]@(dV+y*df)
        return E, D
