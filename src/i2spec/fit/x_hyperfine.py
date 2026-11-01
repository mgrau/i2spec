"""X-state eqQ and C past the core: radial functions with a separated-atom limit, in place of BKT02's polynomial
(eqQ) and pole (C) extrapolations.

The core takes the X hyperfine constants from Bodermann, Knoeckel & Tiemann 2002 (BKT02): eqQ_X a polynomial in
u = v + 1/2 and x = J(J + 1), C_X = c0 + c1 u + c2/(E - 12340 cm-1). Past the core both are extrapolations, and the
pole sends C_X through infinity near v'' = 95. Bacis et al., JCP 73, 2641 (1980) (data/x_hyperfine/bacis1980)
measured eqQ_X up to v'' = 88 and quote precise values of others at v'' = 42, 62, 79, 81, 83: BKT02 holds to
v'' ~ 60 and is off by 180 MHz at v'' = 62 and 300-700 MHz at v'' = 73-88.

Model (levels of the hyperfine stage's continuation, on the fit engine's X grid):
  eqQ_X(R) = P(z) for R <= knot, P a polynomial in the core's X radial variable z (linear in z below INNER, where
             only the inner turning points of high v'' reach), and
             L + (P(z_k) - L) exp(-((R - knot)/lam)^p) beyond: the LCAO -> separated-atom transition Bacis et al.
             describe at 3.5-4.5 A, reaching L = -eQq(2P3/2) = -1146.4 MHz (their Table I).
  P is fitted to BKT02 at the core levels (v'' <= 17, J <= 200) together with lam, p to the measured values above
  the core; inside the core the packaged constants stay BKT02 exactly (constants() switches at v'' = 17).
  Bacis' v'' >= 89 values are effective (1g-perturbed, C held at 0) and are left out; their Table IV values
  (v'' 26-60) rested on eqQ_B(v' = 32) = -533 MHz and are moved by today's -544.75 (Hong et al. 2001) + 533.
  C_X(R) = polynomial in z fitted to BKT02 at the core, continued linearly in z (bounded: z -> 1/scale) with the
  slope fitted to the one measurement above the core, 59 kHz at v'' = 83 (Bacis' ref. 15). UNVALIDATED: near
  dissociation C_X grows through the X-1g interaction (Martin et al. 1983), which a smooth form cannot hold.

"""
from __future__ import annotations

import csv

import numpy as np
from scipy.optimize import least_squares

from ..bspline import BSplineSolver
from ..constants import reduced_mass
from . import rows as training
from . import second_order as so
from .continuation import CORE_J, CORE_V, SolverPotential
from .core import release_spec as core_spec
from .hyperfine_beyond import potential

KNOT, ORDER = 3.4, 5                                # A, and the polynomial order of the release

GRID = dict(rmin=2.1, rmax=40., h=.01, order=10, mesh=((5., .02), (8., .05), (15., .1)))   # the fit engine's X tier
NLEV = 118                                         # (h = 0.02 A puts a spurious level between v'' = 20 and 40)
INNER = 2.45                                       # A: below, the core polynomial continues linearly in z
EQQ_LIMIT = -so.ATOMIC['B32']                       # -eQq(2P3/2), MHz
EQQB_SHIFT_TABLE_IV = -544.751+533.                 # Hong 2001 eqQ_B(v'=32, J'=57) minus Bacis' assumed value
C_V83 = (83, 59., 10.)                              # v'', C_X kHz (BKT02's unit), assumed 1 sigma (Bacis, ref. 15)
ANCHOR_SIGMA = .005                                 # MHz: weight of the BKT02 core anchor


class XLevels:
    def __init__(self):
        curves, values = potential()
        self.solver = BSplineSolver(SolverPotential(curves['X'], values['X']), reduced_mass('127I2'),
                                    joins=curves['X'].joins, nlev=NLEV, **GRID)
        spec = core_spec()
        self.re = np.exp(spec['values']['X']['logRe'])
        self.scale = spec['hyperfine']['radial_scale']['X']
        self.fixed = spec['hyperfine']['fixed_x_constants']
        self.R = self.solver.R
        self.cache = {}

    def z(self, R):
        return (R**2-self.re**2)/(R**2+self.re**2)/self.scale

    def level(self, J):
        if J not in self.cache:
            E, psi = self.solver.wavefunctions(J)
            self.cache[J] = (E, (self.solver.W[:, None]*psi**2).T)        # (nlev,), (nlev, npts)
        return self.cache[J]

    def average(self, f, v, J):
        return float(self.level(J)[1][v]@f)

    def energy_x00(self, v, J):
        """cm-1 above X(0, 0), the energy scale of BKT02's C_X."""
        return float(self.level(J)[0][v]-self.level(0)[0][0])

    def bkt02(self, v, J):
        q, c = self.fixed['eqQ'], self.fixed['C']
        u, x = v+.5, J*(J+1)
        e = self.energy_x00(v, 0)
        return q[0]+q[1]*u+q[2]*u*u+q[3]*x+q[4]*u*x, c[0]+c[1]*u+c[2]/(e-c[3])


class Radial:
    """f(R) = P(z) inside the knot (linear in z below INNER), an outer form beyond."""

    def __init__(self, lv, knot, order):
        self.lv, self.knot, self.order = lv, knot, order
        R = lv.R
        z, zk, zi = lv.z(R), lv.z(np.array(knot)), lv.z(np.array(INNER))
        self.inside = R <= knot
        self.basis = np.array([np.where(R < INNER, zi**k+k*zi**max(k-1, 0)*(z-zi), z**k) for k in range(order)])
        self.knot_values = np.array([zk**k for k in range(order)])
        self.knot_slopes = np.array([k*zk**max(k-1, 0) for k in range(order)])
        self.dz = z-zk

    def eqq(self, coef, lam, p):
        pk = coef@self.knot_values
        with np.errstate(over="ignore"):         # far out exp(-inf) = 0, as wanted
            outer = EQQ_LIMIT+(pk-EQQ_LIMIT)*np.exp(-(np.maximum(self.lv.R-self.knot, 0.)/lam)**p)
        return np.where(self.inside, coef@self.basis, outer)

    def c(self, coef, slope):
        pk, dpk = coef@self.knot_values, coef@self.knot_slopes
        return np.where(self.inside, coef@self.basis, pk+(dpk+slope)*self.dz)


def bacis_rows():
    out = []
    with open(training.ROOT/'data/x_hyperfine/bacis1980/data.csv') as f:
        for r in csv.DictReader(f):
            v = int(r['v'])
            if v <= CORE_V['X'] or v >= 89 or not r['eqQ_unc']:
                continue
            q = float(r['eqQ'])+(EQQB_SHIFT_TABLE_IV if r['table'] == 'IV' else 0.)
            out.append(dict(v=v, J=14, eqQ=q, sigma=float(r['eqQ_unc']), table=r['table'], source=r['source']))
    return out


def fit(lv, knot, order):
    rad = Radial(lv, knot, order)
    core = [(v, J) for J in range(0, CORE_J['X']+1, 10) for v in range(CORE_V['X']+1)]
    target_q = np.array([lv.bkt02(v, J)[0] for v, J in core])
    target_c = np.array([lv.bkt02(v, J)[1] for v, J in core])
    rows = bacis_rows()

    def resid(t):
        coef, lam, p = t[:order], np.exp(t[order]), np.exp(t[order+1])
        f = rad.eqq(coef, lam, p)
        a = np.array([lv.average(f, v, J) for v, J in core])-target_q
        d = np.array([(lv.average(f, r['v'], r['J'])-r['eqQ'])/r['sigma'] for r in rows])
        return np.concatenate([a/ANCHOR_SIGMA, d])

    # start: P from the core alone (outer held at the knot value), lam 0.9 A, p 1.3 (Bacis' r2 values)
    A = np.array([[lv.average(np.where(rad.inside, b, bk), v, J) for b, bk in zip(rad.basis, rad.knot_values)]
                  for v, J in core])
    coef0 = np.linalg.lstsq(A, target_q, rcond=None)[0]
    t0 = np.concatenate([coef0, [np.log(.9), np.log(1.3)]])
    n_anchor = len(core)
    res = least_squares(resid, t0, loss='soft_l1', f_scale=2., x_scale='jac')
    r = resid(res.x)
    coef, lam, p = res.x[:order], float(np.exp(res.x[order])), float(np.exp(res.x[order+1]))
    # C: anchor by linear least squares (the outer slope does not reach the core levels), then the slope
    Ac = np.array([[lv.average(np.where(rad.inside, b, bk+dbk*rad.dz), v, J)
                    for b, bk, dbk in zip(rad.basis, rad.knot_values, rad.knot_slopes)] for v, J in core])
    coef_c = np.linalg.lstsq(Ac, target_c, rcond=None)[0]
    v83 = C_V83[0]
    base = lv.average(rad.c(coef_c, 0.), v83, 14)
    unit = lv.average(rad.c(coef_c, 1.), v83, 14)-base
    slope = (C_V83[1]-base)/unit
    anchor_c = Ac@coef_c-target_c
    return dict(rad=rad, coef=coef, lam=lam, p=p, coef_c=coef_c, slope_c=slope, rows=rows,
                anchor_rms=float(np.sqrt(np.mean((r[:n_anchor]*ANCHOR_SIGMA)**2))),
                anchor_max=float(np.max(np.abs(r[:n_anchor]*ANCHOR_SIGMA))),
                anchor_c_rms=float(np.sqrt(np.mean(anchor_c**2))), cost=float(res.cost))


def constants(lv, model, v, J):
    """(eqQ MHz, C kHz) of X (v, J): BKT02 in the core, the radial model above."""
    if v <= CORE_V['X']:
        return lv.bkt02(v, J)
    rad = model['rad']
    return (lv.average(rad.eqq(model['coef'], model['lam'], model['p']), v, J),
            lv.average(rad.c(model['coef_c'], model['slope_c']), v, J))


def model_from(release, lv):
    """The model dict of `fit` from a release's X record (hyperfine.json "X" "model")."""
    return dict(rad=Radial(lv, release["knot"], release["order"]), coef=np.array(release["eqQ"]["coef"]),
                lam=release["eqQ"]["lam"], p=release["eqQ"]["p"], coef_c=np.array(release["C"]["coef"]),
                slope_c=release["C"]["slope"])
