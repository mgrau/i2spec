"""Fit an MLR B potential to the published levels below v' = 45 and the atlas bands above it.

Usage:  uv run python prototypes/mlr_b_fit.py [n_beta] [Rref/Re]        (about 6 minutes for 14 1.5)
Output: prototypes/out/b_fit_<n>_<rf>.json

The published B potential is out by up to 2.2 cm-1 above v' = 44 (docs/research/spectra-validation.md).
This fits the published levels for v' <= 44 and the levels implied by the Salami & Ross band shifts
above that, walking the atlas weight in from 1 cm-1 to 0.005. Te(B) + De(B) is tied to De(X) plus the
atomic 2P1/2 - 2P3/2 splitting, and C5, C6, C8, C10 carry their published uncertainties.

Needs prototypes/out/band_shifts_salami_ross.txt, which prototypes/band_shifts_salami_ross.py writes.
R_ref > Re matters here: with R_ref = Re the same fit is 188 cm-1 out or numerically unstable
(docs/design/mlr-x.md).
"""
import json, sys
import numpy as np
from scipy.optimize import least_squares, minimize_scalar
from i2spec.bspline import BSplineSolver
from i2spec.constants import MHZ_PER_CM, reduced_mass
from i2spec.model import RovibronicModel
from i2spec.potentials import MLRPotential, load_potentials

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "prototypes/out"
X, B = load_potentials()["X"], load_potentials()["B"]
MU = reduced_mass("127I2")
G = dict(rmin=2.35, rmax=8.0, h=0.01, order=10, nlev=70)
REF = BSplineSolver(B, MU, **G)
ASYM = X.De + 7602.9762
PRIORS = {"C5": (B.C[5], 0.2e5), "C6": (B.C[6], 0.3e6), "C8": (B.C[8], 1.0e7), "C10": (B.C[10], 2.0e8)}
SIGMA_LOW, SIGMA_ATLAS = 0.005, 0.005     # cm-1 both: what the atlas can say about a band position
N = int(sys.argv[1]) if len(sys.argv) > 1 else 14
RF = float(sys.argv[2]) if len(sys.argv) > 2 else 1.5
Re0 = minimize_scalar(lambda r: float(B(np.array([r]))[0]), bracket=(2.9, 3.03, 3.2)).x
Te0 = float(B(np.array([Re0]))[0])


def potential(x):
    Te, Re, C5, C6, C8, C10 = (float(v) for v in x[N:])
    return MLRPotential(De=ASYM - Te, Re=Re, C={5: C5, 6: C6, 8: C8, 10: C10}, beta=tuple(x[:N]),
                        p=6, q=4, Rref=RF * Re, Te=Te, bo=B)


def atlas_targets():
    """E_B(v', J) implied by each atlas band shift, using the published X for the lower level."""
    m = RovibronicModel("127I2", grids={"B": G})
    rows = []
    for vu, vl, s1, c1, n1, s2, c2, n2 in np.loadtxt(OUT / "band_shifts_salami_ross.txt"):
        for shift, contrast, J in ((s1, c1, 25), (s2, c2, 75)):
            if not (np.isfinite(shift) and np.isfinite(contrast)) or contrast < 1.9 or int(vu) <= 44:
                continue
            rows.append((int(vu), J + 1, m.energy("B", int(vu), J + 1) + shift))
    return rows


ATLAS = atlas_targets()
LOW = [(v, J, REF.energy(v, J)) for J in (0, 40, 80, 120, 160) for v in range(45)]
print(f"{len(LOW)} published levels (v' <= 44), {len(ATLAS)} atlas-implied levels "
      f"(v' = {min(a[0] for a in ATLAS)}-{max(a[0] for a in ATLAS)})", flush=True)
JS = sorted({t[1] for t in LOW + ATLAS})


def solve(x):
    s = BSplineSolver(potential(x), MU, **G, joins=())
    return s, {J: s.wavefunctions(J) for J in JS}


def targets(w_atlas):
    return [(v, J, E, SIGMA_LOW) for v, J, E in LOW] + [(v, J, E, w_atlas) for v, J, E in ATLAS]


def resid(x, tgt):
    try:
        _, out = solve(x)
        d = [(out[J][0][v] - E) / s for v, J, E, s in tgt]
        pr = [(x[N + 2 + k] - PRIORS[p][0]) / PRIORS[p][1] for k, p in enumerate(PRIORS)]
        return np.nan_to_num(np.r_[d, pr], nan=1e7, posinf=1e7, neginf=1e7)
    except Exception:
        return np.full(len(tgt) + len(PRIORS), 1e7)


def jac(x, tgt):
    s, out = solve(x)
    dV = []
    for k in range(len(x)):
        h = 1e-6 * max(abs(x[k]), 1.0)
        xp, xm = np.array(x, float), np.array(x, float)
        xp[k] += h; xm[k] -= h
        dV.append((potential(xp)(s.R) - potential(xm)(s.R)) / (2 * h))
    rows = [[(s.W * out[J][1][:, v] ** 2) @ d / sg for d in dV] for v, J, E, sg in tgt]
    pr = np.zeros((len(PRIORS), len(x)))
    for k, p in enumerate(PRIORS):
        pr[k, N + 2 + k] = 1.0 / PRIORS[p][1]
    return np.vstack([np.array(rows), pr])


def report(x, label):
    _, out = solve(x)
    d_low = np.array([out[J][0][v] - E for v, J, E in LOW])
    d_at = np.array([out[J][0][v] - E for v, J, E in ATLAS])
    print(f"{label:<26} v'<=44: rms {np.sqrt(np.mean(d_low**2))*MHZ_PER_CM:9.1f} MHz | "
          f"atlas v'>44: rms {np.sqrt(np.mean(d_at**2)):7.4f} cm-1, max {np.abs(d_at).max():7.4f}", flush=True)


R = np.linspace(2.7, 7.0, 800)
curve = least_squares(lambda b: potential(np.r_[b, Te0, Re0, B.C[5], B.C[6], B.C[8], B.C[10]])(R) - B(R),
                      np.r_[2.0, np.zeros(N - 1)], max_nfev=6000)
x = np.r_[curve.x, Te0, Re0, B.C[5], B.C[6], B.C[8], B.C[10]]
report(x, "curve fit only")
scale = np.r_[np.full(N, 0.1), 1.0, 0.001, 0.2e5, 0.3e6, 1.0e7, 2.0e8]
lo = np.r_[np.full(N, -np.inf), Te0 - 20, 2.9, 1e5, 5e5, 1e6, 1e7]
hi = np.r_[np.full(N, np.inf), Te0 + 20, 3.2, 5e5, 3e6, 1e8, 1e9]
for w in (1.0, 0.3, 0.1, 0.03, SIGMA_ATLAS):
    tgt = targets(w)
    x = least_squares(resid, x, jac=jac, args=(tgt,), x_scale=scale, bounds=(lo, hi), max_nfev=60).x
    report(x, f"  atlas weight {w:g} cm-1")
json.dump(dict(n=N, rf=RF, x=list(map(float, x))), open(OUT / f"b_fit_{N}_{RF}.json", "w"))
