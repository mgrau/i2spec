"""Refit the MLR X potential to the Orsay-atlas levels v'' = 18-25, so their information reaches v'' = 26-47.

Usage:  .venv/bin/python prototypes/mlr_x_atlas.py [--nbeta=12] [--write=mlr_x_2026d]
Output: prototypes/out/mlr_x_atlas.json (and data/potentials/<name>.json with --write)

i2spec2026g corrects X v'' = 18-25 level by level (level_corrections_2026f, from 3 700 assigned lines of
the Orsay atlas part I). Corrections stop at the levels they were fitted on; the potential does not.
This refits mlr_x_2026c - the extended-range X that supplies v'' >= 18 - to reproduce:

  * v'' = 18-25: the corrected levels, E_2026c(v, J) + dX(v, J), at J inside each level's atlas coverage,
    with each level's uncertainty (15-52 MHz);
  * v'' <= 17: mlr_x_2026c's own levels at 3 MHz, as the anchor it already was (inside the published
    range the hybrid takes the published curves and their comb-level corrections, so these only keep
    the MLR continuous with them);
  * v'' = 48, 53, 54: the corrected emission levels at their own uncertainty (7-20 MHz);
  * De, C6, C8, C10: their published priors.

The Jacobian is Hellmann-Feynman: dE/dtheta = <psi| dV/dtheta |psi>, dV by central differences.
Leave-one-level-out over v'' = 18-25 says whether the potential can predict a level it has not seen.
"""
import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

from i2spec.bspline import BSplineSolver
from i2spec.constants import MHZ_PER_CM, reduced_mass
from i2spec.level_corrections import load_level_corrections
from i2spec.potentials import MLRPotential, load_potentials

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "prototypes/out"
X0 = load_potentials()["X"]
MU = reduced_mass("127I2")
GRID = dict(rmin=2.10, rmax=6.5, h=0.01, order=10, nlev=60)
START = json.loads((ROOT / "src/i2spec/data/mlr_x_2026c.json").read_text())
CORR = load_level_corrections("level_corrections_2026f")
PRIORS = {"De": (X0.De, 0.128), "C6": (X0.C[6], 0.12e6), "C8": (X0.C[8], 1.20e7), "C10": (X0.C[10], 0.5e8)}
ATLAS = range(18, 26)
HIGH = (48, 53, 54)
ANCHOR_SIGMA = 3.0      # MHz, v'' <= 17


def potential(x, n):
    De, Re, C6, C8, C10 = (float(v) for v in x[n:])
    return MLRPotential(De=De, Re=Re, C={6: C6, 8: C8, 10: C10}, beta=tuple(x[:n]), p=START["p"], q=START["q"],
                        Rref=START["Rref"], bo=X0)


def x_start(n):
    b = list(START["beta"]) + [0.0] * max(0, n - len(START["beta"]))
    return np.r_[b[:n], START["De"], START["Re"], START["C"]["6"], START["C"]["8"], START["C"]["10"]]


def targets(n):
    base = BSplineSolver(potential(x_start(n), n), MU, **GRID, joins=())
    tg = []
    for v in range(0, 18):
        for J in (0, 60, 120, 180):
            tg.append(dict(v=v, J=J, E=base.energy(v, J), s=ANCHOR_SIGMA / MHZ_PER_CM, kind="anchor"))
    for v in list(ATLAS) + list(HIGH):
        lo, hi = CORR.coverage["X"][v]
        Js = sorted({int(J) for J in np.linspace(lo, hi, 5 if hi - lo > 40 else 2)})
        u = CORR.uncertainty("X", v) or 20.0
        for J in Js:
            tg.append(dict(v=v, J=J, E=base.energy(v, J) + CORR.shift("X", v, J) / MHZ_PER_CM,
                           s=u / MHZ_PER_CM, kind="atlas" if v in ATLAS else "emission"))
    return tg


class Fit:
    def __init__(self, n, tg):
        self.n, self.tg = n, tg
        self.Js = sorted({t["J"] for t in tg})
        self.names = ("De", "C6", "C8", "C10")
        self.pk = {"De": n, "C6": n + 2, "C8": n + 3, "C10": n + 4}

    def solve(self, x):
        s = BSplineSolver(potential(x, self.n), MU, **GRID, joins=())
        return s, {J: s.wavefunctions(J) for J in self.Js}

    def resid(self, x, mask):
        try:
            _, out = self.solve(x)
            d = [(out[t["J"]][0][t["v"]] - t["E"]) / t["s"] for t, m in zip(self.tg, mask) if m]
            pr = [(x[self.pk[p]] - PRIORS[p][0]) / PRIORS[p][1] for p in self.names]
            return np.nan_to_num(np.r_[d, pr], nan=1e7, posinf=1e7, neginf=1e7)
        except Exception:
            return np.full(int(np.sum(mask)) + 4, 1e7)

    def jac(self, x, mask):
        s, out = self.solve(x)
        dV = []
        for k in range(len(x)):
            h = 1e-6 * max(abs(x[k]), 1.0)
            xp, xm = np.array(x, float), np.array(x, float)
            xp[k] += h
            xm[k] -= h
            dV.append((potential(xp, self.n)(s.R) - potential(xm, self.n)(s.R)) / (2 * h))
        dV = np.array(dV)
        rows = []
        for t, m in zip(self.tg, mask):
            if m:
                w = s.W * out[t["J"]][1][:, t["v"]] ** 2
                rows.append((dV @ w) / t["s"])
        pr = np.zeros((4, len(x)))
        for r, p in enumerate(self.names):
            pr[r, self.pk[p]] = 1.0 / PRIORS[p][1]
        return np.vstack([np.array(rows), pr])

    def run(self, x, mask, nfev=30):
        n = self.n
        scale = np.r_[np.full(n, 0.05), 0.128, 0.0005, 0.12e6, 1.2e7, 0.5e8]
        lo = np.r_[np.full(n, -np.inf), X0.De - 5, 2.6, 0.5e6, 1e6, 1e6]
        hi = np.r_[np.full(n, np.inf), X0.De + 5, 2.75, 3.0e6, 1.5e8, 4.0e8]
        return least_squares(self.resid, x, jac=self.jac, args=(mask,), x_scale=scale, bounds=(lo, hi),
                             max_nfev=nfev, loss="soft_l1", f_scale=3.0).x

    def levels(self, x, vs, J):
        s = BSplineSolver(potential(x, self.n), MU, **GRID, joins=())
        return np.array([s.energy(v, J) for v in vs])


def main(argv):
    opts = {a.split("=")[0]: a.split("=", 1)[1] for a in argv if "=" in a}
    n = int(opts.get("--nbeta", len(START["beta"])))
    tg = targets(n)
    f = Fit(n, tg)
    kind = np.array([t["kind"] for t in tg])
    v_of = np.array([t["v"] for t in tg])
    x0 = x_start(n)
    all_ = np.ones(len(tg), bool)

    def report(x, label):
        _, out = f.solve(x)
        d = np.array([(out[t["J"]][0][t["v"]] - t["E"]) * MHZ_PER_CM for t in tg])
        z = np.array([(out[t["J"]][0][t["v"]] - t["E"]) / t["s"] for t in tg])
        rms = lambda a: float(np.sqrt(np.mean(a ** 2)))      # noqa: E731
        print(f"{label:<28} anchor {rms(d[kind == 'anchor']):6.2f} MHz | atlas {rms(d[kind == 'atlas']):6.1f} MHz "
              f"(chi/pt {np.mean(z[kind == 'atlas']**2):5.2f}) | emission {rms(d[kind == 'emission']):6.1f} MHz", flush=True)
        return d

    report(x0, "start (mlr_x_2026c)")
    x = f.run(x0, all_)
    d = report(x, "refit")
    print("per level (MHz, mean of targets):", {int(v): round(float(np.mean(d[v_of == v])), 1) for v in list(ATLAS) + list(HIGH)})

    # leave one atlas level out: fit without it, predict it
    print("\nleave-one-level-out over the atlas levels:")
    loo = {}
    for v in ATLAS:
        m = all_ & (v_of != v)
        xv = f.run(x, m, nfev=15)
        _, out = f.solve(xv)
        pr = np.array([(out[t["J"]][0][t["v"]] - t["E"]) * MHZ_PER_CM for t in tg if t["v"] == v])
        u = CORR.uncertainty("X", v)
        start_err = np.array([CORR.shift("X", v, t["J"]) for t in tg if t["v"] == v])
        loo[v] = float(np.sqrt(np.mean(pr ** 2)))
        print(f"   v''={v}: held-out {loo[v]:6.1f} MHz  (level uncertainty {u:5.1f}; mlr_x_2026c was off by "
              f"{np.sqrt(np.mean(start_err**2)):5.1f})", flush=True)

    # what moves where there are no data
    vs = list(range(18, 48))
    for J in (50, 150):
        e0, e1 = f.levels(x0, vs, J), f.levels(x, vs, J)
        print(f"\nJ={J}: refit - mlr_x_2026c (MHz) by v'':", " ".join(f"{v}:{(b - a) * MHZ_PER_CM:+.0f}" for v, a, b in zip(vs, e0, e1)))

    De, Re, C6, C8, C10 = (float(v) for v in x[n:])
    res = dict(START, De=De, Re=Re, C={"6": C6, "8": C8, "10": C10}, beta=list(map(float, x[:n])),
               id="mlr_x_atlas", loo_MHz={str(k): v for k, v in loo.items()},
               fitted_to=("mlr_x_2026c refitted to the Orsay-atlas levels v'' = 18-25 (i2spec2026g), anchored on its own "
                          "levels v'' <= 17 at 3 MHz and the emission levels v'' = 48, 53, 54; prototypes/mlr_x_atlas.py"))
    OUT.mkdir(exist_ok=True)
    (OUT / "mlr_x_atlas.json").write_text(json.dumps(res, indent=1))
    if "--write" in opts:
        name = opts["--write"]
        res.update(id=name, citation=f"i2spec, {name}: {res['fitted_to']}. Not a published parameter set.",
                   status="research result (docs/research/orsay-atlas-11000-14000.md)")
        res.pop("loo_MHz")
        (ROOT / "src/i2spec/data" / f"{name}.json").write_text(json.dumps(res, indent=1))
        print("wrote", name)


if __name__ == "__main__":
    main(sys.argv[1:])
