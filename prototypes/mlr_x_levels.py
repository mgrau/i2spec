"""Block X of the alternating fit: the X MLR fitted to X levels, with B held at the level-fit B.

Usage:  uv run python prototypes/mlr_x_levels.py [b_json]       (about 15 minutes)
Output: prototypes/out/x_levels.json   (the schema of data/potentials/mlr_x_*.json)

The mirror of prototypes/mlr_b_levels.py. Every 127I2 absolute frequency -- or interval between two
lines, through its reference component -- with v' <= 50 implies an X level once the B level (from the
level-fit B) and the frozen hyperfine offset are removed, at hypot(sigma_obs, 1 MHz for v' <= 30,
3 MHz above, which is what the B is good to there). The published X levels for v'' <= 17 guide the
shape at 3 MHz where the sets are sparse. De, C6, C8, C10 carry their published priors.
"""

import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

from i2spec.bspline import BSplineSolver
from i2spec.constants import MHZ_PER_CM, reduced_mass
from i2spec.observations import Predictor, load_all
from i2spec.potentials import MLRPotential, load_potentials

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "prototypes/out"
X0, B0 = load_potentials()["X"], load_potentials()["B"]
MU = reduced_mass("127I2")
GRID_X = dict(rmin=2.10, rmax=6.0, h=0.01, order=10, nlev=70)
GRID_B = dict(rmin=2.35, rmax=8.0, h=0.01, order=10, nlev=70)
ASYM = X0.De + 7602.9762
PRIORS = {"De": (X0.De, 0.128), "C6": (X0.C[6], 0.12e6), "C8": (X0.C[8], 1.20e7), "C10": (X0.C[10], 0.5e8)}
START = json.loads((ROOT / "data/potentials/mlr_x_2026b.json").read_text())
N = len(START["beta"])


def load_b(path):
    d = json.loads(Path(path).read_text())
    n, rf, x = d["n"], d["rf"], d["x"]
    Te, Re, C5, C6, C8, C10 = x[n:]
    return MLRPotential(De=ASYM - Te, Re=Re, C={5: C5, 6: C6, 8: C8, 10: C10}, beta=tuple(x[:n]), p=6, q=4,
                        Rref=rf * Re, Te=Te, bo=B0)


def potential(x):
    De, Re, C6, C8, C10 = (float(v) for v in x[N:])
    return MLRPotential(De=De, Re=Re, C={6: C6, 8: C8, 10: C10}, beta=tuple(x[:N]), p=START["p"], q=START["q"],
                        Rref=START["Rref"], bo=X0)


def targets(b_path):
    pub = BSplineSolver(X0, MU, **GRID_X)
    e0 = pub.energy(0, 0)
    out = [(v, J, pub.energy(v, J) - e0, 1e-4, "published") for J in (0, 60, 120) for v in range(18)]
    bm = BSplineSolver(load_b(b_path), MU, **GRID_B, joins=())
    pred = Predictor()
    for ds in load_all():
        scale = MHZ_PER_CM if ds.meta["unit"] == "cm-1" else 1.0
        absolute = {(o.line, o.component): (o.value * scale, o.uncertainty * scale)
                    for o in ds.observations if o.kind == "frequency"}
        for o in ds.observations:
            L = o.line
            if L.isotopologue != "127I2" or L.v_upper > 50 or L.v_lower > 60:
                continue
            if o.kind == "frequency":
                value, sigma = o.value * scale, o.uncertainty * scale
            elif o.ref_line is not None and o.ref_line != L and (o.ref_line, o.ref_component) in absolute:
                rv, rs = absolute[(o.ref_line, o.ref_component)]
                value, sigma = rv + o.value * scale, np.hypot(o.uncertainty * scale, rs)
            else:
                continue
            offset = pred.position(L, o.component) - pred.position(L)
            Ju = L.J_lower + (1 if L.branch == "R" else -1)
            E = bm.energy(L.v_upper, Ju) - (value - offset) / MHZ_PER_CM          # X level, potential-minimum origin
            model_error = 1.0 if L.v_upper <= 30 else 3.0
            out.append((L.v_lower, L.J_lower, E, np.hypot(sigma, model_error) / MHZ_PER_CM, "measured"))
    return out, e0


def main():
    b_path = sys.argv[1] if len(sys.argv) > 1 else OUT / "b_levels_14_1.5.json"
    tg, e0_pub = targets(b_path)
    kinds = np.array([t[4] for t in tg])
    print(f"{len(tg)} targets: {np.sum(kinds == 'published')} published levels, {np.sum(kinds == 'measured')} measured", flush=True)
    Js = sorted({t[1] for t in tg})
    x = np.r_[START["beta"], START["De"], START["Re"], START["C"]["6"], START["C"]["8"], START["C"]["10"]]

    def solve(x):
        s = BSplineSolver(potential(x), MU, **GRID_X, joins=())
        return s, {J: s.wavefunctions(J) for J in Js}

    def energies(out):
        """Published-level targets are relative to X(0,0); measured ones are absolute above the minimum."""
        return out

    def sig(weight):
        return np.array([t[3] if t[4] != "measured" else np.hypot(t[3], weight) for t in tg])

    def value(out, t):
        v, J, E, _, kind = t
        return out[J][0][v] - (out[0][0][0] if kind == "published" else 0.0)

    def resid(x, s_arr):
        try:
            _, out = solve(x)
            d = [(value(out, t) - t[2]) / sg for t, sg in zip(tg, s_arr)]
            pr = [(x[N + k] - PRIORS[p][0]) / PRIORS[p][1] for k, p in ((0, "De"), (2, "C6"), (3, "C8"), (4, "C10"))]
            return np.nan_to_num(np.r_[d, pr], nan=1e7, posinf=1e7, neginf=1e7)
        except Exception:
            return np.full(len(tg) + 4, 1e7)

    def jac(x, s_arr):
        s, out = solve(x)
        dV = []
        for k in range(len(x)):
            h = 1e-6 * max(abs(x[k]), 1.0)
            xp, xm = np.array(x, float), np.array(x, float)
            xp[k] += h
            xm[k] -= h
            dV.append((potential(xp)(s.R) - potential(xm)(s.R)) / (2 * h))
        w0 = s.W * out[0][1][:, 0] ** 2
        rows = []
        for t, sg in zip(tg, s_arr):
            v, J, _, _, kind = t
            w = s.W * out[J][1][:, v] ** 2
            rows.append([((w @ d) - (w0 @ d if kind == "published" else 0.0)) / sg for d in dV])
        pr = np.zeros((4, len(x)))
        for row, (k, p) in enumerate(((0, "De"), (2, "C6"), (3, "C8"), (4, "C10"))):
            pr[row, N + k] = 1.0 / PRIORS[p][1]
        return np.vstack([np.array(rows), pr])

    def report(x, label):
        _, out = solve(x)
        d = np.array([value(out, t) - t[2] for t in tg]) * MHZ_PER_CM
        m = kinds == "measured"
        vl = np.array([t[0] for t in tg])
        by = {"v″≤5": d[m & (vl <= 5)], "6-17": d[m & (vl > 5) & (vl <= 17)], "≥48": d[m & (vl >= 48)]}
        print(f"{label:<24} published {np.sqrt(np.mean(d[~m]**2)):7.2f} MHz | measured {np.sqrt(np.mean(d[m]**2)):8.2f} MHz | by v″: "
              + ", ".join(f"{k} {np.sqrt(np.mean(v**2)):.2f}" for k, v in by.items() if v.size), flush=True)

    report(x, "start (mlr_x_2026b)")
    scale = np.r_[np.full(N, 0.1), 0.128, 0.001, 0.12e6, 1.20e7, 0.5e8]
    lo = np.r_[np.full(N, -np.inf), X0.De - 5, 2.5, 0.5e6, 1e6, 1e6]
    hi = np.r_[np.full(N, np.inf), X0.De + 5, 2.8, 3.0e6, 1.5e8, 4.0e8]
    for w in (10.0, 3.0, 1.0, 0.0):
        s_arr = sig(w / MHZ_PER_CM)
        x = least_squares(resid, x, jac=jac, args=(s_arr,), x_scale=scale, bounds=(lo, hi), max_nfev=40,
                          loss="soft_l1", f_scale=5.0).x
        report(x, f"  measured weight {w:g} MHz")
    De, Re, C6, C8, C10 = (float(v) for v in x[N:])
    json.dump(dict(START, De=De, Re=Re, C={"6": C6, "8": C8, "10": C10}, beta=list(map(float, x[:N])),
                   id="x_levels", fitted_to=f"X levels implied by the observations with B from {b_path}",
                   status="alternating fit, block X"),
              open(OUT / "x_levels.json", "w"), indent=1)


if __name__ == "__main__":
    main()
