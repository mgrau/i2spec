"""Block B of the alternating fit: the B MLR fitted to B levels, with X held at the MLR X.

Usage:  uv run python prototypes/mlr_b_levels.py [n_beta] [Rref/Re]       (about 15 minutes)
Output: prototypes/out/b_levels_<n>_<rf>.json

The joint fit stalls (docs/design/mlr-x.md). Here the problem is cut down to what converged before:
a level fit. Every 127I2 absolute frequency with v'' <= 17 implies a B level once the X level (from
mlr_x_2026b) and the frozen hyperfine offset are removed, at hypot(sigma_obs, 1 MHz); an interval
between two different lines counts too, as an absolute frequency through its reference component
(the BIPM 532 nm lines through R(56) 32-0 a10, R(98) 58-1 through P(13) 43-0 a3); the published B
levels for v' <= 30 guide the shape where data are sparse, at 3 MHz; the atlas bands give v' = 45-62 at
150 MHz. The measured levels are walked in from a 30 MHz weight. Te + De is tied to De(X) plus the
atomic splitting; C5, C6, C8, C10 carry their published uncertainties.
"""

import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

from i2spec.bspline import BSplineSolver
from i2spec.constants import MHZ_PER_CM, reduced_mass
from i2spec.model import RovibronicModel
from i2spec.observations import Predictor, load_all
from i2spec.potentials import MLRPotential, load_potentials

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "prototypes/out"
X0, B0 = load_potentials()["X"], load_potentials()["B"]
MU = reduced_mass("127I2")
GRID = dict(rmin=2.35, rmax=8.0, h=0.01, order=10, nlev=70)
GRID_X = dict(rmin=2.10, rmax=6.0, h=0.01, order=10, nlev=70)
ASYM = X0.De + 7602.9762
PRIORS = {"C5": (B0.C[5], 0.2e5), "C6": (B0.C[6], 0.3e6), "C8": (B0.C[8], 1.0e7), "C10": (B0.C[10], 2.0e8)}
N = int(sys.argv[1]) if len(sys.argv) > 1 else 14
RF = float(sys.argv[2]) if len(sys.argv) > 2 else 1.5
X_PATH = Path(sys.argv[3]) if len(sys.argv) > 3 else ROOT / "data/potentials/mlr_x_2026b.json"
#: Above this v' the B state is perturbed by the 1g state (Chen 2004, v' = 57-60), which a single-channel
#: potential cannot follow; measured levels there carry this floor (MHz) instead of their own precision.
V_PERTURBED, PERTURBED_FLOOR = 50, 300.0


def load_mlr_x():
    d = json.loads(X_PATH.read_text())
    return MLRPotential(De=d["De"], Re=d["Re"], C={int(k): v for k, v in d["C"].items()}, beta=tuple(d["beta"]),
                        p=d["p"], q=d["q"], Rref=d["Rref"], bo=X0)


def potential(x):
    Te, Re, C5, C6, C8, C10 = (float(v) for v in x[N:])
    return MLRPotential(De=ASYM - Te, Re=Re, C={5: C5, 6: C6, 8: C8, 10: C10}, beta=tuple(x[:N]),
                        p=6, q=4, Rref=RF * Re, Te=Te, bo=B0)


def targets():
    """(v', J', E_B, sigma cm-1, kind) from published levels, measured lines and atlas bands."""
    pub = BSplineSolver(B0, MU, **GRID)
    out = [(v, J, pub.energy(v, J), 1e-4, "published") for J in (0, 40, 80, 120, 160) for v in range(31)]
    xm = BSplineSolver(load_mlr_x(), MU, **GRID_X, joins=())
    pred = Predictor()                                  # frozen hyperfine offsets from the default model
    for ds in load_all():
        scale = MHZ_PER_CM if ds.meta["unit"] == "cm-1" else 1.0
        # absolute frequencies of this set, by (line, component), so that an interval between two
        # different lines becomes an absolute frequency through its reference component
        absolute = {(o.line, o.component): (o.value * scale, o.uncertainty * scale)
                    for o in ds.observations if o.kind == "frequency"}
        for o in ds.observations:
            L = o.line
            if L.isotopologue != "127I2" or L.v_lower > 17 or L.v_upper > 62:
                continue
            if o.kind == "frequency":
                value, sigma = o.value * scale, o.uncertainty * scale
            elif o.ref_line is not None and o.ref_line != L and (o.ref_line, o.ref_component) in absolute:
                ref_value, ref_sigma = absolute[(o.ref_line, o.ref_component)]
                value, sigma = ref_value + o.value * scale, np.hypot(o.uncertainty * scale, ref_sigma)
            else:
                continue
            offset = pred.position(L, o.component) - pred.position(L)                 # MHz
            nu = value - offset
            Ju = L.J_lower + (1 if L.branch == "R" else -1)
            E = nu / MHZ_PER_CM + xm.energy(L.v_lower, L.J_lower)
            floor = PERTURBED_FLOOR if L.v_upper > V_PERTURBED else 1.0
            out.append((L.v_upper, Ju, E, np.hypot(sigma, floor) / MHZ_PER_CM, "measured"))
    model = RovibronicModel("127I2", grids={"B": GRID})
    path = OUT / "band_shifts_salami_ross.txt"
    if path.exists():
        for vu, vl, s1, c1, n1, s2, c2, n2 in np.loadtxt(path):
            for shift, contrast, J in ((s1, c1, 25), (s2, c2, 75)):
                if np.isfinite(shift) and np.isfinite(contrast) and contrast >= 1.9 and int(vu) > 44:
                    out.append((int(vu), J + 1, model.energy("B", int(vu), J + 1) + shift, 0.005, "atlas"))
    return out


def main():
    tg = targets()
    kinds = np.array([t[4] for t in tg])
    print(f"{len(tg)} targets: {np.sum(kinds == 'published')} published levels, {np.sum(kinds == 'measured')} measured, "
          f"{np.sum(kinds == 'atlas')} atlas", flush=True)
    Js = sorted({t[1] for t in tg})
    start_file = next((a.split("=", 1)[1] for a in sys.argv if a.startswith("--start=")), OUT / f"b_fit_{N}_{RF}.json")
    start = json.loads(Path(start_file).read_text())["x"]
    x = np.array(start, float)

    def solve(x):
        s = BSplineSolver(potential(x), MU, **GRID, joins=())
        return s, {J: s.wavefunctions(J) for J in Js}

    def sig(weight):
        return np.array([t[3] if t[4] != "measured" else np.hypot(t[3], weight) for t in tg])

    def resid(x, s_arr):
        try:
            _, out = solve(x)
            d = [(out[J][0][v] - E) / sg for (v, J, E, _, _), sg in zip(tg, s_arr)]
            pr = [(x[N + 2 + k] - PRIORS[p][0]) / PRIORS[p][1] for k, p in enumerate(PRIORS)]
            return np.nan_to_num(np.r_[d, pr], nan=1e7, posinf=1e7, neginf=1e7)
        except Exception:
            return np.full(len(tg) + len(PRIORS), 1e7)

    def jac(x, s_arr):
        s, out = solve(x)
        dV = []
        for k in range(len(x)):
            h = 1e-6 * max(abs(x[k]), 1.0)
            xp, xm_ = np.array(x, float), np.array(x, float)
            xp[k] += h
            xm_[k] -= h
            dV.append((potential(xp)(s.R) - potential(xm_)(s.R)) / (2 * h))
        rows = [[(s.W * out[J][1][:, v] ** 2) @ d / sg for d in dV] for (v, J, E, _, _), sg in zip(tg, s_arr)]
        pr = np.zeros((len(PRIORS), len(x)))
        for k, p in enumerate(PRIORS):
            pr[k, N + 2 + k] = 1.0 / PRIORS[p][1]
        return np.vstack([np.array(rows), pr])

    def report(x, label):
        _, out = solve(x)
        d = np.array([out[J][0][v] - E for v, J, E, _, _ in tg]) * MHZ_PER_CM
        msg = " | ".join(f"{k} {np.sqrt(np.mean(d[kinds == k] ** 2)):8.2f} MHz" for k in ("published", "measured", "atlas"))
        m = kinds == "measured"
        vs = np.array([t[0] for t in tg])
        by_v = {"v'≤30": d[m & (vs <= 30)], "31-44": d[m & (vs > 30) & (vs <= 44)], "45-50": d[m & (vs > 44) & (vs <= 50)],
                ">50": d[m & (vs > 50)]}
        print(f"{label:<24} {msg} | measured by v′: " +
              ", ".join(f"{k} {np.sqrt(np.mean(v**2)):.1f}" for k, v in by_v.items() if v.size), flush=True)

    report(x, "start")
    scale = np.r_[np.full(N, 0.1), 1.0, 0.001, 0.2e5, 0.3e6, 1.0e7, 2.0e8]
    lo = np.r_[np.full(N, -np.inf), x[N] - 20, 2.9, 1e5, 5e5, 1e6, 1e7]
    hi = np.r_[np.full(N, np.inf), x[N] + 20, 3.2, 5e5, 3e6, 1e8, 1e9]
    for w in (30.0, 10.0, 3.0, 1.0, 0.3, 0.0):
        s_arr = sig(w / MHZ_PER_CM)
        x = least_squares(resid, x, jac=jac, args=(s_arr,), x_scale=scale, bounds=(lo, hi), max_nfev=40,
                          loss="soft_l1", f_scale=5.0).x
        report(x, f"  measured weight {w:g} MHz")
    tag = next((a.split("=", 1)[1] for a in sys.argv if a.startswith("--tag=")), "")
    json.dump(dict(n=N, rf=RF, x=list(map(float, x)), x_potential=str(X_PATH)),
              open(OUT / f"b_levels_{N}_{RF}{tag}.json", "w"))


if __name__ == "__main__":
    main()
