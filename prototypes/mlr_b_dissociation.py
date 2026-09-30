"""The B MLR refitted to the dissociation limit with the Orsay Partie IV atlas (docs/research/orsay-atlas-19700-20035.md).

Usage:  uv run python prototypes/mlr_b_dissociation.py [n_beta] [--start=FILE] [--from=V] [--vmax=79] [--nfev=30] [--tag=]
Output: prototypes/out/b_dissociation_<n><tag>.json

mlr_b_levels.py fitted the B MLR to levels up to v' = 62; above that it had no data and extrapolated to
errors of 10 cm-1. Here the targets are, all as B levels (cm-1 above the X minimum):
  - published B levels v' <= 30, at 3 MHz, as there (the shape inside the well);
  - every 127I2 absolute frequency with v'' <= 17, now at every v', less its frozen hyperfine offset, plus
    the model's X level, at hypot(sigma, 1 MHz), or hypot(sigma, 30 MHz) above v' = 50 where a single
    channel cannot follow the 1g perturbation at the kHz level;
  - the unblended lines of the Orsay Partie IV atlas (v' = 51-79, v'' = 0-1) plus the model's X level, at
    hypot(eps, 1.5 mK), with one offset for the atlas scale; the three atlas lines that were also measured
    against a comb (R(26) 62-0, Goncharov 2007; P(40) 52-0 and P(52) 53-0, Sakamoto 2024) put it at
    +82, +150 and +138 MHz, and the offset carries that as a prior, 123(10) MHz (a looser prior lets the offset trade against the
    potential: it went to -186 and +296 MHz at 20 MHz);
    levels within 0.3 cm-1 of the asymptote or above it (quasi-bound, behind the centrifugal barrier) are
    left out.
The centrifugal term carries one more parameter, q_far: a constant added to the published alpha(R)
beyond 5 A (MLRPotential.q_far). Without it every band shows the same J trend, +500 MHz at J' = 0 to
-500 at J' = 80: alpha(R) was fitted in the well, and its extrapolation (-8e-4 out to 20 A) sets the
rotational energy of the last levels. Levels come from a 40 A box with a graded mesh (0.01 A to 5 A, coarser beyond; exact to 1e-4 MHz for
every level bound by more than 0.3 cm-1). Levels are added by v' in stages up to --vmax, each stage
starting from the last, because the starting potential binds only v' <= 75.
"""

import csv
import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

from i2spec.bspline import BSplineSolver
from i2spec.constants import HBAR2_2U, MHZ_PER_CM, reduced_mass
from i2spec.model import RovibronicModel
from i2spec.observations import Predictor, load_all
from i2spec.potentials import MLRPotential, load_potentials

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "prototypes/out"
X0, B0 = load_potentials()["X"], load_potentials()["B"]
MU = reduced_mass("127I2")
GRID = dict(rmin=2.35, rmax=40.0, h=0.01, order=10, nlev=95, mesh=((5.0, 0.02), (8.0, 0.05), (15.0, 0.1)))
ASYM = X0.De + 7602.9762
PRIORS = {"C5": (B0.C[5], 0.2e5), "C6": (B0.C[6], 0.3e6), "C8": (B0.C[8], 1.0e7), "C10": (B0.C[10], 2.0e8)}
V_SINGLE_CHANNEL, SINGLE_CHANNEL_FLOOR = 50, 30.0       # MHz
ATLAS_FLOOR = 1.5e-3                                    # cm-1
ATLAS_OFFSET = (123.0 / MHZ_PER_CM, 10.0 / MHZ_PER_CM)  # comb - atlas on the three lines both hold
ARGS = {a.split("=", 1)[0]: a.split("=", 1)[1] for a in sys.argv[1:] if a.startswith("--")}
N = int(next((a for a in sys.argv[1:] if not a.startswith("--")), 14))
VMAX = int(ARGS.get("--vmax", 79))


def potential(x):
    Te, Re, C5, C6, C8, C10 = (float(v) for v in x[N:N + 6])
    q_far = float(x[N + 7]) if len(x) > N + 7 else 0.0
    return MLRPotential(De=ASYM - Te, Re=Re, C={5: C5, 6: C6, 8: C8, 10: C10}, beta=tuple(x[:N]),
                        p=6, q=4, Rref=1.5 * Re, Te=Te, bo=B0, q_far=q_far)


def targets():
    """(v', J', E_B, sigma cm-1, kind): kind 'published', 'measured' or 'atlas' (E_B without the atlas offset)."""
    pub = BSplineSolver(B0, MU, **{**GRID, "rmax": 8.0, "mesh": ()})
    out = [(v, J, pub.energy(v, J), 1e-4, "published") for J in (0, 40, 80, 120, 160) for v in range(31)]
    model, pred = RovibronicModel(), Predictor()
    for ds in load_all():
        scale = MHZ_PER_CM if ds.meta["unit"] == "cm-1" else 1.0
        absolute = {(o.line, o.component): (o.value * scale, o.uncertainty * scale)
                    for o in ds.observations if o.kind == "frequency"}
        for o in ds.observations:
            L = o.line
            if L.isotopologue != "127I2" or L.v_lower > 17:
                continue
            if o.kind == "frequency":
                value, sigma = o.value * scale, o.uncertainty * scale
            elif o.ref_line is not None and o.ref_line != L and (o.ref_line, o.ref_component) in absolute:
                ref_value, ref_sigma = absolute[(o.ref_line, o.ref_component)]
                value, sigma = ref_value + o.value * scale, np.hypot(o.uncertainty * scale, ref_sigma)
            else:
                continue
            nu = value - (pred.position(L, o.component) - pred.position(L))
            Ju = L.J_lower + (1 if L.branch == "R" else -1)
            floor = SINGLE_CHANNEL_FLOOR if L.v_upper > V_SINGLE_CHANNEL else 1.0
            out.append((L.v_upper, Ju, nu / MHZ_PER_CM + model.energy("X", L.v_lower, L.J_lower),
                        np.hypot(sigma, floor) / MHZ_PER_CM, "measured"))
    with open(ROOT / "data/atlas_lines/orsay1983_part4_assigned.csv", newline="") as f:
        for r in csv.DictReader(f):
            if r["n_assignments"] != "1":
                continue
            vl, J = int(r["v_lower"]), int(r["J_lower"])
            E = float(r["sigma_cm1"]) + model.energy("X", vl, J)
            if E > ASYM - 0.3:
                continue
            eps = float(r["eps_mk"]) * 1e-3 if r["eps_mk"] else 5e-3
            out.append((int(r["v_upper"]), J + (1 if r["branch"] == "R" else -1), E, np.hypot(eps, ATLAS_FLOOR), "atlas"))
    return out


def main():
    everything = targets()
    kinds_all = np.array([t[4] for t in everything])
    print(f"{len(everything)} targets: " + ", ".join(f"{np.sum(kinds_all == k)} {k}" for k in ("published", "measured", "atlas")),
          flush=True)
    start = ARGS.get("--start")
    if start:
        d = json.loads(Path(start).read_text())
        n0, x0 = d["n"], np.array(d["x"], float)
        x = np.r_[x0[:n0], np.zeros(max(0, N - n0)), x0[n0:]] if N >= n0 else None     # extra beta start at 0
    else:
        d = json.loads((ROOT / "src/i2spec/data/mlr_b_2026c.json").read_text())
        x = np.r_[d["beta"], d["Te"], d["Re"], d["C"]["5"], d["C"]["6"], d["C"]["8"], d["C"]["10"], 0.0]
    n_par = N + 8                                           # beta, Te, Re, C5, C6, C8, C10, atlas offset, q_far

    def solve(x, Js):
        s = BSplineSolver(potential(x), MU, **GRID, joins=())
        return s, {J: s.wavefunctions(J) for J in Js}

    def fit(tg, x, stage):
        kinds = np.array([t[4] for t in tg])
        Js = sorted({t[1] for t in tg})
        sg = np.array([t[3] for t in tg])
        atlas = kinds == "atlas"

        def diffs(x, out):
            return np.array([out[J][0][v] - E for v, J, E, _, _ in tg]) - np.where(atlas, x[N + 6], 0.0)

        def resid(x):
            try:
                _, out = solve(x, Js)
                pr = [(x[N + 2 + k] - PRIORS[p][0]) / PRIORS[p][1] for k, p in enumerate(PRIORS)]
                pr.append((x[N + 6] - ATLAS_OFFSET[0]) / ATLAS_OFFSET[1])
                return np.nan_to_num(np.r_[diffs(x, out) / sg, pr], nan=1e7, posinf=1e7, neginf=1e7)
            except Exception:
                return np.full(len(tg) + len(PRIORS) + 1, 1e7)

        def jac(x):
            s, out = solve(x, Js)
            dV = []
            for k in range(N + 6):
                h = 1e-6 * max(abs(x[k]), 1.0)
                xp, xm = np.array(x, float), np.array(x, float)
                xp[k] += h
                xm[k] -= h
                dV.append((potential(xp)(s.R) - potential(xm)(s.R)) / (2 * h))
            dV = np.array(dV).T
            rows = np.array([(s.W * out[J][1][:, v] ** 2) @ dV for v, J, _, _, _ in tg])
            drot = HBAR2_2U / MU * potential(x).far_switch(s.R) / s.R ** 2          # dE/dq_far per J(J+1)
            dq = np.array([J * (J + 1) * (s.W * out[J][1][:, v] ** 2) @ drot for v, J, _, _, _ in tg])
            rows = np.c_[rows, -atlas.astype(float), dq] / sg[:, None]
            pr = np.zeros((len(PRIORS) + 1, n_par))
            for k, p in enumerate(PRIORS):
                pr[k, N + 2 + k] = 1.0 / PRIORS[p][1]
            pr[-1, N + 6] = 1.0 / ATLAS_OFFSET[1]
            return np.vstack([rows, pr])

        def report(x, label):
            _, out = solve(x, Js)
            d = diffs(x, out) * MHZ_PER_CM
            vs = np.array([t[0] for t in tg])
            parts = [f"{k} {np.sqrt(np.mean(d[kinds == k] ** 2)):.1f}" for k in ("published", "measured", "atlas")]
            bands = [f"{lo}-{hi} {np.median(np.abs(d[atlas & (vs >= lo) & (vs <= hi)])):.0f}"
                     for lo, hi in ((51, 56), (57, 61), (62, 66), (67, 71), (72, 76), (77, 79)) if np.any(atlas & (vs >= lo) & (vs <= hi))]
            print(f"{label:<28} rms MHz: {' | '.join(parts)} | atlas median |d| by v′: {', '.join(bands)} | "
                  f"offset {x[N + 6] * MHZ_PER_CM:+.1f} MHz | q_far {x[N + 7]:+.2e}", flush=True)

        report(x, f"stage v′≤{stage} start")
        scale = np.r_[np.full(N, 0.1), 1.0, 0.001, 0.2e5, 0.3e6, 1.0e7, 2.0e8, 0.003, 1e-3]
        lo = np.r_[np.full(N, -np.inf), x[N] - 20, 2.9, 1e5, 5e5, 1e6, 1e7, -0.02, -0.05]
        hi = np.r_[np.full(N, np.inf), x[N] + 20, 3.2, 5e5, 3e6, 1e8, 1e9, 0.02, 0.05]
        for loss in ("linear", "soft_l1"):       # linear first: the start can be tens of GHz off
            x = least_squares(resid, x, jac=jac, x_scale=scale, bounds=(lo, hi), max_nfev=int(ARGS.get("--nfev", 30)),
                              loss=loss, f_scale=3.0).x
            report(x, f"stage v′≤{stage} {loss}")
        return x

    x = np.r_[x, np.zeros(N + 8 - len(x))]        # a start without the atlas offset or q_far
    stages = [s for s in (66, 71, 75, 77, VMAX) if s <= VMAX and s >= int(ARGS.get("--from", 0))]
    for stage in sorted(set(stages)):
        x = fit([t for t in everything if t[4] != "atlas" or t[0] <= stage], x, stage)
        json.dump(dict(n=N, x=list(map(float, x)), stage=stage, grid={k: v for k, v in GRID.items()}),
                  open(OUT / f"b_dissociation_{N}{ARGS.get('--tag', '')}.json", "w"))


if __name__ == "__main__":
    main()
