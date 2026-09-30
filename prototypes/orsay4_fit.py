"""B levels v' = 51-79 from the Orsay Partie IV atlas and the comb lines there, on a B potential.

Model for an atlas line (v', v'' = 0-1, J''), in MHz:
    obs + offset - model = dB(v', J')        [offset: the atlas scale; dB: correction to the B term value]
and for a comb-referenced line (hyperfine-free, frozen hyperfine offsets as in the potential fit):
    obs - model = dB(v', J').
X v'' = 0 and 1 are the published levels, uncorrected and good to MHz. dB(v', J') = sum_k c_k y^k with
y = J'(J'+1)/1e4, of the degree (0-3) that predicts the level's own lines best when each is left out;
near dissociation the levels bend with J far more than the X levels of part I. The comb lines at
v' = 52, 53, 62 and the atlas lines of the same levels fix the offset. Robust (Huber) weights; atlas line
sigma hypot(eps, FLOOR), comb lines hypot(sigma, 1 MHz). Levels bound by less than 0.3 cm-1 are left
out (quasi-bound above the asymptote, or too close to it for the bound-state solver).

Held-out figure per level: the robust spread (1.4826 MAD) of its lines left out one at a time, with the
atlas's per-line scatter (LINE_SCATTER, the spread of the best-fitting levels) removed in quadrature,
floored at FLOOR_HELD_OUT.

usage: orsay4_fit.py --potential=prototypes/out/b_dissociation_20.json [--floor=1.5] [--maxdeg=3]
                     [--write=level_corrections_2026j --base=level_corrections_2026i]
"""
import csv
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from mlr_b_dissociation import ASYM, GRID, MU, potential as mlr        # noqa: E402

from i2spec.bspline import BSplineSolver
from i2spec.constants import MHZ_PER_CM
from i2spec.model import RovibronicModel
from i2spec.observations import Predictor, load_all

ROOT = Path(__file__).resolve().parents[1]
V_MIN = 51
FLOOR_HELD_OUT = 20.0          # MHz
LINE_SCATTER = 25.0            # MHz: atlas residual spread of v' = 54-56 after the fit, where one line's
                               # scatter is all that is left


def huber(X, r, s, prior=1e4, iters=10, mask=None, cov=False):
    mask = np.ones(len(r), bool) if mask is None else mask
    w = 1 / s
    for _ in range(iters):
        Xw, rw = X[mask] * w[mask, None], r[mask] * w[mask]
        c = np.linalg.solve(Xw.T @ Xw + np.eye(X.shape[1]) / prior**2, Xw.T @ rw)
        z = np.abs(r - X @ c) / s
        w = np.minimum(1.0, 2.5 / np.maximum(z, 1e-9)) / s
    if cov:            # covariance at the final robust weights, scaled by the reduced chi^2
        Xw = X[mask] * w[mask, None]
        C = np.linalg.inv(Xw.T @ Xw + np.eye(X.shape[1]) / prior**2)
        chi2 = float(np.sum(((r[mask] - X[mask] @ c) * w[mask]) ** 2)) / max(int(mask.sum()) - X.shape[1], 1)
        return c, C * max(chi2, 1.0)
    return c


def solve_all(vu, Ju, res0, s, atlas, levels, deg):
    """The joint robust fit: the atlas offset and every level's polynomial."""
    cols = ["offset"] + [(v, k) for v in levels for k in range(deg[v] + 1)]
    ix = {c: i for i, c in enumerate(cols)}
    y = Ju * (Ju + 1) / 1e4
    X = np.zeros((len(vu), len(cols)))
    X[:, 0] = -atlas.astype(float)              # an atlas line: obs - model = dB - offset
    for i in range(len(vu)):
        if vu[i] in deg:
            for k in range(deg[vu[i]] + 1):
                X[i, ix[(vu[i], k)]] = y[i] ** k
    keep = np.isin(vu, levels)
    return huber(X[keep], res0[keep], s[keep]), X, cols, ix


def main(argv):
    opts = {a.split("=")[0]: a.split("=", 1)[1] for a in argv if "=" in a}
    d = json.loads(Path(opts["--potential"]).read_text())
    import mlr_b_dissociation as mb
    mb.N = d["n"]
    B = BSplineSolver(mlr(np.array(d["x"], float)), MU, **GRID, joins=())
    floor = float(opts.get("--floor", 1.5)) * 1e-3 * MHZ_PER_CM
    model, pred = RovibronicModel(), Predictor()
    rows = []                    # (v', J', residual MHz, sigma MHz, is_atlas)

    def base(vu, Ju):
        return B.levels(Ju)[vu]

    with open(ROOT / "data/atlas_lines/orsay1983_part4_assigned.csv", newline="") as f:
        for r in csv.DictReader(f):
            if r["n_assignments"] != "1":
                continue
            vu, vl, J = int(r["v_upper"]), int(r["v_lower"]), int(r["J_lower"])
            Ju = J + (1 if r["branch"] == "R" else -1)
            EX = model.energy("X", vl, J)
            if float(r["sigma_cm1"]) + EX > ASYM - 0.3:
                continue
            eps = float(r["eps_mk"]) * 1e-3 * MHZ_PER_CM if r["eps_mk"] else 150.0
            rows.append((vu, Ju, (float(r["sigma_cm1"]) + EX - base(vu, Ju)) * MHZ_PER_CM, np.hypot(eps, floor), True))
    n_atlas = len(rows)
    sets = load_all()
    # an interval between two lines becomes absolute through its reference component, measured in any set
    # (a compilation's own reference may have given way to its sources: bipm2005a's P(13) 43-0 a3)
    absolute = {}
    for ds in sets:
        scale = MHZ_PER_CM if ds.meta["unit"] == "cm-1" else 1.0
        for o in ds.observations:
            if o.kind == "frequency":
                v = (o.value * scale, o.uncertainty * scale)
                if (o.line, o.component) not in absolute or v[1] < absolute[(o.line, o.component)][1]:
                    absolute[(o.line, o.component)] = v
    for ds in sets:
        scale = MHZ_PER_CM if ds.meta["unit"] == "cm-1" else 1.0
        for o in ds.observations:
            L = o.line
            if L.isotopologue != "127I2" or L.v_upper < V_MIN or L.v_lower > 17:
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
            E = nu / MHZ_PER_CM + model.energy("X", L.v_lower, L.J_lower)
            rows.append((L.v_upper, Ju, (E - base(L.v_upper, Ju)) * MHZ_PER_CM, np.hypot(sigma, 1.0), False))
    print(f"{n_atlas} atlas lines, {len(rows) - n_atlas} comb-referenced rows at v' >= {V_MIN}")

    vu = np.array([r[0] for r in rows]); Ju = np.array([r[1] for r in rows])
    res0 = np.array([r[2] for r in rows]); s = np.array([r[3] for r in rows]); atlas = np.array([r[4] for r in rows])
    levels = [v for v in sorted(set(vu)) if np.sum(vu == v) >= 3]
    maxdeg = int(opts.get("--maxdeg", 3))
    deg = {v: 1 for v in levels}
    for _ in range(2):                  # the offset from a first pass, then each level's degree by leave-one-out
        c, X, cols, ix = solve_all(vu, Ju, res0, s, atlas, levels, deg)
        off = c[0]
        for v in levels:
            m = vu == v
            r_v = res0[m] + np.where(atlas[m], off, 0.0)
            yv = (Ju[m] * (Ju[m] + 1) / 1e4)
            best = None
            for dgr in range(min(maxdeg, max(0, m.sum() // 6)) + 1):
                A = np.vander(yv, dgr + 1, increasing=True)
                loo = [r_v[i] - A[i] @ huber(A, r_v, s[m], iters=4, mask=np.arange(m.sum()) != i) for i in range(m.sum())]
                score = 1.4826 * np.median(np.abs(loo))
                if best is None or score < 0.97 * best[0]:
                    best = (score, dgr)
            deg[v] = best[1]
    c, X, cols, ix = solve_all(vu, Ju, res0, s, atlas, levels, deg)
    keep = np.isin(vu, levels)
    X, res0, s, atlas, vu, Ju = X[keep], res0[keep], s[keep], atlas[keep], vu[keep], Ju[keep]
    c, C = huber(X, res0, s, cov=True)
    fitted = res0 - X @ c
    offset = c[0]
    print(f"atlas offset {offset:+.1f} MHz (comb - atlas on the three shared lines: +82, +150, +138)")
    loo = fitted.copy()
    for i in range(len(res0)):
        m = np.ones(len(res0), bool); m[i] = False
        loo[i] = res0[i] - X[i] @ huber(X, res0, s, iters=4, mask=m)
    mad = lambda x: float(1.4826 * np.median(np.abs(x - np.median(x))))  # noqa: E731
    atlas_scatter = mad(fitted[atlas])
    print(f"atlas residual MAD after the fit {atlas_scatter:.1f} MHz, comb rows rms "
          f"{np.sqrt(np.mean(fitted[~atlas] ** 2)):.2f} MHz")
    coef, coverage, held, covs, disc = {}, {}, {}, {}, {}
    print(f"{'v':>3} {'n':>4} {'deg':>3} {'J′':>9} {'before med':>11} {'after MAD':>10} {'held-out':>9}  dB at J′ = 0 / 40 / 80 (MHz)")
    for v in levels:
        m = vu == v
        cc = [float(c[ix[(v, k)]]) for k in range(deg[v] + 1)]
        coef[v] = cc
        coverage[v] = [int(Ju[m].min()), int(Ju[m].max())]
        spread = 1.4826 * np.median(np.abs(loo[m]))
        held[v] = max(float(np.sqrt(max(spread ** 2 - LINE_SCATTER ** 2, 0.0))), FLOOR_HELD_OUT)
        ids = [ix[(v, k)] for k in range(deg[v] + 1)]
        covs[v] = C[np.ix_(ids, ids)]
        yv = Ju[m] * (Ju[m] + 1) / 1e4
        pv = np.mean([np.vander([y], deg[v] + 1, increasing=True)[0] @ covs[v] @ np.vander([y], deg[v] + 1, increasing=True)[0]
                      for y in yv])
        disc[v] = float(np.sqrt(max(held[v] ** 2 - pv, 0.0)))     # what the held-out figure has beyond the covariance
        at = [np.polyval(cc[::-1], J * (J + 1) / 1e4) for J in (0, 40, 80)]
        print(f"{v:3d} {m.sum():4d} {deg[v]:3d} {coverage[v][0]:4d}-{coverage[v][1]:<4d} {np.median(res0[m]):+11.0f} "
              f"{mad(fitted[m]):10.1f} {held[v]:9.1f}  " + " / ".join(f"{a:+.0f}" for a in at))
    if "--write" in opts:
        name, basename = opts["--write"], opts.get("--base", "level_corrections_2026i")
        out = json.loads((ROOT / f"src/i2spec/data/{basename}.json").read_text())
        out["id"] = name
        for st_key in ("B",):
            out[st_key] = {k: v for k, v in out[st_key].items() if int(k) < V_MIN}
            out["coverage"][st_key] = {k: v for k, v in out["coverage"][st_key].items() if int(k) < V_MIN}
            out["held_out_MHz"][st_key] = {k: v for k, v in out["held_out_MHz"][st_key].items() if int(k) < V_MIN}
        for key in ("covariance_MHz2", "discrepancy_MHz"):
            out.setdefault(key, {"X": {}, "B": {}})
            out[key]["B"] = {k: v for k, v in out[key].get("B", {}).items() if int(k) < V_MIN}
        for v in levels:
            out["B"][str(v)] = coef[v]
            out["coverage"]["B"][str(v)] = coverage[v]
            out["held_out_MHz"]["B"][str(v)] = held[v]
            out["covariance_MHz2"]["B"][str(v)] = covs[v].tolist()
            out["discrepancy_MHz"]["B"][str(v)] = disc[v]
        out["orsay4"] = dict(
            note=(f"B v' = {V_MIN}-79 from the unblended lines of the Orsay atlas Partie IV (Gerstenkorn & Luc 1983) "
                  "and the comb-referenced lines there, on the potential of " + str(opts["--potential"]) +
                  "; prototypes/orsay4_fit.py. Held-out figures: leave-one-line-out rms with the atlas lines' own "
                  f"scatter removed, floored at {FLOOR_HELD_OUT:g} MHz."),
            offset_MHz=float(offset), levels=[int(v) for v in levels])
        (ROOT / f"src/i2spec/data/{name}.json").write_text(json.dumps(out, indent=1) + "\n")
        print(f"wrote src/i2spec/data/{name}.json")


if __name__ == "__main__":
    main(sys.argv[1:])
