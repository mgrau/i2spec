"""Level corrections for the near infrared: δX(v'', J'') and δB(v', J') on top of the potentials.

The potentials miss the comb-referenced NIR lines (755-815 nm, bands (0-7)-(11-17)) by 2-6 MHz where
the data are 0.015-0.3 MHz. IodineSpec5 answers that with a local Dunham model fitted to those data.
Here the same information goes into corrections to the term values the potentials give,
    δE(state, v, J) = Σ_k c_k(state, v) y^k,   y = J(J+1)/1e4,
which is linear in the c_k, so the fit is a weighted ridge regression on the residuals, and every band
sharing a level gets the correction. δB(v' = 0) is held at zero: with nearly every line from v' = 0,
X and B(0) corrections are degenerate, and the convention puts the sum on X.

usage: nir_corrections.py [--floor=0.1] [--prior=5] [--write]
Prints in-sample and leave-one-set-out rms per set, and with --write saves the corrections to
src/i2spec/data/level_corrections_2026a.json.
"""
from __future__ import annotations

import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

from i2spec.constants import MHZ_PER_CM
from i2spec.observations import Predictor, load_all

ROOT = Path(__file__).resolve().parents[1]
SETS = ("liao2010a", "bodermann2000a", "bodermann1998b", "reinhardt2007a", "cornish2000a")
BASE = "i2spec2026a"          # the set the corrections sit on; never the default, which may already carry them
X_V = range(11, 18)          # v'' with corrections
B_V = range(1, 8)            # v' > 0 with corrections (constants only); v' = 0 is the reference
DEGREE_X = {11: 0, 12: 2, 13: 2, 14: 3, 15: 1, 16: 1, 17: 0}   # polynomial degree in y per v''
DEGREE_B = {1: 1, 2: 0, 3: 0, 5: 0, 7: 0}                          # v' = 1 needs a J term: band 1-14 spans 12 MHz


#: Dunham-like form: δX = Σ a_lk u^l y^k with u = v'' − 14, y = J''(J''+1)/1e4 (smooth in v'' and J'');
#: δB(v' ≥ 1) = Σ b_l (v')^l, l = 1, 2. FORM = "perv" is the per-v'' polynomial instead.
X_TERMS = [(0, 0), (1, 0), (2, 0), (3, 0), (0, 1), (1, 1), (2, 1), (0, 2), (1, 2)]
B_TERMS = [1, 2]
FORM = "dunham"


def basis(line):
    """Column names and values of the design matrix for one line: δν = δB(v', J') − δX(v'', J'')."""
    cols = {}
    y = line.J_lower * (line.J_lower + 1) / 1e4
    if line.v_lower in X_V:
        if FORM == "perv":
            for k in range(DEGREE_X[line.v_lower] + 1):
                cols[("X", line.v_lower, k)] = -(y ** k)
        else:
            u = line.v_lower - 14
            for l, k in X_TERMS:
                cols[("X", l, k)] = -(u ** l) * (y ** k)
    if line.v_upper in B_V:
        if FORM == "perv":
            Ju = line.J_lower + (1 if line.branch == "R" else -1)
            for k in range(DEGREE_B.get(line.v_upper, 0) + 1):
                cols[("B", line.v_upper, k)] = (Ju * (Ju + 1) / 1e4) ** k
        else:
            for l in B_TERMS:
                cols[("B", l, 0)] = float(line.v_upper ** l)
    return cols


def main(argv):
    global FORM, X_TERMS
    opts = {a.split("=")[0]: a.split("=", 1)[1] if "=" in a else True for a in argv}
    FORM = opts.get("--form", FORM)
    if "--xterms" in opts:                        # e.g. --xterms=00,10,20,01,11,02
        X_TERMS = [(int(t[0]), int(t[1])) for t in opts["--xterms"].split(",")]
    floor = float(opts.get("--floor", 0.1))
    prior = float(opts.get("--prior", 5.0))       # MHz, ridge prior on every coefficient
    pred = Predictor(BASE)                     # BASE carries no corrections; see level_corrections_fit.BarePredictor
    rows = []                                      # (set, line, residual MHz, sigma MHz, cols)
    for ds in load_all():
        if ds.id not in SETS:
            continue
        scale = MHZ_PER_CM if ds.unit == "cm-1" else 1.0
        for o in ds.observations:
            if o.kind == "interval" and (o.ref_line is None or o.ref_line == o.line):
                continue                           # intra-line intervals carry no level information
            r = (o.value - pred(o, ds.unit)) * scale
            cols = basis(o.line)
            if o.kind == "interval":                   # f(line) − f(ref_line): both levels' corrections enter
                for c, v in basis(o.ref_line).items():
                    cols[c] = cols.get(c, 0.0) - v
            rows.append((ds.id, o.line, r, np.hypot(o.uncertainty * scale, floor), cols))
    names = sorted({c for *_, cols in rows for c in cols}, key=str)
    idx = {n: i for i, n in enumerate(names)}
    A = np.zeros((len(rows), len(names)))
    for i, (*_, cols) in enumerate(rows):
        for c, v in cols.items():
            A[i, idx[c]] = v
    r = np.array([x[2] for x in rows]); w = 1 / np.array([x[3] for x in rows]); sets = np.array([x[0] for x in rows])
    print(f"{len(rows)} lines, {len(names)} coefficients, floor {floor} MHz, prior {prior} MHz")

    def fit(mask, robust=True):
        """Weighted ridge; with robust=True, iteratively reweighted (Huber, 3 sigma_eff) so a few
        discordant lines do not steer the polynomial."""
        ww = w.copy()
        for _ in range(8 if robust else 1):
            Aw, rw = A[mask] * ww[mask, None], r[mask] * ww[mask]
            M = Aw.T @ Aw + np.eye(len(names)) / prior**2
            c = np.linalg.solve(M, Aw.T @ rw)
            z = np.abs((r - A @ c) * w)
            ww = w * np.minimum(1.0, 3.0 / np.maximum(z, 1e-9))
        return c

    rms = lambda v: np.sqrt(np.mean(np.square(v))) if len(v) else float("nan")   # noqa: E731
    c_all = fit(np.ones(len(rows), bool))
    print(f"\n{'set':16s} {'n':>3s} {'before':>8s} {'in-sample':>10s} {'held-out':>9s}")
    for s in SETS:
        m = sets == s
        c_out = fit(~m)
        print(f"{s:16s} {m.sum():3d} {rms(r[m]):8.3f} {rms(r[m] - A[m] @ c_all):10.3f} {rms(r[m] - A[m] @ c_out):9.3f}")
    m = np.ones(len(rows), bool)
    print(f"{'all':16s} {m.sum():3d} {rms(r):8.3f} {rms(r - A @ c_all):10.3f}")
    # leave-one-line-out: the test of interpolation within measured bands, which is what the model is for
    lines_ = np.array([str(x[1]) for x in rows])
    loo = np.zeros(len(rows))
    for ln in set(lines_):
        m = lines_ == ln
        loo[m] = r[m] - A[m] @ fit(~m)
    print(f"\nleave-one-line-out: {'set':16s} {'rms':>6s} {'median |r|':>10s}")
    for s_ in SETS:
        m = sets == s_
        print(f"{'':20s}{s_:16s} {rms(loo[m]):6.3f} {np.median(np.abs(loo[m])):10.3f}")
    print(f"{'':20s}{'all':16s} {rms(loo):6.3f} {np.median(np.abs(loo)):10.3f}")
    z = np.abs((r - A @ c_all) * w)
    print("\nlines more than 3 sigma_eff from the fit (in-sample):")
    for i in np.argsort(-z)[:8]:
        if z[i] > 3:
            print(f"  {rows[i][0]:15s} {str(rows[i][1]):22s} {rows[i][1] and ''}residual {r[i] - A[i] @ c_all:+7.2f} MHz  (sigma_eff {1/w[i]:.2f})")
    print("\ncoefficients (MHz):")
    for n, v in zip(names, c_all):
        print(f"  {n[0]} {'v' if FORM == 'perv' else 'u^'}{n[1]:2d} y^{n[2]}: {v:+8.3f}")
    if "--write" in opts:
        out = {"id": "level_corrections_2026a", "base": BASE, "unit": "MHz", "y": "J(J+1)/1e4",
               "note": ("Corrections to the term values of hannover2008/i2spec2026a for 127I2, fitted to the comb-"
                        "referenced NIR lines (liao2010a, bodermann2000a, bodermann1998b, reinhardt2007a, cornish2000a); "
                        "delta B(v'=0) = 0 by convention. prototypes/nir_corrections.py"),
               "floor_MHz": floor, "prior_MHz": prior,
               "X": {str(v): [float(c_all[idx[("X", v, k)]]) for k in range(DEGREE_X[v] + 1)] for v in X_V if ("X", v, 0) in idx},
               "B": {str(v): [float(c_all[idx[("B", v, k)]]) for k in range(DEGREE_B.get(v, 0) + 1)] for v in B_V if ("B", v, 0) in idx}}
        cov = defaultdict(lambda: [10**6, -1])
        for _, line, *_ in rows:
            Ju = line.J_lower + (1 if line.branch == "R" else -1)
            for st, v, J in (("X", line.v_lower, line.J_lower), ("B", line.v_upper, Ju)):
                if (st == "X" and v in X_V) or (st == "B" and v in B_V):
                    cov[(st, v)] = [min(cov[(st, v)][0], J), max(cov[(st, v)][1], J)]
        out["coverage"] = {st: {str(v): cov[(st, v)] for (s_, v) in sorted(cov) if s_ == st} for st in ("X", "B")}
        out["fit"] = {"lines": len(rows), "rms_before_MHz": float(rms(r)), "rms_in_sample_MHz": float(rms(r - A @ c_all)),
                      "leave_one_line_out_median_MHz": float(np.median(np.abs(loo))),
                      "leave_one_line_out_median_by_set_MHz": {s_: float(np.median(np.abs(loo[sets == s_]))) for s_ in SETS}}
        path = ROOT / "src/i2spec/data/level_corrections_2026a.json"
        path.write_text(json.dumps(out, indent=1))
        print("wrote", path)


if __name__ == "__main__":
    main(sys.argv[1:])
