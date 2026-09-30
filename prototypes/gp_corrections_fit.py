"""The Gaussian-process B-state correction that i2spec2026n applies to visible B levels without their own.

Usage:  uv run --group research python prototypes/gp_corrections_fit.py [--validate] [--write=gp_b_2026n]
Output: src/i2spec/data/<name>.npz and <name>.json (with --write); prototypes/out/gp_corrections_fit.json

The bake-off (docs/research/model-bakeoff.md, prototypes/bakeoff.py scheme B) found that Gaussian-process
correction surfaces f_B(v', y), f_X(v'', y), y = J(J+1)/1e4, predict a held-out visible B level (v' = 3-43)
to 1.66 MHz rms, against 2.75 for no correction, and the 76 absolute rows of the later sets to 1.55 MHz
(4.85 for the per-level polynomials). Here the same model is fitted in the gauge of the level corrections:
X v'' <= 10 and B v' = 0 carry no correction (their columns are removed), so f_B at an uncorrected level adds
to the published curves exactly as a level correction would. --validate repeats the leave-one-level-out
test in this gauge (every B level v' = 3-50 that >= 2 rows touch; hyperparameters refitted per fold).

The package keeps only what the prediction at a new B level needs: the B training levels (v', J'),
beta = S_B^T alpha (alpha = K^-1 r over the rows) and M = S_B^T K^-1 S_B, with the hyperparameters, so that
    mean(v, J)  = k_B((v, y), L) . beta
    var(v, J)   = k_B((v, y), (v, y)) - k_B((v, y), L) M k_B(L, (v, y))
(src/i2spec/gp_corrections.py). Rows: prototypes/out/bakeoff_rows.pkl (run `bakeoff.py rows` first).
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

for _k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_k, "2")

import numpy as np
import scipy.sparse as sp

sys.path.insert(0, str(Path(__file__).parent))
import bakeoff  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
V_MIN_B, V_MAX_B = 3, 50


def gauged_data():
    D = bakeoff.Data()
    fixed = ((D.state == "X") & (D.v <= 10)) | ((D.state == "B") & (D.v == 0))
    D.S = (D.S @ sp.diags((~fixed).astype(float))).tocsr()
    D.Q = (D.Q @ sp.diags((~fixed).astype(float))).tocsr()
    return D


def validate(D, theta):
    """Leave one B level out: fit on the rows that do not touch it, predict those that do."""
    out = []
    levels = sorted({(s, v) for t in D.touch for (s, v) in t if s == "B" and V_MIN_B <= v <= V_MAX_B})
    for s, v in levels:
        test = np.array([(s, v) in t for t in D.touch])
        if test.sum() < 2:
            continue
        g = bakeoff.SchemeB(D).fit(~test, theta0=theta)
        mu, sd = g.predict(np.flatnonzero(test))
        res = D.r[test] - mu
        out.append(dict(v=v, n=int(test.sum()), rms=float(np.sqrt(np.mean(res**2))),
                        bare=float(np.sqrt(np.mean(D.r[test] ** 2))),
                        within1=float(np.mean(np.abs(res) <= sd)), z_rms=float(np.sqrt(np.mean((res / sd) ** 2)))))
        print(f"  B v' = {v:2d}: {out[-1]['n']:3d} rows, bare {out[-1]['bare']:6.2f}, GP {out[-1]['rms']:6.2f} MHz, "
              f"within 1 sigma {out[-1]['within1']:.2f}", flush=True)
    return out


def main(argv):
    opts = {a.split("=")[0]: (a.split("=", 1)[1] if "=" in a else True) for a in argv}
    D = gauged_data()
    g = bakeoff.SchemeB(D).fit(np.ones(D.N, bool))
    th = g.th
    print("hyperparameters:", dict(zip(bakeoff.SchemeB.NAMES, np.round(np.exp(th), 4))), flush=True)
    report = dict(theta=dict(zip(bakeoff.SchemeB.NAMES, map(float, np.exp(th)))))
    if "--validate" in opts:
        val = validate(D, th)
        vis = [x for x in val if x["v"] <= 43]
        pool = lambda xs, k: float(np.sqrt(sum(x[k] ** 2 * x["n"] for x in xs) / sum(x["n"] for x in xs)))  # noqa: E731
        report["validation"] = dict(levels=val, visible_rms=pool(vis, "rms"), visible_bare=pool(vis, "bare"),
                                    visible_within1=float(sum(x["within1"] * x["n"] for x in vis) / sum(x["n"] for x in vis)))
        print(f"visible B v' = 3-43, held out by level: bare {report['validation']['visible_bare']:.2f}, "
              f"GP {report['validation']['visible_rms']:.2f} MHz rms, "
              f"{report['validation']['visible_within1'] * 100:.0f} % within 1 sigma", flush=True)
    # the B block of the posterior
    from scipy.linalg import cho_solve
    P = g.st["B"]
    iiB = P["ii"]
    SB = P["S"][g.tr]                                   # rows x B levels
    used = np.flatnonzero(np.asarray(abs(SB).sum(axis=0)).ravel() > 0)
    SBu = SB[:, used].toarray()
    beta = SBu.T @ g.al
    M = SBu.T @ cho_solve(g.cf, SBu)
    vB, JB = D.v[iiB][used], D.J[iiB][used]
    o = 0
    a1, lv, ly, a2, ly2 = map(float, np.exp(th[o:o + 5]))
    meta = dict(description="Gaussian-process B-state correction (MHz) to the bare potentials of i2spec2026n, in the gauge "
                            "of the level corrections (X v'' <= 10 and B v' = 0 fixed at zero). Kernel: a1^2 Matern52(dv/len_v) "
                            "SE(dy/len_y) + a2^2 [same v] SE(dy/len_y2), y = J(J+1)/1e4. prototypes/gp_corrections_fit.py.",
                a1=a1, len_v=lv, len_y=ly, a2=a2, len_y2=ly2, tau=float(np.exp(th[-1])),
                v_min=3, v_max=35,
                apply_to="B v' = 3-35 levels with no level correction of their own (held out, v' = 37-43 miss by 2.5-15 MHz)",
                n_rows=int(g.tr.size), n_levels=int(used.size), validation=report.get("validation"))
    report["n_levels"] = int(used.size)
    (ROOT / "prototypes/out").mkdir(exist_ok=True)
    (ROOT / "prototypes/out/gp_corrections_fit.json").write_text(json.dumps(report, indent=1))
    if "--write" in opts:
        name = opts["--write"]
        np.savez_compressed(ROOT / f"src/i2spec/data/{name}.npz", v=vB.astype(np.int16), J=JB.astype(np.int16),
                            beta=beta, M=M)
        (ROOT / f"src/i2spec/data/{name}.json").write_text(json.dumps(meta, indent=1))
        print(f"wrote {name}: {used.size} B levels")


if __name__ == "__main__":
    main(sys.argv[1:])
