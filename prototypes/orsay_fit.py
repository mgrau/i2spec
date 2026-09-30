"""X-state levels v'' = 17-25 from the assigned Orsay part I lines, with the atlas calibration.

Model for the residual r = obs - model (MHz) of an assigned line (v', v'', J''):
    r = delta * sigma - dX(v'', J'')          [delta: atlas scale error; dX: correction to the X term value]
The atlas's own page IV: its wavenumbers are absolute, Fourier-transform errors are proportional to
frequency, and calibration plus pointing should stay within +-0.005 cm-1 for strong lines. So the
calibration is one scale factor, not an offset and a slope.
    dX(v'', J'') = sum_k c_{v''k} y^k,   y = J''(J''+1)/1e4,   for v'' in FREE; zero below.
X v'' <= 16 is held, and a line there enters only if both its levels lie inside the J range the
comb-referenced corrections were fitted on (level_corrections_2026e coverage +- 15; B v' = 0 counts as
covered): outside it those levels are extrapolations, not anchors. Robust (Huber) weights; line sigma = hypot(eps, FLOOR).

Validation: fit on the v' = 0 lines only and predict the v' = 1 lines, whose upper levels the fit never
saw; and leave-one-line-out per level.

--freeB also frees the J-dependence of the upper levels: B v' = 0 (y'^1, y'^2; its constant is the
convention) and B v' = 1 (constant, y', y'^2). Needed because these lines reach J' ~ 265, where B v' = 1 is
known only to J' = 147 and B v' = 0's J-dependence lives in the NIR band correction, which is confined to
v'' = 12-17.

usage: orsay_fit.py [--free=17-25] [--floor=20] [--maxdeg=3] [--freeB] [--write=level_corrections_2026f]
"""
import csv, json, sys
from pathlib import Path
import numpy as np

MHZ = 29979.2458
SYSTEMATIC = 15.0      # MHz: spread of a level's v' = 0 / 1 / 2 sub-samples after the fit
ROOT = Path(__file__).resolve().parents[1]


def degree(Js, maxdeg=3):
    span = max(Js) - min(Js)
    d = 0 if len(Js) < 5 or span < 40 else 1 if len(Js) < 15 or span < 100 else 2 if len(Js) < 60 or span < 180 else 3
    return min(d, maxdeg)


def design(A, free, deg, freeB=False):
    cols = ["delta"] + [(v, k) for v in free for k in range(deg[v] + 1)]
    if freeB:
        cols += [("B", 0, 1), ("B", 0, 2), ("B", 1, 0), ("B", 1, 1), ("B", 1, 2)]
    X = np.zeros((len(A), len(cols)))
    ix = {c: i for i, c in enumerate(cols)}
    for i, a in enumerate(A):
        X[i, 0] = a["sigma"] * MHZ * 1e-9          # delta in ppb
        if a["vl"] in deg:
            for k in range(deg[a["vl"]] + 1):
                X[i, ix[(a["vl"], k)]] = -(a["y"] ** k)
        Ju = a["J"] + (1 if a["br"] == "R" else -1)
        for k in range(3):
            if ("B", a["vu"], k) in ix:
                X[i, ix[("B", a["vu"], k)]] = (Ju * (Ju + 1) / 1e4) ** k
    return X, cols


def fit(X, r, s, mask, prior=1e4, iters=10):
    w = 1 / s
    for _ in range(iters):
        Xw, rw = X[mask] * w[mask, None], r[mask] * w[mask]
        c = np.linalg.solve(Xw.T @ Xw + np.eye(X.shape[1]) / prior**2, Xw.T @ rw)
        z = np.abs(r - X @ c) / s
        w = np.minimum(1.0, 2.5 / np.maximum(z, 1e-9)) / s
    return c


def main(argv):
    opts = {a.split("=")[0]: a.split("=", 1)[1] for a in argv if "=" in a}
    lo, hi = (int(x) for x in opts.get("--free", "17-25").split("-"))
    floor = float(opts.get("--floor", 20))
    A = []
    for a in csv.DictReader(open(ROOT / "data/atlas_lines/orsay1982_part1_assigned.csv")):
        eps = float(a["eps_mk"]) * 1e-3 * MHZ if a["eps_mk"] else 150.0
        J = int(a["J_lower"])
        A.append(dict(N=int(a["N"]), sigma=float(a["sigma_cm1"]), vu=int(a["v_upper"]), vl=int(a["v_lower"]), br=a["branch"],
                      J=J, y=J * (J + 1) / 1e4, r=float(a["obs_minus_model_MHz"]), s=np.hypot(eps, floor)))
    A = [a for a in A if a["vl"] <= hi]
    cov = json.loads((ROOT / "src/i2spec/data/level_corrections_2026e.json").read_text())["coverage"]
    def covered(st, v, J):
        if st == "B" and v == 0:
            return True
        r_ = cov[st].get(str(v))
        return r_ is not None and r_[0] - 15 <= J <= r_[1] + 15
    n0 = len(A)
    freeB_ = "--freeB" in argv
    A = [a for a in A if a["vl"] >= lo or (covered("X", a["vl"], a["J"]) and
         (freeB_ and a["vu"] <= 1 or covered("B", a["vu"], a["J"] + (1 if a.get("br", "R") == "R" else -1))))]
    print(f"anchor filter: {n0 - len(A)} lines at v'' < {lo} dropped as outside the comb-measured J ranges")
    free = [v for v in range(lo, hi + 1) if sum(a["vl"] == v for a in A) >= 5]
    maxdeg = int(opts.get("--maxdeg", 3))
    deg = {v: degree([a["J"] for a in A if a["vl"] == v], maxdeg) for v in free}
    freeB = "--freeB" in argv
    X, cols = design(A, free, deg, freeB)
    r = np.array([a["r"] for a in A]); s = np.array([a["s"] for a in A])
    vl = np.array([a["vl"] for a in A]); vu = np.array([a["vu"] for a in A])
    c = fit(X, r, s, np.ones(len(A), bool))
    res = r - X @ c
    rms = lambda x: float(np.sqrt(np.mean(np.square(x))))     # noqa: E731
    mad = lambda x: float(1.4826 * np.median(np.abs(x - np.median(x))))  # noqa: E731
    print(f"{len(A)} lines, free v'' {free}, degrees {deg}, floor {floor} MHz")
    ppb = float(c[0])
    print(f"atlas scale error delta = {ppb:+.1f} ppb  ({ppb*1e-9*12000*MHZ:+.1f} MHz at 12 000 cm-1, "
          f"{ppb*1e-9*12000*1e3:+.2f} mk)")
    # held-out: fit v'=0 only, predict v'=1
    m0 = vu == 0
    c0 = fit(X, r, s, m0)
    pred = r - X @ c0
    print(f"\n{'v':>3} {'n':>4} {'before med/MAD':>16} {'after med/MAD':>15} {'v1 from v0-only fit (n, med, MAD)':>36} {'LOO rms':>8}")
    loo = res.copy()
    Ns = np.array([a["N"] for a in A])
    for i in range(len(A)):
        m = np.ones(len(A), bool); m[i] = False
        if vl[i] in deg:
            loo[i] = r[i] - X[i] @ fit(X, r, s, m, iters=4)
    for v in sorted(set(vl)):
        m = vl == v
        m1 = m & (vu == 1)
        h = pred[m1]
        extra = f"{m1.sum():4d} {np.median(h):+6.1f} {mad(h):5.1f}" if m1.sum() >= 5 else ""
        print(f"{v:3d} {m.sum():4d} {np.median(r[m]):+8.1f} {mad(r[m]):6.1f} {np.median(res[m]):+8.1f} {mad(res[m]):6.1f}   {extra:>36} {rms(np.clip(loo[m], -300, 300)):8.1f}")
    print(f"all: before MAD {mad(r):.1f}, after {mad(res):.1f} MHz")
    print("\nresidual median by (v'', v') after the fit, MHz (n):")
    for v in sorted(set(vl)):
        print(f"   v''={v:2d}: " + "  ".join(f"v'={u}: {np.median(res[(vl==v)&(vu==u)]):+6.1f} ({((vl==v)&(vu==u)).sum()})"
                                          for u in (0, 1, 2) if ((vl == v) & (vu == u)).sum() >= 8))
    if freeB:
        print("B terms:", {c_: round(float(c[cols.index(c_)]), 1) for c_ in cols if isinstance(c_, tuple) and c_[0] == "B"})
    fr = np.isin(vl, free)
    print(f"LOO rms over the free levels: {rms(np.clip(loo[fr], -300, 300)):.2f} MHz  (maxdeg {maxdeg})")
    # dX values at a few J for the record
    print("\ncorrection to the X term values, dX (MHz):")
    for v in free:
        cc = [c[cols.index((v, k))] for k in range(deg[v] + 1)]
        Js = [a["J"] for a in A if a["vl"] == v]
        print(f"   v''={v}: " + "  ".join(f"J={J}: {np.polyval(cc[::-1], J*(J+1)/1e4):+6.1f}" for J in (min(Js), 100, 150, max(Js))) + f"   (J {min(Js)}-{max(Js)})")
    if "--write" in opts:
        name = opts["--write"]
        base = json.loads((ROOT / "src/i2spec/data/level_corrections_2026e.json").read_text())
        base["id"] = name
        for v in free:
            if v <= 16:
                continue
            cc = [float(c[cols.index((v, k))]) for k in range(deg[v] + 1)]
            old = base["X"].get(str(v))
            if old is not None:           # an existing correction (v'' = 17): add, keep the comb-fixed part
                cc = [a + (old[k] if k < len(old) else 0.0) for k, a in enumerate(cc)] + old[len(cc):]
            Js = [a["J"] for a in A if a["vl"] == v]
            base["X"][str(v)] = cc
            base["coverage"]["X"][str(v)] = [int(min(Js)), int(max(Js))]
            # the model's uncertainty, not the atlas's: held-out rms with the atlas lines' own scatter
            # taken out, floored at the 15 MHz by which the v' = 0, 1, 2 sub-samples of a level disagree
            lo_ = rms(np.clip(loo[vl == v], -300, 300)); sa = float(np.median(s[vl == v]))
            base["held_out_MHz"]["X"][str(v)] = max(float(np.sqrt(max(lo_**2 - sa**2, 0.0))), SYSTEMATIC)
        base["orsay"] = {"note": ("X v'' = 18-25 from the assigned lines of the Orsay atlas part I (Gerstenkorn, "
                                  "Verges & Chevillard 1982), prototypes/orsay_fit.py " + " ".join(a for a in argv if not a.startswith("--write"))
                                  + ". Held-out figures for these levels are the leave-one-line-out rms with the atlas lines' own "
                                  "scatter removed in quadrature, floored at 15 MHz (the disagreement between a level's v' = 0, 1 "
                                  "and 2 sub-samples)."),
                         "calibration_scale_ppb": ppb,
                         "free": free}
        (ROOT / "src/i2spec/data" / f"{name}.json").write_text(json.dumps(base, indent=1))
        print("wrote", name)


if __name__ == "__main__":
    main(sys.argv[1:])
