"""A local correction for the NIR bands v' = 0 -> v'' = 12-17, on top of the global level corrections.

Why local, and why not more level corrections: the residual that remains on the comb-referenced NIR
lines after level_corrections_2026c is R-high / P-low at the same J'' (R(106) 0-13 +0.70, R(100) 0-12
+0.59, P(82) 0-13 -0.39 MHz on liao2010a). No function of the lower level can absorb that; it needs a
J'-dependence of B v' = 0. Fitting one as a *level* correction works on the NIR lines and moves every
band that shares only one side of it - 0-9 and 0-10 by +13 MHz at J = 100 and +47 at J = 200, the v' >= 1
bands to v'' = 11-17 by -10 to -130 MHz - because the common J-dependence of B(v'=0) and X(v''=12-17)
is invisible to the NIR lines and there are no data elsewhere to fix it. IodineSpec5 has the same
structure and confines it the same way: its local NIR Dunham replaces the potential only for these
bands (docs/research/iodinespec5.md, "The NIR loss").

So this fits, to the residuals the default model leaves on lines inside the domain,

    delta nu = sum_{k>=1} b_k y'^k  -  sum_{l,k} a_lk (v'' - CENTRE)^l y''^k,      y = J(J+1)/1e4

and the model adds it only to v' = 0, v'' = 12-17 lines within J_MARGIN of the J range measured at that
v''. Outside, nothing changes.

Liao 2010 is treated as the paper says: +114 kHz on every row (the cell pressure shift, Table 2
footnote b), and the printed per-row measurement uncertainty with LIAO_RESIDUAL_SHIFT for the
transition dependence the paper did not study - not the 200 kHz overall accuracy, which is mostly that
one common offset.

usage: nir_band_fit.py [--parameters=SET] [--degree=L,K,KB] [--floor=0.05] [--prior=5] [--write=NAME]
       Fits against SET (default: the package default) with any band section of its corrections removed.
       --write copies SET's level corrections and adds this as their "band" section.
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import numpy as np

from dataclasses import replace

from i2spec.constants import DEFAULT_PARAMETERS, MHZ_PER_CM
from i2spec.level_corrections import J_MARGIN
from i2spec.observations import Predictor, load_all
from i2spec.potentials import parameter_set

ROOT = Path(__file__).resolve().parents[1]
V_UPPER, V_LOWER = 0, (12, 17)
CENTRE = 14
SIGMA_PRECISE = 0.3
LIAO_SHIFT = 0.0              # MHz: the loader now applies Liao 2010's +114 kHz (Table 2 footnote b) from its [shift_correction]
LIAO_RESIDUAL_SHIFT = 0.03    # MHz: the three calibration lines' shifts scatter over 102-121 kHz


def in_domain(line):
    return (line.isotopologue == "127I2" and line.v_upper == V_UPPER
            and V_LOWER[0] <= line.v_lower <= V_LOWER[1])


def features(line, L, K, KB):
    """d(nu)/d(parameter) for one line: B(v'=0) terms, then the X surface a_lk."""
    Ju = line.J_lower + (1 if line.branch == "R" else -1)
    yu, yl, dv = Ju * (Ju + 1) / 1e4, line.J_lower * (line.J_lower + 1) / 1e4, line.v_lower - CENTRE
    return np.array([yu**k for k in range(1, KB + 1)]
                    + [-(dv**l) * yl**k for l in range(L + 1) for k in range(K + 1)])


def rows(pred, L, K, KB, floor, liao=True):
    out = []
    for ds in load_all():
        scale = MHZ_PER_CM if ds.unit == "cm-1" else 1.0
        for o in ds.observations:
            ref = o.ref_line or o.line
            if o.kind == "interval" and ref == o.line:
                continue                                        # intra-line: hyperfine only
            if not in_domain(o.line) or (o.kind == "interval" and not in_domain(ref)):
                continue
            r = (o.value - pred(o, ds.unit)) * scale
            sigma = o.uncertainty * scale
            if liao and ds.id == "liao2010a":
                r += LIAO_SHIFT
                m = re.search(r"uncertainty (\d+) kHz", o.note)
                sigma = np.hypot(float(m.group(1)) / 1e3, LIAO_RESIDUAL_SHIFT) if m else sigma
            f = features(o.line, L, K, KB)
            if o.kind == "interval":
                f = f - features(ref, L, K, KB)
            out.append(dict(set=ds.id, line=str(o.line), obj=o.line, r=r, sigma=sigma,
                            w=1 / np.hypot(sigma, floor), f=f))
    return out


def coverage_of(R):
    """Where the correction may be applied: the J'' range of the precise lines at each v''."""
    coverage = {}
    for x in R:
        if x["sigma"] <= SIGMA_PRECISE:
            v, J = x["obj"].v_lower, x["obj"].J_lower
            lo, hi = coverage.get(v, (J, J))
            coverage[v] = (min(lo, J), max(hi, J))
    return coverage


def fit(A, r, w, ridge, mask, robust=True):
    ww = w.copy()
    for _ in range(8 if robust else 1):
        Aw, rw = A[mask] * ww[mask, None], r[mask] * ww[mask]
        c = np.linalg.solve(Aw.T @ Aw + ridge, Aw.T @ rw)
        z = np.abs((r - A @ c) * w)
        ww = w * np.minimum(1.0, 3.0 / np.maximum(z, 1e-9))
    return c


def main(argv):
    opts = {a.split("=")[0]: a.split("=", 1)[1] if "=" in a else True for a in argv}
    L, K, KB = (int(x) for x in opts.get("--degree", "3,3,3").split(","))
    floor, prior = float(opts.get("--floor", 0.05)), float(opts.get("--prior", 5.0))
    pset = opts.get("--parameters", DEFAULT_PARAMETERS)
    pred = Predictor(pset)                               # global corrections included, any band section removed
    m127 = pred.model("127I2")
    if m127.corrections is not None and m127.corrections.band is not None:
        m127.corrections = replace(m127.corrections, band=None)
    R = rows(pred, L, K, KB, floor)
    A = np.array([x["f"] for x in R]); r = np.array([x["r"] for x in R]); w = np.array([x["w"] for x in R])
    sets = np.array([x["set"] for x in R]); lines = np.array([x["line"] for x in R])
    sig = np.array([x["sigma"] for x in R])
    ridge = np.eye(A.shape[1]) / prior**2
    c = fit(A, r, w, ridge, np.ones(len(R), bool))
    # held out a line at a time, under the rule the model ships with: the correction reaches a line only
    # inside the J range the *remaining* precise lines cover at its v'' (+- J_MARGIN); outside, none
    loo = r - A @ c
    for ln in set(lines):
        m = lines == ln
        cov = coverage_of([x for x, keep in zip(R, ~m) if keep])
        held = R[int(np.nonzero(m)[0][0])]["obj"]
        lo, hi = cov.get(held.v_lower, (None, None))
        inside = lo is not None and lo - J_MARGIN <= held.J_lower <= hi + J_MARGIN
        loo[m] = r[m] - (A[m] @ fit(A, r, w, ridge, ~m) if inside else 0.0)
    rms = lambda v: float(np.sqrt(np.mean(np.square(v)))) if len(v) else float("nan")   # noqa: E731
    print(f"{len(R)} rows in v'={V_UPPER}, v''={V_LOWER[0]}-{V_LOWER[1]} ({len(set(lines))} lines); "
          f"degree L,K,KB = {L},{K},{KB} ({A.shape[1]} parameters); floor {floor}, prior {prior} MHz")
    print(f"{'set':16s} {'rows':>4s} {'before':>7s} {'in-sample':>9s} {'held-out rms / median':>24s}")
    for s in sorted(set(sets)):
        m = sets == s
        print(f"{s:16s} {m.sum():4d} {rms(r[m]):7.3f} {rms(r[m] - A[m] @ c):9.3f} {rms(loo[m]):14.3f} / {np.median(np.abs(loo[m])):6.3f}")
    m = sig <= SIGMA_PRECISE
    print(f"{'NIR precise':16s} {m.sum():4d} {rms(r[m]):7.3f} {rms(r[m] - A[m] @ c):9.3f} {rms(loo[m]):14.3f} / "
          f"{np.median(np.abs(loo[m])):6.3f}   (max |held-out| {np.abs(loo[m]).max():.2f})")

    coverage = coverage_of(R)
    print("coverage (J''):", dict(sorted(coverage.items())), f"+- {J_MARGIN}")

    if "--write" in opts:
        name = opts["--write"]
        base_name = parameter_set(pset)["level_corrections"]
        out = json.loads((ROOT / "src/i2spec/data" / f"{base_name}.json").read_text())
        out.pop("band", None)
        out["id"] = name
        out["band"] = {
            "note": ("Local correction for the NIR bands, fitted to the residuals the global corrections leave "
                     "(prototypes/nir_band_fit.py " + " ".join(a for a in argv if not a.startswith("--write"))
                     + "). Applied only to v' = 0 -> v'' in v_lower lines whose J'' lies within J_MARGIN of the "
                     "coverage at that v''. delta nu (MHz) = sum_k B[k-1] y'^k - sum_lk X[l][k] (v'' - centre)^l y''^k, "
                     "y = J(J+1)/1e4. liao2010a enters with its +114 kHz pressure shift applied (by the loader, [shift_correction])."),
            "based_on": base_name, "v_upper": V_UPPER, "v_lower": list(V_LOWER), "centre": CENTRE,
            "B": [float(x) for x in c[:KB]],
            "X": [[float(c[KB + l * (K + 1) + k]) for k in range(K + 1)] for l in range(L + 1)],
            "coverage": {str(v): list(map(int, j)) for v, j in sorted(coverage.items())},
            "held_out_MHz": rms(loo[m]),
            "held_out_median_MHz": float(np.median(np.abs(loo[m]))),
        }
        path = ROOT / "src/i2spec/data" / f"{name}.json"
        path.write_text(json.dumps(out, indent=1))
        print("wrote", path)


if __name__ == "__main__":
    main(sys.argv[1:])
