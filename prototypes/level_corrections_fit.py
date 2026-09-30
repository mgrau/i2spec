"""Level corrections from every comb-referenced ¹²⁷I₂ line: δX(v'' >= 11, J) and δB(v', J) on the potentials.

Generalises prototypes/nir_corrections.py to the whole range. Each corrected level gets a polynomial in
y = J(J+1)/1e4 whose degree follows the J coverage of the precise lines (sigma <= SIGMA_PRECISE) that
reach it; Doppler-limited lines (sansonetti, velchev, xu) enter with their weights but create no
parameters. Conventions: δX = 0 for v'' <= 10 and δB(v' = 0) = 0, so B carries the v'' = 0 bands'
error and X the NIR's. Robust ridge regression on the residuals against BASE, leave-one-line-out
validation, and a held-out rms per level that lookup.uncertainty uses.

usage: level_corrections_fit.py [--floor=0.3] [--prior=5] [--base=i2spec2026a|pair] [--extended]
                                [--write=level_corrections_2026b]
--base=pair puts the corrections on the MLR pair of data/potentials/*_2026c.json instead of the published
potentials; --extended lets the fit reach v' <= 62 and v'' <= 54 (only sensible on the pair).
"""
from __future__ import annotations

import json
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

import numpy as np

from i2spec.constants import MHZ_PER_CM
from i2spec.observations import Predictor, load_all


class BarePredictor(Predictor):
    """A Predictor whose models carry no level corrections, whatever the parameter set names: the
    corrections are what this script fits, and a base that already has them fits nothing."""

    def model(self, isotopologue):
        m = super().model(isotopologue)
        if m.corrections is not None:
            m.corrections = None
            m._levels.clear()
        return m

ROOT = Path(__file__).resolve().parents[1]
BASE = "i2spec2026a"          # bare potentials + the V_ad constant; never the default, which may carry corrections
SIGMA_PRECISE = 0.3           # MHz: lines that may create a correction
SIGMA_PERTURBED = 5.0         # above v' = 50 (the 1g region, GHz errors) any line measured to this creates one
V_PERTURBED = 50
X_MIN_V = 11                  # X corrections only from here up (below, the potentials are the reference)
V_MAX_B, V_MAX_X = 43, 17     # the region where the potentials are right to MHz; beyond it the MLR pair takes over


def degree(Js):
    """Polynomial degree in y a level can support: from the distinct J of its precise lines."""
    Js = sorted(set(Js))
    if len(Js) < 3 or Js[-1] - Js[0] < 40:
        return 0
    if len(Js) < 5 or Js[-1] - Js[0] < 100:
        return 1
    if len(Js) < 10 or Js[-1] - Js[0] < 150:
        return 2
    return 3


def main(argv):
    opts = {a.split("=")[0]: a.split("=", 1)[1] if "=" in a else True for a in argv}
    floor, prior = float(opts.get("--floor", 0.3)), float(opts.get("--prior", 5.0))
    base = opts.get("--base", BASE)
    if base == "pair":
        from mlr_evaluate import load_x, load_b, GRID_B
        pred = BarePredictor(BASE, potentials={"X": load_x(ROOT / "data/potentials/mlr_x_2026c.json"),
                                               "B": load_b(ROOT / "data/potentials/mlr_b_2026c.json")}, grids=GRID_B)
    else:
        pred = BarePredictor(base)
    v_max_b, v_max_x = (62, 54) if "--extended" in opts else (V_MAX_B, V_MAX_X)
    rows = []                 # (set, line, residual MHz, sigma_eff MHz, {level: (state, v, J)} with signs)
    for ds in load_all():
        scale = MHZ_PER_CM if ds.unit == "cm-1" else 1.0
        for o in ds.observations:
            ref = o.ref_line or o.line
            if o.line.isotopologue != "127I2" or ref.isotopologue != "127I2":
                continue
            if o.kind == "interval" and ref == o.line:
                continue                                   # intra-line: hyperfine only
            if max(o.line.v_upper, ref.v_upper) > v_max_b or max(o.line.v_lower, ref.v_lower) > v_max_x:
                continue
            r = (o.value - pred(o, ds.unit)) * scale
            levels = defaultdict(float)                    # (state, v, J) -> sign
            for line, sign in ((o.line, +1.0), (ref, -1.0 if o.kind == "interval" else 0.0)):
                if sign == 0.0:
                    continue
                Ju = line.J_lower + (1 if line.branch == "R" else -1)
                levels[("B", line.v_upper, Ju)] += sign
                levels[("X", line.v_lower, line.J_lower)] -= sign
            rows.append((ds.id, o.line, r, np.hypot(o.uncertainty * scale, floor), o.uncertainty * scale, dict(levels)))
    # which levels get parameters, and of what degree
    precise = defaultdict(list)
    for _, line, _, _, sig, levels in rows:
        if sig <= SIGMA_PRECISE or (sig <= SIGMA_PERTURBED and line.v_upper > V_PERTURBED):
            for (st, v, J), s in levels.items():
                if s != 0 and ((st == "X" and v >= X_MIN_V) or (st == "B" and v >= 1)):
                    precise[(st, v)].append(J)
    deg = {k: degree(Js) for k, Js in precise.items()}
    names = sorted({(st, v, k) for (st, v), d in deg.items() for k in range(d + 1)}, key=str)
    idx = {n: i for i, n in enumerate(names)}
    A = np.zeros((len(rows), len(names)))
    for i, (_, _, _, _, _, levels) in enumerate(rows):
        for (st, v, J), s in levels.items():
            if (st, v) in deg:
                y = J * (J + 1) / 1e4
                for k in range(deg[(st, v)] + 1):
                    A[i, idx[(st, v, k)]] += s * y**k
    r = np.array([x[2] for x in rows]); w = 1 / np.array([x[3] for x in rows]); sets = np.array([x[0] for x in rows])
    used = A.any(axis=1)
    print(f"{len(rows)} rows, {used.sum()} touch a corrected level; {len(deg)} levels, {len(names)} coefficients; "
          f"floor {floor} MHz, prior {prior} MHz")

    # the ridge prior per coefficient: `prior` where the potentials are right to MHz, wide where they
    # are known to be off by tens of MHz (the extended region) or GHz (the 1g region)
    def prior_of(name):
        st, v, _ = name
        if st == "B" and v > V_PERTURBED:
            return 1e5
        if (st == "B" and v > V_MAX_B) or (st == "X" and v > V_MAX_X):
            return 100.0
        return prior
    ridge = np.diag([1 / prior_of(n)**2 for n in names])

    def fit(mask, robust=True, cov=False):
        ww = w.copy()
        for _ in range(8 if robust else 1):
            Aw, rw = A[mask] * ww[mask, None], r[mask] * ww[mask]
            c = np.linalg.solve(Aw.T @ Aw + ridge, Aw.T @ rw)
            z = np.abs((r - A @ c) * w)
            ww = w * np.minimum(1.0, 3.0 / np.maximum(z, 1e-9))
        if cov:        # the coefficients' covariance at the final robust weights, scaled by the reduced chi^2
            Aw = A[mask] * ww[mask, None]
            C = np.linalg.inv(Aw.T @ Aw + ridge)
            dof = max(int(mask.sum()) - A.shape[1], 1)
            chi2 = float(np.sum(((r[mask] - A[mask] @ c) * ww[mask]) ** 2)) / dof
            return c, C * max(chi2, 1.0)
        return c

    rms = lambda v: float(np.sqrt(np.mean(np.square(v)))) if len(v) else float("nan")   # noqa: E731
    c_all, C_all = fit(np.ones(len(rows), bool), cov=True)
    lines_ = np.array([str(x[1]) for x in rows])
    loo = r - A @ c_all
    for ln in set(lines_[used]):
        m = lines_ == ln
        loo[m] = r[m] - A[m] @ fit(~m)
    print(f"\n{'set':16s} {'rows':>4s} {'before':>7s} {'in-sample':>9s} {'leave-one-line-out rms / median':>32s}")
    for s_ in sorted(set(sets)):
        m = (sets == s_) & used
        if m.any():
            print(f"{s_:16s} {m.sum():4d} {rms(r[m]):7.2f} {rms(r[m] - A[m] @ c_all):9.2f} {rms(loo[m]):16.2f} / {np.median(np.abs(loo[m])):6.2f}")
    m = used
    print(f"{'all corrected':16s} {m.sum():4d} {rms(r[m]):7.2f} {rms(r[m] - A[m] @ c_all):9.2f} {rms(loo[m]):16.2f} / {np.median(np.abs(loo[m])):6.2f}")
    # per level: coefficients, coverage, held-out rms of the lines that reach it
    per_level = {}
    sig = np.array([x[4] for x in rows])
    # a level's held-out figure comes from the precise lines that reach it whose other level is inside the
    # well-determined region: a B level is not judged by lines from X v'' >= 18, whose residual is the X
    # hyperfine model there, nor an X level by lines to B v' >= 44
    def partner_ok(lv, st):
        return all((k0 != "X" or k1 <= V_MAX_X) if st == "B" else (k0 != "B" or k1 <= V_MAX_B) for (k0, k1, _), s_ in lv.items() if s_ != 0)
    for (st, v), d in sorted(deg.items()):
        touch = np.array([any(k == (st, v) for (k0, k1, _), s in lv.items() for k in [(k0, k1)] if s != 0) and partner_ok(lv, st) for *_, lv in rows]) \
            & (sig <= (SIGMA_PERTURBED if (st == "B" and v > V_PERTURBED) else SIGMA_PRECISE))
        Js = precise[(st, v)]
        n_lines = int(len(set(lines_[touch])))
        ids = [idx[(st, v, k)] for k in range(d + 1)]
        cov_v = C_all[np.ix_(ids, ids)]
        # the part of the lines' held-out scatter that neither their own uncertainty nor the coefficients'
        # covariance explains: the model-discrepancy term of this level (roadmap item 3)
        if n_lines >= 2:
            pv = np.array([sum(sg * (J * (J + 1) / 1e4) ** k * (J * (J + 1) / 1e4) ** l * cov_v[k, l]
                               for (k0, k1, J), sg in lv.items() if (k0, k1) == (st, v) for k in range(d + 1) for l in range(d + 1))
                           for *_, lv in (rows[i] for i in np.flatnonzero(touch))])
            md = float(np.sqrt(max(rms(loo[touch]) ** 2 - np.mean(sig[touch] ** 2 + np.abs(pv)), 0.0)))
        else:
            md = None
        per_level[(st, v)] = dict(covariance=cov_v.tolist(), discrepancy=md,
                                  coefficients=[float(c_all[idx[(st, v, k)]]) for k in range(d + 1)],
                                  coverage=[int(min(Js)), int(max(Js))], lines=n_lines,
                                  loo_rms_MHz=rms(loo[touch]) if n_lines >= 2 else None,   # one line: no validation
                                  loo_median_MHz=float(np.median(np.abs(loo[touch]))))
    print("\nlevel   deg  lines  J range      c0      LOO rms / median (MHz), precise lines only")
    for (st, v), d in per_level.items():
        print(f"{st} {v:2d}   {len(d['coefficients'])-1}   {d['lines']:4d}  {d['coverage'][0]:3d}-{d['coverage'][1]:3d}  {d['coefficients'][0]:+7.2f}   "
              + (f"{d['loo_rms_MHz']:6.2f} / {d['loo_median_MHz']:5.2f}" if d['loo_rms_MHz'] is not None else "  (one line)"))
    if "--write" in opts:
        name = opts["--write"] if isinstance(opts["--write"], str) else "level_corrections_2026b"
        out = {"id": name, "base": base, "unit": "MHz", "y": "J(J+1)/1e4",
               "note": ("Corrections to the 127I2 term values of the base set, fitted to every comb-referenced line "
                        "(sigma <= 0.3 MHz creates a level's parameters; Doppler-limited sets weigh in). delta X = 0 "
                        "below v'' = 11 and delta B(v' = 0) = 0 by convention. prototypes/level_corrections_fit.py"),
               "floor_MHz": floor, "prior_MHz": prior,
               "X": {str(v): d["coefficients"] for (st, v), d in per_level.items() if st == "X"},
               "B": {str(v): d["coefficients"] for (st, v), d in per_level.items() if st == "B"},
               "coverage": {st: {str(v): d["coverage"] for (s_, v), d in per_level.items() if s_ == st} for st in ("X", "B")},
               "held_out_MHz": {st: {str(v): d["loo_rms_MHz"] for (s_, v), d in per_level.items() if s_ == st and d["loo_rms_MHz"] is not None} for st in ("X", "B")},
               "covariance_MHz2": {st: {str(v): d["covariance"] for (s_, v), d in per_level.items() if s_ == st} for st in ("X", "B")},
               "discrepancy_MHz": {st: {str(v): d["discrepancy"] for (s_, v), d in per_level.items() if s_ == st and d["discrepancy"] is not None} for st in ("X", "B")},
               "fit": {"rows": int(used.sum()), "rms_before_MHz": rms(r[used]), "rms_in_sample_MHz": rms(r[used] - A[used] @ c_all),
                       "leave_one_line_out_rms_MHz": rms(loo[used]), "leave_one_line_out_median_MHz": float(np.median(np.abs(loo[used])))}}
        path = ROOT / "src/i2spec/data" / f"{name}.json"
        path.write_text(json.dumps(out, indent=1))
        print("wrote", path)


if __name__ == "__main__":
    main(sys.argv[1:])
