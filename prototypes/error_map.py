"""Lever 5, the discrepancy half: what the current default model actually achieves, by region.

Usage:  uv run python prototypes/error_map.py          (about 5 minutes)
Output: prototypes/out/error_map.json

For every observation that constrains a line position (frequencies and inter-line intervals) the
residual obs - model is computed with the default model, and the residuals are pooled by isotopologue,
v' band and v'' band. The rms in each cell, with the count and the data floor, is what a per-line
uncertainty can honestly quote today, and it is compared with the rule in lookup.uncertainty.
"""

import json
from collections import defaultdict
from pathlib import Path

import numpy as np

from i2spec.constants import MHZ_PER_CM
from i2spec.lookup import Line, uncertainty
from i2spec.observations import Predictor, load_all

ROOT = Path(__file__).resolve().parents[1]
V_UPPER_BANDS = ((0, 17), (18, 30), (31, 44), (45, 53), (54, 70))
V_LOWER_BANDS = ((0, 5), (6, 12), (13, 17), (18, 47), (48, 60))
GRID_B = {"B": dict(rmin=2.35, rmax=8.0, h=0.01, order=10, nlev=70)}


def band(v, bands):
    for lo, hi in bands:
        if lo <= v <= hi:
            return f"{lo}-{hi}"
    return "beyond"


def main():
    pred = Predictor(grids=GRID_B)
    rows = []
    for ds in load_all():
        scale = MHZ_PER_CM if ds.meta["unit"] == "cm-1" else 1.0
        for o in ds.observations:
            ref = o.ref_line or o.line
            if o.kind == "interval" and ref == o.line:
                continue
            r = (o.value - pred(o, "MHz")) * scale if ds.meta["unit"] == "MHz" else (o.value - pred(o, "cm-1")) * MHZ_PER_CM
            L = o.line
            rule = uncertainty(Line(L.isotopologue, L.branch, L.J_lower, L.v_upper, L.v_lower,
                                    pred.position(L) / MHZ_PER_CM, 1e-20, 0.0, 293.15))[0]
            rows.append(dict(dataset=ds.id, line=str(L), iso=L.isotopologue, v_upper=L.v_upper, v_lower=L.v_lower,
                             J=L.J_lower, residual=float(r), sigma=float(o.uncertainty * scale), rule=float(rule),
                             kind=o.kind))
    cells = defaultdict(list)
    for r in rows:
        cells[(r["iso"], band(r["v_upper"], V_UPPER_BANDS), band(r["v_lower"], V_LOWER_BANDS))].append(r)
    table = []
    for (iso, vu, vl), rs in sorted(cells.items()):
        res = np.array([r["residual"] for r in rs])
        table.append(dict(iso=iso, v_upper=vu, v_lower=vl, n=len(rs), lines=len({r["line"] for r in rs}),
                          sets=sorted({r["dataset"] for r in rs}), rms=float(np.sqrt(np.mean(res ** 2))),
                          median_abs=float(np.median(np.abs(res))), p90_abs=float(np.percentile(np.abs(res), 90)),
                          floor=float(np.sqrt(np.mean(np.square([r["sigma"] for r in rs])))),
                          rule_median=float(np.median([r["rule"] for r in rs]))))
    (ROOT / "prototypes/out/error_map.json").write_text(json.dumps(dict(rows=rows, cells=table), indent=1))
    print(f"{'iso':<9}{'v′':>7}{'v″':>7}{'n':>5}{'lines':>6}{'rms MHz':>10}{'median':>9}{'90%':>9}{'floor':>8}{'rule':>10}")
    for c in table:
        print(f"{c['iso']:<9}{c['v_upper']:>7}{c['v_lower']:>7}{c['n']:>5}{c['lines']:>6}{c['rms']:>10.2f}{c['median_abs']:>9.2f}"
              f"{c['p90_abs']:>9.2f}{c['floor']:>8.2f}{c['rule_median']:>10.1f}")
    # the rule's coverage: fraction of residuals within 1 and 2 rule-sigmas
    z = np.array([abs(r["residual"]) / r["rule"] for r in rows])
    print(f"\nrule coverage over {len(rows)} rows: within 1σ {np.mean(z < 1):.2f}, within 2σ {np.mean(z < 2):.2f}, "
          f"beyond 5σ {np.mean(z > 5):.3f}  (Gaussian: 0.68, 0.95, 0.000)")


if __name__ == "__main__":
    main()
