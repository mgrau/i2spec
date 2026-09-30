"""Does the measured-parameter table predict hyperfine splittings it has not seen? (lever 3)

Usage:  uv run python prototypes/hfs_table_validate.py        (about 5 minutes)
Output: prototypes/out/hfs_table_validate.json

For every precisely measured line, the table is rebuilt without that line and its splittings are
predicted from the remaining lines at the same v' (or, above v' = 53, from the neighbouring v'). That
is compared with the published formulae on the same splittings. A line whose v' has no other measured
line falls back to the formula, so it scores 1:1; the gain is on the rest.
"""

import json
from pathlib import Path

import numpy as np

from i2spec.hfs_table import HyperfineTable
from i2spec.hyperfine_fit import HyperfineFit, rms
from i2spec.observations import load_dataset

ROOT = Path(__file__).resolve().parents[1]
SETS = ("bipm2003a", "bipm2003b", "bipm2003c", "bipm2003d", "bipm2003e", "bipm2005a", "bipm2012a", "reinhardt2006a",
        "reinhardt2007a", "bodermann1998b", "yoshiki2023a", "matsunaga2024a", "kobayashi2016a", "nishiyama2024a")


def main():
    h = HyperfineFit(datasets=[load_dataset(ROOT / "data/observations" / s) for s in SETS], max_v_upper=70)
    full = HyperfineTable.load()
    results = []
    for key in h.lines:
        iso, br, Jl, vu, vl = key
        rows = [d for d in h.data if h.line_of(d) == key]
        if iso != "127I2" or len(rows) < 8 or np.median([d.uncertainty for d in rows]) > 0.025:
            continue
        label = f"{br}({Jl}) {vu}-{vl}"
        held_out = HyperfineTable.load(exclude=(label,))
        model = h._model(iso)
        devs = {}
        for name, table in (("formula", None), ("table (held out)", held_out)):     # never the default table
            _, comps = model.hyperfine_components(vu, vl, Jl, br, table=table)
            a = np.array(sorted(q.offset for q in comps if q.label))
            devs[name] = np.array([(a[d.rank - 1] - a[d.ref_rank - 1]) - d.value for d in rows if max(d.rank, d.ref_rank) <= len(a)]) * 1e3
        others = [r["J"] for r in full.rows if r["v"] == vu and r["line"] != label]
        results.append(dict(line=label, dataset=rows[0].dataset, v=vu, J=Jl + (1 if br == "R" else -1), n=len(rows),
                            sigma_kHz=float(np.median([d.uncertainty for d in rows])) * 1e3,
                            others_at_v=len(others), formula=rms(devs["formula"]), table=rms(devs["table (held out)"])))
        r = results[-1]
        print(f"  {r['dataset']:<15} {label:<13} v′={vu:2d} J′={r['J']:3d}  others at v′: {len(others):2d}  "
              f"formula {r['formula']:8.1f} kHz  table {r['table']:8.1f} kHz", flush=True)
    (ROOT / "prototypes/out/hfs_table_validate.json").write_text(json.dumps(results, indent=1))
    print("\nby data set (rms over the set's lines, kHz):")
    for s in SETS:
        rs = [r for r in results if r["dataset"] == s]
        if rs:
            f = np.sqrt(np.mean([r["formula"] ** 2 for r in rs])); t = np.sqrt(np.mean([r["table"] ** 2 for r in rs]))
            print(f"  {s:<15} {len(rs):2d} lines  formula {f:8.1f}  table {t:8.1f}  {'better' if t < f - 0.05 else ('same' if abs(t - f) <= 0.05 else 'WORSE')}")
    with_others = [r for r in results if r["others_at_v"] > 0]
    print(f"\nlines with another measured line at their v′: {len(with_others)}; "
          f"formula {np.sqrt(np.mean([r['formula']**2 for r in with_others])):.1f} kHz, "
          f"table {np.sqrt(np.mean([r['table']**2 for r in with_others])):.1f} kHz; "
          f"worse on {sum(r['table'] > r['formula'] for r in with_others)} of them")


if __name__ == "__main__":
    main()
