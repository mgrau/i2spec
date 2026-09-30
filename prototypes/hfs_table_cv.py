"""Does a hyperfine table predict lines it has not seen? Leave-one-line-out against the formulae.

Usage:  uv run python prototypes/hfs_table_cv.py [--new=b_state_lines_2026n] [--old=b_state_lines]
Output: prototypes/out/hfs_table_cv.json, and a per-set summary on stdout

Every 127I2 intra-line splitting of the data in use (observations.load_all), v' <= 70. For each line the
splittings are predicted three ways: by the published formulae alone (no table), and by the old and the
new table each rebuilt without that line (HyperfineTable(rows, exclude=...)), so a line never predicts
itself. Lines that are in neither table are predicted by the full tables, which is what the model does
for them. The score is z = (model - observed) / sigma, sigma the stated uncertainty; kHz rms as well.
"""
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

from i2spec.hfs_table import DATA, HyperfineTable
from i2spec.hyperfine_fit import HyperfineFit
from i2spec.observations import load_all

ROOT = Path(__file__).resolve().parents[1]


def main(argv):
    opts = {a.split("=")[0]: a.split("=", 1)[1] for a in argv if "=" in a}
    tables = {"old": opts.get("--old", "b_state_lines"), "new": opts.get("--new", "b_state_lines_2026n")}
    rows = {k: json.loads((DATA / f"{v}.json").read_text())["rows"] for k, v in tables.items()}
    full = {k: HyperfineTable(r) for k, r in rows.items()}
    in_table = {k: {r["line"] for r in r_} for k, r_ in rows.items()}

    h = HyperfineFit(datasets=load_all(), max_v_upper=70, model_floor=0.0)
    model = h._model("127I2")
    data = [d for d in h.data if d.isotopologue == "127I2"]
    by_line = defaultdict(list)
    for d in data:
        by_line[(d.branch, d.J_lower, d.v_upper, d.v_lower)].append(d)

    def predict(key, table):
        br, J, vu, vl = key
        _, comps = model.hyperfine_components(vu, vl, J, br, table=table)
        main = sorted([q for q in comps if q.label], key=lambda q: q.offset)
        return np.array([q.offset for q in main])

    out = []
    for n, (key, ds) in enumerate(sorted(by_line.items())):
        br, J, vu, vl = key
        name = f"{br}({J}) {vu}-{vl}"
        preds = {"formula": predict(key, None)}
        for k in ("old", "new"):
            t = HyperfineTable(rows[k], exclude=(name,)) if name in in_table[k] else full[k]
            preds[k] = predict(key, t)
        for d in ds:
            rec = dict(set=d.dataset, line=name, v=vu, J=J, label=d.key, sigma_kHz=d.uncertainty * 1e3,
                       in_old=name in in_table["old"], in_new=name in in_table["new"])
            for k, a in preds.items():
                ok = max(d.rank, d.ref_rank) <= len(a)
                rec[k] = float(((a[d.rank - 1] - a[d.ref_rank - 1]) - d.value) * 1e3) if ok else None
            out.append(rec)
        if n % 10 == 0:
            print(f"  {n}/{len(by_line)} lines", flush=True)
    (ROOT / "prototypes/out").mkdir(exist_ok=True)
    (ROOT / "prototypes/out/hfs_table_cv.json").write_text(json.dumps(dict(tables=tables, rows=out), indent=1))

    def summary(sel, label):
        z = {k: np.array([r[k] / r["sigma_kHz"] for r in sel if r[k] is not None]) for k in ("formula", "old", "new")}
        e = {k: np.array([r[k] for r in sel if r[k] is not None]) for k in ("formula", "old", "new")}
        if not len(z["new"]):
            return
        cells = " | ".join(f"{k} {np.mean(np.abs(z[k]) <= 1) * 100:4.0f}% {np.mean(np.abs(z[k]) <= 3) * 100:4.0f}% "
                           f"{np.sqrt(np.mean(e[k] ** 2)):8.1f}" for k in ("formula", "old", "new"))
        print(f"{label:<22} {len(sel):5d}  {cells}")

    print(f"\n{'set':<22} {'n':>5}  per model: within 1 sigma, within 3 sigma, rms kHz (each line held out)")
    summary(out, "all")
    summary([r for r in out if r["sigma_kHz"] <= 25], "sigma <= 25 kHz")
    for s in sorted({r["set"] for r in out}):
        summary([r for r in out if r["set"] == s], s)


if __name__ == "__main__":
    main(sys.argv[1:])
