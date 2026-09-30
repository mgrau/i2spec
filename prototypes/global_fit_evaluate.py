"""Score a global-fit result against every observation, next to the published model and IodineSpec5.

usage: global_fit_evaluate.py [result.json] [--nx=12] [--nb=16]
The result is prototypes/out/global_fit<tag>.json (from prototypes/global_fit.py). Atlas rows are
scored after removing each atlas's fitted calibration offset, since that offset is a property of the
atlas and not of the model.
"""
from __future__ import annotations

import csv
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from mlr_from_rkr import mlr                                                   # noqa: E402

from i2spec.constants import MHZ_PER_CM                                        # noqa: E402
from i2spec.observations import Line, Predictor, load_all                      # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
GRID_B = {"B": dict(rmin=2.35, rmax=8.0, h=0.01, order=10, nlev=85)}


def load(path, nx=None, nb=None):
    """The fitted pair. The long-range corrections are part of the answer: dropped, the X well is
    9 cm-1 too deep and every line moves by tens of GHz."""
    d = json.loads(Path(path).read_text())
    nx = d.get("nx", nx); nb = d.get("nb", nb)
    lr = d.get("lr") or None
    return mlr("X", d["X"]["x"], nx, lr=lr), mlr("B", d["B"]["x"], nb, lr=lr), d.get("offsets", [0.0, 0.0])


def main(argv):
    opts = {a.split("=")[0]: a.split("=", 1)[1] for a in argv if "=" in a}
    path = next((a for a in argv if not a.startswith("--")), ROOT / "prototypes/out/global_fit_v1.json")
    px, pb, offsets = load(path, int(opts.get("--nx", 12)), int(opts.get("--nb", 16)))
    ours = Predictor("hannover2008", potentials={"X": px, "B": pb}, grids=GRID_B)
    published = Predictor("hannover2008")
    default = Predictor()
    rms = lambda v: float(np.sqrt(np.mean(np.square(v)))) if len(v) else float("nan")   # noqa: E731

    print("atlas lines (each atlas's fitted calibration offset removed), MHz:")
    print(f"  {'atlas':18s} {'n':>6s} {'published':>10s} {'default':>10s} {'ours':>10s}")
    for k, name in enumerate(("salami_ross_2005", "apo_nist_2009")):
        rows = list(csv.DictReader(open(ROOT / "data/atlas_lines" / f"{name}.csv")))
        r = {"published": [], "default": [], "ours": []}
        for row in rows[::3]:
            line = Line.parse(row["line"])
            nu = float(row["value"]) * MHZ_PER_CM
            for tag, p in (("published", published), ("default", default), ("ours", ours)):
                try:
                    r[tag].append(nu - p.position(line))
                except Exception:
                    r[tag].append(np.nan)
        for tag in r:
            v = np.array(r[tag]); v = v[np.isfinite(v)]
            r[tag] = v - np.median(v)
        print(f"  {name:18s} {len(rows[::3]):6d} " + " ".join(f"{rms(r[t]):10.1f}" for t in ("published", "default", "ours")))

    print("\nprecision sets (absolute, hyperfine offsets from the default model), MHz:")
    print(f"  {'set':16s} {'n':>5s} {'published':>10s} {'default':>10s} {'ours':>10s}")
    tot = defaultdict(list)
    for ds in load_all():
        scale = MHZ_PER_CM if ds.unit == "cm-1" else 1.0
        r = defaultdict(list)
        for o in ds.observations:
            ref = o.ref_line or o.line
            if o.line.isotopologue != "127I2" or ref.isotopologue != "127I2":
                continue
            if o.kind == "interval" and ref == o.line:
                continue
            for tag, p in (("published", published), ("default", default), ("ours", ours)):
                try:
                    model = p(o, ds.unit) * scale
                except Exception:
                    continue
                r[tag].append(o.value * scale - model)
                tot[tag].append(o.value * scale - model)
        if r["ours"]:
            print(f"  {ds.id:16s} {len(r['ours']):5d} " + " ".join(f"{rms(r[t]):10.2f}" for t in ("published", "default", "ours")))
    print(f"  {'all':16s} {len(tot['ours']):5d} " + " ".join(f"{rms(tot[t]):10.2f}" for t in ("published", "default", "ours")))
    print(f"\ncalibration offsets fitted: " + ", ".join(f"{n} {o * MHZ_PER_CM:+.1f} MHz"
                                                        for n, o in zip(("salami_ross", "apo"), offsets)))


if __name__ == "__main__":
    main(sys.argv[1:])
