"""Every observation against a candidate pair of MLR potentials, next to the published model.

Usage:  uv run python prototypes/mlr_evaluate.py <b_json> [x_json]      (about 8 minutes)
        b_json: a B fit from prototypes/mlr_b_levels.py (prototypes/out/b_levels_<n>_<rf>.json) or a
                data/potentials/mlr_b_*.json file; x_json defaults to data/potentials/mlr_x_2026b.json.
Output: prototypes/out/mlr_evaluate.json
"""

import json
import sys
from pathlib import Path

import numpy as np

from i2spec.constants import MHZ_PER_CM
from i2spec.observations import Predictor, load_all, residuals
from i2spec.potentials import MLRPotential, load_potentials

ROOT = Path(__file__).resolve().parents[1]
X0, B0 = load_potentials()["X"], load_potentials()["B"]
ASYM = X0.De + 7602.9762
GRID_B = {"B": dict(rmin=2.35, rmax=8.0, h=0.01, order=10, nlev=70)}


def load_b(path):
    d = json.loads(Path(path).read_text())
    if "beta" in d:
        return MLRPotential(De=d["De"], Re=d["Re"], C={int(k): v for k, v in d["C"].items()}, beta=tuple(d["beta"]),
                            p=d["p"], q=d["q"], Rref=d["Rref_factor"] * d["Re"], Te=d["Te"], bo=B0)
    n, rf, x = d["n"], d["rf"], d["x"]
    Te, Re, C5, C6, C8, C10 = x[n:]
    return MLRPotential(De=ASYM - Te, Re=Re, C={5: C5, 6: C6, 8: C8, 10: C10}, beta=tuple(x[:n]), p=6, q=4,
                        Rref=rf * Re, Te=Te, bo=B0)


def load_x(path):
    d = json.loads(Path(path).read_text())
    return MLRPotential(De=d["De"], Re=d["Re"], C={int(k): v for k, v in d["C"].items()}, beta=tuple(d["beta"]),
                        p=d["p"], q=d["q"], Rref=d["Rref"], bo=X0)


def main():
    b = load_b(sys.argv[1])
    x = load_x(sys.argv[2] if len(sys.argv) > 2 else ROOT / "data/potentials/mlr_x_2026b.json")
    preds = {"published": Predictor(grids=GRID_B), "MLR X+B": Predictor(potentials={"X": x, "B": b}, grids=GRID_B)}
    rms = lambda v: float(np.sqrt(np.mean(np.square(v))))          # noqa: E731
    out, tot = [], {k: [] for k in preds}
    print(f"{'set':<16}{'rows':>5}{'v′':>8}{'published':>12}{'MLR X+B':>10}")
    for ds in load_all():
        scale = MHZ_PER_CM if ds.meta["unit"] == "cm-1" else 1.0
        obs = ds.observations
        keep = np.array([not (o.kind == "interval" and (o.ref_line or o.line) == o.line) for o in obs])
        if not keep.any():
            continue
        row = dict(id=ds.id, n=int(keep.sum()), v_upper=(min(o.line.v_upper for o in obs), max(o.line.v_upper for o in obs)))
        for k, p in preds.items():
            r = residuals(p, ds) * scale
            row[k] = rms(r[keep]); tot[k] += list(r[keep])
        out.append(row)
        print(f"{ds.id:<16}{row['n']:>5}{row['v_upper'][0]:>4}-{row['v_upper'][1]:<3}{row['published']:>12.2f}{row['MLR X+B']:>10.2f}")
    print(f"{'ALL':<16}{len(tot['published']):>5}{'':>8}{rms(tot['published']):>12.2f}{rms(tot['MLR X+B']):>10.2f}")
    (ROOT / "prototypes/out/mlr_evaluate.json").write_text(json.dumps(dict(b=str(sys.argv[1]), sets=out), indent=1))


if __name__ == "__main__":
    main()
