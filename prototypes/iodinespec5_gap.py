"""IodineSpec5's Dunham column at v'' > 17 (Martin 1986 X levels via Gerstenkorn-Luc B) vs our potentials,
and the BKT02 hyperfine parameters it prints vs ours, by (v, J)."""
import sys
from collections import defaultdict
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).parent))
from iodinespec5_compare import DEFAULT_DIR, parse_L, parse_H
from mlr_evaluate import load_x, load_b, GRID_B, ROOT
from i2spec.constants import MHZ_PER_CM
from i2spec.model import RovibronicModel
from i2spec import hfs_params

out_dir = DEFAULT_DIR
pub = RovibronicModel("127I2", "hannover2008", grids=GRID_B)
mlr = RovibronicModel("127I2", "hannover2008", grids=GRID_B,
                      potentials={"X": load_x(ROOT / "data/potentials/mlr_x_2026c.json"),
                                  "B": load_b(ROOT / "data/potentials/mlr_b_2026c.json")})
# --- the gap: Dun - ours, by v'' ---
byv = defaultdict(list)
for p in sorted(out_dir.glob("*_L.OUT")):
    for line, pot, nir, dun in parse_L(p):
        if line.v_lower <= 17 or line.v_upper > 43 or dun is None or line.v_lower > 54 or line.J_lower > 150:
            continue
        a = pub.transition(line.v_upper, line.v_lower, line.J_lower, line.branch) * MHZ_PER_CM
        b = mlr.transition(line.v_upper, line.v_lower, line.J_lower, line.branch) * MHZ_PER_CM
        byv[line.v_lower].append((line, dun - a, dun - b, pot - a if pot else None))
print("IodineSpec5 FREQofDun (Martin 1986 X + Gerstenkorn-Luc B) minus ours, MHz, lines with v'' > 17, v' <= 43")
print(f"  {'v″':>3} {'n':>4} {'Dun-hannover2008':>18} {'Dun-MLR 2026c':>16}   (mean, and spread)   v' range")
for v in sorted(byv):
    d1 = np.array([x[1] for x in byv[v]]); d2 = np.array([x[2] for x in byv[v]])
    vps = sorted({x[0].v_upper for x in byv[v]})
    print(f"  {v:3d} {len(d1):4d} {d1.mean():+12.1f} ±{d1.std():7.1f} {d2.mean():+10.1f} ±{d2.std():7.1f}    v'={vps[0]}-{vps[-1]}")
# the Gerstenkorn-Luc B part alone: v'' <= 17 lines where Dun and Pot both exist, Dun - Pot by v'
byvp = defaultdict(list)
for p in sorted(out_dir.glob("*_L.OUT")):
    for line, pot, nir, dun in parse_L(p):
        if pot is not None and dun is not None and line.v_lower <= 17:
            byvp[line.v_upper].append(dun - pot)
print("\nDun - Pot (Gerstenkorn-Luc Dunham vs the potential, both IodineSpec5) for v'' <= 17, by v':")
for v in sorted(byvp):
    d = np.array(byvp[v]); print(f"  v'={v:2d} n={len(d):4d} mean {d.mean():+9.1f} MHz  std {d.std():7.1f}")

# --- hyperfine parameters: IodineSpec5 - ours by v and J ---
seen = {}
for p in sorted(out_dir.glob("*_H.OUT")):
    for b in parse_H(p):
        if b["line"] is not None and b["line"].v_upper <= 43:
            seen[b["line"]] = b
res = defaultdict(list)
for line, b in seen.items():
    xp, bp = hfs_params.line_states("127I2", line.v_upper, line.v_lower, pub._reference_term("B", line.v_upper),
                                    pub._reference_term("X", line.v_lower))
    Jl = line.J_lower; Ju = Jl + (1 if line.branch == "R" else -1)
    x, bb = xp(Jl), bp(Ju)
    res["eqQ_X"].append((line.v_lower, Jl, b["params"]["eqQ_X"] - x.eqQ))
    res["eqQ_B"].append((line.v_upper, Ju, b["params"]["eqQ_B"] - bb.eqQ))
    res["C_B"].append((line.v_upper, Ju, b["params"]["Csr_B"] - bb.C))
    res["d_B"].append((line.v_upper, Ju, b["params"]["d_B"] - bb.d))
    res["delta_B"].append((line.v_upper, Ju, b["params"]["del_B"] - bb.delta))
for name in res:
    arr = np.array(res[name])
    v, J, d = arr[:, 0], arr[:, 1], arr[:, 2]
    y = J * (J + 1) / 1e4
    A = np.c_[np.ones_like(v), v, y]
    coef, *_ = np.linalg.lstsq(A, d, rcond=None)
    r = d - A @ coef
    print(f"\n{name}: IodineSpec5 - ours = {coef[0]:+.4f} {coef[1]:+.4f}*v {coef[2]:+.4f}*J(J+1)/1e4   residual rms {np.sqrt((r**2).mean()):.4f}  (raw rms {np.sqrt((d**2).mean()):.4f})")
    for lo, hi in ((0, 5), (6, 12), (13, 17)) if name.endswith("X") else ((0, 10), (11, 20), (21, 30), (31, 43)):
        m = (v >= lo) & (v <= hi)
        if m.any():
            print(f"   v={lo:2d}-{hi:2d} n={m.sum():4d} mean {d[m].mean():+8.4f}  rms {np.sqrt((d[m]**2).mean()):8.4f}")
