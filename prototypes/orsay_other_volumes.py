"""What v'' and v' would the other Orsay volumes measure? Model lines thinned per 15 cm-1 window.

Usage: python prototypes/orsay_other_volumes.py LO HI T_KELVIN B_NLEV
"""
import numpy as np, sys
from collections import Counter
from i2spec.model import RovibronicModel
from i2spec.intensity import line_list, SHARED_GRID, X_RMIN

g = dict(SHARED_GRID)
n_inner = max(int(round((g["rmin"] - X_RMIN) / g["step"])), 0)
gx = dict(g, rmin=g["rmin"] - n_inner * g["step"], nlev=48)
m = RovibronicModel("127I2", "i2spec2026d", grids={"X": gx, "B": dict(g, nlev=int(sys.argv[4]))}, solver="dvr")
LO, HI, T = float(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3])
W, PER = 15.0, 36
ll = line_list(m, T, LO, HI, S_min=1e-30)
nu, S, vl, vu = (np.asarray(x) for x in (ll.nu, ll.S, ll.v_lower, ll.v_upper))
print(f"{LO:.0f}-{HI:.0f} at {T:.0f} K: {len(nu)} model lines, B nlev={sys.argv[4]}")
kv, ku = [], []
for lo in np.arange(LO, HI, W):
    sel = np.where((nu >= lo) & (nu < lo + W))[0]
    k = sel[np.argsort(S[sel])[::-1][:PER]]
    if len(k) == 0:
        print(f'  {lo:7.0f}: no model lines'); continue
    kv += list(vl[k]); ku += list(vu[k])
    if True:
        print(f"  {lo:7.0f}: v'' {vl[k].min()}-{vl[k].max()} (mode {Counter(vl[k].tolist()).most_common(1)[0][0]})"
              f"   v' {vu[k].min()}-{vu[k].max()}")
kv, ku = np.array(kv), np.array(ku)
print(f"  tabulated-density set: {len(kv)} lines; v'' > 17: {np.sum(kv>17)}; v' > 50: {np.sum(ku>50)};"
      f" v' range {ku.min()}-{ku.max()}")
