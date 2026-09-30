"""Calibration of each part of the Orsay 11 000-14 000 atlas against the model where its levels are known.

Parts II-IV (500 C) sit at v'' <= 12, where the model's X and B levels are fixed to ~1 MHz by comb data, so
each part's wavenumber scale can be fitted on its own. Part I was calibrated inside the level fit
(prototypes/orsay_fit.py, +24 ppb). The overlaps check the result: part II sits +2.0 mk above part I at
12 908-13 009 cm-1 (orsay-atlas-11000-14000.md), which the two scales must reproduce.

Assignment as in orsay_assign.py (model list thinned to 1.5x each part's own density, 120 MHz window);
only lines whose X and B levels both lie inside their measured J ranges are used.
"""
import csv, sys
import numpy as np
sys.path.insert(0, "prototypes")
import orsay_assign as A
from i2spec.model import RovibronicModel
from i2spec.intensity import line_list, SHARED_GRID, X_RMIN
from i2spec.level_corrections import load_level_corrections, J_MARGIN

MHZ = 29979.2458
PARTS = {1: 1063.15, 2: 773.15, 3: 773.15, 4: 773.15}


def model_lines(T, lo, hi):
    g = dict(SHARED_GRID)
    n_inner = max(int(round((g["rmin"] - X_RMIN) / g["step"])), 0)
    gx = dict(g, rmin=g["rmin"] - n_inner * g["step"], nlev=48)
    m = RovibronicModel("127I2", "i2spec2026i", grids={"X": gx, "B": dict(g, nlev=70)}, solver="dvr")
    ll = line_list(m, T, lo, hi, S_min=1e-28)
    d = dict(nu=np.asarray(ll.nu), S=np.asarray(ll.S), vu=np.asarray(ll.v_upper), vl=np.asarray(ll.v_lower),
             J=np.asarray(ll.J_lower), br=np.asarray(ll.branch))
    o = np.argsort(d["nu"])
    return {k: v[o] for k, v in d.items()}


def main():
    C = load_level_corrections("level_corrections_2026h")

    def covered(st, v, J):
        if (st == "X" and v <= 10) or (st == "B" and v == 0):
            return True
        r = C.coverage[st].get(v)
        return r is not None and r[0] - J_MARGIN <= J <= r[1] + J_MARGIN

    out = {}
    for part in (2, 3, 4, 1):
        rows = list(csv.DictReader(open(f"data/atlas_lines/orsay1982_part{part}.csv")))
        obs = np.array([float(r["sigma_cm1"]) for r in rows])
        eps = np.array([float(r["eps_mk"]) if r["eps_mk"] else 9.0 for r in rows]) * 1e-3 * MHZ
        M = model_lines(PARTS[part], obs.min() - 2, obs.max() + 2)
        M = A.thin(M, obs, 1.5, W=30.0) if part == 1 else A.thin_range(M, obs, 1.5)
        a = A.assign(obs, M, 0.004)
        sel = [(k, ai) for k, ai in enumerate(a) if ai is not None and M["vl"][ai] <= 17
               and covered("X", int(M["vl"][ai]), int(M["J"][ai]))
               and covered("B", int(M["vu"][ai]), int(M["J"][ai]) + (1 if M["br"][ai] > 0 else -1))]
        k = np.array([s[0] for s in sel]); ai = np.array([s[1] for s in sel])
        r = (obs[k] - M["nu"][ai]) * MHZ                   # MHz
        x = obs[k] * MHZ * 1e-9                            # MHz per ppb
        w = 1 / np.hypot(eps[k], 20.0) ** 2
        good = np.ones(len(r), bool)
        for _ in range(6):
            ppb = np.sum(w[good] * x[good] * r[good]) / np.sum(w[good] * x[good] ** 2)
            res = r - ppb * x
            s = 1.4826 * np.median(np.abs(res[good]))
            good = np.abs(res) < 3 * s
        err = 1 / np.sqrt(np.sum(w[good] * x[good] ** 2)) * np.sqrt(np.sum(w[good] * res[good] ** 2) / (good.sum() - 1))
        out[part] = (ppb, err)
        print(f"part {part} ({PARTS[part] - 273.15:.0f} C): {len(obs)} lines, {len(a) - a.count(None)} assigned, "
              f"{good.sum()} anchors -> scale {ppb:+.1f} +- {err:.1f} ppb; residual MAD {s:.0f} MHz", flush=True)
    d12 = (out[2][0] - out[1][0]) * 1e-9 * 12958 * 1e3
    print(f"\npart II - part I at 12 958 cm-1 from the scales: {d12:+.2f} mk (overlap measured +2.0 +- 0.3 mk)")


if __name__ == "__main__":
    main()
