"""Compare the Nölleke 2018 raw absorption scans (915-985 nm) with modelled transmission for their cell.

The earlier attempt (prototypes/nolleke_assign.py) matched their peak list against a model at 300 K. But
absorption at 10 150-10 930 cm-1 needs lower levels 5 000+ cm-1 up, which are populated only because
their cell body was at 300 C (catalog broadband.yaml: 573 K, cold finger 39 C, 127 Pa, 225 cm path).
This uses the right temperature and the raw scans (data/external/nolleke_2018/scan_data.csv), so it
does not depend on their peak-picking threshold or on how they split hyperfine blends.

For each 10 cm-1 window the observed absorbance -ln T and the modelled sigma*N are both high-passed
(a running median of 0.5 cm-1 removed), and cross-correlated over shifts of +-0.3 cm-1. A model that has
the X levels right registers at a shift near zero with a correlation peak far above the off-peak spread.

usage: nolleke_transmission.py [--parameters=i2spec2026i,hannover2008] [--window=10]
"""
import sys
from pathlib import Path

import numpy as np
from scipy.ndimage import median_filter

from i2spec.intensity import intensity_model, line_list
from i2spec.model import RovibronicModel
from i2spec.spectrum import cross_section

ROOT = Path(__file__).resolve().parents[1]
SCAN = ROOT / "data/external/nolleke_2018/scan_data.csv"
T_CELL = 573.15
P_PA, PATH_CM = 127.0, 225.0
KB = 1.380649e-23
N_COL = P_PA / (KB * T_CELL) * 1e-6 * PATH_CM          # molecules / cm2
MHZ = 29979.2458
WATER_HALF = 0.05      # cm-1 blanked each side of a flagged water line
STRONG = 0.03          # absorbance: no iodine line here reaches this (the list's median is ~0.005)


def highpass(y, step, width_cm=0.5):
    k = max(3, int(width_cm / step) | 1)
    return y - median_filter(y, size=k, mode="nearest")


def main(argv):
    opts = {a.split("=")[0]: a.split("=", 1)[1] for a in argv if "=" in a}
    sets = opts.get("--parameters", "i2spec2026i,hannover2008").split(",")
    W = float(opts.get("--window", 10))
    d = np.genfromtxt(SCAN, delimiter=";", skip_header=1)
    nu = 1e7 / d[:, 0]
    o = np.argsort(nu)
    nu, tr = nu[o], d[o, 1]
    step = float(np.median(np.diff(nu)))
    grid = np.arange(nu[0], nu[-1], step)
    A_obs = np.interp(grid, nu, -np.log(np.clip(tr, 1e-3, None)))
    A_obs = highpass(A_obs, step)
    # the open beam path crosses the 940 nm water band; the list flags the water lines (column 3): blank
    # +-WATER_HALF around each, and anything stronger than the strongest iodine line could be, in both spectra
    L = np.genfromtxt(ROOT / "data/external/nolleke_2018/iodine_atlas.csv", delimiter=",", skip_header=1)
    water = 1e7 / L[np.isfinite(L[:, 2]), 2]
    mask = np.zeros(grid.size, bool)
    for w in water:
        mask[(grid > w - WATER_HALF) & (grid < w + WATER_HALF)] = True
    mask |= A_obs > STRONG
    print(f"masked {mask.mean()*100:.1f} % of the scan: {len(water)} flagged water lines and features above {STRONG}")
    print(f"scan {grid[0]:.1f}-{grid[-1]:.1f} cm-1, step {step*MHZ:.0f} MHz; column {N_COL:.2e} cm-2 at {T_CELL:.0f} K")

    shifts = np.arange(-0.3, 0.3 + step / 2, step)
    ks = np.round(shifts / step).astype(int)
    results = {}
    for ps in sets:
        m = intensity_model("127I2")
        if ps != m.parameters:
            g = m.grids
            m = RovibronicModel("127I2", ps, grids=g, solver="dvr")
        ll = line_list(m, T_CELL, grid[0] - 1, grid[-1] + 1, S_min=1e-27)
        A_mod = highpass(cross_section(grid, ll, T_CELL) * N_COL, step)
        A_mod[mask] = 0.0
        A_o = np.where(mask, 0.0, A_obs)
        vl = np.asarray(ll.v_lower); S = np.asarray(ll.S); lnu = np.asarray(ll.nu)
        rows = []
        for lo in np.arange(grid[0], grid[-1] - W / 2, W):
            sel = (grid >= lo) & (grid < lo + W)
            a, b = A_o[sel], A_mod[sel]
            if b.std() == 0 or a.std() == 0:
                continue
            a = (a - a.mean()) / a.std()
            b = (b - b.mean()) / b.std()
            cc = np.array([np.mean(a[max(0, k):len(a) + min(0, k)] * b[max(0, -k):len(b) - max(0, k)]) for k in ks])
            i = int(np.argmax(cc))
            off = np.abs(shifts - shifts[i]) > 0.03
            z = (cc[i] - cc[off].mean()) / cc[off].std()
            inw = (lnu >= lo) & (lnu < lo + W)
            top = vl[inw][np.argsort(S[inw])[::-1][:40]] if inw.any() else np.array([0])
            rows.append((lo, shifts[i], cc[i], z, int(np.bincount(top).argmax()), top.min(), top.max()))
        results[ps] = np.array(rows)
        r = results[ps]
        good = r[:, 3] > 5
        print(f"\n{ps}: {len(r)} windows, {good.sum()} register at z > 5; "
              f"median |shift| of those {np.median(np.abs(r[good, 1]))*MHZ if good.any() else float('nan'):.0f} MHz")
        print(f"   {'window':>8} {'shift MHz':>10} {'peak r':>7} {'z':>6}  v'' (strong lines)")
        for row in r[::4]:
            print(f"   {row[0]:8.0f} {row[1]*MHZ:+10.0f} {row[2]:7.3f} {row[3]:6.1f}  {row[4]:.0f} ({row[5]:.0f}-{row[6]:.0f})")
    np.save(ROOT / "prototypes/out/nolleke_transmission.npy", results, allow_pickle=True)


if __name__ == "__main__":
    main(sys.argv[1:])
