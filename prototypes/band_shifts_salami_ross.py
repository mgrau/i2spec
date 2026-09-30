"""Per-band line-position offsets between i2spec and the Salami & Ross atlas (20 °C segments).

Usage:  uv run python prototypes/band_shifts_salami_ross.py

For each band (v', v'') the model lines of a J'' range are shifted rigidly against the atlas: the
score of a shift s is Σ S_k · depth(ν_k + s), with depth = (local maximum − transmission) of the
atlas. The best s (atlas − model) and its contrast (peak score over the median score) are printed.
With the X levels known far better than the B levels, s(v', v'') ≈ ΔE_B(v', J') for all v''.
Output: prototypes/out/band_shifts_salami_ross.txt.
"""

from pathlib import Path

import numpy as np
from scipy.ndimage import maximum_filter1d

from i2spec.constants import MHZ_PER_CM
from i2spec.intensity import intensity_model, master_line_list

ATLAS = Path("data/external/salami_ross_2005/salami_ross_2005_mmc1.txt")
T = 293.15
RANGE = (16600.0, 20050.0)  # 20 °C segments
SHIFTS = np.arange(-3.0, 3.0 + 1e-9, 0.001)  # cm⁻¹
J_RANGES = ((0, 50), (51, 100))


def main():
    data = np.loadtxt(ATLAS, skiprows=5)
    data = data[(data[:, 0] >= RANGE[0]) & (data[:, 0] <= RANGE[1])]
    nu_a, t = data[:, 0], data[:, 1]
    depth = maximum_filter1d(t, size=801) - t  # 4 cm⁻¹ running maximum as the local baseline
    master = master_line_list(intensity_model("127I2"), 11000.0, 20100.0, T_range=(200.0, 600.0), S_min=1e-27)
    lines = master.at(T, *RANGE, S_min=1e-24)
    rows = []
    for v_lower in (0, 1, 2):
        for v_upper in range(20, 63):
            out = [v_upper, v_lower]
            for j0, j1 in J_RANGES:
                m = (lines.v_upper == v_upper) & (lines.v_lower == v_lower) & (lines.J_lower >= j0) & (lines.J_lower <= j1)
                if m.sum() < 10:
                    out += [np.nan, np.nan, 0]
                    continue
                nu, S = lines.nu[m], lines.S[m]
                score = np.array([np.dot(S, np.interp(nu + s, nu_a, depth, left=0, right=0)) for s in SHIFTS])
                k = int(np.argmax(score))
                out += [SHIFTS[k], score[k] / np.median(score), int(m.sum())]
            rows.append(out)
            print(f"v'={v_upper:2d} v''={v_lower}: " + "   ".join(
                f"J''{j0}-{j1}: shift {out[2 + 3 * i]:+.4f} cm-1 ({out[2 + 3 * i] * MHZ_PER_CM:+8.0f} MHz), "
                f"contrast {out[3 + 3 * i]:4.1f}, {out[4 + 3 * i]:3d} lines" for i, (j0, j1) in enumerate(J_RANGES)))
    np.savetxt(Path("prototypes/out/band_shifts_salami_ross.txt"), np.array(rows, dtype=float), fmt="%.5g",
               header="v_upper v_lower shift_J0-50_cm-1 contrast n shift_J51-100_cm-1 contrast n")


if __name__ == "__main__":
    main()
