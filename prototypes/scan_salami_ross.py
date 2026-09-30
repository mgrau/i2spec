"""Scan the Salami & Ross (2005) atlas with narrow-window fits across 14 300-20 000 cm⁻¹.

Usage:  uv run python prototypes/scan_salami_ross.py [STEP WIDTH]          full scan (default 50 4 cm⁻¹)
        uv run python prototypes/scan_salami_ross.py --one CENTRE [T_K]    one window, printed

Each window is fitted (compare_salami_ross.fit) at the cell temperature of every atlas segment that
covers it: column density, linear baseline, zero offset (fixed at 0 when the deepest line is
< 30 %), wavenumber shift and Gaussian instrument width. The model is the master line list plus
the continuum; lines with S >= S_HYPERFINE get hyperfine structure (ΔJ = 0 blocks).

Within a segment the I₂ pressure was fixed by the sidearm, so the fitted column density should not
drift with wavenumber: its variation measures errors in the relative band intensities (μ(R),
Franck-Condon factors, populations). The shift maps line-position errors against the atlas
calibration (±0.003 cm⁻¹). Output: prototypes/out/scan_salami_ross.txt.
"""

import sys
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

sys.path.insert(0, str(Path(__file__).resolve().parent))
from compare_salami_ross import ATLAS, CELL_LENGTH, fit  # noqa: E402

from i2spec.constants import MHZ_PER_CM  # noqa: E402
from i2spec.continuum import continuum_model  # noqa: E402
from i2spec.intensity import intensity_model, master_line_list  # noqa: E402
from i2spec.spectrum import cross_section, hyperfine_patterns, number_density, vapor_pressure  # noqa: E402

#: Atlas segments (Salami & Ross 2005, Table 1): range (cm⁻¹) and cell-body temperature. The first
#: range is extended down to the start of the data.
SEGMENTS = ((19100.0, 20500.0, 293.15), (16500.0, 19500.0, 293.15), (15400.0, 17400.0, 323.15),
            (14250.0, 16500.0, 463.15))
S_HYPERFINE = 1e-22  # cm: weaker lines change σN by < 0.5 % through their hyperfine structure
OUT = Path(__file__).resolve().parent / "out" / "scan_salami_ross.txt"
HEADER = "centre_cm-1 T_K lines depth rms noise log10N cold_point_C shift_MHz width_cm-1 b0 b1 offset"
_state = {}


def _init():
    model = intensity_model("127I2")
    master = master_line_list(model, 11000.0, 20100.0, T_range=(200.0, 600.0), S_min=1e-27)
    _state.update(data=np.loadtxt(ATLAS, skiprows=5), model=model, master=master,
                  continuum=continuum_model("127I2", master.upper_cut, T_max=600.0))


def fit_window(job):
    centre, width, T = job
    s = _state
    lo, hi = centre - width / 2, centre + width / 2
    sel = (s["data"][:, 0] >= lo) & (s["data"][:, 0] <= hi)
    nu_obs, t_obs = s["data"][sel, 0], s["data"][sel, 1] / 100
    lines = s["master"].at(T, lo - 0.5, hi + 0.5, S_min=1e-27)
    pattern, single = hyperfine_patterns(s["model"], lines, dJ=0), (np.zeros(1), np.ones(1))
    grid = np.arange(lo - 0.3, hi + 0.3, 0.0005)
    sigma = (cross_section(grid, lines, T, hyperfine=lambda k: pattern(k) if lines.S[k] >= S_HYPERFINE else single)
             + s["continuum"].cross_section(grid, T))
    fix_offset = t_obs.min() > 0.7 * t_obs.max()
    (log_n, b0, b1, c0, shift, w), t_model = fit(grid, sigma, nu_obs, t_obs, number_density(295.0, T) * CELL_LENGTH,
                                                 fix_offset)
    p_cell = 10**log_n / CELL_LENGTH * 1e6 * 1.380649e-23 * T  # Pa
    try:
        cold = brentq(lambda x: vapor_pressure(x) - p_cell, 150, 600) - 273.15
    except ValueError:
        cold = np.nan
    noise = np.std(np.diff(t_obs, 2)) / np.sqrt(6)
    return (centre, T, len(lines), 1 - t_obs.min() / t_obs.max(), np.std(t_obs - t_model), noise, log_n, cold,
            shift * MHZ_PER_CM, abs(w), b0, b1, 0.0 if fix_offset else c0)


def jobs(step, width):
    for centre in np.arange(14300.0, 20000.0 + 1e-9, step):
        for T in sorted({T for a, b, T in SEGMENTS if a <= centre - width / 2 and centre + width / 2 <= b}):
            yield float(centre), width, T


def main(step=50.0, width=4.0):
    todo = list(jobs(step, width))
    OUT.parent.mkdir(exist_ok=True)
    t0, rows = time.time(), []
    with Pool(initializer=_init) as pool, open(OUT, "w") as f:
        f.write(f"# {HEADER}\n")
        for n, row in enumerate(pool.imap_unordered(fit_window, todo), start=1):
            rows.append(row)
            f.write(" ".join(f"{x:.6g}" for x in row) + "\n")
            f.flush()
            print(f"{n}/{len(todo)} {time.time() - t0:.0f} s: {row[0]:.0f} cm-1 at {row[1]:.0f} K, "
                  f"rms {100 * row[4]:.2f}% (noise {100 * row[5]:.2f}%), cold point {row[7]:.1f} °C, "
                  f"shift {row[8]:+.0f} MHz", flush=True)
    rows.sort()
    np.savetxt(OUT, rows, fmt="%.6g", header=HEADER)


if __name__ == "__main__":
    if "--one" in sys.argv:
        i = sys.argv.index("--one")
        centre = float(sys.argv[i + 1])
        temps = ([float(sys.argv[i + 2])] if len(sys.argv) > i + 2 else
                 sorted({T for a, b, T in SEGMENTS if a <= centre - 2 and centre + 2 <= b}))
        t = time.time()
        _init()
        print(f"setup {time.time() - t:.1f} s")
        for T in temps:
            t = time.time()
            row = fit_window((centre, 4.0, T))
            print(f"{time.time() - t:.1f} s: " + ", ".join(f"{k}={v:.5g}" for k, v in zip(HEADER.split(), row)))
    else:
        main(*map(float, sys.argv[1:3]))
