"""Compare i2spec transmission with the Salami & Ross (2005) FTS atlas in one window.

Usage:  uv run python prototypes/compare_salami_ross.py NU_MIN NU_MAX [T_CELL_K] [--fix-offset]

The atlas records neither the I₂ pressure nor an absolute zero of transmission. Each fit therefore
adjusts:
  - the column density,
  - a linear baseline,
  - an additive zero offset (unless --fix-offset, which is needed where every line is weak,
    because then offset, baseline and column density are degenerate),
  - a wavenumber shift,
  - a Gaussian instrument width.
The fit runs with and without hyperfine structure. Outputs go to prototypes/out/: the
observed and model transmission with residuals, the strongest lines, and a line-list cache.
"""

import sys
import time
from pathlib import Path

import numpy as np
from scipy.optimize import brentq, least_squares

from i2spec.constants import MHZ_PER_CM
from i2spec.intensity import LineList, intensity_model, line_list
from i2spec.spectrum import (convolve_gaussian, cross_section, hyperfine_patterns, ils_kernel, number_density,
                             vapor_pressure)

ATLAS = Path("data/external/salami_ross_2005/salami_ross_2005_mmc1.txt")
CELL_LENGTH = 50.0  # cm (Salami & Ross 2005)
FIELDS = ("nu", "S", "v_upper", "v_lower", "J_lower", "branch", "E_lower")


ILS_START = {"gauss": 0.015, "sinc": 0.012, "sinc2": 0.008}  # width parameter (see i2spec.spectrum.ils_kernel)


def convolve(grid, y, kind, width):
    if kind == "gauss":
        return convolve_gaussian(grid, y, width)
    k = ils_kernel(grid[1] - grid[0], kind, width)
    return np.convolve(np.pad(y, k.size // 2, mode="edge"), k, mode="valid")


def fit(grid, sigma, nu_obs, t_obs, N0, fix_offset, ils="gauss"):
    x = nu_obs - nu_obs.mean()

    def model(p):
        log_n, b0, b1, c0, shift, width = p
        tr = convolve(grid, np.exp(-sigma * 10**log_n), ils, max(abs(width), 1e-4))
        return (b0 + b1 * x) * np.interp(nu_obs + shift, grid, tr) + (0.0 if fix_offset else c0)

    res = least_squares(lambda p: model(p) - t_obs, x0=[np.log10(N0), t_obs.max(), 0.0, 0.0, 0.0, ILS_START[ils]],
                        x_scale=[0.1, 0.1, 0.01, 0.02, 0.002, 0.005])
    p = res.x.copy()
    if fix_offset:
        p[3] = 0.0
    return p, model(p)


def cached_line_list(model, T, lo, hi, out):
    path = out / f"linelist_{lo:.0f}_{hi:.0f}_{T:.0f}K.npz"
    if path.exists():
        z = np.load(path)
        return LineList(**{k: z[k] for k in FIELDS}, T=T)
    lines = line_list(model, T, lo - 0.5, hi + 0.5, S_min=1e-26)
    np.savez(path, **{k: getattr(lines, k) for k in FIELDS})
    return lines


def main(lo, hi, T, fix_offset=False, ils="gauss"):
    data = np.loadtxt(ATLAS, skiprows=5)
    sel = (data[:, 0] >= lo) & (data[:, 0] <= hi)
    nu_obs, t_obs = data[sel, 0], data[sel, 1] / 100
    out = Path("prototypes/out")
    out.mkdir(exist_ok=True)

    t0 = time.time()
    model = intensity_model("127I2")
    lines = cached_line_list(model, T, lo, hi, out)
    t1 = time.time()
    grid = np.arange(lo - 0.3, hi + 0.3, 0.0005)
    sigmas = {"no hyperfine": cross_section(grid, lines, T),
              "hyperfine": cross_section(grid, lines, T, hyperfine=hyperfine_patterns(model, lines, dJ=0))}
    print(f"{len(lines)} lines with S >= 1e-26 cm; line list {t1 - t0:.0f} s, cross sections {time.time() - t1:.0f} s")

    N0 = number_density(293.15, T) * CELL_LENGTH
    stem = out / f"salami_ross_{lo:.0f}_{hi:.0f}"
    for name, sigma in sigmas.items():
        (log_n, b0, b1, c0, shift, fwhm), t_model = fit(grid, sigma, nu_obs, t_obs, N0, fix_offset, ils)
        resid = t_obs - t_model
        p_cell = 10**log_n / CELL_LENGTH * 1e6 * 1.380649e-23 * T  # Pa
        t_cold = brentq(lambda x: vapor_pressure(x) - p_cell, 150, 500)
        depth = 1 - t_model / t_model.max()
        bins = ((0, .02), (.02, .2), (.2, .5), (.5, 1))
        bias = ", ".join(f"{100 * resid[m].mean():+.1f}%" if m.any() else "n/a"
                         for m in ((depth >= a) & (depth < b) for a, b in bins))
        print(f"{name:13s}: rms residual {100 * resid.std():.2f}% of transmission (max {100 * np.abs(resid).max():.1f}%; "
              f"mean at depth <2/2-20/20-50/>50%: {bias}); max depth {100 * depth.max():.0f}%")
        print(f"{'':15s}N = {10**log_n:.3e} cm^-2 (I2 {p_cell:.1f} Pa, cold point {t_cold - 273.15:.1f} °C); "
              f"shift {shift * MHZ_PER_CM:+.1f} MHz; instrument {ils} width {abs(fwhm):.4f} cm^-1; baseline {b0:.3f}; "
              f"zero offset {'fixed at 0' if fix_offset else f'{c0:+.3f}'}")
        tag = name.replace(" ", "_") + ("" if ils == "gauss" else f"_{ils}")
        np.savetxt(f"{stem}_{tag}.txt", np.column_stack([nu_obs, t_obs, t_model, resid]),
                   header="nu_cm-1 t_obs t_model residual", fmt="%.6f")

    in_window = np.nonzero((lines.nu >= lo) & (lines.nu <= hi))[0]
    strongest = sorted(in_window[np.argsort(lines.S[in_window])[::-1][:10]], key=lambda k: lines.nu[k])
    with open(f"{stem}_lines.txt", "w") as f:
        f.writelines(f"{lines.nu[k]:.4f}\t{lines.label(k)}\n" for k in strongest)
    print("strongest lines:", ", ".join(f"{lines.label(k)} {lines.nu[k]:.3f}" for k in strongest))


if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    ils = next((a.split("=", 1)[1] for a in sys.argv if a.startswith("--ils=")), "gauss")
    main(float(args[0]), float(args[1]), float(args[2]) if len(args) > 2 else 293.15, "--fix-offset" in sys.argv, ils)
