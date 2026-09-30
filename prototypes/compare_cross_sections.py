"""Absolute I₂ cross sections: i2spec vs. Spietz et al. (2006) and Saiz-Lopez et al. (2004), with no free parameters.

Usage:  uv run python prototypes/compare_cross_sections.py DATASET      (spietz025 | spietz059 | saizlopez)

Model: B-X lines from the master line list at temperature T, with Voigt profiles for ~1 bar of
pressure broadening, plus the continuum (A←X, C←X, and B←X above the line list; i2spec.continuum).
It is turned into what each spectrometer reports: exp(-σN) for the stated I₂ column N, convolved
with the instrument function, then -ln(...)/N.
"""

import sys
from pathlib import Path

import numpy as np

from i2spec.continuum import continuum_model
from i2spec.intensity import intensity_model, master_line_list
from i2spec.spectrum import apparent_cross_section

EXT = Path(__file__).resolve().parents[1] / "data" / "external"
DATASETS = {  # conditions: docs/research/cross-section-data.md
    "spietz025": dict(file=EXT / "spietz_2006/I2DOASref_1200.TXT", T=298.0, N=6.86e15, ils=("nm", 0.25)),
    "spietz059": dict(file=EXT / "spietz_2006/I2DOASref_0300.TXT", T=298.0, N=1.42e16, ils=("nm", 0.59)),
    "saizlopez": dict(file=EXT / "mpi_mainz_i2/I2_Saiz-Lopez(2004)_295K_182-750nm.txt", T=295.0, N=3.1e16,
                      ils=("cm", 4.0)),
}
LORENTZ_FWHM = 0.24  # cm⁻¹: about 9.5 MHz/Torr for air or N₂ (Fletcher & McDaniel 1995; Wolf 2009) at ~750 Torr
RANGE_NM = (440.0, 640.0)


def two_columns(path):
    """(x, y) from a text file, skipping any line that does not start with two numbers."""
    rows = []
    for line in open(path, encoding="latin-1"):
        parts = line.replace(",", " ").split()
        try:
            rows.append((float(parts[0]), float(parts[1])))
        except (ValueError, IndexError):
            continue
    return np.array(rows).T


def main(name):
    d = DATASETS[name]
    lam, sigma_obs = two_columns(d["file"])
    sel = (lam >= max(RANGE_NM[0], lam.min() + 2)) & (lam <= min(RANGE_NM[1], lam.max() - 2))
    lam, sigma_obs = lam[sel], sigma_obs[sel]
    master = master_line_list(intensity_model("127I2"), 11000.0, 20100.0, T_range=(200.0, 600.0), S_min=1e-27)
    nu = 1e7 / lam
    lines = master.at(d["T"], nu.min() - 40, nu.max() + 40, S_min=1e-27)
    continuum = continuum_model("127I2", master.upper_cut, T_max=600.0).at(d["T"])

    kind, width = d["ils"]
    model, discrete = np.empty_like(lam), np.empty_like(lam)
    chunks = np.array_split(np.arange(lam.size), max(1, int(np.ceil(np.ptp(lam) / 5))))  # ~5 nm pieces
    for idx in chunks:
        w_cm = 1e7 * width / lam[idx].mean() ** 2 if kind == "nm" else width
        args = dict(ils=("gauss", w_cm), lorentz_fwhm=LORENTZ_FWHM, step=0.03)
        model[idx] = apparent_cross_section(nu[idx], lines, d["T"], d["N"], continuum=continuum, **args)
        discrete[idx] = apparent_cross_section(nu[idx], lines, d["T"], d["N"], **args)

    ratio = model / sigma_obs
    print(f"{name}: {lam.min():.1f}-{lam.max():.1f} nm, {lam.size} points, T = {d['T']} K, N = {d['N']:.3g} cm^-2; "
          f"{len(lines)} lines")
    print(f"  model/data: mean {ratio.mean():.3f}, rms scatter {ratio.std():.3f}; "
          f"continuum is {100 * np.mean((model - discrete) / model):.1f}% of the model on average")
    for lo in np.arange(np.floor(lam.min() / 10) * 10, lam.max(), 10):
        m = (lam >= lo) & (lam < lo + 10)
        if m.sum() > 3:
            print(f"    {lo:.0f}-{lo + 10:.0f} nm: model/data {ratio[m].mean():.3f} ± {ratio[m].std():.3f}")
    out = Path(__file__).resolve().parent / "out"
    out.mkdir(exist_ok=True)
    np.savetxt(out / f"xsec_{name}.txt", np.column_stack([lam, sigma_obs, model, discrete]), fmt="%.6e",
               header="lambda_vac_nm sigma_obs sigma_model sigma_model_lines_only  (cm^2)")


if __name__ == "__main__":
    main(sys.argv[1])
