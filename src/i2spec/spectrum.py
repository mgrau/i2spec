"""Absorption cross sections and cell transmission from a line list."""

from __future__ import annotations

from math import log, pi, sqrt

import numpy as np
from scipy.signal import fftconvolve
from scipy.special import voigt_profile

from .constants import ATOMIC_MASS, MHZ_PER_CM

_KB, _C, _U = 1.380649e-23, 299792458.0, 1.66053906892e-27
_TORR = 133.322368  # Pa


def doppler_fwhm(nu, T, isotopologue_mass=2 * ATOMIC_MASS[127]):
    """Doppler FWHM in cm⁻¹ (≈ 4.494e-8 ν √T for ¹²⁷I₂)."""
    return nu * np.sqrt(8 * _KB * T * log(2) / (isotopologue_mass * _U * _C**2))


def vapor_pressure(T):
    """I₂ vapour pressure over the solid in Pa (Tellinghuisen 2011 eq. 8, T in K)."""
    return 10 ** (12.1891 - 0.001301 * T - 0.3523 * np.log10(T) - 3410.71 / T) * _TORR


def number_density(T_coldfinger, T_cell):
    """I₂ molecules per cm³ in a cell whose pressure is set by its cold finger."""
    return vapor_pressure(T_coldfinger) / (_KB * T_cell) * 1e-6


def air_refractive_index(lam_vac_nm):
    """Refractive index of standard air at a vacuum wavelength in nm (Edlén 1966, as given by Morton 2000)."""
    s2 = (1e3 / np.asarray(lam_vac_nm, dtype=float)) ** 2
    return 1 + 8.34254e-5 + 2.406147e-2 / (130 - s2) + 1.5998e-4 / (38.9 - s2)


def air_to_vacuum(lam_air_nm):
    lam = np.asarray(lam_air_nm, dtype=float)
    return lam * air_refractive_index(lam * air_refractive_index(lam))


def vacuum_to_air(lam_vac_nm):
    lam = np.asarray(lam_vac_nm, dtype=float)
    return lam / air_refractive_index(lam)


def hyperfine_patterns(model, lines, dJ=0, j_step=None):
    """Callable for cross_section(): hyperfine offsets (cm⁻¹) and relative strengths of line k.

    j_step: for J'' >= 20, reuse the pattern of the nearest J'' on a grid of that step, keeping the
    parity of J'' (which fixes the number of components). Patterns change slowly with J, so this is
    fine for Doppler-limited and lower-resolution spectra.
    """
    cache = {}

    def pattern(k):
        J = int(lines.J_lower[k])
        if j_step and J >= 20:
            J_ref = j_step * round(J / j_step)
            J = J_ref + 1 if (J_ref - J) % 2 else J_ref
        key = (int(lines.v_upper[k]), int(lines.v_lower[k]), J, "R" if lines.branch[k] > 0 else "P")
        if key not in cache:
            _, comps = model.hyperfine_components(*key, dJ=dJ)
            cache[key] = (np.array([c.offset for c in comps]) / MHZ_PER_CM, np.array([c.strength for c in comps]))
        return cache[key]

    return pattern


def cross_section(nu_grid, lines, T, *, hyperfine=None, lorentz_fwhm=0.0, reach=6.0):
    """Absorption cross section σ(ν) in cm² on nu_grid (cm⁻¹).

    Each line (or hyperfine component) gets a Gaussian Doppler profile at temperature T, or a
    Voigt profile if lorentz_fwhm (cm⁻¹) > 0. hyperfine: None, or a callable k -> (offsets, weights),
    e.g. hyperfine_patterns().
    """
    nu_grid = np.asarray(nu_grid, dtype=float)
    sigma = np.zeros_like(nu_grid)
    for k in range(len(lines)):
        offsets, weights = (np.zeros(1), np.ones(1)) if hyperfine is None else hyperfine(k)
        g = doppler_fwhm(lines.nu[k], T)
        half = reach * (g + lorentz_fwhm) + np.abs(offsets).max()
        lo, hi = np.searchsorted(nu_grid, [lines.nu[k] - half, lines.nu[k] + half])
        if lo == hi:
            continue
        x = nu_grid[lo:hi, None] - (lines.nu[k] + offsets)[None, :]
        if lorentz_fwhm > 0:
            profile = voigt_profile(x, g / (2 * sqrt(2 * log(2))), lorentz_fwhm / 2)
        else:
            profile = sqrt(4 * log(2) / pi) / g * np.exp(-4 * log(2) * (x / g) ** 2)
        if 0 < lo and hi < nu_grid.size:
            # keep each line's integrated strength despite truncating the profile (Lorentzian wings)
            profile = profile / (profile.sum(axis=0) * (nu_grid[1] - nu_grid[0]))
        sigma[lo:hi] += lines.S[k] * profile @ weights
    return sigma


def transmission(sigma, column_density):
    """exp(-σ N) for column density N in molecules/cm²."""
    return np.exp(-sigma * column_density)


def convolve_gaussian(nu_grid, y, fwhm):
    """Convolve y, sampled on a uniform grid, with a unit-area Gaussian of the given FWHM (cm⁻¹)."""
    step = nu_grid[1] - nu_grid[0]
    n = int(np.ceil(3 * fwhm / step))
    x = np.arange(-n, n + 1) * step
    kernel = np.exp(-4 * log(2) * (x / fwhm) ** 2)
    kernel /= kernel.sum()
    return np.convolve(np.pad(y, n, mode="edge"), kernel, mode="valid")


def ils_kernel(step, kind="gauss", width=0.1):
    """Instrument function sampled at ``step`` (cm⁻¹), normalized to unit sum.

    kind:
      "gauss"   Gaussian; width = FWHM.
      "boxcar"  rectangular slit function; width = full width.
      "sinc"    unapodized FTS; width = 1/(2·MOPD), the usual nominal resolution. The FWHM is 1.207·width.
      "sinc2"   FTS with triangular apodization, same width convention. The FWHM is 1.772·width.
    The sinc kernels are truncated at ±40 widths.
    """
    half = {"gauss": 3 * width, "boxcar": width / 2, "sinc": 40 * width, "sinc2": 40 * width}[kind]
    n = int(np.ceil(half / step))
    x = np.arange(-n, n + 1) * step
    if kind == "gauss":
        k = np.exp(-4 * log(2) * (x / width) ** 2)
    elif kind == "boxcar":
        k = (np.abs(x) <= width / 2 + 1e-12).astype(float)
    elif kind == "sinc":
        k = np.sinc(x / width)
    else:
        k = np.sinc(x / (2 * width)) ** 2
    return k / k.sum()


def apparent_cross_section(nu_out, lines, T, column_density, ils=("gauss", 0.1), *, hyperfine=None,
                           continuum=None, step=None, lorentz_fwhm=0.0):
    """Cross section a spectrometer would report: -ln(ILS ⊗ exp(-σN)) / N, on nu_out (cm⁻¹).

    Strong lines saturate at column density N (cm⁻²), so at finite resolution the apparent cross
    section depends on N (Spietz et al. 2006). ils = (kind, width) as in ils_kernel(). continuum:
    optional callable ν -> σ_c(ν) in cm², added to the line cross section. lorentz_fwhm (cm⁻¹):
    pressure broadening. The fine grid step defaults to 1/8 of the line width.
    """
    nu_out = np.asarray(nu_out, dtype=float)
    step = (float(doppler_fwhm(nu_out.min(), T)) + lorentz_fwhm) / 8 if step is None else step
    kernel = ils_kernel(step, *ils)
    pad = (kernel.size // 2) * step + 1.0
    grid = np.arange(nu_out.min() - pad, nu_out.max() + pad, step)
    sigma = cross_section(grid, lines, T, hyperfine=hyperfine, lorentz_fwhm=lorentz_fwhm)
    if continuum is not None:
        sigma = sigma + continuum(grid)
    transmitted = fftconvolve(np.exp(-sigma * column_density), kernel, mode="same")
    return -np.log(np.interp(nu_out, grid, transmitted)) / column_density
