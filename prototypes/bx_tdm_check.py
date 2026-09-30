"""Is the discrete B<-X band strength too high at 540-630 nm, and what B-X moment change would fix it?

Question: docs/research/tellinghuisen-components.md found the model's banded B<-X 2-5 % above Tellinghuisen's
pseudocontinuum and his 600-625 nm spectra, and fig. 9 of the paper shows 1.01-1.11 against Spietz 2006.
This script redoes the comparison per data set, at each set's temperature, column and resolution, splits
the excess into what saturation, wavelength scale, temperature and absolute anchors explain, maps it onto
the R axis of the moment, and fits a minimal correction mu(R) = mu_T11(R) * (1 + a z + b z^2), z = R - R0,
to the absolute cross sections. Note: docs/research/bx-band-strength.md.

Method
* Line strengths for any such correction, without rebuilding the line list: for each line,
  M_k = <v'J'| mu_T11 z^k |v''J''> (k = 0, 1, 2), from the same DVR eigenvectors master_line_list uses
  (checked to reproduce strength0 to 1e-12), so S(a, b) = S0 (1 + a r1 + b r2)^2 with r_k = M_k/M_0, exactly.
  R_bar = R0 + r1 is the line's R-centroid.
* The B<-X continuum (continuum.py uses the same mu) for the same correction: six continuum builds with
  mu * (1, 1+z, 1-z, 1+z^2, 1-z^2, 1+z+z^2) give the six products <E|mu z^i|v><E|mu z^j|v>, so the
  continuum is also exact in (a, b). The builds patch i2spec.continuum.mu_tellinghuisen2011 in this
  process only; no package code changes.
* Weak-absorption ("linear") model at each data point: every line's strength histogrammed on a 0.1 cm-1
  grid, plus the continuum, averaged over a Gaussian instrument function of the data set's FWHM. The
  optical-depth effect is then a multiplicative factor F = sigma_app(N)/sigma_weak computed with the full
  line-by-line -ln(ILS (x) exp(-sigma N))/N (spectrum.py) at nominal strengths (a few-% change of S
  changes F by a few % of 1 - F, negligible here). Tellinghuisen's column is unknown: F is reported for
  trial columns and not applied.
* Comparison: means over 10 nm bins (and 4 nm Gaussian smooths for the figure) of model and data at the
  data's own sampling points.

Run from the repository root (single process; about 6 minutes cold, seconds with the caches in
~/.cache/i2spec/bx_tdm_check/):

    I2SPEC_WORKERS=1 uv run --group research python prototypes/bx_tdm_check.py
"""
from __future__ import annotations

import hashlib
import json
import re
import sys
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import least_squares
from scipy.signal import fftconvolve

import i2spec.continuum as C
from i2spec.constants import DEFAULT_PARAMETERS
from i2spec.intensity import (S_UNIT, cache_dir, intensity_model, master_line_list, mu_tellinghuisen2011,
                              nuclear_spin_weight, shared_offset)
from i2spec.model import RovibronicModel
from i2spec.spectrum import air_to_vacuum, cross_section, doppler_fwhm, ils_kernel

ROOT = Path(__file__).resolve().parents[1]
EXT = ROOT / "data/external"
OUT_JSON = ROOT / "prototypes/out/bx_tdm_check.json"
OUT_PNG = ROOT / "docs/figures/bx_tdm_check.png"
CACHE = cache_dir() / "bx_tdm_check"
N_A = 6.02214076e23
EPS = 1000 * np.log(10) / N_A          # cm^2 per (l mol-1 cm-1)
R0 = 2.85                              # Å, origin of z (Tellinghuisen 2011 eq. 9)
LORENTZ = 0.24                         # cm-1 FWHM at ~1 bar N2/air (compare_cross_sections.py)
FINE = (11400.0, 25600.0, 0.1)         # cm-1 grid for the weak-limit model
BINS = np.arange(420.0, 781.0, 10.0)   # nm, comparison bins
T_SENS = (293.0, 303.0)                # sensitivity to the Spietz room temperature
#: Tellinghuisen's spectra are nominally 2 nm, but Table IS shows none of the ~4 nm band structure a 2 nm
#: Gaussian leaves at 600-630 nm (model/data alternates 0.86/1.20 point to point); 4 nm is used, 2 nm is a variant.
TELL_ILS = 4.0
TELL_N = (1e16, 3e16, 1e17, 3e17)      # cm-2, trial columns for Tellinghuisen's pure-vapour spectra
MU_VARIANTS = {"1": (1, 0, 0), "1+z": (1, 1, 0), "1-z": (1, -1, 0), "1+z2": (1, 0, 1), "1-z2": (1, 0, -1),
               "1+z+z2": (1, 1, 1)}


def log(*a):
    print(*a, flush=True)


def key(*args):
    return hashlib.sha1(repr((DEFAULT_PARAMETERS, R0) + args).encode()).hexdigest()[:12]


# ------------------------------------------------------------------------------------------------ inputs
def columns(path, n=2):
    rows = []
    for line in open(path, encoding="latin-1"):
        parts = line.replace(",", " ").split()
        try:
            rows.append([float(p) for p in parts[:n]])
        except (ValueError, IndexError):
            continue
    return np.array([r for r in rows if len(r) == n]).T


def table_IS():
    """Tellinghuisen 2011 [T11C] supplement Table IS: λ (air), ε and s.d. at 35.4 and 64.0 °C."""
    text = (EXT / "tellinghuisen_2011_supplement/JT-Supplement-JCP-I2.txt").read_text(encoding="utf-8")
    rows = [[float(s) for s in line.split()] for line in text.split("Table IIS.")[0].splitlines()
            if len(line.split()) == 5 and re.fullmatch(r"\d{3}", line.split()[0])]
    return np.array(rows)


def table_IIS_B():
    """Table IIS B<-X column at 35 °C (λ air, ε, LS s.e.), the pseudocontinuum used above 500 nm."""
    text = (EXT / "tellinghuisen_2011_supplement/JT-Supplement-JCP-I2.txt").read_text(encoding="utf-8")
    t2 = text.split("Table IIS.")[1].split("Least-squares normalized residuals")[0]
    out = []
    for line in t2.splitlines():
        v = line.split()
        if not v or not re.fullmatch(r"\d{3}", v[0]):
            continue
        x = [float(s) for s in v[1:]]
        if len(x) == 12:
            out.append((float(v[0]), x[4], x[5]))
        elif len(x) == 9:
            out.append((float(v[0]), x[1], x[2]))
        if v[0] == "850":
            break
    return np.array(out)


def datasets():
    """Absolute cross sections with their conditions (data/external/README.md, cross-section-data.md)."""
    d = {}
    lam, s = columns(EXT / "spietz_2006/I2DOASref_0300.TXT")
    d["spietz059"] = dict(lam_file=lam, obs=s, T=298.0, N=1.42e16, fwhm_nm=0.59, bath=True, scale_sd=0.027,
                          label="Spietz 2006, 0.59 nm", scale_note="pinned to sigma(500.0 nm air) = 2.191e-18")
    lam, s = columns(EXT / "spietz_2006/I2DOASref_1200.TXT")
    d["spietz025"] = dict(lam_file=lam, obs=s, T=298.0, N=6.86e15, fwhm_nm=0.25, bath=True, scale_sd=0.04,
                          label="Spietz 2006, 0.25 nm", scale_note="column transferred from 500 nm (2.191e-18)")
    lam, s = columns(EXT / "mpi_mainz_i2/I2_Saiz-Lopez(2004)_295K_182-750nm.txt")
    m = (lam >= 420) & (lam <= 760)
    d["saizlopez"] = dict(lam_file=lam[m], obs=s[m], T=295.0, N=3.1e16, fwhm_cm=4.0, bath=True, scale_sd=0.12,
                          label="Saiz-Lopez 2004, 4 cm-1", scale_note="vapour pressure, +-12 %")
    t = table_IS()
    for col, T, name in ((1, 308.55, "tell308"), (3, 337.15, "tell337")):
        d[name] = dict(lam_file=t[:, 0], obs=t[:, col] * EPS, sd=t[:, col + 1] * EPS, T=T, N=None, fwhm_nm=TELL_ILS,
                       bath=False, scale_sd=0.005, air=True,
                       label=f"Tellinghuisen 2011, 2 nm, {T - 273.15:.1f} °C", scale_note="stated ~0.5 %")
    return d


# ------------------------------------------------------------------------------------ moments per line
def line_moments(master):
    """r1 = M1/M0 and r2 = M2/M0 for every line of ``master`` (z = R - R0), from the DVR eigenvectors."""
    path = CACHE / f"moments_{key(len(master), float(master.nu.sum()))}.npz"
    if path.exists():
        z = np.load(path)
        return z["r1"], z["r2"], float(z["check"])
    model = intensity_model("127I2")
    X, B = model.states["X"], model.states["B"]
    k0 = shared_offset(X, B)
    R = B.R
    mu = mu_tellinghuisen2011(R)
    zz = R - R0
    order = np.lexsort((master.v_lower, master.v_upper, master.branch, master.J_lower))
    J_sorted = master.J_lower[order]
    starts = np.searchsorted(J_sorted, np.arange(master.J_lower.max() + 2))
    r1 = np.full(len(master), np.nan)
    r2 = np.full(len(master), np.nan)
    worst = 0.0
    upper = {}
    for J in range(master.J_lower.max() + 1):
        idx = order[starts[J]:starts[J + 1]]
        if idx.size == 0:
            continue
        _, cx = X.states(J)
        cx = cx[k0:]
        for br in (+1, -1):
            sel = idx[master.branch[idx] == br]
            if sel.size == 0:
                continue
            Jp = J + br
            if Jp not in upper:
                upper[Jp] = B.states(Jp)[1]
            cb = upper[Jp]
            vb, vx = master.v_upper[sel], master.v_lower[sel]
            M = [np.einsum("ri,ri->i", cb[:, vb], (mu * zz**k)[:, None] * cx[:, vx]) for k in range(3)]
            r1[sel], r2[sel] = M[1] / M[0], M[2] / M[0]
            s = J + 1 if br > 0 else J
            s0 = S_UNIT * master.nu[sel] * s * M[0] ** 2 * nuclear_spin_weight(J, "127I2")
            worst = max(worst, float(np.max(np.abs(s0 / master.strength0[sel] - 1))))
        upper.pop(J - 1, None)
        if J % 50 == 0:
            log(f"  moments J''={J}")
    CACHE.mkdir(parents=True, exist_ok=True)
    np.savez(path, r1=r1, r2=r2, check=worst)
    return r1, r2, worst


# --------------------------------------------------------------------------------- continuum products
def continuum_basis(master, temps):
    """σ_B basis on a 1 cm-1 grid: products P00, P01, P02, P11, P12, P22 (cm²) at each T."""
    nu = np.arange(FINE[0], FINE[1] + 1, 1.0)
    path = CACHE / f"contB_{key(tuple(temps), float(np.nansum(master.upper_cut)))}.npz"
    if path.exists():
        z = np.load(path)
        return nu, {T: z[f"T{T}"] for T in temps}
    base = C.mu_tellinghuisen2011
    raw = {}
    for name, (c0, c1, c2) in MU_VARIANTS.items():
        C.mu_tellinghuisen2011 = (lambda c0, c1, c2: lambda R: base(R) * (c0 + c1 * (np.asarray(R) - R0)
                                                                           + c2 * (np.asarray(R) - R0) ** 2))(c0, c1, c2)
        try:
            cont = C.continuum_model("127I2", master.upper_cut, T_max=600.0, cache=False)
        finally:
            C.mu_tellinghuisen2011 = base
        raw[name] = {T: cont.cross_section(nu, T, states=("B",)) for T in temps}
        log(f"  continuum mu*({name}) built")
    out = {}
    for T in temps:
        u = {k: raw[k][T] for k in raw}
        P00 = u["1"]
        P01 = (u["1+z"] - u["1-z"]) / 4
        P11 = (u["1+z"] + u["1-z"]) / 2 - P00
        P02 = (u["1+z2"] - u["1-z2"]) / 4
        P22 = (u["1+z2"] + u["1-z2"]) / 2 - P00
        P12 = (u["1+z+z2"] - P00 - P11 - P22 - 2 * P01 - 2 * P02) / 2
        out[T] = np.array([P00, P01, P02, P11, P12, P22])
    CACHE.mkdir(parents=True, exist_ok=True)
    np.savez(path, **{f"T{T}": v for T, v in out.items()})
    return nu, out


def coeffs(a, b, g=1.0):
    """Weights of the six products (00, 01, 02, 11, 12, 22) in |<mu g (1 + a z + b z²)>|²."""
    return g * g * np.array([1.0, 2 * a, 2 * b, a * a, 2 * a * b, b * b])


def fcorr(R, par):
    """mu_new / mu_T11 at R for par = (g, a, b)."""
    g, a, b = par
    z = np.asarray(R) - R0
    return g * (1 + a * z + b * z * z)


# ------------------------------------------------------------------------------------- weak-limit model
class WeakModel:
    """Line + continuum cross section on the fine grid, as six basis spectra, at one temperature."""

    def __init__(self, master, r1, r2, cont, cB_nu, cB, T):
        lo, hi, st = FINE
        self.edges = np.arange(lo, hi + st, st)
        self.nu = (self.edges[1:] + self.edges[:-1]) / 2
        L = master.at(T, S_min=0.0)
        S = L.S
        w = np.array([np.ones_like(r1), r1, r2, r1 * r1, r1 * r2, r2 * r2])
        self.lines = np.array([np.histogram(L.nu, self.edges, weights=S * wk)[0] / st for wk in w])
        self.cB = np.array([np.interp(self.nu, cB_nu, p) for p in cB[T]])
        self.cAC = cont.cross_section(self.nu, T, states=("A", "C"))
        self.T = T

    def sigma(self, a=0.0, b=0.0, parts=False):
        c = coeffs(a, b)
        lines, cb = c @ self.lines, c @ self.cB
        return (lines, cb, self.cAC) if parts else lines + cb + self.cAC

    def basis(self):
        """(8, n): six line+B-continuum product spectra, A+C continuum, nothing else."""
        return np.vstack([self.lines + self.cB, self.cAC[None, :]])


def ils_average(nu_fine, y, lam_vac, fwhm_nm=None, fwhm_cm=None):
    """Gaussian ILS average of y (rows) at vacuum wavelengths lam_vac; FWHM in nm or cm-1."""
    y = np.atleast_2d(y)
    out = np.empty((y.shape[0], lam_vac.size))
    for i, lv in enumerate(lam_vac):
        n0 = 1e7 / lv
        fw = fwhm_cm if fwhm_cm is not None else 1e7 * fwhm_nm / lv**2
        a, b = np.searchsorted(nu_fine, [n0 - 3 * fw, n0 + 3 * fw])
        g = np.exp(-4 * np.log(2) * ((nu_fine[a:b] - n0) / fw) ** 2)
        out[:, i] = y[:, a:b] @ g / g.sum()
    return out


# --------------------------------------------------------------------------- optical-depth (saturation)
_HFS = {}


def hfs_function(model, j_step=20):
    """Factory: sub-list -> callable k -> (offsets cm-1, weights), with one pattern cache for the whole run."""
    from i2spec.constants import MHZ_PER_CM

    def make(lines):
        def pattern(k):
            J = int(lines.J_lower[k])
            if J >= 20:
                Jr = j_step * round(J / j_step)
                J = Jr + 1 if (Jr - J) % 2 else Jr
            kk = (int(lines.v_upper[k]), int(lines.v_lower[k]), J, "R" if lines.branch[k] > 0 else "P")
            if kk not in _HFS:
                _, comps = model.hyperfine_components(*kk, dJ=0)
                _HFS[kk] = (np.array([c.offset for c in comps]) / MHZ_PER_CM, np.array([c.strength for c in comps]))
            return _HFS[kk]
        return pattern
    return make


def _sub(lines, sel, T):
    return type(lines)(lines.nu[sel], lines.S[sel], lines.v_upper[sel], lines.v_lower[sel], lines.J_lower[sel],
                       lines.branch[sel], lines.E_lower[sel], T=T)


def saturation_factor(lines, cont_fn, T, lam_vac, Ns, fwhm_nm=None, fwhm_cm=None, lorentz=LORENTZ,
                      hyperfine=None, S_hfs=0.0, block_nm=5.0):
    """{N: F} with F = sigma_app(N) / sigma_weak at lam_vac: full line-by-line, Gaussian ILS, blocks of block_nm.
    With ``hyperfine`` (a hfs_function), lines with S >= S_hfs get their hyperfine components."""
    F = {N: np.full(lam_vac.size, np.nan) for N in Ns}
    order = np.argsort(lam_vac)
    lam_s = lam_vac[order]
    for blk in np.array_split(np.arange(lam_s.size), max(1, int(np.ceil(np.ptp(lam_s) / block_nm)))):
        lv = lam_s[blk]
        fw = fwhm_cm if fwhm_cm is not None else 1e7 * fwhm_nm / lv.mean() ** 2
        nu_out = 1e7 / lv
        step = (float(doppler_fwhm(nu_out.min(), T)) + lorentz) / 8
        kern = ils_kernel(step, "gauss", fw)
        pad = (kern.size // 2) * step + 1.0
        grid = np.arange(nu_out.min() - pad, nu_out.max() + pad, step)
        inside = (lines.nu > grid[0] - 2) & (lines.nu < grid[-1] + 2)
        if hyperfine is None:
            sig = cross_section(grid, _sub(lines, inside, T), T, lorentz_fwhm=lorentz)
        else:
            strong = _sub(lines, inside & (lines.S >= S_hfs), T)
            sig = (cross_section(grid, strong, T, hyperfine=hyperfine(strong), lorentz_fwhm=lorentz)
                   + cross_section(grid, _sub(lines, inside & (lines.S < S_hfs), T), T, lorentz_fwhm=lorentz))
        sig = sig + cont_fn(grid)
        weak = np.interp(nu_out, grid, fftconvolve(sig, kern, mode="same"))
        for N in Ns:
            app = -np.log(np.interp(nu_out, grid, fftconvolve(np.exp(-sig * N), kern, mode="same"))) / N
            F[N][order[blk]] = app / weak
        log(f"    block {lv.min():.1f}-{lv.max():.1f} nm done")
    return F


# ------------------------------------------------------------------------------------------- helpers
def smooth(lam, y, fwhm):
    out = np.empty_like(y)
    for i, l0 in enumerate(lam):
        m = np.abs(lam - l0) < 1.5 * fwhm
        g = np.exp(-4 * np.log(2) * ((lam[m] - l0) / fwhm) ** 2)
        out[i] = (g * y[m]).sum() / g.sum()
    return out


def binned(lam, y, bins=BINS, min_n=3):
    idx = np.digitize(lam, bins) - 1
    cen, val = [], []
    for k in range(bins.size - 1):
        m = idx == k
        if m.sum() >= min_n:
            cen.append((bins[k] + bins[k + 1]) / 2)
            val.append(y[m].mean())
    return np.array(cen), np.array(val)


def bin_ratio(lam, model, obs, bins=BINS, min_n=3):
    """Ratio of bin means (model/obs) and the bin centres."""
    c, m = binned(lam, model, bins, min_n)
    _, o = binned(lam, obs, bins, min_n)
    return c, m / o


def highpass(lam, y, w=3.0):
    return y - smooth(lam, y, w)


# ----------------------------------------------------------------------------------------------- main
def main():
    t0 = time.time()
    CACHE.mkdir(parents=True, exist_ok=True)
    model_i = intensity_model("127I2")
    master = master_line_list(model_i, 11000.0, 20100.0, T_range=(200.0, 600.0), S_min=1e-27, workers=1)
    cont = C.continuum_model("127I2", master.upper_cut, T_max=600.0)
    log(f"{len(master)} lines, parameters {DEFAULT_PARAMETERS}")

    log("moments ...")
    r1, r2, check = line_moments(master)
    bad = ~np.isfinite(r1) | (np.abs(r1) > 1.0)      # lines at an FC node: centroid undefined; keep S0
    r1 = np.where(bad, 0.0, r1)
    r2 = np.where(bad, 0.0, r2)
    log(f"  M0 reproduces strength0 to {check:.1e}; {bad.sum()} lines with |r1| > 1 Å kept at S0 "
        f"({master.strength0[bad].sum() / master.strength0.sum():.1e} of the total strength0)")

    D = datasets()
    temps = sorted({d["T"] for d in D.values()} | set(T_SENS) | {273.15, 308.15})
    log("continuum basis ...")
    cB_nu, cB = continuum_basis(master, temps)
    W = {T: WeakModel(master, r1, r2, cont, cB_nu, cB, T) for T in temps}

    # --- wavelength scales ---------------------------------------------------------------------------
    for name, d in D.items():
        d["lam"] = air_to_vacuum(d["lam_file"]) if d.get("air") else d["lam_file"].copy()
    # Spietz 0.59 nm: shift against the model (cross-section-data.md §4 found -0.09 to -0.17 nm vs the atlas)
    d = D["spietz059"]
    lam_f = np.arange(500.0, 590.0, 0.02)
    mod_f = ils_average(W[298.0].nu, W[298.0].sigma(), lam_f, fwhm_nm=0.59)[0]
    shifts = np.arange(-0.30, 0.301, 0.005)
    win = (d["lam_file"] > 505) & (d["lam_file"] < 586)
    hp_obs = highpass(d["lam_file"][win], d["obs"][win])
    cc = []
    for s in shifts:
        m = np.interp(d["lam_file"][win] + s, lam_f, mod_f)
        cc.append(np.corrcoef(highpass(d["lam_file"][win], m), hp_obs)[0, 1])
    shift059 = float(shifts[int(np.argmax(cc))])
    d["lam"] = d["lam_file"] + shift059
    log(f"  Spietz 0.59 nm: best shift {shift059:+.3f} nm (corr {max(cc):.3f})")
    shift_check = {}
    for name in ("spietz025", "saizlopez"):
        d = D[name]
        win = (d["lam_file"] > 505) & (d["lam_file"] < 586) & (np.gradient(d["lam_file"]) < 0.25)
        lf = np.arange(d["lam_file"][win].min() - 1, d["lam_file"][win].max() + 1, 0.01)
        kw = dict(fwhm_nm=d["fwhm_nm"]) if "fwhm_nm" in d else dict(fwhm_cm=d["fwhm_cm"])
        mf = ils_average(W[d["T"]].nu, W[d["T"]].sigma(), lf, **kw)[0]
        hp = highpass(d["lam_file"][win], d["obs"][win])
        cc2 = [np.corrcoef(highpass(d["lam_file"][win], np.interp(d["lam_file"][win] + s, lf, mf)), hp)[0, 1]
               for s in shifts]
        shift_check[name] = float(shifts[int(np.argmax(cc2))])
    log(f"  shifts found (not applied): {shift_check}")

    # --- weak-limit basis at every data point ----------------------------------------------------------
    for name, d in D.items():
        kw = dict(fwhm_nm=d["fwhm_nm"]) if "fwhm_nm" in d else dict(fwhm_cm=d["fwhm_cm"])
        w = W[d["T"]]
        d["basis"] = ils_average(w.nu, w.basis(), d["lam"], **kw)
        if name.startswith("tell"):      # vacuum reading and 2 nm ILS, as sensitivities
            d["basis_vacread"] = ils_average(w.nu, w.basis(), d["lam_file"], **kw)
            d["basis_2nm"] = ils_average(w.nu, w.basis(), d["lam"], fwhm_nm=2.0)
        if name.startswith("spietz"):
            d["basis_T"] = {T: ils_average(W[T].nu, W[T].basis(), d["lam"], **kw) for T in T_SENS}
        if name == "spietz059":
            d["basis_noshift"] = ils_average(w.nu, w.basis(), d["lam_file"], **kw)

    for name in ("tell308", "tell337"):
        w = W[D[name]["T"]]
        D[name]["Bfrac"] = (ils_average(w.nu, w.lines[0], D[name]["lam"], fwhm_nm=TELL_ILS)[0]
                            / ils_average(w.nu, w.sigma(), D[name]["lam"], fwhm_nm=TELL_ILS)[0])

    def model_at(d, par=(1.0, 0.0, 0.0), which="basis"):
        B = d[which] if not isinstance(which, np.ndarray) else which
        g, a, b = par
        return coeffs(a, b, g) @ B[:6] + B[6]

    # --- optical-depth factors -------------------------------------------------------------------------
    log("optical-depth factors ...")
    sat_path = CACHE / f"sat_{key(shift059)}.npz"
    if sat_path.exists():
        z = np.load(sat_path, allow_pickle=True)
        sat = z["sat"].item()
    else:
        sat = {}
        for name in ("spietz059", "spietz025", "saizlopez"):
            d = D[name]
            m = d["lam"] > 498.0            # no lines below the B limit
            lines = master.at(d["T"], S_min=1e-27)
            kw = dict(fwhm_nm=d["fwhm_nm"]) if "fwhm_nm" in d else dict(fwhm_cm=d["fwhm_cm"])
            F = np.ones(d["lam"].size)
            F[m] = saturation_factor(lines, cont.at(d["T"]), d["T"], d["lam"][m], [d["N"]], **kw)[d["N"]]
            sat[name] = F
            log(f"  {name}: F {np.nanmin(F):.4f}-{np.nanmax(F):.4f}")
        # Tellinghuisen: pure I2 vapour, Doppler + hyperfine, 2 nm; column unknown -> a range of N
        hfs = hfs_function(RovibronicModel("127I2"))
        for name in ("tell308", "tell337"):
            d = D[name]
            m = (d["lam"] > 599.5) & (d["lam"] < 641)
            T = d["T"]
            lines = master.at(T, 1e7 / 646, 1e7 / 594, S_min=1e-30)
            Fs = saturation_factor(lines, cont.at(T), T, d["lam"][m], TELL_N, fwhm_nm=2.0, lorentz=0.0,
                                   hyperfine=hfs, S_hfs=3e-26, block_nm=10.0)
            for N in TELL_N:
                F = np.ones(d["lam"].size)
                F[m] = Fs[N]
                sat[f"{name}_N{N:.0e}"] = F
                log(f"  {name} N={N:.0e}: F(600-628) {np.mean(F[(d['lam'] > 600) & (d['lam'] < 628.5)]):.4f}")
        np.savez(sat_path, sat=np.array(sat, dtype=object))
    for name in ("spietz059", "spietz025", "saizlopez"):
        D[name]["F"] = sat[name]
    for name in ("tell308", "tell337"):
        D[name]["F"] = np.ones(D[name]["lam"].size)      # weak limit is the default for his data (see note)

    # --- nominal comparison --------------------------------------------------------------------------
    res = {"parameters": DEFAULT_PARAMETERS, "R0_A": R0, "moment_check_rel": check,
           "spietz059_shift_nm": shift059, "shift_check_nm_not_applied": shift_check, "datasets": {}}
    anchor = {}
    # σ(500.0 nm air) at 0.58 nm, 298 K (Spietz's anchor, 2.186(21)e-18; weighted mean 2.191(20)e-18)
    lam500 = np.array([float(air_to_vacuum(500.0))])

    B500 = ils_average(W[298.0].nu, W[298.0].basis(), lam500, fwhm_nm=0.58)[:, 0]

    def sigma500(par=(1.0, 0.0, 0.0)):
        g, a, b = par
        return float(coeffs(a, b, g) @ B500[:6] + B500[6])

    s500 = sigma500()
    anchor["model_sigma500air_298K"] = s500
    anchor["spietz_sigma500"] = 2.186e-18
    anchor["spietz_weighted_mean"] = 2.191e-18
    anchor["tellinghuisen_room_T_eps500"] = 590.0 * EPS
    log(f"  model sigma(500 air, 0.58 nm, 298 K) = {s500:.4e}  (Spietz 2.186e-18: {s500 / 2.186e-18:.4f}; "
        f"Tellinghuisen 590(4) eps = {590 * EPS:.4e}: {s500 / (590 * EPS):.4f})")
    res["anchor_500nm"] = anchor

    # R-centroid map: strength-weighted mean R_bar of the lines in each bin, at each T
    rmap = {}
    for T in (298.0, 308.55):
        L = master.at(T, S_min=0.0)
        lam_l = 1e7 / L.nu
        Rb = R0 + r1
        cen, mR, sdR = [], [], []
        for lo in np.arange(500.0, 700.0, 10.0):
            m = (lam_l >= lo) & (lam_l < lo + 10)
            if m.sum() < 10:
                continue
            w = L.S[m]
            mu_ = np.average(Rb[m], weights=w)
            cen.append(lo + 5)
            mR.append(mu_)
            sdR.append(np.sqrt(np.average((Rb[m] - mu_) ** 2, weights=w)))
        rmap[f"{T}"] = dict(lam_nm=cen, R_bar_A=mR, R_spread_A=sdR)
    res["R_centroid_map"] = rmap
    # continuum B<-X: effective R at 420-500 nm (share-weighted from the z products)
    wB = W[308.55]
    cb0, cb1 = wB.cB[0], wB.cB[1]
    lamc = np.arange(420.0, 505.0, 10.0)
    Rc = [R0 + float(np.interp(1e7 / l, wB.nu, cb1) / np.interp(1e7 / l, wB.nu, cb0)) for l in lamc]
    res["R_centroid_continuum_B"] = dict(lam_nm=lamc.tolist(), R_bar_A=Rc)
    log("  R-centroid (298 K): " + ", ".join(f"{c:.0f} nm {r:.3f}" for c, r in zip(rmap['298.0']['lam_nm'], rmap['298.0']['R_bar_A'])))
    log("  continuum B R-centroid: " + ", ".join(f"{c:.0f} {r:.3f}" for c, r in zip(lamc, Rc)))

    def summarize(par=(1.0, 0.0, 0.0)):
        out = {}
        for name, d in D.items():
            mod = model_at(d, par) * d["F"]
            lam, obs = d["lam"], d["obs"]
            ok = np.isfinite(mod) & (obs > 0)
            if name.startswith("tell"):
                ok &= obs > 5 * d["sd"]
            c, r = bin_ratio(lam[ok], mod[ok], obs[ok])
            out[name] = dict(bin_nm=c.tolist(), ratio=r.tolist())
        return out

    nominal = summarize()
    res["nominal_bins"] = nominal
    # Spietz spectra carry sigma(500.0 nm air) = 2.191e-18; on the model's own 500 nm value they become:
    fac = 2.191e-18 / s500
    res["spietz_on_model_500nm_anchor"] = {n: dict(bin_nm=nominal[n]["bin_nm"],
                                                   ratio=(np.array(nominal[n]["ratio"]) * fac).tolist())
                                           for n in ("spietz059", "spietz025")}
    res["spietz_anchor_factor"] = fac
    for name, v in nominal.items():
        log(f"  {name}: " + " ".join(f"{c:.0f}:{r:.3f}" for c, r in zip(v["bin_nm"], v["ratio"])))

    # --- sensitivities ---------------------------------------------------------------------------------
    sens = {}
    for name in ("spietz059", "spietz025", "saizlopez"):
        d = D[name]
        m0 = model_at(d)
        sens[name] = dict(optical_depth=dict(zip(*[np.round(x, 4).tolist() for x in binned(d["lam"], d["F"])])))
        if name.startswith("spietz"):
            for T in T_SENS:
                mT = model_at(d, which=d["basis_T"][T])
                c, r = bin_ratio(d["lam"], mT, m0)
                sens[name][f"T{T:.0f}_over_298"] = dict(zip(np.round(c).tolist(), np.round(r, 4).tolist()))
        if name == "spietz059":
            mns = model_at(d, which=d["basis_noshift"])
            c, r = bin_ratio(d["lam"], m0, d["obs"])
            c2, r2_ = bin_ratio(d["lam_file"], mns, d["obs"])
            sens[name]["ratio_without_shift"] = dict(zip(np.round(c2).tolist(), np.round(r2_, 4).tolist()))
    for name in ("tell308", "tell337"):
        d = D[name]
        mv = model_at(d, which=d["basis_vacread"])
        ok = d["obs"] > 5 * d["sd"]
        c, r = bin_ratio(d["lam_file"][ok], mv[ok], d["obs"][ok])
        sens[name] = dict(ratio_if_vacuum=dict(zip(np.round(c).tolist(), np.round(r, 4).tolist())))
        for lo, hi in ((600, 628.5), (629, 660), (660, 700), (420, 500)):
            m = (d["lam_file"] >= lo) & (d["lam_file"] < hi) & ok
            sens[name][f"window_{lo:.0f}-{hi:.0f}"] = dict(
                ils4nm=float(model_at(d)[m].mean() / d["obs"][m].mean()),
                ils2nm=float(model_at(d, which=d["basis_2nm"])[m].mean() / d["obs"][m].mean()),
                vacuum_reading=float(mv[m].mean() / d["obs"][m].mean()),
                B_line_fraction=float(np.mean(d["Bfrac"][m])))
        blk = (d["lam"] > 600) & (d["lam"] < 629)
        sens[name]["optical_depth_600_628"] = {
            f"{N:.0e}": float(np.mean(sat[f"{name}_N{N:.0e}"][blk])) for N in TELL_N}
        blk2 = (d["lam"] > 629) & (d["lam"] < 641)
        sens[name]["optical_depth_630_640"] = {
            f"{N:.0e}": float(np.mean(sat[f"{name}_N{N:.0e}"][blk2])) for N in TELL_N}
        # the absorbance at 500 nm and at the 600-628 nm band that each N would give
        sens[name]["absorbance10_500nm_per_N"] = {f"{N:.0e}": float(s500 * N / np.log(10)) for N in TELL_N}
    res["sensitivity"] = sens
    log(json.dumps(sens, indent=0)[:3000])

    # --- fits ------------------------------------------------------------------------------------------
    # Each data set enters as its 10 nm bin means; the per-bin error is the larger of 1 % and the stated
    # per-point error / sqrt(n) (Tellinghuisen: his s.d.). Scale factors: one per data set, with a Gaussian
    # prior of the stated scale uncertainty ("absolute" fit) or free ("shape" fit). Windows:
    WIN = {"spietz059": [(490, 588)], "spietz025": [(545, 576)], "saizlopez": [(500, 640)],
           "tell308": [(420, 500), (600, 700)], "tell337": [(420, 500), (600, 700)]}
    WIN_BLUE = {"spietz059": [(440, 588)]}      # adds the Spietz window-deposit region, as a variant

    def prep(win):
        P = {}
        for name, d in D.items():
            m = np.zeros(d["lam"].size, bool)
            for lo, hi in win.get(name, WIN[name]):
                m |= (d["lam"] >= lo) & (d["lam"] < hi)
            if name.startswith("tell"):
                m &= d["obs"] > 5 * d["sd"]
            idx = np.digitize(d["lam"], BINS) - 1
            groups = [np.flatnonzero(m & (idx == k)) for k in range(BINS.size - 1)]
            groups = [g for g in groups if g.size >= 3]
            P[name] = groups
        return P

    FORMS = {"scale": ("g",), "linear": ("g", "a"), "quadratic": ("g", "a", "b"), "slope_fixed_2.85": ("a",)}

    def unpack(p, form):
        v = dict(g=1.0, a=0.0, b=0.0)
        v.update(zip(FORMS[form], p))
        return (v["g"], v["a"], v["b"])

    def residuals(p, names, P, prior, form):
        k = len(FORMS[form])
        par = unpack(p[:k], form)
        out = []
        for i, name in enumerate(names):
            d = D[name]
            s = p[k + i]
            mod = model_at(d, par) * d["F"]
            for g in P[name]:
                o, mm = d["obs"][g].mean(), mod[g].mean()
                err = 0.01
                if name.startswith("tell"):
                    err = max(0.005, float(np.sqrt(np.mean(d["sd"][g] ** 2) / g.size) / o))
                out.append((mm / (s * o) - 1) / err)
            if prior:
                out.append((s - 1) / d["scale_sd"])
        return np.array(out)

    Rg = np.array([2.45, 2.5, 2.55, 2.6, 2.63, 2.66, 2.7, 2.75, 2.8, 2.85, 2.9, 2.93, 2.95, 3.0, 3.1, 3.3])
    fits = {}
    names = list(D)
    noSL = [n for n in names if n != "saizlopez"]
    TELL = ["tell308", "tell337"]
    WIN_TELL_BAND = {"tell308": [(420, 500), (600, 629)], "tell337": [(420, 500), (600, 629)]}
    runs = [("scale", "scale", True, {}, noSL), ("linear", "linear", True, {}, noSL),
            ("quadratic", "quadratic", True, {}, noSL), ("linear_shape_only", "linear", False, {}, noSL),
            ("linear_with_saizlopez", "linear", True, {}, names),
            ("linear_with_blue_spietz", "linear", True, WIN_BLUE, noSL),
            ("slope_fixed_2.85", "slope_fixed_2.85", True, {}, noSL),
            ("scale_tell", "scale", True, {}, TELL), ("linear_tell", "linear", True, {}, TELL),
            ("linear_tell_600-629_420-500", "linear", True, WIN_TELL_BAND, TELL),
            ("quadratic_tell", "quadratic", True, {}, TELL),
            ("scale_spietz", "scale", True, {}, ["spietz059", "spietz025"]),
            ("linear_spietz", "linear", True, {}, ["spietz059", "spietz025"]),
            ("linear_spietz059", "linear", True, {}, ["spietz059"])]
    for fname, form, prior, win, use in runs:
        P = prep(win)
        k = len(FORMS[form])
        p0 = np.r_[[1.0 if n == "g" else 0.0 for n in FORMS[form]], np.ones(len(use))]
        r = least_squares(residuals, p0, args=(use, P, prior, form))
        dof = max(1, r.fun.size - r.x.size)
        chi2 = float((r.fun ** 2).sum())
        cov = np.linalg.pinv(r.jac.T @ r.jac) * max(1.0, chi2 / dof)
        chi2_nom = float((least_squares(lambda q: residuals(np.r_[p0[:k], q], use, P, prior, form),
                                        np.ones(len(use))).fun ** 2).sum())
        par = unpack(r.x[:k], form)
        # uncertainty of mu_new/mu at each R from the covariance of the shape parameters
        Jm = np.array([[(fcorr(R, unpack(r.x[:k] + dp, form)) - fcorr(R, par)) / 1e-6
                        for dp in np.eye(k) * 1e-6] for R in Rg])
        f_se = np.sqrt(np.einsum("ri,ij,rj->r", Jm, cov[:k, :k], Jm))
        fits[fname] = dict(form=form, datasets=use, params=dict(zip(FORMS[form], r.x[:k].tolist())),
                           params_se=dict(zip(FORMS[form], np.sqrt(np.diag(cov)[:k]).tolist())),
                           cov_params=cov[:k, :k].tolist(),
                           scales={n: float(v) for n, v in zip(use, r.x[k:])},
                           scales_se={n: float(np.sqrt(cov[k + i, k + i])) for i, n in enumerate(use)},
                           chi2=chi2, chi2_nominal=chi2_nom, dof=dof,
                           mu_ratio=dict(R_A=Rg.tolist(), value=fcorr(Rg, par).tolist(), se=f_se.tolist()),
                           sigma500_ratio=sigma500(par) / s500 if "spietz059" in D else None,
                           bins=summarize(par))
        log(f"  fit {fname:32s} " + " ".join(f"{n}={v:+.4f}±{e:.4f}" for n, v, e in
                                             zip(FORMS[form], r.x[:k], np.sqrt(np.diag(cov)[:k])))
            + f" chi2 {chi2:.1f} (nominal {chi2_nom:.1f}) dof {dof}; mu ratio at 2.55/2.70/2.80/2.93 Å: "
            + "/".join(f"{fcorr(R, par):.4f}" for R in (2.55, 2.70, 2.80, 2.93))
            + "; scales " + " ".join(f"{n}:{v:.3f}" for n, v in zip(use, r.x[k:])))
    res["fits"] = fits

    # --- effect of each candidate elsewhere --------------------------------------------------------------
    tb = table_IIS_B()
    lv_b = air_to_vacuum(tb[:, 0])
    w = W[308.15]
    Bn = ils_average(w.nu, w.lines + w.cB, lv_b, fwhm_nm=4.0)
    eff = {}
    for fname in ("scale", "linear", "quadratic", "linear_tell", "scale_tell", "linear_spietz059", "quadratic_tell"):
        par = unpack(list(fits[fname]["params"].values()), fits[fname]["form"])
        e = dict(par=par, sigma500air_298K_ratio=sigma500(par) / s500,
                 dube2004=dict(R_A=2.93, mu_D=1.10, mu_se=0.03, model_nominal=float(mu_tellinghuisen2011(2.93)),
                               model_new=float(mu_tellinghuisen2011(2.93) * fcorr(2.93, par))))
        for lab in TELL:
            d = D[lab]
            for lo, hi in ((420, 500), (600, 629), (629, 660), (660, 700), (700, 780)):
                m = (d["lam"] >= lo) & (d["lam"] < hi) & (d["obs"] > 5 * d["sd"])
                e[f"{lab}_{lo}-{hi}"] = dict(nominal=float(model_at(d)[m].mean() / d["obs"][m].mean()),
                                            new=float(model_at(d, par)[m].mean() / d["obs"][m].mean()))
        for lab, (lo, hi) in (("spietz059", (490, 520)), ("spietz059", (520, 560)), ("spietz059", (560, 588)),
                              ("spietz059", (440, 490)), ("spietz025", (545, 578))):
            d = D[lab]
            m = (d["lam"] >= lo) & (d["lam"] < hi)
            e[f"{lab}_{lo}-{hi}"] = dict(nominal=float((model_at(d) * d["F"])[m].mean() / d["obs"][m].mean()),
                                        new=float((model_at(d, par) * d["F"])[m].mean() / d["obs"][m].mean()))
        for lo, hi in ((420, 495), (525, 560), (565, 625)):
            m = (tb[:, 0] >= lo) & (tb[:, 0] <= hi)
            e[f"T11C_TableIIS_B_{lo}-{hi}"] = dict(
                nominal=float(np.mean((coeffs(0, 0) @ Bn)[m] / (tb[m, 1] * EPS))),
                new=float(np.mean((coeffs(par[1], par[2], par[0]) @ Bn)[m] / (tb[m, 1] * EPS))))
        eff[fname] = e
    eff["mu_T11_eq10_D"] = dict(zip(map(str, Rg.tolist()), mu_tellinghuisen2011(Rg).tolist()))
    eff["mu_T11_eq9_over_eq10_squared"] = dict(zip(map(str, Rg.tolist()),
                                                   ((1.1123 + 0.712 * (Rg - 2.85)) / mu_tellinghuisen2011(Rg)) ** 2))
    res["effects"] = eff
    for fname in ("scale", "linear", "quadratic", "linear_tell"):
        log(f"  effects {fname}: " + ", ".join(f"{k} {v['nominal']:.3f}->{v['new']:.3f}"
                                             for k, v in eff[fname].items() if isinstance(v, dict) and "new" in v)
            + f"; sigma500 x{eff[fname]['sigma500air_298K_ratio']:.4f}; Dubé mu {eff[fname]['dube2004']['model_new']:.3f}")

    # --- spot check: strengths from a rebuilt mu equal the moment prediction -----------------------------
    par_c = unpack(list(fits["quadratic"]["params"].values()), "quadratic")
    X, Bst = model_i.states["X"], model_i.states["B"]
    k0 = shared_offset(X, Bst)
    Rgrid = Bst.R
    mu_new = mu_tellinghuisen2011(Rgrid) * fcorr(Rgrid, par_c)
    worst = 0.0
    for J in (10, 60, 120):
        cx = X.states(J)[1][k0:]
        cb = Bst.states(J + 1)[1]
        Mn = cb.T @ (mu_new[:, None] * cx)
        sel = (master.J_lower == J) & (master.branch == 1)
        vb, vx = master.v_upper[sel], master.v_lower[sel]
        pred = master.strength0[sel] * (par_c[0] * (1 + par_c[1] * r1[sel] + par_c[2] * r2[sel])) ** 2
        direct = S_UNIT * master.nu[sel] * (J + 1) * Mn[vb, vx] ** 2 * nuclear_spin_weight(J, "127I2")
        big = master.strength0[sel] > 1e-3 * master.strength0[sel].max()
        worst = max(worst, float(np.max(np.abs(pred[big] / direct[big] - 1))))
    res["spot_check_rebuilt_strengths_max_rel_diff"] = worst
    log(f"  rebuilt-mu spot check: max rel diff {worst:.1e}")
    best_name = "linear"
    best = fits[best_name]
    par_best = unpack(list(best["params"].values()), best["form"])

    # --- figure ----------------------------------------------------------------------------------------
    fig, axs = plt.subplots(3, 1, figsize=(8.5, 10.5), gridspec_kw=dict(height_ratios=[1.3, 1.3, 1]))
    sty = {"spietz059": ("k", "-", None), "spietz025": ("#0072B2", "-", None), "saizlopez": ("0.6", "-", None),
           "tell308": ("#D55E00", "none", "o"), "tell337": ("#009E73", "none", "^")}
    for ax, (par, title) in zip(axs[:2], (((1.0, 0.0, 0.0), f"{DEFAULT_PARAMETERS} as is (μ: Tellinghuisen 2011 eq. 10)"),
                                        (par_best, f"candidate: μ = μ_T11 × {par_best[0]:.4f} (1 {par_best[1]:+.4f} (R − 2.85 Å))  "
                                                   "(fit 'linear', Spietz + Tellinghuisen)"))):
        ax.axhline(1, color="0.5", lw=0.6)
        for name, d in D.items():
            col, ls, mk = sty[name]
            mod = model_at(d, par) * d["F"]
            if name.startswith("tell"):
                ok = d["obs"] > 5 * d["sd"]
                for blk in (d["lam"] < 550, d["lam"] > 550):
                    l = d["lam"][ok & blk]
                    r = smooth(l, mod[ok & blk], 4.0) / smooth(l, d["obs"][ok & blk], 4.0)
                    ax.plot(l, r, ls="none", marker=mk, ms=3, mfc="white", mec=col,
                            label=d["label"] if blk[0] else None)
            else:
                l = d["lam"]
                r = smooth(l, mod, 4.0) / smooth(l, d["obs"], 4.0)
                e = (l > l.min() + 4) & (l < l.max() - 4)
                ax.plot(l[e], r[e], color=col, ls=ls, lw=1.1, label=d["label"])
                if name == "spietz059":
                    ax.plot(l[e], r[e] * 2.191e-18 / sigma500(par), color=col, ls="--", lw=0.9,
                            label="Spietz 0.59 nm on the model's σ(500)")
        ax.axhline(sigma500(par) / 2.191e-18, color="k", ls=":", lw=0.8,
                   label="model / Spietz σ(500 nm air)")
        ax.set_xlim(420, 780)
        ax.set_ylim(0.8, 1.2)
        ax.set_ylabel("model / measured (4 nm avg)")
        ax.set_title(title, fontsize=9)
        ax.legend(fontsize=7, ncol=3, loc="lower center", frameon=False)
    ax = axs[2]
    Rg2 = np.linspace(2.45, 3.1, 200)
    ax.axhline(1, color="0.5", lw=0.6)
    ax.plot(Rg2, ((1.1123 + 0.712 * (Rg2 - 2.85)) / mu_tellinghuisen2011(Rg2)) ** 2, "k--", lw=0.8,
            label="T11 eq. 9 (linear fit to his 520–640 nm spectra)")
    for fname, col in (("linear", "#D55E00"), ("quadratic", "#CC79A7"), ("linear_tell", "#009E73"),
                       ("linear_spietz059", "#0072B2")):
        f = fits[fname]
        par = unpack(list(f["params"].values()), f["form"])
        se = np.interp(Rg2, Rg, f["mu_ratio"]["se"])
        fc = fcorr(Rg2, par)
        ax.plot(Rg2, fc ** 2, color=col, lw=1.1, label=f"fit '{fname}'")
        ax.fill_between(Rg2, (fc - se) ** 2, (fc + se) ** 2, color=col, alpha=0.18, lw=0)
    mu_d = float(mu_tellinghuisen2011(2.93))
    ax.errorbar([2.93], [(1.10 / mu_d) ** 2], yerr=[2 * 1.10 * 0.03 / mu_d ** 2], fmt="s", color="k", ms=4,
                label="Dubé & Trinczek 2004 (718 nm lines)")
    rm = rmap["298.0"]
    ax.set_xlim(2.45, 3.1)
    ax.set_ylim(0.8, 1.1)
    ax2 = ax.twiny()
    ax2.set_xlim(2.45, 3.1)
    ticks = [l for l in (505, 540, 580, 620, 660, 690) if l <= max(rm["lam_nm"])]
    ax2.set_xticks(np.interp(ticks, rm["lam_nm"], rm["R_bar_A"]))
    ax2.set_xticklabels([str(t) for t in ticks], fontsize=7)
    ax2.set_xlabel("λ (nm) whose bands have that strength-weighted R-centroid (298 K); B continuum 420–500 nm ↔ 2.49–2.63 Å",
                   fontsize=8)
    ax.set_xlabel("R (Å)")
    ax.set_ylabel("μ² / μ²(T11 eq. 10)")
    ax.legend(fontsize=7, frameon=False, ncol=2, loc="lower left")
    fig.tight_layout()
    OUT_PNG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PNG, dpi=150)

    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    OUT_JSON.write_text(json.dumps(res, indent=1, default=float))
    log(f"wrote {OUT_JSON} and {OUT_PNG} in {time.time() - t0:.0f} s")


if __name__ == "__main__":
    sys.exit(main())
