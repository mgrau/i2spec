"""Line positions from a Doppler-limited FTS atlas, by fitting per-line shifts of the model template.

For each window the model transmission (master line list at the cell temperature, hyperfine patterns,
continuum, instrument function; compare_salami_ross.fit for column density, baseline, zero offset and
width) is the template, and every line deep enough to matter gets its own wavenumber shift δ_k, fitted
jointly by linearised least squares with a weak ridge prior. Blends come out with large σ(δ) from the
covariance instead of a wrong number. The result is the hyperfine-free position of each line on the
atlas's own wavenumber scale, with an uncertainty from the fit.

usage: atlas_lines.py --one CENTRE [T]          one 2 cm-1 window, printed
       atlas_lines.py [STEP WIDTH]              the whole atlas (default 2 2.4), to prototypes/out/atlas_lines_salami_ross.txt
"""
from __future__ import annotations

import sys
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from compare_salami_ross import CELL_LENGTH, convolve, fit                       # noqa: E402
from scan_salami_ross import SEGMENTS, S_HYPERFINE, _init, _state                # noqa: E402
from i2spec.constants import MHZ_PER_CM                                          # noqa: E402
from i2spec.model import RovibronicModel                                          # noqa: E402
from i2spec.spectrum import cross_section, hyperfine_patterns, number_density    # noqa: E402

MIN_DEPTH = 0.02          # a line gets a shift parameter if its own absorption exceeds this at line centre
PRIOR = 0.01              # cm-1, ridge prior on every shift
EPS = 0.001               # cm-1, finite-difference step for the Jacobian
OUT = Path(__file__).resolve().parent / "out" / "atlas_lines_salami_ross.txt"
HEADER = ("nu_model_cm-1 delta_cm-1 sigma_cm-1 depth S_cm v_upper v_lower J_lower branch T_K centre_cm-1 "
          "window_rms")


ATLASES = {
    # name: (loader, segments (lo, hi, T_K), output file)
    "salami_ross": (None, SEGMENTS, "atlas_lines_salami_ross.txt"),
    "apo": ("apo", ((15380.0, 20085.0, None),), "atlas_lines_apo.txt"),        # T fitted (see --temperature)
}
_atlas = {"name": "salami_ross", "T": None}


def load_apo():
    """The NIST 2-m FTS scan of the APO cell (data/external/apo_pyodine_atlas): air wavelengths (Å)
    and continuum-normalised flux -> (vacuum wavenumber, transmission), ascending in wavenumber."""
    import h5py
    from i2spec.lookup import air_to_vacuum
    with h5py.File(Path("data/external/apo_pyodine_atlas/APO_I2_from_FITS_cleaned_newscan_vnorm.h5")) as h:
        lam_air, flux = h["wavelength_air"][...] / 10.0, h["flux_normalized"][...]
    nu = 1e7 / air_to_vacuum(lam_air)
    ok = np.isfinite(nu) & np.isfinite(flux) & (flux > -0.5) & (flux < 2.0)   # the scan has NaN and spikes at its edges
    nu, flux = nu[ok], flux[ok]
    order = np.argsort(nu)
    return np.column_stack([nu[order], 100.0 * flux[order]])          # the scan code expects per cent


def init(atlas="salami_ross", T=None):
    """The scan's state plus a B-spline model for the hyperfine patterns: its levels are cached per J
    across lines and windows, where the DVR intensity model re-diagonalises for every line."""
    _init()
    _state["hfs"] = RovibronicModel("127I2", grids={"B": dict(rmin=2.35, rmax=8.0, h=0.01, order=10, nlev=70)})   # v' to 62
    _atlas.update(name=atlas, T=T)
    if atlas == "apo":
        _state["data"] = load_apo()


def extract(job):
    centre, width, T = job
    s = _state
    lo, hi = centre - width / 2, centre + width / 2
    sel = (s["data"][:, 0] >= lo) & (s["data"][:, 0] <= hi)
    nu_obs, t_obs = s["data"][sel, 0], s["data"][sel, 1] / 100
    if len(nu_obs) < 200 or not np.isfinite(t_obs).all() or t_obs.max() <= 0:
        return []                                    # a gap in the scan, or an unusable stretch
    lines = s["master"].at(T, lo - 0.5, hi + 0.5, S_min=1e-27)
    pattern, single = hyperfine_patterns(s["hfs"], lines, dJ=0), (np.zeros(1), np.ones(1))
    grid = np.arange(lo - 0.3, hi + 0.3, 0.0005)
    N0 = number_density(295.0, T) * CELL_LENGTH
    # which lines get a parameter: own peak absorption above MIN_DEPTH at the starting column density
    peak = lines.S * np.sqrt(4 * np.log(2) / np.pi) / (7.16e-7 * lines.nu * np.sqrt(T / 253.8))   # Doppler peak of a unit line
    strong = np.nonzero((1 - np.exp(-N0 * peak) > MIN_DEPTH) & (lines.nu >= lo) & (lines.nu <= hi))[0]
    hf = lambda k: pattern(k) if (lines.S[k] >= S_HYPERFINE or k in set(strong)) else single     # noqa: E731
    sig_line = {k: cross_section(grid, lines.at_index(k) if hasattr(lines, "at_index") else _one(lines, k), T,
                                 hyperfine=lambda _: hf(k)) for k in strong}
    sig_total = cross_section(grid, lines, T, hyperfine=hf) + s["continuum"].cross_section(grid, T)
    fix_offset = t_obs.min() > 0.7 * t_obs.max()
    (log_n, b0, b1, c0, shift, w), _ = fit(grid, sig_total, nu_obs, t_obs, N0, fix_offset)
    N, w, c0 = 10**log_n, max(abs(w), 1e-4), (0.0 if fix_offset else c0)
    x = nu_obs - nu_obs.mean()

    base = b0 + b1 * x

    def shifted(k, d):
        return np.interp(grid - d, grid, sig_line[k], left=0.0, right=0.0) if d != 0.0 else sig_line[k]

    def model(delta):
        sig = sig_total.copy()
        for k, d in zip(strong, delta):
            sig += shifted(k, d) - sig_line[k]
        e = np.exp(-N * sig)
        return base * np.interp(nu_obs, grid, convolve(grid, e, "gauss", w)) + c0, e

    delta = np.zeros(len(strong))
    for _ in range(2):
        t0, e = model(delta)
        # dT/d(delta_k) = base * ILS (x) [N sigma_k'(nu - delta_k) exp(-N sigma)]: one convolution per line
        J = np.empty((len(nu_obs), len(strong)))
        for i, k in enumerate(strong):
            g = np.gradient(shifted(k, delta[i]), grid)
            J[:, i] = base * np.interp(nu_obs, grid, convolve(grid, N * g * e, "gauss", w))
        r = t_obs - t0
        noise = max(np.std(np.diff(t_obs, 2)) / np.sqrt(6), 1e-4)
        A = J / noise
        M = A.T @ A + np.eye(len(strong)) / PRIOR**2
        step = np.linalg.solve(M, A.T @ (r / noise) - delta / PRIOR**2)
        delta = delta + step
    t0, _ = model(delta)
    r = t_obs - t0
    scale = max(1.0, np.sqrt(np.mean((r / noise) ** 2)))
    cov = np.linalg.inv(M) * scale**2
    rows = []
    for i, k in enumerate(strong):
        depth = 1 - np.exp(-N * peak[k])
        rows.append((lines.nu[k], delta[i], np.sqrt(cov[i, i]), depth, lines.S[k], lines.v_upper[k], lines.v_lower[k],
                     lines.J_lower[k], lines.branch[k], T, centre, np.std(r)))
    return rows


def _one(lines, k):
    from i2spec.intensity import LineList
    return LineList(lines.nu[k:k + 1], lines.S[k:k + 1], lines.v_upper[k:k + 1], lines.v_lower[k:k + 1],
                    lines.J_lower[k:k + 1], lines.branch[k:k + 1], lines.E_lower[k:k + 1], lines.T)


def safe_extract(job):
    """One bad window must not lose the whole run; its centre is reported and skipped."""
    try:
        return extract(job)
    except Exception as e:                                    # noqa: BLE001
        print(f"window {job[0]:.1f} failed: {type(e).__name__}: {e}", flush=True)
        return []


def jobs(step, width, atlas="salami_ross", T_fixed=None):
    segments = ATLASES[atlas][1]
    lo_all, hi_all = min(a for a, *_ in segments), max(b for _, b, _ in segments)
    for centre in np.arange(np.ceil(lo_all) + width, hi_all - width, step):
        for T in sorted({T_fixed or T for a, b, T in segments if a <= centre - width / 2 and centre + width / 2 <= b}):
            yield float(centre), width, T


def main(argv):
    opts = {a.split("=")[0]: a.split("=", 1)[1] for a in argv if a.startswith("--") and "=" in a}
    argv = [a for a in argv if not (a.startswith("--") and "=" in a)]
    atlas = opts.get("--atlas", "salami_ross")
    T_fixed = float(opts["--temperature"]) if "--temperature" in opts else None
    if argv[:1] == ["--one"]:
        centre = float(argv[1]); T = float(argv[2]) if len(argv) > 2 else (T_fixed or 293.15)
        init(atlas, T_fixed)
        t = time.time(); rows = extract((centre, 2.4, T))
        print(f"first window {time.time() - t:.1f} s (includes the level cache warm-up)")
        t = time.time(); rows = extract((centre + 2.0, 2.4, T))
        print(f"{len(rows)} lines in {time.time() - t:.1f} s; window rms {rows[0][-1]:.4f}" if rows else "no lines")
        for nu, d, sd, depth, S, vu, vl, J, br, *_ in sorted(rows):
            print(f"  {'R' if br > 0 else 'P'}({J:3d}) {vu:2d}-{vl:<2d} {nu:11.4f}  δ {d * MHZ_PER_CM:+8.1f} ± {sd * MHZ_PER_CM:6.1f} MHz  depth {depth:4.2f}")
        return
    step, width = (float(argv[0]), float(argv[1])) if len(argv) >= 2 else (2.0, 2.4)
    todo = list(jobs(step, width, atlas, T_fixed))
    out = OUT.parent / ATLASES[atlas][2]
    t = time.time()
    with Pool(initializer=init, initargs=(atlas, T_fixed)) as pool, open(out, "w") as f:
        f.write("# " + HEADER + "\n")
        for n, rows in enumerate(pool.imap_unordered(safe_extract, todo, chunksize=4), 1):
            for row in rows:
                # nu needs every digit: 6 significant figures would round it to 0.1 cm-1
                f.write("%.6f %.7f %.7f %.4f %.4g %d %d %d %d %.2f %.3f %.5f\n" % tuple(row))
            if n % 100 == 0:
                print(f"{n}/{len(todo)} windows, {time.time() - t:.0f} s", flush=True)
    print("done", time.time() - t)


if __name__ == "__main__":
    main(sys.argv[1:])
