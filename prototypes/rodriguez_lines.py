"""Line centres from the Rodríguez Fernández / Lefrán Torres diode-laser iodine spectra, 14 400-14 710 cm-1.

The spectra (data/external/rodriguez_fernandez_marcassa/, sent by the authors, not for redistribution)
are Doppler-limited transmission scans read on a HighFinesse WS7 wavemeter: 0.5 cm-1 scans joined into
two files. The first file's axis steps backwards at 207 of the joins; it is sorted here, and those
backward steps mark scan boundaries (the other joins happen to be monotonic and are invisible).

The fit is prototypes/atlas_lines.py extract, as for the FTS atlases, in 1.0 cm-1 windows stepped by
0.8 cm-1 (short, because each 0.5 cm-1 scan may have its own baseline and offset). There is no instrument
function (a free-running diode laser); the fitted Gaussian width should come out near zero.
The axis is shifted by PRESHIFT_MHZ before the fits so the ridge prior does not pull against the
wavemeter offset; the dataset step puts it back.

Notes: docs/research/kitt-peak-and-rodriguez-spectra.md.

usage: rodriguez_lines.py --tscan                      the cell temperature, from window rms
       rodriguez_lines.py [--T=...] [--workers=2]       window fits -> prototypes/out/atlas_lines_rodriguez.txt
       rodriguez_lines.py --dataset                     data/atlas_lines/rodriguez_fernandez_2023.{csv,toml}
"""
from __future__ import annotations

import os
import sys
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np

os.environ.setdefault("I2SPEC_WORKERS", "2")
sys.path.insert(0, str(Path(__file__).resolve().parent))
import atlas_lines                                                           # noqa: E402
from kitt_peak_lines import (ATLAS, OUT, ROW_FMT, build, calibrate_sigma, clipped_fit, compare,  # noqa: E402
                             model_centres, reference_atlases, robust, drop_blends)
from i2spec.constants import MHZ_PER_CM                                     # noqa: E402

RAW = Path(__file__).resolve().parents[1] / "data" / "external" / "rodriguez_fernandez_marcassa"
FILES = ("Iodine spectrum 14400-14600 cm-1.txt", "Iodine spectrum 14600-14710 cm-1.txt")
PRESHIFT_MHZ = 131.0          # the first comparison (data/external/README.md); taken out, then put back
T_CELL = 350.0                # K: lowest mean window rms in --tscan (295-500 K; 0.79 % at 350 K)
OUTFILE = OUT / "atlas_lines_rodriguez.txt"
CUTS = dict(sigma=30.0, rms=0.04, depth=0.03)
#: blends: about half the Doppler FWHM, neighbours down to 0.3 of the line's strength (kept lines with such a
#: neighbour read +6 to +9 MHz against -1 MHz for the others; weaker neighbours made no difference)
BLEND = (0.01, 0.3)
MODEL_SIGMA_MHZ = 1.0         # i2spec2026m at 14 400-14 710 cm-1 (data/external/README.md)


def read(name):
    return np.loadtxt(RAW / name, skiprows=3)


def scan_edges():
    """Wavenumbers at which the first file's axis steps backwards: known joins between 0.5 cm-1 scans."""
    d = read(FILES[0])
    back = np.nonzero(np.diff(d[:, 0]) <= 0)[0]
    return np.sort(d[back + 1, 0])


def load(preshift_mhz=PRESHIFT_MHZ):
    parts = []
    for name in FILES:
        d = read(name)
        parts.append(d[np.argsort(d[:, 0], kind="stable")])
    d = np.vstack(parts)
    d = d[np.argsort(d[:, 0], kind="stable")]
    nu = d[:, 0] - preshift_mhz / MHZ_PER_CM
    return np.column_stack([nu, 100.0 * d[:, 1]])


def init(T):
    atlas_lines.init("salami_ross", T)
    atlas_lines._state["data"] = load()


def windows(T, step=0.8, width=1.0):
    return [(float(c), width, T) for c in np.arange(14400.0 + width / 2, 14710.0 - width / 2 + 1e-9, step)]


def tscan(temps=(295.0, 325.0, 350.0, 375.0, 400.0, 450.0, 500.0), centres=(14420.3, 14470.3, 14530.3, 14580.3,
                                                                               14620.3, 14650.3, 14690.3)):
    """Window rms (and the fitted Gaussian width) against the assumed cell temperature: the Doppler width
    and the hot-band populations both depend on it."""
    atlas_lines.init("salami_ross", 300.0)
    atlas_lines._state["data"] = load()
    for T in temps:
        res = []
        for c in centres:
            rows = atlas_lines.extract((c, 1.0, T))
            res.append(rows[0][-1] if rows else np.nan)
        print(f"T = {T:5.0f} K: window rms " + " ".join(f"{100 * r:5.2f}" for r in res)
              + f"  mean {100 * np.nanmean(res):.3f}%", flush=True)


def run(T, workers):
    todo = windows(T)
    t = time.time()
    print(f"{len(todo)} windows at {T} K, {workers} workers", flush=True)
    with Pool(workers, initializer=init, initargs=(T,)) as pool, open(OUTFILE, "w") as f:
        f.write("# " + atlas_lines.HEADER + "\n")
        for n, rows in enumerate(pool.imap_unordered(atlas_lines.safe_extract, todo, chunksize=2), 1):
            for row in rows:
                f.write(ROW_FMT % tuple(row))
            f.flush()
            if n % 50 == 0:
                print(f"  {n}/{len(todo)} windows, {time.time() - t:.0f} s", flush=True)
    print(f"done in {time.time() - t:.0f} s -> {OUTFILE}")


def dataset():
    refs = reference_atlases()
    raw = np.loadtxt(OUTFILE)
    raw[:, 1] += PRESHIFT_MHZ / MHZ_PER_CM                  # back on the spectra's own axis
    lines0 = build(raw, CUTS, 1.0, 0.0, prescale_ppb=0.0)
    n0 = len(lines0)
    lines0 = drop_blends(lines0, T_CELL, *BLEND)
    print(f"{n0} lines after cuts, {len(lines0)} without blends")
    Ks, fs, _ = calibrate_sigma(lines0, refs["salami_ross_2005"])
    # Salami & Ross scatter 33 MHz about the model here, more than their own sigma, so they overstate the
    # error of this set (floor Ks, fs above). The model is known to ~1 MHz in this range (its precision
    # sets), so the scatter about it, after the constant offset, is the calibration used.
    mc0 = model_centres(lines0)
    model_ref = {k: (v / MHZ_PER_CM, MODEL_SIGMA_MHZ / MHZ_PER_CM) for k, v in mc0.items()}
    K, floor, table = calibrate_sigma(lines0, model_ref)
    print(f"sigma model: vs Salami & Ross K = {Ks:.2f}, floor = {fs:.1f} MHz; vs the model K = {K:.2f}, floor = "
          f"{floor:.1f} MHz; bins (sigma_fit, excess spread): " + ", ".join(f"({a:.1f}, {b:.1f})" for a, b in table))
    lines = {k: v for k, v in build(raw, CUTS, K, floor, prescale_ppb=0.0).items() if k in lines0}
    mc = model_centres(lines)
    keys = sorted((k for k in lines if k in mc), key=lambda k: lines[k]["nu"])
    nu = np.array([lines[k]["nu"] for k in keys])
    r = np.array([lines[k]["nu"] * MHZ_PER_CM - mc[k] for k in keys])
    s = np.array([lines[k]["sigma"] * MHZ_PER_CM for k in keys])
    med, spread = robust(r)
    print(f"\n{len(keys)} lines vs model: median {med:+.1f} MHz, robust spread {spread:.1f} MHz")
    for lo, hi in ((14400, 14500), (14500, 14600), (14600, 14710)):
        m = (nu >= lo) & (nu < hi)
        a, b = robust(r[m])
        print(f"   {lo}-{hi}: n={m.sum():4d}  median {a:+7.1f} MHz  spread {b:5.1f}")
    p, cov, keep, chi = clipped_fit(np.ones((len(r), 1)), r, s)
    rng = np.random.default_rng(1)
    blocks = np.floor(nu / 5).astype(int); ub = np.unique(blocks)
    boot = [np.median(np.concatenate([r[(blocks == b) & keep] for b in rng.choice(ub, len(ub))])) for _ in range(300)]
    offset, offset_err = float(np.median(r[keep])), float(np.std(boot))
    print(f"offset: median {offset:+.1f} MHz (block bootstrap over 5 cm-1: {offset_err:.1f}); weighted "
          f"{p[0]:+.1f} ± {np.sqrt(cov[0, 0]):.1f}; sqrt(chi2/dof) {chi:.2f}")
    sl, sc, *_ = clipped_fit(np.column_stack([np.ones_like(nu), nu - nu.mean()]), r, s)
    print(f"offset + slope: {sl[0]:+.1f} MHz at {nu.mean():.0f}, slope {sl[1]:+.2f} MHz per cm-1 "
          f"(± {np.sqrt(sc[1, 1]):.2f})")

    # per-scan offsets: do the 0.5 cm-1 scans differ by more than their lines scatter?
    edges = np.r_[14399.0, scan_edges(), 14600.2]
    lo_file = nu < 14600.2
    seg = np.searchsorted(edges, nu[lo_file]) - 1
    groups = [np.nonzero((seg == g) & keep[lo_file])[0] for g in np.unique(seg)]
    groups = [g for g in groups if len(g) >= 3]
    means = np.array([np.average(r[g], weights=1 / s[g]**2) for g in groups])
    errs = np.array([1 / np.sqrt(np.sum(1 / s[g]**2)) for g in groups])
    within = np.concatenate([(r[g] - np.average(r[g], weights=1 / s[g]**2)) / s[g] for g in groups])
    chi_between = np.sum(((means - np.average(means, weights=1 / errs**2)) / errs) ** 2) / (len(means) - 1)
    extra = np.sqrt(max(np.var(means) - np.mean(errs**2), 0.0))
    print(f"per-scan (first file, {len(groups)} scans with >= 3 lines): scan means spread {np.std(means):.1f} MHz "
          f"(robust {robust(means)[1]:.1f}), mean error {np.mean(errs):.1f}; chi2/dof between scans {chi_between:.2f};"
          f" within-scan z rms {np.sqrt(np.mean(within**2)):.2f}; extra scan-to-scan scatter {extra:.1f} MHz")
    # the same test on 0.5 cm-1 blocks at arbitrary phase, as a control, and on the second file
    for label, sel, e0 in (("first file, 0.5 cm-1 blocks", lo_file, 14400.0),
                           ("second file, 0.5 cm-1 blocks", ~lo_file, 14600.0)):
        b = np.floor((nu[sel] - e0) / 0.5).astype(int)
        idx = np.nonzero(sel)[0]
        gs = [idx[(b == g) & keep[idx]] for g in np.unique(b)]
        gs = [g for g in gs if len(g) >= 3]
        mm = np.array([np.average(r[g], weights=1 / s[g]**2) for g in gs])
        ee = np.array([1 / np.sqrt(np.sum(1 / s[g]**2)) for g in gs])
        c2 = np.sum(((mm - np.average(mm, weights=1 / ee**2)) / ee) ** 2) / (len(mm) - 1)
        print(f"   {label}: {len(gs)} blocks, means spread {np.std(mm):.1f} MHz, chi2/dof {c2:.2f}, "
              f"extra {np.sqrt(max(np.var(mm) - np.mean(ee**2), 0.0)):.1f} MHz")

    print("\nagainst the FTS atlases (offset left in):")
    compare(lines, refs, 0.0, "Rodriguez")
    for name in ("salami_ross_2005",):
        mref = model_centres([k for k in refs[name] if k in lines])
        rr = [refs[name][k][0] * MHZ_PER_CM - mref[k] for k in mref]
        print(f"  {name} - model on the same lines: n={len(rr)} median {np.median(rr):+.1f} MHz, "
              f"spread {robust(rr)[1]:.1f}")
    sr = refs["salami_ross_2005"]
    new = [k for k in lines if k not in sr]
    print(f"  lines not in salami_ross_2005: {len(new)} of {len(lines)}")
    write(lines, K, floor, offset, offset_err, med, spread, extra, fs)


def write(lines, K, floor, offset, offset_err, med, spread, extra, floor_sr):
    name = "rodriguez_fernandez_2023"
    edges = np.r_[14399.0, scan_edges()]

    def where(nu):
        """The file, and in the first file the segment between backward steps of its axis (one or more
        0.5 cm-1 scans sharing a wavemeter offset), so a fit can give each segment its own offset."""
        return f"file 1 segment {int(np.searchsorted(edges, nu))}" if nu < 14600.0 else "file 2"
    with open(ATLAS / f"{name}.csv", "w") as f:
        f.write("line,component,kind,value,uncertainty,ref_line,ref_component,group,note\n")
        for k, r in sorted(lines.items(), key=lambda kv: (kv[1]["vu"], kv[1]["vl"], kv[0])):
            f.write(f"{k},,frequency,{r['nu']:.6f},{r['sigma']:.6f},,,{name},"
                    f"depth {r['depth']:.2f}; {r['n']} window{'s' if r['n'] > 1 else ''}; {where(r['nu'])}\n")
    nus = [r["nu"] for r in lines.values()]
    (ATLAS / f"{name}.toml").write_text(
        f'id = "{name}"\nunit = "cm-1"\nisotopologue = "127I2"\ndoi = "10.1016/j.jms.2023.111789"\n'
        f'kind = "atlas line positions"\n'
        f'citation = """Transmission spectra behind D. Rodriguez Fernandez et al., J. Mol. Spectrosc. 395, 111789\n'
        f'(2023), doi:10.1016/j.jms.2023.111789 (14 400-14 600 cm-1), and M. A. Lefran Torres et al., J. Mol.\n'
        f'Spectrosc. 387, 111668 (2022), doi:10.1016/j.jms.2022.111668 (14 600-14 710 cm-1); sent by\n'
        f'D. Rodriguez Fernandez on request, 2026-09-30. The spectra are not redistributed; these line\n'
        f'centres are our fit to them, not the authors\' line lists."""\n'
        f'instrument = """Doppler-limited diode-laser transmission, 0.5 cm-1 scans in about 0.0006 cm-1 steps,\n'
        f'HighFinesse WS7 wavemeter calibrated to 85Rb 5S1/2 F=3 -> 5P3/2 F\'=4 (stated absolute uncertainty\n'
        f'0.002 cm-1, 60 MHz)."""\n'
        f'calibration = """Wavemeter. Against i2spec2026m the lines read {med:+.0f} MHz (median; robust spread\n'
        f'{spread:.0f} MHz); as a constant offset {offset:+.1f} +/- {offset_err:.1f} MHz (block bootstrap over 5 cm-1).\n'
        f'Subtract it to put the values on the model\'s scale; it is NOT in the per-line uncertainty and\n'
        f'belongs in the fit as a group parameter. In the first file (14 400-14 600 cm-1) the scans also differ\n'
        f'from each other: {extra:.0f} MHz of offset scatter between the segments delimited by the backward steps\n'
        f'of its axis (chi2/dof between segments 6.5), none in the second file. That scatter is in the per-line\n'
        f'floor, not corrected; the note column names each line\'s segment for a fit that wants per-segment\n'
        f'offsets."""\n'
        f'extraction = """prototypes/rodriguez_lines.py: both files sorted and joined; axis shifted by\n'
        f'-{PRESHIFT_MHZ:.0f} MHz before the fits and restored after. prototypes/atlas_lines.py extract at\n'
        f'T = {T_CELL:.0f} K (from the window rms, --tscan) in 1.0 cm-1 windows stepped by 0.8 cm-1: per-line\n'
        f'wavenumber shifts of the i2spec template, with column density, linear baseline, zero offset and a\n'
        f'Gaussian width per window, hyperfine structure in the template, ridge prior 0.01 cm-1. Repeats\n'
        f'combined by inverse variance. Kept: fit sigma < {CUTS["sigma"]:.0f} MHz, window rms < {100 * CUTS["rms"]:.0f}%,\n'
        f'depth > {CUTS["depth"]:.2f}, no other model line within +/-{BLEND[0]} cm-1 at >= {BLEND[1]} of its strength.\nAssignments are the template\'s; blends fail the sigma cut."""\n'
        f'uncertainty = """sqrt(({K:.2f} sigma_fit)^2 + {floor:.1f}^2) MHz, fitted to the scatter of the lines about\n'
        f'i2spec2026m (after the constant offset; the model is known to ~1 MHz here from its precision sets) in\n'
        f'bins of sigma_fit. The floor holds the template error and the scan-to-scan wavemeter offsets. Against\n'
        f'salami_ross_2005 the same calibration gives {floor_sr:.0f} MHz, but Salami & Ross themselves scatter 33 MHz\n'
        f'about the model in this range. Excludes the offset."""\n'
        f'offset_mhz = {offset:.1f}\noffset_mhz_uncertainty = {offset_err:.1f}\n'
        f'range_cm-1 = [{min(nus):.1f}, {max(nus):.1f}]\nn_lines = {len(lines)}\n')
    print(f"wrote {name}: {len(lines)} lines")


def main(argv):
    opts = {a.split("=")[0]: a.split("=", 1)[1] for a in argv if "=" in a}
    if "--tscan" in argv:
        return tscan()
    if "--dataset" in argv:
        return dataset()
    run(float(opts.get("--T", T_CELL)), int(opts.get("--workers", 2)))


if __name__ == "__main__":
    main(sys.argv[1:])
