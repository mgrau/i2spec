"""Condense the spectra we may not redistribute into line positions and intensities.

The raw spectra stay in data/external/ (not in the repository). What this writes goes to the private
companion repository i2spec_data ($I2SPEC_DATA, default ../i2spec_data), until the terms of each source
allow more.

Each spectrum is fitted window by window with prototypes/atlas_lines.py extract, amplitudes on: every line
deep enough to matter gets a wavenumber shift and a strength factor on the i2spec template (hyperfine
structure, Doppler width at the cell temperature, a Gaussian instrument width, column density, linear
baseline and zero offset per window). Per line it keeps

  position   the hyperfine-free line centre on the spectrum's own wavenumber scale (cm-1), with
             sigma = sqrt((K sigma_fit)^2 + floor^2), K and floor fitted to the scatter against a
             reference (another atlas, or the model where it is known to ~1 MHz) after a smooth scale curve;
  intensity  the integrated optical depth A = N S (cm-1): column density times line strength, the area
             under -ln(transmission) of that line before the instrument function. It depends on the model
             only through the template's line shape and the split of a blend. A / N_window is the strength
             on the scale of the window's fitted column density.

Repeats (overlapping windows) are combined by inverse variance, intensities per cell temperature.

Per source, in the private repository:
  lines/<id>.csv               positions, in the observations format (i2spec.observations.read_observations)
  lines/<id>.intensities.csv   line, T_K, A_cm-1, A_unc_cm-1, N_window_cm-2, depth, windows, blend
  lines/<id>.toml              provenance: source, conditions, terms, method, uncertainty calibration

usage: condense_spectra.py --list
       condense_spectra.py --calibrate ID ...     scan scale (and temperature where unknown), Kitt Peak only
       condense_spectra.py --fit ID ... [--workers=6]          window fits -> prototypes/out/condense/ID.txt
       condense_spectra.py --build ID ...                       the files above, and a report
       condense_spectra.py --index                              lines/README.md, a table of what is there
       condense_spectra.py --coadd                              co-add the Debus et al. 2023 spectra (before the rest)
       ID may also be 'all' or 'kitt_peak'
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
import time
from collections import defaultdict
from multiprocessing import Pool
from pathlib import Path

import numpy as np

os.environ.setdefault("I2SPEC_WORKERS", "2")          # keep any line-list rebuild modest
sys.path.insert(0, str(Path(__file__).resolve().parent))
import atlas_lines                                                           # noqa: E402
from atlas_dataset import ATLASES as PUBLISHED, combine                     # noqa: E402
from kitt_peak_lines import calibrate_sigma, curve, drop_blends, envelope, model_centres, reference_atlases, robust  # noqa: E402
from scan_salami_ross import SEGMENTS                                       # noqa: E402
from i2spec.constants import MHZ_PER_CM                                     # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
RAW_KP = ROOT / "data" / "external" / "kitt_peak_fts"
OUT = ROOT / "prototypes" / "out" / "condense"
CALIBRATION = OUT / "calibration"          # one <id>.json per source, so calibrations can run in parallel
DEST = Path(os.environ.get("I2SPEC_DATA", ROOT.parent / "i2spec_data")) / "lines"
#: sources whose terms allow it (public=True) go to the repository itself instead
PUBLIC = ROOT / "data" / "atlas_lines"
CUTS = dict(sigma=30.0, rms=0.04, depth=0.03, sigma_eta=0.25)
#: kept from the published extraction, where a source has too few lines in common with any reference
K_DEFAULT, FLOOR_DEFAULT = 1.55, 7.5
ROW_FMT = "%.6f %.7f %.7f %.4f %.4g %d %d %d %d %.2f %.3f %.5f %.5f %.5f %.5f\n"

SOURCES = {}


def source(id, **kw):
    SOURCES[id] = dict(dict(id=id, cuts=dict(CUTS), blend=(0.01, 0.1), drop_blends=False, reference="apo_nist_2009",
                            prescale_ppb=0.0, preshift_mhz=0.0), **kw)


source("salami_ross_2005", kind="salami_ross", title="Salami & Ross (2005) FTS atlas, 50 cm cell", jobs=lambda: list(atlas_lines.jobs(2.0, 2.4, "salami_ross")),
       reference="apo_nist_2009", url="https://doi.org/10.1016/j.jms.2005.06.002",
       citation=PUBLISHED["salami_ross_2005"]["citation"], instrument=PUBLISHED["salami_ross_2005"]["instrument"],
       conditions="50 cm cell; four segments with the cell body at 20, 20, 50 and 190 C (Salami & Ross 2005, Table 1): "
                  "19 100-20 500, 16 500-19 500, 15 400-17 400 and 14 250-16 500 cm-1",
       terms="Elsevier supplementary file mmc1 of the article; no licence stated",
       intensity_note="Against the 1993 and 1995 Kitt Peak cells at the same temperature (323 K), the "
                      "intensities of the 50 C segment scatter 17-22 % per line about a constant ratio, about three "
                      "times the stated sigma, while those cells agree with each other to 3-6 %: treat the stated "
                      "sigma of this source's intensities as a lower limit")
source("apo_nist_2009", kind="apo", title="NIST 2-m FTS scan of the APO cell (2009)", jobs=lambda: list(atlas_lines.jobs(2.0, 2.4, "apo", 333.15)),
       reference="salami_ross_2005", url="https://github.com/jaklusmeyer/APO_pyodine",
       citation=PUBLISHED["apo_nist_2009"]["citation"], instrument=PUBLISHED["apo_nist_2009"]["instrument"],
       conditions="APO cell, sample AS-17; temperature and path not recorded. The template is at 333.15 K, as in "
                  "the published extraction",
       terms="no licence file in the APO_pyodine repository: research use only, not redistributed")
source("rodriguez_fernandez_2023", kind="rodriguez", title="Rodríguez Fernández / Lefrán Torres diode-laser scans (2022-23)", reference="model", blend=(0.01, 0.3), drop_blends=True,
       preshift_mhz=131.0, T=350.0, url="https://doi.org/10.1016/j.jms.2023.111789",
       jobs=lambda: [(float(c), 1.0, 350.0) for c in np.arange(14400.5, 14709.5 + 1e-9, 0.8)],
       citation="Transmission spectra behind D. Rodriguez Fernandez et al., J. Mol. Spectrosc. 395, 111789 (2023), "
                "doi:10.1016/j.jms.2023.111789 (14 400-14 600 cm-1), and M. A. Lefran Torres et al., J. Mol. "
                "Spectrosc. 387, 111668 (2022), doi:10.1016/j.jms.2022.111668 (14 600-14 710 cm-1); sent by "
                "D. Rodriguez Fernandez on request, 2026-09-30",
       instrument="Doppler-limited diode-laser transmission, 0.5 cm-1 scans in about 0.0006 cm-1 steps, HighFinesse "
                  "WS7 wavemeter calibrated to 85Rb (stated absolute uncertainty 0.002 cm-1, 60 MHz)",
       conditions="room-temperature cell; the template is at 350 K, the lowest window rms (rodriguez_lines.py --tscan)",
       terms="shared by the authors for this work; not for redistribution")

#: Kitt Peak McMath-Pierce 1-m FTS laboratory iodine-cell spectra (NSO archive): file -> (cell temperature
#: in C or None to fit it, band in cm-1, description). docs/research/kitt-peak-and-rodriguez-spectra.md lists
#: all 54; the others are sunlight or a laser through the cell, emission or calibration runs, other gases,
#: or coarser than 0.04 cm-1. The June 1995 long and short cells (950626R0.001, .002) are left out too: a ~10 %
#: ripple runs through both, and their lines are shallower and broader than the header's 0.020 cm-1 allows.
#: So are the November 1990 scans (901126R0.001, .002): their lines are far broader than the 1993 scans' at the
#: same nominal 0.037 cm-1, and no position comes out better than 300 MHz. The April 1996 scans of a 10 cm cell
#: behind 2.4 m of O2 (960405R0.016, .017, .022) fit well but give no position better than 30 MHz at 0.030-0.035
#: cm-1 and strengths to only ~30 %: kept out as adding nothing the 1993 scans do not.
#: the fitted column densities of 12-14 October 1988 are 0.94-1.00 of the vapour pressure over 25 cm at the logged
#: temperatures (3.7-23.6 C); those of 15 October agree with each other but are 0.85-0.875 of it at the logged 23.9 C
OCT15 = "25 cm cell, logged at 23.9 C"
NOTES = {n: ("All five scans of 15 October fit a column density 0.85-0.875 of the vapour pressure at the logged "
             "23.9 C over 25 cm (as if the cold point were 22.2 C), where the scans of 12-14 October agree with their "
             "logged temperatures to 0.94-1.00. ") for n in ("881015R0.014", "881015R0.015", "881015R0.016",
                                                            "881015R0.017", "881015R0.018")}
KITT_PEAK = {
    "840531R0.001": (None, None, "1.5 inch cell, DZE lamp, 10 scans"),
    "840531R0.002": (None, None, "1.5 inch cell, DZE lamp, 40 scans"),
    "840611R0.002": (None, None, "1.5 inch cell, quartz lamp, 20 scans"),
    "881012R0.002": (23.6, None, "25 cm cell at 23.6 C"),
    "881012R0.003": (3.95, None, "25 cm cell at 3.6-4.3 C"),
    "881012R0.004": (3.65, None, "25 cm cell at 3.5-3.8 C, narrow band"),
    "881013R0.006": (14.05, None, "25 cm cell at 13.9-14.2 C"),
    "881014R0.007": (7.9, None, "25 cm cell at 7.8-8.0 C"),
    "881014R0.008": (7.9, None, "25 cm cell at 7.8-8.0 C"),
    "881014R0.009": (7.9, None, "25 cm cell at 7.8-8.0 C"),
    "881014R0.010": (7.8, None, "25 cm cell at 7.7-7.9 C"),
    "881014R0.011": (16.25, None, "25 cm cell at 16.1-16.4 C"),
    "881014R0.012": (16.25, None, "25 cm cell at 16.1-16.4 C"),
    "881014R0.013": (16.25, None, "25 cm cell at 16.1-16.4 C"),
    "881015R0.014": (23.9, None, OCT15),
    "881015R0.015": (23.9, None, OCT15),
    "881015R0.016": (23.9, None, OCT15),
    "881015R0.017": (23.9, None, OCT15),
    "881015R0.018": (23.9, None, OCT15),
    "930317R0.001": (50.0, None, "Lick backup cell, 100 mm, 50.0 C"),
    "930317R0.002": (50.1, None, "McDonald cell 2, 160 mm, 50.1 C"),
    "930317R0.003": (50.0, None, "Keck cell, 100 mm, 50.0 C"),
    "930317R0.004": (60.0, None, "Keck cell, 100 mm, 60.0 C"),
    "930317R0.005": (70.0, None, "Keck cell, 100 mm, 70.0 C"),
    "930317R0.006": (50.1, None, "McDonald cell 1, 150 mm, 50.1 C"),
    "930317R0.007": (50.0, None, "cell 4A, 100 mm, 50.0 C"),
    "930317R0.008": (50.0, None, "ESO cell, 100 mm, 50.0 C"),
    "930317R0.009": (70.0, None, "ESO cell, 100 mm, 70.0 C"),
    "930317R0.011": (50.0, None, "cell 1A, 100 mm, 50.0 C"),
    "930317R0.013": (50.0, None, "cell 3A, 100 mm, 50.0 C"),
    "951009R0.001": (50.0, None, "Lick cell, 10 cm, 50 C"),
    "951009R0.002": (50.0, None, "Lick cell, 10 cm, 50 C"),
}


def kp_id(name):
    return "kitt_peak_" + name.replace("R0.", "_")


def read_fits(path):
    """(header, wavenumber, signal) of an NSO FTS FITS file (BITPIX -32; x = WSTART + i DELW whatever
    XAXIS_IS says; byte order from BOCODE)."""
    raw = Path(path).read_bytes()
    hdr, i = {}, 0
    while True:
        card = raw[i:i + 80].decode("ascii", "replace"); i += 80
        key = card[:8].strip()
        if card[8:10] == "= ":
            hdr[key] = card[10:].split("/")[0].strip().strip("'").strip()
        if key == "END":
            break
    assert hdr["BITPIX"] == "-32" and hdr["DATA_IS"].startswith("Real")
    off, n = -(-i // 2880) * 2880, int(hdr["NPO"])
    y = np.frombuffer(raw[off:off + 4 * n], dtype=">f4" if hdr.get("BOCODE", "1") == "1" else "<f4").astype(float)
    nu = float(hdr["WSTART"]) + np.arange(n) * float(hdr["DELW"])
    # WSTOP is the last point (1993) or one step past it (the older files; their CRVAL1/CDELT1 agree with
    # WSTART + i DELW): either way the axis is the same
    assert min(abs(nu[-1] - float(hdr["WSTOP"])), abs(nu[-1] + float(hdr["DELW"]) - float(hdr["WSTOP"]))) < 1e-3
    return hdr, nu, y


def kp_band(name, hdr):
    band = KITT_PEAK[name][1]
    lo, hi = band or (float(hdr.get("BANDLO") or 0), float(hdr.get("BANDHI") or 0))
    if hi <= lo:
        lo, hi = float(hdr["WSTART"]), float(hdr["WSTOP"])
    return max(lo, float(hdr["WSTART"]) + 2, 14250.0), min(hi, float(hdr["WSTOP"]) - 2, 20050.0)


def strong_band(nu, y, lo, hi, thr=0.2):
    """The longest stretch of lo-hi where the lamp-times-filter signal (upper envelope in 10 cm-1 blocks) is at
    least thr of its 95th percentile: beyond it the normalisation, and with it the intensities, fail (pairs of
    narrow-band 1988 scans overlapping only there disagree by 20-25 % in A)."""
    m = (nu >= lo) & (nu <= hi)
    nu, y = nu[m], y[m]
    edges = np.arange(lo, hi + 10.0, 10.0)
    c, v = [], []
    for a, b in zip(edges[:-1], edges[1:]):
        k = (nu >= a) & (nu < b)
        if k.sum() > 50:
            c.append(0.5 * (a + b)); v.append(np.percentile(y[k], 99.5))
    c, v = np.array(c), np.array(v)
    ok = v >= thr * np.percentile(v, 95)
    best, start = (0, 0), None
    for i, flag in enumerate(np.append(ok, False)):
        if flag and start is None:
            start = i
        elif not flag and start is not None:
            if i - start > best[1] - best[0]:
                best = (start, i)
            start = None
    return max(lo, c[best[0]] - 5.0), min(hi, c[best[1] - 1] + 5.0)


def usable_band(name, hdr, thr=0.2):
    """strong_band of a Kitt Peak scan within kp_band."""
    _, nu, y = read_fits(RAW_KP / f"{name}.fits")
    return strong_band(nu, y, *kp_band(name, hdr), thr)


for _name, (_T, _band, _desc) in KITT_PEAK.items():
    source(kp_id(_name), kind="kitt_peak", file=_name, T=None if _T is None else round(_T + 273.15, 2), description=_desc,
           blend=None, drop_blends=True, url="https://nispdata.nso.edu/ftp/FTS_cdrom/",
           cuts=dict(CUTS, depth=0.05), terms="NSO/NOAO archive data; public, no licence stated")


#: Debus, Schäfer & Reiners 2023: 222 spectra of a Thorlabs GC19100-I cell (halogen lamp) taken with a Fabry-Perot
#: etalon over about 90 h from 20 January 2023 on the IAG Göttingen IFS 125HR at 0.02 cm-1. GRO.data
#: doi:10.25625/PONFBL, CC0. The co-added single-channel spectrum is the input (--coadd); the iodine channel
#: covers 498-659 nm, less +/-3 nm about the 633 nm HeNe artifact (the paper), which the windows skip.
RAW_DEBUS = ROOT / "data" / "external" / "debus2023"
DEBUS_COADD = OUT / "debus_2023_coadd.npz"
DEBUS_RANGE = (15180.0, 20050.0)
HENE_GAP = (15720.0, 15876.0)
source("debus_2023", kind="debus", title="IAG Göttingen FTS, January 2023: Thorlabs iodine cell (Debus et al. 2023)",
       T=None, blend=None, drop_blends=True, url="https://doi.org/10.25625/PONFBL", cuts=dict(CUTS, depth=0.05),
       band_threshold=0.05,      # the halogen lamp falls smoothly by 10x to the blue, no filter edge; the co-add holds up
       public=True,
       citation="M. Debus, S. Schäfer, A. Reiners, Toward 10 cm/s radial velocity accuracy on the Sun using a Fourier "
                "transform spectrometer, J. Astron. Telesc. Instrum. Syst. 9, 045003 (2023), doi:10.1117/1.JATIS.9.4.045003; "
                "replication data on GRO.data, doi:10.25625/PONFBL (FP_iodine_spectra.7z)",
       terms="CC0 1.0 (GRO.data): no restriction on use or redistribution; cite the paper and the data set")


#: The two FTS iodine-cell spectra distributed with the pyodine radial-velocity code (gitlab.com/Heeren/pyodine,
#: iodine_atlas/, MIT licence): continuum-normalised transmission against vacuum wavenumber. Their scale was
#: resampled for RV use (about -2e-7 against Salami & Ross), which --calibrate measures and the fits take out.
RAW_PYODINE = ROOT / "data" / "external" / "heeren_pyodine_atlas"
PYODINE = {
    "pyodine_song_65c": ("song_iodine_cell_01_65C.h5", 338.15, "SONG iodine cell at 65 C",
                         "the SONG telescope's (Teide) iodine cell, 65 C"),
    "pyodine_fischer_2022": ("Fischer_Cell_May2022_downsampled3.h5", None, "D. Fischer's iodine cell, May 2022",
                             "D. Fischer's (Lick/Yale) iodine cell, scanned May 2022, downsampled 3x by the pyodine authors"),
}
for _id, (_file, _T, _title, _desc) in PYODINE.items():
    source(_id, kind="pyodine", file=_file, T=_T, title=f"{_title} (pyodine FTS atlas)", description=_desc, public=True,
           blend=(0.01, 0.1), drop_blends=True, cuts=dict(CUTS, depth=0.05), url="https://gitlab.com/Heeren/pyodine",
           citation="FTS iodine-cell spectrum distributed in iodine_atlas/ of the pyodine code, https://gitlab.com/Heeren/pyodine "
                    "(file " + _file + "); P. Heeren, R. J. Tronsgaard, F. Grundahl, S. Reffert, A. Quirrenbach, T. Marquart, "
                    "Astron. Astrophys. 674, A164 (2023), doi:10.1051/0004-6361/202244441",
           instrument="Fourier transform spectrum of the cell, continuum-normalised by the pyodine authors; instrument and "
                      "resolution are not recorded in the file",
           terms="MIT (the pyodine repository's licence); cite Heeren et al. 2023")


def pyodine_read(name):
    """(ascending vacuum wavenumber, normalised transmission, raw flux) of a pyodine atlas file."""
    import h5py
    with h5py.File(RAW_PYODINE / name) as h:
        nu, fl, raw = h["wavenumber"][...], h["flux_normalized"][...].astype(float), h["flux"][...].astype(float)
    o = np.argsort(nu)
    ok = np.isfinite(nu[o]) & np.isfinite(fl[o])
    return nu[o][ok], fl[o][ok], raw[o][ok]


def debus_files():
    return sorted(p for p in RAW_DEBUS.rglob("*.0") if "Iodine" in p.name)


def debus_coadd():
    """Sum the single-channel (ScSm) spectra of every Debus file on the first file's axis (they share it; a file
    that does not is interpolated) and keep 14 500-20 500 cm-1."""
    from opus import read
    files = debus_files()
    if not files:
        raise SystemExit(f"no spectra under {RAW_DEBUS}")
    xf, total, status, params = read(files[0])
    m = (xf >= 14500) & (xf <= 20500)
    x0, total = xf[m], total[m].copy()
    times, regridded = [status.get("DAT", "") + " " + status.get("TIM", "")], 0
    for f in files[1:]:
        x, y, st, _ = read(f)
        if len(x) == len(xf) and abs(x[0] - xf[0]) < 1e-9 and abs(x[-1] - xf[-1]) < 1e-9:
            total += y[m]
        else:
            total += np.interp(x0, x, y); regridded += 1
        times.append(st.get("DAT", "") + " " + st.get("TIM", ""))
    info = {f"{k[0]}_{k[1]}_{k[2]}": v for k, v in params.items()}
    OUT.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(DEBUS_COADD, nu=x0, y=total, n=len(files), first=times[0], last=times[-1],
                        params=json.dumps(info, default=str))
    print(f"{len(files)} spectra co-added ({regridded} regridded), {times[0]} to {times[-1]}, {len(x0)} points, "
          f"step {(x0[-1] - x0[0]) / (len(x0) - 1):.6f} cm-1 -> {DEBUS_COADD}")


def debus_params():
    d = np.load(DEBUS_COADD)
    return d, json.loads(str(d["params"]))


# --- loading and window fits -----------------------------------------------------------------------

def calibration():
    return {p.stem: json.loads(p.read_text()) for p in sorted(CALIBRATION.glob("*.json"))}


def settings(src):
    """The source with its scale and temperature from --calibrate filled in (Kitt Peak)."""
    s = dict(src)
    if s["kind"] == "kitt_peak":
        hdr = read_fits(RAW_KP / f"{s['file']}.fits")[0]
        s["header"] = hdr
        s["resolution"] = float(hdr.get("RESOLUTN") or 0)
        s["band"] = usable_band(s["file"], hdr)
        s["width"] = max(2.4, round(230 * float(hdr["DELW"]), 1))
        s["step"] = round(0.8 * s["width"], 1)
        if s["blend"] is None:
            s["blend"] = (max(s["resolution"] / 2, 0.01), 0.1)
        s.update(calibration().get(s["id"], {}))
        if s.get("T") is not None:
            lo, hi = s["band"]
            s["jobs"] = (lambda lo=lo, hi=hi, w=s["width"], st=s["step"], T=s["T"]:
                         [(float(c), w, T) for c in np.arange(np.ceil(lo) + w / 2, hi - w / 2, st)])
    elif s["kind"] == "pyodine":
        nu, fl, raw = pyodine_read(s["file"])
        step = (nu[-1] - nu[0]) / (len(nu) - 1)
        s["resolution"] = 0.0
        s["band"] = strong_band(nu, raw, max(nu[0] + 2, 15400.0), min(nu[-1] - 2, 20050.0), 0.1)
        s["width"] = max(2.4, round(230 * step, 1))
        s["step"] = round(0.8 * s["width"], 1)
        s.update(calibration().get(s["id"], {}))
        if s.get("T") is not None:
            lo, hi = s["band"]
            s["jobs"] = (lambda lo=lo, hi=hi, w=s["width"], st=s["step"], T=s["T"]:
                         [(float(c), w, T) for c in np.arange(np.ceil(lo) + w / 2, hi - w / 2, st)])
    elif s["kind"] == "debus":
        d, _ = debus_params()
        nu, y = d["nu"], d["y"]
        step = (nu[-1] - nu[0]) / (len(nu) - 1)
        s["resolution"] = 0.02
        s["band"] = strong_band(nu, y, *DEBUS_RANGE, s["band_threshold"])
        s["width"] = max(2.4, round(230 * step, 1))
        s["step"] = round(0.8 * s["width"], 1)
        s["blend"] = (0.01, 0.1)
        s.update(calibration().get(s["id"], {}))
        if s.get("T") is not None:
            lo, hi = s["band"]
            gap = lambda c, w: c + w / 2 > HENE_GAP[0] and c - w / 2 < HENE_GAP[1]          # noqa: E731
            s["jobs"] = (lambda lo=lo, hi=hi, w=s["width"], st=s["step"], T=s["T"]:
                         [(float(c), w, T) for c in np.arange(np.ceil(lo) + w / 2, hi - w / 2, st) if not gap(c, w)])
    return s


def load(s):
    """(nu, per cent transmission) on the axis the window fits see: the spectrum's own axis divided by
    (1 + prescale) and shifted by -preshift, so the 0.01 cm-1 ridge prior on each shift does not pull
    against a calibration offset; build() puts both back."""
    if s["kind"] == "kitt_peak":
        hdr, nu, y = read_fits(RAW_KP / f"{s['file']}.fits")
        nu = nu / (1 + s["prescale_ppb"] * 1e-9)
        lo, hi = s["band"]
        m = (nu >= lo) & (nu <= hi)
        nu, y = nu[m], y[m]
        return np.column_stack([nu, 100.0 * y / envelope(nu, y)])
    if s["kind"] == "debus":
        d, _ = debus_params()
        nu, y = d["nu"] / (1 + s["prescale_ppb"] * 1e-9), d["y"]
        lo, hi = s["band"]
        m = (nu >= lo) & (nu <= hi) & ~((nu > HENE_GAP[0]) & (nu < HENE_GAP[1]))
        nu, y = nu[m], y[m]
        return np.column_stack([nu, 100.0 * y / envelope(nu, y)])
    if s["kind"] == "pyodine":
        nu, fl, _ = pyodine_read(s["file"])
        nu = nu / (1 + s["prescale_ppb"] * 1e-9)
        lo, hi = s["band"]
        m = (nu >= lo) & (nu <= hi)
        return np.column_stack([nu[m], 100.0 * fl[m]])
    if s["kind"] == "rodriguez":
        import rodriguez_lines
        return rodriguez_lines.load(s["preshift_mhz"])
    return None                                                     # atlas_lines loads its own


def init(sid):
    s = settings(SOURCES[sid])
    if s["kind"] == "apo":
        atlas_lines.init("apo", 333.15, amplitudes=True)
    else:
        atlas_lines.init("salami_ross", s.get("T"), amplitudes=True)
        data = load(s)
        if data is not None:
            atlas_lines._state["data"] = data


def fit(sid, workers):
    s = settings(SOURCES[sid])
    if s["kind"] == "kitt_peak" and s.get("T") is None:
        raise SystemExit(f"{sid}: no temperature yet; run --calibrate {sid}")
    todo = s["jobs"]()
    OUT.mkdir(parents=True, exist_ok=True)
    out = OUT / f"{sid}.txt"
    t = time.time()
    print(f"{sid}: {len(todo)} windows, {workers} workers", flush=True)
    with Pool(workers, initializer=init, initargs=(sid,)) as pool, open(out, "w") as f:
        f.write("# " + atlas_lines.HEADER_AMPLITUDES + "\n")
        for n, rows in enumerate(pool.imap_unordered(atlas_lines.safe_extract, todo, chunksize=2), 1):
            for row in rows:
                f.write(ROW_FMT % tuple(row))
            if n % 200 == 0:
                f.flush()
                print(f"  {sid}: {n}/{len(todo)} windows, {time.time() - t:.0f} s", flush=True)
    print(f"{sid}: done in {time.time() - t:.0f} s -> {out}", flush=True)


def xcorr_prescale(s, half=20.0):
    """A model-free first scale (ppb): the lag that best aligns the spectrum with the Salami & Ross atlas in
    40 cm-1 stretches across the band, median over stretches. Scales of several 1000 ppb (the 1990 scans) are
    beyond what the per-line shifts can take up."""
    from scan_salami_ross import ATLAS
    sr = np.loadtxt(ATLAS, skiprows=5)
    d = load(dict(s, prescale_ppb=0.0))
    lo, hi = max(s["band"][0], 14400.0) + half, min(s["band"][1], 20000.0) - half
    lags, ppb = np.arange(-0.3, 0.3, 0.001), []
    for c in np.linspace(lo, hi, 6):
        grid = np.arange(c - half, c + half, 0.002)
        a = np.interp(grid, sr[:, 0], sr[:, 1]); a = a - a.mean()
        cc = []
        for L in lags:
            b = np.interp(grid + L, d[:, 0], d[:, 1])
            cc.append(np.dot(a, b - b.mean()))
        ppb.append(lags[int(np.argmax(cc))] / c * 1e9)          # obs = true (1 + s): a feature at x sits at x (1 + s)
    print(f"{s['id']}: cross-correlation with Salami & Ross: " + ", ".join(f"{x:+.0f}" for x in ppb) + " ppb", flush=True)
    return float(np.median(ppb))


def calibrate(sid, n=36):
    """Kitt Peak: the scale (ppb) from the median shift of well-measured lines in n windows, iterated, and
    the cell temperature, where unknown, from the lowest mean window rms."""
    s = settings(SOURCES[sid])
    lo, hi = s["band"]
    centres = np.linspace(lo + s["width"], hi - s["width"], n)
    cal = calibration()
    entry = dict(cal.get(sid, {}))
    known_T = SOURCES[sid].get("T")                           # from the header; otherwise fitted here
    T = known_T or 296.0
    atlas_lines.init("salami_ross", T, amplitudes=False)     # the template temperature comes with each window

    def pass_(T, prescale):
        if s.get("_loaded") != prescale:
            s["prescale_ppb"] = prescale
            atlas_lines._state["data"] = load(s)
            s["_loaded"] = prescale
        rows = [r for c in centres for r in atlas_lines.safe_extract((float(c), s["width"], T))]
        return np.array(rows) if rows else np.zeros((0, 12))

    def scale_passes(T, prescale, sigma_mhz, rms, passes):
        for it in range(passes):
            rows = pass_(T, prescale)
            good = (rows[:, 2] * MHZ_PER_CM < sigma_mhz) & (rows[:, 3] > 0.15) & (rows[:, 11] < rms)
            if good.sum() < 20:
                print(f"{sid}: only {good.sum()} lines with sigma < {sigma_mhz} MHz at {T:.1f} K; scale left at "
                      f"{prescale:+.0f} ppb", flush=True)
                break
            step = float(np.median(rows[good, 1] / rows[good, 0])) * 1e9
            prescale += step
            print(f"{sid}: pass at {T:.1f} K: {good.sum()} lines, scale {prescale:+.0f} ppb (step {step:+.0f})", flush=True)
            if abs(step) < 5:
                break
        return prescale

    prescale = entry.get("prescale_ppb") or xcorr_prescale(s)
    # loose passes first: an offset of a few 100 MHz inflates every line's sigma past the strict cut
    prescale = scale_passes(T, prescale, 150.0, 0.10, 2)
    if known_T is None:
        # the temperature from the lowest mean window rms (Doppler width and hot-band populations both move it)
        temps = (280.0, 296.0, 310.0, 323.0, 340.0, 360.0, 390.0)
        rms = []
        for Tc in temps:
            win = {row[10]: row[11] for row in pass_(Tc, prescale)}
            rms.append(np.mean(list(win.values())) if win else np.inf)
            print(f"   T = {Tc:.0f} K: mean window rms {100 * rms[-1]:.3f} %", flush=True)
        i = int(np.argmin(rms))
        T = temps[i]
        if 0 < i < len(temps) - 1:                        # parabola through the minimum and its neighbours
            c = np.polyfit(np.array(temps[i - 1:i + 2]), np.array(rms[i - 1:i + 2]), 2)
            if c[0] > 0:
                T = float(-c[1] / (2 * c[0]))
        entry["T_fitted"] = True
        print(f"{sid}: template temperature {T:.1f} K", flush=True)
    prescale = scale_passes(T, prescale, 60.0, 0.04, 4)
    entry.update(prescale_ppb=round(prescale, 1), T=round(T, 2))
    CALIBRATION.mkdir(parents=True, exist_ok=True)
    (CALIBRATION / f"{sid}.json").write_text(json.dumps(entry, indent=1, sort_keys=True) + "\n")


# --- the data set ----------------------------------------------------------------------------------

def key(vu, vl, J, br):
    return f"127I2 {'R' if br > 0 else 'P'}({J}) {vu}-{vl}"


def references(s):
    if s["reference"] == "model":
        return None
    return reference_atlases()[s["reference"]]


def build(sid):
    s = settings(SOURCES[sid])
    raw = np.loadtxt(OUT / f"{sid}.txt")
    if "band" in s:
        raw = raw[(raw[:, 0] >= s["band"][0]) & (raw[:, 0] <= s["band"][1])]
    scale = 1 + s["prescale_ppb"] * 1e-9
    shift = s["preshift_mhz"] / MHZ_PER_CM
    cuts = s["cuts"]

    # positions: repeats of a line combined (both temperature segments of Salami & Ross count as repeats)
    keep = (raw[:, 2] * MHZ_PER_CM < cuts["sigma"]) & (raw[:, 11] < cuts["rms"]) & (raw[:, 3] > cuts["depth"])
    groups = defaultdict(list)
    for row in raw[keep]:
        groups[(int(row[5]), int(row[6]), int(row[7]), int(row[8]))].append(row)
    lines = {}
    for (vu, vl, J, br), rows in groups.items():
        rows = np.array(rows)
        delta, sigma = combine(rows)
        lines[key(vu, vl, J, br)] = dict(nu=(float(rows[0, 0]) + delta) * scale + shift, sigma_fit=sigma, n=len(rows),
                                         depth=float(rows[:, 3].max()), vu=vu, vl=vl, T=float(rows[0, 9]))
    # the blend flag, for every line that has a position or an intensity
    if not lines:
        raise ValueError(f"no line passes the cuts ({len(raw)} window-lines fitted)")
    T_ref = float(np.median(raw[:, 9]))
    every = {key(int(r[5]), int(r[6]), int(r[7]), int(r[8])): dict(nu=float(r[0]), vu=int(r[5]), vl=int(r[6]))
             for r in raw}
    clean = drop_blends(every, T_ref, *s["blend"])
    blend = {k: k not in clean for k in every}
    if s["drop_blends"]:
        lines = {k: v for k, v in lines.items() if not blend[k]}

    # the uncertainty: K and floor from the scatter against the reference, after a quadratic scale curve
    ref = references(s)
    if ref is None:
        mc = model_centres(lines)
        ref = {k: (v / MHZ_PER_CM, 1.0 / MHZ_PER_CM) for k, v in mc.items()}
        ref_name = "the i2spec model (known to ~1 MHz in this range)"
    else:
        ref_name = s["reference"]
    both = [k for k in lines if k in ref]
    trend, scale_curve = None, None
    if len(both) >= 300:
        nu_b = np.array([lines[k]["nu"] for k in both])
        d_b = np.array([(lines[k]["nu"] - ref[k][0]) * MHZ_PER_CM for k in both])
        span = nu_b.max() - nu_b.min()
        deg = 2 if span > 2000 else 1 if span > 500 else 0      # a curve the range can determine
        p, e, m, chi, trend = curve(nu_b, d_b, np.ones_like(d_b), deg)
        scale_curve = (p, e)
        K, floor, table = calibrate_sigma(lines, ref, trend=trend)
        how = (f"fitted to the scatter of the positions about {ref_name} ({len(both)} common lines), after a "
               f"degree-{deg} scale curve, in bins of sigma_fit (K = {K:.2f}, floor = {floor:.1f} MHz)")
        if floor < FLOOR_DEFAULT:
            # errors of the template common to both spectra cancel in their difference, so the fitted floor is
            # a lower limit; the published extraction's floor was checked against the precision measurements
            how += (f"; the floor is raised to {FLOOR_DEFAULT} MHz, the published extraction's (checked against "
                    f"precision measurements), because template errors common to both spectra cancel here")
            floor = FLOOR_DEFAULT
    else:
        K, floor, table = K_DEFAULT, FLOOR_DEFAULT, []
        how = (f"the published extraction's values: only {len(both)} lines in common with {ref_name}, too few "
               f"to fit them")
    for r in lines.values():
        r["sigma"] = float(np.hypot(K * r["sigma_fit"] * MHZ_PER_CM, floor)) / MHZ_PER_CM

    # intensities: per (line, cell temperature); ln A combined by inverse variance
    ok = (raw[:, 11] < cuts["rms"]) & (raw[:, 3] > cuts["depth"]) & (raw[:, 13] < cuts["sigma_eta"])
    igroups = defaultdict(list)
    for row in raw[ok]:
        igroups[(int(row[5]), int(row[6]), int(row[7]), int(row[8]), round(float(row[9]), 2))].append(row)
    intens = []
    for (vu, vl, J, br, T), rows in igroups.items():
        rows = np.array(rows)
        lnA = rows[:, 14] * np.log(10) + rows[:, 12] + np.log(rows[:, 4])
        mean, sig = combine(np.column_stack([np.zeros(len(rows)), lnA, rows[:, 13]]))
        k = key(vu, vl, J, br)
        intens.append(dict(line=k, T=T, A=float(np.exp(mean)), sA=float(np.exp(mean) * sig),
                           N=float(10 ** np.mean(rows[:, 14])), depth=float(rows[:, 3].max()), n=len(rows),
                           blend=int(blend.get(k, True))))

    report(s, lines, intens, ref, ref_name, K, floor, table, scale_curve)
    write(s, lines, intens, K, floor, how, scale_curve, raw)


def report(s, lines, intens, ref, ref_name, K, floor, table, scale_curve):
    nus = [r["nu"] for r in lines.values()]
    print(f"\n{s['id']}: {len(lines)} positions, {min(nus):.0f}-{max(nus):.0f} cm-1, median sigma "
          f"{MHZ_PER_CM * np.median([r['sigma'] for r in lines.values()]):.1f} MHz; {len(intens)} intensities "
          f"(median relative sigma {np.median([i['sA'] / i['A'] for i in intens]):.3f})")
    print(f"  sigma = sqrt(({K:.2f} sigma_fit)^2 + {floor:.1f}^2) MHz, against {ref_name}"
          + (": bins " + ", ".join(f"({a:.1f}, {b:.1f})" for a, b in table) if table else ""))
    if scale_curve is not None:
        p, e = scale_curve
        print("  scale against it (ppb, x = (nu - 17900)/1000): "
              + ", ".join(f"c{len(p) - 1 - j} = {c:+.1f} ± {x:.1f}" for j, (c, x) in enumerate(zip(p, e))))
    # against the published positions of the same spectrum, where there are some
    published = {"salami_ross_2005": "salami_ross_2005", "apo_nist_2009": "apo_nist_2009",
                 "kitt_peak_930317_003": "kitt_peak_1993"}.get(s["id"])
    old = None
    if published:
        from i2spec.observations import read_observations
        old = {str(o.line): (o.value, o.uncertainty) for o in read_observations(ROOT / "data/atlas_lines" / f"{published}.csv")}
    elif s["id"] == "rodriguez_fernandez_2023":
        from i2spec.observations import read_observations
        old = {str(o.line): (o.value, o.uncertainty) for o in read_observations(
            ROOT / "data/external/rodriguez_fernandez_marcassa/rodriguez_fernandez_2023.csv")}
    if old:
        both = [k for k in lines if k in old]
        d = np.array([(lines[k]["nu"] - old[k][0]) * MHZ_PER_CM for k in both])
        z = d / np.array([np.hypot(lines[k]["sigma"], old[k][1]) * MHZ_PER_CM for k in both])
        med, spread = robust(d)
        print(f"  against the positions-only extraction: {len(both)} of its {len(old)} lines here, {len(lines) - len(both)} new; "
              f"new - old median {med:+.2f} MHz, robust spread {spread:.1f} MHz, |z| > 1 for {np.mean(np.abs(z) > 1):.1%}")
    # intensities: their ratio to the model, by band, shows what the template's strengths get wrong
    by_T = defaultdict(list)
    for i in intens:
        by_T[i["T"]].append(i)
    for T, rows in sorted(by_T.items()):
        clean = [r for r in rows if not r["blend"]]
        print(f"  {T:.1f} K: {len(rows)} intensities ({len(clean)} unblended), A {min(r['A'] for r in rows):.2e}-"
              f"{max(r['A'] for r in rows):.2e} cm-1")


def git_state():
    head = subprocess.run(["git", "-C", str(ROOT), "rev-parse", "--short", "HEAD"], capture_output=True, text=True).stdout.strip()
    dirty = subprocess.run(["git", "-C", str(ROOT), "status", "--porcelain", "--untracked-files=no"],
                           capture_output=True, text=True).stdout.strip()
    return head + (" with uncommitted changes" if dirty else "")


def kp_meta(s):
    hdr = s["header"]
    name = s["file"]
    date = f"19{name[:2]}-{name[2:4]}-{name[4:6]}"
    citation = (f"Kitt Peak McMath-Pierce 1-m FTS spectrum {name}, {date}, NSO FTS archive "
                f"(https://nispdata.nso.edu/ftp/FTS_cdrom/). Header ID: {hdr.get('ID', '').strip()}")
    instrument = (f"McMath-Pierce 1-m FTS, RESOLUTN {hdr.get('RESOLUTN', '?')} cm-1, {float(hdr['DELW']):.5f} cm-1 per "
                  f"point, {hdr.get('NSCAN', '?')} co-added scans, source {hdr.get('SOURCE', '?').strip()}, "
                  f"filter {hdr.get('FILTER1', '?').strip()}; vacuum wavenumbers (AIRCORR = {hdr.get('AIRCORR', '?').strip()}), "
                  f"HeNe-referenced (REFWAVNO {hdr.get('REFWAVNO', '?')})")
    T = s["T"] - 273.15
    s["title"] = f"Kitt Peak 1-m FTS, {date}: {s['description']}"
    conditions = (f"{s['description']}. {NOTES.get(name, '')}The template is at {s['T']:.2f} K"
                  + (" (fitted: the lowest mean window rms)" if s.get("T_fitted") else f" ({T:.2f} C, from the header)")
                  + f"; band used {s['band'][0]:.0f}-{s['band'][1]:.0f} cm-1")
    return citation, instrument, conditions


def debus_meta(s):
    d, params = debus_params()
    flat = {k: v for p in params.values() for k, v in p.items()}
    pick = ", ".join(f"{k} {flat[k]}" for k in ("RES", "APF", "ZFF", "LWN", "NSS", "SRC", "BMS", "DTC", "APT", "VEL")
                     if k in flat)
    instrument = (f"Bruker IFS 125HR at IAG Göttingen, evacuated; {int(d['n'])} single-channel spectra of 10 scans each, "
                  f"co-added here; OPUS parameters of the first file: {pick}. HeNe-referenced, not comb-calibrated")
    conditions = (f"Thorlabs GC19100-I cell (19 mm x 100 mm by its part number) in a GCH25-75 heater, halogen lamp; "
                  f"spectra from {d['first']} to {d['last']}. The deposit records no cell temperature or column; the "
                  f"template is at {s['T']:.2f} K, the lowest mean window rms. Band used {s['band'][0]:.0f}-"
                  f"{s['band'][1]:.0f} cm-1 less {HENE_GAP[0]:.0f}-{HENE_GAP[1]:.0f} cm-1 (the 633 nm HeNe artifact)")
    return s["citation"], instrument, conditions


def scale_fields(s, scale_curve):
    """calibration and scale_ppb for a spectrum whose scale was fitted against the model (--calibrate):
    observed = true (1 + scale_ppb 1e-9). The explorer divides it out (webapp._atlas_rows)."""
    if s["kind"] not in ("kitt_peak", "debus", "pyodine"):
        return ""
    ref = (f"; {scale_curve[0][-1]:+.1f} ppb against {s['reference']}" if scale_curve is not None else "")
    what = ("the FTS scale as resampled by the pyodine authors for radial-velocity use, not comb-calibrated"
            if s["kind"] == "pyodine" else "HeNe-referenced FTS scale, not comb-calibrated")
    return (f'calibration = """{what[0].upper() + what[1:]}: {s["prescale_ppb"]:+.1f} ppb against the '
            f'i2spec model (scale_ppb{ref}).\nNot in the per-line uncertainty; the explorer divides it out."""\n'
            f'scale_ppb = {s["prescale_ppb"]:.1f}\n')


def write(s, lines, intens, K, floor, how, scale_curve, raw):
    sid = s["id"]
    DEST = PUBLIC if s.get("public") else globals()["DEST"]
    DEST.mkdir(parents=True, exist_ok=True)
    order = sorted(lines, key=lambda k: (lines[k]["vu"], lines[k]["vl"], k))
    with open(DEST / f"{sid}.csv", "w") as f:
        f.write("line,component,kind,value,uncertainty,ref_line,ref_component,group,note\n")
        for k in order:
            r = lines[k]
            f.write(f"{k},,frequency,{r['nu']:.6f},{r['sigma']:.6f},,,{sid},"
                    f"depth {r['depth']:.2f}; {r['n']} window{'s' if r['n'] > 1 else ''}\n")
    with open(DEST / f"{sid}.intensities.csv", "w") as f:
        f.write("line,T_K,A_cm-1,A_unc_cm-1,N_window_cm-2,depth,windows,blend\n")
        for i in sorted(intens, key=lambda i: (i["T"], i["line"])):
            f.write(f"{i['line']},{i['T']:.2f},{i['A']:.5g},{i['sA']:.3g},{i['N']:.4g},{i['depth']:.3f},{i['n']},{i['blend']}\n")
    if s["kind"] == "kitt_peak":
        citation, instrument, conditions = kp_meta(s)
    elif s["kind"] == "debus":
        citation, instrument, conditions = debus_meta(s)
    elif s["kind"] == "pyodine":
        citation, instrument = s["citation"], s["instrument"]
        conditions = (f"{s['description']}. The template is at {s['T']:.2f} K"
                      + (" (fitted: the lowest mean window rms)" if s.get("T_fitted") else " (the cell temperature)")
                      + f"; band used {s['band'][0]:.0f}-{s['band'][1]:.0f} cm-1")
    else:
        citation, instrument, conditions = s["citation"], s["instrument"], s["conditions"]
    cuts = s["cuts"]
    nus = [r["nu"] for r in lines.values()]
    scale_text = ""
    if scale_curve is not None:
        p, e = scale_curve
        scale_text = ("Against the reference the scale is sum_j c_j x^j ppb, x = (nu - 17900)/1000 cm-1: "
                      + ", ".join(f"c{len(p) - 1 - j} = {c:+.1f} +/- {x:.1f}" for j, (c, x) in enumerate(zip(p, e)))
                      + " (bootstrap over 50 cm-1 blocks). ")
    blend_text = (f"no other model line within +/-{s['blend'][0]:.3f} cm-1 at >= {s['blend'][1]} of its strength"
                  if s["drop_blends"] else "blends kept, as in the published extraction")
    toml = f'''id = "{sid}"
unit = "cm-1"
isotopologue = "127I2"
url = "{s['url']}"
kind = "line positions and intensities condensed from a spectrum"
title = "{s['title']}"
citation = """{citation}"""
instrument = """{instrument}"""
conditions = """{conditions}"""
terms = """{s['terms']}. The spectrum is not in this repository; these are our fits to it."""
extraction = """prototypes/condense_spectra.py (i2spec {git_state()}, {time.strftime('%Y-%m-%d')}): prototypes/atlas_lines.py
extract with amplitudes on, {len(np.unique(raw[:, 10]))} windows. Every model line whose own absorption exceeds
{atlas_lines.MIN_DEPTH} at line centre gets a wavenumber shift (ridge prior {atlas_lines.PRIOR} cm-1) and a strength factor
exp(eta) (ridge prior {atlas_lines.ETA_PRIOR} on eta) on the i2spec template: hyperfine structure, Doppler width at the
template temperature, a Gaussian instrument width, column density, linear baseline and zero offset per window.
The fit axis was divided by 1 + s, s = {s['prescale_ppb']:+.1f} ppb, and shifted by {-s['preshift_mhz']:+.1f} MHz so that the
shift prior does not pull against the calibration; the positions here are back on the spectrum's own axis."""
positions = """{sid}.csv: hyperfine-free line centres (cm-1) on the spectrum's own wavenumber scale, in the
observations format. Repeats combined by inverse variance, sigma inflated by sqrt(chi2/dof). Kept: fit sigma
< {cuts['sigma']:.0f} MHz, window rms < {100 * cuts['rms']:.0f}% of transmission, depth > {cuts['depth']:.2f}; {blend_text}.
The scale carries a calibration offset (and here possibly a curve) that belongs in a fit as a group
parameter and is NOT in the per-line uncertainty. {scale_text}"""
uncertainty = """sqrt(({K:.2f} sigma_fit)^2 + {floor:.1f}^2) MHz, {how}."""
intensities = """{sid}.intensities.csv: per line and template temperature, A = N S (cm-1), the integrated
optical depth of the line (column density times line strength: the area under -ln(transmission) before the
instrument function), and its uncertainty from the window fits (statistical: it does not include what the
template's line shape or a blend's split gets wrong; the blend column flags lines with a model neighbour within
+/-{s['blend'][0]:.3f} cm-1 at >= {s['blend'][1]} of their strength). N_window_cm-2 is the window's fitted column
density, on the scale of the template's strengths, so A / N_window is a strength in cm/molecule on that scale;
with an independent column density (path times number density) A gives absolute strengths. Kept: window rms
< {100 * cuts['rms']:.0f}%, depth > {cuts['depth']:.2f}, sigma(ln A) < {cuts['sigma_eta']}. Repeats combined in ln A.{(" " + s["intensity_note"] + ".") if s.get("intensity_note") else ""}"""
{scale_fields(s, scale_curve)}range_cm-1 = [{min(nus):.1f}, {max(nus):.1f}]
n_lines = {len(lines)}
n_intensities = {len(intens)}
'''
    (DEST / f"{sid}.toml").write_text(toml)
    print(f"  -> {DEST / sid}.{{csv,intensities.csv,toml}}")


def index():
    """lines/README.md: one row per source, from the .toml files."""
    import tomllib
    rows = []
    for p in sorted(DEST.glob("*.toml")):
        m = tomllib.loads(p.read_text())
        lo, hi = m["range_cm-1"]
        rows.append(f"| `{m['id']}` | {m['title']} | {lo:.0f}–{hi:.0f} | "
                    f"{m['n_lines']} | {m['n_intensities']} | {m['terms'].split('.')[0]} |")
    (DEST / "README.md").write_text(
        "# Line positions and intensities condensed from spectra\n\n"
        "Made by `prototypes/condense_spectra.py` in the i2spec repository; each `<id>.toml` says how. Per source:\n\n"
        "- `<id>.csv`: hyperfine-free line positions (cm⁻¹) on the spectrum's own scale, in i2spec's observations\n"
        "  format (`i2spec.observations.read_observations`). The scale's offset is not in the uncertainties.\n"
        "- `<id>.intensities.csv`: `line, T_K, A_cm-1, A_unc_cm-1, N_window_cm-2, depth, windows, blend`, where A = N S\n"
        "  is the line's integrated optical depth (column density times line strength).\n"
        "- `<id>.toml`: the source, conditions, terms, method and uncertainty calibration.\n\n"
        "| id | source | range (cm⁻¹) | positions | intensities | terms |\n|---|---|---|---|---|---|\n"
        + "\n".join(rows) + "\n")
    print(f"{len(rows)} sources -> {DEST / 'README.md'}")


def expand(ids):
    out = []
    for i in ids:
        if i == "all":
            out += list(SOURCES)
        elif i == "kitt_peak":
            out += [k for k in SOURCES if k.startswith("kitt_peak_")]
        elif i in SOURCES:
            out.append(i)
        else:
            raise SystemExit(f"unknown source {i!r}; --list shows them")
    return out


def main(argv):
    opts = {a.split("=")[0]: a.split("=", 1)[1] for a in argv if a.startswith("--") and "=" in a}
    args = [a for a in argv if "=" not in a]
    if not args or args[0] == "--list":
        cal = calibration()
        for sid, s in SOURCES.items():
            done = (OUT / f"{sid}.txt").exists()
            built = (DEST / f"{sid}.csv").exists()
            T = s.get("T") or cal.get(sid, {}).get("T")
            print(f"{sid:28s} {s['kind']:12s} T {T or '?':>7}  fit {'yes' if done else 'no ':3s}  built {'yes' if built else 'no'}")
        return
    if args[0] == "--coadd":
        debus_coadd()
        return
    if args[0] == "--index":
        index()
        return
    ids = expand(args[1:])
    if args[0] == "--calibrate":
        for sid in ids:
            if SOURCES[sid]["kind"] in ("kitt_peak", "debus", "pyodine"):
                calibrate(sid)
    elif args[0] in ("--fit", "--build"):
        for sid in ids:                  # one failure must not stop the others; it is reported at the end
            try:
                fit(sid, int(opts.get("--workers", 6))) if args[0] == "--fit" else build(sid)
            except (Exception, SystemExit) as e:                               # noqa: BLE001
                print(f"{sid}: FAILED {type(e).__name__}: {e}", flush=True)


if __name__ == "__main__":
    main(sys.argv[1:])
