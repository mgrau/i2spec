"""Line centres from the 1993 Kitt Peak McMath-Pierce 1-m FTS scans of RV iodine cells, and their scale.

The scans (data/external/kitt_peak_fts/, NSO archive, volume FTS33) are raw FITS spectra: vacuum
wavenumbers WSTART + i DELW (AIRCORR = No, tank evacuated), WAVCORR = 0 (no scale correction applied),
0.037 cm-1 resolution sampled every 0.0141 cm-1, DZE lamp continuum not divided out. Useful band
15 875-20 000 cm-1 (BANDLO/BANDHI).

Each scan is normalised by an upper envelope (the 99.5th percentile in 4 cm-1 blocks, interpolated),
then fitted as the NIST/APO scan was (prototypes/atlas_lines.py extract: per-line shifts of the i2spec
template), in 3.2 cm-1 windows stepped by 2.4 cm-1 (the extractor needs 200 points per window). The
dataset step combines repeats, leaves out blends, fits the scan's scale against the model (a quadratic
in wavenumber: a single scale does not fit), and writes data/atlas_lines/kitt_peak_1993.{csv,toml}.
Notes: docs/research/kitt-peak-and-rodriguez-spectra.md.

reproduce: kitt_peak_lines.py --file=930317R0.003 --workers=3
           kitt_peak_lines.py --file=930317R0.005 --T=343.15 --step=6    (and .008, .001 at 323.15 K)
           kitt_peak_lines.py --dataset

usage: kitt_peak_lines.py --one CENTRE [--file=930317R0.003] [--T=323.15]
       kitt_peak_lines.py [--file=930317R0.003] [--T=323.15] [--workers=2] [--step=2.4] window fits
       kitt_peak_lines.py --dataset [--files=930317R0.003,...]                            csv/toml + report
"""
from __future__ import annotations

import os
import sys
import time
from collections import defaultdict
from multiprocessing import Pool
from pathlib import Path

import numpy as np

os.environ.setdefault("I2SPEC_WORKERS", "2")          # keep any line-list rebuild modest
sys.path.insert(0, str(Path(__file__).resolve().parent))
import atlas_lines                                                           # noqa: E402
from atlas_dataset import combine                                            # noqa: E402
from i2spec.constants import MHZ_PER_CM                                     # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "external" / "kitt_peak_fts"
OUT = ROOT / "prototypes" / "out"
ATLAS = ROOT / "data" / "atlas_lines"
#: A first estimate of the scan's scale (obs = true (1 + s)), taken out of the axis before the window
#: fits so that the 0.01 cm-1 ridge prior on each shift does not pull against a -700 MHz offset; the
#: dataset step puts it back, so the published positions are on the scan's own scale.
PRESCALE_PPB = -1330.0
ROW_FMT = "%.6f %.7f %.7f %.4f %.4g %d %d %d %d %.2f %.3f %.5f\n"


# --- reading ---------------------------------------------------------------------------------------

def read_fits(path):
    """(header dict, wavenumber, signal) of an NSO FTS FITS file: 2880-byte header blocks, BITPIX -32,
    byte order from BOCODE (1 = big-endian), x axis WSTART + i DELW."""
    raw = Path(path).read_bytes()
    hdr, i = {}, 0
    while True:
        card = raw[i:i + 80].decode("ascii", "replace"); i += 80
        key = card[:8].strip()
        if card[8:10] == "= ":
            hdr[key] = card[10:].split("/")[0].strip().strip("'").strip()
        if key == "END":
            break
    assert hdr["BITPIX"] == "-32" and hdr["XAXIS_IS"].startswith("Index") and hdr["DATA_IS"].startswith("Real")
    off, n = -(-i // 2880) * 2880, int(hdr["NPO"])
    y = np.frombuffer(raw[off:off + 4 * n], dtype=">f4" if hdr["BOCODE"] == "1" else "<f4").astype(float)
    nu = float(hdr["WSTART"]) + np.arange(n) * float(hdr["DELW"])
    assert abs(nu[-1] - float(hdr["WSTOP"])) < 1e-4
    return hdr, nu, y


def envelope(nu, y, block=4.0, q=99.5):
    """Upper envelope of a lamp spectrum: a high percentile per block, interpolated linearly."""
    edges = np.arange(nu[0], nu[-1] + block, block)
    c, v = [], []
    for a, b in zip(edges[:-1], edges[1:]):
        m = (nu >= a) & (nu < b)
        if m.sum() > 50:
            c.append(nu[m].mean()); v.append(np.percentile(y[m], q))
    return np.interp(nu, c, v)


def load(name, prescale_ppb=PRESCALE_PPB):
    """(nu, per cent transmission) of one scan over its useful band, as atlas_lines expects, with the
    axis divided by (1 + prescale)."""
    hdr, nu, y = read_fits(RAW / f"{name}.fits")
    nu = nu / (1 + prescale_ppb * 1e-9)
    if "Iodine" not in hdr["ID"]:
        raise ValueError(f"{name} is not an iodine cell: {hdr['ID']}")
    m = (nu >= float(hdr.get("BANDLO", 15875))) & (nu <= float(hdr.get("BANDHI", 20000)))
    nu, y = nu[m], y[m]
    return hdr, np.column_stack([nu, 100.0 * y / envelope(nu, y)])


# --- window fits -----------------------------------------------------------------------------------

def init(name, T):
    atlas_lines.init("salami_ross", T)           # model, master list, continuum, hyperfine model
    atlas_lines._state["data"] = load(name)[1]


def windows(data, step, width, T):
    lo, hi = data[0, 0], data[-1, 0]
    return [(float(c), width, T) for c in np.arange(np.ceil(lo) + width / 2, hi - width / 2, step)]


def run(name, T, workers, step, width=3.2):
    hdr, data = load(name)
    todo = windows(data, step, width, T)
    out = OUT / f"atlas_lines_kitt_peak_{name}.txt"
    t = time.time()
    print(f"{name}: {hdr['ID']}; {len(todo)} windows at {T} K, {workers} workers", flush=True)
    with Pool(workers, initializer=init, initargs=(name, T)) as pool, open(out, "w") as f:
        f.write("# " + atlas_lines.HEADER + "\n")
        for n, rows in enumerate(pool.imap_unordered(atlas_lines.safe_extract, todo, chunksize=2), 1):
            for row in rows:
                f.write(ROW_FMT % tuple(row))
            f.flush()
            if n % 50 == 0:
                print(f"  {n}/{len(todo)} windows, {time.time() - t:.0f} s", flush=True)
    print(f"done in {time.time() - t:.0f} s -> {out}")


# --- the data set ----------------------------------------------------------------------------------

def build(raw, cuts, K, floor, prescale_ppb=PRESCALE_PPB):
    """Repeats of one line (overlapping windows, several scans) combined by inverse variance; sigma
    = sqrt((K sigma_fit)^2 + floor^2), MHz. The positions are put back on the scan's own scale
    (the window fits ran on an axis divided by 1 + prescale). Returns {line key: dict}."""
    keep = (raw[:, 2] * MHZ_PER_CM < cuts["sigma"]) & (raw[:, 11] < cuts["rms"]) & (raw[:, 3] > cuts["depth"])
    groups = defaultdict(list)
    for row in raw[keep]:
        groups[(int(row[5]), int(row[6]), int(row[7]), int(row[8]))].append(row)
    out = {}
    for (vu, vl, J, br), rows in sorted(groups.items()):
        rows = np.array(rows)
        delta, sigma = combine(rows)
        key = f"127I2 {'R' if br > 0 else 'P'}({J}) {vu}-{vl}"
        out[key] = dict(nu=(float(rows[0, 0]) + delta) * (1 + prescale_ppb * 1e-9), delta=delta, sigma_fit=sigma,
                        sigma=float(np.hypot(K * sigma * MHZ_PER_CM, floor)) / MHZ_PER_CM,
                        n=len(rows), depth=float(np.max(rows[:, 3])), vu=vu, vl=vl)
    return out


def drop_blends(lines, T, half, ratio):
    """Leave out lines with another model line within +/- half cm-1 at least ratio times as strong (master
    line list at T): the window fit gives only the strong lines a parameter, so a neighbour below its
    threshold is held at its model position and a blend can still come out with a small sigma."""
    from scan_salami_ross import _init, _state
    if "master" not in _state:
        _init()
    nus = [r["nu"] for r in lines.values()]
    ll = _state["master"].at(T, min(nus) - 1, max(nus) + 1, S_min=1e-27)
    index = {(int(ll.v_upper[i]), int(ll.v_lower[i]), int(ll.J_lower[i]), 1 if ll.branch[i] > 0 else -1): i
             for i in range(len(ll.nu))}
    keep = {}
    for k, r in lines.items():
        br, J = (1 if k.split()[1][0] == "R" else -1), int(k.split()[1][2:-1])
        i = index.get((r["vu"], r["vl"], J, br))
        if i is None:
            continue
        m = np.abs(ll.nu - ll.nu[i]) < half
        m[i] = False
        if not m.any() or ll.S[m].max() < ratio * ll.S[i]:
            keep[k] = r
    return keep


def model_centres(keys, parameters=None):
    """Hyperfine-free model line centres (MHz), as paper/figures/rows.py atlas_rows computes them."""
    from i2spec.constants import DEFAULT_PARAMETERS
    from i2spec.observations import Predictor
    model = Predictor(parameters or DEFAULT_PARAMETERS).model("127I2")
    out = {}
    for k in keys:
        br, rest = k.split()[1][0], k.split()[1]
        J = int(rest[2:-1]); vu, vl = map(int, k.split()[2].split("-"))
        try:
            out[k] = model.transition(vu, vl, J, br) * MHZ_PER_CM
        except Exception:                                            # noqa: BLE001
            pass
    return out


def robust(x):
    """median and 1.4826 * MAD."""
    x = np.asarray(x)
    m = np.median(x)
    return float(m), float(1.4826 * np.median(np.abs(x - m)))


def clipped_fit(X, y, s, nsig=4.0, iters=6):
    """Weighted least squares with iterative clipping at nsig robust sigma; returns (p, cov, mask, scale)."""
    m = np.ones(len(y), bool)
    for _ in range(iters):
        w = 1 / s[m]
        p, *_ = np.linalg.lstsq(X[m] * w[:, None], y[m] * w, rcond=None)
        r = (y - X @ p) / s
        _, spread = robust(r[m])
        m_new = np.abs(r - np.median(r[m])) < nsig * max(spread, 1e-9)
        if (m_new == m).all():
            break
        m = m_new
    w = 1 / s[m]
    A = X[m] * w[:, None]
    chi2 = np.sum(((y[m] - X[m] @ p) * w) ** 2) / max(m.sum() - X.shape[1], 1)
    cov = np.linalg.inv(A.T @ A) * max(chi2, 1.0)
    return p, cov, m, float(np.sqrt(chi2))


def reference_atlases():
    from i2spec.observations import read_observations
    return {n: {str(o.line): (o.value, o.uncertainty) for o in read_observations(ATLAS / f"{n}.csv")}
            for n in ("salami_ross_2005", "apo_nist_2009")}


def compare(lines, refs, correction=0.0, label=""):
    """This set minus each reference atlas, MHz, on the lines both hold, after subtracting a calibration:
    a constant scale in ppb, or a function of wavenumber returning MHz."""
    corr = correction if callable(correction) else (lambda nu: correction * 1e-9 * nu * MHZ_PER_CM)
    for name, ref in refs.items():
        both = [k for k in lines if k in ref]
        if not both:
            print(f"  {label} vs {name}: no common lines"); continue
        d = np.array([(lines[k]["nu"] - ref[k][0]) * MHZ_PER_CM - corr(lines[k]["nu"]) for k in both])
        s = np.array([np.hypot(lines[k]["sigma"], ref[k][1]) for k in both]) * MHZ_PER_CM
        med, spread = robust(d)
        z = (d - med) / s
        print(f"  {label} - {name}: n={len(both)}  median {med:+.1f} MHz  robust spread {spread:.1f} MHz  "
              f"z rms (4-sigma clipped) {np.sqrt(np.mean(z[np.abs(z) < 4]**2)):.2f}")
        for lo in range(14000, 20000, 1000):
            m = np.array([lo <= lines[k]["nu"] < lo + 1000 for k in both])
            if m.sum() >= 20:
                a, b = robust(d[m])
                print(f"      {lo}-{lo + 1000}: n={m.sum():5d}  median {a:+7.1f}  spread {b:5.1f}")


def calibrate_sigma(lines, ref, trend=None):
    """K, floor such that d^2 ~ sigma_ref^2 + (K sigma_fit)^2 + floor^2, fitted to the scatter of
    (this - reference) in bins of sigma_fit after removing the median (and trend(nu), MHz, if given)."""
    both = [k for k in lines if k in ref]
    d = np.array([lines[k]["nu"] - ref[k][0] for k in both]) * MHZ_PER_CM
    if trend is not None:
        d = d - np.array([trend(lines[k]["nu"]) for k in both])
    sf = np.array([lines[k]["sigma_fit"] for k in both]) * MHZ_PER_CM
    sr = np.array([ref[k][1] for k in both]) * MHZ_PER_CM
    d = d - np.median(d)
    edges = np.percentile(sf, np.linspace(0, 100, 9))
    x, y = [], []
    for a, b in zip(edges[:-1], edges[1:]):
        m = (sf >= a) & (sf <= b)
        if m.sum() > 20:
            spread = robust(d[m])[1]
            x.append(np.median(sf[m]) ** 2); y.append(max(spread**2 - np.median(sr[m]) ** 2, 0.0))
    A = np.column_stack([x, np.ones(len(x))])
    (k2, f2), *_ = np.linalg.lstsq(A, np.array(y), rcond=None)
    return float(np.sqrt(max(k2, 1.0))), float(np.sqrt(max(f2, 1.0))), list(zip(np.sqrt(x), np.sqrt(y)))


CUTS = dict(sigma=30.0, rms=0.04, depth=0.05)
KITT_FILES = ("930317R0.003", "930317R0.005", "930317R0.008", "930317R0.001")   # the first is published; the rest (6 cm-1 steps) test the curve


#: blends (drop_blends): half the resolution, and neighbours down to 0.1 of the line's strength. Kept lines
#: with such a neighbour scattered 28-34 MHz about the scale curve against 21 MHz for the others.
BLEND = (0.02, 0.1)
NU0 = 17900.0     # cm-1, centre of the calibration polynomial


def curve(nu, r, s, deg=2, block=50.0, nboot=200):
    """The scan's scale as a polynomial in x = (nu - NU0)/1000 cm-1, in ppb: r(MHz) = nu f(nu) 1e-9 sum c_j x^j.
    Weighted, 4-sigma clipped; uncertainties from a bootstrap over 50 cm-1 blocks (lines in one block share
    template and model errors). Returns coefficients (highest power first), their errors, mask, fn(nu)->MHz."""
    f = nu * MHZ_PER_CM
    X = np.vander((nu - NU0) / 1000, deg + 1) * f[:, None] * 1e-9
    p, _, m, chi = clipped_fit(X, r, s)
    rng = np.random.default_rng(1)
    blocks = np.floor(nu / block).astype(int)
    ub = np.unique(blocks)
    boot = []
    for _ in range(nboot):
        idx = np.concatenate([np.nonzero((blocks == b) & m)[0] for b in rng.choice(ub, len(ub))])
        w = 1 / s[idx]
        boot.append(np.linalg.lstsq(X[idx] * w[:, None], r[idx] * w, rcond=None)[0])
    err = np.std(boot, axis=0)

    def fn(x, p=p):
        return float(np.polyval(p, (x - NU0) / 1000) * x * MHZ_PER_CM * 1e-9)
    return p, err, m, chi, fn


def dataset(files):
    refs = reference_atlases()
    apo = refs["apo_nist_2009"]
    per_file = {}
    for name in files:
        per_file[name] = np.loadtxt(OUT / f"atlas_lines_kitt_peak_{name}.txt")
    main_file = files[0]
    raw = per_file[main_file]
    lines0 = build(raw, CUTS, 1.0, 0.0)
    n0 = len(lines0)
    lines0 = drop_blends(lines0, 323.15, *BLEND)
    print(f"{main_file}: {len(raw)} window-lines, {n0} after cuts, {len(lines0)} without blends "
          f"(a neighbour >= {BLEND[1]} x as strong within {BLEND[0]} cm-1)")

    # the scale against the APO scan first (model independent), to detrend the sigma calibration
    ka = [k for k in lines0 if k in apo]
    nua = np.array([lines0[k]["nu"] for k in ka])
    da = np.array([(lines0[k]["nu"] - apo[k][0]) * MHZ_PER_CM for k in ka])
    pa, ea, *_, fa = curve(nua, da, np.ones_like(da))
    print("scale vs apo_nist_2009 (ppb, x = (nu - 17900)/1000): "
          + ", ".join(f"c{len(pa) - 1 - j} = {c:+.1f} ± {e:.1f}" for j, (c, e) in enumerate(zip(pa, ea))))
    K, floor, table = calibrate_sigma(lines0, apo, trend=fa)
    print(f"sigma model vs APO (detrended): K = {K:.2f}, floor = {floor:.1f} MHz; bins (sigma_fit, excess spread): "
          + ", ".join(f"({a:.1f}, {b:.1f})" for a, b in table))
    lines = {k: v for k, v in build(raw, CUTS, K, floor).items() if k in lines0}

    # the scale against the model
    mc = model_centres(lines)
    keys = [k for k in lines if k in mc]
    nu = np.array([lines[k]["nu"] for k in keys])
    r = np.array([lines[k]["nu"] * MHZ_PER_CM - mc[k] for k in keys])
    s = np.array([lines[k]["sigma"] * MHZ_PER_CM for k in keys])
    med, spread = robust(r)
    print(f"\n{len(keys)} lines vs model: median {med:+.1f} MHz, robust spread {spread:.1f} MHz")
    for lo in range(15500, 20000, 500):
        m = (nu >= lo) & (nu < lo + 500)
        if m.sum() > 20:
            a, b = robust(r[m]); rr = robust(r[m] / (nu[m] * MHZ_PER_CM) * 1e9)[0]
            print(f"   {lo}-{lo + 500}: n={m.sum():5d}  median {a:+7.1f} MHz ({rr:+6.1f} ppb)  spread {b:5.1f}")
    fits = {}
    for deg in (0, 1, 2, 3):
        p, e, m, chi, fn = curve(nu, r, s, deg)
        res = r - np.array([fn(x) for x in nu])
        bins = [np.median(res[(nu >= a) & (nu < a + 250) & m]) for a in range(15750, 20000, 250)
                if ((nu >= a) & (nu < a + 250) & m).sum() > 20]
        fits[deg] = (p, e, m, chi, fn, res)
        print(f"deg {deg}: " + ", ".join(f"c{len(p) - 1 - j} = {c:+.1f} ± {ee:.1f}" for j, (c, ee) in enumerate(zip(p, e)))
              + f" ppb; residual spread {robust(res[m])[1]:.1f} MHz, 250 cm-1 bin medians rms "
              f"{np.sqrt(np.mean(np.square(bins))):.1f} MHz (max {np.max(np.abs(bins)):.1f}), z rms {np.sqrt(np.mean((res[m] / s[m])**2)):.2f}")
    p, e, m, chi, fn, res = fits[2]
    print("after the quadratic scale, by region (MHz):")
    for lo in range(15500, 20000, 500):
        mm = (nu >= lo) & (nu < lo + 500)
        if mm.sum() > 20:
            a, b = robust(res[mm])
            print(f"   {lo}-{lo + 500}: n={mm.sum():5d}  median {a:+6.1f}  spread {b:5.1f}")

    # the other cells (coarser window step): is the curve common to the session?
    for name in files[1:]:
        ln = build(per_file[name], CUTS, K, floor)
        mo = model_centres(ln)
        kk = [k for k in ln if k in mo]
        nn = np.array([ln[k]["nu"] for k in kk])
        rr = np.array([ln[k]["nu"] * MHZ_PER_CM - mo[k] for k in kk])
        ss = np.array([ln[k]["sigma"] * MHZ_PER_CM for k in kk])
        pp, ee, mm, _, ff = curve(nn, rr, ss, 2)
        diff = rr - np.array([fn(x) for x in nn])
        print(f"   {name} ({read_fits(RAW / f'{name}.fits')[0]['ID'].split(',  ')[0].strip()}): {len(kk)} lines, "
              + ", ".join(f"c{2 - j} = {c:+.1f} ± {x:.1f}" for j, (c, x) in enumerate(zip(pp, ee)))
              + f" ppb; minus the {main_file} curve: median {np.median(diff):+.1f} MHz, spread {robust(diff)[1]:.1f}")

    print("\nagainst the other atlases, before the scale correction:")
    compare(lines, refs, 0.0, "Kitt Peak")
    print("after removing the quadratic scale fitted against the model:")
    compare(lines, refs, fn, "Kitt Peak")
    for name in ("salami_ross_2005", "apo_nist_2009"):
        mref = model_centres([k for k in refs[name] if k in lines])
        rr = [refs[name][k][0] * MHZ_PER_CM - mref[k] for k in mref]
        print(f"  {name} - model on the same lines: n={len(rr)} median {np.median(rr):+.1f} MHz, "
              f"spread {robust(rr)[1]:.1f}")
    print(f"  lines in neither atlas: {sum(1 for k in lines if k not in refs['apo_nist_2009'] and k not in refs['salami_ross_2005'])}"
          f"; not in APO: {sum(1 for k in lines if k not in refs['apo_nist_2009'])}")
    write(lines, [main_file], K, floor, fits, med, spread, pa, ea)


def write(lines, files, K, floor, fits, med, spread, pa, ea):
    hdrs = {n: read_fits(RAW / f"{n}.fits")[0] for n in files}
    name = "kitt_peak_1993"
    p, e, m, chi, fn, res = fits[2]
    p0 = fits[0][0][0]; e0 = fits[0][1][0]
    with open(ATLAS / f"{name}.csv", "w") as f:
        f.write("line,component,kind,value,uncertainty,ref_line,ref_component,group,note\n")
        for k, r in sorted(lines.items(), key=lambda kv: (kv[1]["vu"], kv[1]["vl"], kv[0])):
            f.write(f"{k},,frequency,{r['nu']:.6f},{r['sigma']:.6f},,,{name},"
                    f"depth {r['depth']:.2f}; {r['n']} window{'s' if r['n'] > 1 else ''}\n")
    ids = "; ".join(f"{n}: {hdrs[n]['ID'].split(',  ')[0].strip()}" for n in files)
    nus = [r["nu"] for r in lines.values()]
    c2, c1, c0 = p; e2, e1, e0q = e
    (ATLAS / f"{name}.toml").write_text(
        f'id = "{name}"\nunit = "cm-1"\nisotopologue = "127I2"\n'
        f'url = "https://nispdata.nso.edu/ftp/FTS_cdrom/"\nkind = "atlas line positions"\n'
        f'citation = """Kitt Peak McMath-Pierce 1-m FTS, 17 March 1993 (observer Marcy), NSO Digital Library\n'
        f'volume FTS33, file {ids}. The RV cell templates of Butler et al., PASP 108, 500\n'
        f'(1996), doi:10.1086/133755."""\n'
        f'instrument = """McMath-Pierce 1-m FTS, RESOLUTN 0.037 cm-1, 0.01409 cm-1 per point, 10 co-added scans,\n'
        f'DZE lamp, 8 mm aperture, evacuated tank; vacuum wavenumbers (AIRCORR = No), WAVCORR = 0 (no scale\n'
        f'correction applied in the archive). Keck cell, 100 mm, 50.0 C; useful band 15 875-20 000 cm-1."""\n'
        f'calibration = """HeNe-referenced (REFWAVNO 15798.0025 cm-1), uncorrected, and NOT a single scale.\n'
        f'Against i2spec2026m the lines read {med:+.0f} MHz (median; robust spread {spread:.0f} MHz), from\n'
        f'-1210 ppb at 15 750 cm-1 to -1700 ppb at 19 750 cm-1. A constant scale fits {p0:+.0f} +/- {e0:.0f} ppb but\n'
        f'leaves +/-100 MHz of structure. The scan\'s scale, obs = model (1 + s(nu)), is\n'
        f'  s(nu) = c0 + c1 x + c2 x^2 ppb, x = (nu - {NU0:.0f} cm-1)/1000 cm-1,\n'
        f'  c0 = {c0:+.1f} +/- {e0q:.1f}, c1 = {c1:+.1f} +/- {e1:.1f}, c2 = {c2:+.1f} +/- {e2:.1f}\n'
        f'(block bootstrap over 50 cm-1), which leaves 250 cm-1 bin medians within a few MHz. The same curve\n'
        f'against apo_nist_2009 is c0 = {pa[2]:+.1f}, c1 = {pa[1]:+.1f}, c2 = {pa[0]:+.1f} ppb (APO itself reads -18 MHz\n'
        f'from the model). Divide the values by (1 + s(nu)) to put them on the model\'s scale. The curve is NOT\n'
        f'in the per-line uncertainty; a fit using this set needs these three group parameters, not one."""\n'
        f'extraction = """prototypes/kitt_peak_lines.py: the scan divided by its upper envelope (99.5th percentile\n'
        f'in 4 cm-1 blocks) and its axis by (1 - 1.33e-6) (restored afterwards, so the 0.01 cm-1 ridge prior\n'
        f'does not pull against a -700 MHz offset), then fitted by prototypes/atlas_lines.py extract at 323.15 K:\n'
        f'per-line wavenumber shifts of the i2spec template in 3.2 cm-1 windows (the extractor needs 200 points)\n'
        f'stepped by 2.4 cm-1; column density, linear baseline, zero offset and Gaussian instrument width per\n'
        f'window; hyperfine structure in the template; ridge prior 0.01 cm-1. Repeats combined by inverse\n'
        f'variance. Kept: fit sigma < {CUTS["sigma"]:.0f} MHz, window rms < {100 * CUTS["rms"]:.0f}% of transmission, line depth >\n'
        f'{CUTS["depth"]:.2f}, and no other model line within +/-{BLEND[0]} cm-1 at >= {BLEND[1]} of its strength. Assignments are the template\'s, one line per parameter; blends come out with large\n'
        f'sigma and fail the cut. Positions are hyperfine-free line centres on the scan\'s own axis."""\n'
        f'uncertainty = """sqrt(({K:.2f} sigma_fit)^2 + {floor:.1f}^2) MHz, the scaling fitted to the scatter of\n'
        f'(Kitt Peak - apo_nist_2009), after removing the scale curve, in bins of sigma_fit and after removing\n'
        f'the APO uncertainty. Excludes the scale curve."""\n'
        f'scale_ppb_x0_cm-1 = {NU0:.1f}\nscale_ppb_coefficients = [{c0:.2f}, {c1:.2f}, {c2:.2f}]\n'
        f'scale_ppb_coefficient_uncertainties = [{e0q:.2f}, {e1:.2f}, {e2:.2f}]\n'
        f'range_cm-1 = [{min(nus):.1f}, {max(nus):.1f}]\nn_lines = {len(lines)}\n')
    print(f"wrote {name}: {len(lines)} lines")


def main(argv):
    opts = {a.split("=")[0]: a.split("=", 1)[1] for a in argv if a.startswith("--") and "=" in a}
    flags = [a for a in argv if not ("=" in a)]
    name, T = opts.get("--file", "930317R0.003"), float(opts.get("--T", 323.15))
    if "--dataset" in flags:
        return dataset(opts.get("--files", ",".join(KITT_FILES)).split(","))
    if "--one" in flags:
        centre = float(flags[flags.index("--one") + 1])
        init(name, T)
        for c in (centre, centre + 2.0):
            t = time.time(); rows = atlas_lines.extract((c, 3.2, T))
            print(f"window {c}: {len(rows)} lines in {time.time() - t:.1f} s; rms {rows[0][-1] if rows else 0:.4f}")
        for nu, d, sd, depth, S, vu, vl, J, br, *_ in sorted(rows):
            print(f"  {'R' if br > 0 else 'P'}({J:3d}) {vu:2d}-{vl:<2d} {nu:11.4f}  δ {d * MHZ_PER_CM:+8.1f} ± "
                  f"{sd * MHZ_PER_CM:6.1f} MHz  depth {depth:4.2f}")
        return
    run(name, T, int(opts.get("--workers", 2)), float(opts.get("--step", 2.4)))


if __name__ == "__main__":
    main(sys.argv[1:])
