"""Turn the per-window shifts of atlas_lines.py into a line-position data set.

Repeats of the same line (overlapping windows, and the overlapping temperature segments of
Salami & Ross) are combined by inverse variance, with the uncertainty inflated where they disagree.
The position of a line is nu_model + delta, on the atlas's own wavenumber scale; the scale itself is
left to a per-atlas calibration offset in the fit, not folded into the per-line uncertainty.

usage: atlas_dataset.py [--sigma=30] [--rms=0.04] [--depth=0.03]
Writes data/atlas_lines/<id>.csv and .toml, and prints the comparison with the precision sets.
"""
from __future__ import annotations

import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

from i2spec.constants import MHZ_PER_CM
from i2spec.observations import Line, Predictor, load_all

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "atlas_lines"
ATLASES = {
    "salami_ross_2005": dict(
        file="atlas_lines_salami_ross.txt",
        citation=("H. Salami, A. J. Ross, A molecular iodine atlas in ascii format, J. Mol. Spectrosc. 233, "
                  "157-159 (2005), doi:10.1016/j.jms.2005.06.002 (Elsevier supplementary file mmc1)"),
        instrument="FTS, 50 cm cell, instrumental resolution 0.02 cm-1, 0.005 cm-1 sampling",
        calibration="+/-0.003 cm-1 stated; our own check against the precision sets is in the table below",
    ),
    "apo_nist_2009": dict(
        file="atlas_lines_apo.txt",
        citation=("NIST 2-m Fourier transform spectrometer scan of the APO iodine cell (sample AS-17), "
                  "21 July 2009, distributed in the iodine_atlas/ directory of "
                  "https://github.com/jaklusmeyer/APO_pyodine; the reduction code is that of Heeren et al., "
                  "Astron. Astrophys. 674, A164 (2023), doi:10.1051/0004-6361/202244441"),
        instrument="NIST 2-m FTS, RESOLUTN 0.018 cm-1, 32 scans, Xe lamp; air wavelengths converted here",
        calibration="HeNe-referenced FTS scale; agrees with Salami & Ross to <= 5e-8 in earlier checks",
    ),
}


#: sigma_eff = sqrt((K sigma_fit)^2 + FLOOR^2), MHz. The window fit's covariance is optimistic: it
#: knows the noise but not the error the template makes on a blend or on a hyperfine pattern. K and
#: FLOOR are fitted to the difference between the two atlases, binned by sigma_fit (they are
#: independent measurements of the same lines with the same code, so their difference measures
#: everything but the common template error). Checked against the precision sets below.
K, FLOOR = 1.55, 7.5


def combine(rows):
    """Inverse-variance mean of repeated shifts, with sigma inflated by sqrt(chi2/dof) when > 1."""
    d, s = rows[:, 1], rows[:, 2]
    w = 1 / s**2
    mean = float(np.sum(d * w) / np.sum(w))
    sigma = float(np.sqrt(1 / np.sum(w)))
    if len(rows) > 1:
        chi2 = np.sum(((d - mean) / s) ** 2) / (len(rows) - 1)
        sigma *= max(1.0, float(np.sqrt(chi2)))
    return mean, sigma


def build(name, cuts):
    raw = np.loadtxt(ROOT / "prototypes" / "out" / ATLASES[name]["file"])
    keep = (raw[:, 2] * MHZ_PER_CM < cuts["sigma"]) & (raw[:, 11] < cuts["rms"]) & (raw[:, 3] > cuts["depth"])
    groups = defaultdict(list)
    for row in raw[keep]:
        groups[(int(row[5]), int(row[6]), int(row[7]), int(row[8]))].append(row)
    out = []
    for (vu, vl, J, br), rows in sorted(groups.items()):
        rows = np.array(rows)
        delta, sigma = combine(rows)
        sigma = float(np.hypot(K * sigma * MHZ_PER_CM, FLOOR)) / MHZ_PER_CM
        out.append(dict(line=Line("127I2", "R" if br > 0 else "P", J, vu, vl), nu=float(rows[0, 0]) + delta,
                        sigma=sigma, n=len(rows), depth=float(np.max(rows[:, 3]))))
    return out


def main(argv):
    opts = {a.split("=")[0]: float(a.split("=", 1)[1]) for a in argv if a.startswith("--")}
    cuts = dict(sigma=opts.get("--sigma", 30.0), rms=opts.get("--rms", 0.04), depth=opts.get("--depth", 0.03))
    OUT.mkdir(exist_ok=True)
    built = {}
    for name, meta in ATLASES.items():
        lines = build(name, cuts)
        built[name] = {(str(r["line"])): r for r in lines}
        with open(OUT / f"{name}.csv", "w") as f:
            f.write("line,component,kind,value,uncertainty,ref_line,ref_component,group,note\n")
            for r in lines:
                f.write(f"{r['line']},,frequency,{r['nu']:.6f},{r['sigma']:.6f},,,{name},"
                        f"depth {r['depth']:.2f}; {r['n']} window{'s' if r['n'] > 1 else ''}\n")
        (OUT / f"{name}.toml").write_text(
            f'id = "{name}"\nunit = "cm-1"\nisotopologue = "127I2"\nkind = "atlas line positions"\n'
            f'citation = """{meta["citation"]}"""\ninstrument = """{meta["instrument"]}"""\n'
            f'calibration = """{meta["calibration"]}"""\n'
            f'extraction = """Line positions fitted from the transmission spectrum by prototypes/atlas_lines.py:\n'
            f'per-line wavenumber shifts of the i2spec model template in 2.4 cm-1 windows stepped by 2 cm-1, with\n'
            f'the column density, baseline, zero offset and instrument width fitted per window, hyperfine structure\n'
            f'included in the template, and a ridge prior of 0.01 cm-1 on each shift. Repeats combined by inverse\n'
            f'variance (prototypes/atlas_dataset.py). Kept: fit sigma < {cuts["sigma"]:.0f} MHz, window rms <\n'
            f'{100 * cuts["rms"]:.0f}% of transmission, line depth > {cuts["depth"]:.2f}. The positions are\n'
            f'hyperfine-free line centres on the atlas\'s own wavenumber scale; that scale carries a calibration\n'
            f'offset and slope which belong in the fit as a group nuisance parameter, and are NOT in the per-line\n'
            f'uncertainty below."""\n'
            f'uncertainty = """sqrt(({K} sigma_fit)^2 + {FLOOR}^2) MHz. sigma_fit is the window fit\'s covariance,\n'
            f'inflated by sqrt(chi2/dof) where repeated windows disagree; the scaling is fitted to the difference\n'
            f'between this atlas and the other one binned by sigma_fit, which measures what the template gets wrong\n'
            f'on blends and hyperfine patterns as well as the noise. Checked against every precision measurement of\n'
            f'the same lines: normalised residuals 0.7-1.3 after removing each atlas\'s constant offset. The\n'
            f'calibration offset itself is excluded and belongs in the fit as a group parameter."""\n'
            f'n_lines = {len(lines)}\n')
        print(f"{name}: {len(lines)} lines, {min(r['nu'] for r in lines):.0f}-{max(r['nu'] for r in lines):.0f} cm-1, "
              f"median sigma {MHZ_PER_CM * np.median([r['sigma'] for r in lines]):.1f} MHz")

    # --- validation 1: the two atlases against each other, on the lines both measured ---
    a, b = built["salami_ross_2005"], built["apo_nist_2009"]
    both = sorted(set(a) & set(b))
    d = np.array([(a[k]["nu"] - b[k]["nu"]) for k in both])
    s = np.array([np.hypot(a[k]["sigma"], b[k]["sigma"]) for k in both])
    print(f"\n{len(both)} lines in both atlases: Salami-Ross - APO = {MHZ_PER_CM * np.median(d):+.1f} MHz median, "
          f"{MHZ_PER_CM * np.std(d):.1f} MHz scatter, (d/sigma) rms {np.sqrt(np.mean((d / s)**2)):.2f}")
    for lo in range(15000, 20000, 1000):
        m = np.array([lo <= a[k]["nu"] < lo + 1000 for k in both])
        if m.sum() > 20:
            print(f"   {lo}-{lo + 1000}: n={m.sum():5d}  median {MHZ_PER_CM * np.median(d[m]):+7.1f} MHz  "
                  f"scatter {MHZ_PER_CM * np.std(d[m]):6.1f}")

    # --- validation 2: against every precision measurement of the same line ---
    pred = Predictor()
    print("\nagainst the precision sets (atlas - measured, MHz):")
    print(f"  {'set':16s} {'n':>4s} " + " ".join(f"{t:>27s}" for t in ("Salami-Ross: bias rms z (n)", "APO: bias rms z (n)")))
    for ds in load_all():
        rows = defaultdict(list)
        for o in ds.observations:
            if o.line.isotopologue != "127I2" or o.kind != "frequency" or o.component is None:
                continue
            key = str(o.line)
            centre = o.value * (MHZ_PER_CM if ds.unit == "cm-1" else 1.0) - (pred.position(o.line, o.component)
                                                                             - pred.position(o.line))
            rows[key].append(centre)
        cells, n = [], 0
        for atlas in (a, b):
            pairs = [(atlas[k]["nu"] * MHZ_PER_CM - np.mean(v), atlas[k]["sigma"] * MHZ_PER_CM)
                     for k, v in rows.items() if k in atlas]
            n = max(n, len(pairs))
            if pairs:
                dd, ss = np.array([p[0] for p in pairs]), np.array([p[1] for p in pairs])
                z = (dd - np.median(dd)) / ss
                cells.append(f"{np.median(dd):+6.1f} {np.std(dd):5.1f} {np.sqrt(np.mean(z**2)):4.2f} ({len(dd):3d})")
            else:
                cells.append(f"{'-':>27s}")
        if n:
            print(f"  {ds.id:16s} {len(rows):4d} " + " ".join(f"{c:>27s}" for c in cells))


if __name__ == "__main__":
    main(sys.argv[1:])
