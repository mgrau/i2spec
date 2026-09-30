"""The model's A-X, B-X and C-X absorption against Tellinghuisen's fitted component bands.

Reference: J. Tellinghuisen, J. Chem. Phys. 135, 054301 (2011) [T11C], supplement Table IIS (molar
absorptivity of A-X, B-X and C-X and their total at 0 and 35 °C, every 5 nm from 400 nm, with LS standard
errors at 35 °C) and Table IS (measured ε at 35.4 and 64.0 °C). Read from
data/external/tellinghuisen_2011_supplement/JT-Supplement-JCP-I2.txt (gitignored, third-party).

Conventions
* ε is decadic molar absorptivity (A = ε l C, [T11B] eq. 2), so σ = ε · 1000 ln10 / N_A = 3.8235e-21 ε cm².
* Table wavelengths are read as air and converted to vacuum (continuum-model.md §10: every check of the
  model against [T11C] Table II favours air). The vacuum reading is reported as a sensitivity.
* Model (default parameter set): A and C are continuum.py's A<-X and C<-X (Table I potentials and moments).
  B is what the model has for B<-X at that wavelength: the bound-free continuum (Hannover B with the
  [T11C] inner wall, μ_B of [T11B] eq. 10) and, at λ >= 500 nm, also the discrete line list; both are
  averaged over a Gaussian of SMOOTH_NM FWHM there so they can be set against Tellinghuisen's smooth
  pseudocontinuum B. At λ <= 495 nm (above the B dissociation limit for every X level) B is purely
  bound-free and point values are used.

    uv run --group research python prototypes/tellinghuisen_components.py
"""
import json
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from i2spec.continuum import continuum_model
from i2spec.intensity import intensity_model, master_line_list
from i2spec.spectrum import air_to_vacuum

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "data/external/tellinghuisen_2011_supplement/JT-Supplement-JCP-I2.txt"
OUT_JSON = ROOT / "prototypes/out/tellinghuisen_components.json"
OUT_PNG = ROOT / "docs/figures/tellinghuisen_components.png"
N_A = 6.02214076e23
EPS = 1000 * np.log(10) / N_A          # cm² per (l mol⁻¹ cm⁻¹) = 3.8235e-21
TEMPS = {"0": 273.15, "35": 308.15}
SMOOTH_NM = 4.0                        # Gaussian FWHM for the banded B-X region (λ >= 500 nm)
SMOOTH_ALT = (2.0, 8.0)                # sensitivity to that choice
COMPS = ("A", "B", "C", "total")


def parse():
    text = SRC.read_text(encoding="utf-8")
    t1, t2 = text.split("Table IIS.")
    t2 = t2.split("Least-squares normalized residuals")[0]
    rows = {}
    for line in t2.splitlines():
        v = line.split()
        if not v or not re.fullmatch(r"\d{3}", v[0]):
            continue
        lam, x = int(v[0]), [float(s) for s in v[1:]]
        r = {c: (np.nan, np.nan, np.nan) for c in COMPS}
        if len(x) == 12:
            r.update(A=tuple(x[0:3]), B=tuple(x[3:6]), C=tuple(x[6:9]), total=tuple(x[9:12]))
        elif len(x) == 9:          # 400-480 nm: no A-X column
            r.update(A=(0.0, 0.0, 0.0), B=tuple(x[0:3]), C=tuple(x[3:6]), total=tuple(x[6:9]))
        elif len(x) == 6:          # 730-850 nm: A-X and total only
            r.update(A=tuple(x[0:3]), B=(0.0, 0.0, 0.0), C=(0.0, 0.0, 0.0), total=tuple(x[3:6]))
        else:
            raise ValueError(line)
        rows[lam] = r
        if lam == 850:             # last row; the residual plot's axis labels follow
            break
    meas = []
    for line in t1.splitlines():
        v = line.split()
        if len(v) == 5 and re.fullmatch(r"\d{3}", v[0]):
            meas.append([float(s) for s in v])
    return rows, np.array(meas)


def gauss_avg(nu_c, sig, lam_vac, fwhm_nm):
    out = []
    for lv in np.atleast_1d(lam_vac):
        n0, fw = 1e7 / lv, 1e7 * fwhm_nm / lv**2
        m = np.abs(nu_c - n0) < 3 * fw
        g = np.exp(-4 * np.log(2) * ((nu_c[m] - n0) / fw) ** 2)
        out.append((g * sig[m]).sum() / g.sum())
    return np.array(out)


def model_components(master, cont, lam_vac, T, fwhm=SMOOTH_NM):
    """ε (l mol⁻¹ cm⁻¹) of A, B, C, total at vacuum wavelengths lam_vac and temperature T."""
    nu = 1e7 / lam_vac
    A = cont.cross_section(nu, T, states=("A",)) / EPS
    C = cont.cross_section(nu, T, states=("C",)) / EPS
    Bpt = cont.cross_section(nu, T, states=("B",)) / EPS
    lines = master.at(T, S_min=0.0)
    edges = np.arange(1e7 / 780.0, 1e7 / 480.0, 0.05)
    cen = (edges[1:] + edges[:-1]) / 2
    sig_b = np.histogram(lines.nu, edges, weights=lines.S)[0] / 0.05 + cont.cross_section(cen, T, states=("B",))
    banded = lam_vac >= 499.5
    B = Bpt.copy()
    B_lines_only = np.zeros_like(Bpt)
    sel = banded & (lam_vac < 760)
    B[sel] = gauss_avg(cen, sig_b, lam_vac[sel], fwhm) / EPS
    sig_l = np.histogram(lines.nu, edges, weights=lines.S)[0] / 0.05
    B_lines_only[sel] = gauss_avg(cen, sig_l, lam_vac[sel], fwhm) / EPS
    return dict(A=A, B=B, C=C, total=A + B + C, B_boundfree=Bpt, B_lines=B_lines_only)


def main():
    rows, meas = parse()
    lam = np.array(sorted(rows), dtype=float)
    ref = {c: {t: np.array([rows[int(l)][c][i] for l in lam]) for i, t in enumerate(("0", "35"))} for c in COMPS}
    err = {c: np.array([rows[int(l)][c][2] for l in lam]) for c in COMPS}

    master = master_line_list(intensity_model("127I2"), 11000.0, 20100.0, T_range=(200.0, 600.0), S_min=1e-27)
    cont = continuum_model("127I2", master.upper_cut, T_max=600.0)

    lam_air_vac = air_to_vacuum(lam)
    mod = {t: model_components(master, cont, lam_air_vac, T) for t, T in TEMPS.items()}
    mod_vac = {t: model_components(master, cont, lam, T) for t, T in TEMPS.items()}
    mod_alt = {w: model_components(master, cont, lam_air_vac, TEMPS["35"], w) for w in SMOOTH_ALT}

    regions = {"400-495": (400, 495), "500-650": (500, 650), "655-725": (655, 725), "730-850": (730, 850)}
    summary = {}
    for c in COMPS:
        summary[c] = {}
        for rname, (lo, hi) in regions.items():
            m = (lam >= lo) & (lam <= hi) & (ref[c]["35"] > 1.0)
            if not m.any():
                continue
            d = {}
            for t in TEMPS:
                r = mod[t][c][m] / ref[c][t][m]
                rv = mod_vac[t][c][m] / ref[c][t][m]
                d[t] = dict(ratio_min=r.min(), ratio_max=r.max(), ratio_mean=r.mean(),
                            ratio_mean_if_vacuum=rv.mean(),
                            diff_eps_mean=float(np.mean(mod[t][c][m] - ref[c][t][m])),
                            diff_sigma_cm2_mean=float(np.mean(mod[t][c][m] - ref[c][t][m]) * EPS))
            z = (mod["35"][c][m] - ref[c]["35"][m]) / err[c][m]
            d["z35_min"], d["z35_max"], d["z35_rms"] = z.min(), z.max(), float(np.sqrt(np.mean(z**2)))
            d["n"] = int(m.sum())
            summary[c][rname] = {k: (float(v) if not isinstance(v, dict) else {kk: float(vv) for kk, vv in v.items()})
                                 for k, v in d.items()}

    # measured ε (Table IS) against the model total, lines + continuum at 2 nm (weak-absorption limit)
    lam_m = meas[:, 0]
    lv = air_to_vacuum(lam_m)
    tab_is = {}
    for col, T in ((1, 308.55), (3, 337.15)):
        lines = master.at(T, S_min=0.0)
        edges = np.arange(1e7 / 870.0, 1e7 / 390.0, 0.05)
        cen = (edges[1:] + edges[:-1]) / 2
        sig = np.histogram(lines.nu, edges, weights=lines.S)[0] / 0.05 + cont.cross_section(cen, T)
        mt = gauss_avg(cen, sig, lv, 2.0) / EPS
        o, s = meas[:, col], meas[:, col + 1]
        blk = {}
        for lo, hi in ((420, 500), (600, 650), (650, 750), (750, 850)):
            m = (lam_m >= lo) & (lam_m <= hi)
            z = (mt[m] - o[m]) / s[m]
            blk[f"{lo}-{hi}"] = dict(mean_diff_eps=float(np.mean(mt[m] - o[m])), mean_ratio=float(np.mean(mt[m] / o[m])),
                                     z_mean=float(z.mean()), z_rms=float(np.sqrt(np.mean(z**2))))
        tab_is[f"{T:.2f} K"] = dict(blocks=blk, lam=lam_m.tolist(), meas=o.tolist(), sd=s.tolist(), model=mt.tolist())

    out = dict(
        source=str(SRC.relative_to(ROOT)), eps_to_sigma_cm2=EPS, temperatures_K=TEMPS,
        wavelengths="Table IIS read as air, converted to vacuum; *_if_vacuum reads them as vacuum",
        smooth_fwhm_nm_banded_B=SMOOTH_NM, summary=summary,
        table=dict(lam_nm=lam.tolist(), lam_vac_nm=lam_air_vac.tolist(),
                   tellinghuisen={c: {t: ref[c][t].tolist() for t in TEMPS} | {"se35": err[c].tolist()} for c in COMPS},
                   model={t: {k: v.tolist() for k, v in mod[t].items()} for t in TEMPS},
                   model_vacuum_reading={t: {k: v.tolist() for k, v in mod_vac[t].items()} for t in TEMPS},
                   model35_B_smoothing={f"{w} nm": mod_alt[w]["B"].tolist() for w in SMOOTH_ALT}),
        table_IS_vs_model_total=tab_is,
    )
    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    OUT_JSON.write_text(json.dumps(out, indent=1, default=lambda x: None if np.isnan(x) else x))

    # printed table
    print(f"{'λ':>4} | " + " | ".join(f"{c:>5}: T35  model  r    z " for c in COMPS))
    for i, l in enumerate(lam):
        cells = []
        for c in COMPS:
            t, m, e = ref[c]["35"][i], mod["35"][c][i], err[c][i]
            r = m / t if t > 0.05 else np.nan
            z = (m - t) / e if e > 0 else np.nan
            cells.append(f"{t:7.2f} {m:7.2f} {r:5.3f} {z:+6.1f}")
        print(f"{l:4.0f} | " + " | ".join(cells))
    print(json.dumps(summary, indent=1))
    for k, v in tab_is.items():
        print(k, json.dumps(v["blocks"]))
    for w in SMOOTH_ALT:
        m = (lam >= 505) & (lam <= 640)
        print(f"B smoothing {w} nm, 505-640 nm ratio 35°C:", np.round(mod_alt[w]["B"][m] / ref["B"]["35"][m], 3).tolist())

    # figure
    col = {"0": "#0072B2", "35": "#D55E00"}
    fig, axs = plt.subplots(2, 4, figsize=(12, 5.6), sharex="col", gridspec_kw=dict(height_ratios=[1.6, 1]))
    titles = dict(A="A←X", B="B←X (≥500 nm: lines + continuum, 4 nm avg)", C="C←X", total="Total")
    xr = dict(A=(480, 855), B=(400, 700), C=(400, 700), total=(400, 855))
    for j, c in enumerate(COMPS):
        ax, axr = axs[0, j], axs[1, j]
        for t in TEMPS:
            ok = ref[c][t] > 0
            ax.plot(lam[ok], ref[c][t][ok] * EPS, "o", ms=3.5, mfc="white", mec=col[t], mew=0.9,
                    label=f"Tellinghuisen {t} °C")
            ax.plot(lam, mod[t][c] * EPS, "-", color=col[t], lw=1.2, label=f"i2spec {t} °C")
            good = ref[c][t] > 0.5
            axr.plot(lam[good], mod[t][c][good] / ref[c][t][good], "-o", ms=2.5, color=col[t], lw=1)
        good = ref[c]["35"] > 0.5
        band = err[c][good] / ref[c]["35"][good]
        axr.fill_between(lam[good], 1 - band, 1 + band, color="0.8", lw=0, label="±1 LS s.e. (35 °C)")
        axr.axhline(1, color="0.4", lw=0.6)
        if c in ("B", "total"):
            for a in (ax, axr):
                a.axvspan(500, 650, color="0.93", lw=0, zorder=0)
        ax.set_title(titles[c], fontsize=9)
        ax.set_xlim(*xr[c])
        axr.set_ylim(0.9, 1.15)
        axr.set_xlabel("λ, air (nm)")
        ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    axs[0, 0].set_ylabel("σ (cm²)")
    axs[1, 0].set_ylabel("i2spec / Tellinghuisen")
    axs[0, 3].legend(fontsize=7, frameon=False)
    axs[1, 3].legend(fontsize=7, frameon=False, loc="lower right")
    fig.suptitle("i2spec (i2spec2026m) against Tellinghuisen 2011 Table IIS component bands; "
                 "grey band 500–650 nm: his pseudocontinuum B", fontsize=9)
    fig.tight_layout()
    OUT_PNG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PNG, dpi=150)
    print("wrote", OUT_JSON, OUT_PNG)


if __name__ == "__main__":
    main()
