"""Stage 2 of the refit (src/i2spec/hyperfine_fit.py): every number in docs/design/hyperfine-fit.md.

Usage:  uv run --with matplotlib python prototypes/hyperfine_stage2.py
Output: prototypes/out/hyperfine_stage2.json, docs/figures/hyperfine_stage2.png   (about 2 minutes)
"""

import json
import time
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from i2spec.hyperfine_fit import HyperfineFit, rms  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
SURFACE, INK, INK2, GRID, AXIS = "#fcfcfb", "#1f1e1c", "#52514e", "#e1e0d9", "#c3c2b7"
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"

C1 = ["C_b", "C_b_s", "C_b_y"]
C2 = C1 + ["C_b_s2", "C_b_sy", "C_b_y2"]
SS = ["d_b", "d_b_s", "d_b_y", "delta_b", "delta_b_s", "delta_b_y"]
QB = ["eqQ_b", "eqQ_b_s", "eqQ_b_y"]
QX = ["eqQ_x", "eqQ_x_s", "C_x", "C_x_s"]
MODELS = {"published": [], "C_B linear": C1, "C_B quadratic": C2, "C_B quadratic + spin-spin": C2 + SS,
          "+ eqQ_B": C2 + SS + QB, "everything": C2 + SS + QB + QX, "spin-spin only": SS}
CHOSEN = "C_B quadratic + spin-spin"


def label(key):
    iso, br, J, vu, vl = key
    return f"{'' if iso == '127I2' else iso + ' '}{br}({J}) {vu}-{vl}"


def main():
    t0 = time.perf_counter()
    out = {}
    h = HyperfineFit(model_floor=0.005)
    dev0 = h.deviations(h.x0)
    ds = np.array([d.dataset for d in h.data])
    iso = np.array([d.isotopologue for d in h.data])
    precise = (h.sigma <= 0.025) & (iso == "127I2")

    # 1. who owns chi-squared when weighted by measurement uncertainty alone
    chi2 = (dev0 / h.sigma) ** 2
    share = defaultdict(float)
    for d, c in zip(h.data, chi2):
        share[f"{d.dataset} {label(h.line_of(d))}"] += c / chi2.sum()
    out["chi2_share"] = dict(sorted(share.items(), key=lambda kv: -kv[1])[:5])

    # 2. per-line fits: is the Hamiltonian right, and which parameter carries the error
    per_line = []
    for key in h.lines:
        n = sum(h.line_of(d) == key for d in h.data)
        if n < 8:
            continue
        four = h.per_line_fit(key)
        c_only = h.per_line_fit(key, terms=("C_b",))
        per_line.append(dict(line=label(key), key=list(key), J_upper=key[2] + (1 if key[1] == "R" else -1),
                             dataset=next(d.dataset for d in h.data if h.line_of(d) == key), **four,
                             rms_C_only=c_only["rms_after"], normalised_C_only=c_only["normalised_after"]))
    out["per_line"] = per_line

    # 3. cross-validated correction models, linearised
    J = h.jacobian()
    lines, vgroups = h.groups("line"), h.groups("v_upper")
    cv = {}
    for name, terms in MODELS.items():
        coef, fit = h.linear_fit(terms, J, dev0) if terms else ({}, dev0)
        lo = h.cross_validate(terms, lines, J, dev0)
        lv = h.cross_validate(terms, vgroups, J, dev0)
        cv[name] = dict(terms=terms, in_sample=rms(fit[precise]) * 1e3, held_out_line=rms(lo[precise]) * 1e3,
                        held_out_v=rms(lv[precise]) * 1e3,
                        by_dataset={k: rms(lo[precise & (ds == k)]) * 1e3 for k in sorted(set(ds[precise]))})
        if name == CHOSEN:
            out["chosen_coefficients"] = coef
    out["cross_validation"] = cv

    # the same with the 25 kHz floor, to show the conclusion does not hinge on the weighting
    h25 = HyperfineFit(model_floor=0.025)
    out["cross_validation_25kHz"] = {name: dict(held_out_line=rms(h25.cross_validate(t, lines, J, dev0)[precise]) * 1e3)
                                     for name, t in MODELS.items()}

    # 4. within the 532 nm lines: predict each v' from the other five
    vu = np.array([d.v_upper for d in h.data])
    region = ds == "bipm2012a"
    rows = []
    for g in sorted(set(vu[region])):
        test, train = region & (vu == g), region & (vu != g)
        row = dict(v_upper=int(g), n=int(test.sum()), published=rms(dev0[test]) * 1e3)
        for name in ("C_B quadratic", CHOSEN):
            cols = [h.free.index(t) for t in MODELS[name]]
            pred = dev0[test] + J[np.ix_(test, cols)] @ h._ridge(J, dev0, cols, train)
            row[name] = rms(pred) * 1e3
        rows.append(row)
    out["bipm2012a_leave_one_v_out"] = rows
    out["seconds"] = time.perf_counter() - t0

    (ROOT / "prototypes/out").mkdir(exist_ok=True)
    (ROOT / "prototypes/out/hyperfine_stage2.json").write_text(json.dumps(out, indent=1, default=float))
    figure(out)
    print(f"done in {out['seconds']:.0f} s")


def style(ax):
    ax.set_facecolor(SURFACE)
    ax.grid(True, color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(AXIS)
    ax.tick_params(colors=INK2, labelsize=8.5)


def figure(out):
    fig, (a, b, c) = plt.subplots(3, 1, figsize=(11, 12.5), facecolor=SURFACE,
                                  gridspec_kw=dict(height_ratios=[1.15, 1, 0.85]))
    precise = sorted((r for r in out["per_line"] if r["median_uncertainty"] <= 25 and r["key"][0] == "127I2"),
                     key=lambda r: (r["key"][3], r["J_upper"]))

    # (a) per-line normalised rms: published, C_B alone, all four B parameters
    style(a)
    x = np.arange(len(precise))
    a.axhline(1.0, color=INK2, lw=1.0, ls=(0, (4, 3)), zorder=2)
    a.text(len(precise) - 0.4, 1.12, "measurement noise", ha="right", va="bottom", fontsize=8, color=INK2)
    a.scatter(x, [r["normalised_before"] for r in precise], s=34, marker="o", color=BLUE, label="published formulae", zorder=3)
    a.scatter(x, [r["normalised_C_only"] for r in precise], s=34, marker="s", color=ORANGE, label="C_B freed for this line", zorder=3)
    a.scatter(x, [r["normalised_after"] for r in precise], s=40, marker="^", color=AQUA, label="eqQ_B, C_B, d_B, δ_B freed", zorder=3)
    a.set_yscale("log")
    a.set_xticks(x, [r["line"] for r in precise], rotation=60, ha="right", fontsize=8)
    a.set_ylabel("rms of (model − measured) / σ", color=INK)
    a.set_title("(a) The Hamiltonian is right: per-line parameters reach the measurement noise", loc="left", color=INK, fontsize=11)
    a.legend(frameon=False, fontsize=8.5, loc="upper left", ncol=3, labelcolor=INK)

    # (b) the fitted C_B shift per line, against J'
    style(b)
    clusters, labelled = defaultdict(list), set()
    for r in precise:
        bipm = r["dataset"] == "bipm2012a"
        b.errorbar(r["J_upper"], r["x"]["C_b"], yerr=r["error"]["C_b"], fmt="o" if bipm else "D", ms=6,
                   color=BLUE if bipm else ORANGE, ecolor=INK2, elinewidth=0.8, zorder=3)
        if bipm:
            clusters[r["key"][3]].append((r["J_upper"], r["x"]["C_b"]))
        elif r["key"][3] not in labelled:   # P(62) 17-1 sits inside the v' = 32 cluster: label it below
            labelled.add(r["key"][3])
            below = r["key"][3] == 17
            b.annotate(f"v′={r['key'][3]}", (r["J_upper"], r["x"]["C_b"]), textcoords="offset points",
                       xytext=(6, -12 if below else 5), fontsize=7.5, color=INK2)
    for v, pts in clusters.items():   # one label per v' group of the 532 nm lines; v' = 34 goes below, clear of v' = 28
        J_mid = np.mean([q[0] for q in pts])
        y, dy = (min(q[1] for q in pts), -16) if v == 34 else (max(q[1] for q in pts), 9)
        b.annotate(f"v′={v}", (J_mid, y), textcoords="offset points", xytext=(0, dy), ha="center",
                   fontsize=7.5, color=INK2)
    b.axhline(0, color=AXIS, lw=0.8)
    b.scatter([], [], marker="o", color=BLUE, label="BIPM 532 nm lines (bipm2012a)")
    b.scatter([], [], marker="D", color=ORANGE, label="other data sets")
    b.legend(frameon=False, fontsize=8.5, loc="upper left", labelcolor=INK)
    b.set_xlabel("J′", color=INK)
    b.set_ylabel("ΔC_B, fitted − published (kHz)", color=INK)
    b.set_title("(b) The error is in C_B; within the 532 nm lines v′ and J′ rise together", loc="left", color=INK, fontsize=11)

    # (c) predicting a held-out v' within the 532 nm lines
    style(c)
    rows = out["bipm2012a_leave_one_v_out"]
    x = np.arange(len(rows))
    wbar = 0.38
    c.bar(x - wbar / 2 - 0.01, [r["published"] for r in rows], wbar, color=BLUE, label="published formulae", zorder=3)
    c.bar(x + wbar / 2 + 0.01, [r[CHOSEN] for r in rows], wbar, color=ORANGE,
          label="correction fitted to the other five v′ (C_B quadratic + spin-spin)", zorder=3)
    for i, r in enumerate(rows):
        for dx, v in ((-wbar / 2 - 0.01, r["published"]), (wbar / 2 + 0.01, r[CHOSEN])):
            c.text(i + dx, v * 1.08, f"{v:.0f}", ha="center", va="bottom", fontsize=8, color=INK2)
    c.set_yscale("log")
    c.set_ylim(10, 3000)
    c.set_xticks(x, [f"v′ = {r['v_upper']}\n({r['n']} splittings)" for r in rows], fontsize=8.5)
    c.set_ylabel("held-out rms (kHz)", color=INK)
    c.set_title("(c) …but a correction does not predict an unseen v′: better for 34, 36, 37, worse for 32, 33, 35",
                loc="left", color=INK, fontsize=11)
    c.legend(frameon=False, fontsize=8.5, loc="upper left", labelcolor=INK)

    fig.tight_layout(h_pad=2.0)
    fig.savefig(ROOT / "docs/figures/hyperfine_stage2.png", dpi=130, facecolor=SURFACE)


if __name__ == "__main__":
    main()
