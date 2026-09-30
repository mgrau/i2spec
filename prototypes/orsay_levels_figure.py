"""Figure: the X levels v'' = 18-25 from the Orsay atlas against the potentials. -> docs/figures/orsay_x_levels.png"""
import csv, json, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
d = json.load(open("src/i2spec/data/level_corrections_2026f.json"))
A = list(csv.DictReader(open("data/atlas_lines/orsay1982_part1_assigned.csv")))
MHZ = 29979.2458; ppb = d["orsay"]["calibration_scale_ppb"]
fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
vs = list(range(13, 26))
before = [np.median([float(a["obs_minus_model_MHz"]) - ppb*1e-9*float(a["sigma_cm1"])*MHZ for a in A if int(a["v_lower"]) == v]) for v in vs]
n = [sum(int(a["v_lower"]) == v for a in A) for v in vs]
ax[0].bar(vs, before, color=["#999999" if v <= 17 else "#4C72B0" for v in vs])
for v, b, k in zip(vs, before, n):
    ax[0].text(v, b + (2 if b >= 0 else -6), str(k), ha="center", fontsize=7)
ax[0].axhline(0, color="k", lw=0.8); ax[0].axvline(17.5, color="k", ls="--", lw=1)
ax[0].set_xlabel("v$''$"); ax[0].set_ylabel("median atlas − i2spec2026f (MHz)")
ax[0].set_title("Assigned lines per level (count on bar), after the\n+24 ppb atlas scale: grey = held, blue = fitted", fontsize=9)
v2 = [v for v in range(18, 26)]
y100 = [np.polyval(d["X"][str(v)][::-1], 1.0) for v in v2]
err = [d["held_out_MHz"]["X"][str(v)] for v in v2]
ax[1].errorbar(v2, y100, yerr=err, fmt="o", color="#4C72B0", capsize=3)
ax[1].axhline(0, color="k", lw=0.8)
ax[1].set_xlabel("v$''$"); ax[1].set_ylabel("correction to X term value at J$''$=100 (MHz)")
ax[1].set_title("X levels measured above the published range\n(error bars: level uncertainty; previously 300 MHz)", fontsize=9)
fig.tight_layout(); fig.savefig("docs/figures/orsay_x_levels.png", dpi=150)
