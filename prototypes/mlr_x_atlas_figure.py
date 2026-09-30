"""Figure for the atlas refit of the X potential: what moved, and the held-out test. -> docs/figures/mlr_x_atlas.png"""
import json, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from i2spec.model import RovibronicModel
from i2spec.constants import MHZ_PER_CM
from i2spec.level_corrections import load_level_corrections
g, h = RovibronicModel("127I2", "i2spec2026g"), RovibronicModel("127I2", "i2spec2026h")
C = load_level_corrections("level_corrections_2026f")
vs = np.arange(18, 48)
fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
for J, ls in ((60, "-"), (150, "--")):
    ax[0].plot(vs, [(h.energy("X", v, J) - g.energy("X", v, J)) * MHZ_PER_CM for v in vs], ls, color="#4C72B0", label=f"J''={J}")
ax[0].axvspan(17.5, 25.5, color="#dddddd", zorder=0, label="atlas-measured")
ax[0].axvline(47.5, color="k", ls=":", lw=1); ax[0].text(46.8, -120, "emission\nv''=48", ha="right", fontsize=8)
ax[0].axhline(0, color="k", lw=0.8)
ax[0].set_xlabel("v''"); ax[0].set_ylabel("i2spec2026h − i2spec2026g (MHz)")
ax[0].set_title("Where the refit X potential moves the levels", fontsize=10); ax[0].legend(fontsize=8, frameon=False)
loo = {18: 7.2, 19: 11.3, 20: 19.0, 21: 17.4, 22: 22.3, 23: 20.6, 24: 11.8, 25: 38.1}
old = {18: 10.6, 19: 16.8, 20: 26.9, 21: 28.7, 22: 28.8, 23: 31.4, 24: 41.8, 25: 81.0}
v2 = sorted(loo)
ax[1].bar(np.array(v2) - 0.2, [old[v] for v in v2], 0.4, color="#bbbbbb", label="mlr_x_2026c (before)")
ax[1].bar(np.array(v2) + 0.2, [loo[v] for v in v2], 0.4, color="#4C72B0", label="refit, level held out")
ax[1].plot(v2, [C.uncertainty("X", v) for v in v2], "k_", ms=18, mew=2, label="level uncertainty")
ax[1].set_xlabel("v''"); ax[1].set_ylabel("rms error on the atlas level (MHz)")
ax[1].set_title("Predicting each atlas level without it", fontsize=10); ax[1].legend(fontsize=8, frameon=False)
fig.tight_layout(); fig.savefig("docs/figures/mlr_x_atlas.png", dpi=150)
