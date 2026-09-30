"""Plot the Salami & Ross atlas against the i2spec model, from the output of compare_salami_ross.py.

Usage:  uv run --with matplotlib python prototypes/plot_salami_ross.py STEM OUT.png [T_CELL_C] [--fixed-offset]
        (STEM e.g. prototypes/out/salami_ross_18799_18806; use --fixed-offset when the fit was run with it)
"""

import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import FuncFormatter, MultipleLocator

# Reference palette (dataviz skill): light surface, ink tokens, categorical slots 1 and 2.
SURFACE, INK, INK2, MUTED, GRID, AXIS = "#fcfcfb", "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"
ATLAS_COLOR, MODEL_COLOR = "#2a78d6", "#eb6834"

args = [a for a in sys.argv[1:] if not a.startswith("--")]
stem, out = args[0], args[1]
t_cell = float(args[2]) if len(args) > 2 else 20.0
fixed_offset = "--fixed-offset" in sys.argv
nu, obs, mod, res = np.loadtxt(f"{stem}_hyperfine.txt").T
labels = [line.rstrip("\n").split("\t") for line in open(f"{stem}_lines.txt")]

# Full 0-100% scale when lines are deep; otherwise zoom to the data so weak features stay visible.
lo_t, hi_t = 100 * min(obs.min(), mod.min()), 100 * max(obs.max(), mod.max())
if hi_t - lo_t > 40:
    y_lo, y_hi = 0.0, 100.0
else:
    pad = 0.25 * (hi_t - lo_t)
    y_lo, y_hi = max(0.0, np.floor(lo_t - pad)), min(100.0, np.ceil(hi_t + pad))
r_lim = float(np.ceil(110 * max(4 * res.std(), np.abs(res).max())))  # percent, 10% headroom

plt.rcParams.update({"font.family": ["Helvetica Neue", "Helvetica", "Arial", "DejaVu Sans"], "font.size": 9})
fig, (ax, axr) = plt.subplots(2, 1, figsize=(11, 6.0), sharex=True, facecolor=SURFACE,
                              gridspec_kw=dict(height_ratios=[3, 1], hspace=0.08))
for a in (ax, axr):
    a.set_facecolor(SURFACE)
    a.grid(True, color=GRID, lw=0.75)
    a.set_axisbelow(True)
    for side in ("top", "right"):
        a.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        a.spines[side].set_color(AXIS)
    a.tick_params(colors=AXIS, labelcolor=INK2)

ax.plot(nu, 100 * obs, color=ATLAS_COLOR, lw=1.5, solid_joinstyle="round", label="Salami & Ross (2005) FTS atlas")
ax.plot(nu, 100 * mod, color=MODEL_COLOR, lw=1.5, solid_joinstyle="round", label="i2spec model")
ax.set_ylabel("Transmission (%)", color=INK2)
ax.set_ylim(y_lo, y_hi)
ax.legend(frameon=False, loc="lower left", labelcolor=INK, ncols=2, bbox_to_anchor=(0, 1.0), borderaxespad=0.2)

span = y_hi - y_lo
last = -np.inf
for x, text in ((float(a), b) for a, b in labels):
    if x - last > 0.6:
        ax.plot([x, x], [y_hi - 0.07 * span, y_hi - 0.045 * span], color=MUTED, lw=0.75)
        ax.annotate(text, xy=(x, y_hi - 0.04 * span), ha="center", va="bottom", color=MUTED, fontsize=8)
        last = x

axr.plot(nu, 100 * res, color=INK2, lw=1.0)
axr.axhline(0, color=AXIS, lw=0.75)
axr.set_ylabel("Atlas − model (%)", color=INK2)
axr.set_xlabel("Wavenumber (cm⁻¹)", color=INK2)
axr.set_ylim(-r_lim, r_lim)
axr.xaxis.set_major_locator(MultipleLocator(1))
axr.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:,.0f}"))

lam = 1e7 / nu.mean()
offset_note = "zero offset fixed at 0" if fixed_offset else "zero offset"
ax.set_title(f"I₂ B–X absorption near {lam:.0f} nm: i2spec model vs. FTS atlas", loc="left", color=INK, fontsize=12, pad=52)
ax.text(0, 1.09, f"50 cm cell at {t_cell:.0f} °C, 0.02 cm⁻¹ resolution. Model: Hannover 2008 potentials, Bodermann 2002 hyperfine, "
        f"Tellinghuisen 2011 transition moment. Fitted: column density, baseline, shift, instrument width; {offset_note}. "
        f"RMS residual {100 * res.std():.2f}% of transmission.", transform=ax.transAxes, color=INK2, fontsize=8, wrap=True)
fig.savefig(out, dpi=160, facecolor=SURFACE, bbox_inches="tight")
print("wrote", out)
