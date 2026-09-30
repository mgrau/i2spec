"""Figure for the Salami & Ross atlas scan (prototypes/scan_salami_ross.py).

Usage:  uv run --with matplotlib python prototypes/plot_scan_salami_ross.py
Output: docs/figures/salami_ross_scan.png
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

SURFACE, INK, INK2, GRID, AXIS = "#fcfcfb", "#1f1e1c", "#52514e", "#e1e0d9", "#c3c2b7"
SEGMENT_COLOURS = {293.15: ("#2a78d6", "cell 20 °C"), 323.15: ("#eb6834", "cell 50 °C"), 463.15: ("#1baf7a", "cell 190 °C")}
ROOT = Path(__file__).resolve().parents[1]


def main():
    d = np.loadtxt(ROOT / "prototypes/out/scan_salami_ross.txt")
    centre, T, rms, noise, cold, shift = d[:, 0], d[:, 1], d[:, 4], d[:, 5], d[:, 7], d[:, 8]
    best = np.array([rms[k] <= rms[centre == centre[k]].min() for k in range(len(d))])
    fig, axes = plt.subplots(3, 1, figsize=(11, 9), sharex=True, facecolor=SURFACE)
    panels = [(cold, "Implied I₂ cold point (°C)", (0, 50)), (shift, "Line shift, model − atlas (MHz)", (-150, 150)),
              (rms / noise, "RMS residual / atlas noise", (0, 22))]
    for ax, (y, label, (lo, hi)) in zip(axes, panels):
        clipped = (y < lo) | (y > hi)
        y = np.clip(y, lo, hi)
        ax.set_ylim(lo - 0.03 * (hi - lo), hi + 0.03 * (hi - lo))
        ax.set_facecolor(SURFACE)
        ax.grid(True, color=GRID, lw=0.75)
        ax.set_axisbelow(True)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        for side in ("left", "bottom"):
            ax.spines[side].set_color(AXIS)
        ax.tick_params(colors=AXIS, labelcolor=INK2)
        ax.set_ylabel(label, color=INK2)
        for temp, (colour, name) in SEGMENT_COLOURS.items():
            m = (T == temp) & ~clipped
            ax.plot(centre[m & best], y[m & best], "o", ms=5, color=colour, label=name)
            ax.plot(centre[m & ~best], y[m & ~best], "o", ms=5, mfc="none", mec=colour, mew=1.2)
            m = (T == temp) & clipped & best
            ax.plot(centre[m], y[m], "^" if hi == y[m].max(initial=lo) else "v", ms=6, color=colour)
        if clipped.any():
            ax.text(0.995, 0.97, "▲▼ off scale", transform=ax.transAxes, ha="right", va="top", color=INK2, fontsize=8)
    axes[1].axhspan(-90, 90, color=GRID, alpha=0.6, lw=0, zorder=0)
    axes[1].text(0.005, 0.03, "grey band: atlas calibration ±0.003 cm⁻¹", transform=axes[1].transAxes, color=INK2,
                 fontsize=8)
    axes[2].axhline(1, color=AXIS, lw=1)
    axes[0].legend(frameon=False, fontsize=8, loc="best", labelcolor=INK2,
                   title="filled: better fit where segments overlap", title_fontsize=8)
    axes[-1].set_xlabel("Wavenumber (cm⁻¹)", color=INK2)
    top = axes[0].secondary_xaxis("top", functions=(lambda x: 1e7 / np.maximum(x, 1), lambda x: 1e7 / np.maximum(x, 1)))
    top.set_xlabel("Vacuum wavelength (nm)", color=INK2)
    top.tick_params(colors=AXIS, labelcolor=INK2)
    fig.tight_layout()
    fig.savefig(ROOT / "docs/figures/salami_ross_scan.png", dpi=120, facecolor=SURFACE)


if __name__ == "__main__":
    main()
