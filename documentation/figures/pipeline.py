"""The calculation at a glance, for documentation/calculation.md: six steps, each with a small picture drawn
from the model, and what each step takes in.

    uv run --group research python documentation/figures/pipeline.py

Writes documentation/figures/pipeline.svg with the same placeholder colours as make_figures.py, so it follows
the site's light and dark themes.
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("svg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
from matplotlib.path import Path as MPath

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from i2spec.constants import DEFAULT_PARAMETERS  # noqa: E402
from i2spec.level_corrections import corrections_for  # noqa: E402
from i2spec.lookup import Catalog  # noqa: E402
from i2spec.model import RovibronicModel  # noqa: E402
from i2spec.potentials import load_extended, load_potentials  # noqa: E402
from i2spec.spectrum import cross_section, number_density, transmission  # noqa: E402

# the helpers and colour tokens of make_figures.py, without running it
C = dict(ink="#010101", ink2="#020202", ink3="#030303", rule="#040404", trace="#050505", pick="#060606",
         mark="#070707", surface="#080808")
TOKENS = {C["ink"]: "var(--ink)", C["ink2"]: "var(--ink-2)", C["ink3"]: "var(--ink-3)", C["rule"]: "var(--rule)",
          C["trace"]: "var(--trace)", C["pick"]: "var(--pick)", C["mark"]: "var(--mark)", C["surface"]: "var(--surface)"}
plt.rcParams.update({
    "svg.fonttype": "none", "font.family": "sans-serif", "font.size": 9,
    "axes.edgecolor": C["rule"], "axes.labelcolor": C["ink3"], "xtick.color": C["ink3"], "ytick.color": C["ink3"],
    "text.color": C["ink"], "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": 0.7,
    "xtick.major.size": 2, "ytick.major.size": 2, "legend.frameon": False, "figure.facecolor": "none",
    "axes.facecolor": "none", "savefig.facecolor": "none",
})


def save(fig, name):
    import re
    path = HERE / f"{name}.svg"
    fig.savefig(path, pad_inches=0.02)
    plt.close(fig)
    svg = path.read_text()
    for hexcode, var in TOKENS.items():
        svg = svg.replace(hexcode, var)
    svg = re.sub(r"font-family:[^;\"]*", "font-family:var(--body)", svg)
    svg = re.sub(r"font: ([\d.]+)px '[^']*'", r"font-size:\1px;font-family:var(--body)", svg)
    svg = re.sub(r'<svg([^>]*?) width="[^"]*" height="[^"]*"', r'<svg\1', svg, count=1)
    svg = re.sub(r"<metadata>.*?</metadata>", "", svg, flags=re.S)
    path.write_text(svg)
    print(f"{path.name}: {len(svg) / 1e3:.0f} kB")


def turning(V, R, E):
    inside = np.flatnonzero(V(R) <= E)
    return R[inside[0]], R[inside[-1]]


model = RovibronicModel("127I2")
ext, from_v = load_extended()
pub = load_potentials()
corr = corrections_for(DEFAULT_PARAMETERS)

# ---------------------------------------------------------------------------------------------- layout
W, H = 7.4, 5.6
fig = plt.figure(figsize=(W, H))
CARD_W, CARD_H = 0.295, 0.335
XS = [0.025, 0.3525, 0.68]
YS = [0.585, 0.135]
STEPS = [
    ("1", "Level energies", "radial equation in V(R)"),
    ("2", "Level corrections", "B levels, fitted to precision lines"),
    ("3", "Line positions", "ν = E′ − E″, with a 1σ each"),
    ("4", "Hyperfine components", "15 or 21 per line"),
    ("5", "Line strengths", "⟨v′|μₑ(R)|v″⟩², per band"),
    ("6", "Spectrum", "at T, for a given cell"),
]
INPUTS = {
    0: "X and B potentials: 2008 curves,\nMLR curves beyond their range",
    1: "comb lines → per-level polynomials;\na Gaussian process fills the rest",
    3: "hyperfine formulae + a table of\nmeasured B-state parameters",
    4: "transition moment μₑ(R)\n(Tellinghuisen 2011)",
    5: "temperature, path length,\ncold-finger temperature",
}
plot_axes = []
for k, (num, title, sub) in enumerate(STEPS):
    x, y = XS[k % 3], YS[k // 3]
    fig.patches.append(FancyBboxPatch((x, y), CARD_W, CARD_H, boxstyle="round,pad=0.004,rounding_size=0.012",
                                      transform=fig.transFigure, fc="none", ec=C["rule"], lw=0.9))
    fig.text(x + 0.012, y + CARD_H - 0.012, num, fontsize=14, fontweight="bold", color=C["trace"], va="top")
    fig.text(x + 0.045, y + CARD_H - 0.014, title, fontsize=9.5, fontweight="bold", color=C["ink"], va="top")
    fig.text(x + 0.045, y + CARD_H - 0.050, sub, fontsize=7.5, color=C["ink3"], va="top")
    pad = 0.035 if k in (1, 2) else 0.0
    plot_axes.append(fig.add_axes([x + 0.035 + pad, y + 0.075, CARD_W - 0.05 - pad, CARD_H - 0.16]))
    if k in INPUTS:
        fig.text(x + CARD_W / 2, y - 0.014, INPUTS[k], fontsize=7, color=C["ink2"], ha="center", va="top",
                 linespacing=1.15, style="italic")
# arrows: along each row, and from the end of the first row to the start of the second
for row in (0, 1):
    for i in (0, 1):
        y = YS[row] + CARD_H / 2
        fig.patches.append(FancyArrowPatch((XS[i] + CARD_W + 0.004, y), (XS[i + 1] - 0.004, y),
                                           transform=fig.transFigure, arrowstyle="-|>", mutation_scale=9,
                                           color=C["ink2"], lw=1.0))
# from the right side of step 3, round the gap between the rows, into the left side of step 4: straight runs
# with rounded right-angle bends, drawn in inches so the bends are circular
def elbow(points, r=0.09):
    verts, codes = [points[0]], [MPath.MOVETO]
    for a, b, c in zip(points, points[1:], points[2:]):
        a, b, c = map(np.asarray, (a, b, c))
        u, w = (a - b) / np.linalg.norm(a - b), (c - b) / np.linalg.norm(c - b)
        verts += [tuple(b + r * u), tuple(b), tuple(b + r * w)]
        codes += [MPath.LINETO, MPath.CURVE3, MPath.CURVE3]
    verts.append(points[-1]); codes.append(MPath.LINETO)
    return MPath(verts, codes)


x_out, x_in = (XS[2] + CARD_W + 0.013) * W, (XS[0] - 0.013) * W
y_top, y_gap, y_bot = (YS[0] + CARD_H / 2) * H, (YS[1] + CARD_H + 0.025) * H, (YS[1] + CARD_H / 2) * H
fig.patches.append(FancyArrowPatch(path=elbow([((XS[2] + CARD_W + 0.004) * W, y_top), (x_out, y_top),
                                               (x_out, y_gap), (x_in, y_gap), (x_in, y_bot),
                                               ((XS[0] - 0.004) * W, y_bot)]),
                                   transform=fig.dpi_scale_trans, arrowstyle="-|>", mutation_scale=9,
                                   color=C["ink2"], lw=1.0))
fig.text(0.5, 0.985, "computed once per isotopologue, independent of temperature: steps 1–5", ha="center",
         va="top", fontsize=8, color=C["ink3"])


def clean(ax, xlabel=None):
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    if xlabel:
        ax.set_xlabel(xlabel, fontsize=7, labelpad=1)
    ax.tick_params(labelsize=6.5, pad=1)


# 1. potentials
ax = plot_axes[0]
R = np.linspace(2.25, 6.5, 800)
for s, col in (("X", C["trace"]), ("B", C["pick"])):
    a0, b0 = turning(pub[s], np.linspace(2.0, 9.0, 20000), model.energy(s, from_v[s] - 1, 0))
    V = np.where((R >= a0) & (R <= b0), pub[s](R), ext[s](R))
    ax.plot(R, V, color=col, lw=1.4)
    for v in range(0, 80 if s == "B" else 100, 6):
        E = model.energy(s, v, 0)
        a, b = turning(lambda r: np.where((r >= a0) & (r <= b0), pub[s](r), ext[s](r)), np.linspace(2.2, 12, 6000), E)
        ax.plot([a, min(b, 6.5)], [E, E], color=col, lw=0.5, alpha=0.6)
ax.text(2.95, 1500, "X", color=C["trace"], fontsize=9, fontweight="bold")
ax.text(3.45, 15800, "B", color=C["pick"], fontsize=9, fontweight="bold")
ax.set_ylim(-500, 21500)
clean(ax, "R (Å)")

# 2. corrections: B levels at J' = 40, from comb data (points) and from the Gaussian process (band)
ax = plot_axes[1]
JS = 40
vs, mid, err = [], [], []
for v in sorted(corr.coefficients["B"]):
    if v > 50:
        continue
    lo, hi = corr.coverage["B"].get(v, (JS, JS))
    j = int(np.clip(JS, lo, hi))
    vs.append(v); mid.append(corr.shift("B", v, j)); err.append(corr.uncertainty_at("B", v, j) or 1.0)
ax.errorbar(vs, mid, err, fmt="o", ms=2.6, color=C["trace"], elinewidth=0.7, capsize=0)
vv = np.arange(1, 36)
ms = np.array([corr.gp(int(v), JS) for v in vv])
ax.fill_between(vv, ms[:, 0] - ms[:, 1], ms[:, 0] + ms[:, 1], color=C["pick"], alpha=0.2, lw=0)
ax.plot(vv, ms[:, 0], color=C["pick"], lw=1.0)
ax.axhline(0, color=C["rule"], lw=0.6)
ax.set_yscale("symlog", linthresh=10)
ax.set_yticks([-100, -10, 0, 10, 100], ["−100", "−10", "0", "10", "100"])
ax.set_xlim(0, 50)
ax.text(0.97, 0.95, "comb data", color=C["trace"], fontsize=6.5, ha="right", va="top", transform=ax.transAxes)
ax.text(0.97, 0.83, "Gaussian process", color=C["pick"], fontsize=6.5, ha="right", va="top", transform=ax.transAxes)
ax.set_ylabel("ΔE (MHz)", fontsize=7, labelpad=1)
ax.set_xlabel("v′  (J′ = 40)", fontsize=7, labelpad=1)
ax.tick_params(labelsize=6.5, pad=1)

# 3. line positions: the 32-0 band near its head (Fortrat)
ax = plot_axes[2]
J = np.arange(1, 81)
for br, col in (("R", C["trace"]), ("P", C["pick"])):
    nu = np.array([model.transition(32, 0, j, br) for j in J])
    ax.plot(nu, J, "o", ms=1.3, color=col)
    i = 30 if br == "P" else 20
    ax.text(nu[i] + (-4 if br == "P" else 4), J[i] - 8 * (br == "R"), br, color=col, fontsize=8,
            fontweight="bold", ha="right" if br == "P" else "left", va="center")
ax.set_ylim(0, 92)
ax.set_ylabel("J″", fontsize=7, labelpad=1)
ax.tick_params(labelsize=6.5, pad=1)
ax.set_xlabel("ν (cm⁻¹), band 32–0", fontsize=7, labelpad=1)
ax.xaxis.get_major_formatter().set_useOffset(False)
ax.locator_params(axis="x", nbins=3)

# 4. hyperfine: R(56) 32-0
ax = plot_axes[3]
nu0, comps = model.hyperfine_components(32, 0, 56, "R")
main = [c for c in comps if c.label]
off = np.array([c.offset for c in main]); st = np.array([c.strength for c in main])
x = np.linspace(-1000, 900, 1200)
prof = sum(s * np.exp(-0.5 * ((x - o) / (434 / 2.3548)) ** 2) for o, s in zip(off, st))
ax.fill_between(x, 0, prof / prof.max(), color=C["ink3"], alpha=0.12, lw=0)
ax.plot(x, prof / prof.max(), color=C["ink3"], lw=0.8)
ax.vlines(off, 0, st / st.max() * 0.85, color=C["trace"], lw=1.2)
ax.text(-950, 0.95, "R(56) 32–0", fontsize=7, color=C["ink2"])
ax.set_ylim(0, 1.1)
clean(ax, "offset (MHz)")

# 5. line strengths: band strengths of the v'' = 0 progression at 300 K
ax = plot_axes[4]
catalog = Catalog(temperature=300.0)
master = catalog.master("127I2")
L = master.at(300.0, 14000.0, 20100.0, S_min=1e-24)
vu, vl, S = np.asarray(L.v_upper), np.asarray(L.v_lower), np.asarray(L.S)
for v2, col, al in ((0, C["trace"], 1.0), (1, C["pick"], 0.8), (2, C["ink3"], 0.7)):
    band = np.array([S[(vu == v) & (vl == v2)].sum() for v in range(0, 60)])
    ax.bar(np.arange(60) + (v2 - 1) * 0.28, band / 1e-17, width=0.3, color=col, alpha=al, lw=0, label=f"v″ = {v2}")
ax.legend(fontsize=6.5, loc="upper right", handlelength=0.8, borderaxespad=0.1)
ax.set_xlim(0, 60)
clean(ax, "v′")

# 6. spectrum: cross section and 10 cm cell transmission at 532 nm
ax = plot_axes[5]
lines = master.at(293.15, 18787.7, 18789.0, S_min=1e-26)
grid = np.linspace(18787.75, 18789.0, 1500)
sigma = cross_section(grid, lines, 293.15)
tr = transmission(sigma, number_density(293.15, 293.15) * 10.0)
lam = 1e7 / grid
ax.fill_between(lam, 0, sigma / sigma.max(), color=C["trace"], alpha=0.15, lw=0)
ax.plot(lam, sigma / sigma.max(), color=C["trace"], lw=1.0)
ax.plot(lam, 1.05 + 0.6 * (tr - tr.min()) / (1 - tr.min()), color=C["pick"], lw=1.0)
ax.text(lam.max(), 0.85, "σ(ν)", color=C["trace"], fontsize=7, ha="left")
ax.text(lam.max(), 1.92, "transmission, 10 cm", color=C["pick"], fontsize=7, ha="left", va="top")
ax.invert_xaxis()
ax.xaxis.get_major_formatter().set_useOffset(False)
ax.locator_params(axis="x", nbins=3)
ax.set_ylim(0, 1.95)
clean(ax, "λ (nm), near 532 nm")

save(fig, "pipeline")
