"""The calculation at a glance, for documentation/calculation.md: six steps, each with a small picture drawn
from the model, and what each step takes in.

    uv run --group research python documentation/figures/pipeline.py

Writes documentation/figures/pipeline.svg with the same placeholder colors as make_figures.py, so it follows
the site's light and dark themes.
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("svg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from i2spec import preferred  # noqa: E402
from i2spec.lookup import Catalog  # noqa: E402
from i2spec.spectrum import cross_section, number_density, transmission  # noqa: E402

# the helpers and color tokens of make_figures.py, without running it
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


model = preferred.model()                       # balanced-expanded (research-20261020)
pot = model._levels_x.potentials
table = preferred.preferred()


def energy(state, v, J):
    try:
        return model._energy(state, v, J)
    except ValueError:                          # dissociated at this J
        return np.nan

# ---------------------------------------------------------------------------------------------- layout
W, H = 7.4, 5.6
fig = plt.figure(figsize=(W, H))
CARD_W, CARD_H = 0.295, 0.335
XS = [0.025, 0.3525, 0.68]
YS = [0.585, 0.135]
STEPS = [
    ("1", "Level energies", "radial equation in V(R)"),
    ("2", "Line positions", "ν = E′ − E″, P and R branches"),
    ("3", "Hyperfine components", "15 or 21 per line"),
    ("4", "Line strengths", "⟨v′|μₑ(R)|v″⟩², per band"),
    ("5", "Spectrum", "at T, for a given cell"),
    ("6", "Uncertainties", "1σ for every line and component"),
]
INPUTS = {
    0: "shared MLR potentials, fitted to all\n¹²⁷I₂ data and continued to both limits",
    2: "radial hyperfine functions; second-\norder sums beyond v′ = 43",
    3: "transition moment μₑ(R)\n(Tellinghuisen 2011)",
    4: "temperature, path length,\ncold-finger temperature",
    5: "fit covariance; errors of held-out\nlevels and hyperfine intervals",
}
plot_axes = []
for k, (num, title, sub) in enumerate(STEPS):
    x, y = XS[k % 3], YS[k // 3]
    fig.patches.append(FancyBboxPatch((x, y), CARD_W, CARD_H, boxstyle="round,pad=0.004,rounding_size=0.012",
                                      transform=fig.transFigure, fc="none", ec=C["rule"], lw=0.9))
    fig.text(x + 0.012, y + CARD_H - 0.012, num, fontsize=14, fontweight="bold", color=C["trace"], va="top")
    fig.text(x + 0.045, y + CARD_H - 0.014, title, fontsize=9.5, fontweight="bold", color=C["ink"], va="top")
    fig.text(x + 0.045, y + CARD_H - 0.050, sub, fontsize=7.5, color=C["ink3"], va="top")
    pad = 0.035 if k in (1, 5) else 0.0
    plot_axes.append(fig.add_axes([x + 0.035 + pad, y + 0.075, CARD_W - 0.05 - pad, CARD_H - 0.16]))
    if k in INPUTS:
        fig.text(x + CARD_W / 2, y - 0.014, INPUTS[k], fontsize=7, color=C["ink2"], ha="center", va="top",
                 linespacing=1.15, style="italic")
# arrows along each row; the reader carries on from the end of the first row to the start of the second.
# None into step 6: the uncertainties belong to the line list of steps 1-4, not to the spectrum.
for row, i in ((0, 0), (0, 1), (1, 0)):
    y = YS[row] + CARD_H / 2
    fig.patches.append(FancyArrowPatch((XS[i] + CARD_W + 0.004, y), (XS[i + 1] - 0.004, y),
                                       transform=fig.transFigure, arrowstyle="-|>", mutation_scale=9,
                                       color=C["ink2"], lw=1.0))
fig.text(0.5, 0.985, "computed once per isotopologue, independent of temperature: steps 1–4 and 6", ha="center",
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
    ax.plot(R, pot[s](R), color=col, lw=1.4)
    for v in range(0, 80 if s == "B" else 108, 6):
        E = energy(s, v, 0)
        a, b = turning(pot[s], np.linspace(2.2, 12, 6000), E)
        ax.plot([a, min(b, 6.5)], [E, E], color=col, lw=0.5, alpha=0.6)
ax.text(2.95, 1500, "X", color=C["trace"], fontsize=9, fontweight="bold")
ax.text(3.45, 15800, "B", color=C["pick"], fontsize=9, fontweight="bold")
ax.set_ylim(-500, 21500)
clean(ax, "R (Å)")

# 6. uncertainties: 1σ of R(20) v'-0 and of R(20) 5-v'' (the center of each line)
ax = plot_axes[5]
vu = np.arange(0, 81)
u_b, _, _ = table.uncertainties(np.full(vu.size, "R"), np.full(vu.size, 20), vu, np.zeros(vu.size, int))
vl = np.arange(0, 109)
u_x, _, _ = table.uncertainties(np.full(vl.size, "R"), np.full(vl.size, 20), np.full(vl.size, 5), vl)
ax.semilogy(vu, u_b, color=C["pick"], lw=1.2)
ax.semilogy(vl, u_x, color=C["trace"], lw=1.2)
ax.set_xlim(0, 108)
ax.set_ylim(0.05, 500)
ax.text(0.03, 0.95, "R(20) 5–v″", color=C["trace"], fontsize=6.5, ha="left", va="top", transform=ax.transAxes)
ax.text(0.03, 0.83, "R(20) v′–0", color=C["pick"], fontsize=6.5, ha="left", va="top", transform=ax.transAxes)
ax.set_ylabel("1σ (MHz)", fontsize=7, labelpad=1)
ax.set_xlabel("v′ or v″", fontsize=7, labelpad=1)
ax.tick_params(labelsize=6.5, pad=1)

# 2. line positions: the 32-0 band near its head (Fortrat)
ax = plot_axes[1]
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

# 3. hyperfine: R(56) 32-0
ax = plot_axes[2]
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

# 4. line strengths: band strengths of the v'' = 0 progression at 300 K
ax = plot_axes[3]
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

# 5. spectrum: cross section and 10 cm cell transmission at 532 nm
ax = plot_axes[4]
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
