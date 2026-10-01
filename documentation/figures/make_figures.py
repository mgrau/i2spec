"""Explanatory figures for the documentation site (documentation/), drawn from the model itself.

    uv run --group research python documentation/figures/make_figures.py

They are illustrations, not results: the levels, wavefunctions and spectra are the model's, but the
choice of what to draw is for explanation. The SVGs are written with placeholder colours that are
then replaced by CSS variables (documentation/stylesheets/extra.css), and the pages inline them, so
they follow the light and dark themes like the rest of the page.
"""
import re
from pathlib import Path

import matplotlib
matplotlib.use("svg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import PathPatch, Patch
from matplotlib.path import Path as MPath

from i2spec.constants import reduced_mass
from i2spec.intensity import mu_tellinghuisen2011
from i2spec.lookup import Catalog, Line, uncertainty
from i2spec.model import RovibronicModel
from i2spec.potentials import load_extended, load_potentials
from i2spec.solver import RadialSolver
from i2spec.spectrum import cross_section, number_density, transmission

OUT = Path(__file__).resolve().parent
# placeholder colours -> the site's tokens (web/site.css)
C = dict(ink="#010101", ink2="#020202", ink3="#030303", rule="#040404", trace="#050505", pick="#060606",
         mark="#070707", surface="#080808", u1="#0a0a01", u2="#0a0a02", u3="#0a0a03", u4="#0a0a04", u5="#0a0a05")
TOKENS = {**{C[f"u{k}"]: f"var(--unc-{k})" for k in range(1, 6)},C["ink"]: "var(--ink)", C["ink2"]: "var(--ink-2)", C["ink3"]: "var(--ink-3)", C["rule"]: "var(--rule)",
          C["trace"]: "var(--trace)", C["pick"]: "var(--pick)", C["mark"]: "var(--mark)", C["surface"]: "var(--surface)"}

plt.rcParams.update({
    "svg.fonttype": "none", "font.family": "sans-serif", "font.size": 10.5,
    "axes.edgecolor": C["rule"], "axes.labelcolor": C["ink2"], "xtick.color": C["ink3"], "ytick.color": C["ink3"],
    "text.color": C["ink"], "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": 0.8,
    "xtick.major.size": 3, "ytick.major.size": 3, "legend.frameon": False, "figure.facecolor": "none",
    "axes.facecolor": "none", "savefig.facecolor": "none", "lines.solid_capstyle": "round",
})


def save(fig, name):
    path = OUT / f"{name}.svg"
    fig.savefig(path, bbox_inches="tight", pad_inches=0.05)
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


def turning_points(V, R, E):
    inside = np.flatnonzero(V(R) <= E)
    return R[inside[0]], R[inside[-1]]


model = RovibronicModel("127I2")
ext, from_v = load_extended()
pub = load_potentials()
e00 = model.energy("X", 0, 0)
mu = reduced_mass("127I2")

# --- 1. the two potentials, their levels, and the 532 nm transition -------------------------------
R = np.linspace(2.3, 7.0, 1500)
fig, ax = plt.subplots(figsize=(7.2, 4.6))
for s, colour in (("X", C["trace"]), ("B", C["pick"])):
    V = ext[s]
    ax.plot(R, V(R), color=colour, lw=2)
    levels = model.levels(s, 0)
    top = {"X": 47, "B": 58}[s]
    fitted = {"X": 17, "B": 43}[s]
    for v in range(0, top + 1):
        if v > fitted and v % 3:            # beyond the fitted range, every third level is enough
            continue
        if s == "B" and v <= fitted and v % 2:
            continue
        E = levels[v]
        a, b = turning_points(V, R, E)
        ax.plot([a, b], [E, E], color=colour, lw=0.7, alpha=0.9 if v <= fitted else 0.35)
    E = levels[fitted]
    a, b = turning_points(V, R, E)
    ax.annotate(f"{'v″' if s == 'X' else 'v′'} = {fitted}: highest level in the 2008 fit", xy=(b, E),
                xytext=(b + 0.45, E - (900 if s == "X" else 1300)),
                fontsize=8.5, color=C["ink3"], va="center", arrowprops=dict(arrowstyle="-", color=C["ink3"], lw=0.6))
# dissociation limits
ax.axhline(ext["X"].De if hasattr(ext["X"], "De") else 12547.3, color=C["ink3"], lw=0.7, ls=(0, (4, 3)))
ax.axhline(20150.3, color=C["ink3"], lw=0.7, ls=(0, (4, 3)))
ax.text(6.95, 12547.3 + 250, "I ²P3/2 + I ²P3/2", ha="right", fontsize=8.5, color=C["ink3"])
ax.text(6.95, 20150.3 + 250, "I ²P3/2 + I ²P1/2", ha="right", fontsize=8.5, color=C["ink3"])
ax.text(3.25, 700, "X ¹Σg⁺", color=C["trace"], fontsize=12, fontweight="bold")
ax.text(3.55, 15600, "B ³Π(0u⁺)", color=C["pick"], fontsize=12, fontweight="bold")
# the 532 nm line, R(56) 32-0, drawn vertically at the equilibrium distance of X (Franck-Condon)
Ex, Eb = model.energy("X", 0, 56), model.energy("B", 32, 57)
ax.annotate("", xy=(2.70, Eb), xytext=(2.70, Ex), arrowprops=dict(arrowstyle="-|>", color=C["ink"], lw=1.4))
ax.text(2.76, 13600, "532 nm\nR(56) 32–0", ha="left", va="center", fontsize=9, color=C["ink"])
ax.set_xlim(2.3, 7.0)
ax.set_ylim(-800, 22500)
ax.set_xlabel("internuclear distance R (Å)")
ax.set_ylabel("energy above the X minimum (cm⁻¹)")
save(fig, "potentials")

# --- 2. published curve vs MLR: where each is trusted ---------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(7.2, 2.7), sharey=False)
for ax, s, colour, span in ((axes[0], "X", C["trace"], (2.427, 3.079)), (axes[1], "B", C["pick"], (2.65, 4.589))):
    R = np.linspace(span[0] - 0.1, 8.0, 800)
    d = pub[s](R) - ext[s](R)
    ax.axvspan(*span, color=colour, alpha=0.12, lw=0)
    ax.axhline(0, color=C["rule"], lw=0.8)
    ax.plot(R, d, color=colour, lw=1.8)
    ax.text(span[1] + 0.1, 0.12, "← range of the\n    fitted levels", transform=ax.get_xaxis_transform(),
            ha="left", va="bottom", fontsize=8, color=C["ink3"])
    ax.set_title(f"{s} state: V(2008) − V(MLR)", fontsize=10, color=C["ink2"], loc="left")
    ax.set_xlabel("R (Å)")
    lim = np.abs(d[R > span[0]]).max() * 1.15
    ax.set_ylim(-lim, lim)
axes[0].set_ylabel("ΔV (cm⁻¹)")
fig.tight_layout()
save(fig, "published_vs_mlr")

# --- 3. wavefunctions and the transition moment: Franck-Condon ------------------------------------
sx = RadialSolver(ext["X"], mu, rmin=2.3, rmax=4.0, step=0.004, nlev=6)
sb = RadialSolver(ext["B"], mu, rmin=2.4, rmax=7.0, step=0.004, nlev=40)
ex, wx = sx.states(0)
eb, wb = sb.states(0)
fig, (ax, ax2) = plt.subplots(2, 1, figsize=(7.2, 4.8), sharex=True, gridspec_kw=dict(height_ratios=[3, 1.2], hspace=0.08))
R = np.linspace(2.3, 7.0, 1200)
ax.plot(R, ext["X"](R), color=C["trace"], lw=1.4)
ax.plot(R, ext["B"](R), color=C["pick"], lw=1.4)
for solver, e, w, v, colour, height, V in ((sx, ex, wx, 0, C["trace"], 2600, ext["X"]),
                                          (sb, eb, wb, 5, C["pick"], 1100, ext["B"]),
                                          (sb, eb, wb, 32, C["pick"], 1100, ext["B"])):
    psi = w[:, v] / np.abs(w[:, v]).max() * height
    psi *= np.sign(psi[np.argmax(np.abs(psi))])
    keep = np.abs(psi) > 1e-3 * height
    a, b = turning_points(V, solver.R, e[v])
    ax.plot([a, b], [e[v], e[v]], color=colour, lw=0.6, alpha=0.6)
    ax.fill_between(solver.R[keep], e[v], e[v] + psi[keep], color=colour, alpha=0.18, lw=0)
    ax.plot(solver.R[keep], e[v] + psi[keep], color=colour, lw=1.1)
    ax.text(b + 0.06, e[v] + (500 if colour == C["trace"] else -700), f"{'v″' if colour == C['trace'] else 'v′'} = {v}", fontsize=9, color=colour, va="center")
ax.axvspan(2.55, 2.80, color=C["ink3"], alpha=0.08, lw=0)
ax.text(2.675, 23000, "Franck–Condon\nregion", ha="center", va="top", fontsize=8.5, color=C["ink3"])
ax.set_ylim(-500, 23200)
ax.set_ylabel("energy (cm⁻¹)")
ax2.plot(R, mu_tellinghuisen2011(R), color=C["ink"], lw=1.6)
ax2.axvspan(2.55, 2.80, color=C["ink3"], alpha=0.08, lw=0)
ax2.set_ylabel("|μₑ| (D)")
ax2.set_xlabel("internuclear distance R (Å)")
ax2.set_xlim(2.3, 5.5)
save(fig, "wavefunctions")

# --- 4. rotational structure: the Fortrat diagram of 32-0 ------------------------------------------
J = np.arange(1, 41)
nuR = np.array([model.transition(32, 0, j, "R") for j in J])
nuP = np.array([model.transition(32, 0, j, "P") for j in J])
fig, ax = plt.subplots(figsize=(7.2, 3.3))
ax.plot(nuR, J, "o", ms=2.4, color=C["trace"], label="R branch  (J′ = J″ + 1)")
ax.plot(nuP, J, "o", ms=2.4, color=C["pick"], label="P branch  (J′ = J″ − 1)")
k = np.argmax(nuR)
ax.annotate(f"band head, J″ = {J[k]}", xy=(nuR[k], J[k]), xytext=(nuR[k] - 6, J[k] + 3), fontsize=9, ha="right",
            color=C["ink2"], arrowprops=dict(arrowstyle="-", color=C["ink3"], lw=0.6))
nu0 = model.energy("B", 32, 0) - model.energy("X", 0, 0)
ax.axvline(nu0, color=C["ink3"], lw=0.7, ls=(0, (4, 3)))
ax.text(nu0 - 0.1, 38, "band origin", ha="right", fontsize=8.5, color=C["ink3"])
ax.set_xlim(right=nu0 + 1.5)
ax.set_xlabel("wavenumber (cm⁻¹)")
ax.set_ylabel("J″")
ax.legend(loc="lower left", fontsize=9)
ax.set_title("32–0 band, J″ ≤ 40 (the band extends to J″ ≈ 170)", fontsize=10, color=C["ink2"], loc="left")
save(fig, "fortrat")

# --- 5. hyperfine structure: 15 and 21 components ----------------------------------------------------
fig, axes = plt.subplots(2, 1, figsize=(7.2, 3.9), sharex=True)
for ax, (branch, j, name) in zip(axes, (("R", 56, "R(56) 32–0, J″ even: 15 components"),
                                       ("P", 53, "P(53) 32–0, J″ odd: 21 components"))):
    nu0, comps = model.hyperfine_components(32, 0, j, branch)
    main = [c for c in comps if c.label]
    off = np.array([c.offset for c in main])
    st = np.array([c.strength for c in main])
    x = np.linspace(-1100, 1100, 1500)
    sig = 434 / 2.3548
    prof = sum(s * np.exp(-0.5 * ((x - o) / sig) ** 2) for o, s in zip(off, st))
    ax.fill_between(x, 0, prof / prof.max(), color=C["ink3"], alpha=0.12, lw=0)
    ax.plot(x, prof / prof.max(), color=C["ink3"], lw=1)
    ax.vlines(off, 0, st / st.max() * 0.85, color=C["trace"], lw=1.6)
    for c in main:
        if c.label in ("a1", f"a{len(main)}") or (len(main) == 15 and c.label == "a10"):
            ax.text(c.offset, c.strength / st.max() * 0.85 + 0.04, c.label, ha="center", fontsize=8, color=C["ink2"])
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.set_title(name, fontsize=10, color=C["ink2"], loc="left")
    ax.set_ylim(0, 1.15)
axes[1].set_xlabel("offset from the hyperfine-free position (MHz)")
axes[0].text(1080, 0.95, "Doppler profile at 20 °C", ha="right", fontsize=8.5, color=C["ink3"])
fig.tight_layout()
save(fig, "hyperfine")

# --- 6. a spectrum: cross section and cell transmission at 532 nm ---------------------------------
catalog = Catalog(temperature=293.15)
master = catalog.master("127I2")
T = 293.15
lines = master.at(T, 18787.7, 18789.0, S_min=1e-26)
grid = np.linspace(18787.75, 18789.0, 3000)
sigma = cross_section(grid, lines, T)
N = number_density(293.15, T) * 10.0
fig, (ax, ax2) = plt.subplots(2, 1, figsize=(7.2, 3.8), sharex=True, gridspec_kw=dict(hspace=0.1))
lam = 1e7 / grid
sigma16 = sigma * 1e16
ax.fill_between(lam, 0, sigma16, color=C["trace"], alpha=0.15, lw=0)
ax.plot(lam, sigma16, color=C["trace"], lw=1.4)
placed = []
for k in np.argsort(lines.S)[::-1]:
    x = 1e7 / lines.nu[k]
    if any(abs(x - p) < 0.0015 for p in placed) or len(placed) == 4:
        continue
    near = any(abs(x - p) < 0.005 for p in placed)
    placed.append(x)
    lab = f"{'R' if lines.branch[k] > 0 else 'P'}({lines.J_lower[k]}) {lines.v_upper[k]}–{lines.v_lower[k]}"
    y = np.interp(x, lam[::-1], sigma16[::-1])
    ax.text(x, y + 0.04 * sigma16.max(), lab, ha="left" if near else "center", fontsize=8, color=C["ink2"])
ax.set_ylabel("σ (10⁻¹⁶ cm²)")
ax.set_ylim(0, sigma16.max() * 1.25)
ax2.xaxis.get_major_formatter().set_useOffset(False)
ax2.plot(lam, transmission(sigma, N), color=C["pick"], lw=1.4)
ax2.set_ylabel("transmission")
ax2.set_xlabel("vacuum wavelength (nm)")
ax2.text(0.99, 0.08, "10 cm cell, cold finger 20 °C", transform=ax2.transAxes, ha="right", fontsize=8.5, color=C["ink3"])
ax2.invert_xaxis()
save(fig, "spectrum")

# --- 7. how well each band is known ------------------------------------------------------------------
VU, VL = 60, 54
classes = [(1, "≤ 1 MHz", C["u1"]), (10, "1–10 MHz", C["u2"]), (100, "10–100 MHz", C["u3"]),
           (1000, "0.1–1 GHz", C["u4"]), (np.inf, "> 1 GHz", C["u5"])]
cells = {k: [] for k in range(len(classes))}
for vu in range(VU):
    for vl in range(VL + 1):
        nu = model.energy("B", vu, 51) - model.energy("X", vl, 50)
        if nu <= 0:
            continue
        u, _ = uncertainty(Line("127I2", "R", 50, vu, vl, nu, 0.0, 0.0, 300.0))
        k = next(i for i, (lim, _, _) in enumerate(classes) if u <= lim)
        cells[k].append((vl, vu))
fig, ax = plt.subplots(figsize=(7.2, 4.4))
for k, pts in cells.items():
    if not pts:
        continue
    verts, codes = [], []
    for x, y in pts:
        verts += [(x - 0.5, y - 0.5), (x + 0.5, y - 0.5), (x + 0.5, y + 0.5), (x - 0.5, y + 0.5), (0, 0)]
        codes += [MPath.MOVETO, MPath.LINETO, MPath.LINETO, MPath.LINETO, MPath.CLOSEPOLY]
    ax.add_patch(PathPatch(MPath(verts, codes), facecolor=classes[k][2], lw=0))
ax.axvline(17.5, color=C["ink3"], lw=0.7, ls=(0, (4, 3)))
ax.axhline(43.5, color=C["ink3"], lw=0.7, ls=(0, (4, 3)))
ax.text(17.8, VU - 1.5, "v″ > 17: beyond the 2008 fit", fontsize=8.5, color=C["ink2"], va="top")
ax.text(VL, 44.2, "v′ > 43: beyond the 2008 fit", fontsize=8.5, color=C["ink2"], ha="right")
ax.set_xlim(-0.5, VL + 0.5)
ax.set_ylim(-0.5, VU - 0.5)
ax.set_xlabel("v″ (X)")
ax.set_ylabel("v′ (B)")
ax.legend(handles=[Patch(facecolor=col, label=l) for k, (_, l, col) in enumerate(classes) if cells[k]], loc="upper left",
          bbox_to_anchor=(1.01, 1.0), fontsize=9, title="1σ at J″ = 50", title_fontsize=9)
save(fig, "uncertainty_map")


# --- 8. hyperfine levels of R(56) 32-0 and the 15 main transitions -----------------------------------
from i2spec import hfs_params
from i2spec.constants import MHZ_PER_CM, NUCLEAR_SPIN
from i2spec.hfs_table import default_table
from i2spec.hyperfine import level_structure, line_components

vu, vl, Jl = 32, 0, 56
Ju = Jl + 1
xp, bp = hfs_params.line_states("127I2", vu, vl, model._reference_term("B", vu), model._reference_term("X", vl),
                                None, default_table())
ratios = dict(zip(("eqQ_ratio", "C_ratio"), hfs_params.nucleus_ratios("127I2")))
spins = dict(i1=NUCLEAR_SPIN[127], i2=NUCLEAR_SPIN[127], dJ=2, **ratios)
lower = level_structure(Jl, xp, lambda J: model.energy("X", vl, J) * MHZ_PER_CM, symmetry="g", **spins)
upper = level_structure(Ju, bp, lambda J: model.energy("B", vu, J) * MHZ_PER_CM, symmetry="u", **spins)
comps = [c for c in line_components(upper, lower, Ju, Jl) if c.label]

fig, ax = plt.subplots(figsize=(7.2, 4.6))
eu = np.array([u.energy for u in upper]); el = np.array([l.energy for l in lower])
eu0, el0 = eu.mean(), el.mean()
SU = 3.0                                        # the upper group is drawn 3x enlarged, with its own scale bar
top_lo = (el.max() - el0)
base_up = top_lo + 700 + SU * (eu0 - eu.min())  # vertical gap between the groups: not to scale
y_up = lambda e: base_up + SU * (e - eu0)       # noqa: E731
y_lo = lambda e: (e - el0)                      # noqa: E731
X0, X1 = 0.24, 0.99
for u in upper:
    ax.plot([X0, X1], [y_up(u.energy)] * 2, color=C["pick"], lw=1.0)
for l in lower:
    ax.plot([X0, X1], [y_lo(l.energy)] * 2, color=C["trace"], lw=1.0)
ordered = sorted(comps, key=lambda c: c.offset)
xs = np.linspace(X0 + 0.03, X1 - 0.03, len(ordered))
for x, c in zip(xs, ordered):
    # the pair of levels this component joins: same F and I, energy difference equal to its offset
    # the pair of levels joined: right F on each side and an energy difference equal to the offset (I is
    # only approximately good at high J, so it is not used to match)
    u, l = next((u, l) for u in upper for l in lower
                if u.F == c.F_upper and l.F == c.F_lower and abs((u.energy - l.energy) - c.offset) < 1e-6)
    ax.annotate("", xy=(x, y_up(u.energy)), xytext=(x, y_lo(l.energy)),
                arrowprops=dict(arrowstyle="-|>", color=C["ink2"], lw=0.8, shrinkA=0, shrinkB=0, mutation_scale=7))
    ax.text(x, y_lo(el.min()) - 70, c.label, ha="center", va="top", fontsize=8, color=C["ink2"])
ax.text(0.0, y_up(eu.max()), f"B³Π(0u⁺)\nv′ = {vu}, J′ = {Ju}\n{len(upper)} sublevels", va="top", fontsize=9.5, color=C["pick"])
ax.text(0.0, y_lo(el.max()), f"X¹Σg⁺\nv″ = {vl}, J″ = {Jl}\n{len(lower)} sublevels", va="top", fontsize=9.5, color=C["trace"])
for y0, length, text in ((y_up(eu.min()), SU * 100, "100 MHz"), (y_lo(el.min()), 500, "500 MHz")):
    ax.plot([0.195, 0.195], [y0, y0 + length], color=C["ink"], lw=1.5)
    ax.text(0.185, y0 + length / 2, text, ha="right", va="center", fontsize=8.5, color=C["ink"])
ax.set_xlim(0, 1)
ax.set_ylim(y_lo(el.min()) - 260, y_up(eu.max()) + 80)
ax.axis("off")
save(fig, "hyperfine_levels")
print("hyperfine levels: upper spread %.0f MHz, lower spread %.0f MHz, %d main components" % (np.ptp(eu), np.ptp(el), len(comps)))
