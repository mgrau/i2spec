"""Register Orsay part I (4 826 lines, 11 000-13 009 cm-1, 790 C) against our model and Hannover 2008.

Per 30 cm-1 window (about one plate): median distance from each atlas line to the nearest model line
among the strongest model lines at the atlas's own density, against a null that shifts the window's
atlas lines by a common offset in +-1 cm-1. Writes docs/figures/orsay_part1.png.
"""
import csv, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from i2spec.model import RovibronicModel
from i2spec.intensity import line_list, SHARED_GRID, X_RMIN

MHZ, T, W = 29979.2458, 1063.15, 30.0
obs = np.array([float(r["sigma_cm1"]) for r in csv.DictReader(open("data/atlas_lines/orsay1982_part1.csv"))])
g = dict(SHARED_GRID)
n_inner = max(int(round((g["rmin"] - X_RMIN) / g["step"])), 0)
gx = dict(g, rmin=g["rmin"] - n_inner * g["step"], nlev=48)
rng = np.random.default_rng(0)
offs = rng.uniform(-1, 1, 200)

def nearest(nu, x):
    i = np.clip(np.searchsorted(nu, x), 1, len(nu) - 1)
    return np.minimum(np.abs(x - nu[i - 1]), np.abs(nu[i] - x))

res = {}
for tag in ("i2spec2026f", "hannover2008"):
    m = RovibronicModel("127I2", tag, grids={"X": gx, "B": dict(g, nlev=70)}, solver="dvr")
    ll = line_list(m, T, 10995, 13015, S_min=3e-28)
    nu_all, S, vl = np.asarray(ll.nu), np.asarray(ll.S), np.asarray(ll.v_lower)
    out = []
    for lo in np.arange(11000, 13009, W):
        o = obs[(obs >= lo) & (obs < lo + W)]
        sel = (nu_all >= lo - 1.5) & (nu_all < lo + W + 1.5)
        k = np.argsort(S[sel])[::-1][: int(len(o) * (W + 3) / W * 3)]   # 3x the atlas density: the atlas misses blends
        nu = np.sort(nu_all[sel][k])
        v_mode = np.bincount(vl[sel][k]).argmax()
        real = np.median(nearest(nu, o)) * MHZ
        null = np.array([np.median(nearest(nu, o + d)) for d in offs]) * MHZ
        out.append((lo, len(o), v_mode, real, null.mean(), null.std(), (null.mean() - real) / null.std()))
    res[tag] = np.array(out)
    a = res[tag]
    print(f"{tag}: windows {len(a)}, median z {np.median(a[:,6]):.2f}, windows with z>3: {np.sum(a[:,6]>3)}, "
          f"median |r| {np.median(a[:,3]):.0f} MHz (null {np.median(a[:,4]):.0f})")
    for row in a[::7]:
        print(f"   {row[0]:.0f} n={row[1]:.0f} v''~{row[2]:.0f}  real {row[3]:5.0f}  null {row[4]:5.0f}+-{row[5]:3.0f}  z {row[6]:5.1f}")

fig, ax = plt.subplots(2, 1, figsize=(10, 6.5), sharex=True)
for tag, c, lab in (("i2spec2026f", "#4C72B0", "i2spec (this work)"), ("hannover2008", "#C44E52", "Hannover 2008 (published)")):
    a = res[tag]
    ax[0].plot(a[:, 0] + W / 2, a[:, 6], "o-", color=c, ms=3, lw=1, label=lab)
    ax[1].plot(a[:, 0] + W / 2, a[:, 3], "o-", color=c, ms=3, lw=1, label=lab)
ax[1].plot(res["i2spec2026f"][:, 0] + W / 2, res["i2spec2026f"][:, 4], color="gray", lw=1, ls="--", label="shifted-offset null")
ax[0].axhline(3, color="k", ls=":", lw=1); ax[0].set_ylabel("registration z (per plate)")
ax[1].set_ylabel("median |atlas - model| (MHz)"); ax[1].set_xlabel("wavenumber (cm$^{-1}$)")
ax2 = ax[0].twiny(); ax2.set_xlim(ax[0].get_xlim())
t = res["i2spec2026f"][::8]; ax2.set_xticks(t[:, 0] + W / 2); ax2.set_xticklabels([f"{v:.0f}" for v in t[:, 2]])
ax2.set_xlabel("dominant v$''$ of the strongest model lines")
for a_ in ax: a_.legend(fontsize=8, frameon=False)
fig.suptitle("Orsay atlas part I (4 826 lines, cell at 790 °C) against both X potentials", fontsize=11)
fig.tight_layout(); fig.savefig("docs/figures/orsay_part1.png", dpi=150)
np.save("/tmp/orsay_part1_res.npy", np.array([res["i2spec2026f"], res["hannover2008"]]))
