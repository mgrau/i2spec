"""The Phase C trial: an MLR X potential against the published X-representation one.

Usage:  uv run --with matplotlib python prototypes/mlr_x.py
Output: prototypes/out/mlr_x.json, docs/figures/mlr_x.png   (about 5 minutes)

Loads the fitted parameters from data/potentials/mlr_x_2026a.json; it does not refit. The fit itself is
in docs/design/mlr-x.md, and the script that produced it is recorded there.
"""

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from i2spec.bspline import BSplineSolver  # noqa: E402
from i2spec.constants import MHZ_PER_CM, reduced_mass  # noqa: E402
from i2spec.observations import Predictor, load_all, residuals  # noqa: E402
from i2spec.potentials import MLRPotential, load_potentials  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
SURFACE, INK, INK2, GRID, AXIS = "#fcfcfb", "#1f1e1c", "#52514e", "#e1e0d9", "#c3c2b7"
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
SOLVER = dict(rmin=2.10, rmax=6.0, h=0.01, order=10, nlev=70)
MEASURED = (48, 53, 54)


def load_mlr(path=ROOT / "data/potentials/mlr_x_2026a.json", bo=None):
    """The fitted MLR X potential; its Born-Oppenheimer corrections come from ``bo``."""
    d = json.loads(Path(path).read_text())
    return MLRPotential(De=d["De"], Re=d["Re"], C={int(k): v for k, v in d["C"].items()},
                        beta=tuple(d["beta"]), p=d["p"], q=d["q"], Rref=d["Rref"],
                        bo=bo if bo is not None else load_potentials()["X"])


def level_shifts(published, mlr, v_max=64):
    a = BSplineSolver(published, reduced_mass("127I2"), **SOLVER)
    b = BSplineSolver(mlr, reduced_mass("127I2"), **SOLVER, joins=())
    e0a, e0b = a.energy(0, 0), b.energy(0, 0)
    out = []
    for v in range(v_max + 1):
        out.append(dict(v=v, published=a.energy(v, 0) - e0a, mlr=b.energy(v, 0) - e0b,
                        shift=(b.energy(v, 0) - e0b) - (a.energy(v, 0) - e0a)))
    return out


def dataset_residuals(mlr):
    pub, new = Predictor(), Predictor(potentials={"X": mlr})
    rows = []
    for ds in load_all():
        scale = MHZ_PER_CM if ds.meta["unit"] == "cm-1" else 1.0
        a, b = residuals(pub, ds) * scale, residuals(new, ds) * scale
        vs = [o.line.v_lower for o in ds.observations]
        rows.append(dict(id=ds.id, n=len(a), v_min=min(vs), v_max=max(vs),
                         published=float(np.sqrt(np.mean(a ** 2))), mlr=float(np.sqrt(np.mean(b ** 2)))))
    return rows


def main(figure_only=False):
    published = load_potentials()["X"]
    mlr = load_mlr(bo=published)
    cache = ROOT / "prototypes/out/mlr_x.json"
    if figure_only and cache.exists():
        result = json.loads(cache.read_text())
    else:
        result = dict(shifts=level_shifts(published, mlr), datasets=dataset_residuals(mlr))
        R = np.linspace(2.4, 6.0, 400)
        result["curve"] = dict(R=R.tolist(), difference=(mlr(R) - published(R)).tolist())
        cache.parent.mkdir(exist_ok=True)
        cache.write_text(json.dumps(result, indent=1))
    figure(result, published, mlr)
    s = {r["v"]: r["shift"] for r in result["shifts"]}
    print("level shifts (cm-1):", {v: round(s[v], 3) for v in (10, 17, 30, 42, 48, 53, 54, 60)})
    for r in result["datasets"]:
        if r["mlr"] > 10 or r["published"] > 10:
            print(f"  {r['id']:<16} v'' {r['v_min']}-{r['v_max']}: {r['published']:12.1f} -> {r['mlr']:10.1f} MHz")


def style(ax):
    ax.set_facecolor(SURFACE)
    ax.grid(True, color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(AXIS)
    ax.tick_params(colors=INK2, labelsize=8.5)


def figure(result, published, mlr):
    fig, (a, b, c) = plt.subplots(3, 1, figsize=(11, 12), facecolor=SURFACE)
    R = np.array(result["curve"]["R"])

    style(a)
    a.plot(R, np.array(result["curve"]["difference"]), color=BLUE, lw=1.8, label="MLR − published")
    outer = R[R >= published.RO]                       # the published tail only exists beyond R_O = 3.3 Å
    tail = published.De - sum(cn / outer ** n for n, cn in published.C.items())
    a.plot(outer, tail - published(outer), color=ORANGE, lw=1.4, ls=(0, (5, 3)),
           label="the published outer exponential, −A_O exp(−B_O (R − R_O))")
    for v, Rt in ((17, 3.02), (48, 3.49), (54, 3.61)):
        a.axvline(Rt, color=AXIS, lw=0.9, ls=(0, (2, 3)), zorder=1)
        a.annotate(f"v″={v} turns here", (Rt, -105), textcoords="offset points", xytext=(4, 0), fontsize=7.5,
                   color=INK2, rotation=90, va="bottom")
    a.axhline(0, color=AXIS, lw=0.8)
    a.set_xlim(2.4, 6.0)
    a.set_ylim(-120, 360)
    a.set_xlabel("R (Å)", color=INK)
    a.set_ylabel("ΔV (cm⁻¹)", color=INK)
    a.set_title("(a) The two potentials differ where the published one splices on its outer branch",
                loc="left", color=INK, fontsize=11)
    a.legend(frameon=False, fontsize=8.5, loc="upper left", labelcolor=INK)

    style(b)
    shifts = result["shifts"]
    v = np.array([r["v"] for r in shifts])
    s = np.array([r["shift"] for r in shifts])
    b.axvspan(17.5, 47.5, color=GRID, alpha=0.6, zorder=0)
    b.text(32.5, s.min() * 0.92, "v″ = 18–47: never measured", ha="center", fontsize=9, color=INK2)
    b.plot(v, s, color=BLUE, lw=2, zorder=3, label="MLR − published")
    b.scatter([q for q in MEASURED], [s[q] for q in MEASURED], s=70, marker="o", color=ORANGE, zorder=4,
              label="required by matyugin2012 and nesterenko2019")
    b.axhline(0, color=AXIS, lw=0.8)
    b.set_xlabel("v″", color=INK)
    b.set_ylabel("level shift (cm⁻¹)", color=INK)
    b.set_title("(b) Fitting the three measured levels bends v″ = 18–47 by up to −21 cm⁻¹", loc="left",
                color=INK, fontsize=11)
    b.legend(frameon=False, fontsize=8.5, loc="upper left", labelcolor=INK)

    style(c)
    rows = sorted(result["datasets"], key=lambda r: -r["published"])
    x = np.arange(len(rows))
    w = 0.38
    c.bar(x - w / 2 - 0.01, [r["published"] for r in rows], w, color=BLUE, label="published X", zorder=3)
    c.bar(x + w / 2 + 0.01, [r["mlr"] for r in rows], w, color=ORANGE, label="MLR X", zorder=3)
    c.set_yscale("log")
    c.set_xticks(x, [f"{r['id']}\nv″ {r['v_min']}–{r['v_max']}" for r in rows], rotation=60, ha="right", fontsize=7.5)
    c.set_ylabel("rms residual (MHz)", color=INK)
    c.set_title("(c) Against every measurement: the high-v″ sets fall by 10⁴, the rest hold", loc="left",
                color=INK, fontsize=11)
    c.legend(frameon=False, fontsize=8.5, loc="upper right", labelcolor=INK)

    fig.tight_layout(h_pad=2.0)
    fig.savefig(ROOT / "docs/figures/mlr_x.png", dpi=130, facecolor=SURFACE)


if __name__ == "__main__":
    import sys

    main(figure_only="--figure-only" in sys.argv)
