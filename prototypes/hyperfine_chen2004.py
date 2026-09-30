"""Chen 2004 B-state hyperfine parameters against the published formulae: the numbers in
docs/design/hyperfine-fit.md, section "Chen 2004".

Usage:  uv run --with matplotlib python prototypes/hyperfine_chen2004.py
        (after prototypes/hyperfine_stage2.py, whose per-line fits it reuses)
Output: prototypes/out/hyperfine_chen2004.json, docs/figures/hyperfine_chen2004.png   (about 1 minute)
"""

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from i2spec import hfs_params as hp  # noqa: E402
from i2spec.hfs_params import Corrections  # noqa: E402
from i2spec.model import RovibronicModel  # noqa: E402
from i2spec.observations import load_hyperfine_parameters  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
SURFACE, INK, INK2, GRID, AXIS = "#fcfcfb", "#1f1e1c", "#52514e", "#e1e0d9", "#c3c2b7"
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
PARAMS = ("eqQ", "C", "d", "delta")
STAGE2 = {"eqQ": "eqQ_b", "C": "C_b", "d": "d_b", "delta": "delta_b"}
FLOOR = {"eqQ": 0.010, "C": 0.02, "d": 0.5, "delta": 0.5}          # eqQ in MHz, the others in kHz
REGIONS = ("BKT02", "S06", "S06 frozen")
#: The B grid reaches v' = 79; the default stops at v' = 59. Level energies are converged at rmax = 12.
BIG_B = {"B": dict(rmin=2.35, rmax=12.0, h=0.01, order=10, nlev=80)}


def y(J):
    return J * (J + 1) / 1e4


def region(E_b):
    return "BKT02" if E_b <= hp.E_B_MAX else ("S06" if E_b <= hp.E_B_MAX_S06 else "S06 frozen")


def build_table(model):
    """measured - formula for every line: Chen's parameters, and the Stage 2 per-line fits of splittings."""
    rows = []
    for r in load_hyperfine_parameters(ROOT / "data/hyperfine_parameters/chen2004a").rows:
        L = r.line
        E_b, E_x = model._reference_term("B", L.v_upper), model._reference_term("X", L.v_lower)
        p = hp.line_states("127I2", L.v_upper, L.v_lower, E_b, E_x)[1](r.J_upper)
        formula = {"eqQ": p.eqQ, "C": p.C, "d": p.d, "delta": p.delta}
        measured = {"eqQ": r.eqQ, "C": r.C * 1e3, "d": r.d * 1e3, "delta": r.delta * 1e3}
        unc = {"eqQ": r.eqQ_unc, "C": r.C_unc * 1e3, "d": r.d_unc * 1e3, "delta": r.delta_unc * 1e3}
        rows.append(dict(source="chen2004a", line=f"{L.branch}({L.J_lower}) {L.v_upper}-{L.v_lower}", v=L.v_upper,
                         v_lower=L.v_lower, branch=L.branch, J_lower=L.J_lower, J=r.J_upper, region=region(E_b),
                         measured=measured, formula=formula,
                         delta={k: measured[k] - formula[k] for k in PARAMS}, unc=unc))
    stage2 = json.loads((ROOT / "prototypes/out/hyperfine_stage2.json").read_text())
    for r in stage2["per_line"]:
        iso, br, J_lower, vu, vl = r["key"]
        if iso != "127I2" or r["median_uncertainty"] > 25:
            continue
        rows.append(dict(source=r["dataset"], line=r["line"], v=vu, v_lower=vl, branch=br, J_lower=J_lower,
                         J=r["J_upper"], region=region(model._reference_term("B", vu)),
                         delta={k: r["x"][STAGE2[k]] for k in PARAMS}, unc={k: r["error"][STAGE2[k]] for k in PARAMS}))
    return rows


def wfit(X, t, w):
    return np.linalg.lstsq(X * w[:, None], t * w, rcond=None)[0]


def same_v(rows, i, k, degree):
    """Predict row i's correction from the other lines at its v', as a polynomial in y; None if too few."""
    r = rows[i]
    g = [q for j, q in enumerate(rows) if q["v"] == r["v"] and j != i]
    if len(g) < degree + 1:
        return None
    Y = np.array([y(q["J"]) for q in g])
    X = np.column_stack([Y ** n for n in range(degree + 1)])
    c = wfit(X, np.array([q["delta"][k] for q in g]), 1 / np.hypot([q["unc"][k] for q in g], FLOOR[k]))
    return float(sum(c[n] * y(r["J"]) ** n for n in range(degree + 1)))


def neighbouring_v(rows, i, k, degree):
    """Interpolate linearly in v' between the per-v' polynomials of the nearest measured v' either side."""
    r = rows[i]
    fits = {}
    for v in sorted({q["v"] for q in rows if q["v"] != r["v"]}):
        g = [q for q in rows if q["v"] == v]
        deg = min(degree, len(g) - 1)
        Y = np.array([y(q["J"]) for q in g])
        X = np.column_stack([Y ** n for n in range(deg + 1)])
        c = wfit(X, np.array([q["delta"][k] for q in g]), 1 / np.hypot([q["unc"][k] for q in g], FLOOR[k]))
        fits[v] = float(sum(c[n] * y(r["J"]) ** n for n in range(deg + 1)))
    below, above = [v for v in fits if v < r["v"]], [v for v in fits if v > r["v"]]
    if not below or not above:
        return None
    v0, v1 = max(below), min(above)
    return fits[v0] + (fits[v1] - fits[v0]) * (r["v"] - v0) / (v1 - v0)


def cross_validate(table):
    """Held-out rms by region and parameter, on the lines with at least three partners at their v'."""
    out = {}
    for reg in REGIONS:
        rows = [r for r in table if r["region"] == reg]
        for k in PARAMS:
            errs = {"published": [], "same v' line": [], "same v' quadratic": [], "neighbouring v'": []}
            for i, r in enumerate(rows):
                if sum(q["v"] == r["v"] for q in rows) < 4:
                    continue
                errs["published"].append(r["delta"][k])
                errs["same v' line"].append(r["delta"][k] - same_v(rows, i, k, 1))
                errs["same v' quadratic"].append(r["delta"][k] - same_v(rows, i, k, 2))
                nv = neighbouring_v(rows, i, k, 2 if k == "C" and reg != "BKT02" else 1)
                if nv is not None:
                    errs["neighbouring v'"].append(r["delta"][k] - nv)
            out[f"{reg}/{k}"] = {name: dict(rms=float(np.sqrt(np.mean(np.square(e)))) if e else None, n=len(e))
                                 for name, e in errs.items()}
    return out


def splitting_impact(table, model):
    """How far apart the model's hyperfine intervals are with the formulae and with Chen's parameters."""
    out = []
    for r in table:
        if r["source"] != "chen2004a":
            continue
        c = Corrections(eqQ_b=r["delta"]["eqQ"], C_b=r["delta"]["C"], d_b=r["delta"]["d"], delta_b=r["delta"]["delta"])
        _, pub = model.hyperfine_components(r["v"], r["v_lower"], r["J_lower"], r["branch"], table=None)
        _, meas = model.hyperfine_components(r["v"], r["v_lower"], r["J_lower"], r["branch"], corrections=c, table=None)
        a = np.array(sorted(q.offset for q in pub if q.label))
        b = np.array(sorted(q.offset for q in meas if q.label))
        d = (a - a[0]) - (b - b[0])
        out.append(dict(line=r["line"], v=r["v"], J=r["J"], region=r["region"],
                        rms_MHz=float(np.sqrt(np.mean(d ** 2))), max_MHz=float(np.abs(d).max())))
    return out


def main():
    model = RovibronicModel("127I2", grids=BIG_B)
    table = build_table(model)
    result = dict(table=table, cross_validation=cross_validate(table), splitting_impact=splitting_impact(table, model))
    (ROOT / "prototypes/out/hyperfine_chen2004.json").write_text(json.dumps(result, indent=1))
    figure(result)
    for key, v in result["cross_validation"].items():
        print(f"{key:<18}" + "  ".join(f"{n} {x['rms']:.3g} ({x['n']})" for n, x in v.items() if x["rms"] is not None))
    for reg in REGIONS:
        s = [x for x in result["splitting_impact"] if x["region"] == reg]
        print(f"{reg:<11} interval differences, median rms {np.median([x['rms_MHz'] for x in s]):.3f} MHz, "
              f"worst {max(x['max_MHz'] for x in s):.1f} MHz")


def style(ax):
    ax.set_facecolor(SURFACE)
    ax.grid(True, color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(AXIS)
    ax.tick_params(colors=INK2, labelsize=8.5)


def figure(result):
    table = result["table"]
    chen = [r for r in table if r["source"] == "chen2004a"]
    fig, (a, b) = plt.subplots(2, 1, figsize=(11, 10), facecolor=SURFACE, gridspec_kw=dict(height_ratios=[1, 1]))

    # (a) measured C_B against v', with the formula the model uses for the same lines
    style(a)
    a.axvspan(53.5, 71, color=GRID, alpha=0.55, zorder=0)
    a.scatter([r["v"] for r in chen], [r["measured"]["C"] for r in chen], s=34, marker="o", color=BLUE, zorder=3,
              label="measured, Chen 2004 (4–7 J′ per v′)")
    a.scatter([r["v"] for r in chen], [r["formula"]["C"] for r in chen], s=34, marker="s", facecolor="none",
              edgecolor=ORANGE, linewidths=1.3, zorder=3, label="formula the model uses, same lines")
    a.text(54.0, 2150, "above v′ = 53 the model holds the formulae at the v′ = 53 energy", fontsize=8.5, color=INK2)
    a.set_xlabel("v′", color=INK)
    a.set_ylabel("C_B (kHz)", color=INK)
    a.set_title("(a) C_B keeps rising toward dissociation, to 2.2 MHz; the frozen formula stays at 500–650 kHz", loc="left",
                color=INK, fontsize=11)
    a.legend(frameon=False, fontsize=8.5, loc="upper left", labelcolor=INK)

    # (b) measured - formula, C_B, for v' <= 53: smooth in J within a v', not across v'
    style(b)
    groups = sorted({r["v"] for r in table if r["region"] != "S06 frozen" and sum(q["v"] == r["v"] for q in table) >= 4})
    cmap = plt.get_cmap("viridis")
    for n, v in enumerate(groups):
        g = sorted([r for r in table if r["v"] == v], key=lambda r: r["J"])
        colour = cmap(n / max(len(groups) - 1, 1))
        b.errorbar([r["J"] for r in g], [r["delta"]["C"] for r in g], yerr=[r["unc"]["C"] for r in g], fmt="-o", ms=4, lw=1,
                   color=colour, ecolor=colour, elinewidth=0.8, zorder=3)
        b.annotate(f"v′={v}", (g[-1]["J"], g[-1]["delta"]["C"]), textcoords="offset points", xytext=(5, -3), fontsize=7.5,
                   color=INK2)
    b.axhline(0, color=AXIS, lw=0.8)
    b.set_xlabel("J′", color=INK)
    b.set_ylabel("C_B, measured − formula (kHz)", color=INK)
    b.set_title("(b) Within a v′ the error is smooth in J′ and predictable; between v′ it jumps", loc="left", color=INK, fontsize=11)
    fig.tight_layout(h_pad=2.0)
    fig.savefig(ROOT / "docs/figures/hyperfine_chen2004.png", dpi=130, facecolor=SURFACE)


if __name__ == "__main__":
    main()
