"""Refit the extended X potential to Martin et al. 1986's level constants, so X reaches its dissociation limit.

Usage:  uv run python prototypes/mlr_x_martin.py [--nbeta=N] [--write=mlr_x_2026e]
Output: prototypes/out/mlr_x_martin.json (and src/i2spec/data/<name>.json with --write)

mlr_x_2026d (prototypes/mlr_x_atlas.py) is measured to v'' = 25 (Orsay atlas part I) and at v'' = 48, 53, 54
(emission lines), and interpolated between; above v'' ~ 60 it is off by up to 70 cm-1 against Martin et al.,
J. Mol. Spectrosc. 116, 71 (1986), Table I (data/x_levels/martin1986/): G, B, D, H, L, M of 93 levels,
v'' = 8-108, from 14 820 B->X fluorescence lines. This refit reproduces:

  * v'' <= 17: mlr_x_2026d's own levels at 3 MHz (the hybrid takes the published curve there; these only keep
    the MLR continuous with it);
  * v'' = 18-25 and 48, 53, 54: the corrected levels of level_corrections_2026j at their own uncertainty, as
    in mlr_x_atlas.py;
  * Martin's levels at every other v'' from 26, at J inside each level's coverage (0-120 for v'' <= 85, 0-60 for
    86-96, 0-20 above, where the higher J approach the asymptote), with the paper's own precision of lines
    recomputed from the constants (Sect. V: 0.005 cm-1 to J = 80, 0.01 to 100, 0.02 to 120), and one free
    offset for Martin's origin (G(9) fixed to Luc's value, ~1 mK from the comb-referenced levels);
  * De, C6, C8, C10: the priors of Bacis, Cerny & Martin 1986, as before.

Levels come from a 40 A box graded beyond 5 A (BSplineSolver(mesh=...)), which holds every level to v'' ~ 113.
The Jacobian is Hellmann-Feynman: dE/dtheta = <psi| dV/dtheta |psi>, dV by central differences.
"""
import csv
import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

from i2spec.bspline import BSplineSolver
from i2spec.constants import MHZ_PER_CM, reduced_mass
from i2spec.level_corrections import load_level_corrections
from i2spec.potentials import MLRPotential, load_potentials

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "prototypes/out"
X0 = load_potentials()["X"]
MU = reduced_mass("127I2")
GRID = dict(rmin=2.10, rmax=40.0, h=0.01, order=10, nlev=118, mesh=((5.0, 0.02), (8.0, 0.05), (15.0, 0.1)))
START = json.loads((ROOT / "src/i2spec/data/mlr_x_2026d.json").read_text())
CORR = load_level_corrections("level_corrections_2026j")
PRIORS = {"De": (X0.De, 0.128), "C6": (X0.C[6], 0.12e6), "C8": (X0.C[8], 1.20e7), "C10": (X0.C[10], 0.5e8)}
#: Martin's origin is Luc's G(9), about 1 mK from the comb-referenced levels: the offset carries 0(2) mK, so it
#: cannot take up a shape error of the potential (free, it ran to -1.2 GHz)
OFFSET_PRIOR = 0.002
STAGES = (60, 75, 89, 108)      # Martin's levels are added in stages: the start is 70 cm-1 off near v'' = 95
ATLAS, HIGH = range(18, 26), (48, 53, 54)
ANCHOR_SIGMA = 3.0      # MHz, v'' <= 17
SKIP_MARTIN = {24}      # its G sits 0.009 cm-1 from the Orsay-atlas level while its neighbours agree to 0.001-0.002


def martin_levels():
    rows = list(csv.DictReader(open(ROOT / "data/x_levels/martin1986/data.csv")))
    out = {}
    for r in rows:
        f = lambda k, s=1.0: float(r[k]) * s if r.get(k) else 0.0      # noqa: E731
        out[int(r["v"])] = (f("G"), f("B_x10", 0.1), f("D_x1e8", 1e-8), -f("negH_x1e14", 1e-14), -f("negL_x1e20", 1e-20),
                            -f("negM_x1e23", 1e-23))
    return out


def martin_term(c, J):
    x = J * (J + 1.0)
    G, B, D, H, L, M = c
    return G + B * x - D * x**2 + H * x**3 + L * x**4 + M * x**5


def martin_sigma(J):
    return 0.005 if J <= 80 else 0.01 if J <= 100 else 0.02


def potential(x, n):
    De, Re, C6, C8, C10 = (float(v) for v in x[n:n + 5])
    return MLRPotential(De=De, Re=Re, C={6: C6, 8: C8, 10: C10}, beta=tuple(x[:n]), p=START["p"], q=START["q"],
                        Rref=START["Rref"], bo=X0)


def x_start(n):
    b = list(START["beta"]) + [0.0] * max(0, n - len(START["beta"]))
    return np.r_[b[:n], START["De"], START["Re"], START["C"]["6"], START["C"]["8"], START["C"]["10"], 0.0]


def targets(n):
    base = BSplineSolver(potential(x_start(n), n), MU, **GRID, joins=())
    e00 = base.energy(0, 0)
    tg = []
    for v in range(0, 18):
        for J in (0, 60, 120, 180):
            tg.append(dict(v=v, J=J, E=base.energy(v, J), s=ANCHOR_SIGMA / MHZ_PER_CM, kind="anchor"))
    for v in list(ATLAS) + list(HIGH):
        lo, hi = CORR.coverage["X"][v]
        Js = sorted({int(J) for J in np.linspace(lo, hi, 5 if hi - lo > 40 else 2)})
        u = CORR.uncertainty("X", v) or 20.0
        for J in Js:
            tg.append(dict(v=v, J=J, E=base.energy(v, J) + CORR.shift("X", v, J) / MHZ_PER_CM,
                           s=u / MHZ_PER_CM, kind="atlas" if v in ATLAS else "emission"))
    for v, c in sorted(martin_levels().items()):
        if v < 26 or v in HIGH or v in SKIP_MARTIN:
            continue
        Jmax = 120 if v <= 85 else 60 if v <= 96 else 20
        for J in (0, 20, 40, 60, 80, 100, 120):
            if J <= Jmax:
                tg.append(dict(v=v, J=J, E=e00 + martin_term(c, J), s=martin_sigma(J), kind="martin"))
    return tg, e00


class Fit:
    def __init__(self, n, tg):
        self.n, self.tg = n, tg
        self.mask = np.ones(len(tg), bool)
        self.Js = sorted({t["J"] for t in tg})
        self.names = ("De", "C6", "C8", "C10")
        self.pk = {"De": n, "C6": n + 2, "C8": n + 3, "C10": n + 4}
        self.martin = np.array([t["kind"] == "martin" for t in tg])

    def solve(self, x):
        s = BSplineSolver(potential(x, self.n), MU, **GRID, joins=())
        return s, {J: s.wavefunctions(J) for J in self.Js}

    def diffs(self, x, out):
        return np.array([out[t["J"]][0][t["v"]] - t["E"] for t in self.tg]) - np.where(self.martin, x[-1], 0.0)

    def resid(self, x):
        try:
            _, out = self.solve(x)
            d = (self.diffs(x, out) / np.array([t["s"] for t in self.tg]))[self.mask]
            pr = [(x[self.pk[p]] - PRIORS[p][0]) / PRIORS[p][1] for p in self.names] + [x[-1] / OFFSET_PRIOR]
            return np.nan_to_num(np.r_[d, pr], nan=1e7, posinf=1e7, neginf=1e7)
        except Exception:
            return np.full(int(self.mask.sum()) + 5, 1e7)

    def jac(self, x):
        s, out = self.solve(x)
        n5 = self.n + 5
        dV = []
        for k in range(n5):
            h = 1e-6 * max(abs(x[k]), 1.0)
            xp, xm = np.array(x, float), np.array(x, float)
            xp[k] += h
            xm[k] -= h
            dV.append((potential(xp, self.n)(s.R) - potential(xm, self.n)(s.R)) / (2 * h))
        dV = np.array(dV).T
        sg = np.array([t["s"] for t in self.tg])
        rows = np.array([(s.W * out[t["J"]][1][:, t["v"]] ** 2) @ dV for t in self.tg])
        rows = (np.c_[rows, -self.martin.astype(float)] / sg[:, None])[self.mask]
        pr = np.zeros((5, len(x)))
        for r, p in enumerate(self.names):
            pr[r, self.pk[p]] = 1.0 / PRIORS[p][1]
        pr[4, -1] = 1.0 / OFFSET_PRIOR
        return np.vstack([rows, pr])

    def run(self, x, nfev=40):
        n = self.n
        scale = np.r_[np.full(n, 0.05), 0.128, 0.0005, 0.12e6, 1.2e7, 0.5e8, 0.001]
        lo = np.r_[np.full(n, -np.inf), X0.De - 5, 2.6, 0.5e6, 1e6, 1e6, -0.05]
        hi = np.r_[np.full(n, np.inf), X0.De + 5, 2.75, 3.0e6, 1.5e8, 4.0e8, 0.05]
        for loss in ("linear", "soft_l1"):          # linear first: the start is up to 70 cm-1 off at high v''
            x = least_squares(self.resid, x, jac=self.jac, x_scale=scale, bounds=(lo, hi), max_nfev=nfev,
                              loss=loss, f_scale=3.0).x
        return x


def main(argv):
    opts = {a.split("=")[0]: a.split("=", 1)[1] for a in argv if "=" in a}
    n = int(opts.get("--nbeta", len(START["beta"])))
    tg, e00 = targets(n)
    if "--high" in opts:        # cm-1: sigma floor above v'' = 89, where the X-a'-a interactions near 5 A perturb
        for t in tg:            # the last levels and a single-channel potential cannot follow Martin's precision
            if t["kind"] == "martin" and t["v"] > 89:
                t["s"] = max(t["s"], float(opts["--high"]))
    f = Fit(n, tg)
    kind = np.array([t["kind"] for t in tg])
    v_of = np.array([t["v"] for t in tg])
    print(f"{len(tg)} targets: " + ", ".join(f"{np.sum(kind == k)} {k}" for k in ("anchor", "atlas", "emission", "martin")),
          flush=True)

    def report(x, label):
        _, out = f.solve(x)
        d = f.diffs(x, out) * MHZ_PER_CM
        rms = lambda a: float(np.sqrt(np.mean(a ** 2))) if a.size else float("nan")     # noqa: E731
        bands = " ".join(f"{lo}-{hi}:{rms(d[(kind == 'martin') & (v_of >= lo) & (v_of <= hi)]):.0f}"
                         for lo, hi in ((26, 47), (49, 60), (61, 75), (76, 89), (91, 108)))
        print(f"{label:<14} anchor {rms(d[kind == 'anchor']):.2f} | atlas {rms(d[kind == 'atlas']):.1f} | emission "
              f"{rms(d[kind == 'emission']):.1f} | Martin rms by v'' (MHz) {bands} | offset {x[-1] * MHZ_PER_CM:+.1f} MHz",
              flush=True)
        return d

    x0 = x_start(n)
    report(x0, "start (2026d)")
    x = x0
    stages = STAGES
    if "--from" in opts:        # continue from a saved stage, e.g. --from=89 --stages=92,95,98,101,104,108
        saved = json.load(open(OUT / f"mlr_x_martin_stage{opts['--from']}.json"))
        x = np.array(saved["x"])
        if saved["n"] < n:          # more beta: the new ones start at zero
            x = np.r_[x[:saved["n"]], np.zeros(n - saved["n"]), x[saved["n"]:]]
        stages = tuple(int(v) for v in opts["--stages"].split(","))
    soft = float(opts.get("--soft", 0))      # cm-1: first pass with the levels above v'' = 89 at this sigma
    if soft:
        hi = [k for k, t in enumerate(tg) if t["kind"] == "martin" and t["v"] > 89]
        keep = [tg[k]["s"] for k in hi]
        for k in hi:
            tg[k]["s"] = soft
        f.mask = np.ones(len(tg), bool)
        x = f.run(x)
        report(x, f"soft {soft} cm-1")
        for k, s0 in zip(hi, keep):
            tg[k]["s"] = s0
        json.dump(dict(n=n, x=list(map(float, x)), stage="soft"), open(OUT / "mlr_x_martin_stage_soft.json", "w"))
    for vmax in stages:
        f.mask = (kind != "martin") | (v_of <= vmax)
        x = f.run(x)
        report(x, f"stage v''<={vmax}")
        json.dump(dict(n=n, x=list(map(float, x)), stage=vmax), open(OUT / f"mlr_x_martin_stage{vmax}.json", "w"))
    f.mask = np.ones(len(tg), bool)
    d = report(x, "refit")
    by_v = {int(v): round(float(np.sqrt(np.mean(d[(v_of == v) & (kind == "martin")] ** 2))), 1)
            for v in sorted(set(v_of[kind == "martin"]))}
    print("Martin rms per level (MHz):", by_v)
    De, Re, C6, C8, C10 = (float(v) for v in x[n:n + 5])
    res = dict(START, De=De, Re=Re, C={"6": C6, "8": C8, "10": C10}, beta=list(map(float, x[:n])), id="mlr_x_martin",
               martin_offset_cm1=float(x[-1]), martin_rms_MHz={str(k): v for k, v in by_v.items()},
               fitted_to=("mlr_x_2026d refitted to Martin et al. 1986 Table I (v'' = 26-108, at the paper's precision of "
                          "recomputed lines, one origin offset), the Orsay-atlas levels v'' = 18-25 and the emission levels "
                          "v'' = 48, 53, 54 of level_corrections_2026j, anchored on its own levels v'' <= 17 at 3 MHz; "
                          "prototypes/mlr_x_martin.py"))
    OUT.mkdir(exist_ok=True)
    (OUT / "mlr_x_martin.json").write_text(json.dumps(res, indent=1))
    if "--write" in opts:
        name = opts["--write"]
        res.update(id=name, citation=f"i2spec, {name}: {res['fitted_to']}. Not a published parameter set.",
                   status="research result (docs/research/x-levels-martin1986.md)")
        for k in ("martin_rms_MHz",):
            res.pop(k)
        (ROOT / "src/i2spec/data" / f"{name}.json").write_text(json.dumps(res, indent=1))
        print("wrote", name)


if __name__ == "__main__":
    main(sys.argv[1:])
