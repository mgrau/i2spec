"""Bridge the X potential across the unmeasured gap: v'' <= 19 from the atlas, v'' = 48-54 from measurement.

The Gerstenkorn & Luc constants stop at v'' = 19 (R <= 3.07 A) and the only measurements above that are
36 emission lines at v'' = 48, 53, 54. An MLR fitted to the first alone has an exponent nobody has
constrained beyond 3.07 A: with Rref = 1.25 Re its levels come out 408 cm-1 low at v'' = 48, and with
Rref = Re it diverges outright. Adding the 36 levels at their own weight is no better -- they outweigh
the 260 atlas levels 10^8 to one and the fit throws the well away to chase them.

So the weight comes down in steps: the levels enter at 100 cm-1 and end at their own uncertainty, and
the potential walks across the gap instead of jumping. Every stage is reported, so the walk is visible.

usage: x_bridge.py [--nbeta=14] [--rref=1.25] [--start=prototypes/out/global_fit_dunham_v1.json]
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

sys.path.insert(0, str(Path(__file__).parent))
from global_fit import GRIDS, precision_rows                                   # noqa: E402
from mlr_from_rkr import X_LONG_RANGE, mlr                                     # noqa: E402

from i2spec import dunham                                                      # noqa: E402
from i2spec.bspline import BSplineSolver                                       # noqa: E402
from i2spec.constants import MHZ_PER_CM, reduced_mass                          # noqa: E402
from i2spec.observations import Predictor                                      # noqa: E402
from i2spec.potentials import MLRPotential                                     # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
MU = reduced_mass("127I2")
STEPS = (100.0, 30.0, 10.0, 3.0, 1.0, 0.3, 0.1, 0.05)     # cm-1, the implied levels' weight


def implied_levels(b_params, nb):
    """E_X(v'', J'') = E_B(v', J') - nu for every measured line with v'' > 17."""
    pb = mlr("B", b_params, nb)
    sb = BSplineSolver(pb, MU, **GRIDS["B"], joins=())
    out = []
    for line, nu, sig, _ in precision_rows(Predictor()):
        if line.v_lower <= 17:
            continue
        Ju = line.J_lower + (1 if line.branch == "R" else -1)
        out.append((line.v_lower, line.J_lower, sb.energy(line.v_upper, Ju) - nu, float(np.hypot(sig, 0.05))))
    return out


def main(argv):
    opts = {a.split("=")[0]: a.split("=", 1)[1] for a in argv if "=" in a}
    n_beta = int(opts.get("--nbeta", 14))
    rref = float(opts.get("--rref", 1.25))
    start = json.loads(Path(opts.get("--start", ROOT / "prototypes/out/global_fit_dunham_v1.json")).read_text())
    nb = len(start["B"]["x"]) - 2
    x0 = np.array(start["X"]["x"][:1 + n_beta] + [0.0] * max(0, n_beta - (len(start["X"]["x"]) - 1)))

    def build(x):
        return MLRPotential(De=X_LONG_RANGE["De"], Re=x[0], C=X_LONG_RANGE["C"], beta=tuple(x[1:]),
                            p=5, q=3, Rref=rref * x[0])

    js = list(range(0, 121, 10))
    dun = [(v, J, dunham.energy("X", v, J) - dunham.energy("X", 0, 0)) for v in range(0, 20) for J in js]
    imp = implied_levels(start["B"]["x"], nb)
    print(f"{len(dun)} atlas levels (v'' <= 19) and {len(imp)} implied levels (v'' = "
          f"{sorted({a[0] for a in imp})}), {n_beta} beta, Rref = {rref} Re")
    all_js = sorted(set(js) | {a[1] for a in imp})

    def solve(x, want_psi=False):
        s = BSplineSolver(build(x), MU, **GRIDS["X"], joins=())
        lev, psi = {}, {}
        for J in all_js:
            if want_psi:
                e, v = s.wavefunctions(J)
                psi[J] = v
            else:
                e = s.levels(J)
            lev[J] = e
        return (lev, s, psi) if want_psi else lev

    v_all = np.array([t[0] for t in dun] + [a[0] for a in imp])
    j_all = np.array([t[1] for t in dun] + [a[1] for a in imp])
    val = np.array([t[2] for t in dun] + [a[2] for a in imp])
    n_rel = len(dun)
    rel = np.arange(len(val)) < n_rel

    def model(x, lev=None):
        lev = lev if lev is not None else solve(x)
        e = np.array([lev[J][v] for v, J in zip(v_all, j_all)])
        return e - np.where(rel, lev[0][0], 0.0)

    def make(sigma_imp):
        sig = np.where(rel, 0.002, sigma_imp)

        def resid(x):
            p = build(x)
            with np.errstate(over="ignore", invalid="ignore"):
                probe = p(np.arange(2.2, 6.0, 0.05))
            if not np.isfinite(probe).all() or np.abs(probe).max() > 1e6 or not (2.5 < x[0] < 2.8):
                return np.full(len(val), 1e4)
            try:
                return (model(x) - val) / sig
            except (ValueError, np.linalg.LinAlgError):
                return np.full(len(val), 1e4)

        def jac(x):
            lev, s, psi = solve(x, want_psi=True)
            J = np.zeros((len(val), len(x)))
            for i in range(len(x)):
                h = 1e-4 * max(abs(x[i]), 1.0)
                xp, xm = x.copy(), x.copy()
                xp[i] += h; xm[i] -= h
                dV = (build(xp)(s.R) - build(xm)(s.R)) / (2 * h)
                wdv = s.W * dV
                d = {k: (wdv[:, None] * psi[k] ** 2).sum(axis=0) for k in all_js}
                col = np.array([d[k][v] for v, k in zip(v_all, j_all)])
                J[:, i] = (col - np.where(rel, d[0][0], 0.0)) / sig
            return J

        return resid, jac

    x = x0
    print(f"\n{'sigma_imp':>10} {'atlas rms':>12} {'implied rms':>14}")
    for sigma_imp in STEPS:
        resid, jac = make(sigma_imp)
        res = least_squares(resid, x, jac=jac, x_scale=[0.005] + [0.2] * n_beta, max_nfev=60)
        x = res.x
        r = model(x) - val
        print(f"{sigma_imp:10.2f} {MHZ_PER_CM * np.sqrt(np.mean(r[rel] ** 2)):9.1f} MHz "
              f"{MHZ_PER_CM * np.sqrt(np.mean(r[~rel] ** 2)):11.1f} MHz", flush=True)
    out = ROOT / "prototypes/out/x_bridge.json"
    out.write_text(json.dumps({"x": list(map(float, x)), "n_beta": n_beta, "rref": rref}, indent=1))
    lev = solve(x)
    print("\nlevels above X(0,0), cm-1:")
    ref = None
    try:
        from mlr_evaluate import load_x
        sg = BSplineSolver(load_x(ROOT / "data/potentials/mlr_x_2026c.json"), MU, **GRIDS["X"], joins=())
        ref = {v: sg.energy(v, 0) - sg.energy(0, 0) for v in (19, 30, 42, 48, 54)}
    except Exception:
        pass
    for v in (19, 30, 42, 48, 54):
        e = lev[0][v] - lev[0][0]
        print(f"  v''={v:2d}: {e:10.3f}" + (f"   mlr_x_2026c {ref[v]:10.3f}  diff {e - ref[v]:+8.3f}" if ref else ""))
    print("wrote", out)


if __name__ == "__main__":
    main(sys.argv[1:])
