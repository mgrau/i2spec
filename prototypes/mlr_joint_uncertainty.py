"""Lever 5, the parameter half: the joint fit's covariance propagated to any line.

Usage:  uv run python prototypes/mlr_joint_uncertainty.py        (about 10 minutes)
Output: prototypes/out/mlr_joint_uncertainty.json

Takes the joint fit's end point (prototypes/out/mlr_joint.json), rebuilds its weighted Jacobian at the
final floor, forms the Gauss-Newton covariance (J^T J)^-1 scaled by the reduced chi-squared of the rows
inside the robust loss's linear range, and propagates it to a grid of lines through
d nu / d theta = <psi_B|dV_B/d theta|psi_B> - <psi_X|dV_X/d theta|psi_X>. That is the part of a
line's uncertainty the error map cannot see: how far the fitted potentials are constrained where no
observation sits. It is meaningful only once the fit has converged; the report says whether it has.
"""

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from mlr_joint_fit import GRID, MU, Joint, band_shifts, observations  # noqa: E402

from i2spec.bspline import BSplineSolver  # noqa: E402
from i2spec.constants import MHZ_PER_CM  # noqa: E402
from i2spec.model import RovibronicModel  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "prototypes/out"
#: Lines to report: (v', v'', J'') across the measured region and beyond it
PROBES = [(v, vl, J) for v in (0, 10, 20, 32, 40, 44, 50, 58, 62) for vl in (0, 5, 12, 17, 30, 42, 54) for J in (10, 60, 120)]


def main():
    fit_result = json.loads((OUT / "mlr_joint.json").read_text())
    nx, nb, x = fit_result["nx"], fit_result["nb"], np.array(fit_result["x"])
    model = RovibronicModel("127I2", grids=GRID)
    rows = observations() + band_shifts(model)
    fit = Joint(rows, nx, nb)
    fit.set_floor(1.0)
    r = fit.residuals(x)
    J = fit.jacobian(x)
    n_data = len(rows)
    linear = np.abs(r[:n_data]) < 5.0                       # inside soft_l1's linear range at f_scale 5
    chi2_red = float(np.mean(r[:n_data][linear] ** 2))
    print(f"{n_data} rows; {linear.sum()} inside the linear range; reduced chi2 there {chi2_red:.2f} "
          f"(1 means the 1 MHz floor is the model's error)")
    JtJ = J.T @ J
    cov = np.linalg.pinv(JtJ) * chi2_red
    sig = np.sqrt(np.diag(cov))
    names = [f"bX{i}" for i in range(nx)] + ["X_De", "X_Re", "X_C6", "X_C8", "X_C10"] + \
            [f"bB{i}" for i in range(nb)] + ["B_Te", "B_Re", "B_C5", "B_C6", "B_C8", "B_C10"]
    print("parameter sigmas (De, Re, C in their own units):",
          {n: f"{s:.3g}" for n, s in zip(names, sig) if not n.startswith("b")})
    # propagate to probe lines
    X, B = fit.potentials(x)
    sx = BSplineSolver(X, MU, **GRID["X"], joins=())
    sb = BSplineSolver(B, MU, **GRID["B"], joins=())
    dVx, dVb = [], []
    for k in range(len(x)):
        h = 1e-6 * max(abs(x[k]), 1.0)
        xp, xm = x.copy(), x.copy()
        xp[k] += h
        xm[k] -= h
        Xp, Bp = fit.potentials(xp)
        Xm, Bm = fit.potentials(xm)
        dVx.append((Xp(sx.R) - Xm(sx.R)) / (2 * h))
        dVb.append((Bp(sb.R) - Bm(sb.R)) / (2 * h))
    dVx, dVb = np.array(dVx), np.array(dVb)
    probes = []
    cache_x, cache_b = {}, {}
    for v, vl, Jl in PROBES:
        Ju = Jl + 1
        if Jl not in cache_x:
            cache_x[Jl] = sx.wavefunctions(Jl)
        if Ju not in cache_b:
            cache_b[Ju] = sb.wavefunctions(Ju)
        ex, px = cache_x[Jl]
        eb, pb = cache_b[Ju]
        if vl >= len(ex) or v >= len(eb):
            continue
        g = dVb @ (sb.W * pb[:, v] ** 2) - dVx @ (sx.W * px[:, vl] ** 2)       # d nu / d theta, cm-1
        s = float(np.sqrt(g @ cov @ g)) * MHZ_PER_CM
        probes.append(dict(v_upper=v, v_lower=vl, J_lower=Jl, sigma_MHz=s))
    (OUT / "mlr_joint_uncertainty.json").write_text(json.dumps(dict(chi2_red=chi2_red, sigma=dict(zip(names, sig.tolist())),
                                                                   probes=probes), indent=1))
    print(f"\n{'v′':>4}{'v″':>4}{'J″':>5}{'σ_param MHz':>13}")
    for p in probes:
        if p["J_lower"] == 60:
            print(f"{p['v_upper']:>4}{p['v_lower']:>4}{p['J_lower']:>5}{p['sigma_MHz']:>13.2f}")


if __name__ == "__main__":
    main()
