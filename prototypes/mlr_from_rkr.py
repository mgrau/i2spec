"""Fit MLR potentials to the RKR curves of the Gerstenkorn & Luc constants: the start of our own fit.

The RKR turning points (i2spec.rkr) are the potential that produced the atlas's Dunham constants. An
MLR fitted to them is the same curve in a form that extrapolates on theory rather than on a series:
V(R) = De [1 - u(R)/u(Re) exp(-beta(R) y_p^eq)]^2 with the long-range u(R) = sum C_n/R^n held at its
measured values. Nothing here comes from the Hannover parameter sets.

X: De and C6, C8, C10 from Bacis, Cerny & Martin 1986 (the source Knoeckel 2004 Table 4 itself cites).
B: the asymptote is De(X) + 7602.9762 cm-1 (the I(2P3/2) + I(2P1/2) limit), so Te + De(B) is fixed;
   C5 is the quadrupole-quadrupole term of the 0u+ state and C6, C8 follow Gerstenkorn, Luc & Amiot 1985.

usage: mlr_from_rkr.py [--n-beta-x=10] [--n-beta-b=12] [--write]
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

from i2spec import dunham, rkr
from i2spec.constants import ATOMIC_SPLITTING
from i2spec.potentials import MLRPotential

ROOT = Path(__file__).resolve().parents[1]

#: X state long range, Bacis, Cerny & Martin, J. Mol. Spectrosc. 118, 434 (1986): De = 12547.335(128),
#: C6 = 1.48(12)e6, C8 = 3.86(1.20)e7, C10 = 1.0(0.5)e8 cm-1 Angstrom^n.
X_LONG_RANGE = dict(De=12547.335, C={6: 1.48e6, 8: 3.86e7, 10: 1.0e8})
#: B state long range, Gerstenkorn, Luc & Amiot, J. Phys. (Paris) 46, 355 (1985), "Analysis of the long
#: range potential of iodine in the B 3Pi(0u+) state": the C5 (quadrupole-quadrupole), C6, C8 and C10 of
#: that analysis, as tabulated by Salumbides et al. 2008. Literature values, not fitted here.
B_LONG_RANGE = dict(C={5: 3.161e5, 6: 1.506e6, 8: 2.48e7, 10: 4.20e8})
#: RKR near dissociation is unreliable (first order, and the constants themselves are extrapolating):
#: the B fit uses v' <= 70, which reaches 4360 of the 4381 cm-1 well depth.
B_V_MAX = 70


#: Long-range parameters may be fitted as fractional corrections to the published values, with those
#: papers' own uncertainties as priors: De 128/12547, C6 8%, C8 31%, C10 50% (Bacis 1986). They are the
#: only handle that reaches v'' = 48-54 while leaving v'' <= 17 alone -- the orthogonal direction
#: docs/design/fitting.md found the series coefficients could not supply.
X_LR_PRIOR = {"De": 128.0 / 12547.335, "C6": 0.081, "C8": 0.311, "C10": 0.50}
#: R_ref of each state as a multiple of R_e. It is not a detail: with 1.25 the X potential extrapolates
#: smoothly but low, with 1.0 it fits the well better and extrapolates wildly, and mlr_x_2026c -- which
#: reaches 2-4 MHz on the precision sets where this fit reaches 16 -- uses 1.0. Settable to test that.
RREF = {"X": float(__import__("os").environ.get("I2SPEC_X_RREF", "1.25")),
        "B": float(__import__("os").environ.get("I2SPEC_B_RREF", "1.5"))}
#: The polynomial basis of the MLR exponent. Powers of y on [-1, 1] are nearly collinear and leave the
#: fit with a Jacobian of condition number 1e17; Chebyshev spans the same functions and is resolvable.
BASIS = __import__("os").environ.get("I2SPEC_BASIS", "power")
X_LR_KEYS = ("De", "C6", "C8", "C10")


def mlr(state, params, n_beta, bo=None, lr=None):
    """Build an MLR from the fitted vector: (Re, beta_0..beta_n) for X, (Te, Re, beta...) for B.

    ``lr``: fractional corrections (De, C6, C8, C10) to the X long-range parameters, or None.
    """
    if state == "X":
        Re, beta = params[0], params[1:]
        f = dict(zip(X_LR_KEYS, lr)) if lr is not None else {}
        De = X_LONG_RANGE["De"] * (1 + f.get("De", 0.0))
        C = {n: X_LONG_RANGE["C"][n] * (1 + f.get(f"C{n}", 0.0)) for n in X_LONG_RANGE["C"]}
        return MLRPotential(De=De, Re=Re, C=C, beta=tuple(beta), p=5, q=3, Rref=RREF["X"] * Re, bo=bo,
                            basis=BASIS)
    Te, Re, beta = params[0], params[1], params[2:]
    De_x = X_LONG_RANGE["De"] * (1 + (lr[0] if lr is not None else 0.0))
    return MLRPotential(De=De_x + ATOMIC_SPLITTING - Te, Re=Re, C=B_LONG_RANGE["C"],
                        beta=tuple(beta), p=6, q=4, Rref=RREF["B"] * Re, Te=Te, bo=bo, basis=BASIS)


#: With R_ref = R_e the exponent is centred on the well and the fit is far better there, but nothing
#: constrains it beyond the RKR range and it can diverge. These points hold it to the theoretical tail,
#: V(R) -> De - u(R), which is what the MLR is built to do anyway: they cost nothing where the RKR
#: curve already speaks and keep the extrapolation physical where it does not.
#: only where the asymptotic form is actually the potential: at 5.5 A the dispersion sum is 0.6 % of
#: De for X, at 3.2 A it is 11 % and the true curve is nowhere near it -- points there wreck the fit.
TAIL_R = np.arange(5.5, 12.01, 0.5)
TAIL_WEIGHT = 5.0                      # cm-1: loose enough to bend, tight enough to stop a divergence


def tail_points(state, De):
    C = X_LONG_RANGE["C"] if state == "X" else B_LONG_RANGE["C"]
    u = sum(c / TAIL_R ** n for n, c in C.items())
    return TAIL_R, De - u


def fit_state(state, n_beta, v_max=None):
    R, V = rkr.curve(state, v_max=v_max)
    Re0 = rkr.re(state)
    Te0 = dunham.constants()["T00"] - dunham.G("B", -0.5) + dunham.G("X", -0.5) if state == "B" else 0.0
    x0 = ([Re0] if state == "X" else [Te0, Re0]) + [0.0] * n_beta
    x0[len(x0) - n_beta] = -1.0                       # beta_0: a sensible scale for both states
    w = 1.0 / np.maximum(np.sqrt(np.maximum(V, 1.0)), 1.0)     # relative weight: the well bottom matters most

    De = X_LONG_RANGE["De"] if state == "X" else None

    def resid(x):
        p = mlr(state, x, n_beta)
        model = p(R) - (0.0 if state == "X" else x[0])
        rt, vt = tail_points(state, p.De)
        with np.errstate(over="ignore", invalid="ignore"):
            tail = p(rt) - (0.0 if state == "X" else x[0])
        tail = np.where(np.isfinite(tail), tail, 1e6)
        return np.concatenate([(model - V) * w, (tail - vt) / TAIL_WEIGHT])

    res = least_squares(resid, x0, x_scale=[0.01] + ([1.0] if state == "B" else []) + [0.5] * n_beta,
                        max_nfev=2000)
    p = mlr(state, res.x, n_beta)
    d = p(R) - (0.0 if state == "X" else res.x[0]) - V
    return res.x, p, R, V, d


def main(argv):
    opts = {a.split("=")[0]: a.split("=", 1)[1] for a in argv if "=" in a}
    out = {}
    for state, key in (("X", "--n-beta-x"), ("B", "--n-beta-b")):
        n_beta = int(opts.get(key, 10 if state == "X" else 12))
        x, p, R, V, d = fit_state(state, n_beta, v_max=(B_V_MAX if state == "B" else None))
        rms = float(np.sqrt(np.mean(d**2)))
        print(f"{state}: {n_beta} beta, {len(R)} RKR points over {R.min():.3f}-{R.max():.3f} A, "
              f"V up to {V.max():.0f} cm-1")
        print(f"   rms {rms:.4f} cm-1, max |{np.abs(d).max():.4f}|; Re = {p.Re:.6f} A"
              + (f", Te = {x[0]:.4f} cm-1, De = {p.De:.3f}" if state == "B" else f", De = {p.De:.3f} (fixed)"))
        for lo, hi in ((0, 500), (500, 2000), (2000, 4000), (4000, 12600)):
            m = (V >= lo) & (V < hi)
            if m.any():
                print(f"      V {lo:5d}-{hi:5d} cm-1: n={m.sum():3d}  rms {np.sqrt(np.mean(d[m]**2)):.4f} cm-1")
        out[state] = dict(x=list(map(float, x)), n_beta=n_beta, rms=rms, Re=p.Re, De=p.De,
                          Te=(float(x[0]) if state == "B" else 0.0))
    if "--write" in argv:
        path = ROOT / "prototypes/out/mlr_from_rkr.json"
        path.write_text(json.dumps(out, indent=1))
        print("wrote", path)


if __name__ == "__main__":
    main(sys.argv[1:])
