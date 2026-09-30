"""Fit an MLR X potential to the well-determined levels plus the measured high-v'' ones.

Usage:  uv run python prototypes/mlr_x_fit.py [n_beta] [p] [q]      (about 4 minutes for 12 5 3)
Output: prototypes/out/mlr_x_fit_<n>_<pq>.json

The result of `12 5 3` is what data/potentials/mlr_x_2026a.json holds, and docs/design/mlr-x.md
reports. Targets: the published X levels for v'' <= 17 (3 MHz each), the X levels implied by
matyugin2012 (v'' = 48) and nesterenko2019 (v'' = 53, 54) through the published B levels (90 MHz
each), and the published long-range parameters with their own uncertainties. The Jacobian is analytic
by Hellmann-Feynman, so one eigensolve per J gives every derivative. The measured levels are walked in
from a 10 cm^-1 weight to their real one, with De held until the last stage: released early, the fit
falls into a local minimum with De pinned at a bound.
"""
import json, sys, time
import numpy as np
from scipy.optimize import least_squares, minimize_scalar
from i2spec.bspline import BSplineSolver
from i2spec.constants import MHZ_PER_CM, reduced_mass
from i2spec.potentials import MLRPotential, load_potentials

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "prototypes/out"
X = load_potentials()["X"]
MU = reduced_mass("127I2")
GRID = dict(rmin=2.10, rmax=6.0, h=0.01, order=10, nlev=70)
REF = BSplineSolver(X, MU, **GRID)
E00 = REF.energy(0, 0)
LOW_J = (0, 40, 80, 120, 160)
SIGMA_LOW, SIGMA_HIGH = 1e-4, 3e-3
def measured_x_levels():
    """X level energies implied by the high-v'' emission data: E_X = E_B(published) - nu(measured).

    The upper levels are v' = 32 and 33, which the BIPM 532 nm tables anchor, so the published B state
    carries them. Components of one line agree to 26-45 MHz, which is the hyperfine model's error there.
    """
    from i2spec.model import RovibronicModel
    from i2spec.observations import Predictor, load_dataset

    pred, model, by = Predictor(), RovibronicModel("127I2"), {}
    for name in ("nesterenko2019", "matyugin2012"):
        for o in load_dataset(ROOT / "data/observations" / name).observations:
            if o.kind != "frequency":
                continue
            E = model.energy("X", o.line.v_lower, o.line.J_lower) - model.energy("X", 0, 0)
            implied = E - (o.value - pred(o, "MHz")) / MHZ_PER_CM
            by.setdefault((o.line.v_lower, o.line.J_lower), []).append(implied)
    return [dict(v=v, J=J, E=float(np.mean(e)), n=len(e)) for (v, J), e in sorted(by.items())]


HIGH = measured_x_levels()
#: Published values and 1 sigma priors of the long-range parameters (docs/design/fitting.md:
#: Knoeckel 2004 Table 4). De is known to 0.128 cm-1, so it is not a free shape parameter.
PRIORS = {"De": (X.De, 0.128), "C6": (X.C[6], 0.12e6), "C8": (X.C[8], 1.20e7), "C10": (X.C[10], 0.5e8)}


class Fit:
    #: parameters after the beta coefficients
    EXTRA = ("De", "Re", "C6", "C8", "C10")

    def __init__(self, n_beta, p, q, Rref_factor=None, use_high=True):
        self.n, self.p, self.q = n_beta, p, q
        self.Rref_factor, self.use_high = Rref_factor, use_high
        self.targets = [dict(v=v, J=J, E=REF.energy(v, J) - E00, sigma=SIGMA_LOW, high=False)
                        for J in LOW_J for v in range(18)]
        if use_high:
            self.targets += [dict(v=h["v"], J=h["J"], E=h["E"], sigma=SIGMA_HIGH, high=True) for h in HIGH]
        self.Js = sorted({t["J"] for t in self.targets})
        self.sigma = np.array([t["sigma"] for t in self.targets])
        self.high = np.array([t["high"] for t in self.targets])

    def potential(self, x):
        De, Re, C6, C8, C10 = (float(v) for v in x[self.n:])
        Rref = None if self.Rref_factor is None else self.Rref_factor * Re
        return MLRPotential(De=De, Re=Re, C={6: C6, 8: C8, 10: C10}, beta=tuple(x[:self.n]),
                            p=self.p, q=self.q, Rref=Rref, bo=X)

    def prior_residuals(self, x):
        """The long-range parameters carry their published uncertainties, in units of sigma."""
        return np.array([(x[self.n + k] - PRIORS[name][0]) / PRIORS[name][1]
                         for k, name in enumerate(self.EXTRA) if name in PRIORS])

    def solve(self, x):
        pot = self.potential(x)
        s = BSplineSolver(pot, MU, **GRID, joins=())
        out = {J: s.wavefunctions(J) for J in self.Js}
        return pot, s, out

    def deviations(self, x, solved=None):
        pot, s, out = solved or self.solve(x)
        e0 = out[0][0][0]
        return np.array([out[t["J"]][0][t["v"]] - e0 - t["E"] for t in self.targets])

    def residuals(self, x):
        try:
            return np.concatenate([self.deviations(x) / self.sigma, self.prior_residuals(x)])
        except Exception:
            return np.full(len(self.targets) + len(PRIORS), 1e6)

    def jacobian(self, x, step=1e-6):
        """d(residual)/dx by Hellmann-Feynman: <psi|dV/dx|psi>, with dV/dx on the grid."""
        pot, s, out = self.solve(x)
        R, W = s.R, s.W
        dV = []
        for k in range(len(x)):
            h = step * max(abs(x[k]), 1.0)
            xp, xm = np.array(x, float), np.array(x, float)
            xp[k] += h
            xm[k] -= h
            dV.append((self.potential(xp)(R) - self.potential(xm)(R)) / (2 * h))
        rows = []
        for t in self.targets:
            e, psi = out[t["J"]]
            w = W * psi[:, t["v"]] ** 2
            w0 = W * out[0][1][:, 0] ** 2                     # the v=0, J=0 zero point is subtracted
            rows.append([(w @ d - w0 @ d) / t["sigma"] for d in dV])
        prior = np.zeros((len(PRIORS), len(x)))
        for row, (k, name) in enumerate((k, n) for k, n in enumerate(self.EXTRA) if n in PRIORS):
            prior[row, self.n + k] = 1.0 / PRIORS[name][1]
        return np.vstack([np.array(rows), prior])

    def report(self, x, label):
        d = self.deviations(x)
        pr = self.prior_residuals(x)
        low = ~self.high
        msg = (f"{label:<34} v''<=17: rms {np.sqrt(np.mean(d[low]**2))*MHZ_PER_CM:8.2f} MHz, "
               f"max {np.abs(d[low]).max()*MHZ_PER_CM:9.2f} MHz")
        if self.high.any():
            msg += (f" | measured: rms {np.sqrt(np.mean(d[self.high]**2)):8.4f} cm-1, "
                    f"max {np.abs(d[self.high]).max():7.4f}")
        msg += "  | priors " + " ".join(f"{n} {v:+.1f}σ" for n, v in zip([k for k in self.EXTRA if k in PRIORS], pr))
        print(msg, flush=True)
        return d

    def x0(self):
        R = np.linspace(2.35, 5.0, 600)
        Re0 = minimize_scalar(lambda r: float(X(np.array([r]))[0]), bracket=(2.5, 2.7, 2.9)).x
        def curve(b):
            Rref = None if self.Rref_factor is None else self.Rref_factor * Re0
            return MLRPotential(De=X.De, Re=Re0, C=X.C, beta=tuple(b), p=self.p, q=self.q, Rref=Rref, bo=X)(R) - X(R)
        b = least_squares(curve, np.r_[2.0, np.zeros(self.n - 1)], max_nfev=4000).x
        return np.r_[b, X.De, Re0, X.C[6], X.C[8], X.C[10]]

    hold_De = False

    def bounds(self):
        """Physical bounds: the long-range parameters stay positive and within about 5 sigma.

        While the measured levels are walked in, De is held: it is known to 0.128 cm-1, and letting it
        drift early sends the fit into a local minimum with De pinned at a bound.
        """
        w = 1e-7 if self.hold_De else 5.0
        lo = np.r_[np.full(self.n, -np.inf), X.De - w, 2.5, 0.5e6, 1e6, 1e6]
        hi = np.r_[np.full(self.n, np.inf), X.De + w, 2.8, 3.0e6, 1.5e8, 4.0e8]
        return lo, hi

    def run(self, x0=None, max_nfev=120, verbose=0):
        scale = np.r_[np.full(self.n, 0.1), 0.128, 0.001, 0.12e6, 1.20e7, 0.5e8]
        res = least_squares(self.residuals, self.x0() if x0 is None else x0, jac=self.jacobian,
                            x_scale=scale, max_nfev=max_nfev, verbose=verbose, bounds=self.bounds())
        return res.x


if __name__ == "__main__":
    n, p, q = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    rf = float(sys.argv[4]) if len(sys.argv) > 4 else None
    t0 = time.perf_counter()
    low = Fit(n, p, q, rf, use_high=False)
    x = low.run()
    both = Fit(n, p, q, rf, use_high=True)
    both.report(x, f"N={n} p={p} q={q} published only")
    # walk the measured levels in: from 10 cm-1 down to their real weight, each stage starting from the last
    x2 = x
    for sigma in (10.0, 3.0, 1.0, 0.3, 0.1, 0.03, 0.01, SIGMA_HIGH):
        stage = Fit(n, p, q, rf, use_high=True)
        for t in stage.targets:
            if t["high"]:
                t["sigma"] = sigma
        stage.sigma = np.array([t["sigma"] for t in stage.targets])
        stage.hold_De = sigma > SIGMA_HIGH
        x2[n] = np.clip(x2[n], *(X.De + np.r_[-1, 1] * (1e-7 if stage.hold_De else 5.0)))
        x2 = stage.run(x0=x2)
        both.report(x2, f"   measured at sigma {sigma:g} cm-1")
    pot = both.potential(x2)
    print(f"   De = {pot.De:.4f} (published {X.De}), Re = {pot.Re:.6f}, Rref = {pot.reference:.4f}, "
          f"C6 = {pot.C[6]:.4g}, C8 = {pot.C[8]:.4g}, C10 = {pot.C[10]:.4g}, {time.perf_counter()-t0:.0f} s", flush=True)
    OUT.mkdir(exist_ok=True)
    json.dump(dict(n=n, p=p, q=q, Rref_factor=rf, x=list(map(float, x2)), x_low=list(map(float, x))),
              open(OUT / f"mlr_x_fit_{n}_{p}{q}_{'Re' if rf is None else rf}.json", "w"), indent=1)
