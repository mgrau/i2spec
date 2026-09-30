"""Refit the MLR X potential to the measurements themselves, not to the published levels.

Usage:  uv run python prototypes/mlr_x_refit.py [n_beta] [max_v_upper]     (about 10 minutes)
Output: prototypes/out/mlr_x_refit.json

prototypes/mlr_x_fit.py fitted the published X levels for v'' <= 17, so the MLR inherited the published
errors there and added its own. This fits the 890 observations that constrain the potentials -- every
frequency and every interval between different lines -- with the B state, the hyperfine constants and
each line's hyperfine offset held fixed.

Lines above v' = MAX_V_UPPER are dropped: the published B potential is out by up to 2.2 cm-1 there
(bipm2005a's R(98) 58-1 is 5.4 GHz out), and an X-only fit would bend the X potential to absorb it.
"""

import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

from i2spec.bspline import BSplineSolver
from i2spec.constants import MHZ_PER_CM, reduced_mass
from i2spec.model import RovibronicModel
from i2spec.observations import Predictor, load_all
from i2spec.potentials import MLRPotential, load_potentials

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "prototypes/out"
X = load_potentials()["X"]
MU = reduced_mass("127I2")
GRID = dict(rmin=2.10, rmax=6.0, h=0.01, order=10, nlev=70)
PRIORS = {"De": (X.De, 0.128), "C6": (X.C[6], 0.12e6), "C8": (X.C[8], 1.20e7), "C10": (X.C[10], 0.5e8)}
MAX_V_UPPER = 44
#: Floor added to each measurement uncertainty, MHz. The published model reaches a few MHz at best on
#: these lines, and without a floor the 5 kHz BIPM rows would set the whole fit.
FLOOR = 0.5


def constraints(max_v_upper=MAX_V_UPPER):
    """Observations that constrain the potentials, with their hyperfine offsets frozen."""
    pred = Predictor()
    model = RovibronicModel("127I2")
    rows, skipped = [], 0
    for ds in load_all():
        scale = MHZ_PER_CM if ds.meta["unit"] == "cm-1" else 1.0
        for o in ds.observations:
            ref = o.ref_line or o.line
            if o.kind == "interval" and ref == o.line:
                continue                                   # intra-line: hyperfine only (Stage 2's data)
            lines = [o.line] + ([ref] if o.kind == "interval" else [])
            if any(L.isotopologue != "127I2" for L in lines) or any(L.v_upper > max_v_upper for L in lines):
                skipped += 1
                continue
            offset = pred.position(o.line, o.component) - pred.position(o.line)
            if o.kind == "interval":
                offset -= pred.position(ref, o.ref_component) - pred.position(ref)
            rows.append(dict(dataset=ds.id, value=o.value * scale, sigma=o.uncertainty * scale, offset=offset,
                             line=(o.line.v_upper, o.line.J_lower + (1 if o.line.branch == "R" else -1),
                                   o.line.v_lower, o.line.J_lower),
                             ref=None if o.kind == "frequency" else
                             (ref.v_upper, ref.J_lower + (1 if ref.branch == "R" else -1), ref.v_lower, ref.J_lower)))
    e_b = {}
    for r in rows:
        for key in (r["line"], r["ref"]):
            if key is not None:
                e_b.setdefault((key[0], key[1]), model.energy("B", key[0], key[1]))
    return rows, e_b, skipped


class Refit:
    def __init__(self, rows, e_b, n_beta, p=5, q=3):
        self.rows, self.e_b, self.n, self.p, self.q = rows, e_b, n_beta, p, q
        self.Js = sorted({k[3] for r in rows for k in (r["line"], r["ref"]) if k is not None})
        self.index = {J: i for i, J in enumerate(self.Js)}
        self.sigma = np.hypot([r["sigma"] for r in rows], FLOOR)

    def potential(self, x):
        De, Re, C6, C8, C10 = (float(v) for v in x[self.n:])
        return MLRPotential(De=De, Re=Re, C={6: C6, 8: C8, 10: C10}, beta=tuple(x[:self.n]),
                            p=self.p, q=self.q, bo=X)

    def solve(self, x):
        s = BSplineSolver(self.potential(x), MU, **GRID, joins=())
        return s, {J: s.wavefunctions(J) for J in self.Js}

    def model_values(self, out):
        """Model frequency (MHz) of every row, with the frozen hyperfine offsets added."""
        v = np.empty(len(self.rows))
        for i, r in enumerate(self.rows):
            vu, Ju, vl, Jl = r["line"]
            nu = (self.e_b[(vu, Ju)] - out[Jl][0][vl]) * MHZ_PER_CM
            if r["ref"] is not None:
                rvu, rJu, rvl, rJl = r["ref"]
                nu -= (self.e_b[(rvu, rJu)] - out[rJl][0][rvl]) * MHZ_PER_CM
            v[i] = nu + r["offset"]
        return v

    def prior_residuals(self, x):
        return np.array([(x[self.n + k] - PRIORS[n][0]) / PRIORS[n][1]
                         for k, n in enumerate(("De", "Re", "C6", "C8", "C10")) if n in PRIORS])

    def residuals(self, x):
        try:
            _, out = self.solve(x)
            d = (self.model_values(out) - np.array([r["value"] for r in self.rows])) / self.sigma
            return np.concatenate([np.nan_to_num(d, nan=1e6, posinf=1e6, neginf=1e6), self.prior_residuals(x)])
        except Exception:
            return np.full(len(self.rows) + len(PRIORS), 1e6)

    def jacobian(self, x, step=1e-6):
        s, out = self.solve(x)
        dV = []
        for k in range(len(x)):
            h = step * max(abs(x[k]), 1.0)
            xp, xm = np.array(x, float), np.array(x, float)
            xp[k] += h
            xm[k] -= h
            dV.append((self.potential(xp)(s.R) - self.potential(xm)(s.R)) / (2 * h))
        dE = {J: np.array([[(s.W * out[J][1][:, v] ** 2) @ d for d in dV] for v in range(GRID["nlev"])])
              for J in self.Js}                                   # dE_X(v, J)/dtheta, cm-1 per unit
        rows = []
        for i, r in enumerate(self.rows):
            _, _, vl, Jl = r["line"]
            g = -dE[Jl][vl]
            if r["ref"] is not None:
                g = g + dE[r["ref"][3]][r["ref"][2]]
            rows.append(g * MHZ_PER_CM / self.sigma[i])
        prior = np.zeros((len(PRIORS), len(x)))
        for row, k in enumerate(k for k, n in enumerate(("De", "Re", "C6", "C8", "C10")) if n in PRIORS):
            prior[row, self.n + k] = 1.0 / PRIORS[list(PRIORS)[row]][1]
        return np.vstack([np.array(rows), prior])

    def report(self, x, label):
        _, out = self.solve(x)
        d = self.model_values(out) - np.array([r["value"] for r in self.rows])
        by = {}
        for r, e in zip(self.rows, d):
            by.setdefault(r["dataset"], []).append(e)
        print(f"{label}: {len(d)} rows, rms {np.sqrt(np.mean(d ** 2)):.3f} MHz, "
              f"max {np.abs(d).max():.1f} | priors " +
              " ".join(f"{n} {v:+.1f}σ" for n, v in zip(PRIORS, self.prior_residuals(x))), flush=True)
        return {k: float(np.sqrt(np.mean(np.square(v)))) for k, v in by.items()}


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 12
    max_v = int(sys.argv[2]) if len(sys.argv) > 2 else MAX_V_UPPER
    t0 = time.perf_counter()
    rows, e_b, skipped = constraints(max_v)
    print(f"{len(rows)} constraints ({skipped} dropped: not ¹²⁷I₂, or v' > {max_v}), "
          f"{len({k[3] for r in rows for k in (r['line'], r['ref']) if k})} distinct J″", flush=True)
    fit = Refit(rows, e_b, n)
    d = json.loads((ROOT / "data/potentials/mlr_x_2026a.json").read_text())
    x0 = np.r_[d["beta"], d["De"], d["Re"], d["C"]["6"], d["C"]["8"], d["C"]["10"]]
    before = fit.report(x0, "start (fitted to published levels)")
    scale = np.r_[np.full(n, 0.1), 0.128, 0.001, 0.12e6, 1.20e7, 0.5e8]
    res = least_squares(fit.residuals, x0, jac=fit.jacobian, x_scale=scale, max_nfev=60,
                        loss="soft_l1", f_scale=3.0)
    after = fit.report(res.x, "refitted to the observations")
    OUT.mkdir(exist_ok=True)
    json.dump(dict(n=n, max_v_upper=max_v, x=list(map(float, res.x)), before=before, after=after),
              open(OUT / "mlr_x_refit.json", "w"), indent=1)
    print(f"\n{'dataset':<16}{'start':>10}{'refit':>10}   (rms, MHz)")
    for k in sorted(before, key=lambda k: -before[k]):
        print(f"{k:<16}{before[k]:>10.3f}{after[k]:>10.3f}")
    print(f"{time.perf_counter() - t0:.0f} s")


if __name__ == "__main__":
    main()
