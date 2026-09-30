"""Fit MLR potentials for X and B together, to the observations and the atlas band shifts.

Usage:  uv run python prototypes/mlr_joint_fit.py [n_beta_X] [n_beta_B]      (about 30 minutes)
Output: prototypes/out/mlr_joint.json

This is Stage 1 with the Phase C forms. prototypes/mlr_x_refit.py had to drop every line above
v' = 44, because the published B potential is out by up to 2.2 cm-1 there and an X-only fit would bend
X to absorb it. Here B is fitted too, so those lines come back, and the Salami & Ross band shifts
(prototypes/out/band_shifts_salami_ross.txt, from band_shifts_salami_ross.py) add constraints up to
v' = 62 where no frequency measurement exists.

Held fixed: the hyperfine constants, each line's hyperfine offset, and the Born-Oppenheimer
corrections of both states. Tied: Te(B) + De(B) = De(X) + 7602.9762 cm-1, the atomic splitting.
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
X0, B0 = load_potentials()["X"], load_potentials()["B"]
MU = reduced_mass("127I2")
GRID = {"X": dict(rmin=2.10, rmax=6.0, h=0.01, order=10, nlev=70),
        "B": dict(rmin=2.35, rmax=8.0, h=0.01, order=10, nlev=70)}
SPLITTING = 7602.9762          # 2P(1/2) - 2P(3/2) of atomic iodine: Te(B) + De(B) = De(X) + this
PRIORS = {"X_De": (X0.De, 0.128), "X_C6": (X0.C[6], 0.12e6), "X_C8": (X0.C[8], 1.20e7),
          "X_C10": (X0.C[10], 0.5e8), "B_C5": (B0.C[5], 0.2e5), "B_C6": (B0.C[6], 0.3e6),
          "B_C8": (B0.C[8], 1.0e7), "B_C10": (B0.C[10], 2.0e8)}
#: Model-error floor, MHz, added to every measurement uncertainty in quadrature, and walked down.
#: The data run from 5 kHz to 150 MHz in stated precision while the model starts tens to hundreds of
#: MHz out, so with the true uncertainties the fit is owned by whichever 5 kHz row is furthest off --
#: R(98) 58-1 at 265 GHz was 5e5 sigma -- and the 150 MHz atlas rows that carry v' > 44 count for
#: nothing. A floor that starts at the model's own error and shrinks as the fit improves lets every
#: row vote in proportion to what the model can currently do (docs/design/hyperfine-fit.md has the
#: same lesson).
FLOORS = (50.0, 10.0, 3.0, 1.0)
BAND_SIGMA = 0.005             # cm-1: the band-shift method's resolution
MIN_CONTRAST = 1.9


def observations():
    """Frequencies and inter-line intervals of 127I2, with hyperfine offsets frozen."""
    pred = Predictor()
    rows = []
    for ds in load_all():
        scale = MHZ_PER_CM if ds.meta["unit"] == "cm-1" else 1.0
        for o in ds.observations:
            ref = o.ref_line or o.line
            lines = [o.line] + ([ref] if o.kind == "interval" else [])
            if (o.kind == "interval" and ref == o.line) or any(L.isotopologue != "127I2" for L in lines):
                continue
            offset = pred.position(o.line, o.component) - pred.position(o.line)
            if o.kind == "interval":
                offset -= pred.position(ref, o.ref_component) - pred.position(ref)
            key = lambda L: (L.v_upper, L.J_lower + (1 if L.branch == "R" else -1), L.v_lower, L.J_lower)  # noqa: E731
            rows.append(dict(kind=ds.id, value=o.value * scale, sigma=o.uncertainty * scale, offset=offset,
                             line=key(o.line), ref=None if o.kind == "frequency" else key(ref)))
    return rows


def band_shifts(model):
    """Atlas band shifts as constraints on one representative line per band and J range."""
    path = OUT / "band_shifts_salami_ross.txt"
    if not path.exists():
        print(f"  (no {path.name}; run prototypes/band_shifts_salami_ross.py for the v' > 44 constraints)")
        return []
    rows = []
    for vu, vl, s1, c1, n1, s2, c2, n2 in np.loadtxt(path):
        for shift, contrast, J in ((s1, c1, 25), (s2, c2, 75)):
            if not (np.isfinite(shift) and np.isfinite(contrast)) or contrast < MIN_CONTRAST:
                continue        # a NaN contrast is a band the shift search could not fit at all
            vu, vl = int(vu), int(vl)
            nu = (model.energy("B", vu, J + 1) - model.energy("X", vl, J)) * MHZ_PER_CM
            rows.append(dict(kind="atlas", value=nu + shift * MHZ_PER_CM, sigma=BAND_SIGMA * MHZ_PER_CM,
                             offset=0.0, line=(vu, J + 1, vl, J), ref=None))
    return rows


class Joint:
    def __init__(self, rows, nx, nb, rf_x=1.25, rf_b=1.5):
        self.rows, self.nx, self.nb = rows, nx, nb
        self.rf_x, self.rf_b = rf_x, rf_b
        self.JX = sorted({k[3] for r in rows for k in (r["line"], r["ref"]) if k is not None})
        self.JB = sorted({k[1] for r in rows for k in (r["line"], r["ref"]) if k is not None})
        self.set_floor(FLOORS[0])

    def set_floor(self, floor):
        self.floor = floor
        self.sigma = np.hypot([r["sigma"] for r in self.rows], floor)
        self.names = (["X_De", "X_Re", "X_C6", "X_C8", "X_C10"], ["B_Te", "B_Re", "B_C5", "B_C6", "B_C8", "B_C10"])

    def split(self, x):
        nx, nb = self.nx, self.nb
        return x[:nx], x[nx:nx + 5], x[nx + 5:nx + 5 + nb], x[nx + 5 + nb:]

    def potentials(self, x):
        bx, px, bb, pb = self.split(x)
        De, Re, C6, C8, C10 = (float(v) for v in px)
        Te, ReB, C5, C6b, C8b, C10b = (float(v) for v in pb)
        X = MLRPotential(De=De, Re=Re, C={6: C6, 8: C8, 10: C10}, beta=tuple(bx), p=5, q=3,
                         Rref=self.rf_x * Re, bo=X0)
        B = MLRPotential(De=De + SPLITTING - Te, Re=ReB, C={5: C5, 6: C6b, 8: C8b, 10: C10b}, beta=tuple(bb),
                         p=6, q=4, Rref=self.rf_b * ReB, Te=Te, bo=B0)
        return X, B

    def solve(self, x):
        X, B = self.potentials(x)
        sx = BSplineSolver(X, MU, **GRID["X"], joins=())
        sb = BSplineSolver(B, MU, **GRID["B"], joins=())
        return sx, sb, {J: sx.wavefunctions(J) for J in self.JX}, {J: sb.wavefunctions(J) for J in self.JB}

    def model_values(self, ox, ob):
        out = np.empty(len(self.rows))
        for i, r in enumerate(self.rows):
            vu, Ju, vl, Jl = r["line"]
            nu = (ob[Ju][0][vu] - ox[Jl][0][vl]) * MHZ_PER_CM
            if r["ref"] is not None:
                rvu, rJu, rvl, rJl = r["ref"]
                nu -= (ob[rJu][0][rvu] - ox[rJl][0][rvl]) * MHZ_PER_CM
            out[i] = nu + r["offset"]
        return out

    def prior_residuals(self, x):
        bx, px, bb, pb = self.split(x)
        values = dict(zip(self.names[0], px)) | dict(zip(self.names[1], pb))
        return np.array([(values[k] - m) / s for k, (m, s) in PRIORS.items()])

    def residuals(self, x):
        try:
            _, _, ox, ob = self.solve(x)
            d = (self.model_values(ox, ob) - np.array([r["value"] for r in self.rows])) / self.sigma
            return np.concatenate([np.nan_to_num(d, nan=1e7, posinf=1e7, neginf=1e7), self.prior_residuals(x)])
        except Exception:
            return np.full(len(self.rows) + len(PRIORS), 1e7)

    def jacobian(self, x, step=1e-6):
        sx, sb, ox, ob = self.solve(x)
        dVx, dVb = [], []
        for k in range(len(x)):
            h = step * max(abs(x[k]), 1.0)
            xp, xm = np.array(x, float), np.array(x, float)
            xp[k] += h
            xm[k] -= h
            Xp, Bp = self.potentials(xp)
            Xm, Bm = self.potentials(xm)
            dVx.append((Xp(sx.R) - Xm(sx.R)) / (2 * h))
            dVb.append((Bp(sb.R) - Bm(sb.R)) / (2 * h))
        dEx = {J: np.array([[(sx.W * ox[J][1][:, v] ** 2) @ d for d in dVx] for v in range(GRID["X"]["nlev"])])
               for J in self.JX}
        dEb = {J: np.array([[(sb.W * ob[J][1][:, v] ** 2) @ d for d in dVb] for v in range(GRID["B"]["nlev"])])
               for J in self.JB}
        rows = []
        for i, r in enumerate(self.rows):
            vu, Ju, vl, Jl = r["line"]
            g = dEb[Ju][vu] - dEx[Jl][vl]
            if r["ref"] is not None:
                rvu, rJu, rvl, rJl = r["ref"]
                g = g - dEb[rJu][rvu] + dEx[rJl][rvl]
            rows.append(g * MHZ_PER_CM / self.sigma[i])
        prior = np.zeros((len(PRIORS), len(x)))
        order = list(self.names[0]) + list(self.names[1])
        for row, k in enumerate(PRIORS):
            j = order.index(k)
            prior[row, (self.nx + j) if j < 5 else (self.nx + 5 + self.nb + j - 5)] = 1.0 / PRIORS[k][1]
        return np.vstack([np.array(rows), prior])

    def report(self, x, label):
        _, _, ox, ob = self.solve(x)
        d = self.model_values(ox, ob) - np.array([r["value"] for r in self.rows])
        by = {}
        for r, e in zip(self.rows, d):
            by.setdefault(r["kind"], []).append(e)
        rms = lambda v: float(np.sqrt(np.mean(np.square(v))))          # noqa: E731
        print(f"{label}: {len(d)} rows, rms {rms(d):.3f} MHz | priors " +
              " ".join(f"{k.split('_')[0]}{k.split('_')[1]} {v:+.1f}σ" for k, v in zip(PRIORS, self.prior_residuals(x))),
              flush=True)
        return {k: rms(v) for k, v in by.items()}


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    nx = int(args[0]) if args else 12
    nb = int(args[1]) if len(args) > 1 else 14
    t0 = time.perf_counter()
    model = RovibronicModel("127I2", grids=GRID)   # the atlas bands reach v' = 62
    rows = observations() + band_shifts(model)
    print(f"{len(rows)} constraints: {sum(r['kind'] != 'atlas' for r in rows)} measured frequencies and intervals, "
          f"{sum(r['kind'] == 'atlas' for r in rows)} atlas band shifts", flush=True)
    fit = Joint(rows, nx, nb)
    resume = "--resume" in sys.argv
    b_only = "--b-only" in sys.argv        # hold X at its start: B has 20 parameters against 989 rows, X is already good
    tag = next((a.split("=", 1)[1] for a in sys.argv if a.startswith("--tag=")), "")
    if resume or b_only:      # continue from the last run's end point, at the smallest floors, with more iterations
        x0 = np.array(json.loads((OUT / "mlr_joint.json").read_text())["x"])
        floors, nfev = (3.0, 1.0, 0.5), 60
    else:
        x = json.loads((OUT / f"mlr_x_fit_{nx}_53_1.25.json").read_text())["x"]
        b = json.loads((OUT / f"b_fit_{nb}_1.5.json").read_text())["x"]
        x0 = np.r_[x[:nx], x[nx:], b[:nb], b[nb:]]
        floors, nfev = FLOORS, 25
    before = fit.report(x0, "resumed" if resume else "start (X refit, B from the published levels)")
    scale = np.r_[np.full(nx, 0.1), 0.128, 0.001, 0.12e6, 1.20e7, 0.5e8,
                  np.full(nb, 0.1), 0.1, 0.001, 0.2e5, 0.3e6, 1.0e7, 2.0e8]
    # The floor walks down as the fit improves; soft_l1 at 5 sigma keeps a few far-off rows from
    # owning a stage without writing them off (docs/design/fitting.md, "Two ways the fit went wrong").
    x = x0
    free = np.ones(len(x0), bool)
    if b_only:
        free[:nx + 5] = False
    def with_fixed(xf):
        full = x0.copy(); full[free] = xf; return full
    for floor in floors:
        fit.set_floor(floor)
        res = least_squares(lambda xf: fit.residuals(with_fixed(xf)), x[free],
                            jac=lambda xf: fit.jacobian(with_fixed(xf))[:, free], x_scale=scale[free],
                            max_nfev=nfev, loss="soft_l1", f_scale=5.0)
        x = with_fixed(res.x)
        fit.report(x, f"  floor {floor:g} MHz" + (" (B only)" if b_only else ""))
    after = fit.report(x, "joint fit" if not b_only else "B-only fit")
    res = type("R", (), {"x": x})
    OUT.mkdir(exist_ok=True)
    json.dump(dict(nx=nx, nb=nb, x=list(map(float, x)), before=before, after=after),
              open(OUT / f"mlr_joint{tag}.json", "w"), indent=1)
    print(f"\n{'set':<18}{'start':>10}{'joint':>10}   (rms, MHz)")
    for k in sorted(before, key=lambda k: -before[k]):
        print(f"{k:<18}{before[k]:>10.3f}{after[k]:>10.3f}")
    print(f"{time.perf_counter() - t0:.0f} s")


if __name__ == "__main__":
    main()
