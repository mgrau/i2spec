"""Which X levels do the Nölleke 2018 lines come from? A ground-state combination-difference census.

Two lines from one upper level (v', J') to X (v'', J'-1) and (v'', J'+1) - R(J'-1) and P(J'+1) - are
separated by D2F''(v'', J') = X(v'', J'+1) - X(v'', J'-1), whatever the upper state is. So counting pairs
of listed lines separated by D2F'' for each candidate v'' (over J), against the same count with D2F''
scaled by 1 +- 3 % (a null with the same spacing distribution), says which lower levels the spectrum
comes from without any model of the upper state.
"""
import numpy as np
from pathlib import Path
from i2spec.model import RovibronicModel

ROOT = Path(__file__).resolve().parents[1]
MHZ = 29979.2458


def main():
    L = np.genfromtxt(ROOT / "data/external/nolleke_2018/iodine_atlas.csv", delimiter=",", skip_header=1)
    nu = np.sort(1e7 / L[~np.isfinite(L[:, 2]), 0])          # drop the lines flagged as near water lines
    tol = 60 / MHZ                                               # their ~50 MHz, per line
    m = RovibronicModel("127I2")
    print(f"{len(nu)} lines; tolerance {tol*MHZ:.0f} MHz")
    print(f"{'v':>3} {'pairs':>6} {'null':>7} {'excess':>7} {'z':>6}   J with most pairs")
    for v in range(0, 36):
        Js = np.arange(6, 220)
        d2 = np.array([m.energy("X", v, J + 1) - m.energy("X", v, J - 1) for J in Js])

        def count(scale):
            n, byJ = 0, np.zeros(len(Js), int)
            for k, d in enumerate(d2 * scale):
                i = np.searchsorted(nu, nu + d)
                i = np.clip(i, 1, len(nu) - 1)
                near = np.minimum(np.abs(nu[i] - (nu + d)), np.abs(nu[i - 1] - (nu + d)))
                c = int(np.sum(near < tol))
                n += c
                byJ[k] = c
            return n, byJ
        real, byJ = count(1.0)
        nulls = np.array([count(s)[0] for s in (0.97, 0.98, 1.02, 1.03)])
        z = (real - nulls.mean()) / max(nulls.std(), np.sqrt(nulls.mean()))
        top = Js[np.argsort(byJ)[::-1][:3]]
        print(f"{v:3d} {real:6d} {nulls.mean():7.1f} {real - nulls.mean():+7.1f} {z:6.1f}   {top}")


if __name__ == "__main__":
    main()
