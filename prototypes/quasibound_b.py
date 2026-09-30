"""The B levels at or above the asymptote that the Orsay atlas Partie IV measured: are they narrow
quasi-bound resonances the box solver already gets right?

Usage:  uv run python prototypes/quasibound_b.py
Output: prototypes/out/quasibound_b.json, and a table on stdout

orsay4_fit.py leaves out every line whose upper level lies within 0.3 cm-1 of the B asymptote or above it:
70 unique lines, v' = 63-77, J' = 86-96. Behind the centrifugal barrier of such J' those levels are
resonances. For each one this computes, on the bare extended B potential of the default set:
  * the eigenvalues of the model's 40 A graded box and of an 80 A box graded the same way, keeping the
    states localised inside the barrier (> 90 % of the probability at R < R_barrier) and numbering v by
    them, so box states of the continuum outside the barrier cannot shift the count;
  * the barrier top, and a WKB estimate of the tunnelling width,
      Gamma = (dG / 2 pi) exp(-2 theta),  theta = int sqrt(2 mu (V_eff - E)) / hbar dR over the barrier,
    with dG the local vibrational spacing;
  * the model's own level energy (model.energy, which counts every box state) for comparison.
"""
import csv
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "prototypes"))

from i2spec.bspline import BSplineSolver                      # noqa: E402
from i2spec.constants import HBAR2_2U, MHZ_PER_CM, reduced_mass  # noqa: E402
from i2spec.model import grid_for                            # noqa: E402
from i2spec.observations import Predictor                    # noqa: E402
from mlr_b_dissociation import ASYM                          # noqa: E402

GRID40 = grid_for("B", 90)
GRID80 = dict(GRID40, rmax=80.0, nlev=GRID40["nlev"] + 40)


def main():
    m = Predictor().model("127I2")
    pot = m.extended_potentials["B"]
    mu = reduced_mass("127I2")
    rows = [r for r in csv.DictReader(open(ROOT / "data/atlas_lines/orsay1983_part4_assigned.csv"))
            if r["n_assignments"] == "1"]
    targets = {}
    for r in rows:
        vl, J, br, vu = int(r["v_lower"]), int(r["J_lower"]), r["branch"], int(r["v_upper"])
        Eu = float(r["sigma_cm1"]) + m.energy("X", vl, J)
        if Eu > ASYM - 0.3:
            targets.setdefault((vu, J + (1 if br == "R" else -1)), []).append((r, Eu))
    Js = sorted({J for _, J in targets})
    solvers = {k: BSplineSolver(pot, mu, **g) for k, g in (("40", GRID40), ("80", GRID80))}
    R = np.linspace(3.0, 30.0, 27001)

    def veff(J):
        return pot(R) + pot.adiabatic(R, 1.0) + HBAR2_2U / mu * (1 + pot.nonadiabatic(R, 1.0)) * J * (J + 1) / R**2

    def inner_levels(s, J, Rb):
        E, psi = s.wavefunctions(J)
        w = s.W[:, None] * psi**2
        inside = (w * (s.R[:, None] < Rb)).sum(axis=0) / w.sum(axis=0)
        return np.asarray(E)[inside > 0.9], np.asarray(E)

    out = []
    for J in Js:
        V = veff(J)
        i = int(np.argmax(np.where(R > 4.5, V, -np.inf)))
        Rb, Vb = float(R[i]), float(V[i])
        lev = {k: inner_levels(s, J, Rb) for k, s in solvers.items()}
        for (vu, Jp), lines in sorted(targets.items()):
            if Jp != J:
                continue
            inner40, all40 = lev["40"]
            inner80, _ = lev["80"]
            if vu >= len(inner40) or vu >= len(inner80):
                out.append(dict(v=vu, J=J, note="not localised inside the barrier"))
                continue
            E40, E80 = float(inner40[vu]), float(inner80[vu])
            dG = float(inner40[vu] - inner40[vu - 1])
            under = (R > 4.5) & (V > E40)
            seg = R[under]
            theta = float(np.trapezoid(np.sqrt(np.maximum(V[under] - E40, 0) * mu / HBAR2_2U), seg)) if seg.size else 0.0
            gamma = dG / (2 * np.pi) * np.exp(-2 * theta) * MHZ_PER_CM
            model_E = m.energy("B", vu, J)
            box_index = int(np.argmin(np.abs(all40 - E40)))
            obs = [(float(r["sigma_cm1"]), Eu) for r, Eu in lines]
            out.append(dict(v=vu, J=J, E_minus_asym_cm1=E40 - ASYM, barrier_minus_E_cm1=Vb - E40, R_barrier=Rb,
                            box_40_vs_80_MHz=(E40 - E80) * MHZ_PER_CM, wkb_width_MHz=gamma,
                            model_index_matches=box_index == vu,
                            model_minus_inner_MHz=(model_E - E40) * MHZ_PER_CM,
                            n_lines=len(lines)))
        print(f"  J' = {J} done", flush=True)
    (ROOT / "prototypes/out").mkdir(exist_ok=True)
    (ROOT / "prototypes/out/quasibound_b.json").write_text(json.dumps(out, indent=1))
    print(f"\n{'v':>3} {'J':>4} {'E-asym':>8} {'Vb-E':>8} {'40-80 Å':>10} {'width':>10} {'model idx':>9}  (cm-1, MHz)")
    for o in out:
        if "note" in o:
            print(f"{o['v']:3d} {o['J']:4d}  {o['note']}")
            continue
        print(f"{o['v']:3d} {o['J']:4d} {o['E_minus_asym_cm1']:8.3f} {o['barrier_minus_E_cm1']:8.2f} "
              f"{o['box_40_vs_80_MHz']:10.4f} {o['wkb_width_MHz']:10.3g} {str(o['model_index_matches']):>9}")


if __name__ == "__main__":
    main()
