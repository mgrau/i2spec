"""Assign the Nölleke et al. (2018) 915-985 nm line list to B-X bands.

The list (data/external/nolleke_2018, CC BY 4.0) is 10 162 absorption lines at 50 MHz absolute
accuracy with no quantum numbers. It lies at 10 152-10 929 cm-1, where the absorbing levels are
v'' = 25-45: the gap no other measurement reaches, and the region where the published X potential
and the MLR X disagree by GHz.

A single line cannot be assigned -- the model's X levels are uncertain by 0.01-0.03 cm-1 there and the
observed lines are 0.07 cm-1 apart -- but a whole band can: its lines follow a parabola in J that no
accidental coincidence reproduces. For each candidate (v', v'') this slides the model's predicted
band over a range of origin shifts and counts intensity-weighted matches, and compares the best count
with what the same band gets against shifts far from it (the null). Bands that stand out get per-line
assignments at their fitted shift.

usage: nolleke_assign.py [--parameters=i2spec2026d] [--tol=0.01] [--max-shift=0.6]
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

from i2spec.constants import MHZ_PER_CM
from i2spec.intensity import intensity_model, master_line_list

ROOT = Path(__file__).resolve().parents[1]
LIST = ROOT / "data/external/nolleke_2018/iodine_atlas.csv"
OUT = ROOT / "prototypes/out/nolleke_bands.txt"
RANGE = (10152.0, 10929.0)


def observed():
    import csv
    rows = list(csv.DictReader(open(LIST)))
    nu = np.array([1e7 / float(r["Vacuum-wavelength identified line (nm)"]) for r in rows])
    a = np.array([float(r["Absorption"]) if r["Absorption"] else np.nan for r in rows])
    order = np.argsort(nu)
    return nu[order], a[order]


def main(argv):
    opts = {a.split("=")[0]: a.split("=", 1)[1] for a in argv if "=" in a}
    tol = float(opts.get("--tol", 0.01))
    max_shift = float(opts.get("--max-shift", 0.6))
    nu_obs, a_obs = observed()
    print(f"{len(nu_obs)} observed lines, {nu_obs[0]:.1f}-{nu_obs[-1]:.1f} cm-1")

    model = intensity_model("127I2")
    master = master_line_list(model, RANGE[0] - 5, RANGE[1] + 5, T_range=(290.0, 320.0), S_min=1e-30)
    lines = master.at(300.0, RANGE[0] - 5, RANGE[1] + 5, S_min=0.0)
    print(f"{len(lines)} model lines in range; bands: "
          f"{len(set(zip(lines.v_upper.tolist(), lines.v_lower.tolist())))}")

    shifts = np.arange(-max_shift, max_shift + 1e-9, tol / 4)
    rows = []
    for (vu, vl) in sorted(set(zip(lines.v_upper.tolist(), lines.v_lower.tolist()))):
        m = (lines.v_upper == vu) & (lines.v_lower == vl) & (lines.nu > RANGE[0]) & (lines.nu < RANGE[1])
        if m.sum() < 20:
            continue
        nu_b, S_b = lines.nu[m], lines.S[m]
        w = S_b / S_b.max()
        score = np.empty(len(shifts))
        for i, s in enumerate(shifts):
            idx = np.searchsorted(nu_obs, nu_b + s)
            near = np.minimum(np.abs(nu_obs[np.clip(idx, 0, len(nu_obs) - 1)] - (nu_b + s)),
                              np.abs(nu_obs[np.clip(idx - 1, 0, len(nu_obs) - 1)] - (nu_b + s)))
            score[i] = np.sum(w * (near < tol))
        k = int(np.argmax(score))
        # the null: the same band's score at shifts more than 0.1 cm-1 away from the peak
        far = np.abs(shifts - shifts[k]) > 0.1
        mu, sd = score[far].mean(), score[far].std()
        rows.append((vu, vl, shifts[k], score[k], mu, sd, (score[k] - mu) / max(sd, 1e-9), int(m.sum()),
                     float(np.median(nu_b))))
    rows.sort(key=lambda r: -r[6])
    with open(OUT, "w") as f:
        f.write("# v_upper v_lower shift_cm-1 score null_mean null_sd z n_lines nu_median\n")
        for r in rows:
            f.write(" ".join(f"{v:.6g}" for v in r) + "\n")
    print("\n  v'  v''     shift   score    null       z     n        nu")
    for r in rows[:25]:
        print(f"{r[0]:4d} {r[1]:4d} {r[2]:+9.4f} {r[3]:7.1f} {r[4]:7.1f} {r[6]:7.1f} {r[7]:5d} {r[8]:9.1f}")
    strong = [r for r in rows if r[6] > 5]
    print(f"\n{len(strong)} bands with z > 5, {len([r for r in rows if r[6] > 3])} with z > 3, out of {len(rows)}")


if __name__ == "__main__":
    main(sys.argv[1:])
