"""Build the Orsay Partie IV line lists from the transcribed pages (docs/research/orsay-atlas-19700-20035.md).

  data/atlas_lines/orsay1983_part4.csv           every plate line: N, sigma, eps, depth, plate
  data/atlas_lines/orsay1983_part4_assigned.csv  every assignment of the authors' classification with an
                                                 observed line: N, sigma, eps, depth, v', v'', J'', branch,
                                                 the authors' sigma_cal, the number of assignments of the line
                                                 (1 = unblended), and the i2spec model value and obs - model

usage: orsay4_dataset.py
"""
import csv
from pathlib import Path

from i2spec.constants import MHZ_PER_CM
from i2spec.model import RovibronicModel

ROOT = Path(__file__).resolve().parents[1] / "data" / "atlas_lines"


def main():
    plate = {}
    for p in sorted((ROOT / "orsay1983_pages").glob("p*.tsv")):
        for line in p.read_text().splitlines():
            n, s, e, i = line.split("\t")
            plate[int(n)] = (s, e, i, int(p.stem[1:]))
    with open(ROOT / "orsay1983_part4.csv", "w", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["N", "sigma_cm1", "eps_mk", "depth", "plate"])
        for n in sorted(plate):
            w.writerow([n, *plate[n]])

    rows = []
    for p in sorted((ROOT / "orsay1983_classification").glob("c*.tsv")):
        for line in p.read_text().splitlines():
            n, s, i, J, br, vu, vl, sc, omc = line.split("\t")
            if n and J:
                rows.append((int(n), int(vu), int(vl), int(J), br, sc))
    blend = {}
    for r in rows:
        blend[r[0]] = blend.get(r[0], 0) + 1
    model = RovibronicModel()
    with open(ROOT / "orsay1983_part4_assigned.csv", "w", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["N", "sigma_cm1", "eps_mk", "depth", "v_upper", "v_lower", "J_lower", "branch", "sigma_cal_GL",
                    "n_assignments", "model_cm1", "obs_minus_model_MHz"])
        for n, vu, vl, J, br, sc in rows:
            s, e, i, _ = plate.get(n, (None, "", "", None))
            if s is None:           # N = 3614, printed only in the classification
                continue
            try:
                m = model.transition(vu, vl, J, br)
                mm, d = f"{m:.5f}", f"{(float(s) - m) * MHZ_PER_CM:.1f}"
            except (ValueError, IndexError):
                mm = d = ""
            w.writerow([n, s, e, i, vu, vl, J, br, sc, blend[n], mm, d])
    print(f"{len(plate)} lines, {len(rows)} assignments of observed lines")


if __name__ == "__main__":
    main()
