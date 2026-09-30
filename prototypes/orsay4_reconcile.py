"""Cross-check the two tables of the Orsay atlas Partie IV (19 700-20 035 cm-1, Gerstenkorn & Luc 1983).

The plates (data/atlas_lines/orsay1983_pages/p*.tsv: N sigma eps I) and the classification
(data/atlas_lines/orsay1983_classification/c*.tsv: N sigma_mes I J branch v_up v_lo sigma_cal omc_mk) were
read independently from photographs. Every line of the classification is also on a plate, and its
sigma_cal must be the Gerstenkorn-Luc Dunham value (i2spec.dunham.transition). Lists every disagreement:

  plate/class  sigma differs by more than 1 in the last digit, or I differs
  omc          sigma_mes - sigma_cal differs from the printed o-c by more than 1.5e-4
  dunham       sigma_cal differs from the Dunham value by more than 3e-4 cm-1 (a misread digit or a
               J label paired with the wrong row); the value for J +- 1 is shown when it fits
  duplicate    the same transition printed on two rows with different sigma_cal
  missing      a plate N with no classification row

usage: orsay4_reconcile.py
"""
from collections import defaultdict
from pathlib import Path

from i2spec.dunham import transition

ROOT = Path(__file__).resolve().parents[1] / "data" / "atlas_lines"


def plates():
    out = {}
    for p in sorted((ROOT / "orsay1983_pages").glob("p*.tsv")):
        for line in p.read_text().splitlines():
            n, s, e, i = line.split("\t")
            out[int(n)] = (float(s), int(i), p.name)
    return out


def classification():
    rows = []
    for p in sorted((ROOT / "orsay1983_classification").glob("c*.tsv")):
        for k, line in enumerate(p.read_text().splitlines(), 1):
            f = line.split("\t")
            if len(f) != 9:
                raise SystemExit(f"{p.name}:{k}: {len(f)} fields: {line!r}")
            rows.append((f"{p.name}:{k}", f))
    return rows


def main():
    plate, rows = plates(), classification()
    seen, trans = set(), defaultdict(set)
    for where, (n, s, i, J, br, vu, vl, sc, omc) in rows:
        if n:
            n = int(n); seen.add(n)
            if n not in plate:
                print(f"noplate {where} N={n} {s}: not on any plate")
            elif abs(float(s) - plate[n][0]) > 1.5e-4 or int(i) != plate[n][1]:
                print(f"plate/class {where} N={n}: class {s} I={i}  plate {plate[n][0]:.4f} I={plate[n][1]} ({plate[n][2]})")
        if J:
            key = (br, int(vu), int(vl), int(J))
            trans[key].add((sc, where))
            d = transition(int(vu), int(vl), int(J), br) - float(sc)
            if abs(d) > 3e-4:
                alt = [f"J={j}: {transition(int(vu), int(vl), j, br) - float(sc):+.4f}" for j in (int(J) - 1, int(J) + 1)
                       if j >= 0 and abs(transition(int(vu), int(vl), j, br) - float(sc)) < 3e-4]
                print(f"dunham {where} {J} {br}({vu}-{vl}) {sc}: dunham - cal = {d:+.4f} {' '.join(alt)}")
            if n and omc and abs(float(s) - float(sc) - float(omc) / 1000) > 1.5e-4:
                print(f"omc {where} N={n}: {s} - {sc} vs printed {omc}")
    for key, v in trans.items():
        if len({sc for sc, _ in v}) > 1:
            print(f"duplicate {key[0]}({key[3]}) {key[1]}-{key[2]}: " + ", ".join(f"{sc} @{w}" for sc, w in sorted(v)))
    missing = sorted(set(plate) - seen - set(range(1, 7)))
    print(f"missing: {len(missing)} plate lines without a classification row: {missing[:40]}")
    print(f"{len(plate)} plate lines, {len(rows)} classification rows, {len(seen)} classified N")


if __name__ == "__main__":
    main()
