"""Condense the IodineSpec5 output we hold into line positions and intensities, for the private repository.

IodineSpec5 (Knöckel, Bodermann & Tiemann, Hannover; closed source) writes ``*_L.OUT`` (rotational lines:
FREQofPOT / FREQofNIR / FREQofDun, 3-4 decimals, and an intensity) and ``*_H.OUT`` (per rotational line at
6 decimals: the four hyperfine parameters of each state and every component with its intensity). Its output
is not redistributed. This writes, to $I2SPEC_DATA/predictions/iodinespec5 (default ../i2spec_data):

  lines.csv                     one row per rotational line: the three position columns, the 6-decimal
                                position and intensity of the H file where there is one, the L file's
                                intensity, and the hyperfine parameters
  hyperfine_components.csv.gz   every component: its offset from the rotational line (MHz) and intensity
  README.md                     what the columns are, and where the files came from

A line computed more than once (overlapping requests) is kept once; later files win, as in
iodinespec5_compare.py, and the spread between repeats is reported.

usage: condense_iodinespec5.py [output_dir]      default data/external/iodinespec5/output
"""
from __future__ import annotations

import gzip
import os
import re
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from iodinespec5_compare import ISO, _line, _unit, parse_L   # noqa: E402
from i2spec.constants import MHZ_PER_CM                       # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "data/external/iodinespec5/output"
DEST = Path(os.environ.get("I2SPEC_DATA", ROOT.parent / "i2spec_data")) / "predictions" / "iodinespec5"
PARAMS = ("eqQ_X", "Csr_X", "del_X", "d_X", "eqQ_B", "Csr_B", "del_B", "d_B",
          "eqQ127_X", "Csr127_X", "eqQ129_X", "Csr129_X", "eqQ127_B", "Csr127_B", "eqQ129_B", "Csr129_B")   # 127I129I


def temperature(text):
    m = re.search(r"@T=\s*([\d.]+)\s*K", text)
    return float(m.group(1)) if m else None


def intensities_L(path):
    """{line: intensity} of an L file (the INTENSITY column, read by its position under the header)."""
    text = path.read_text(errors="replace").splitlines()
    iso = ISO[text[0].split()[1]]
    header = next(l for l in text if "FREQofPOT" in l)
    edge = header.index("INTENSITY") + len("INTENSITY")
    out = {}
    for l in text[text.index(header) + 3:]:
        m = re.match(r"\s*(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s", l)
        if not m:
            continue
        line = _line(iso, *map(int, m.groups()))
        for tok in re.finditer(r"\S+", l):
            if abs(tok.end() - edge) <= 2 and re.fullmatch(r"-?[\d.]+E[-+]\d+", tok.group()):
                out[line] = float(tok.group())
    return out


def parse_H(path):
    """Blocks: dict(line, freq MHz, intensity, params, comps [(label, F"-J", MHz, intensity)])."""
    text = path.read_text(errors="replace").splitlines()
    iso = ISO[text[0].split()[1]]
    scale = _unit(next(l for l in text if "@T=" in l))
    blocks, i = [], 0
    while i < len(text):
        if text[i].strip() == "rotational line:":
            Jp, vp, Js, vs, f, _, inten = text[i + 1].split()[:7]
            block = dict(line=_line(iso, int(Jp), int(vp), int(Js), int(vs)), freq=float(f) * scale,
                         intensity=float(inten), params={}, comps=[])
            for l in text[i + 3: i + 5]:
                for m in re.finditer(r"(eqQ|Csr|del|d)(127|129)?\s*\(([XB])\)\s*=\s*(-?[\d.]+)\s*(MHz|kHz)", l):
                    block["params"][f"{m[1]}{m[2] or ''}_{m[3]}"] = (float(m[4]), m[5])
            j = i + 7
            while j < len(text) and re.match(r"\s+[a-z]\s*\d+\s+-?\d+\s+", text[j]):
                m = re.match(r"\s+([a-z])\s*(\d+)\s+(-?\d+)\s+(-?[\d.]+)\s+[\d.]+\s+(-?[\d.]+E[-+]\d+|NaN)", text[j])
                block["comps"].append((f"{m[1]}{m[2]}", int(m[3]), float(m[4]) * scale, float(m[5])))
                j += 1
            blocks.append(block)
            i = j
        else:
            i += 1
    return blocks


def main(argv):
    src = Path(argv[0]) if argv else SRC
    rows = defaultdict(dict)                # line -> columns
    comps = {}
    repeats = defaultdict(list)             # line -> H positions (MHz) of every repeat
    units = {}
    for p in sorted(src.glob("*_L.OUT")):
        T = temperature(p.read_text(errors="replace"))
        inten = intensities_L(p)
        for line, pot, nir, dun in parse_L(p):
            if line is None:
                continue
            r = rows[line]
            r.update(pot=pot, nir=nir, dun=dun, L_file=p.name)
            if line in inten:
                r.update(L_intensity=inten[line], L_T=T)
    for p in sorted(src.glob("*_H.OUT")):
        T = temperature(p.read_text(errors="replace"))
        for b in parse_H(p):
            if b["line"] is None:
                continue
            repeats[b["line"]].append(b["freq"])
            r = rows[b["line"]]
            r.update(H=b["freq"], H_intensity=b["intensity"], H_T=T, H_file=p.name)
            for k, (v, u) in b["params"].items():
                r[k] = v
                units.setdefault(k, u)
            comps[b["line"]] = (b["freq"], b["comps"], T)

    spread = [max(v) - min(v) for v in repeats.values() if len(v) > 1]
    print(f"{len(rows)} rotational lines ({sum(1 for r in rows.values() if 'H' in r)} with hyperfine structure), "
          f"{sum(len(c[1]) for c in comps.values())} components; {len(spread)} lines computed more than once, "
          f"largest difference between repeats {max(spread, default=0):.4f} MHz")
    DEST.mkdir(parents=True, exist_ok=True)

    def f(x, fmt):
        return "" if x is None else fmt % x

    order = sorted(rows, key=lambda l: (l.isotopologue, l.v_upper, l.v_lower, l.branch, l.J_lower))
    with open(DEST / "lines.csv", "w") as out:
        out.write("line,FREQofPOT_cm-1,FREQofNIR_cm-1,FREQofDun_cm-1,H_cm-1,H_intensity,H_T_K,L_intensity,L_T_K,"
                  + ",".join(f"{k}_{units.get(k, '')}" for k in PARAMS) + "\n")
        for line in order:
            r = rows[line]
            cm = lambda k: None if r.get(k) is None else r[k] / MHZ_PER_CM          # noqa: E731
            out.write(",".join([str(line), f(cm("pot"), "%.4f"), f(cm("nir"), "%.4f"), f(cm("dun"), "%.4f"),
                                f(cm("H"), "%.6f"), f(r.get("H_intensity"), "%.3g"), f(r.get("H_T"), "%.1f"),
                                f(r.get("L_intensity"), "%.3g"), f(r.get("L_T"), "%.1f")]
                               + [f(r.get(k), "%.6g") for k in PARAMS]) + "\n")
    with gzip.open(DEST / "hyperfine_components.csv.gz", "wt", compresslevel=9) as out:
        out.write("line,component,F_minus_J_lower,offset_MHz,intensity,T_K\n")
        for line in order:
            if line not in comps:
                continue
            freq, cs, T = comps[line]
            for label, dF, nu, inten in cs:
                out.write(f"{line},{label},{dF},{nu - freq:.3f},{inten:.3g},{T:.1f}\n")
    neg = [(line.v_upper, r["H_intensity"]) for line, r in rows.items() if (r.get("H_intensity") or 0) < 0]
    n_neg = len(neg)
    neg_range = f"{min(v for v, _ in neg)}–{max(v for v, _ in neg)}" if neg else "–"
    neg_max = max((abs(x) for _, x in neg), default=0.0)
    print(f"{n_neg} rotational lines with a negative intensity, v' {neg_range}")
    (DEST / "README.md").write_text(f"""# IodineSpec5 predictions

Output of IodineSpec5 (H. Knöckel, B. Bodermann and E. Tiemann, Leibniz Universität Hannover; closed-source
software), run locally for comparison with i2spec in September 2026. The original `*_L.OUT` and `*_H.OUT` files
({len(list(src.glob('*_L.OUT')))} and {len(list(src.glob('*_H.OUT')))}) stay outside this repository; these
tables hold every position and intensity in them. Not for redistribution: the program's licence covers its use,
not the publication of its output.

## `lines.csv`

One row per rotational line, {len(rows)} in all, isotopologue in the line label.

- `FREQofPOT_cm-1`, `FREQofNIR_cm-1`, `FREQofDun_cm-1`: the L file's three positions (vacuum cm⁻¹): from the
  potentials, from the near-infrared correction, and from the Dunham expansion. Blank where it printed none.
- `H_cm-1`: the H file's position of the rotational line (6 decimals; the hyperfine-free centre that its
  components are measured from). `H_intensity` at `H_T_K`.
- `L_intensity` at `L_T_K`: the L file's intensity (two L requests were run at 1000 K).
- The hyperfine parameters the H file printed for the line: eqQ, C (spin-rotation, `Csr`), δ (scalar spin-spin,
  `del`) and d (tensor spin-spin) of X and of B, in the units in the column name.

## `hyperfine_components.csv.gz`

{sum(len(c[1]) for c in comps.values())} components: line, IodineSpec's component label, F″ − J″, the offset
from `H_cm-1` in MHz (the source has 10⁻⁶ cm⁻¹, 0.03 MHz, so three decimals keep every digit), the intensity, and
its temperature.

## Intensities as printed

- IodineSpec5 prints NaN for the intensities of {sum(1 for c in comps.values() for x in c[1] if x[3] != x[3])}
  ¹²⁷I¹²⁹I components.
- It prints a negative intensity for every line it lists with B v′ ≥ 44, in all three isotopologues (v′ =
  {neg_range} here: {n_neg} rotational lines, {sum(1 for c in comps.values() for x in c[1] if x[3] < 0)} components),
  and for none with lower v′.
  Their magnitudes run up to {neg_max:.1e}, as large as the strongest positive lines, so the sign is not
  rounding. Presumably these lines lie outside the range its intensity model was built for; that has not
  been checked. The positions of these lines are unaffected.

Made by `prototypes/condense_iodinespec5.py` in the i2spec repository.
""")
    print(f"-> {DEST}")


if __name__ == "__main__":
    main(sys.argv[1:])
