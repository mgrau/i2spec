"""Compare IodineSpec5 output files with i2spec.

IodineSpec5 (Knöckel & Tiemann, Hannover; closed source) writes two kinds of line list: ``*_L.OUT``
(rotational lines: FREQofPOT / FREQofNIR / FREQofDun, 3-4 decimals) and ``*_H.OUT`` (one block per
rotational line at 6 decimals, with the hyperfine parameters used and every component). This reads
whatever is in an IodineSpec5 ``output/`` folder and compares it with a parameter set of ours.

usage: iodinespec5_compare.py [output_dir ...] [--parameters=hannover2008]
       default directories: $IODINESPEC5_OUTPUT if set, and data/external/iodinespec5/output (not in the
       repository: IodineSpec5's output is not redistributed). Later files win on a repeat.
"""

from __future__ import annotations

import re
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from i2spec.constants import MHZ_PER_CM

from i2spec.model import RovibronicModel
from i2spec.observations import Line, component_rank

import os
DEFAULT_DIR = Path(os.environ.get("IODINESPEC5_OUTPUT", Path(__file__).resolve().parents[1] / "data/external/iodinespec5/output"))
BATCH_DIR = Path(__file__).resolve().parents[1] / "data/external/iodinespec5/output"


def files(dirs, suffix):
    return [p for d in dirs for p in sorted(Path(d).glob(f"*{suffix}"))]
ISO = {"127I2": "127I2", "129I2": "129I2", "127I129I": "127I129I"}


def _unit(header):
    return 1e3 if "GHz" in header else MHZ_PER_CM


def _line(iso, Jp, vp, Js, vs):
    if Jp == Js + 1:
        branch = "R"
    elif Jp == Js - 1:
        branch = "P"
    else:
        return None
    return Line(iso, branch, Js, vp, vs)


def parse_L(path):
    """Rows (line, pot, nir, dun) in MHz, None where IodineSpec5 printed nothing."""
    text = path.read_text(errors="replace").splitlines()
    iso = ISO[text[0].split()[1]]
    header = next(l for l in text if "FREQofPOT" in l)
    scale = _unit(text[text.index(header) + 1])
    cols = {name: header.index(name) + len(name) for name in ("FREQofPOT", "FREQofNIR", "FREQofDun")}
    rows = []
    for l in text[text.index(header) + 3:]:
        m = re.match(r"\s*(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s", l)
        if not m:
            continue
        line = _line(iso, *map(int, m.groups()))
        if line is None:
            continue
        values = {}
        for tok in re.finditer(r"\S+", l[m.end():]):
            end = m.end() + tok.end()
            for name, edge in cols.items():
                if abs(end - edge) <= 2 and name not in values:
                    values[name] = float(tok.group()) * scale
        rows.append((line, values.get("FREQofPOT"), values.get("FREQofNIR"), values.get("FREQofDun")))
    return rows


def parse_H(path):
    """Blocks: dict(line, freq [MHz], params {eqQ_X, C_X, delta_X, d_X, eqQ_B, ...}, comps {label: MHz})."""
    text = path.read_text(errors="replace").splitlines()
    iso = ISO[text[0].split()[1]]
    scale = _unit(next(l for l in text if "@T=" in l))
    blocks, i = [], 0
    while i < len(text):
        if text[i].strip() == "rotational line:":
            Jp, vp, Js, vs, f = text[i + 1].split()[:5]
            block = dict(line=_line(iso, int(Jp), int(vp), int(Js), int(vs)), freq=float(f) * scale,
                         params={}, comps={})
            for l in text[i + 3: i + 5]:
                for m in re.finditer(r"(eqQ|Csr|del|d)\s*\(([XB])\)\s*=\s*(-?[\d.]+)\s*(MHz|kHz)", l):
                    block["params"][f"{m[1]}_{m[2]}"] = float(m[3])
            j = i + 7
            while j < len(text) and re.match(r"\s+[a-z]\s*\d+\s+-?\d+\s+", text[j]):
                m = re.match(r"\s+([a-z])\s*(\d+)\s+(-?\d+)\s+(-?[\d.]+)", text[j])
                block["comps"][f"{m[1]}{m[2]}"] = float(m[4]) * scale
                j += 1
            blocks.append(block)
            i = j
        else:
            i += 1
    return blocks


def main(argv):
    dirs = [d for d in (DEFAULT_DIR, BATCH_DIR) if d.exists()]
    parameters = "hannover2008"
    for a in argv:
        if a.startswith("--parameters="):
            parameters = a.split("=", 1)[1]
        else:
            dirs = [Path(a)] if dirs in ([DEFAULT_DIR, BATCH_DIR], [DEFAULT_DIR], [BATCH_DIR]) else dirs + [Path(a)]
    models = {}

    def model(iso):
        if iso not in models:
            models[iso] = RovibronicModel(iso, parameters)
        return models[iso]

    # --- rotational lines, from the H files (6 decimals) and the L files (3-4 decimals) ---
    import random
    from i2spec.observations import load_all
    seen, nir_lines = {}, set()
    for p in files(dirs, "_H.OUT"):
        for b in parse_H(p):
            if b["line"] is not None and b["line"].v_upper <= 59 and b["line"].v_lower <= 54:
                seen[b["line"]] = b
    for p in files(dirs, "_L.OUT"):
        for line, pot, nir, dun in parse_L(p):
            if nir is not None:
                nir_lines.add(line)          # IodineSpec5 replaced the potential value by its NIR local Dunham
    obs_lines = set()
    for ds in load_all():
        for o in ds.observations:
            obs_lines.update(l for l in (o.line, o.ref_line) if l is not None)
    sample = set(random.Random(0).sample(sorted(seen, key=str), min(2000, len(seen)))) | obs_lines
    print(f"{len(seen)} distinct rotational lines in {len(files(dirs, '_H.OUT'))} H files; ours: {parameters}; "
          f"hyperfine compared on {len(sample & set(seen))} of them")
    groups = defaultdict(list)
    comp_res = defaultdict(list)
    param_res = defaultdict(list)
    worst = []
    from i2spec import hfs_params
    for line, b in seen.items():
        m = model(line.isotopologue)
        nu0 = m.transition(line.v_upper, line.v_lower, line.J_lower, line.branch) * MHZ_PER_CM
        d = b["freq"] - nu0
        key = line.isotopologue + (" v'<=43" if line.v_upper <= 43 else " v'>43") + (" v''<=17" if line.v_lower <= 17 else " v''>17") \
            + (" NIR-replaced" if line in nir_lines else "") + (" J''>200" if line.J_lower > 200 else "")
        groups[key].append((line, d))
        if line.v_upper <= 43 and line.v_lower <= 17 and line not in nir_lines and line.isotopologue == "127I2":
            worst.append((abs(d), str(line), d, b["freq"]))
        if line not in sample:
            continue
        try:
            _, comps = m.hyperfine_components(line.v_upper, line.v_lower, line.J_lower, line.branch, table=None)
        except Exception as e:
            print("skip", line, e)
            continue
        ours = {component_rank(c.label): c.offset for c in comps if c.label}
        for lab, f in b["comps"].items():
            r = component_rank(lab)
            if r in ours:
                comp_res[key].append((line, (f - b["freq"]) - ours[r], lab))
        xp, bp = hfs_params.line_states(line.isotopologue, line.v_upper, line.v_lower,
                                        m._reference_term("B", line.v_upper), m._reference_term("X", line.v_lower))
        Jl, Ju = line.J_lower, line.J_lower + (1 if line.branch == "R" else -1)
        x, bb = xp(Jl), bp(Ju)
        for name, ours_v, theirs_key in (("eqQ_X", x.eqQ, "eqQ_X"), ("C_X", x.C, "Csr_X"), ("delta_X", x.delta, "del_X"),
                                         ("d_X", x.d, "d_X"), ("eqQ_B", bb.eqQ, "eqQ_B"), ("C_B", bb.C, "Csr_B"),
                                         ("delta_B", bb.delta, "del_B"), ("d_B", bb.d, "d_B")):
            if theirs_key in b["params"]:
                param_res[(key, name)].append(b["params"][theirs_key] - ours_v)
    print("\nrotational line frequency, IodineSpec5 (H file, 6 decimals) - ours, MHz:")
    for key in sorted(groups):
        d = np.array([x[1] for x in groups[key]])
        print(f"  {key:34s} n={len(d):6d}  mean {d.mean():+9.4f}  rms {np.sqrt((d**2).mean()):9.4f}  max |{np.abs(d).max():9.4f}|")
    print("\nhyperfine component offsets (component - rotational line), IodineSpec5 - ours (bare formulae), kHz:")
    for key in sorted(comp_res):
        d = np.array([x[1] for x in comp_res[key]]) * 1e3
        pc = np.percentile(np.abs(d), [50, 90, 99])
        print(f"  {key:34s} n={len(d):6d}  mean {d.mean():+8.2f}  rms {np.sqrt((d**2).mean()):8.2f}  |50/90/99%| {pc[0]:.1f}/{pc[1]:.1f}/{pc[2]:.1f}  max |{np.abs(d).max():8.2f}|")
    allc = sorted(((abs(x[1]), str(x[0]), x[2], x[1]) for k in comp_res for x in comp_res[k]), reverse=True)
    print("  worst components:")
    for a, line, lab, d in allc[:8]:
        print(f"    {line:26s} {lab:4s} {d*1e3:+9.1f} kHz")
    print("\nhyperfine parameters printed by IodineSpec5 - ours (eqQ MHz, others kHz):")
    for (key, name) in sorted(param_res):
        d = np.array(param_res[(key, name)])
        print(f"  {key:34s} {name:8s} n={len(d):5d}  mean {d.mean():+10.4f}  rms {np.sqrt((d**2).mean()):10.4f}  max |{np.abs(d).max():10.4f}|")
    worst.sort(reverse=True)
    print("\nworst 12 rotational lines (127I2, v' <= 43, v'' <= 17, not NIR-replaced):")
    for a, line, d, f in worst[:12]:
        print(f"  {line:28s} IodineSpec5 {f:16.3f} MHz   theirs-ours {d:+10.3f} MHz")

    # --- the observation sets, on the lines IodineSpec5 has computed: theirs vs ours ---
    from i2spec.observations import load_all
    from i2spec.constants import DEFAULT_PARAMETERS
    from mlr_evaluate import load_x, load_b, GRID_B, ROOT
    pair = {"X": load_x(ROOT / "data/potentials/mlr_x_2026c.json"), "B": load_b(ROOT / "data/potentials/mlr_b_2026c.json")}
    ours = {"hannover2008 bare": (RovibronicModel, "hannover2008", None),
            f"{DEFAULT_PARAMETERS}+table": (RovibronicModel, DEFAULT_PARAMETERS, "default"),
            "2026c pair+table": (RovibronicModel, DEFAULT_PARAMETERS, "default")}
    cache = {}

    def our_position(tag, line, comp):
        if (tag, line) not in cache:
            _, pset, table = ours[tag]
            if (tag, line.isotopologue) not in cache:
                mlr = tag.startswith("2026c") and line.isotopologue == "127I2"
                cache[(tag, line.isotopologue)] = RovibronicModel(line.isotopologue, pset, grids=GRID_B if mlr else None,
                                                                  potentials=pair if mlr else None)
            m = cache[(tag, line.isotopologue)]
            nu0, comps = m.hyperfine_components(line.v_upper, line.v_lower, line.J_lower, line.branch, table=table)
            cache[(tag, line)] = (nu0, {component_rank(c.label): c.offset for c in comps if c.label})
        nu0, offs = cache[(tag, line)]
        return nu0 if comp is None else nu0 + offs[component_rank(comp)]

    def their_position(line, comp):
        """IodineSpec5 labels every component a1..a21 by frequency; the data sets use the BIPM letters
        (a, b, ... for the line, the number for the rank), so match by rank."""
        b = seen[line]
        if comp is None:
            return b["freq"]
        by_rank = {component_rank(k): v for k, v in b["comps"].items()}
        return by_rank[component_rank(comp)]

    print("\nobservation sets on lines present in these H files: rms of observed - model, MHz")
    print(f"  {'set':16s} {'rows':>9s} {'IodineSpec5':>12s} " + " ".join(f"{t:>22s}" for t in ours))
    for ds in load_all():
        res = defaultdict(list)
        for o in ds.observations:
            ref = o.ref_line or o.line
            if o.line not in seen or (o.kind == "interval" and ref not in seen):
                continue
            try:
                t = their_position(o.line, o.component) - (their_position(ref, o.ref_component) if o.kind == "interval" else 0)
                res["IodineSpec5"].append(o.value * {"MHz": 1, "kHz": 1e-3, "GHz": 1e3, "cm-1": MHZ_PER_CM}[ds.unit] - t)
                for tag in ours:
                    u = our_position(tag, o.line, o.component) - (our_position(tag, ref, o.ref_component) if o.kind == "interval" else 0)
                    res[tag].append(o.value * {"MHz": 1, "kHz": 1e-3, "GHz": 1e3, "cm-1": MHZ_PER_CM}[ds.unit] - u)
            except KeyError:
                continue
        if res["IodineSpec5"]:
            rms = lambda v: np.sqrt(np.mean(np.square(v)))  # noqa: E731
            print(f"  {ds.id:16s} {len(res['IodineSpec5']):4d}/{len(ds.observations):<4d} {rms(res['IodineSpec5']):12.3f} " + " ".join(f"{rms(res[t]):22.3f}" for t in ours))

    # --- L files: which column IodineSpec5 uses where, and how the three agree with ours ---
    print("\nL files, theirs - ours, MHz (Pot = potential, NIR = local Dunham, Dun = Gerstenkorn/Martin Dunham):")
    lgroups = defaultdict(list)
    for p in files(dirs, "_L.OUT"):
        for line, pot, nir, dun in parse_L(p):
            if line.v_upper > 59 or line.v_lower > 54:   # beyond the solver's 60 levels
                continue
            m = model(line.isotopologue)
            nu0 = m.transition(line.v_upper, line.v_lower, line.J_lower, line.branch) * MHZ_PER_CM
            key = ("v'<=43" if line.v_upper <= 43 else "v'>43") + (" v''<=17" if line.v_lower <= 17 else " v''>17")
            for name, val in (("Pot", pot), ("NIR", nir), ("Dun", dun)):
                if val is not None:
                    lgroups[(key, name)].append(val - nu0)
    for key in sorted(lgroups):
        d = np.array(lgroups[key])
        print(f"  {key[0]:18s} {key[1]}  n={len(d):5d}  mean {d.mean():+10.3f}  rms {np.sqrt((d**2).mean()):10.3f}  max |{np.abs(d).max():10.3f}|")


if __name__ == "__main__":
    main(sys.argv[1:])
