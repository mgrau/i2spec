"""Does IodineSpec5 use the S06 hyperfine formulae below v' = 44 (where we use BKT02)? And which fits the data?"""
import sys
from collections import defaultdict
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).parent))
from iodinespec5_compare import DEFAULT_DIR, parse_H
from i2spec import hfs_params
from i2spec.constants import MHZ_PER_CM
from i2spec.model import RovibronicModel
from i2spec.observations import load_all, component_rank

pub = RovibronicModel("127I2", "hannover2008")
seen = {}
for p in sorted(DEFAULT_DIR.glob("*_H.OUT")):
    for b in parse_H(p):
        if b["line"] is not None and b["line"].v_upper <= 43:
            seen[b["line"]] = b
res = defaultdict(list)
for line, b in seen.items():
    Eb, Ex = pub._reference_term("B", line.v_upper), pub._reference_term("X", line.v_lower)
    Jl = line.J_lower; Ju = Jl + (1 if line.branch == "R" else -1)
    x = hfs_params.x_state_s06(line.v_lower, Jl, Ex); bb = hfs_params.b_state_s06(line.v_upper, Ju, Eb)
    for name, ours, key in (("eqQ_X", x.eqQ, "eqQ_X"), ("C_X", x.C, "Csr_X"), ("eqQ_B", bb.eqQ, "eqQ_B"), ("C_B", bb.C, "Csr_B"),
                            ("d_B", bb.d, "d_B"), ("delta_B", bb.delta, "del_B")):
        res[name].append(b["params"][key] - ours)
print("IodineSpec5 printed parameters minus S06 formulae evaluated by us, v' <= 43 (eqQ MHz, others kHz):")
for name in res:
    d = np.array(res[name]); print(f"  {name:8s} n={len(d)} mean {d.mean():+9.4f} rms {np.sqrt((d**2).mean()):9.4f} max |{np.abs(d).max():9.4f}|")

# --- the data: BKT02 (ours) vs S06 everywhere, bare formulae, hannover2008 potentials ---
def positions(model, line, table, s06):
    old = hfs_params.E_B_MAX
    if s06:
        hfs_params.E_B_MAX = -1.0
    try:
        nu0, comps = model.hyperfine_components(line.v_upper, line.v_lower, line.J_lower, line.branch, table=table)
    finally:
        hfs_params.E_B_MAX = old
    return nu0, {component_rank(c.label): c.offset for c in comps if c.label}

models, cache = {}, {}
def pos(tag, line, comp):
    if line.isotopologue not in models:
        models[line.isotopologue] = RovibronicModel(line.isotopologue, "hannover2008")
    if (tag, line) not in cache:
        cache[(tag, line)] = positions(models[line.isotopologue], line, None, s06=(tag == "S06"))
    nu0, offs = cache[(tag, line)]
    return nu0 if comp is None else nu0 + offs[component_rank(comp)]

print("\nobserved - model, rms MHz, 127I2 rows with v' <= 43 only; 'hfs' = intra-line intervals, 'pos' = the rest")
print(f"  {'set':16s} {'hfs rows':>8s} {'BKT02':>8s} {'S06':>8s}   {'pos rows':>8s} {'BKT02':>8s} {'S06':>8s}")
tot = defaultdict(list)
for ds in load_all():
    r = defaultdict(list)
    scale = MHZ_PER_CM if ds.unit == "cm-1" else 1.0
    for o in ds.observations:
        ref = o.ref_line or o.line
        if o.line.isotopologue != "127I2" or ref.isotopologue != "127I2" or o.line.v_upper > 43 or ref.v_upper > 43:
            continue
        kind = "hfs" if (o.kind == "interval" and ref == o.line) else "pos"
        for tag in ("BKT02", "S06"):
            m = pos(tag, o.line, o.component) - (pos(tag, ref, o.ref_component) if o.kind == "interval" else 0)
            r[(kind, tag)].append(o.value * scale - m)
            tot[(kind, tag)].append(o.value * scale - m)
    rms = lambda v: np.sqrt(np.mean(np.square(v))) if v else float("nan")  # noqa: E731
    if r:
        print(f"  {ds.id:16s} {len(r[('hfs','BKT02')]):8d} {rms(r[('hfs','BKT02')]):8.3f} {rms(r[('hfs','S06')]):8.3f}   "
              f"{len(r[('pos','BKT02')]):8d} {rms(r[('pos','BKT02')]):8.3f} {rms(r[('pos','S06')]):8.3f}")
rms = lambda v: np.sqrt(np.mean(np.square(v)))  # noqa: E731
print(f"  {'all':16s} {len(tot[('hfs','BKT02')]):8d} {rms(tot[('hfs','BKT02')]):8.3f} {rms(tot[('hfs','S06')]):8.3f}   "
      f"{len(tot[('pos','BKT02')]):8d} {rms(tot[('pos','BKT02')]):8.3f} {rms(tot[('pos','S06')]):8.3f}")
