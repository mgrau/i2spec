"""Every measured B-state hyperfine parameter we hold, by line: the input to i2spec.hfs_table.

Usage:  uv run python prototypes/hfs_measured_table.py          (about 3 minutes)
Output: data/hyperfine_parameters/b_state_lines.json

Two kinds of source. Chen 2004 publishes fitted parameters directly (data/hyperfine_parameters/
chen2004a). For every other precise set the parameters come from i2spec's own four-parameter fit to
the line's intra-line splittings (HyperfineFit.per_line_fit), with the X state held at the published
formulae, so the two kinds are on the same footing. Each row also carries the published formula's value
at that (v', J'), so the table can be used as a correction to the formulae.
"""

import json
from pathlib import Path

from i2spec import hfs_params as hp
from i2spec.hyperfine_fit import HyperfineFit
from i2spec.model import RovibronicModel
from i2spec.observations import load_dataset, load_hyperfine_parameters

ROOT = Path(__file__).resolve().parents[1]
PARAMS = ("eqQ", "C", "d", "delta")
SETS = ("bipm2003a", "bipm2003b", "bipm2003c", "bipm2003d", "bipm2003e", "bipm2005a", "bipm2012a", "reinhardt2006a",
        "reinhardt2007a", "bodermann1998b", "yoshiki2023a", "matsunaga2024a", "kobayashi2016a", "nishiyama2024a")


def main():
    model = RovibronicModel("127I2", grids={"B": dict(rmin=2.35, rmax=12.0, h=0.01, order=10, nlev=80)})

    def formula(vu, vl, Jp):
        p = hp.line_states("127I2", vu, vl, model._reference_term("B", vu), model._reference_term("X", vl))[1](Jp)
        return dict(eqQ=p.eqQ, C=p.C, d=p.d, delta=p.delta)

    rows = []
    for r in load_hyperfine_parameters(ROOT / "data/hyperfine_parameters/chen2004a").rows:
        L = r.line
        f = formula(L.v_upper, L.v_lower, r.J_upper)
        measured = dict(eqQ=r.eqQ, C=r.C * 1e3, d=r.d * 1e3, delta=r.delta * 1e3)
        unc = dict(eqQ=r.eqQ_unc, C=r.C_unc * 1e3, d=r.d_unc * 1e3, delta=r.delta_unc * 1e3)
        rows.append(dict(source="chen2004a", line=f"{L.branch}({L.J_lower}) {L.v_upper}-{L.v_lower}", v=L.v_upper,
                         v_lower=L.v_lower, J=r.J_upper, measured=measured, uncertainty=unc, formula=f,
                         fit_sd_kHz=r.fit_sd_kHz))
    for name in SETS:
        h = HyperfineFit(datasets=[load_dataset(ROOT / "data/observations" / name)], max_v_upper=70)
        for key in h.lines:
            iso, br, Jl, vu, vl = key
            n = sum(h.line_of(d) == key for d in h.data)
            if iso != "127I2" or n < 8:
                continue
            fit = h.per_line_fit(key)
            if fit["median_uncertainty"] > 25:
                continue
            Jp = Jl + (1 if br == "R" else -1)
            f = formula(vu, vl, Jp)
            names = dict(eqQ="eqQ_b", C="C_b", d="d_b", delta="delta_b")
            rows.append(dict(source=name, line=f"{br}({Jl}) {vu}-{vl}", v=vu, v_lower=vl, J=Jp,
                             measured={k: f[k] + fit["x"][names[k]] for k in PARAMS},
                             uncertainty={k: fit["error"][names[k]] for k in PARAMS}, formula=f,
                             fit_sd_kHz=fit["rms_after"], n_splittings=n, normalised=fit["normalised_after"]))
            print(f"  {name:<15} {br}({Jl}) {vu}-{vl}: C_B {rows[-1]['measured']['C']:.3f} kHz ({fit['normalised_after']:.1f} sigma)", flush=True)
    out = dict(
        description="Measured B-state hyperfine parameters of 127I2 lines, one row per line. eqQ in MHz; C, d, delta in kHz. "
                    "'formula' is the published BKT02/S06 value at the same (v', J'), so measured - formula is the "
                    "correction the formulae need there.",
        generated_by="prototypes/hfs_measured_table.py", x_state="held at the published formulae in every fit",
        rows=sorted(rows, key=lambda r: (r["v"], r["J"])))
    (ROOT / "data/hyperfine_parameters/b_state_lines.json").write_text(json.dumps(out, indent=1))
    print(f"{len(rows)} lines written")


if __name__ == "__main__":
    main()
