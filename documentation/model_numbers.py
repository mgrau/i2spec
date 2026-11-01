"""The numbers the documentation shows, computed when the site is built from the released model and the files of
its release, so that no page carries a typed value of the model.

documentation/hooks.py replaces each ``{{ key }}`` in a page with ``numbers()[key]``; a key that does not exist stops
the build. Three groups:

  we.*    the worked example, R(56) 32-0 of ¹²⁷I₂ (the 532 nm standard), followed through every step of the model
  fit.*   what the fit adjusts and the data it uses (the release's manifest, model files and training inputs)
  unc.*   the frozen uncertainty model and its tests

``worked_example()``, ``fit_counts()`` and ``uncertainty_numbers()`` return the raw values; tests/test_docs.py checks
them against the package's public interface.
"""
from __future__ import annotations

import csv
import gzip
import json
import math
from collections import Counter
from functools import lru_cache
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
TRAINING = ROOT / "data" / "training"

#: the worked example: R(56) 32-0 of 127I2, at 20 °C
LINE = "R(56) 32-0"
BRANCH, J_LOWER, V_UPPER, V_LOWER = "R", 56, 32, 0
J_UPPER = J_LOWER + 1
TEMPERATURE = 293.15
COMPONENT = "a10"

NNBSP = "\u202f"          # narrow no-break space: the thousands separator of every grouped number
#: the isotopologues as the pages write them (documentation/hooks.py uses it too)
ISO = {"127I2": "¹²⁷I₂", "129I2": "¹²⁹I₂", "127I129I": "¹²⁷I¹²⁹I"}


# ------------------------------------------------------------------------------------------------- formatting

def group(x, decimals=0):
    """x with its integer digits grouped in threes (from five digits on), e.g. 18 788.3355."""
    text = f"{abs(x):.{decimals}f}"
    whole, _, frac = text.partition(".")
    if len(whole) > 4:
        whole = NNBSP.join(whole[max(0, i - 3):i] for i in range(len(whole), 0, -3)[::-1])
    return ("−" if x < 0 and float(text) != 0 else "") + whole + (f".{frac}" if frac else "")


def signed(x, decimals=0):
    return ("+" if x >= 0 else "") + group(x, decimals)


def count(n):
    return group(int(n))


def percent(x):
    return f"{100 * x:.0f} %"


def sci(x, digits=2):
    """x as a × 10ⁿ with Unicode superscripts."""
    mantissa, exponent = f"{x:.{digits}e}".split("e")
    sup = str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹")
    return f"{mantissa} × 10{str(int(exponent)).translate(sup)}"


def _turning_points(potential, energy, r_min, r_max):
    """The classical turning points (Å) of a rotationless level: where V(R) = E, on either side of the well."""
    R = np.linspace(r_min, r_max, 60001)
    V = np.asarray(potential(R), float) - energy
    crossings = np.flatnonzero(np.sign(V[:-1]) != np.sign(V[1:]))
    out = [R[k] - V[k] * (R[k + 1] - R[k]) / (V[k + 1] - V[k]) for k in crossings]
    return out[0], out[-1]


# ------------------------------------------------------------------------------------------------- the example

def _strength(nu):
    """(S at TEMPERATURE in cm/molecule, S0, Q(T), E'' cm-1) of the worked line, as the package's line list has it.

    Catalog().line gives the same from the cached line list of every line (minutes to build the first time); this
    builds only the X levels of the partition function (with the line list's convergence rule) and the line's own
    two levels (about 450 X eigensolves, some 20 s), so the site builds quickly from an empty cache.
    tests/test_docs.py compares it with Catalog().line."""
    from i2spec.intensity import C2, _Build, intensity_model, mu_tellinghuisen2011, nuclear_spin_weight

    model = intensity_model()
    temps = np.unique([200.0, math.sqrt(200.0 * 600.0), 600.0])     # the line list's temperature range
    build = _Build(model.isotopologue, model.parameters, model.grids, model.solver, nu - 0.01, nu + 0.01, temps,
                   0.0, mu_tellinghuisen2011)
    levels, matched, Q, J = [], [], np.zeros_like(temps), 0
    while True:
        e, ok = build.x_levels(J)
        term = nuclear_spin_weight(J, model.isotopologue) * (2 * J + 1) * np.exp(-C2 * e[None, :] / temps[:, None]).sum(axis=1)
        levels.append(e)
        matched.append(ok)
        Q += term
        if J > 10 and term[-1] < 1e-9 * Q[-1]:
            break
        J += 1
    cols, _, _ = build.lines([J_LOWER], {J_LOWER: levels[J_LOWER]}, {J_LOWER: matched[J_LOWER]}, Q)
    k = int(np.flatnonzero((cols["v_upper"] == V_UPPER) & (cols["v_lower"] == V_LOWER) & (cols["branch"] == 1))[0])
    x = np.array(levels)
    weight = np.array([nuclear_spin_weight(j, model.isotopologue) for j in range(len(x))]) * (2 * np.arange(len(x)) + 1)
    Q_T = float((weight[:, None] * np.exp(-C2 * x / TEMPERATURE)).sum())
    S0, E = float(cols["strength0"][k]), float(cols["E_lower"][k])
    return S0 * math.exp(-C2 * E / TEMPERATURE) / Q_T, S0, Q_T, E


def _measured():
    """The most precise measured frequency of the worked component among the fit's data sets: (value MHz, sigma MHz,
    set id, short name)."""
    from i2spec.fit.data import fit_datasets
    from i2spec.observations import UNITS
    from i2spec.webapp import _short

    best = None
    for ds in fit_datasets():
        scale = UNITS[ds.unit]
        for o in ds.observations:
            line = o.line
            if (o.kind == "frequency" and o.component == COMPONENT and line.isotopologue == "127I2"
                    and (line.branch, line.J_lower, line.v_upper, line.v_lower) == (BRANCH, J_LOWER, V_UPPER, V_LOWER)):
                u = o.uncertainty * scale
                if best is None or u < best[1]:
                    best = (o.value * scale, u, ds.id, _short(ds.id))
    return best


def _training_rows():
    """How the fit's training table reaches the worked line and its two levels."""
    from i2spec.fit import rows as R

    rows, extra = R.read(TRAINING / R.RELEASE / "rows.csv.gz")
    c, role = rows.cols, extra["continuation"]
    train = np.char.startswith(role, "train:")
    precision = np.isin(c["category"], ["comb", "doppler"])
    this = ((c["branch"] == BRANCH) & (c["J_lower"] == J_LOWER) & (c["v_upper"] == V_UPPER)
            & (c["v_lower"] == V_LOWER) & (c["kind"] == "abs"))
    on = {state: np.array([any(lv[1] == state and lv[2] == v for (lv, _) in levels) for levels in rows.levels])
          for state, v in (("B", V_UPPER), ("X", V_LOWER))}
    with gzip.open(TRAINING / R.RELEASE / "hyperfine.csv.gz", "rt", newline="", encoding="utf-8") as f:
        hf = [r for r in csv.DictReader(f) if r["kind"] == "interval" and r["line"] == f"127I2 {LINE}"]
    return dict(
        train_atlas=int((this & train & ~precision).sum()),
        held_out_precision=int((this & ~train & precision).sum()),
        held_out_sources=sorted(set(c["source"][this & ~train & precision])),
        hf_held_out=sum(r["future"] == "1" for r in hf),
        b_precision=int((on["B"] & train & precision).sum()), b_all=int((on["B"] & train).sum()),
        b_precision_held_out=int((on["B"] & ~train & precision).sum()),
        x_precision=int((on["X"] & train & precision).sum()), x_all=int((on["X"] & train).sum()))


@lru_cache(maxsize=None)
def worked_example():
    """Every number of the worked example, raw (cm-1, MHz, Å, cm/molecule)."""
    from i2spec import MHZ_PER_CM, RovibronicModel
    from i2spec.constants import ATOMIC_MASS
    from i2spec.intensity import C2, nuclear_spin_weight
    from i2spec.lookup import from_wavenumber
    from i2spec.preferred import model
    from i2spec.spectrum import doppler_fwhm

    m = model("127I2")                      # the released model, ResearchModel("i2spec")
    levels = RovibronicModel("127I2")       # its levels and potentials (the released model is the default)
    out = dict(model_id=m.info.id)

    # potentials (core-fit constants) and the two levels
    values = m.coefficients["base_coefficients"]["values"]
    out.update(x_re=math.exp(values["X"]["logRe"]), x_de=math.exp(values["X"]["logDe"]),
               b_re=math.exp(values["B"]["logRe"]), b_te=values["B"]["Te"], b_de=math.exp(values["B"]["logDe"]),
               iso_a0=values["B"]["iso0"], iso_a1=values["B"]["iso1"])
    ex0, ex = levels.energy("X", V_LOWER, 0), levels.energy("X", V_LOWER, J_LOWER)
    eb0, eb = levels.energy("B", V_UPPER, 0), levels.energy("B", V_UPPER, J_UPPER)
    out.update(ex0=ex0, ex=ex, eb0=eb0, eb=eb, ex_rot=ex - ex0, eb_rot=eb - eb0,
               x_bv=(levels.energy("X", V_LOWER, 2) - ex0) / 6, b_bv=(levels.energy("B", V_UPPER, 2) - eb0) / 6,
               e_lower=ex - levels.energy("X", 0, 0))
    out["x_turn"] = _turning_points(levels.potentials["X"], ex0, 2.2, 3.4)
    out["b_turn"] = _turning_points(levels.potentials["B"], eb0, 2.5, 6.0)

    # line position
    nu = m.transition(V_UPPER, V_LOWER, J_LOWER, BRANCH)
    out.update(nu=nu, nu_MHz=nu * MHZ_PER_CM, nm=1e7 / nu, nm_air=from_wavenumber(nu, "nm-air"),
               origin=eb0 - ex0, below_origin=(eb0 - ex0) - nu)

    # hyperfine
    hb, hx = m.level_hyperfine("B", V_UPPER, J_UPPER), m.level_hyperfine("X", V_LOWER, J_LOWER)
    out.update(b_eqq=hb.eqQ, b_c=hb.C, b_d=hb.d, b_delta=hb.delta, x_eqq=hx.eqQ, x_c=hx.C, x_d=hx.d, x_delta=hx.delta)

    def spread(state, v, J):                # the width of a level's hyperfine structure, MHz
        energies = [h.energy for h in m.hyperfine_levels(state, v, J)]
        return max(energies) - min(energies)
    out.update(x_spread=spread("X", V_LOWER, J_LOWER), b_spread=spread("B", V_UPPER, J_UPPER))
    _, components = m.hyperfine_components(V_UPPER, V_LOWER, J_LOWER, BRANCH)
    main = {c.label: c for c in components if c.label}
    out.update(n_main=len(main), main_strength=sum(c.strength for c in main.values()),
               offsets={k: c.offset for k, c in main.items()}, comp_strength=main[COMPONENT].strength,
               comp_MHz=m.frequency(LINE, COMPONENT), span=main["a15"].offset - main["a1"].offset)

    # strength and Doppler width at 20 °C
    S, S0, Q, E = _strength(nu)
    g = nuclear_spin_weight(J_LOWER, "127I2")
    out.update(S=S, S0=S0, Q=Q, E_strength=E, g_ns=g, honl_london=J_LOWER + 1 if BRANCH == "R" else J_LOWER,
               population=g * (2 * J_LOWER + 1) * math.exp(-C2 * E / TEMPERATURE) / Q,
               doppler_MHz=float(doppler_fwhm(nu, TEMPERATURE, 2 * ATOMIC_MASS[127])) * MHZ_PER_CM)

    # uncertainty
    center_u, comp_u = m.uncertainty(LINE), m.uncertainty(LINE, COMPONENT)
    interval_u = m.uncertainty(LINE, "a15", "a1")
    a, kappa = m.hyperfine_terms(V_UPPER, J_UPPER, V_LOWER, J_LOWER)
    out.update(sigma=center_u.sigma_MHz, sigma_stat=center_u.statistical_MHz, sigma_model=center_u.model_MHz,
               sigma_level=center_u.model_MHz / math.sqrt(2), sigma_comp=comp_u.sigma_MHz,
               sigma_hf=comp_u.hyperfine_MHz, sigma_interval=interval_u.sigma_MHz, hf_a=a, hf_kappa=kappa,
               status=center_u.status)

    # the measurement and the fit
    value, u, set_id, short = _measured()
    out.update(measured=value, measured_sigma=u, measured_set=set_id, measured_short=short,
               residual=out["comp_MHz"] - value, z=(out["comp_MHz"] - value) / math.hypot(out["sigma_comp"], u))
    out.update(_training_rows())
    return out


def _we_strings(w):
    o = w["offsets"]
    return {
        "we.T": f"{TEMPERATURE:.2f}",
        # potentials and levels
        "we.x_re": f"{w['x_re']:.5f}", "we.x_de": group(w["x_de"], 2), "we.b_re": f"{w['b_re']:.5f}",
        "we.b_te": group(w["b_te"], 3), "we.b_de": group(w["b_de"], 2),
        "we.ex0": group(w["ex0"], 4), "we.ex": group(w["ex"], 4), "we.eb0": group(w["eb0"], 4),
        "we.eb": group(w["eb"], 4), "we.ex_rot": group(w["ex_rot"], 4), "we.eb_rot": group(w["eb_rot"], 4),
        "we.x_bv": f"{w['x_bv']:.6f}", "we.b_bv": f"{w['b_bv']:.6f}",
        "we.x_turn": f"{w['x_turn'][0]:.3f}–{w['x_turn'][1]:.3f}", "we.b_turn": f"{w['b_turn'][0]:.3f}–{w['b_turn'][1]:.3f}",
        "we.e_lower": group(w["e_lower"], 4), "we.ex_raw": f"{w['ex']:.4f}", "we.eb_raw": f"{w['eb']:.4f}",
        "we.iso_a0": group(w["iso_a0"], 4), "we.iso_a1": group(w["iso_a1"], 3),
        # line position
        "we.nu": group(w["nu"], 5), "we.nu_raw": f"{w['nu']:.5f}", "we.nu_MHz": group(w["nu_MHz"], 1), "we.nu_MHz_raw": f"{w['nu_MHz']:.1f}",
        "we.nm": f"{w['nm']:.6f}", "we.nm_air": f"{w['nm_air']:.6f}", "we.origin": group(w["origin"], 4),
        "we.below_origin": f"{w['below_origin']:.2f}",
        # hyperfine
        "we.b_eqq": group(w["b_eqq"], 2), "we.b_eqq_raw": f"{w['b_eqq']:.2f}", "we.b_c": f"{w['b_c']:.2f}", "we.b_d": group(w["b_d"], 2),
        "we.b_delta": group(w["b_delta"], 2), "we.x_eqq": group(w["x_eqq"], 2), "we.x_c": f"{w['x_c']:.3f}",
        "we.x_d": f"{w['x_d']:.3f}", "we.x_delta": f"{w['x_delta']:.3f}",
        "we.x_spread": group(w["x_spread"], 0), "we.b_spread": group(w["b_spread"], 0),
        "we.eqq_ratio": f"{w['x_eqq'] / w['b_eqq']:.1f}",
        "we.n_main": str(w["n_main"]),
        "we.main_strength": f"{100 * w['main_strength']:.2f} %",
        "we.a1": group(o["a1"], 3), "we.a15": signed(o["a15"], 3), "we.a10": signed(o[COMPONENT], 3),
        "we.span": group(w["span"], 3),
        "we.comp_MHz": group(w["comp_MHz"], 3), "we.comp_MHz_raw": f"{w['comp_MHz']:.3f}",
        "we.comp_strength": f"{w['comp_strength']:.4f}",
        "we.offsets_table": "\n".join(
            f"| {k} | {signed(o[k], 3)} |" for k in sorted(o, key=lambda k: int(k[1:]))),
        # strength
        "we.S": sci(w["S"], 3), "we.S_raw": f"{w['S']:.3e}", "we.E_strength_raw": f"{w['E_strength']:.4f}", "we.Q": group(w["Q"], 0), "we.g_ns": str(w["g_ns"]),
        "we.honl_london": str(w["honl_london"]), "we.population": sci(w["population"], 3),
        "we.doppler": f"{w['doppler_MHz']:.0f}", "we.doppler_raw": f"{w['doppler_MHz']:.1f}",
        # uncertainty
        "we.sigma": f"{w['sigma']:.3f}", "we.sigma_stat": f"{w['sigma_stat']:.3f}",
        "we.sigma_model": f"{w['sigma_model']:.3f}", "we.sigma_level": f"{w['sigma_level']:.3f}",
        "we.sigma_comp": f"{w['sigma_comp']:.3f}", "we.sigma_raw": f"{w['sigma']:.3f}", "we.sigma_hf": f"{w['sigma_hf']:.3f}",
        "we.sigma_interval": f"{w['sigma_interval']:.3f}", "we.hf_a": f"{1e3 * w['hf_a']:.0f}",
        "we.hf_kappa": sci(w["hf_kappa"], 1), "we.status": w["status"],
        # the measurement and the fit
        "we.measured": group(w["measured"], 4), "we.measured_sigma": f"{1e3 * w['measured_sigma']:.1f}",
        "we.measured_short": w["measured_short"], "we.residual": signed(w["residual"], 3),
        "we.z": signed(w["z"], 2),
        "we.train_atlas": str(w["train_atlas"]), "we.held_out_precision": str(w["held_out_precision"]),
        "we.hf_held_out": str(w["hf_held_out"]),
        "we.b_precision": count(w["b_precision"]), "we.b_all": count(w["b_all"]),
        "we.b_precision_held_out": count(w["b_precision_held_out"]),
        "we.x_precision": count(w["x_precision"]), "we.x_all": count(w["x_all"]),
    }


# ------------------------------------------------------------------------------------------------- the fit

def _release():
    """The released model's manifest entries: (the model's, the continued model's, the core fit's), with the files
    of each."""
    from i2spec.preferred import model
    from i2spec.research import _entry, _read

    entry = _entry(model("127I2").info.id)
    expanded = _entry(entry["extends"]["id"])
    core_id = _read(expanded, "model.json")["base"]["model"]
    return entry, expanded, _entry(core_id), _read


@lru_cache(maxsize=None)
def fit_counts():
    """What each stage of the fit adjusts and how many rows of which data it was fitted to."""
    from i2spec.fit import rows as R
    from i2spec.fit.data import ROLES

    entry, expanded, core, read = _release()
    core_model, ex_model = read(core, "model.json"), read(expanded, "model.json")
    potential = [p for p in core_model["optimizer_parameters"]["potential"] if not p["name"].startswith("iso")]
    shift = [p for p in core_model["optimizer_parameters"]["potential"] if p["name"].startswith("iso")]
    per_state = {s: Counter("c" if p["name"].startswith("c") else "alpha" if p["name"].startswith("alpha") else "other"
                            for p in potential if p["state"] == s) for s in ("X", "B")}
    hf_core = Counter(p["name"] for p in core_model["optimizer_parameters"]["hyperfine"])
    detail = ex_model["inventory_detail"]
    tails = detail["potential_tails"]["names"]

    rows, extra = R.read(TRAINING / R.RELEASE / "rows.csv.gz")
    c = rows.cols
    core_train = np.char.startswith(extra["core"], "train:")
    cont_train = np.char.startswith(extra["continuation"], "train:")
    benchmark = np.char.startswith(extra["core"], "benchmark:")
    by_cat = Counter(c["category"][core_train])
    sources = {cat: len(set(c["source"][core_train & (c["category"] == cat)])) for cat in by_cat}
    with gzip.open(TRAINING / R.RELEASE / "hyperfine.csv.gz", "rt", newline="", encoding="utf-8") as f:
        hf = list(csv.DictReader(f))
    intervals = [r for r in hf if r["kind"] == "interval"]
    constants = Counter(r["source"] for r in hf if r["kind"] == "constant")
    with gzip.open(TRAINING / R.RELEASE / "cross_validation.json.gz", "rt", encoding="utf-8") as f:
        cv = json.load(f)
    future = Counter(r["kind"] for r in cv["core_folds"] if r["fold"] == "future")
    core_folds = sorted({r["fold"] for r in cv["core_folds"]})
    centers = read(entry, "isotopologues.json")["isotope_shift"]["data"]      # the CIPM centers of the core fit

    roles = json.loads((TRAINING / ROLES / "roles.json").read_text())["sets"]
    outcome = json.loads((TRAINING / ROLES / "outcome.json").read_text())["final_roles"]
    final = Counter(outcome.get(s["set"], s["role"]) for s in roles)

    return dict(
        inventory=entry["inventory"], core_inventory=core["inventory"],
        n_potential=len(potential), n_shift=len(shift), n_hf_core=sum(hf_core.values()), hf_core=dict(hf_core),
        x_c=per_state["X"]["c"], x_alpha=per_state["X"]["alpha"], b_c=per_state["B"]["c"], b_alpha=per_state["B"]["alpha"],
        n_tails=detail["potential_tails"]["fitted"],
        tails=dict(x_outer=sum(t.startswith("X:to") for t in tails), x_inner=sum(t.startswith("X:ti") for t in tails),
                   x_rot=sum(t.startswith("X:ro") for t in tails), x_splines=sum(t.startswith("X:w") for t in tails),
                   b_outer=sum(t.startswith("B:to") for t in tails), b_rot=sum(t.startswith("B:ro") for t in tails)),
        hf_beyond=detail["hyperfine_fitted"]["names"],
        core_rows=int(core_train.sum()), core_by_category=dict(by_cat), core_sources=sources,
        core_kinds=dict(Counter(c["kind"][core_train & np.isin(c["category"], ["comb", "doppler"])])),
        benchmark_rows=int(benchmark.sum()), benchmark_sets=len(R.NEW_SETS),
        benchmark_absolute=future["absolute"], benchmark_intraline=future["intraline"],
        benchmark_absolute_rms=core["absolute_rms_MHz"], benchmark_intraline_rms=core["intraline_rms_MHz"],
        core_folds=len(core_folds), core_level_folds=sum(f.startswith("physical_level") for f in core_folds),
        core_source_folds=sum(f.startswith("physical_source") for f in core_folds),
        continuation_folds=len(cv["continuation_folds"]),
        hyperfine_folds=max(len(v["folds"]) for v in cv["hyperfine_folds"].values()),
        hf_core_trained=sum(r["inside"] == "1" and r["future"] == "0" for r in intervals),
        hf_beyond_rows=sum(r["inside"] == "0" for r in intervals), hf_constants=dict(constants),
        hf_core_sources=len({r["source"] for r in intervals if r["inside"] == "1" and r["future"] == "0"}),
        calibration={k: core_model["calibration"][k] for k in ("atlas_sr", "atlas_apo")},
        orsay4_offset=ex_model["calibration"]["orsay4_offset_MHz"],
        isotopologue_centers=[(d["line"], d["sigma_MHz"]) for d in centers],
        continuation_rows=int(cont_train.sum()), continuation_by_category=dict(Counter(c["category"][cont_train])),
        continuation_regions={k: (v["n"], v["rms_MHz"]) for k, v in ex_model["fit"]["training"].items()},
        roles=dict(final),
        exclusion=entry["exclusion"], orsay1_core_v=(int(c["v_lower"][core_train & (c["category"] == "orsay1")].min()),
                                                     int(c["v_lower"][core_train & (c["category"] == "orsay1")].max())),
        martin_core_v=(int(c["v"][core_train & (c["category"] == "martin")].min()),
                       int(c["v"][core_train & (c["category"] == "martin")].max())),
    )


#: how the training table names its data categories
CATEGORY = {"comb": "comb-referenced and other precision lines", "doppler": "Doppler-free lines", "atlas_sr": "Salami & Ross atlas",
            "atlas_apo": "NIST atlas of the APO cell", "orsay1": "Orsay atlas part I", "orsay4": "Orsay atlas Partie IV",
            "martin": "X level constants of Martin *et al.*", "martin_dissoc": "X levels near the limit (Martin *et al.* 1983)",
            "emission": "comb-referenced emission lines"}


def _fit_strings(f):
    inv, core_inv = f["inventory"], f["core_inventory"]
    t = f["tails"]
    regions = f["continuation_regions"]

    def region_rows(prefix):
        out = []
        for key, (n, rms) in regions.items():
            region, cat = key.split(" | ")
            if region == prefix:
                out.append(f"{CATEGORY[cat]}: {count(n)} rows, {rms:.2g} MHz rms" if rms < 1 else
                           f"{CATEGORY[cat]}: {count(n)} rows, {rms:.0f} MHz rms")
        return "; ".join(out)

    def region_name(region):                  # "X v18-25" -> "X v″ = 18–25", "B v78+" -> "B v′ ≥ 78"
        state, v = region.split(" v")
        prime = "″" if state == "X" else "′"
        return f"{state} v{prime} ≥ {v[:-1]}" if v.endswith("+") else f"{state} v{prime} = {v.replace('-', '–')}"

    region_table = "\n".join(f"| {region_name(region)} | {region_rows(region)} |"
                             for region in dict.fromkeys(k.split(" | ")[0] for k in regions) if region != "core")
    sr, apo = f["calibration"]["atlas_sr"], f["calibration"]["atlas_apo"]
    return {
        "fit.fitted": str(inv["fitted_molecular"]), "fit.inherited": str(inv["inherited_molecular"]),
        "fit.calibration": str(inv["fitted_calibration"]),
        "fit.core_calibration": str(core_inv["fitted_calibration"]),
        "fit.n_potential": str(f["n_potential"]), "fit.n_shift": str(f["n_shift"]), "fit.n_hf_core": str(f["n_hf_core"]),
        "fit.hf_eqq": str(f["hf_core"]["eqQ"]), "fit.hf_c": str(f["hf_core"]["C"]), "fit.hf_d": str(f["hf_core"]["d"]),
        "fit.hf_delta": str(f["hf_core"]["delta"]),
        "fit.x_c": str(f["x_c"]), "fit.x_alpha": str(f["x_alpha"]), "fit.b_c": str(f["b_c"]),
        "fit.b_alpha": str(f["b_alpha"]),
        "fit.n_tails": str(f["n_tails"]), "fit.x_outer": str(t["x_outer"]), "fit.x_inner": str(t["x_inner"]),
        "fit.x_rot": str(t["x_rot"]), "fit.x_splines": str(t["x_splines"]), "fit.b_outer": str(t["b_outer"]),
        "fit.b_rot": str(t["b_rot"]), "fit.n_hf_beyond": str(len(f["hf_beyond"])),
        "fit.core_rows": count(f["core_rows"]),
        "fit.core_comb": count(f["core_by_category"]["comb"]), "fit.core_comb_sets": str(f["core_sources"]["comb"]),
        "fit.core_doppler": count(f["core_by_category"]["doppler"]),
        "fit.core_doppler_sets": str(f["core_sources"]["doppler"]),
        "fit.core_sr": count(f["core_by_category"]["atlas_sr"]), "fit.core_apo": count(f["core_by_category"]["atlas_apo"]),
        "fit.core_orsay1": count(f["core_by_category"]["orsay1"]), "fit.core_martin": count(f["core_by_category"]["martin"]),
        "fit.core_intervals": count(f["core_kinds"].get("interval", 0)),
        "fit.benchmark_sets": str(f["benchmark_sets"]), "fit.benchmark_rows": str(f["benchmark_rows"]),
        "fit.benchmark_absolute": str(f["benchmark_absolute"]), "fit.benchmark_intraline": str(f["benchmark_intraline"]),
        "fit.benchmark_absolute_rms": f"{f['benchmark_absolute_rms']:.2f}",
        "fit.benchmark_intraline_rms": f"{f['benchmark_intraline_rms']:.3f}",
        "fit.core_folds": str(f["core_folds"]), "fit.core_level_folds": str(f["core_level_folds"]),
        "fit.core_source_folds": str(f["core_source_folds"]), "fit.continuation_folds": str(f["continuation_folds"]),
        "fit.hyperfine_folds": str(f["hyperfine_folds"]),
        "fit.hf_core_trained": count(f["hf_core_trained"]),
        "fit.hf_core_sources": str(f["hf_core_sources"]), "fit.hf_beyond_rows": count(f["hf_beyond_rows"]),
        "fit.hf_chen": str(f["hf_constants"].get("chen2004a", 0)), "fit.hf_pique": str(f["hf_constants"].get("pique1986", 0)),
        "fit.sr_offset": signed(sr["offset_MHz"], 1), "fit.sr_slope": signed(sr["slope_MHz_per_1000cm"], 1),
        "fit.apo_offset": signed(apo["offset_MHz"], 1), "fit.apo_slope": signed(apo["slope_MHz_per_1000cm"], 1),
        "fit.sr_shared": str(sr["n_used"]), "fit.apo_shared": str(apo["n_used"]),
        "fit.orsay4_offset": f"{f['orsay4_offset']:.1f}",
        "fit.n_iso_centers": str(len(f["isotopologue_centers"])),
        "fit.iso_centers": ", ".join(f"{ISO[iso]} {name}" for iso, name in
                                     (line.split(" ", 1) for line, _ in f["isotopologue_centers"])),
        "fit.continuation_rows": count(f["continuation_rows"]),
        "fit.cont_orsay1": count(f["continuation_by_category"]["orsay1"]),
        "fit.cont_orsay4": count(f["continuation_by_category"]["orsay4"]),
        "fit.cont_martin": count(f["continuation_by_category"]["martin"]),
        "fit.cont_martin_dissoc": count(f["continuation_by_category"]["martin_dissoc"]),
        "fit.cont_emission": count(f["continuation_by_category"]["emission"]),
        "fit.region_table": region_table,
        "fit.calibrate_sets": str(f["roles"].get("calibrate", 0)), "fit.locked_sets": str(f["roles"].get("locked test", 0)),
        "fit.record_sets": str(f["roles"].get("record", 0)),
        "fit.exclusion": f["exclusion"].split(" ", 1)[1].replace("-", "–"),
        "fit.orsay1_core_v": "–".join(map(str, f["orsay1_core_v"])), "fit.martin_core_v": "–".join(map(str, f["martin_core_v"])),
        "fit.calibrate_rows": count(uncertainty_numbers()["u"]["checks"]["step4"]["calibrate_coverage"]["total all | refit"]["n"]),
    }


# ------------------------------------------------------------------------------------------------- uncertainty

@lru_cache(maxsize=None)
def uncertainty_numbers():
    """The frozen uncertainty model (the release's uncertainty.json) and its locked test."""
    entry, expanded, core, read = _release()
    u = read(expanded, "uncertainty.json")
    import i2spec.research as research
    # the locked test was scored once, on the release that first carried it: the released one or the latest before
    releases = [r for r in research._RELEASES if r <= entry["id"].split("-")[1]]
    root = Path(research.__file__).parent / "data" / "research"
    folder = next(root / r / "validation" for r in reversed(releases)
                  if (root / r / "validation" / "locked_test_results.json").is_file())
    locked = json.loads((folder / "locked_test_results.json").read_text())
    # the rows whose prediction and sigma the released version leaves as scored (its manifest's record)
    changed = json.loads((root / entry["id"].split("-")[1] / "manifest.json").read_text()).get("locked_test", {})
    return dict(u=u, locked=locked,
                locked_unchanged=changed.get("rows_unchanged", locked["summary"]["all"][locked["primary"]]["n"]))


def _unc_strings(n):
    u, locked = n["u"], n["locked"]
    me, hf = u["model_error"], u["hyperfine"]
    st = u["statistical"]
    cov = u["checks"]["step4"]["calibrate_coverage"]
    b, x = me["B beyond"], me["X beyond"]
    hb = hf["B beyond"]
    late = [cov[f"intraline | B beyond v' {a}-{a + 9}"] for a in (50, 60, 70)]
    total, core = cov["total all | refit"], cov["total 127I2 core | refit"]
    primary = locked["primary"]
    strata = {k.split(" | ")[1]: v[primary] for k, v in locked["summary"].items() if k.startswith("stratum")}
    verdict = {k.split(" | ")[1]: v["result"] for k, v in locked["verdict"].items() if k.startswith("stratum")}
    label = {"absolute": "¹²⁷I₂ absolute frequencies", "interline": "¹²⁷I₂ intervals between lines",
             "intraline": "¹²⁷I₂ hyperfine intervals within a line", "isotopologues": "¹²⁹I₂ and ¹²⁷I¹²⁹I"}
    every = locked["summary"]["all"][primary]
    rows = [f"| all | {every['n']} | {percent(every['within_1sigma'])} | {percent(every['within_2sigma'])} | "
            f"{every['rms_z']:.2f} | {locked['verdict']['all']['result']} |"]
    for k in ("absolute", "interline", "intraline", "isotopologues"):
        s = strata[k]
        result = "not judged" if verdict[k].startswith("reported") else verdict[k]
        rows.append(f"| {label[k]} | {s['n']} | {percent(s['within_1sigma'])} | {percent(s['within_2sigma'])} | "
                    f"{s['rms_z']:.2f} | {result} |")
    return {
        "unc.n_parameters": str(len(st["parameters"])),
        "unc.core_precise": f"{me['core precise']['sigma_MHz']:.2f}",
        "unc.core_precise_groups": str(me["core precise"]["n_groups"]),
        "unc.core_precise_cal": str(me["core precise"]["n_calibrate"]),
        "unc.core_atlas": f"{me['core atlas']['sigma_MHz']:.2f}", "unc.core_atlas_groups": str(me["core atlas"]["n_groups"]),
        "unc.core_atlas_cal": str(me["core atlas"]["n_calibrate"]),
        "unc.b0": f"{math.exp(b['log_sigma0']):.2f}", "unc.b_slope": f"{b['slope_per_v']:.3f}",
        "unc.b_quad": f"{b['quadratic_per_v2']:.3f}", "unc.b_j2": f"{b['J2_coefficient']:.2f}",
        "unc.b_clip": str(b["J_clip"]), "unc.b_groups": str(b["n_groups"]), "unc.b_cal": str(b["n_calibrate"]),
        "unc.x0": f"{math.exp(x['log_sigma0']):.0f}", "unc.x_j2": f"{x['J2_coefficient']:.2f}",
        "unc.x_groups": str(x["n_groups"]),
        "unc.lab_offset": f"{me['lab offset (core absolute, not in predictions)']['sigma_MHz']:.2f}",
        "unc.hf_a": f"{1e3 * hf['core']['a_MHz']:.0f}", "unc.hf_kappa": sci(hf["core"]["kappa"], 1),
        "unc.hf_growth": f"{hb['a_growth_per_v']:.2f}",
        "unc.x_kappa51": f"{hf['X beyond']['kappa_51']:.3f}",
        "unc.cal_rows": count(total["n"]), "unc.cal_1": percent(total["within_1"]), "unc.cal_2": percent(total["within_2"]),
        "unc.cal_core_rows": count(core["n"]), "unc.cal_core_1": percent(core["within_1"]),
        "unc.cal_core_2": percent(core["within_2"]),
        "unc.cal_hf_core_rows": count(cov["intraline | core"]["n"]),
        "unc.cal_hf_core_1": percent(cov["intraline | core"]["refit"][1]),
        "unc.cal_hf_core_2": percent(cov["intraline | core"]["refit"][2]),
        "unc.cal_hf_b_rows": count(sum(c["n"] for c in late)),
        "unc.cal_hf_b_1": f"{min(100 * c['refit'][1] for c in late):.0f}–{max(100 * c['refit'][1] for c in late):.0f} %",
        "unc.cal_hf_b_2": f"{min(100 * c['refit'][2] for c in late):.0f}–{max(100 * c['refit'][2] for c in late):.0f} %",
        "unc.locked_rows": str(every["n"]), "unc.locked_2": percent(every["within_2sigma"]),
        "unc.locked_unchanged": str(n["locked_unchanged"]),
        "unc.locked_abs_2": percent(strata["absolute"]["within_2sigma"]),
        "unc.locked_hf_2": percent(strata["intraline"]["within_2sigma"]),
        "unc.locked_table": "\n".join(rows),
    }


#: lines whose uncertainties the uncertainty page lists: (isotopologue, line, note, the interval it gives)
TYPICAL = [("127I2", "R(56) 32-0", "532 nm", ("a15", "a1")), ("127I2", "R(150) 20-2", "", ("a15", "a1")),
           ("127I2", "R(41) 51-0", "", ("a15", "a1")), ("127I2", "R(98) 58-1", "", ("a15", "a1")),
           ("127I2", "R(20) 70-0", "", ("a15", "a1")), ("127I2", "R(50) 10-30", "", ("a15", "a1")),
           ("127I2", "R(10) 2-90", "", ("a15", "a1")),
           ("129I2", "P(54) 8-4", "633 nm", ("a2", "a1")), ("129I2", "R(20) 15-1", "", ("a2", "a1")),
           ("129I2", "R(56) 32-0", "532 nm", ("a2", "a1")), ("129I2", "R(41) 51-0", "", ("a2", "a1")),
           ("127I129I", "P(33) 6-3", "633 nm", ("a2", "a1")), ("127I129I", "R(56) 32-0", "", ("a2", "a1"))]


def mhz(x):
    """A 1σ as the pages give it: two decimals below 1 MHz, one below 10 MHz, whole MHz above."""
    return f"{x:.2f}" if x < 1 else f"{x:.1f}" if x < 10 else f"{x:.0f}"


@lru_cache(maxsize=None)
def typical():
    """(isotopologue, line, note, sigma center, sigma a1, sigma of the interval, status) of the TYPICAL lines."""
    from i2spec.preferred import model

    out = []
    for iso, line, note, interval in TYPICAL:
        m = model(iso)
        center = m.uncertainty(line)
        out.append((iso, line, note, center.sigma_MHz, m.uncertainty(line, "a1").sigma_MHz,
                    m.uncertainty(line, *interval).sigma_MHz, center.status))
    return out


def _typical_strings(rows):
    label = lambda line, note: line.replace("-", "–") + (f" ({note})" if note else "")  # noqa: E731
    main = [f"| {label(line, note)} | {mhz(c)} MHz | {mhz(i)} MHz | {status} |"
            for iso, line, note, c, a1, i, status in rows if iso == "127I2"]
    other = [f"| {ISO[iso]} {label(line, note)} | {mhz(c)} MHz | {mhz(a1)} MHz | {mhz(i)} MHz | {status} |"
             for iso, line, note, c, a1, i, status in rows if iso != "127I2"]
    return {"unc.typical_table": "\n".join(main), "unc.typical_iso_table": "\n".join(other)}


def _edge_strings(u):
    """The release's terms at the edges (release 20261031 on): X near dissociation, B beyond the tables, the
    tunneling width of quasi-bound levels."""
    nd = u["model_error"]["X beyond"].get("near_dissociation")
    edges = u.get("edges") or {}
    if not nd or not edges:
        return {}
    w = edges["B quasi-bound width"]["evidence"]
    return {"unc.x_nd0": f"{math.exp(nd['log_sigma_add']):.0f}", "unc.x_nd_slope": f"{nd['slope_per_v']:.3f}",
            "unc.x_nd_n": str(nd["n_cv"]), "unc.b_edge_growth": f"{edges['B above 80']['growth_per_v']:g}",
            "unc.width_deep": f"{1e3 * w['max_width_more_than_5_cm_below_top_MHz']:.1f} kHz",
            "unc.width_within": str(w["change_within_width_or_1kHz"]), "unc.width_n": str(w["n_levels"]),
            "unc.width_worst": f"{w['largest_change_over_width']:.1f}"}


def _hf_beyond_strings():
    """The B hyperfine floor a and pattern-scale error κ beyond the core at a few v′, from the model's own
    uncertainty (a of a line to X v″ = 0, whose floor is B's; κ of the B level alone)."""
    from i2spec.preferred import model

    m = model("127I2")
    unc, core = m._uncertainty_model(), m._spec["core"]
    a = lambda v: unc.hyperfine_terms(core, v, 0)[0]  # noqa: E731
    return {"unc.hf_a_beyond": ", ".join(f"{mhz(a(v))} MHz at v′ = {v}" for v in (50, 60, 70)),
            "unc.hf_kappa_beyond": ", ".join(f"{sci(unc.kappa(core, 'B', v), 1)} at v′ = {v}" for v in (44, 55, 70))}


@lru_cache(maxsize=None)
def numbers():
    """Every key a page may use, formatted as the page shows it."""
    unc = uncertainty_numbers()
    return {**_we_strings(worked_example()), **_fit_strings(fit_counts()), **_unc_strings(unc),
            **_edge_strings(unc["u"]),
            **_hf_beyond_strings(), **_typical_strings(typical())}
