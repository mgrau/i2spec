"""Rebuild a release from its frozen inputs: every stage in order, then the release files.

  1  core          potentials (core.fit) and B hyperfine functions (core_hyperfine.fit) on the frozen core rows
  2  continuation  tail terms on the frozen continuation rows (continuation.fit)
  3  hyperfine     beyond the core: B background + second-order sums (hyperfine_beyond.fit), X eqQ and C
                   (x_hyperfine.fit), the B~1g(1Pig) resonance (b1g_resonance)
  4  package       levels, hyperfine tables, resonance tables, reference predictions (package)
  5  uncertainty   covariance, sigma_model and sigma_HF (uncertainty.build)
  6  isotopologues the isotope shift (isotopologues.fit_shift) and its record

The core and continuation fits start from the released coefficients (a refit returns them: tests/test_fit.py);
`start` replaces them. Each stage's record is returned, and `write` puts the files of a continued release in a
folder. The hyperfine stage of the release was fitted on an earlier continuation (frozen with it as
hyperfine_potential.json), and the resonance tables now carry consistent wavefunction phases (b1g_resonance's
overlap), so those two records differ from the release by construction; every other number reproduces it.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np

from ..research import ResearchModel, _entry, _read
from . import b1g_resonance as rz, continuation, core, core_hyperfine, hyperfine_beyond as hb, isotopologues
from . import package, reduction, rows as training, uncertainty, x_hyperfine as xh
from .curves import Config, Curve


def stage_core(settings=None):
    prob, result = core.fit(settings=settings)
    spec = core.release_spec()
    values = {s: dict(spec["values"][s]) for s in ("X", "B")}
    for (s, n), z in zip(result.names, result.theta):
        values[s][n] = float(z)
    config = Config.from_spec(spec["config"])
    curves = {s: Curve(s, config, values[s]) for s in ("X", "B")}
    records = [r for r in hyperfine_rows_release() if r["inside"]]
    energies, means = core_hyperfine.moments(curves, values, [r["line"] for r in records])
    hf_problem = core_hyperfine.Problem(energies, means, records)
    hf = core_hyperfine.fit(hf_problem, start=core_hyperfine.release_coefficients(spec, hf_problem.names))
    return dict(potential=result, values=values, hyperfine=hf, hyperfine_names=hf_problem.names)


def hyperfine_rows_release():
    from . import hyperfine_rows
    return hyperfine_rows.release()[0]


def stage_continuation(**kwargs):
    prob, result = continuation.fit(**kwargs)
    return dict(problem=prob, result=result)


def stage_hyperfine(workers=12, cache=None):
    rows, params = hb.load_data()
    constants = hb.state_constants(False, hb.FIT_OVERRIDES)
    levels = hb.Levels()
    model = hb.Model(levels, hb.fill_terms(hb.needed_levels(rows, params), constants, workers, cache=cache))
    b, theta = hb.fit(model, rows, params, dict(u=30., g=30.), False, keys=list(hb.FIT_KEYS))
    xlv = xh.XLevels()
    xmodel = xh.fit(xlv, xh.KNOT, xh.ORDER)
    return dict(constants=constants, b=np.asarray(b), theta=theta, training=model.training, levels=levels,
                xlv=xlv, xmodel=xmodel)


def x_model_record(xmodel):
    return dict(tag="hyperfine stage continuation", knot=xh.KNOT, order=xh.ORDER, inner=xh.INNER,
                eqQ_limit=xh.EQQ_LIMIT,
                eqQ=dict(coef=xmodel["coef"].tolist(), lam=xmodel["lam"], p=xmodel["p"],
                         anchor_rms_MHz=xmodel["anchor_rms"], anchor_max_MHz=xmodel["anchor_max"]),
                C=dict(coef=xmodel["coef_c"].tolist(), slope=xmodel["slope_c"], unit="kHz",
                       anchor_rms_kHz=xmodel["anchor_c_rms"],
                       note="slope fitted to the single value 59 kHz at v83; unvalidated"))


def stage_package(cont, hf, workers=12, cache=None):
    every = {f"{s}:{n}": z for s in "XB" for n, z in continuation.release_spec()["values"][s].items()}
    every.update(cont["result"].parameters)
    every["B:dDe"] = every["X:dDe"]                    # one shared shift
    structure, values = package.structure_spec({k: v for k, v in every.items() if not k.startswith("calibration")})
    E = package.levels(values, workers=workers)
    tables = package.hyperfine_tables(E, hf["b"], hf["theta"], hf["constants"], hf["xmodel"], hf["xlv"],
                                      blv=hf["levels"], workers=workers, cache=cache)
    res = rz.Resonance(hf["constants"], hf["b"], hf["theta"], levels=hf["levels"], cache=cache)
    resonance = package.resonance_table(E, res)
    own = reduction.OwnHyperfine(hf["constants"], hf["b"], hf["theta"], hf["xmodel"], xlv=hf["xlv"], cache=cache)
    res2 = rz.Resonance(hf["constants"], hf["b"], hf["theta"], cache=cache)
    references = package.references(E, tables, resonance, own, res2, ResearchModel("balanced"),
                                    _read(_entry("balanced"), "reference_predictions.json"))
    return dict(structure=structure, values=values, levels=E, tables=tables, resonance=resonance,
                references=references)


def rebuild(workers=12, cache=None):
    """Every stage from the frozen inputs: {stage: record}. About five minutes with a warm second-order cache."""
    out = dict(core=stage_core())
    out["continuation"] = stage_continuation()
    out["hyperfine"] = hf = stage_hyperfine(workers, cache)
    out["package"] = stage_package(out["continuation"], hf, workers, cache)
    out["uncertainty"], _ = uncertainty.build(cache=cache, workers=workers)
    out["isotope_shift"] = isotopologues.fit_shift()
    return out


def records(built):
    """The continued release's file contents (model.json, hyperfine.json, reference_predictions.json,
    uncertainty.json) from `rebuild`'s stages."""
    pk, hf, cont = built["package"], built["hyperfine"], built["continuation"]
    res = cont["result"]
    b_model = dict(result="the release's hyperfine fit", theta_deg=hf["theta"], potential="hyperfine stage continuation",
                   grid=hb.so.GRID, curve_constants=hf["constants"], background=hf["b"].tolist(),
                   training=hf["training"], holdout=[])
    x_model = dict(x_model_record(hf["xmodel"]), potential="hyperfine stage continuation",
                   grid=dict(xh.GRID, nlev=xh.NLEV))
    base_spec = ResearchModel("balanced").coefficients
    hyperfine = package.hyperfine_record(pk["tables"], pk["resonance"], b_model, x_model,
                                         base_spec["hyperfine"]["x_spin_spin"])
    counts = package.reference_counts(pk["references"])
    train = training.release("continuation")
    fit = dict(cost=res.cost, converged=res.converged, nfev=res.nfev,
               training=continuation.statistics(train, train.value-res.pred))
    base_entry = _entry("balanced")
    derived = dict(B_background=len(hf["b"]), X_eqQ_polynomial=len(hf["xmodel"]["coef"]),
                   X_C_polynomial=len(hf["xmodel"]["coef_c"]))
    model = package.model_record(pk["structure"], pk["values"], res.parameters["calibration:orsay4"], fit,
                                 hf["constants"], counts, base_entry, derived)
    return {"model.json": model, "hyperfine.json": hyperfine, "reference_predictions.json": pk["references"],
            "uncertainty.json": built["uncertainty"]}


def write(folder, files):
    """Write the files (compact JSON for the large tables) and return their SHA256."""
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    out = {}
    for name, value in files.items():
        compact = name in ("hyperfine.json", "uncertainty.json")
        text = json.dumps(value, separators=(",", ":") if compact else None, indent=None if compact else 1,
                          allow_nan=False)+"\n"
        (folder/name).write_text(text)
        out[name] = hashlib.sha256(text.encode()).hexdigest()
    return out
