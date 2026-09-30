"""Lever 4: the isotopologue lines against the published Born-Oppenheimer corrections.

Usage:  uv run python prototypes/bo_refit.py            (about 5 minutes)
Output: prototypes/out/bo_refit.json

Every 129I2 line sits 7-14 MHz above the published model and the mixed isotopologue 4.6 MHz below
129I2, each with sub-MHz scatter inside the line. The published set puts all Born-Oppenheimer terms in
the B state (X has none): alpha(R), six terms, in the centrifugal part, and V_ad(R), four terms, as
(1 - mu_ref/mu) (2 Rm/(R + Rm))^5 sum_i vad_i x^i. This fits the leading terms of V_ad (and alpha) to
the 124 isotopologue rows that constrain the potentials, 127I2 untouched by construction.
"""

import json
from dataclasses import replace
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

from i2spec.constants import MHZ_PER_CM
from i2spec.observations import Predictor, load_dataset
from i2spec.potentials import load_potentials

ROOT = Path(__file__).resolve().parents[1]
B0 = load_potentials()["B"]


def rows():
    out = []
    for name in ("bipm2003a", "bipm2003d"):
        ds = load_dataset(ROOT / "data/observations" / name)
        scale = MHZ_PER_CM if ds.meta["unit"] == "cm-1" else 1.0
        for o in ds.observations:
            ref = o.ref_line or o.line
            if {o.line.isotopologue, ref.isotopologue} == {"127I2"} or (o.kind == "interval" and ref == o.line):
                continue
            out.append((name, o, scale))
    return out


def frozen_offsets(data):
    """Hyperfine offsets from the published model, held fixed: they depend on the BO terms only through
    the level energies in the S06 pole terms, far below the MHz at stake here."""
    pred = Predictor()
    return [pred(o, "MHz") - _centre(pred, o) for _, o, _ in data]


def _centre(pred, o):
    c = pred.position(o.line)
    return c - (pred.position(o.ref_line or o.line) if o.kind == "interval" else 0.0)


def residuals(vad_delta, alpha_delta, data, offsets):
    vad = tuple(v + d for v, d in zip(B0.vad, list(vad_delta) + [0.0] * (len(B0.vad) - len(vad_delta))))
    alpha = tuple(a + d for a, d in zip(B0.alpha, list(alpha_delta) + [0.0] * (len(B0.alpha) - len(alpha_delta))))
    pred = Predictor(potentials={"B": replace(B0, vad=vad, alpha=alpha)})
    model = np.array([_centre(pred, o) + off for (_, o, _), off in zip(data, offsets)])
    return np.array([o.value * scale for _, o, scale in data]) - model * np.array([scale for _, _, scale in data]), \
        np.array([o.uncertainty * scale for _, o, scale in data])


def main():
    data = rows()
    offsets = frozen_offsets(data)
    r0, sig = residuals([], [], data, offsets)
    w = 1 / np.hypot(sig, 0.5)
    print(f"{len(data)} rows; published: rms {np.sqrt(np.mean(r0**2)):.2f} MHz")
    results = {"published": dict(rms=float(np.sqrt(np.mean(r0 ** 2))))}
    for label, n_vad, n_alpha in (("vad constant", 1, 0), ("vad constant + slope", 2, 0),
                                  ("vad constant + slope + quadratic", 3, 0), ("vad constant + slope, alpha constant", 2, 1)):
        x0 = np.zeros(n_vad + n_alpha)
        fit = least_squares(lambda x: residuals(x[:n_vad], x[n_vad:], data, offsets)[0] * w, x0, diff_step=1e-4,
                            x_scale=np.r_[np.full(n_vad, 0.02), np.full(n_alpha, 1e-4)], max_nfev=200, xtol=1e-10, ftol=1e-10)
        r, _ = residuals(fit.x[:n_vad], fit.x[n_vad:], data, offsets)
        by = {}
        for (name, o, scale), e in zip(data, r):
            key = f"{o.line}" + ("" if o.kind == "frequency" else f" − {(o.ref_line or o.line).isotopologue}")
            by.setdefault(key, []).append(e)
        results[label] = dict(x=fit.x.tolist(), rms=float(np.sqrt(np.mean(r ** 2))),
                              by_line={k: float(np.mean(v)) for k, v in by.items()})
        print(f"{label}: x = {np.round(fit.x, 5).tolist()}; rms {results[label]['rms']:.2f} MHz")
        for k, v in by.items():
            print(f"     {k:<44} mean {np.mean(v):+7.2f} MHz  (n={len(v)})")
    (ROOT / "prototypes/out").mkdir(exist_ok=True)
    (ROOT / "prototypes/out/bo_refit.json").write_text(json.dumps(results, indent=1))


if __name__ == "__main__":
    main()
