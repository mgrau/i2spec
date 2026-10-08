"""Reference values from the Python model, for web/test_physics.mjs to check the browser code against."""
import json

import numpy as np

from i2spec import RovibronicModel, saturation, spectrum
from i2spec.intensity import C2, intensity_model, master_line_list
from i2spec.lookup import Catalog, uncertainty


def make_reference():
    """Compute reference values only in the parent process."""
    out = {
        "air": [[lam, spectrum.air_refractive_index(lam), spectrum.vacuum_to_air(lam),
                 spectrum.air_to_vacuum(lam)] for lam in (400.0, 532.0, 633.0, 800.0, 1000.0)],
        "vapor": [[T, spectrum.vapor_pressure(T), spectrum.number_density(T, 293.15)]
                  for T in (253.15, 273.15, 293.15, 313.15)],
        "doppler": [[nu, T, spectrum.doppler_fwhm(nu, T)]
                    for nu in (12600.0, 18788.0) for T in (200.0, 293.15, 600.0)],
    }

    master = master_line_list(intensity_model("127I2"), 11000.0, 20100.0, T_range=(200.0, 600.0), S_min=1e-27)
    out["Q"] = [[T, float(master.partition_function(T))] for T in (200.0, 250.0, 293.15, 300.0, 400.0, 600.0)]

    lines = master.at(300.0, S_min=1e-24)
    order = np.argsort(lines.nu)
    s0 = lines.S[order] * master.partition_function(300.0) * np.exp(C2 * lines.E_lower[order] / 300.0)
    nu, el = lines.nu[order], lines.E_lower[order]
    sel = [k for k in range(len(nu)) if 18787.0 <= nu[k] <= 18789.0]
    sel = sorted(sel, key=lambda k: -lines.S[order][k])[:6]
    out["strengths"] = [
        {"nu": float(nu[k]), "el": float(el[k]), "s0": float(s0[k]),
         "S": {str(T): float(s0[k] * np.exp(-C2 * el[k] / T) / master.partition_function(T))
               for T in (200.0, 293.15, 300.0, 500.0)}}
        for k in sel]
    # the per-line uncertainty the export carries, for lines the export always has (they are measured)
    catalog = Catalog()
    out["uncertainty"] = []
    for label in ("R(56) 32-0", "P(13) 43-0", "P(82) 0-13", "P(119) 1-24", "P(52) 53-0"):
        line = catalog.line(label)
        out["uncertainty"].append({"label": label, "nu": line.nu, "u": uncertainty(line)[0], "flags": list(line.flags)})

    # sub-Doppler: the full hyperfine calculation of R(56) 32-0, with the levels the crossovers need (eigenstate
    # numbers, as saturation.shared_level uses them and as the export gives them: webapp.weak_links), and
    # what saturation.py makes of it, for the browser's copy of the same model to reproduce
    nu0, comps = RovibronicModel("127I2").hyperfine_components(32, 0, 56, "R")
    width = saturation.doppler_width("127I2", nu0, 300.0)
    out["sub_doppler"] = {
        "line": "R(56) 32-0", "doppler_width": width,
        "components": [{"offset": c.offset, "strength": c.strength, "lower": c.lower_level,
                        "upper": c.upper_level} for c in comps],
        "resonances": [], "lineshape": [], "signal": [],
    }
    for options in ({}, {"exponent": 1.0}, {"lambda_weight": 0.2, "threshold": 1e-3}):
        found = saturation.resonances(comps, width, **options)
        out["sub_doppler"]["resonances"].append({"options": options, "list": [
            {"offset": r.offset, "amplitude": r.amplitude, "kind": r.kind, "sharing": r.sharing,
             "components": list(r.components)} for r in found]})
    x = [-7.3, -2.0, -0.61, 0.0, 0.25, 1.0, 3.9, 40.0]
    for fwhm, harmonic, modulation in ((2.0, 0, 0.0), (2.0, 1, 1.0), (2.0, 3, 1.0), (0.5, 3, 2.0), (1.0, 2, 0.3)):
        out["sub_doppler"]["lineshape"].append({"fwhm": fwhm, "harmonic": harmonic, "modulation": modulation, "x": x,
            "y": saturation.lineshape(np.array(x), fwhm, harmonic=harmonic, modulation=modulation).tolist()})
    found = saturation.resonances(comps, width)
    grid = np.linspace(-600.0, 600.0, 241).tolist() + [c.offset + d for c in comps[:4] for d in (-0.4, 0.0, 0.3)]
    for fwhm, harmonic, modulation in ((2.0, 0, 0.0), (2.0, 3, 1.0)):
        out["sub_doppler"]["signal"].append({"fwhm": fwhm, "harmonic": harmonic, "modulation": modulation, "nu": grid,
            "y": saturation.signal(np.array(grid), found, fwhm, harmonic=harmonic, modulation=modulation).tolist()})
    return out


def main():
    print(json.dumps(make_reference(), allow_nan=False))


if __name__ == "__main__":
    main()
