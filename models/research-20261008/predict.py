"""Evaluate the frozen 127I2 research candidate without the local search workspace."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from functools import lru_cache
import json
from pathlib import Path
import re

import numpy as np
from numpy.polynomial.chebyshev import chebval
from numpy.polynomial.polynomial import polyvander

from i2spec.bspline import BSplineSolver
from i2spec.constants import MHZ_PER_CM, reduced_mass
from i2spec.hyperfine import HyperfineParameters, level_structure, line_components
from i2spec.observations import component_rank
from i2spec.potentials import MLRPotential

HERE = Path(__file__).resolve().parent


@dataclass(frozen=True)
class Potential(MLRPotential):
    alpha: tuple = ()
    mapped_x: bool = False

    def _exponent(self, y):
        return chebval((y - 0.2) / 0.6 if self.mapped_x else y, self.beta)

    def nonadiabatic(self, R, mass_ratio=1.0):
        y = (R - self.Re) / (R + self.Re)
        return mass_ratio * 2 * self.Re / (R + self.Re) * chebval((y - 0.3) / 0.4, self.alpha)


class SnapshotModel:
    """Shared MLR potentials and radial hyperfine functions, with delta-J=2 mixing."""

    def __init__(self):
        self.spec = json.loads((HERE / "model.json").read_text())
        self.solvers = {}
        for state, curve in self.spec["potentials"].items():
            potential = Potential(
                De=curve["De"], Re=curve["Re"], Te=curve["Te"],
                C={int(k): v for k, v in curve["C"].items()},
                beta=tuple(curve["beta"]), p=curve["p"], q=curve["q"],
                Rref=curve["Rref"], basis="chebyshev", alpha=tuple(curve["alpha"]),
                mapped_x=curve["mapped_x"],
            )
            self.solvers[state] = BSplineSolver(potential, reduced_mass("127I2"), **curve["grid"])

    @lru_cache(maxsize=None)
    def moments(self, state, J):
        solver = self.solvers[state]
        energy, psi = solver.wavefunctions(J)
        curve = self.spec["potentials"][state]
        re = curve["Re"]
        z = (solver.R**2 - re**2) / (solver.R**2 + re**2) / self.spec["hyperfine"]["radial_scale"][state]
        mean = (solver.W[:, None] * psi**2).T @ polyvander(z, 7)
        if np.max(np.abs(mean[:, 0] - 1)) > 1e-10:
            raise ValueError("Radial wavefunction normalization failed")
        return energy, mean

    def energy(self, state, v, J):
        return float(self.moments(state, J)[0][v])

    def parameters(self, state, v, J):
        hfs = self.spec["hyperfine"]
        moments = self.moments(state, J)[1][v]
        coefficients = hfs["coefficients"][state]
        values = {p: np.dot(c, moments[:len(c)]) for p, c in coefficients.items()}
        if state == "X":
            anchor = self.moments("X", 0)[1][0]
            for p in ("eqQ", "C"):
                c = coefficients[p]
                values[p] += hfs["x_ground_anchors"][p] - np.dot(c, anchor[:len(c)])
            values.update(hfs["x_spin_spin"])
        return HyperfineParameters(**values)

    @lru_cache(maxsize=None)
    def structure(self, state, v, J):
        e0 = self.energy(state, v, J)
        return level_structure(
            J, lambda j: self.parameters(state, v, j),
            lambda j: (self.energy(state, v, j) - e0) * MHZ_PER_CM,
            symmetry="g" if state == "X" else "u", dJ=2,
        )

    @lru_cache(maxsize=None)
    def line(self, branch, J, vu, vl):
        if branch not in ("P", "R") or any(int(n) != n or n < 0 for n in (J, vu, vl)):
            raise ValueError("Use a P or R line with nonnegative integer quantum numbers")
        if vu > 43 or vl > 17 or (branch == "P" and J == 0):
            raise ValueError("Snapshot domain: 127I2, B v<=43, X v<=17, P J>=1")
        ju = J + (1 if branch == "R" else -1)
        nu = (self.energy("B", vu, ju) - self.energy("X", vl, J)) * MHZ_PER_CM
        components = line_components(self.structure("B", vu, ju), self.structure("X", vl, J), ju, J)
        return nu, {c.label: c.offset for c in components if c.label}

    def frequency(self, line, component=None):
        nu, offsets = self.line(*line)
        return nu + (offsets[f"a{component_rank(component)}"] if component else 0.0)

    def verify(self):
        references = json.loads((HERE / "reference_predictions.json").read_text())
        maximum = {"absolute": 0.0, "intraline": 0.0}
        counts = {k: 0 for k in maximum}
        for row in references:
            prediction = self.frequency(tuple(row["line"]), row["component"])
            kind = "absolute"
            if "ref_component" in row:
                prediction -= self.frequency(tuple(row["line"]), row["ref_component"])
                kind = "intraline"
            maximum[kind] = max(maximum[kind], abs(prediction - row["predicted_MHz"]) * 1e6)
            counts[kind] += 1
        if counts != {"absolute": 60, "intraline": 592}:
            raise ValueError(f"Incomplete reference benchmark: {counts}")
        for kind, error in maximum.items():
            if error > 100:
                raise ValueError(f"{kind} prediction differs from frozen source by {error:.3f} Hz")
        print(json.dumps(dict(counts=counts, maximum_difference_Hz=maximum), indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("line", nargs="?", help="for example 'P(53) 32-0'")
    parser.add_argument("--component", help="main hyperfine component, e.g. a1")
    parser.add_argument("--verify", action="store_true", help="compare all 60 absolute and 592 intraline references")
    args = parser.parse_args()
    if not args.verify and not args.line:
        parser.error("provide a line or --verify")
    model = SnapshotModel()
    if args.verify:
        model.verify()
    if args.line:
        match = re.fullmatch(r"([PR])\((\d+)\)\s+(\d+)-(\d+)", args.line)
        if not match:
            parser.error("line must have the form 'P(53) 32-0'")
        branch, J, vu, vl = match.groups()
        value = model.frequency((branch, int(J), int(vu), int(vl)), args.component)
        print(json.dumps(dict(line=args.line, component=args.component, frequency_MHz=value,
                              wavenumber_cm1=value / MHZ_PER_CM), indent=2))


if __name__ == "__main__":
    main()
