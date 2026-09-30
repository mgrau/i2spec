"""Milestone 1: the published Hannover 2008 potentials reproduce known results.

Parameters: Salumbides et al., EPJD 47, 171 (2008), Table 1. Checks against
Gerstenkorn, Luc & Le Roy, Can. J. Phys. 69, 1299 (1991) (T0, zero-point energy) and the
local NIR Dunham model of Knöckel, Bodermann & Tiemann, EPJD 28, 199 (2004), Table 1.
"""

import numpy as np
import pytest
from scipy.optimize import minimize_scalar

from i2spec import MHZ_PER_CM, RovibronicModel, load_potentials, reduced_mass
from i2spec.local_nir import load_local_nir
from i2spec.solver import RadialSolver

# Knöckel 2004 Table 1: bands 0-v'' (v'' = 12..17), hyperfine-free, stated 1σ < 200 kHz. The
# parameters live in src/i2spec/data/knockel2004_nir.json, so there is one copy of them; this test
# and tests/test_local_nir.py share it.
LOCAL_NIR = load_local_nir("knockel2004")


def dunham_nir_r_line(v_lower, J_lower):
    return LOCAL_NIR.transition(0, v_lower, J_lower, "R")


@pytest.fixture(scope="module")
def model():
    return RovibronicModel("127I2")


def test_extensions_join_the_series():
    for pot in load_potentials().values():
        for R in (pot.RI, pot.RO):
            left, right = pot(np.array([R - 1e-9, R + 1e-9]))
            assert abs(left - right) < 5e-4


def test_zero_point_energy_and_T0(model):
    e00 = model.energy("X", 0, 0)
    assert e00 == pytest.approx(107.101, abs=0.005)  # De - D0 (GLL 1991)
    assert model.energy("B", 0, 0) - e00 == pytest.approx(15724.587, abs=0.002)  # T0 (GLL 1991)


@pytest.mark.parametrize("v_lower", [14, 15])
@pytest.mark.parametrize("J_lower", [20, 50, 80, 100, 120])
def test_central_nir_bands_match_local_dunham_model(model, v_lower, J_lower):
    diff = (model.transition(0, v_lower, J_lower, "R") - dunham_nir_r_line(v_lower, J_lower)) * MHZ_PER_CM
    assert abs(diff) < 5.0


@pytest.mark.parametrize(
    "v_upper, v_lower, J_lower, f_khz",
    [
        (32, 0, 56, 563_260_223_513),  # BIPM MeP 532 nm, R(56)32-0 a10
        (11, 5, 127, 473_612_353_604),  # BIPM MeP 633 nm, R(127)11-5 a16 (f)
        (16, 1, 37, 518_304_551_833),  # Hong et al. 2009, R(37)16-1 a1
    ],
)
def test_absolute_components_lie_within_hyperfine_spread(model, v_upper, v_lower, J_lower, f_khz):
    diff = model.transition(v_upper, v_lower, J_lower, "R") * MHZ_PER_CM - f_khz / 1e3
    assert abs(diff) < 600.0


def test_v0_levels_converged_to_1_khz():
    potentials, mu = load_potentials(), reduced_mass("127I2")
    grids = {"X": dict(rmin=2.10, rmax=4.0, step=0.004), "B": dict(rmin=2.35, rmax=8.0, step=0.005)}
    for state, g in grids.items():
        coarse = RadialSolver(potentials[state], mu, nlev=2, **g).energy(0, 0)
        fine = RadialSolver(potentials[state], mu, nlev=2, **{**g, "step": 0.8 * g["step"]}).energy(0, 0)
        assert abs(coarse - fine) * MHZ_PER_CM < 1e-3


def test_isotopic_Te_shift_matches_salumbides_2008():
    """Te(127I2) - Te(129I2) = 94(11) MHz (Salumbides et al. 2008)."""
    potentials = load_potentials()

    def te(isotopologue):
        ratio = reduced_mass("127I2") / reduced_mass(isotopologue)
        minima = {}
        for state, bounds in (("B", (2.9, 3.2)), ("X", (2.55, 2.8))):
            pot = potentials[state]
            res = minimize_scalar(lambda r: float(pot(r) + pot.adiabatic(r, ratio)), bounds=bounds,
                                  method="bounded", options={"xatol": 1e-10})
            minima[state] = res.fun
        return minima["B"] - minima["X"]

    assert (te("127I2") - te("129I2")) * MHZ_PER_CM == pytest.approx(94, abs=11)
