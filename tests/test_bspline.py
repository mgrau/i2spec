"""B-spline radial solver: analytic Morse levels, convergence at the Hannover joins."""

from dataclasses import dataclass

import numpy as np
import pytest

from i2spec import MHZ_PER_CM, load_potentials, reduced_mass
from i2spec.bspline import BSplineSolver
from i2spec.constants import HBAR2_2U


@dataclass
class Morse:
    De: float = 12547.0
    a: float = 1.86
    Re: float = 2.666
    RI: float = 0.0  # no joins
    RO: float = 99.0

    def __call__(self, R):
        return self.De * (1 - np.exp(-self.a * (np.asarray(R) - self.Re))) ** 2

    def adiabatic(self, R, mass_ratio=1.0):
        return np.zeros_like(np.asarray(R, dtype=float))

    def nonadiabatic(self, R, mass_ratio=1.0):
        return np.zeros_like(np.asarray(R, dtype=float))


def test_morse_levels_are_analytic():
    pot, mu = Morse(), reduced_mass("127I2")
    solver = BSplineSolver(pot, mu, rmin=2.0, rmax=4.5, h=0.01, nlev=12)
    omega = 2 * pot.a * np.sqrt(pot.De * HBAR2_2U / mu)
    v = np.arange(12) + 0.5
    np.testing.assert_allclose(solver.levels(0), omega * v - omega**2 / (4 * pot.De) * v**2, atol=1e-8)


@pytest.fixture(scope="module")
def B_state():
    return load_potentials()["B"], reduced_mass("127I2")


def test_converges_at_the_joins(B_state):
    """B(v'=32, J'=57) sits near the inner join; the sinc-DVR oscillates by ~0.3 MHz there."""
    pot, mu = B_state
    coarse = BSplineSolver(pot, mu, rmin=2.35, rmax=6.0, h=0.012, nlev=34).energy(32, 57)
    fine = BSplineSolver(pot, mu, rmin=2.35, rmax=6.0, h=0.008, nlev=34).energy(32, 57)
    assert abs(coarse - fine) * MHZ_PER_CM < 1e-3  # < 1 kHz
    # the sinc-DVR converged to 19014.4868211(3) cm^-1 at steps <= 0.002 Å (hannover-model-reproduction.md)
    assert fine == pytest.approx(19014.486821110, abs=1e-6)
