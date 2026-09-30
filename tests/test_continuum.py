"""Continuum absorption against Tellinghuisen, J. Chem. Phys. 135, 054301 (2011), and an exact case."""

import numpy as np
import pytest
from scipy.special import airy

from i2spec.constants import HBAR2_2U, reduced_mass
from i2spec.continuum import Continuum, _overlaps, continuum_model
from i2spec.spectrum import air_to_vacuum

EPS = 3.82353e-21  # cm² per L mol⁻¹ cm⁻¹
#: Table II, ε at 0 °C and 35 °C. Its wavelengths fit best as air wavelengths (0.2-0.6 % at
#: 405-495 nm, against 0.5-2 % read as vacuum; docs/research/continuum-model.md §10).
TABLE2 = {405: (1.38, 1.89), 420: (7.08, 8.83), 435: (26.22, 30.24), 450: (74.08, 80.85), 465: (170.00, 178.24),
          480: (328.57, 332.99), 495: (533.16, 523.97), 700: (38.55, 37.07), 720: (31.81, 30.98),
          740: (24.04, 24.06), 760: (16.80, 17.45), 780: (10.99, 11.95), 800: (6.79, 7.79), 850: (1.69, 2.26)}


@pytest.fixture(scope="module")
def continuum():
    return continuum_model("127I2", T_max=320.0, dJ=10, cache=False)


def test_energy_normalization_linear_potential():
    """For V = F (R0 - R) the energy-normalized solution is c^(-1/3) F^(-1/6) Ai(-(F/c)^(1/3) (R - R_t))."""
    c, F, h, stride = HBAR2_2U / reduced_mass("127I2"), 2e4, 0.001, 4
    R = 2.0 + np.arange(3001) * h
    E = np.array([0.0, 2000.0, 5000.0])
    samples = R[::stride][:600]
    weight = np.exp(-(((samples - 2.6) / 0.1) ** 2))[:, None]
    m2 = _overlaps(F * (2.5 - R), h, E, c, stride, samples.size, (2000, 2500), weight)
    Rt = 2.5 - E / F
    exact = [c ** (-1 / 3) * F ** (-1 / 6) * airy(-((F / c) ** (1 / 3)) * (samples - r))[0] @ weight[:, 0] for r in Rt]
    assert np.sqrt(m2[:, 0]) == pytest.approx(np.abs(exact), rel=5e-4)  # Numerov phase error, h = 0.001 Å


@pytest.mark.parametrize("T, col", [(273.15, 0), (308.15, 1)])
def test_tellinghuisen_table2(continuum, T, col):
    lam = np.array(sorted(TABLE2), dtype=float)
    ratio = continuum.cross_section(1e7 / air_to_vacuum(lam), T) / EPS / [TABLE2[int(l)][col] for l in lam]
    assert np.abs(ratio[lam < 500] - 1).max() < 0.01  # A + C + B bound-free
    assert np.abs(ratio[lam > 690] - 1).max() < 0.015  # A←X; beyond 810 nm the A curve is a stand-in


def test_436nm(continuum):
    """Room-temperature value calculated in the paper: 31.0 ± 0.4 L mol⁻¹ cm⁻¹."""
    assert float(continuum.cross_section(1e7 / air_to_vacuum(436.0), 298.0)) / EPS == pytest.approx(31.0, abs=0.6)


def test_save_load(continuum, tmp_path):
    continuum.save(tmp_path / "c.npz")
    again = Continuum.load(tmp_path / "c.npz")
    nu = np.linspace(15000.0, 24000.0, 50)
    assert np.array_equal(again.cross_section(nu, 300.0), continuum.cross_section(nu, 300.0))
