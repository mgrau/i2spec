"""532 nm hyperfine structure against the BIPM mise en pratique tables.

Model: Hannover 2008 potentials plus the Bodermann et al. (2002) hyperfine formulae.
Reference: BIPM, recommended values of standard frequencies, iodine (≈532 nm) (CIPM 2007, file
updated 2012), Tables 1 and 16. Bodermann et al. state prediction uncertainties of about 20-30 kHz
for F - J = 0 components and 0.1-1 MHz for the others; the tolerances below allow for that.
"""

import numpy as np
import pytest

from i2spec import RovibronicModel

# [f(an) - f(a10)] / MHz, R(56)32-0 (Table 16, u = 1.5 kHz); a3 and a4 are not tabulated.
R56_32_0 = {"a1": -571.542, "a2": -311.844, "a5": -260.176, "a6": -170.064, "a7": -154.548, "a8": -131.916,
            "a9": -116.199, "a10": 0.0, "a11": 126.513, "a12": 131.212, "a13": 154.488, "a14": 160.665, "a15": 286.412}
# [f(an) - f(a1)] / MHz, R(87)33-0 (Table 1, u = 2 kHz).
R87_33_0 = dict(zip([f"a{n}" for n in range(1, 22)],
                    [0.0, 51.5768, 101.4407, 282.4331, 332.2313, 342.2223, 390.3168, 445.6559, 462.0620, 497.5450,
                     511.9546, 582.6721, 622.8375, 663.9140, 730.3226, 752.4797, 778.0522, 799.4548, 893.1211,
                     907.5209, 923.5991]))
F_A10 = 563_260_223.513  # MHz, CIPM recommended R(56)32-0 a10 (u = 5 kHz)


@pytest.fixture(scope="module")
def model():
    return RovibronicModel("127I2")


def main_offsets(model, v_upper, J_lower):
    nu0, components = model.hyperfine_components(v_upper, 0, J_lower, "R")
    return nu0, {c.label: c.offset for c in components if c.label}


def test_r56_32_0_intervals(model):
    _, o = main_offsets(model, 32, 56)
    residuals = np.array([(o[k] - o["a10"]) - v for k, v in R56_32_0.items()]) * 1e3
    assert np.abs(residuals).max() < 60.0  # kHz


def test_r87_33_0_intervals(model):
    _, o = main_offsets(model, 33, 87)
    residuals = np.array([(o[k] - o["a1"]) - v for k, v in R87_33_0.items()]) * 1e3
    assert np.abs(residuals).max() < 200.0  # kHz


def test_a10_absolute_frequency(model):
    nu0, o = main_offsets(model, 32, 56)
    assert abs(nu0 + o["a10"] - F_A10) < 5.0  # MHz; rovibronic model ≈ 1.5 MHz (1σ, Salumbides 2008)
