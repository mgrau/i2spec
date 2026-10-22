"""532 nm hyperfine structure against the BIPM mise en pratique tables.

Model: the released model. Reference: BIPM, recommended values of standard frequencies, iodine (≈532 nm)
(CIPM 2007, file updated 2012), Tables 1 and 16. The model's held-out hyperfine intervals are good to 0.06 MHz rms;
the worst component of these two lines is 0.09 MHz off.
"""

import numpy as np
import pytest

from i2spec.preferred import model

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
def released():
    return model("127I2")


def main_offsets(m, v_upper, J_lower):
    nu0, components = m.hyperfine_components(v_upper, 0, J_lower, "R")
    return nu0, {c.label: c.offset for c in components if c.label}


def test_r56_32_0_intervals(released):
    _, o = main_offsets(released, 32, 56)
    residuals = np.array([(o[k] - o["a10"]) - v for k, v in R56_32_0.items()]) * 1e3
    assert np.abs(residuals).max() < 120.0  # kHz


def test_r87_33_0_intervals(released):
    _, o = main_offsets(released, 33, 87)
    residuals = np.array([(o[k] - o["a1"]) - v for k, v in R87_33_0.items()]) * 1e3
    assert np.abs(residuals).max() < 120.0  # kHz


def test_a10_absolute_frequency(released):
    nu0, o = main_offsets(released, 32, 56)
    assert abs(nu0 + o["a10"] - F_A10) < 0.5  # MHz
