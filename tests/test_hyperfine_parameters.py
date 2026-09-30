"""Measured hyperfine parameters (data/hyperfine_parameters/, observations.load_hyperfine_parameters)."""

import numpy as np
import pytest

from i2spec import hfs_params as hp
from i2spec.model import RovibronicModel
from i2spec.observations import PARAMETERS_DIR, load_all_hyperfine_parameters, load_hyperfine_parameters


@pytest.fixture(scope="module")
def chen():
    return load_hyperfine_parameters(PARAMETERS_DIR / "chen2004a")


def test_chen2004_transcription(chen):
    """Values read off the 200 dpi image of Table 1 (continued), and the sign bookkeeping.

    The text layer renders the minus sign as "2"; the 12 values at v' = 69-70 that really begin with 2
    are positive, and 153 printed values are negative.
    """
    assert "chen2004a" in [s.id for s in load_all_hyperfine_parameters()]
    assert len(chen.rows) == 74
    by = {str(r.line): r for r in chen.rows}
    r = by["127I2 P(53) 61-0"]
    assert (r.eqQ, r.C, r.d, r.delta, r.fit_sd_kHz) == (-570.037, 0.906087, -0.4162, 0.5154, 44.0)
    assert (r.eqQ_unc, r.C_unc, r.d_unc, r.delta_unc) == pytest.approx((0.079, 0.000075, 0.0042, 0.0042))
    r = by["127I2 R(45) 70-0"]
    assert (r.eqQ, r.C, r.d, r.delta) == (-558.80, 2.22994, -0.937, 3.045) and r.J_upper == 46
    assert by["127I2 P(10) 42-0"].eqQ_unc == pytest.approx(0.020)
    assert sum((r.eqQ < 0) + (r.C < 0) + (r.d < 0) + (r.delta < 0) for r in chen.rows) == 153
    assert sum((2 <= r.C < 3) + (2 <= r.delta < 3) for r in chen.rows) == 12      # the genuine leading 2s
    assert sorted(str(r.line) for r in chen.rows if r.line.v_lower) == ["127I2 P(49) 59-1", "127I2 R(45) 60-1"]
    excluded = {"P(69) 58-0", "R(18) 59-0", "P(84) 60-0", "P(77) 60-0", "P(63) 70-0"}   # perturbed, left out by the authors
    assert not excluded & {str(r.line).split(" ", 1)[1] for r in chen.rows}


def test_chen2004_agrees_with_bkt02_where_it_is_valid(chen):
    """At v' = 42-43, inside BKT02's range, C_B agrees with the formula to 0.46 kHz rms and eqQ_B to
    0.31 MHz; P(13) 43-0 agrees with the fit to bipm2005a's splittings (+0.066 +- 0.005 kHz)."""
    m = RovibronicModel("127I2")
    dC, dq = [], []
    for r in chen.rows:
        if r.line.v_upper > 43:
            continue
        E_b, E_x = m._reference_term("B", r.line.v_upper), m._reference_term("X", r.line.v_lower)
        assert E_b <= hp.E_B_MAX
        p = hp.line_states("127I2", r.line.v_upper, r.line.v_lower, E_b, E_x)[1](r.J_upper)
        dC.append(r.C * 1e3 - p.C)
        dq.append(r.eqQ - p.eqQ)
        if str(r.line) == "127I2 P(13) 43-0":
            assert abs(dC[-1] - 0.066) < 2 * r.C_unc * 1e3
    assert len(dC) == 11
    assert np.sqrt(np.mean(np.square(dC))) == pytest.approx(0.46, abs=0.02)
    assert np.sqrt(np.mean(np.square(dq))) == pytest.approx(0.31, abs=0.02)


def test_loader_rejects_bad_files(tmp_path):
    d = tmp_path / "x2000"
    d.mkdir()
    (d / "meta.toml").write_text('id = "x2000"\ncitation = "c"\nsource = "s"\nretrieved = "r"\nunit = "MHz"\ntranscription = "t"\n')
    (d / "data.csv").write_text("line,eqQ\n127I2 R(1) 1-0,1\n")
    with pytest.raises(ValueError, match="columns"):
        load_hyperfine_parameters(d)
    (d / "data.csv").write_text("line,eqQ,eqQ_unc,C,C_unc,d,d_unc,delta,delta_unc,fit_sd_kHz,note\n"
                                "127I2 R(1) 1-0,-500,0,1,1,1,1,1,1,,\n")
    with pytest.raises(ValueError, match="positive"):
        load_hyperfine_parameters(d)
    (d / "meta.toml").write_text('id = "other"\ncitation = "c"\nsource = "s"\nretrieved = "r"\nunit = "MHz"\ntranscription = "t"\n')
    with pytest.raises(ValueError, match="directory name"):
        load_hyperfine_parameters(d)
