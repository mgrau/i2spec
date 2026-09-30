"""Measured B-state hyperfine corrections by (v', J') (src/i2spec/hfs_table.py)."""

import json

import numpy as np
import pytest

from i2spec.hfs_table import PARAMS, TABLE_PATH, V_FROZEN, HyperfineTable, default_table
from i2spec.hyperfine import HyperfineParameters
from i2spec.model import RovibronicModel


@pytest.fixture(scope="module")
def table():
    return default_table()


def test_table_covers_the_measured_lines(table):
    rows = json.loads(TABLE_PATH.read_text())["rows"]
    assert len(rows) >= 120 and len(table.rows) == len(rows)
    assert {"chen2004a", "bipm2012a", "yoshiki2023a", "matsunaga2024a"} <= {r["source"] for r in rows}
    assert 32 in table.corrections and 44 in table.corrections and 70 in table.corrections


def test_a_single_measured_line_is_reproduced_exactly(table):
    """v' = 28 has one line, R(106) 28-0: its correction is a constant equal to measured - formula."""
    row, = [r for r in table.rows if r["v"] == 28]
    for k in PARAMS:
        assert table.correction(k, 28, row["J"]) == pytest.approx(row["measured"][k] - row["formula"][k], rel=1e-9)
        assert table.correction(k, 28, 5) == table.correction(k, 28, 150)    # degree 0


def test_a_well_measured_v_interpolates_in_J(table):
    """v' = 32 has six BIPM lines at J' = 52-59; the correction there is a straight line in J'(J'+1)
    that passes within its own scatter of every measured value."""
    rows = [r for r in table.rows if r["v"] == 32]
    assert len(rows) >= 6
    for r in rows:
        assert abs(table.correction("C", 32, r["J"]) - (r["measured"]["C"] - r["formula"]["C"])) < 0.05   # kHz


def test_unmeasured_v_below_the_freeze_keeps_the_formula(table):
    assert 30 not in table.corrections and 30 <= V_FROZEN
    assert all(table.correction(k, 30, 40) == 0.0 for k in PARAMS)


def test_above_the_freeze_the_corrections_interpolate(table):
    """v' = 58 lies between Chen's v' = 57 and 59: its C_B correction is between theirs, and large,
    since the published formulae are frozen at v' = 53 while C_B keeps rising."""
    c57, c58, c59 = (table.correction("C", v, 60) for v in (57, 58, 59))
    assert min(c57, c59) <= c58 <= max(c57, c59)
    assert abs(c58) > 50.0                                                    # kHz


def test_leaving_a_line_out_removes_it(table):
    without = HyperfineTable.load(exclude=("R(56) 32-0",))
    assert len(without.rows) == len(table.rows) - 1
    assert 32 in without.corrections


def test_apply_and_the_model_hook(table):
    p = HyperfineParameters(eqQ=-560.0, C=90.0, d=-40.0, delta=-6.0)
    q = table.apply(p, 32, 57)
    assert q.C != p.C and abs(q.C - p.C) < 1.0
    m = RovibronicModel("127I2")
    _, plain = m.hyperfine_components(32, 0, 56, "R", table=None)
    _, corrected = m.hyperfine_components(32, 0, 56, "R", table=table)
    moved = max(abs(a.offset - b.offset) for a, b in zip(plain, corrected))
    assert 0.001 < moved < 1.0                                                # MHz: a kHz-scale change
    _, default = m.hyperfine_components(32, 0, 56, "R")                      # the table is the default
    assert np.allclose([c.offset for c in corrected], [c.offset for c in default])
    assert table.correction("C", 32, 150) == table.correction("C", 32, 59)   # clamped beyond the measured J'

