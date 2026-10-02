"""The HITRAN export: 160-character records that read back to the model's lines (hitran.py)."""
import math

import pytest

from i2spec import hitran
from i2spec.lookup import Catalog


@pytest.fixture(scope="module")
def catalog():
    return Catalog()


def parse(line):
    """The fields of a HITRAN 2004 record that the export fills."""
    return dict(M=int(line[0:2]), I=int(line[2]), nu=float(line[3:15]), S=float(line[15:25]), A=float(line[25:35]),
                gamma_air=float(line[35:40]), gamma_self=float(line[40:45]), E=float(line[45:55]),
                gq_up=line[67:82], gq_low=line[82:97], lq_up=line[97:112], lq_low=line[112:127],
                ierr=line[127:133], g_up=float(line[146:153]), g_low=float(line[153:160]))


def test_records_are_160_characters_and_read_back(catalog):
    lines = hitran.records(catalog, "127I2", 18788.30, 18788.40, hyperfine=False)
    r56 = next(r for r in lines if (r.v_upper, r.v_lower, r.J_lower, r.branch) == (32, 0, 56, "R"))
    text = hitran.format_record(r56, molecule=0)
    assert len(text) == 160
    f = parse(text)
    assert (f["M"], f["I"]) == (0, 1)
    assert f["nu"] == pytest.approx(18788.335525, abs=2e-6)
    assert f["S"] == pytest.approx(r56.S, rel=1e-3)
    assert f["gq_up"].split() == ["B", "32"] and f["gq_low"].split() == ["X", "0"]
    assert f["lq_low"].split() == ["R", "56"]
    # g = (2J+1) g_ns: J'' = 56 is even, so I is even (0, 2, 4) and g_ns = 15 for 127I2; J' = 57
    assert (f["g_up"], f["g_low"]) == (115 * 15, 113 * 15)
    assert f["ierr"][:2] == "52"         # 0.38 MHz is 1.3e-5 cm-1, HITRAN's code 5 (1e-5 to 1e-4); S code 2


def test_einstein_a_inverts_the_strength(catalog):
    lines = hitran.records(catalog, "127I2", 18788.30, 18788.40, hyperfine=False)
    Q = catalog.master("127I2").partition_function(hitran.T_REF)
    for r in lines:
        S = r.A * r.g_upper * math.exp(-hitran.C2 * r.E_lower / hitran.T_REF) * -math.expm1(-hitran.C2 * r.nu / hitran.T_REF) \
            / (8 * math.pi * hitran.C_CM * r.nu ** 2 * Q)
        assert S == pytest.approx(r.S, rel=1e-12)


def test_hyperfine_components_keep_the_line_strength(catalog):
    whole = hitran.records(catalog, "127I2", 18788.30, 18788.40, hyperfine=False)
    split = hitran.records(catalog, "127I2", 18788.30, 18788.40, hyperfine=True)
    key = lambda r: (r.v_upper, r.v_lower, r.J_lower, r.branch)  # noqa: E731
    for line in whole:
        comps = [r for r in split if key(r) == key(line)]
        assert len(comps) in (15, 21)
        assert sum(r.S for r in comps) == pytest.approx(line.S, rel=1e-12)
        # 2F + 1 over the components of one line: F'' runs over J'' - I .. J'' + I for each I
        assert all(r.g_lower == 2 * r.F_lower + 1 and r.g_upper == 2 * r.F_upper + 1 for r in comps)


def test_fixed_width_fields():
    assert hitran._fortran_f(0.0, 5, 4) == ".0000" and hitran._fortran_f(0.071, 5, 4) == ".0710"
    assert hitran._fortran_f(-0.0012, 8, 6) == "-.001200"
    assert hitran._err_code(1.3e-5) == 5 and hitran._err_code(9e-6) == 6 and hitran._err_code(2.0) == 0
