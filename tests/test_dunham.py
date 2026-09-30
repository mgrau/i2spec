"""Gerstenkorn & Luc 1985 as transcribed: against its own Table V and against the Hannover levels."""
import pytest

from i2spec import dunham
from i2spec.model import RovibronicModel


def test_reproduces_table_v():
    # Table V rows (E_v, B_v, D_v, -H_v, -L_v, -M_v) for v' = 0, 43, 80 and E_v'', D_v'' for v'' = 19
    assert dunham.energy("B", 0, 0) == pytest.approx(15724.5871, abs=1e-4)
    assert dunham.energy("B", 43, 0) == pytest.approx(19433.6763, abs=2e-4)
    assert dunham.energy("B", 80, 0) == pytest.approx(20041.5795, abs=2e-4)
    assert dunham.B("B", 43) == pytest.approx(0.187449932e-1, rel=1e-8)
    D, H, L, M = dunham.cdc("B", 43)
    assert D == pytest.approx(0.261854525e-7, rel=1e-6) and -H == pytest.approx(0.102372922e-12, rel=1e-6)
    assert -L == pytest.approx(0.671072385e-18, rel=1e-5) and -M == pytest.approx(0.140018365e-22, rel=1e-5)
    assert dunham.energy("X", 19, 0) == pytest.approx(0.3834293183969065e4, abs=2e-4)
    assert dunham.cdc("X", 19)[0] == pytest.approx(0.5196987661294230e-8, rel=1e-8)


def test_agrees_with_hannover_where_both_are_fitted():
    """The atlas fit and the potential fit describe the same levels to the atlas accuracy at low J
    (0.002 cm-1); at high J they part smoothly, -0.004 cm-1 at J = 120 and -0.01 at J = 200, the
    Hutson centrifugal constants of the Dunham fit against the potential's own."""
    m = RovibronicModel("127I2", "hannover2008")
    for v in (0, 5, 10, 17):
        for J in (0, 50, 100):
            assert abs(dunham.energy("X", v, J) - (m.energy("X", v, J) - m.energy("X", 0, 0))) < 0.005
    for v in (0, 10, 30, 43):
        for J in (0, 50, 100):
            assert abs(dunham.energy("B", v, J) - (m.energy("B", v, J) - m.energy("X", 0, 0))) < 0.005
    assert abs(dunham.energy("X", 0, 200) - (m.energy("X", 0, 200) - m.energy("X", 0, 0))) < 0.02
    # and above v' = 43 the published potential leaves the atlas: 0.6 cm-1 too low at v' = 50
    assert dunham.energy("B", 50, 0) - (m.energy("B", 50, 0) - m.energy("X", 0, 0)) > 0.5
