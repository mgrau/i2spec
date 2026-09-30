"""The NIR level corrections of i2spec2026b: loaded, applied, and worth their claim on the data."""
import numpy as np
import pytest

from i2spec.constants import MHZ_PER_CM
from i2spec.level_corrections import load_level_corrections
from i2spec.model import RovibronicModel
from i2spec.observations import DATA_DIR, Predictor, load_dataset, residuals


def test_corrections_load_and_apply():
    c = load_level_corrections("level_corrections_2026c")
    assert c.shift("X", 13, 100) != 0.0 and c.shift("X", 3, 100) == 0.0
    assert c.covers("X", 14, 100) and not c.covers("X", 14, 260)
    assert c.shift("B", 39, 33) == pytest.approx(5.57, abs=0.1)           # nishiyama2024a: a pure offset
    assert c.uncertainty("B", 39) < 0.3 and c.uncertainty("B", 43) is None   # one line fixed v' = 43
    a, b = RovibronicModel("127I2", "i2spec2026a"), RovibronicModel("127I2", "i2spec2026d")
    d = (b.energy("X", 13, 100) - a.energy("X", 13, 100)) * MHZ_PER_CM
    assert d == pytest.approx(c.shift("X", 13, 100), abs=1e-6)
    assert b.energy("X", 3, 50) == a.energy("X", 3, 50)                   # untouched elsewhere
    assert RovibronicModel("129I2", "i2spec2026d").corrections is None    # 127I2 only


def test_extended_range_potentials():
    """i2spec2026d: levels from v'' = 18 and v' = 44 up come from the MLR pair, continuous at the seam."""
    a, d = RovibronicModel("127I2", "i2spec2026a"), RovibronicModel("127I2", "i2spec2026d")
    assert d.extended is not None and d.from_v == {"X": 18, "B": 44}
    assert d.energy("X", 5, 50) == a.energy("X", 5, 50)                                # published below the seam
    assert abs(d.energy("X", 18, 50) - a.energy("X", 18, 50)) < 0.01                  # continuous across it
    assert (d.energy("X", 48, 60) - a.energy("X", 48, 60)) < -5.0                     # the published X fails there
    assert (d.energy("B", 58, 99) - a.energy("B", 58, 99)) > 1.0                      # and the published B
    assert RovibronicModel("129I2", "i2spec2026d").energy("X", 48, 60) > 0            # isotopologues carry the BO terms
    assert RovibronicModel("127I2", "i2spec2026d", potentials={"X": a.potentials["X"]}).extended is None


@pytest.mark.parametrize("name, before, after", [("liao2010a", 2.30, 0.45), ("bodermann2000a", 3.13, 0.35),
                                                 ("bodermann1998b", 3.3, 0.1), ("nishiyama2024a", 1.3, 0.35),
                                                 ("bipm2012a", 0.96, 0.45), ("yoshiki2023a", 30.0, 1.0),
                                                 ("matsunaga2024a", 5000.0, 3.0)])
def test_nir_sets_reproduced(name, before, after):
    ds = load_dataset(DATA_DIR / name)
    scale = MHZ_PER_CM if ds.unit == "cm-1" else 1.0
    rms = lambda p: float(np.sqrt(np.mean(np.square(residuals(p, ds) * scale))))   # noqa: E731
    assert rms(Predictor("i2spec2026a")) > before * 0.9
    assert rms(Predictor("i2spec2026d")) < after


def test_band_correction_is_confined_to_the_nir_bands():
    """i2spec2026e adds a line correction to v' = 0 -> v'' = 12-17 only, inside their measured J'' range;
    every other line is exactly i2spec2026d's."""
    d, e = RovibronicModel("127I2", "i2spec2026d"), RovibronicModel("127I2", "i2spec2026e")
    for vu, vl, J, br in ((0, 9, 100, "P"), (0, 10, 200, "R"), (0, 11, 150, "P"), (1, 13, 100, "R"),
                          (2, 12, 60, "P"), (5, 13, 26, "R"), (0, 18, 100, "P"), (0, 17, 150, "R"),
                          (0, 15, 250, "P"), (32, 0, 60, "R")):
        assert e.transition(vu, vl, J, br) == d.transition(vu, vl, J, br), (vu, vl, J, br)
    moved = abs(e.transition(0, 13, 106, "R") - d.transition(0, 13, 106, "R")) * 29979.2458
    assert 0.2 < moved < 2.0                                  # MHz: the R(106) 0-13 residual it exists for


def test_band_correction_improves_liao_and_leaves_other_sets_alone():
    import numpy as np
    from i2spec.observations import load_dataset, residuals
    rms = lambda ps, sid: float(np.sqrt(np.mean(residuals(Predictor(ps), load_dataset(f"data/observations/{sid}")) ** 2)))  # noqa: E731
    assert rms("i2spec2026e", "liao2010a") < 0.16 < 0.3 < rms("i2spec2026d", "liao2010a")
    for sid in ("cornish2000a", "reinhardt2007a", "bipm2003d"):
        assert rms("i2spec2026e", sid) == rms("i2spec2026d", sid)


def test_line_list_carries_the_band_correction():
    """master_line_list takes positions from levels; the band correction has to reach it separately."""
    from i2spec.intensity import intensity_model, line_list, position_model
    m = intensity_model()
    nu = position_model(m).transition(0, 13, 106, "R")
    ll = line_list(m, 800.0, nu - 0.02, nu + 0.02, S_min=1e-32)
    k = [i for i in range(len(ll.nu)) if (ll.v_upper[i], ll.v_lower[i], ll.J_lower[i], ll.branch[i]) == (0, 13, 106, 1)]
    assert len(k) == 1 and abs(ll.nu[k[0]] - nu) * 29979.2458 < 0.01


def test_orsay_levels_are_confined_to_the_atlas_range():
    """i2spec2026g corrects X v'' = 18-25 from the Orsay atlas and leaves every measured data set untouched."""
    import numpy as np
    from i2spec.observations import load_all
    f, g = RovibronicModel("127I2", "i2spec2026f"), RovibronicModel("127I2", "i2spec2026g")
    for v in (5, 12, 16, 17, 26, 30, 48):
        assert g.energy("X", v, 100) == f.energy("X", v, 100), v
    d20 = (g.energy("X", 20, 100) - f.energy("X", 20, 100)) * MHZ_PER_CM
    assert -40 < d20 < -20                     # the atlas puts X v'' = 20 ~30 MHz below the MLR
    F, G = Predictor("i2spec2026f"), Predictor("i2spec2026g")
    for ds in load_all():
        assert np.allclose(residuals(F, ds), residuals(G, ds), atol=1e-6), ds.id


def test_orsay_part1_transcription_is_complete():
    import csv
    rows = list(csv.DictReader(open(DATA_DIR.parent / "atlas_lines/orsay1982_part1.csv")))
    assert [int(r["N"]) for r in rows] == list(range(1, 4827))          # Table I: 4 826 lines in part I
    s = [float(r["sigma_cm1"]) for r in rows]
    assert all(b > a for a, b in zip(s, s[1:])) and 11000 < s[0] < s[-1] < 13010


def test_atlas_refit_potential_keeps_measured_levels():
    """i2spec2026h: mlr_x_2026d carries the atlas levels into v'' = 26-47; inside the atlas coverage and
    everywhere else measured, levels stay what i2spec2026g gave."""
    g, h = RovibronicModel("127I2", "i2spec2026g"), RovibronicModel("127I2", "i2spec2026h")
    for v, J in ((5, 60), (17, 42), (18, 100), (22, 150), (24, 100)):
        assert abs(h.energy("X", v, J) - g.energy("X", v, J)) * MHZ_PER_CM < 2.0, (v, J)
    d34 = (h.energy("X", 34, 60) - g.energy("X", 34, 60)) * MHZ_PER_CM
    assert -150 < d34 < -50                      # the refit moves the unmeasured middle by ~-100 MHz


def test_correction_curvature_is_held_outside_coverage():
    """level_corrections_2026h: beyond coverage the linear part extrapolates (velchev1998a B v'=16 needs it),
    the curvature terms do not; inside coverage nothing changes."""
    h, i = load_level_corrections("level_corrections_2026g"), load_level_corrections("level_corrections_2026h")
    for st, v, J in (("B", 16, 60), ("X", 14, 100), ("X", 13, 150)):
        assert i.shift(st, v, J) == pytest.approx(h.shift(st, v, J), abs=1e-9)
    c = i.coefficients["B"][16]
    assert len(c) == 2 and i.shift("B", 16, 138) == h.shift("B", 16, 138)       # linear: still extrapolates
    lo, hi = i.coverage["X"][14]
    assert len(i.coefficients["X"][14]) > 2 and i.shift("X", 14, 300) != h.shift("X", 14, 300)


def test_levels_to_the_dissociation_limit():
    """B levels above the default grid come from the 12 A grid, and agree with the default where both reach."""
    m = RovibronicModel("127I2")
    assert m.energy("B", 72, 26) < m.energy("B", 0, 0) + 20043.3 - m.energy("X", 0, 0) + m.energy("X", 0, 0)
    e62 = m.energy("B", 62, 26)
    assert m.energy("B", 61, 26) < e62 < m.energy("B", 63, 26)
    far = m._far_levels("B", 26, 62)
    assert abs(far[59] - m.levels("B", 26)[59]) * MHZ_PER_CM < 0.1
