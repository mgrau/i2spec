"""Line widths, pressure coefficients and level lifetimes (data/widths/, data/lifetimes/, i2spec.lineshape)."""

import pytest

from i2spec.lineshape import (LIFETIME_COLUMNS, WIDTH_COLUMNS, load_all_lifetimes, load_all_widths, load_lifetimes,
                              load_widths)

META = 'id = "{id}"\ncatalog = "broadband:tench1983"\nsource = "s"\nretrieved = "2026-10-28"\ntranscription = "t"\nuncertainty = "u"\n'


def _set(tmp_path, name, columns, rows, extra=""):
    d = tmp_path / name
    d.mkdir()
    (d / "meta.toml").write_text(META.format(id=name) + extra)
    (d / "data.csv").write_text(",".join(columns) + "\n" + "".join(r + "\n" for r in rows))
    return d


def test_widths_round_trip(tmp_path):
    d = _set(tmp_path, "w1", WIDTH_COLUMNS, ['127I2 P(13) 43-0,a3,width,64,6,,295,,w1,"printed 32 kHz HWHM"',
                                             "127I2 R(127) 11-5,,broadening,253,25,,,self,w1,"])
    s = load_widths(d)
    assert [r.unit for r in s.rows] == ["kHz", "kHz/Pa"]
    assert s.rows[0].pressure_Pa is None and s.rows[0].temperature_K == 295 and s.rows[0].component == "a3"
    assert s.rows[1].component is None and s.rows[1].perturber == "self"


def test_lifetimes_round_trip(tmp_path):
    d = _set(tmp_path, "l1", LIFETIME_COLUMNS, ["127I2,B,11,128,,lifetime,410,15,l1,",
                                                "127I2,B,43,12,F=17,hyperfine_predissociation_a2,3.36,0.09,l1,"])
    s = load_lifetimes(d)
    assert [r.unit for r in s.rows] == ["ns", "s-1"] and s.rows[1].sublevel == "F=17"


@pytest.mark.parametrize("row, message", [
    ("127I2 P(13) 43-0,a3,fwhm,64,6,,,,w,", "quantity"),
    ("127I2 P(13) 43-0,a3,width,64,0,,,,w,", "uncertainty"),
    ("127I2 Q(13) 43-0,a3,width,64,6,,,,w,", "bad line"),
])
def test_widths_rejects_bad_rows(tmp_path, row, message):
    with pytest.raises(ValueError, match=message):
        load_widths(_set(tmp_path, "w", WIDTH_COLUMNS, [row]))


def test_pending_sets_are_held_out(tmp_path):
    _set(tmp_path, "a", LIFETIME_COLUMNS, ["127I2,B,6,32,,lifetime,310,10,a,"])
    _set(tmp_path, "b", LIFETIME_COLUMNS, ["127I2,B,6,32,,lifetime,310,10,b,"], extra='exclude = "pending"\n')
    assert [s.id for s in load_all_lifetimes(tmp_path)] == ["a"]
    assert len(load_all_lifetimes(tmp_path, include_excluded=True)) == 2


def test_every_set_in_the_repository_loads():
    load_all_widths(include_excluded=True)
    load_all_lifetimes(include_excluded=True)
