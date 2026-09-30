"""Observation data sets: format, validation and model residuals (docs/design/observations.md)."""

import numpy as np
import pytest

from i2spec.observations import (COLUMNS, DATA_DIR, Line, Predictor, component_rank, load_all, load_dataset,
                                 read_observations, residuals, write_observations)


def test_line_keys():
    line = Line.parse("129I2 P(69) 12-6")
    assert (line.isotopologue, line.branch, line.J_lower, line.v_upper, line.v_lower) == ("129I2", "P", 69, 12, 6)
    assert str(line) == "129I2 P(69) 12-6"
    for bad in ("I2 R(56) 32-0", "127I2 Q(56) 32-0", "127I2 P(0) 5-0", "127I2 R56 32-0"):
        with pytest.raises(ValueError):
            Line.parse(bad)


def test_component_rank():
    assert component_rank("a10") == 10 and component_rank("m48") == 48
    with pytest.raises(ValueError):
        component_rank("f")


def test_all_datasets_load():
    datasets = load_all()
    assert datasets and all(len(ds) > 0 for ds in datasets)
    assert len({ds.id for ds in datasets}) == len(datasets)


def test_excluded_datasets_skipped():
    everything = load_all(include_excluded=True)
    excluded = {ds.id for ds in everything if "exclude" in ds.meta}
    assert "chen2023a" in excluded
    assert not excluded & {ds.id for ds in load_all()}


def test_round_trip(tmp_path):
    ds = load_dataset(DATA_DIR / "bipm2012a")
    write_observations(ds.observations, tmp_path / "data.csv")
    assert read_observations(tmp_path / "data.csv") == ds.observations


@pytest.mark.parametrize("row, message", [
    ("127I2 R(56) 32-0,a10,width,1.0,0.1,,,,", "kind"),
    ("127I2 R(56) 32-0,a10,interval,1.0,0.1,,,,", "needs ref"),
    ("127I2 R(56) 32-0,a10,frequency,1.0,0,,,,", "positive"),
    ("127I2 R(56) 32-0,a10,frequency,1.0,0.1,,a1,,", "no reference"),
    ("127I2 R(56) 32,a10,frequency,1.0,0.1,,,,", "bad line"),
    ("127I2 R(56) 32-0,f,frequency,1.0,0.1,,,,", "component label"),
])
def test_validation(tmp_path, row, message):
    (tmp_path / "data.csv").write_text(",".join(COLUMNS) + "\n" + row + "\n")
    with pytest.raises(ValueError, match=f"line 2: .*{message}"):
        read_observations(tmp_path / "data.csv")


def test_columns_checked(tmp_path):
    (tmp_path / "data.csv").write_text("line,value\n127I2 R(56) 32-0,1.0\n")
    with pytest.raises(ValueError, match="columns"):
        read_observations(tmp_path / "data.csv")


def test_bipm_532_residuals():
    """The published parameters: f(a10) 2.0 MHz high, B-level errors of ±6 MHz over v' = 32-37."""
    ds = load_dataset(DATA_DIR / "bipm2012a")
    pairs = list(zip(ds.observations, residuals(Predictor(), ds)))
    absolute = [r for o, r in pairs if o.kind == "frequency" and str(o.line) == "127I2 R(56) 32-0"]
    hyperfine = [r for o, r in pairs if o.kind == "interval" and o.ref_line is None
                 and str(o.line) in ("127I2 R(56) 32-0", "127I2 R(87) 33-0")]
    between_lines = [r for o, r in pairs if o.kind == "interval" and o.ref_line is not None]
    assert abs(absolute[0]) < 3.0  # MHz
    assert len(hyperfine) == 32 and np.abs(hyperfine).max() < 0.15
    assert len(between_lines) == 19 and np.abs(between_lines).max() < 8.0
