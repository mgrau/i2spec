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


def _determined(observations, path, a, b):
    """f(a) - f(b) from a source's rows along a path that joins the two nodes."""
    from i2spec.observations import _ends

    total, node = 0.0, a
    for k in path:
        o = observations[k]
        x, y = _ends(o)                                  # o.value = f(x) - f(y)
        total, node = (total + o.value, y) if node == x else (total - o.value, x)
    assert node == b
    return total


def test_overlaps_name_sets_and_relations():
    from i2spec.observations import RELATIONS

    names = {p.name for p in DATA_DIR.iterdir() if (p / "meta.toml").exists()}
    declared = 0
    for ds in load_all(include_excluded=True):
        assert "defer_to" not in ds.meta, f"{ds.id}: defer_to is now [overlaps]"
        for sid, spec in ds.meta.get("overlaps", {}).items():
            assert sid in names and sid != ds.id, (ds.id, sid)
            assert (spec if isinstance(spec, str) else spec["relation"]) in RELATIONS
            declared += 1
    assert declared >= 17


def test_a_reprint_is_its_source_printed_again():
    """A compilation row merged into a source's record ("same") must be the source's value, re-referenced and
    rounded at most: within 1 kHz or half the compilation's uncertainty. A larger difference means a mean, a
    misprint or a wrong label."""
    from i2spec.observations import _ends, _graph, _overlaps, _path

    datasets = [load_dataset(p, raw=True) for p in sorted(DATA_DIR.iterdir()) if (p / "meta.toml").exists()]
    by_id = {d.id: d for d in datasets}
    _, merged, _ = _overlaps(datasets)
    assert len(merged) >= 200
    def admits(spec, o):
        spec = {"relation": spec} if isinstance(spec, str) else spec
        return (spec["relation"] != "superseded" and o.group in spec.get("groups", [o.group])
                and str(o.line) in spec.get("lines", [str(o.line)]))

    for cid, k in merged:
        o = by_id[cid].observations[k]
        a, b = _ends(o)
        paths = [(s, p) for s, spec in by_id[cid].meta["overlaps"].items() if admits(spec, o)
                 for p in [_path(_graph(by_id[s].observations), a, b)] if p]
        assert len(paths) == 1, (cid, k)
        sid, path = paths[0]
        value = _determined(by_id[sid].observations, path, a, b)
        assert abs(o.value - value) <= max(1e-3, o.uncertainty / 2) + 1e-9, (cid, str(o.line), o.component, sid)


def test_one_record_per_measurement():
    """The record keeps every row of every set once: a compilation's reprints only as citations on their
    source's rows, cited after it. The fit uses part of the record."""
    record, fit = load_all(record=True), load_all()
    ids = {ds.id for ds in record}
    for ds in record:
        for o in ds.observations:
            assert o.citations[0] == ds.id and len(set(o.citations)) == len(o.citations)
            assert set(o.citations[1:]) <= ids
    rows = lambda sets: {(ds.id, o) for ds in sets for o in ds.observations}   # noqa: E731
    assert rows(fit) <= rows(record)
    by_id = {ds.id: ds for ds in record}
    # MEP 532 nm: f(a10) is the CIPM's mean of the absolute measurements, so it is a record of its own and
    # stays out of the fit; Holzwarth's a10 is one of them, and Hong & Ishikawa 2000 is reprinted as Tables 5-6
    a10 = [o for o in by_id["bipm2012a"].observations if o.kind == "frequency"]
    assert [o.value for o in a10] == [563260223.513] and not any(o.kind == "frequency" for o in
                                                                 next(d for d in fit if d.id == "bipm2012a").observations)
    holzwarth = next(o for o in by_id["holzwarth2001a"].observations if o.component == "a10")
    assert holzwarth.citations == ("holzwarth2001a", "bipm2012a")
    assert all(o.citations == ("hong2000a", "bipm2012a") for o in by_id["hong2000a"].observations
               if o.kind == "interval" and o.ref_line is None)
    # Barwood's P(62) 17-1 table (MEP 576 nm) is superseded by Reinhardt's: kept, not fitted, not cited on his
    assert sum(o.kind == "interval" for o in by_id["bipm2003c"].observations) == 10
    assert all(o.citations == ("reinhardt2006a",) for o in by_id["reinhardt2006a"].observations)


def test_a_set_left_out_has_no_say():
    """borde1981a is pending: until it is admitted, the BIPM 515 nm tables it is averaged into stay in the fit."""
    count = lambda sets: len(next(d for d in sets if d.id == "bipm2005a").observations)   # noqa: E731
    assert count(load_all()) == 56
    assert count(load_all(include_excluded=True)) == 56 - 32


def test_a_paper_that_also_prints_others_keeps_its_own():
    """picard2003a prints hong2001b and arie1993a beside its own lasers, on the same R(56) 32-0 intervals: the
    "same" entries name those columns' groups, so only they merge and the lasers' rows stay records."""
    sets = {d.id: d for d in load_all(include_excluded=True, record=True)}
    groups = {o.group for o in sets["picard2003a"].observations}
    assert "picard2003a_mikes" in groups and "picard2003a_bipm" in groups
    assert not groups & {"picard2003a_jilanist2", "picard2003a_arie"}
    r56 = [o for o in sets["hong2001b"].observations if str(o.line) == "127I2 R(56) 32-0"]
    assert r56 and all("picard2003a" in o.citations for o in r56)


def test_a_quantity_measured_twice_is_cited_twice():
    """ye2000a measures f(a16) of R(127) 11-5 by two chains; the MEP 633 nm value averages both, so both cite it."""
    sets = {d.id: d for d in load_all(include_excluded=True, record=True)}
    assert [o.citations for o in sets["ye2000a"].observations] == [("ye2000a", "bipm2003a")] * 2


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
