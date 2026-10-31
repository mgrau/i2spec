"""The data the fit uses: the overlap rules (a reprint enters once, from its latest printing; a mean gives way to its
primaries; a benchmark's reprints stay out with it) and the role decision (trained sets in, row records out)."""
from i2spec.fit import data as fit_data
from i2spec.observations import load_all

HEADER = "line,component,kind,value,uncertainty,ref_line,ref_component,group,note\n"


def _set(root, name, rows, overlaps=None, exclude=None):
    d = root / name
    d.mkdir()
    meta = f'id = "{name}"\ncatalog = "test"\nsource = "test"\nretrieved = "test"\nunit = "MHz"\ntranscription = "test"\n'
    if exclude:
        meta += f'exclude = "{exclude}"\n'
    if overlaps:
        meta += "[overlaps]\n" + "".join(f'{k} = "{v}"\n' for k, v in overlaps.items())
    (d / "meta.toml").write_text(meta)
    (d / "data.csv").write_text(HEADER + "".join(f"127I2 R(56) 32-0,{c},frequency,{v},0.01,,,{g},\n"
                                                 for c, v, g in rows))


def _rows(datasets):
    return sorted((d.id, o.component, o.value) for d in datasets for o in d.observations)


def test_a_reprint_enters_from_its_latest_printing_and_a_mean_gives_way(tmp_path):
    _set(tmp_path, "smith2001a", [("a10", 100.0, ""), ("a1", 50.0, "")])
    _set(tmp_path, "list2010a", [("a10", 100.0, "")], overlaps={"smith2001a": "same"})          # a reprint
    _set(tmp_path, "list2015a", [("a10", 100.0, "")], overlaps={"smith2001a": "same"})          # a later one
    _set(tmp_path, "mean2012a", [("a1", 50.2, "")], overlaps={"smith2001a": "mean"})            # a mean of it
    # the record and the default fit view: the source's rows, the mean left out of the fit
    assert _rows(load_all(tmp_path)) == [("smith2001a", "a1", 50.0), ("smith2001a", "a10", 100.0)]
    # the fit's rule: the reprint from its latest printing, the primary over the mean
    assert _rows(load_all(tmp_path, later_printing=True)) == [("list2015a", "a10", 100.0), ("smith2001a", "a1", 50.0)]
    # a source held out (a benchmark) keeps its reprints out with it
    assert _rows(load_all(tmp_path, later_printing=True, first_printing=("smith2001a",))) == \
        _rows(load_all(tmp_path))


def test_the_role_decision_admits_its_trained_sets_and_leaves_out_its_row_records(tmp_path, monkeypatch):
    _set(tmp_path, "jones2001a", [("a10", 100.0, "")])
    _set(tmp_path, "brown2002a", [("a1", 50.0, "comb"), ("a2", 60.0, "lambdameter")], exclude="pending: test")
    _set(tmp_path, "green2003a", [("a3", 70.0, "")], exclude="pending: test")
    monkeypatch.setattr(fit_data, "roles", lambda name: ({"brown2002a": "train", "green2003a": "calibrate"},
                                                          [{"set": "brown2002a", "group": "lambdameter"}]))
    assert _rows(fit_data.fit_datasets(tmp_path)) == [("brown2002a", "a1", 50.0), ("jones2001a", "a10", 100.0)]


def test_the_fit_view_of_the_current_decision_is_the_release_s():
    """No pending set trains, the row records belong to sets that do not, and the only reprints in the fit view are
    the CIPM 2007 list's of benchmark sets: the rules leave the release's rows as they are."""
    key = lambda ds: sorted((d.id, str(o.line), o.component, o.kind, o.value) for d in ds for o in d.observations)  # noqa: E731
    assert key(fit_data.fit_datasets()) == key(load_all())
    final, records = fit_data.roles()
    assert "train" not in final.values() and all(final[r["set"]] != "train" for r in records)
