"""The catalog (data/catalog/*.yaml) and its links to the data sets: each fact stored once."""

import re

from i2spec import catalog

STATES = ("not needed:", "not possible:", "to do:", "n/a:")


def test_every_set_points_to_one_entry():
    ds = catalog.datasets()          # raises if a pointer names no entry
    assert sum(map(len, ds.values())) >= 295


def test_one_entry_per_publication():
    seen = {}
    for e in catalog.entries().values():
        d = re.sub(r"^https?://(dx\.)?doi\.org/", "", str(e.get("doi") or "").lower())
        if d and d != "10.1051/0004-6361/202244441":   # three atlas spectra cite one reduction-code paper
            assert d not in seen, f"{e['id']} and {seen[d]} are the same publication"
            seen[d] = e["id"]


def test_transcription_state_matches_the_data():
    ds = catalog.datasets()
    for i, e in catalog.entries().items():
        t = e["transcription"]
        if i in ds:
            assert t is None, f"{i}: has data sets, so its transcription field must be null, not {t!r}"
        else:
            assert t and t.startswith(STATES), f"{i}: no data sets, so it needs a reason, not {t!r}"


def test_status_is_given():
    for i, e in catalog.entries().items():
        assert e["status"] == "held" or e["status"].startswith(("unobtainable", "not useful")), i
