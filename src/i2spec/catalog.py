"""The bibliography of iodine measurements (data/catalog/*.yaml) and its links to the data sets in data/.

Each fact is stored once. A catalog entry holds what a source is (citation, DOI, what it measures, whether we
hold it); a data set's meta.toml points to its entry with ``catalog = "<file>:<id>"``. The data sets of an
entry are found from those pointers, never stored in the entry.
"""

from __future__ import annotations

import tomllib
from collections import defaultdict
from functools import cache
from pathlib import Path

import yaml

DATA = Path(__file__).resolve().parents[2] / "data"
CATALOGS = ("precision", "broadband")


@cache
def entries() -> dict[str, dict]:
    """Every catalog entry by id, with the file it lives in as ``entry["catalog"]``."""
    out = {}
    for name in CATALOGS:
        for e in yaml.safe_load((DATA / "catalog" / f"{name}.yaml").read_text(encoding="utf-8")):
            if e["id"] in out:
                raise ValueError(f"catalog id {e['id']!r} appears twice")
            out[e["id"]] = {**e, "catalog": name}
    return out


def _set_metas():
    """(path of the set, its meta) for every data set: folders with meta.toml, and the atlas .toml files."""
    for meta in sorted(DATA.glob("*/*/meta.toml")):
        if meta.parts[-3] not in ("catalog", "external"):
            yield meta.parent, tomllib.loads(meta.read_text(encoding="utf-8"))
    for toml in sorted((DATA / "atlas_lines").glob("*.toml")):
        yield toml, tomllib.loads(toml.read_text(encoding="utf-8"))


def entry_of(meta: dict) -> dict:
    """The catalog entry a data set's meta.toml points to."""
    name, _, eid = meta["catalog"].partition(":")
    e = entries().get(eid)
    if e is None or e["catalog"] != name:
        raise KeyError(f"{meta.get('id')}: catalog {meta['catalog']!r} names no entry")
    return e


@cache
def datasets() -> dict[str, list[str]]:
    """Catalog id -> the data sets transcribed from it (paths relative to data/)."""
    out = defaultdict(list)
    for path, meta in _set_metas():
        out[entry_of(meta)["id"]].append(str(path.relative_to(DATA)))
    return dict(out)
