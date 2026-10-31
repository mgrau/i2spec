"""The measurements the fit uses: observations.load_all's fit view, with the rules for overlapping data and the role
decision of the pending sets applied.

  - a measurement printed twice (a compilation's "same") enters once, from its latest printing; a mean gives way to
    the primaries it combines ("primary is superior to mean"); a superseded measurement leaves for its successor;
  - a reprint of a set held out as the core fit's benchmark (rows.NEW_SETS) stays out with it;
  - the sets the role decision (data/training/roles_<date>/outcome.json, final_roles) trains on are admitted, and
    its row-level records (row_records: a group of rows of a set kept as the record only) are left out.

With the decision of 2026-10-29 no pending set trains, and these rules give the release's own rows unchanged
(tests/test_fit_data.py).
"""
from __future__ import annotations

import json
from dataclasses import replace

from ..observations import DATA_DIR, load_all
from .rows import NEW_SETS, TRAINING

#: the role decision the fit follows
ROLES = "roles_20261029"


def roles(name=ROLES):
    """(final_roles: set id -> role, row_records: [{set, group, reason}]) of a role decision."""
    outcome = json.loads((TRAINING / name / "outcome.json").read_text())
    return outcome["final_roles"], outcome.get("row_records", [])


def fit_datasets(root=DATA_DIR, raw=False, decision=ROLES):
    """The data sets of the fit (see the module docstring), each with only the rows the fit uses."""
    final, row_records = roles(decision) if decision else ({}, [])
    train = sorted(s for s, role in final.items() if role == "train")
    records = {(r["set"], r["group"]) for r in row_records}
    out = []
    for d in load_all(root, raw=raw, include=train, later_printing=True, first_printing=NEW_SETS):
        kept = [o for o in d.observations if (d.id, o.group) not in records]
        if kept:
            out.append(replace(d, observations=kept) if len(kept) != len(d.observations) else d)
    return out
