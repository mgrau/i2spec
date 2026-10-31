"""The hyperfine data of the fits: intervals between components of one line, and published effective B constants.

  intervals   every 127I2 interval between two labelled components of one line (P(13) 43-0 excluded); `inside`
              marks the core (B v'<=43, X v''<=17), `future` the core's benchmark sets
  constants   Chen et al. 2004's effective B constants (v' 42-70; their 1g-flagged lines and P(13) 43-0 left out;
              each with a floor, their errors being statistical only) and Pique et al. 1986's v'-averaged constants
              at v' 71 and 75. Units of hyperfine.HyperfineParameters: eqQ MHz; C, d, delta kHz.

The release's own rows are frozen in data/training/<release>/hyperfine.csv.gz (`release`).
"""
from __future__ import annotations

import csv
import gzip
import re
from math import sqrt

from ..constants import MHZ_PER_CM
from .data import fit_datasets
from ..observations import Line, component_rank
from . import rows as training

EXCLUDED = "127I2 P(13) 43-0"
CORE_V = dict(X=17, B=43)
CHEN_PERTURBED = {"R(65) 57-0", "R(74) 57-0", "P(17) 59-0", "P(27) 59-0", "R(28) 59-0"}
CHEN_FLOOR = dict(eqQ=.2, C=1., d=3., delta=3.)
PIQUE = {71: dict(eqQ=(-555., 2.), C=(2160., 70.), d=(-900., 100.), delta=(3560., 40.)),
         75: dict(eqQ=(-537., 5.), C=(3650., 180.), d=(-2700., 200.), delta=(5680., 80.))}
PARAMETERS = ("eqQ", "C", "d", "delta")


def intervals(datasets=None):
    rows = []
    for ds in datasets if datasets is not None else fit_datasets():
        scale = MHZ_PER_CM if ds.unit == "cm-1" else 1.
        for index, obs in enumerate(ds.observations):
            line = obs.line
            if line.isotopologue != "127I2" or str(line) == EXCLUDED:
                continue
            if obs.kind != "interval" or (obs.ref_line or line) != line or not obs.component or not obs.ref_component:
                continue
            inside = line.v_upper <= CORE_V["B"] and line.v_lower <= CORE_V["X"]
            rows.append(dict(id=f"{ds.id}:{index}", source=ds.id, line=line, rank=component_rank(obs.component),
                             ref_rank=component_rank(obs.ref_component), observed=obs.value*scale,
                             sigma=obs.uncertainty*scale, inside=inside, future=ds.id in training.NEW_SETS))
    return rows


def chen_rows(path=None):
    path = path or training.ROOT / "data" / "hyperfine_parameters" / "chen2004a" / "data.csv"
    out = []
    for r in csv.DictReader(open(path, newline="", encoding="utf-8")):
        label = r["line"].replace("127I2 ", "")
        if r["line"] == EXCLUDED or label in CHEN_PERTURBED:
            continue
        m = re.match(r"([PR])\((\d+)\) (\d+)-(\d+)", label)
        branch, J, vu = m[1], int(m[2]), int(m[3])
        Jp = J+(1 if branch == "R" else -1)
        for p, unit in (("eqQ", 1.), ("C", 1e3), ("d", 1e3), ("delta", 1e3)):
            sigma = sqrt((float(r[f"{p}_unc"])*unit)**2+CHEN_FLOOR[p]**2)
            out.append(dict(source="chen2004a", line=label, v=vu, J=Jp, parameter=p, observed=float(r[p])*unit,
                            sigma=sigma))
    return out


def pique_rows():
    return [dict(source="pique1986", line=f"v={v}", v=v, J=5, parameter=p, observed=x, sigma=s)
            for v, d in PIQUE.items() for p, (x, s) in d.items()]


FIELDS = ("kind", "id", "source", "line", "rank", "ref_rank", "v", "J", "parameter", "observed", "sigma", "inside",
          "future")


def write(path, interval_rows, constant_rows):
    with gzip.open(path, "wt", newline="", encoding="utf-8") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(FIELDS)
        for r in interval_rows:
            w.writerow(["interval", r["id"], r["source"], str(r["line"]), r["rank"], r["ref_rank"], "", "", "",
                        repr(float(r["observed"])), repr(float(r["sigma"])), int(r["inside"]), int(r["future"])])
        for r in constant_rows:
            w.writerow(["constant", "", r["source"], r["line"], "", "", r["v"], r["J"], r["parameter"],
                        repr(float(r["observed"])), repr(float(r["sigma"])), "", ""])


def release(name=training.RELEASE):
    """(intervals, constants) as the release's hyperfine fits read them."""
    with gzip.open(training.TRAINING / name / "hyperfine.csv.gz", "rt", newline="", encoding="utf-8") as f:
        records = list(csv.DictReader(f))
    out, constants = [], []
    for r in records:
        if r["kind"] == "interval":
            out.append(dict(id=r["id"], source=r["source"], line=Line.parse(r["line"]), rank=int(r["rank"]),
                            ref_rank=int(r["ref_rank"]), observed=float(r["observed"]), sigma=float(r["sigma"]),
                            inside=bool(int(r["inside"])), future=bool(int(r["future"]))))
        else:
            constants.append(dict(source=r["source"], line=r["line"], v=int(r["v"]), J=int(r["J"]),
                                  parameter=r["parameter"], observed=float(r["observed"]), sigma=float(r["sigma"])))
    return out, constants
