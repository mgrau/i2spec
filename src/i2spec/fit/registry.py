"""Blind tests of the model: forecasts of held-out data sets, frozen before their values are read, then scored once.

A data set transcribed with ``exclude = "pending: ..."`` in its meta.toml stays out of every fit (load_all skips
it). ``freeze_registry`` reads only the identity columns of such sets (line, component, kind, ref_line,
ref_component, group: never a value or an uncertainty), classes each row by the [overlaps] declarations, and
writes every model's forecast, the primary model's 1-sigma and the SHA256 of every input. ``score_registry``
checks those hashes, then reads the values and scores the frozen forecasts. A registry is a folder of files that
are written once and never overwritten:

  protocol.json    models, partitions, kinds, metrics and comparisons: fixed before any value is read
  targets.json     one target per data row, id "<set>:<row>": its identity columns, its status (PARTITIONS),
                   whether its line is in the training data, and the diagnostic exclusion
  forecasts.json   per model and target: predicted_MHz (None, with a reason, outside the model's domain); the
                   primary model's forecasts carry its 1-sigma, sigma_model_MHz
  snapshot.json    SHA256 of each set's data.csv (meta.toml for information), of the code and of the releases
  manifest.json    SHA256 of the four files above; written last, it marks a complete registry
  results.json     the scores, written once by score_registry

The exact line 127I2 P(13) 43-0 is kept as a diagnostic target and left out of every headline score. The first two
registries (release 20261020, validation/pending_blind_v1 and v2) were frozen by an earlier version of this code
with four comparators; score_registry reads and reproduces them.

From release 20261029 a registry also holds the pending sets of data/hyperfine_parameters (B hyperfine constants,
target id "hyperfine_parameters/<set>:<row>:<quantity>") and data/x_hyperfine (X hyperfine constants and intervals
within an X level, "x_hyperfine/<set>:<file>:<row>[:<quantity>]"), from their identity columns too, and may forecast
again the sets of earlier registries (``rescore``): their values were read when those were scored, so they form
their own partition and never enter a blind score.
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import tomllib

import numpy as np

from ..constants import MHZ_PER_CM
from ..observations import COLUMNS, DATA_DIR, UNITS, Line, _ends, _graph, _observation, _path, load_dataset

ROOT = Path(__file__).resolve().parents[3]
PACKAGE = ROOT / "src" / "i2spec"
#: the data folders beside data/observations whose pending sets a registry holds
HYPERFINE, X_HYPERFINE = "hyperfine_parameters", "x_hyperfine"
B_LEVELS, X_LEVELS = "b_levels", "x_levels"
HF_QUANTITIES = ("eqQ", "C", "d", "delta")
#: printed names of the four constants (yokozeki1980's meta.toml: C1 = C, C3 = d, C4 = delta); other printed
#: quantities (octupole, the fit's B, sigma) are not the model's and get no target
X_ALIASES = {"C1": "C", "C3": "d", "C4": "delta"}
EXCLUSION = "127I2 P(13) 43-0"
THRESHOLDS_MHZ = (.01, .03, .1, .3, 1., 3.)
FILES = ("protocol.json", "targets.json", "forecasts.json", "snapshot.json")
IDENTITY = ("line", "component", "kind", "ref_line", "ref_component", "group")
#: the modules a forecast depends on, hashed into the snapshot with the manifests of the releases used
CODE = ("research.py", "research_expanded.py", "hyperfine.py", "bspline.py", "wigner.py", "lookup.py",
        "constants.py", "observations.py", "fit/registry.py")
PARTITIONS = dict(
    blind="nothing of this measurement is in any training set or an earlier registry",
    scored_before="an [overlaps] declaration ties it to a set of an earlier registry, whose values were read when "
                  "that registry was scored",
    seen="an active set prints this measurement or a mean that contains it ([overlaps] same/mean)",
    superseded="an active set holds the same quantity re-measured or re-evaluated by the same authors ([overlaps] "
               "superseded)",
    read_before="a set of an earlier registry, whose values were read when that registry was scored; forecast "
                "again, never fitted, and never in a blind score",
    derived="a constant a source derived from its own line data (term values, rotational constants, band origins): "
            "the table is in no training set or earlier registry, but lines it was derived from may be")
KINDS = dict(absolute="a component frequency",
             intraline="component minus component of the same line (hyperfine only)",
             interline="component (or center) of one line minus that of another",
             hf_constant="an effective hyperfine constant (eqQ, C, d or delta) of a line's B level, or of B minus X "
                         "where the set's meta.toml kind says the source printed differences",
             x_constant="an effective hyperfine constant of an X level",
             x_interval="the energy of one hyperfine sublevel of an X level minus another's, E(I, F) - E(I', F'); "
                        "sublevels are labeled (I, F) by a one-to-one assignment within each F (the largest total "
                        "weight of each I on the level's own J); a set that prints magnitudes is scored on |forecast|",
             center="a line's hyperfine-free center (Doppler-limited wavenumbers)",
             level="a term value G (J = 0) or rotational constant B, D, H of an X or B level; the model's B, D, H are "
                   "the coefficients of E(J) = G + B x - D x^2 + H x^3, x = J(J + 1), fitted to its levels J = 0-20",
             band="a band origin, E_B(v', 0) - E_X(v'', 0)",
             fluorescence="a line of a fluorescence series from one upper level to X(v'', J''): forecast -E_X(v'', J'') "
                          "and scored after the series' offset (the upper level's energy) is removed")
#: the target kinds of the level and line tables (LEVEL_TABLES), read from and scored in cm-1 converted to MHz
LEVEL_KINDS = ("center", "level", "band", "fluorescence")
#: how an offset is removed: the unweighted mean residual of the group's printed rows, per model; a group of one row
#: is not scored
OFFSET = "one constant per group, the unweighted mean of its residuals per model; a group of one row is not scored"


def _sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()


def _write_once(path, value):
    """JSON to a new file; an existing file is never overwritten."""
    with open(path, "x", encoding="utf-8") as f:
        json.dump(value, f, indent=1, allow_nan=False)
        f.write("\n")


def _relative(path):
    path = Path(path).resolve()
    return str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path)


def _meta(root, set_id):
    return tomllib.loads((Path(root) / set_id / "meta.toml").read_text(encoding="utf-8"))


def pending_sets(root=DATA_DIR):
    """The ids of the sets held out as pending (``exclude = "pending: ..."``)."""
    return [p.name for p in sorted(Path(root).iterdir()) if (p / "meta.toml").exists()
            and str(_meta(root, p.name).get("exclude", "")).startswith("pending")]


def registered_sets(*registries):
    """The sets the given registries hold, sorted."""
    return sorted({t["set"] for r in registries for t in json.loads((Path(r) / "targets.json").read_text())})


def identities(set_id, root=DATA_DIR):
    """The rows of a set with placeholder values: only the identity columns of its data.csv are read."""
    with open(Path(root) / set_id / "data.csv", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        if tuple(reader.fieldnames or ()) != COLUMNS:
            raise ValueError(f"{set_id}: unexpected columns")
        return [_observation(dict({k: row[k] for k in IDENTITY}, value="0", uncertainty="1", note=""))
                for row in reader]


def classify(pending, prior=(), root=DATA_DIR):
    """({(set, row): (status, related set, relation)}, lines in the training data) of the rows of ``pending``
    ({set: identities}), from the [overlaps] between them and the active sets, and between them and the sets of
    earlier registries (``prior``: their values have been read, so 'scored_before'). The stronger tie wins."""
    active = {}
    for p in sorted(Path(root).iterdir()):
        if (p / "meta.toml").exists() and p.name not in pending and (
                "exclude" not in _meta(root, p.name) or p.name in prior):
            active[p.name] = load_dataset(p, raw=True)
    status, rank = {}, {"seen": 0, "scored_before": 1, "superseded": 2}

    def mark(key, relation, other):
        new = "superseded" if relation == "superseded" else "scored_before" if other in prior else "seen"
        if key not in status or rank[new] < rank[status[key][0]]:
            status[key] = (new, other, relation)

    def matches(spec, o):
        lines = spec.get("lines")
        return (not lines or o.line in {Line.parse(s) for s in lines}) and (not spec.get("groups")
                                                                           or o.group in spec["groups"])

    def specs(meta):
        for sid, spec in meta.get("overlaps", {}).items():
            yield sid, {"relation": spec} if isinstance(spec, str) else spec
    for d in active.values():                   # an active set prints a pending set's measurement
        for sid, spec in specs(d.meta):
            if sid in pending:
                graph = _graph(pending[sid])
                for o in d.observations:
                    if matches(spec, o):
                        for j in _path(graph, *_ends(o)) or ():
                            mark((sid, j), spec["relation"], d.id)
    for pid, obs in pending.items():             # a pending set declares an active source of its rows
        for sid, spec in specs(_meta(root, pid)):
            if sid in active:
                graph = _graph(active[sid].observations)
                for k, o in enumerate(obs):
                    if matches(spec, o) and _path(graph, *_ends(o)):
                        mark((pid, k), spec["relation"], sid)
    fitted = [d for sid, d in active.items() if sid not in prior]
    trained = {o.line for d in fitted for o in d.observations} | {
        o.ref_line for d in fitted for o in d.observations if o.ref_line}
    return status, trained


def build_targets(sets, prior=(), root=DATA_DIR):
    """The targets of the pending ``sets``, one per data row, from their identity columns."""
    pending = {s: identities(s, root) for s in sets}
    status, trained = classify(pending, set(prior), root)
    targets = []
    for s, obs in pending.items():
        for k, o in enumerate(obs):
            st = status.get((s, k), ("blind", None, None))
            kind = ("absolute" if o.kind == "frequency" else
                    "intraline" if (o.ref_line or o.line) == o.line else "interline")
            targets.append(dict(id=f"{s}:{k}", set=s, row=k, kind=kind, line=str(o.line), component=o.component,
                                ref_line=str(o.ref_line) if o.ref_line else None, ref_component=o.ref_component,
                                group=o.group, isotopologue=o.line.isotopologue, status=st[0],
                                overlap_set=st[1], overlap_relation=st[2],
                                line_in_training=bool(o.line in trained or (o.ref_line and o.ref_line in trained)),
                                diagnostic_exclusion=EXCLUSION in (str(o.line), str(o.ref_line))))
    return targets


def _identity_rows(path, columns):
    """The rows of a CSV with only ``columns`` kept: the values are never read."""
    with open(path, newline="", encoding="utf-8") as f:
        return [{k: row[k] for k in columns if k in row} for row in csv.DictReader(f)]


def _header(path):
    with open(path, newline="", encoding="utf-8") as f:
        return next(csv.reader(f))


def _hyperfine_training_lines():
    """The lines whose intervals train the hyperfine fits (release 20261021's frozen hyperfine rows)."""
    import gzip
    with gzip.open(ROOT / "data" / "training" / "20261021" / "hyperfine.csv.gz", "rt", newline="") as f:
        return {r["line"] for r in csv.DictReader(f) if r["future"] in ("False", "0", "")}


def _x_training_levels(root):
    """{(set, v)}: the X levels of the active x_hyperfine sets (bacis1980's v'' 18-88 train the X eqQ beyond the
    core)."""
    return {(p.name, int(r["v"])) for p in sorted(Path(root).iterdir()) if (p / "meta.toml").exists()
            and "exclude" not in _meta(root, p.name) for f in sorted(p.glob("*.csv"))
            if "v" in _header(f) for r in _identity_rows(f, ("v",))}


def _constant_target(tid, folder, s, file, row, kind, quantity=None, **fields):
    return dict(dict(id=tid, folder=folder, set=s, file=file, row=row, kind=kind, quantity=quantity, status="blind",
                     overlap_set=None, overlap_relation=None, isotopologue="127I2", line_in_training=False,
                     diagnostic_exclusion=False), **fields)


def hyperfine_targets(sets, root):
    """Targets of pending data/hyperfine_parameters sets (``root``): eqQ, C, d and delta of every row's line, of
    its B level, or of B minus X where the set's meta.toml kind begins "B-X". Only the line column is read."""
    trained = _hyperfine_training_lines()
    out = []
    for s in sets:
        convention = "B - X" if str(_meta(root, s).get("kind", "")).startswith("B-X") else "B"
        for k, r in enumerate(_identity_rows(Path(root) / s / "data.csv", ("line",))):
            line = Line.parse(r["line"])
            for q in HF_QUANTITIES:
                out.append(_constant_target(f"{HYPERFINE}/{s}:{k}:{q}", HYPERFINE, s, "data.csv", k, "hf_constant", q,
                                            line=str(line), isotopologue=line.isotopologue, convention=convention,
                                            line_in_training=str(line) in trained,
                                            diagnostic_exclusion=str(line) == EXCLUSION))
    return out


def x_hyperfine_targets(sets, root):
    """Targets of pending data/x_hyperfine sets (``root``), from v, J and the labels only: a constant per row of a
    file with a quantity column; eqQ and C per row of a file with eqQ and C columns; an interval per row of a file
    with F_from and F_to (I of the second sublevel: eps_to, or eps; ``interval_sign = "magnitude"`` in meta.toml
    when the source prints magnitudes). An [overlaps] entry
    ``<active set> = { relation = "...", v = [...] }`` marks those levels' constants as printed by that set."""
    trained = _x_training_levels(root)
    out = []
    for s in sets:
        tied = {v: (other, spec["relation"]) for other, spec in _meta(root, s).get("overlaps", {}).items()
                for v in spec.get("v", ())}
        for path in sorted((Path(root) / s).glob("*.csv")):
            head = _header(path)
            if "F_from" in head:
                keys = ("v", "J", "eps", "eps_from", "eps_to", "F_from", "F_to")
                rows = [(k, None, r) for k, r in enumerate(_identity_rows(path, keys))]
                kind = "x_interval"
            elif "quantity" in head:
                rows = [(k, X_ALIASES.get(r["quantity"], r["quantity"]), r)
                        for k, r in enumerate(_identity_rows(path, ("v", "J", "quantity")))]
                rows = [(k, q, r) for k, q, r in rows if q in HF_QUANTITIES]
                kind = "x_constant"
            else:
                rows = [(k, q, r) for k, r in enumerate(_identity_rows(path, ("v", "J"))) for q in ("eqQ", "C")
                        if q in head or f"{q}_kHz" in head]
                kind = "x_constant"
            for k, q, r in rows:
                v, J = int(r["v"]), int(r["J"])
                fields = dict(state="X", v=v, J=J)
                if kind == "x_interval":
                    fields.update(I_from=int(r.get("eps_from") or r["eps"]), F_from=int(r["F_from"]),
                                  I_to=int(r.get("eps_to") or r["eps"]), F_to=int(r["F_to"]),
                                  magnitude=_meta(root, s).get("interval_sign") == "magnitude")
                elif v in tied:
                    fields.update(status="seen", overlap_set=tied[v][0], overlap_relation=tied[v][1])
                fields["line_in_training"] = kind == "x_constant" and any((o, v) in trained for o, _ in tied.values())
                tid = f"{X_HYPERFINE}/{s}:{path.name}:{k}" + (f":{q}" if q else "")
                out.append(_constant_target(tid, X_HYPERFINE, s, path.name, k, kind, q, **fields))
    return out


#: The level and line tables a registry holds, by (folder, set, file). form "level": a column per quantity, as
#: (quantity, scale to cm-1, origin of G: "X min" = above the X potential minimum, "X(0,0)" = above X v = 0, J = 0,
#: "offset" = not stated by the table: scored after one offset per set and state); "state": the state, or None for a
#: state column. "center": the line from the identity columns and the measured column. "band": the origin and the
#: two levels' B. "fluorescence": R and P lines (J'' and wavenumber columns) of a series from one upper level.
LEVEL_TABLES = {
    (B_LEVELS, "barrow1973", "data.csv"): dict(form="level", state="B", columns=dict(
        G=("G", 1., "offset"), B_x100=("B", 1e-2, None), D_x1e8=("D", 1e-8, None))),
    (B_LEVELS, "danyluk1977", "data.csv"): dict(form="level", state="B", columns=dict(
        G=("G", 1., "offset"), B_x100=("B", 1e-2, None))),
    (B_LEVELS, "gerstenkorn1983", "data.csv"): dict(form="level", state="B", columns=dict(
        G=("G", 1., "offset"), B_x100=("B", 1e-2, None))),
    (B_LEVELS, "king1980", "data.csv"): dict(form="level", state="B", columns=dict(
        G=("G", 1., "offset"), B_x100=("B", 1e-2, None))),
    (B_LEVELS, "tromp1983", "data.csv"): dict(form="level", state="B", columns=dict(
        G=("G", 1., "offset"), B_x100=("B", 1e-2, None))),
    (B_LEVELS, "luc1980", "data.csv"): dict(form="level", state="B", columns=dict(
        G=("G", 1., "X(0,0)"), B_x100=("B", 1e-2, None), D_x1e8=("D", 1e-8, None), H_x1e14=("H", 1e-14, None))),
    (B_LEVELS, "luc1980", "x_levels.csv"): dict(form="level", state="X", columns=dict(
        G=("G", 1., "X(0,0)"), B=("B", 1., None), D_x1e9=("D", 1e-9, None))),
    (B_LEVELS, "steinfeld1965", "data.csv"): dict(form="level", state="B", columns=dict(
        T_xmin=("G", 1., "X min"), B_x100=("B", 1e-2, None), B_cd_x100=("B", 1e-2, None))),
    (X_LEVELS, "tellinghuisen1980", "table7_levels.csv"): dict(form="level", state="X", columns=dict(
        G=("G", 1., "X min"), G_minus_G0=("G", 1., "X(0,0)"), B_x100=("B", 1e-2, None), D_x1e9=("D", 1e-9, None),
        negH_x1e16=("H", -1e-16, None))),
    (B_LEVELS, "wei1974", "bands.csv"): dict(form="band"),
    (B_LEVELS, "steinfeld1969", "lines.csv"): dict(form="center", value="sigma"),
    (X_LEVELS, "rank1964", "green_series.csv"): dict(form="center", value="sigma"),
    (B_LEVELS, "yee1976", "excitations.csv"): dict(form="center", value="sigma_meas"),
    (B_LEVELS, "gerstenkorn1991", "lines_calc.csv"): dict(form="center", value="sigma_meas"),
    (X_LEVELS, "verma1960", "data.csv"): dict(form="fluorescence"),
    (X_LEVELS, "tellinghuisen1980", "data.csv"): dict(form="fluorescence"),
}


def _integer(text):
    """An integer identity printed perhaps in parentheses ("(45)"); None if there is none."""
    import re
    m = re.fullmatch(r"\(?\s*(\d+)\s*\)?", str(text or "").strip())
    return int(m.group(1)) if m else None


def level_targets(tables, root, trained=frozenset()):
    """Targets of the LEVEL_TABLES entries ``tables`` under ``root`` (the data folder), from identity columns only
    (isotopologue, state, v, line labels, series, J): the value columns are named, never read. ``trained``: the lines
    of the fit data (a center target's line_in_training)."""
    out = []
    for folder, s, file in tables:
        spec, path = LEVEL_TABLES[(folder, s, file)], Path(root) / folder / s / file
        head, meta_iso = _header(path), str(_meta(Path(root) / folder, s).get("isotopologue", "127I2"))
        default_iso = meta_iso if meta_iso in ("127I2", "129I2", "127I129I") else "127I2"

        def add(k, column, kind, **fields):
            out.append(_constant_target(f"{folder}/{s}:{file}:{k}:{column}", folder, s, file, k, kind,
                                        column=column, **fields))
        if spec["form"] == "level":
            keys = ("isotopologue", "v", "state")
            for k, r in enumerate(_identity_rows(path, keys)):
                v, iso = _integer(r["v"]), r.get("isotopologue") or default_iso
                state = spec["state"] or r.get("state")
                if v is None or state not in ("X", "B"):
                    continue
                for column, (q, scale, origin) in spec["columns"].items():
                    if column not in head:
                        continue
                    unc = next((u for u in (f"{column}_unc",) if u in head), None)
                    add(k, column, "level", quantity=q, state=state, v=v, isotopologue=iso, scale=scale,
                        uncertainty_column=unc, origin=origin, status="derived",
                        offset_group=f"{folder}/{s}:{file}:{state}:{iso}:G" if origin == "offset" else None)
        elif spec["form"] == "band":
            for k, r in enumerate(_identity_rows(path, ("isotopologue", "v_upper", "v_lower"))):
                vu, vl, iso = _integer(r["v_upper"]), _integer(r["v_lower"]), r.get("isotopologue") or default_iso
                common = dict(isotopologue=iso, status="derived", offset_group=None, origin=None)
                add(k, "nu0", "band", quantity="nu0", v_upper=vu, v_lower=vl, scale=1., uncertainty_column="nu0_unc",
                    **common)
                for column, state, v in (("Bu_x100", "B", vu), ("Bl_x100", "X", vl)):
                    add(k, column, "level", quantity="B", state=state, v=v, scale=1e-2,
                        uncertainty_column=f"{column}_unc", **common)
        elif spec["form"] == "center":
            keys = ("isotopologue", "line", "branch", "J_lower", "v_upper", "v_lower")
            for k, r in enumerate(_identity_rows(path, keys)):
                iso = r.get("isotopologue") or default_iso
                try:
                    if r.get("branch") and r.get("J_lower"):
                        line = Line.parse(f"{iso} {r['branch']}({int(r['J_lower'])}) {int(r['v_upper'])}-"
                                          f"{int(r['v_lower'])}")
                    else:
                        line = Line.parse(f"{iso} {r['line']}")
                except (ValueError, TypeError):
                    continue
                add(k, spec["value"], "center", component=None, ref_line=None, ref_component=None, line=str(line),
                    isotopologue=iso, scale=1., uncertainty_column=None, offset_group=None,
                    line_in_training=str(line) in trained, diagnostic_exclusion=str(line) == EXCLUSION)
        else:                               # fluorescence: R and P of a series from one upper level
            keys = ("table", "series", "J_upper", "v", "R_J", "P_J")
            for k, r in enumerate(_identity_rows(path, keys)):
                v, Ju = _integer(r["v"]), _integer(r["J_upper"])
                if v is None or Ju is None:
                    continue
                for branch in ("R", "P"):
                    J = _integer(r.get(f"{branch}_J"))
                    if J is None:
                        continue
                    add(k, branch, "fluorescence", quantity=branch, state="X", v=v, J=J, isotopologue=default_iso,
                        scale=1., uncertainty_column=None, J_upper=Ju, series=f"{r.get('table')}/{r.get('series')}",
                        offset_group=f"{folder}/{s}:{file}:{r.get('table')}/{r.get('series')}:J'={Ju}")
    return out


def _label(line):
    return f"{line.branch}({line.J_lower}) {line.v_upper}-{line.v_lower}"


class Forecaster:
    """A frozen release's forecasts of targets (``name``: a ResearchModel id or role), each isotopologue's model
    built on first use; a target outside the release's domain has no forecast."""

    def __init__(self, name):
        from ..research import ResearchModel
        self.name, self._new, self._models = name, ResearchModel, {}
        m = self.model("127I2")
        core = getattr(m, "_core_entry", m._entry)          # the frozen release a later one extends
        entries = [m._entry] + ([core] if core is not m._entry else [])
        self.info = dict(id=m.info.id, coefficients=hashlib.sha256(_canonical(m.coefficients)).hexdigest(),
                         files_sha256={e["release"]: e["files_sha256"] for e in entries})
        self.releases = [e["release"] for e in entries]

    def model(self, isotopologue):
        if isotopologue not in self._models:
            self._models[isotopologue] = self._new(self.name, isotopologue)
        return self._models[isotopologue]

    def energy(self, isotopologue, state, v, J):
        """E(v, J) of a level, cm-1 above the X potential minimum."""
        return self.model(isotopologue)._levels_x.energy(state, v, J)

    def rotational(self, isotopologue, state, v):
        """(G, B, D, H) of a level, cm-1: E(J) = G + B x - D x^2 + H x^3 fitted to the model's levels J = 0-20."""
        x = np.array([J*(J + 1) for J in range(21)], float)
        c = np.polynomial.polynomial.polyfit(x, [self.energy(isotopologue, state, v, J) for J in range(21)], 3)
        return c[0], c[1], -c[2], c[3]

    def level(self, t):
        """The forecast of a target of a level or line table (LEVEL_KINDS), cm-1."""
        iso = t["isotopologue"]
        if t["kind"] == "band":
            return self.energy(iso, "B", t["v_upper"], 0) - self.energy(iso, "X", t["v_lower"], 0)
        if t["kind"] == "fluorescence":
            return -self.energy(iso, "X", t["v"], t["J"])
        if t["quantity"] == "G":
            zero = 0. if t["origin"] == "X min" else self.energy(iso, "X", 0, 0)
            return self.energy(iso, t["state"], t["v"], 0) - zero
        return dict(zip("GBDH", self.rotational(iso, t["state"], t["v"])))[t["quantity"]]

    def predict(self, t):
        """The forecast of a target, MHz: a component (or center) frequency, or the difference of two; a
        hyperfine constant; or an interval between two hyperfine sublevels of an X level; or a quantity of a level or
        line table, converted from cm-1."""
        if t["kind"] in ("level", "band", "fluorescence"):
            return float(self.level(t))*MHZ_PER_CM
        if t["kind"] == "hf_constant":
            line = Line.parse(t["line"])
            m = self.model(line.isotopologue)
            Ju = line.J_lower + (1 if line.branch == "R" else -1)
            value = _constant(m.level_hyperfine("B", line.v_upper, Ju), t["quantity"])
            if t["convention"] == "B - X":
                value -= _constant(m.level_hyperfine("X", line.v_lower, line.J_lower), t["quantity"])
            return value
        if t["kind"] == "x_constant":
            return _constant(self.model("127I2").level_hyperfine("X", t["v"], t["J"]), t["quantity"])
        if t["kind"] == "x_interval":
            levels = _labeled(self.model("127I2").hyperfine_levels("X", t["v"], t["J"]), t["J"])
            return levels[(t["I_from"], t["F_from"])] - levels[(t["I_to"], t["F_to"])]
        line = Line.parse(t["line"])
        value = self.model(line.isotopologue).frequency(_label(line), t["component"])
        if t["kind"] not in ("absolute", "center"):
            ref = Line.parse(t["ref_line"]) if t["ref_line"] else line
            value -= self.model(ref.isotopologue).frequency(_label(ref), t["ref_component"])
        return float(value)

    def alternatives(self, t):
        """An X interval's forecast with each I the second sublevel may have ({I: MHz}): a source may print the
        second level's I only in its note (borde2002), which is read when the registry is scored."""
        levels = _labeled(self.model("127I2").hyperfine_levels("X", t["v"], t["J"]), t["J"])
        if (t["I_from"], t["F_from"]) not in levels:
            return {}
        return {str(I): levels[(t["I_from"], t["F_from"])] - e for (I, F), e in sorted(levels.items())
                if F == t["F_to"]}

    def sigma(self, t):
        """The 1-sigma of a target's forecast (sigma_model_MHz, its statistical, model and hyperfine parts, and the
        data status), or the reason there is none."""
        if t["kind"] in ("hf_constant", "x_constant", "x_interval", "level", "band", "fluorescence"):
            return dict(sigma_model_MHz=None, sigma_reason="the uncertainty model covers line positions and "
                                                           "components, not constants, X intervals or levels")
        line = Line.parse(t["line"])
        ref = Line.parse(t["ref_line"]) if t["ref_line"] else None
        try:
            if ref is not None and ref.isotopologue != line.isotopologue:
                raise ValueError("no 1-sigma of an interval between isotopologues")
            args = (_label(line), t["component"]) + ((t["ref_component"],) if t["kind"] not in ("absolute", "center")
                                                 else ()) \
                + ((_label(ref),) if t["kind"] == "interline" else ())
            u = self.model(line.isotopologue).uncertainty(*args)
        except (ValueError, KeyError) as e:
            return dict(sigma_model_MHz=None, sigma_reason=str(e)[:200])
        return dict(sigma_model_MHz=u.sigma_MHz, sigma_parts_MHz=[u.statistical_MHz, u.model_MHz, u.hyperfine_MHz],
                    level_status=u.status, hyperfine_status=u.hyperfine_status)


def _constant(p, quantity):
    """One constant of HyperfineParameters in MHz (C, d and delta are kept in kHz)."""
    return float(getattr(p, quantity)) * (1. if quantity == "eqQ" else 1e-3)


def _labeled(levels, J):
    """{(I, F): energy MHz} of the hyperfine sublevels of a level of rotation J: within each F, the eigenvectors are
    matched one to one to the I values of J by the largest total weight (the dominant-I label of strongly mixed
    sublevels may repeat)."""
    from scipy.optimize import linear_sum_assignment
    out = {}
    for F in sorted({lv.F for lv in levels}):
        block = [lv for lv in levels if lv.F == F]
        Is = sorted({b[1] for b in block[0].basis if b[0] == J})
        w = np.array([[sum(float(lv.vector[k])**2 for k, b in enumerate(lv.basis) if b[0] == J and b[1] == I)
                       for I in Is] for lv in block])
        for r, c in zip(*linear_sum_assignment(-w)):
            out[(Is[c], F)] = block[r].energy
    return out


def forecasts(targets, models):
    """{model name: Forecaster} and every model's forecast of every target, in model order; the first model's
    forecasts carry its 1-sigma where it has an uncertainty model."""
    out, built = [], {}
    for n, name in enumerate(models):
        f = built[name] = Forecaster(name)
        with_sigma = n == 0 and hasattr(f.model("127I2"), "uncertainty")
        for t in targets:
            try:
                row = dict(model=name, target=t["id"], predicted_MHz=f.predict(t))
            except (ValueError, KeyError, IndexError) as e:
                row = dict(model=name, target=t["id"], predicted_MHz=None, reason=str(e)[:200])
            if t["kind"] == "x_interval":           # also when the same-I pair has no forecast (a printed
                row["alternatives_MHz"] = f.alternatives(t)     # label may name another I)
            elif row["predicted_MHz"] is None:
                out.append(row)
                continue
            out.append(row | (f.sigma(t) if with_sigma else {}))
    return built, out


def freeze_registry(destination, sets=None, *, after=(), models=("i2spec",), root=DATA_DIR, decision=None,
                    rescore=(), constants=False, constants_root=None, disclosure=None, protocol_extra=None,
                    levels=(), held=()):
    """Freeze the forecasts of pending sets into the new folder ``destination``; returns the manifest.

    ``sets``: pending set ids (default: every pending set the registries ``after`` do not hold). ``models``: the
    releases to forecast with (ResearchModel ids or roles), the first the primary. ``rescore``: registries among
    ``after`` whose sets are forecast again, in the partition read_before. ``constants``: also every pending set of
    data/hyperfine_parameters and data/x_hyperfine (under ``constants_root``, default data/), or a collection of set
    ids: only those. ``disclosure``: values of registered sets seen before freezing, stated in the protocol.
    ``protocol_extra``: further fields of the protocol (a test's criteria, fixed before any value is read).
    ``levels``: (folder, set, file) entries of LEVEL_TABLES to register. ``held``: (folder, set, reason) of pending
    sets registered without targets (the model forecasts nothing of them); their files are hashed for the record.
    Every target and forecast is made before the
    folder is created; the manifest, written last, marks a complete registry. An existing destination, a set that
    is not pending and a set an earlier registry holds (unless rescored) are refused."""
    destination = Path(destination)
    if destination.exists():
        raise FileExistsError(f"{destination} exists: a registry is never overwritten")
    if not models or len(set(models)) != len(models):
        raise ValueError("declare unique model names before forecasting")
    earlier = [Path(r) for r in after]
    for r in earlier:
        read_registry(r)                        # an earlier registry must itself be intact
    prior = registered_sets(*earlier)
    if any(Path(r) not in earlier for r in rescore):
        raise ValueError("a rescored registry must be among the earlier ones")
    again = registered_sets(*rescore) if rescore else []
    pending = pending_sets(root)
    sets = [s for s in pending if s not in prior] if sets is None else list(sets)
    if (not sets and not constants and not levels and not held) or len(set(sets)) != len(sets):
        raise ValueError("declare a nonempty list of distinct sets")
    if bad := [s for s in sets if s not in pending or s in prior]:
        raise ValueError(f"only pending sets no earlier registry holds can be registered: {bad}")
    targets = build_targets(sets, prior, root)
    if again:
        if bad := [s for s in again if s not in pending]:
            raise ValueError(f"a rescored set must still be pending: {bad}")
        targets += [dict(t, status="read_before") for t in build_targets(again, prior, root)]
    other = {}
    if constants:
        croot = Path(constants_root) if constants_root is not None else DATA_DIR.parent
        other = {f: pending_sets(croot / f) for f in (HYPERFINE, X_HYPERFINE)}
        if constants is not True:
            named = set(constants)
            if bad := named - {s for ss in other.values() for s in ss}:
                raise ValueError(f"only pending constant sets can be registered: {sorted(bad)}")
            other = {f: [s for s in ss if s in named] for f, ss in other.items()}
        targets += hyperfine_targets(other[HYPERFINE], croot / HYPERFINE)
        targets += x_hyperfine_targets(other[X_HYPERFINE], croot / X_HYPERFINE)
    lroot = Path(constants_root) if constants_root is not None else DATA_DIR.parent
    if levels:
        if bad := [e for e in levels if tuple(e) not in LEVEL_TABLES or e[1] not in pending_sets(lroot / e[0])]:
            raise ValueError(f"only pending sets of LEVEL_TABLES can be registered: {bad}")
        from .data import fit_datasets
        trained = frozenset(str(o.line) for d in fit_datasets() for o in d.observations)
        targets += level_targets([tuple(e) for e in levels], lroot, trained)
        for f, s_, _ in levels:
            other.setdefault(f, [])
            if s_ not in other[f]:
                other[f].append(s_)
        croot = lroot
    built, rows = forecasts(targets, models)
    primary = models[0]
    releases = sorted({r for f in built.values() for r in f.releases})
    created = datetime.now(timezone.utc).isoformat()
    protocol = dict(
        schema_version=3, created_utc=created, status="frozen",
        decision=decision or "the sets of this registry stay out of every fit and are used only to score these "
                             "frozen predictions and uncertainties",
        primary_model=primary, models={name: f.info for name, f in built.items()},
        comparisons=dict(all=list(models)),
        outcome_rule="values are read only by score_registry, after this registry exists; no refit, no model or "
                     "setting selection, no exclusion of targets after scoring. A data correction found later is "
                     "reported beside, not instead of, these scores",
        exclusion=f"{EXCLUSION}: diagnostic target only, never in a headline score",
        partitions=PARTITIONS, kinds=KINDS,
        **(dict(rescored_registries=[_relative(r) for r in rescore]) if rescore else {}),
        **(dict(constants_scope="eqQ, C, d and delta of every row (a quantity a source did not print is not "
                                "scored); summaries per kind and quantity") if constants else {}),
        **(dict(disclosure=disclosure) if disclosure else {}),
        **(dict(levels_scope=dict(tables=[list(e) for e in levels], offset=OFFSET,
                                  rotational="the model's B, D, H: E(J) = G + B x - D x^2 + H x^3 fitted to its "
                                             "levels J = 0-20, x = J(J + 1); G: E(v, 0) from the origin the table states")) if levels else {}),
        **(dict(held_without_forecast=[dict(folder=f, set=s_, reason=why) for f, s_, why in held]) if held else {}),
        **(protocol_extra or {}),
        metrics=dict(residual="observed - predicted (MHz), observed at zero pressure and power",
                     scores="per model, partition, kind and set: n, rms, median |r|, mean, rms of r/sigma, fraction "
                            "within each threshold; paired on targets every compared model predicts",
                     thresholds_MHz=list(THRESHOLDS_MHZ)),
        uncertainty_note="sigma is the measurement uncertainty alone; |r|/sigma describes disagreement with the "
                         "measurement, not model confidence. The primary forecasts also carry the model's frozen "
                         "1-sigma (sigma_model_MHz); coverage uses z = r / sqrt(sigma_model^2 + sigma^2)")
    sets_all = sets + again
    snapshot = dict(data={s: {"data.csv": _sha(Path(root) / s / "data.csv")} for s in sets_all},
                    meta_informational={s: _sha(Path(root) / s / "meta.toml") for s in sets_all},
                    earlier_registries=[_relative(r) for r in earlier], earlier_sets=prior,
                    earlier_manifests_sha256={_relative(r): _sha(r / "manifest.json") for r in earlier},
                    code={_relative(PACKAGE / c): _sha(PACKAGE / c) for c in CODE} | {
                        _relative(PACKAGE / "data" / "research" / r / "manifest.json"):
                        _sha(PACKAGE / "data" / "research" / r / "manifest.json") for r in releases})
    if other:
        snapshot["data_other"] = {f"{f}/{s}": {p.name: _sha(p) for p in sorted((croot / f / s).glob("*.csv"))}
                                  for f, ss in other.items() for s in ss}
        snapshot["meta_informational_other"] = {f"{f}/{s}": _sha(croot / f / s / "meta.toml")
                                                for f, ss in other.items() for s in ss}
        snapshot["data_other_root"] = _relative(croot)
    if held:
        hroot = Path(constants_root) if constants_root is not None else DATA_DIR.parent
        snapshot["held_informational"] = {f"{f}/{s_}": {q.name: _sha(q) for q in sorted(
            (hroot / f / s_).glob("*") if (hroot / f / s_).is_dir() else (hroot / f).glob(f"{s_}*"))
            if q.is_file()} for f, s_, _ in held}
    destination.mkdir(parents=True)
    for name, value in zip(FILES, (protocol, targets, rows, snapshot)):
        _write_once(destination / name, value)
    manifest = dict(created_utc=created, files_sha256={n: _sha(destination / n) for n in FILES})
    _write_once(destination / "manifest.json", manifest)
    return manifest


def read_registry(folder):
    """(protocol, targets, forecasts, snapshot, manifest) of a frozen registry, after checking each file against
    manifest.json and that every model forecasts every target once. A folder without a manifest is incomplete."""
    folder = Path(folder)
    manifest = json.loads((folder / "manifest.json").read_text())
    if set(manifest["files_sha256"]) != set(FILES):
        raise ValueError(f"{folder.name}: unexpected or incomplete registry manifest")
    for name, digest in manifest["files_sha256"].items():
        if _sha(folder / name) != digest:
            raise ValueError(f"{folder.name}/{name} changed after freezing")
    protocol, targets, rows, snapshot = (json.loads((folder / n).read_text()) for n in FILES)
    ids = [t["id"] for t in targets]
    if protocol.get("status") != "frozen" or protocol.get("primary_model") not in protocol["models"]:
        raise ValueError(f"{folder.name}: not a frozen registry with a declared primary model")
    if len(set(ids)) != len(ids) or sorted((f["model"], f["target"]) for f in rows) != sorted(
            (m, i) for m in protocol["models"] for i in ids):
        raise ValueError(f"{folder.name}: every model must forecast every target exactly once")
    return protocol, targets, rows, snapshot, manifest


def changed_inputs(snapshot, root=DATA_DIR):
    """(data files that differ from the snapshot, sets whose meta.toml differs from its informational hash)."""
    root = Path(root)
    data = sorted(f"{s}/{n}" for s, files in snapshot["data"].items() for n, d in files.items()
                  if _sha(root / s / n) != d)
    meta = sorted(s for s, d in snapshot.get("meta_informational", {}).items() if _sha(root / s / "meta.toml") != d)
    if "data_other" in snapshot:
        other = ROOT / snapshot["data_other_root"]
        data += sorted(f"{s}/{n}" for s, files in snapshot["data_other"].items() for n, d in files.items()
                       if _sha(other / s / n) != d)
        meta += sorted(s for s, d in snapshot["meta_informational_other"].items()
                       if _sha(other / s / "meta.toml") != d)
    return data, meta


def _number(text):
    return float(text) if text not in (None, "") else None


def _observed_constant(t, other):
    """(observed MHz, sigma MHz, the second sublevel's I) of a target of data/hyperfine_parameters or
    data/x_hyperfine, read from its file under ``other`` (the data folder); (None, None, None) when the source
    did not print that quantity. A printed label "a,b -> c,d" in an X interval's note names the second I (c)."""
    import re
    with open(other / t["folder"] / t["set"] / t["file"], newline="", encoding="utf-8") as f:
        row = list(csv.DictReader(f))[t["row"]]
    if t["kind"] in LEVEL_KINDS:             # cm-1, scaled as printed
        value = _number(row[t["column"]])
        sigma = _number(row[t["uncertainty_column"]]) if t["uncertainty_column"] in row else None
        k = t["scale"]*MHZ_PER_CM
        return (value*k if value is not None else None), (abs(sigma*k) if sigma is not None else None), None
    if t["kind"] == "hf_constant":
        return _number(row[t["quantity"]]), _number(row[f"{t['quantity']}_unc"]), None
    if t["kind"] == "x_interval":
        label = re.search(r"printed (\d+),(\d+) -> (\d+),(\d+)", row.get("note", ""))
        I_to = int(label.group(3)) if label else t["I_to"]
        value, sigma = _number(row["value_kHz"]), _number(row["uncertainty_kHz"])
        return (value*1e-3 if value is not None else None), (sigma*1e-3 if sigma is not None else None), I_to
    if "quantity" in row:
        scale = {"MHz": 1., "kHz": 1e-3}[row["unit"]]
        value, sigma = _number(row["value"]), _number(row["uncertainty"])
    else:                                   # eqQ in MHz, C in kHz (C_kHz)
        key, unc, scale = ("eqQ", "eqQ_unc", 1.) if t["quantity"] == "eqQ" else ("C_kHz", "C_unc_kHz", 1e-3)
        value, sigma = _number(row[key]), _number(row[unc])
    return (value*scale if value is not None else None), (sigma*scale if sigma is not None else None), None


def _comparisons(protocol):
    """The model groups a registry pairs its scores on. The first two registries name none: they compare every
    model, and the primary against the comparators that are not research releases."""
    if "comparisons" in protocol:
        return protocol["comparisons"]
    models = protocol["models"]
    return dict(all_four=list(models), primary_vs_production=[protocol["primary_model"]] + [
        m for m, v in models.items() if not v["id"].startswith("research-")])


def _summary(r, z, thresholds):
    r, z = np.asarray(r), np.asarray(z, dtype=float)
    if not len(r):
        return None
    if not np.all(np.isfinite(z)):          # rows without a measurement sigma
        z = z[np.isfinite(z)]
    return dict(n=len(r), rms_MHz=float(np.sqrt(np.mean(r*r))), median_abs_MHz=float(np.median(np.abs(r))),
                mean_MHz=float(np.mean(r)), rms_z=float(np.sqrt(np.mean(z*z))) if len(z) else None,
                within_MHz={str(t): float(np.mean(np.abs(r) <= t)) for t in thresholds})


def score_registry(folder, out=None, *, root=DATA_DIR, strict=True):
    """Score a frozen registry against the measurements and write the results once; returns them.

    The registry's files are checked against its manifest, and the data files against its snapshot: a changed
    data file stops scoring with ``strict``; otherwise the scores are made and the changes listed in
    inputs_changed_since_freeze. A meta.toml hashed for information is only listed (meta_changed_since_freeze): a
    later [overlaps] declaration does not change a measurement. ``out`` defaults to the registry's results.json,
    which must not exist: a registry is scored once. Only the registry's own sets are read."""
    folder = Path(folder)
    out = Path(out) if out is not None else folder / "results.json"
    if out.exists():
        raise FileExistsError(f"{out} exists: a registry is scored once")
    protocol, targets, rows, snapshot, manifest = read_registry(folder)
    changed, meta_changed = changed_inputs(snapshot, root)
    if changed and strict:
        raise ValueError(f"{folder.name}: inputs changed after freezing: {changed}")
    models, primary = list(protocol["models"]), protocol["primary_model"]
    thresholds = protocol["metrics"]["thresholds_MHz"]
    observed = {}
    for s in snapshot["data"]:
        d = load_dataset(Path(root) / s)
        for k, o in enumerate(d.observations):
            observed[f"{s}:{k}"] = o, o.value*UNITS[d.unit], o.uncertainty*UNITS[d.unit]
    pred = {(f["model"], f["target"]): f["predicted_MHz"] for f in rows}
    alternatives = {(f["model"], f["target"]): f["alternatives_MHz"] for f in rows if "alternatives_MHz" in f}
    sig = {f["target"]: f.get("sigma_model_MHz") for f in rows if f["model"] == primary}
    other = ROOT / snapshot["data_other_root"] if "data_other_root" in snapshot else None
    scored, not_printed = [], []
    for t in targets:
        if "folder" in t:                       # a hyperfine constant, an X interval or a level or line table
            value, sigma, I_to = _observed_constant(t, other)
            if value is None:
                not_printed.append(t["id"])
                continue
            forecast = {}
            for m in models:
                p = pred[(m, t["id"])]
                if (m, t["id"]) in alternatives:
                    p = alternatives[(m, t["id"])].get(str(I_to))
                    p = abs(p) if t["magnitude"] and p is not None else p
                forecast[m] = p
            scored.append(dict(t, observed_MHz=value, sigma_MHz=sigma, I_to_printed=I_to, residual_MHz={
                m: value-p if p is not None else None for m, p in forecast.items()}))
            continue
        o, value, sigma = observed[t["id"]]
        if (str(o.line), o.component, str(o.ref_line) if o.ref_line else None, o.ref_component) != (
                t["line"], t["component"], t["ref_line"], t["ref_component"]):
            raise ValueError(f"{t['id']}: the data row is not the registered target")
        scored.append(dict(t, observed_MHz=value, sigma_MHz=sigma, residual_MHz={
            m: value-pred[(m, t["id"])] if pred[(m, t["id"])] is not None else None for m in models}))
    offsets = {}
    for r in scored:                            # OFFSET: one constant per group removed
        if r.get("offset_group"):
            offsets.setdefault(r["offset_group"], []).append(r)
    single = []
    for rs in offsets.values():
        for m in models:
            got = [r["residual_MHz"][m] for r in rs if r["residual_MHz"][m] is not None]
            c = float(np.mean(got)) if len(got) > 1 else None
            for r in rs:
                r.setdefault("raw_residual_MHz", {})[m] = r["residual_MHz"][m]
                r.setdefault("offset_MHz", {})[m] = c
                r["residual_MHz"][m] = r["residual_MHz"][m] - c if c is not None and r["residual_MHz"][m] is not None \
                    else None
        if all(r["offset_MHz"][primary] is None for r in rs):
            single += [r["id"] for r in rs]
    headline = [r for r in scored if not r["diagnostic_exclusion"] and r["id"] not in set(single)]

    def scores(subset, names):
        common = [r for r in subset if all(r["residual_MHz"][m] is not None for m in names)]
        return dict(n_targets=len(common), sets=sorted({r["set"] for r in common}),
                    models={m: _summary([r["residual_MHz"][m] for r in common],
                                        [r["residual_MHz"][m]/r["sigma_MHz"] if r["sigma_MHz"] else np.nan
                                         for r in common], thresholds)
                            for m in names})
    results = dict(scored_utc=datetime.now(timezone.utc).isoformat(), registry_manifest=manifest, summary={},
                   by_set={}, inputs_changed_since_freeze=changed, meta_changed_since_freeze=meta_changed,
                   rows=scored)
    if not_printed:
        results["not_printed"] = not_printed
    if single:
        results["offset_group_of_one"] = single
    if any(v is not None for v in sig.values()):          # blind coverage of the frozen 1-sigma
        cov = {}
        for r in headline:
            res, s_m = r["residual_MHz"][primary], sig.get(r["id"])
            if res is None or s_m is None:
                continue
            r["z_coverage"] = res/float(np.hypot(s_m, r["sigma_MHz"] or 0.))   # no printed sigma: the model's alone
            for key in (r["status"], f"{r['status']} | {r['kind']}", f"{r['status']} | {r['set']}"):
                cov.setdefault(key, []).append(r["z_coverage"])
        results["coverage"] = {k: dict(n=len(z), within_1sigma=float(np.mean(np.abs(z) <= 1)),
                                       within_2sigma=float(np.mean(np.abs(z) <= 2)),
                                       rms_z=float(np.sqrt(np.mean(np.square(z))))) for k, z in sorted(cov.items())}
    def label(r):                               # constants are scored per quantity
        return r["kind"] + (f" {r['quantity']}" if r.get("quantity") else "")

    def group(r):
        return f"{r['folder']}/{r['set']}" + (f" | {r['quantity']}" if r.get("quantity") else "") \
            if "folder" in r else r["set"]
    for partition in PARTITIONS:
        for kind in sorted({label(r) for r in headline}, key=lambda k: (list(KINDS).index(k.split()[0]), k)):
            sub = [r for r in headline if r["status"] == partition and label(r) == kind]
            for name, comparison in _comparisons(protocol).items():
                s = scores(sub, comparison)
                if s["n_targets"]:
                    results["summary"][f"{partition} | {kind} | {name}"] = s
    for s in sorted({group(r) for r in headline}):
        sub = [r for r in headline if group(r) == s]
        avail = [m for m in models if all(r["residual_MHz"][m] is not None for r in sub)]
        results["by_set"][s] = dict(status=sorted({r["status"] for r in sub}), kinds=sorted({r["kind"] for r in sub}),
                                    **scores(sub, avail))
    _write_once(out, results)
    return results


def main(argv=None):
    parser = argparse.ArgumentParser(description="Freeze or score a blind-test registry of the pending data sets.")
    commands = parser.add_subparsers(dest="command", required=True)
    freeze = commands.add_parser("freeze", help="freeze the forecasts of the pending sets no earlier registry holds")
    freeze.add_argument("destination", type=Path)
    freeze.add_argument("--after", nargs="*", type=Path, default=[], help="earlier registries")
    freeze.add_argument("--model", nargs="+", default=["i2spec"], help="releases to forecast with, primary first")
    freeze.add_argument("--rescore", nargs="*", type=Path, default=[], help="earlier registries to forecast again")
    freeze.add_argument("--constants", action="store_true", help="also the pending hyperfine constants and X intervals")
    freeze.add_argument("--disclosure", help="values of registered sets seen before freezing")
    score = commands.add_parser("score", help="score a frozen registry, once")
    score.add_argument("registry", type=Path)
    args = parser.parse_args(argv)
    if args.command == "freeze":
        print(json.dumps(freeze_registry(args.destination, after=args.after, models=tuple(args.model),
                                         rescore=args.rescore, constants=args.constants,
                                         disclosure=args.disclosure), indent=1))
    else:
        result = score_registry(args.registry)
        for key, s in result["summary"].items():
            print(f"{key:52s} n={s['n_targets']:4d}  " + "  ".join(
                f"{m[:12]:12s} rms {v['rms_MHz']:9.4f} med {v['median_abs_MHz']:8.4f}" for m, v in s["models"].items()))


if __name__ == "__main__":
    main()
