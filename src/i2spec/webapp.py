"""The web app: its data, the static site, and a local server.

    i2spec web                 build the site into site/ and serve it at http://localhost:8777
    i2spec web export          rewrite web/data/ from the current model (a few minutes; line lists are cached)
    i2spec web build DIR       the static site GitHub Pages publishes: the explorer (web/ and its data) at
                               DIR, the documentation (documentation/, MkDocs) at DIR/docs

The app never solves a Schrödinger equation. It reads line positions and temperature-independent
strengths from web/data/ and does the rest in the browser (web/physics.js), so re-exporting is how the
app picks up a model change. The export also carries every measurement we hold of each line, with its
source, so the app can mark measured lines and link each one to the paper it comes from.
"""

from __future__ import annotations

import csv
import functools
import http.server
import json
import re
import shutil
import subprocess
import time
import tomllib
import webbrowser
from datetime import date
from pathlib import Path

import numpy as np

from .constants import DEFAULT_PARAMETERS, MHZ_PER_CM

ROOT = Path(__file__).resolve().parents[2]
WEB = ROOT / "web"
DATA = WEB / "data"
ATLAS_DIR = ROOT / "data" / "atlas_lines"
OBSERVATIONS_DIR = ROOT / "data" / "observations"
#: What the site is made of; everything else in web/ is tooling (tests, reference values).
SITE_FILES = ("index.html", "app.js", "physics.js", "site.css")
MKDOCS = ROOT / "mkdocs.yml"
LABELS = {"127I2": "¹²⁷I₂", "129I2": "¹²⁹I₂", "127I129I": "¹²⁷I¹²⁹I"}
#: Hyperfine patterns kept per wavenumber shard, on top of the globally strongest ones
PER_SHARD_HFS = 80
#: Atlas line lists we hold with assignments. The Orsay tables carry their own columns (<name>_assigned.csv);
#: the other two are in the observations format. Each Orsay scale is the one its fit found: part I in
#: level_corrections_2026f, Partie IV in level_corrections_2026j (prototypes/orsay4_fit.py). Partie IV shows
#: its unblended lines only: a blend's one position belongs to none of its assignments.
ATLASES = ("salami_ross_2005", "apo_nist_2009", "orsay1982_part1", "orsay1983_part4")
ORSAY_SCALE_PPB = {"orsay1982_part1": 23.85, "orsay1983_part4": 200.8}   # Partie IV: +119.6 MHz at 19 870 cm-1
#: How the app colours each of lookup.Line.flags
FLAG_CLASS = {"corrected levels": "good", "local NIR model": "good", "v″ 18-25: atlas-measured": "warn",
              "v′ > 43: extended model": "warn", "v″ ≥ 48: extended model": "warn", "beyond 815 nm": "warn",
              "v″ > 17: unmeasured": "bad", "v″ 26-89: Martin 1986": "warn", "v′ 51-79: atlas-measured": "warn", "v′ > 50: extrapolated": "bad", "isotope shift": "none", "few data": "none"}
SHORT_NAMES = {"bipm": "BIPM", "salami_ross": "Salami & Ross", "apo_nist": "NIST FTS (APO cell)",
               "orsay": "Gerstenkorn, Vergès & Chevillard"}
#: Names that the prefix alone gets wrong
SHORT_IDS = {"orsay1983_part4": "Gerstenkorn & Luc 1983"}


# --- sources and measurements -----------------------------------------------------------------------

def _short(source_id):
    """"velchev1998a" -> "Velchev 1998", "salami_ross_2005" -> "Salami & Ross 2005"."""
    if source_id in SHORT_IDS:
        return SHORT_IDS[source_id]
    m = re.match(r"^([a-z_]+?)_?(\d{4})([a-z]?)(?:_part\d)?$", source_id)
    if not m:
        return source_id
    name = SHORT_NAMES.get(m[1].rstrip("_"), m[1].replace("_", " ").title())
    return f"{name} {m[2]}"


def _link(meta):
    doi = meta.get("doi")
    if doi:
        return f"https://doi.org/{doi}"
    for key in ("url", "source"):
        text = meta.get(key, "")
        found = re.search(r"https?://\S+", text)
        if found:
            return found[0].rstrip(").,")
    return None


def _one_line(text):
    return " ".join(str(text).split())


def sources():
    """Every data set of measured line positions the model has used, with what it covers.

    Computed from the data files alone, so the references page can be built without the model.
    """
    from .observations import UNITS, load_all

    out = []
    for ds in load_all(OBSERVATIONS_DIR):
        scale = UNITS[ds.unit]
        lines = {str(o.line) for o in ds.observations} | {str(o.ref_line) for o in ds.observations if o.ref_line}
        freq = [o for o in ds.observations if o.kind == "frequency"]
        short = _short(ds.id)
        band = re.search(r"λ ≈ (\d+) nm", ds.meta["citation"])
        if ds.id.startswith("bipm") and band:          # six BIPM documents, one per wavelength
            short = f"BIPM {band[1]} nm"
        out.append(dict(
            id=ds.id, kind="precision", short=short, citation=_one_line(ds.meta["citation"]),
            doi=ds.meta.get("doi"), link=_link(ds.meta), n_rows=len(ds), n_lines=len(lines),
            n_absolute=len(freq), n_intervals=len(ds) - len(freq),
            isotopologues=sorted({line.split()[0] for line in lines}),
            nu=[round(min(o.value * scale for o in freq) / MHZ_PER_CM, 3),
                round(max(o.value * scale for o in freq) / MHZ_PER_CM, 3)] if freq else None,
            u_MHz=float(np.median([o.uncertainty * scale for o in ds.observations])),
            uncertainty=_one_line(ds.meta.get("uncertainty", ""))))
    for name in ATLASES:
        meta = tomllib.loads((ATLAS_DIR / f"{name}.toml").read_text(encoding="utf-8"))
        rows = _atlas_rows(name)
        nu = [r[4] for r in rows]
        out.append(dict(
            id=name, kind="atlas", short=_short(name), citation=_one_line(meta["citation"]), doi=meta.get("doi"),
            link=_link(meta), n_rows=len(rows), n_lines=len({r[:4] for r in rows}), n_absolute=len(rows),
            n_intervals=0, isotopologues=["127I2"], nu=[round(min(nu), 3), round(max(nu), 3)],
            u_MHz=float(np.median([r[5] for r in rows])), uncertainty=_one_line(meta.get("calibration", ""))))
    # two sources by one author in one year keep their letters: "Bodermann 1998b", "Bodermann 1998c"
    for s in out:
        twins = [t for t in out if t["short"] == s["short"]]
        if len(twins) > 1:
            for t in twins:
                t["short"] += re.match(r".*\d{4}([a-z]?)", t["id"])[1]
    return out


@functools.lru_cache(maxsize=None)
def _atlas_rows(name):
    """(branch, J'', v', v'', ν cm⁻¹ on the corrected scale, u MHz) for each assigned atlas line.

    An Orsay line without a printed precision gets 150 MHz, as in its fit (prototypes/orsay_fit.py)."""
    from .observations import read_observations

    if name.startswith("orsay"):
        scale = 1 + ORSAY_SCALE_PPB[name] * 1e-9
        with open(ATLAS_DIR / f"{name}_assigned.csv", newline="", encoding="utf-8") as f:
            return tuple((r["branch"], int(r["J_lower"]), int(r["v_upper"]), int(r["v_lower"]),
                          float(r["sigma_cm1"]) * scale, (float(r["eps_mk"]) * 1e-3 * MHZ_PER_CM if r["eps_mk"] else 150.0))
                         for r in csv.DictReader(f) if r.get("n_assignments", "1") == "1")
    # an atlas on a wavemeter scale carries its fitted offset (offset_mhz in its .toml), subtracted here
    meta = tomllib.loads((ATLAS_DIR / f"{name}.toml").read_text(encoding="utf-8"))
    shift = float(meta.get("offset_mhz", 0.0)) / MHZ_PER_CM
    return tuple((o.line.branch, o.line.J_lower, o.line.v_upper, o.line.v_lower, o.value - shift,
                  o.uncertainty * MHZ_PER_CM) for o in read_observations(ATLAS_DIR / f"{name}.csv"))


def line_key(branch, J, v_upper, v_lower):
    """The key the app uses for a line, shared with the hyperfine file: "32-0R56"."""
    return f"{v_upper}-{v_lower}{branch}{J}"


def measurements(model_frequency):
    """Every measurement of every line, grouped by isotopologue and line.

    ``model_frequency(iso, branch, J, v', v'', component)`` gives the model's value in MHz (component None
    for the hyperfine-free line centre, else a rank) or None when it cannot, so each measurement carries
    its residual. Returns {iso: {key: [{"s": source index, "f": [[component, value, u, obs - model], ...],
    "i": number of hyperfine intervals}, ...]}}.
    """
    from .observations import UNITS, component_rank, load_all

    index = {s["id"]: k for k, s in enumerate(sources())}
    out: dict[str, dict[str, dict[int, dict]]] = {}

    def entry(line_iso, key, source):
        return out.setdefault(line_iso, {}).setdefault(key, {}).setdefault(source, {"s": source, "f": [], "i": 0})

    for ds in load_all(OBSERVATIONS_DIR):
        scale, src = UNITS[ds.unit], index[ds.id]
        for o in ds.observations:
            ln = o.line
            e = entry(ln.isotopologue, line_key(ln.branch, ln.J_lower, ln.v_upper, ln.v_lower), src)
            if o.kind == "interval":
                e["i"] += 1
                if o.ref_line and o.ref_line != ln:
                    r = o.ref_line
                    entry(r.isotopologue, line_key(r.branch, r.J_lower, r.v_upper, r.v_lower), src)["i"] += 1
                continue
            value = o.value * scale
            rank = component_rank(o.component) if o.component else None
            model = model_frequency(ln.isotopologue, ln.branch, ln.J_lower, ln.v_upper, ln.v_lower, rank)
            e["f"].append([o.component or "", round(value, 4), round(o.uncertainty * scale, 4),
                           None if model is None else round(value - model, 4)])
    for name in ATLASES:
        src = index[name]
        for branch, J, vu, vl, nu, u in _atlas_rows(name):
            model = model_frequency("127I2", branch, J, vu, vl, None)
            value = nu * MHZ_PER_CM
            entry("127I2", line_key(branch, J, vu, vl), src)["f"].append(
                ["", round(value, 1), round(u, 1), None if model is None else round(value - model, 1)])
    return {iso: {key: list(by_source.values()) for key, by_source in lines.items()} for iso, lines in out.items()}


# --- hyperfine patterns for every line -------------------------------------------------------------------

_HFS_MODEL = None


def _hfs_init(iso):
    global _HFS_MODEL
    from .model import RovibronicModel
    _HFS_MODEL = RovibronicModel(iso)


#: Weak (ΔF ≠ ΔJ) components kept for crossovers, as a fraction of the line's main components: a crossover
#: with a main component is then 2 x 1e-4 / (1/15), 0.3% of a Lamb dip, at weak saturation
WEAK_FRACTION = 1e-4


def weak_links(comps, main, total):
    """The weak components of a line that can make crossovers: [offset (MHz), strength / total, upper, lower].

    upper and lower number the hyperfine eigenstates: main component k joins upper level k and lower level k
    (main components never share a level), and a level no main component uses gets the next number after
    them. Together with that rule for the main components, the numbers identify every level the kept
    components share, so the explorer finds the same crossovers as saturation.resonances (physics.resonances).
    """
    level = lambda c, side: (getattr(c, f"{side}_level", None) if getattr(c, f"{side}_level", None) is not None  # noqa: E731
                             else (getattr(c, f"I_{side}"), getattr(c, f"F_{side}")))
    number = {side: {level(c, side): k for k, c in enumerate(main)} for side in ("upper", "lower")}
    chosen = {id(c) for c in main}
    offset = lambda c: c.offset if hasattr(c, "offset") else c.offset_MHz   # noqa: E731 (model or lookup)
    out = []
    for c in sorted(comps, key=offset):
        if id(c) in chosen or c.strength < WEAK_FRACTION * total:
            continue
        ids = [number[side].setdefault(level(c, side), max(number[side].values(), default=-1) + 1)
               for side in ("upper", "lower")]
        out.append([round(float(offset(c)), 1), round(float(c.strength / total), 6), *ids])
    return out


def _hfs_chunk(items):
    """Main components of each line, ΔJ = 0: offsets (0.1 MHz) and strengths (per mille), by frequency, and
    the weak components that share a level with one of them (weak_links), or None when there are none."""
    out = []
    for branch, J, vu, vl in items:
        try:
            _, comps = _HFS_MODEL.hyperfine_components(vu, vl, J, branch, dJ=0, position=False)
        except Exception:                 # a level beyond every validated grid
            out.append(None)
            continue
        main = sorted((c for c in comps if c.label), key=lambda c: c.offset)
        total = sum(c.strength for c in main) or 1.0
        out.append(([round(c.offset, 1) for c in main], [round(1000 * c.strength / total) for c in main],
                    weak_links(comps, main, total) or None))
    return out


def hyperfine_patterns(iso, lines, workers=None):
    """ΔJ = 0 hyperfine patterns for a list of (branch, J'', v', v'') lines, computed in parallel.

    Leaving out the J ± 2 couplings changes the offsets by under 1 MHz (0.9 MHz for R(56) 32-0), far below
    the 400 MHz Doppler width the explorer draws them under, and makes each pattern twenty times cheaper,
    which is what lets every line have one. The detail pane keeps the full calculation where it has it.
    """
    import multiprocessing as mp
    import os
    from .intensity import _THREAD_VARS, _n_workers

    lines = list(lines)
    n = _n_workers() if workers is None else max(1, workers)
    # neighbouring lines share levels, so chunks of consecutive lines reuse each worker's level cache
    order = sorted(range(len(lines)), key=lambda k: (lines[k][2], lines[k][3], lines[k][1]))
    chunks = [[lines[k] for k in order[i:i + 400]] for i in range(0, len(order), 400)]
    main = __import__("__main__")
    if n > 1 and getattr(main, "__file__", None):
        saved = {k: os.environ.get(k) for k in _THREAD_VARS}
        os.environ.update({k: "1" for k in _THREAD_VARS})
        try:
            with mp.get_context("spawn").Pool(n, initializer=_hfs_init, initargs=(iso,)) as pool:
                results = pool.map(_hfs_chunk, chunks)
        finally:
            for k, v in saved.items():
                os.environ.pop(k, None) if v is None else os.environ.__setitem__(k, v)
    else:
        _hfs_init(iso)
        results = [_hfs_chunk(c) for c in chunks]
    flat = [r for part in results for r in part]
    out = [None] * len(lines)
    for k, r in zip(order, flat):
        out[k] = r
    return out


# --- the export --------------------------------------------------------------------------------------

def _git_revision():
    try:
        return subprocess.run(["git", "-C", str(ROOT), "rev-parse", "--short", "HEAD"],
                              capture_output=True, text=True, timeout=10).stdout.strip() or None
    except Exception:
        return None


def export(out=DATA, s_min=1e-24, n_shards=64, n_hfs=1500, isotopologues=("127I2", "129I2", "127I129I")):
    """Write the app's data: line shards, partition tables, hyperfine patterns, measurements, references.

    A line weaker than ``s_min`` (cm, at 300 K) is left out unless it has been measured: every measured
    line the model has is exported, so the app can show all of them.
    """
    from .intensity import C2, intensity_model, master_line_list, nuclear_spin_weight, with_dissociation_lines
    from .lookup import Catalog, Line
    from . import __version__

    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    catalog = Catalog()
    # the lines to the B limit (12 A grid, ~20 min to build) for 127I2, the one isotopologue measured there
    masters = {iso: master_line_list(intensity_model(iso), 11000.0, 20100.0, T_range=(200.0, 600.0), S_min=1e-27)
               for iso in isotopologues}
    if "127I2" in masters:
        masters["127I2"] = with_dissociation_lines(masters["127I2"], S_min=1e-27)
    centres = {}
    for iso, master in masters.items():
        centres[iso] = {line_key("R" if b > 0 else "P", int(j), int(u), int(v)): float(nu) * MHZ_PER_CM
                        for b, j, u, v, nu in zip(master.branch, master.J_lower, master.v_upper, master.v_lower,
                                                  master.nu)}

    @functools.lru_cache(maxsize=None)
    def components(iso, branch, J, vu, vl):
        line = Line(iso, branch, J, vu, vl, 0.0, 0.0, 0.0, 300.0)
        return {int(c.label[1:]): c.frequency_MHz for c in catalog.components(line) if c.label}

    def model_frequency(iso, branch, J, vu, vl, rank):
        try:
            if rank is not None:
                return components(iso, branch, J, vu, vl).get(rank)
            key = line_key(branch, J, vu, vl)
            if key in centres.get(iso, {}):
                return centres[iso][key]
            return catalog.model(iso).transition(vu, vl, J, branch) * MHZ_PER_CM
        except Exception:            # a level beyond every validated grid
            return None

    print("  measurements and their residuals ...", flush=True)
    measured = measurements(model_frequency)
    refs = sources()
    missing = {iso: sorted(set(lines) - set(centres.get(iso, {}))) for iso, lines in measured.items()}
    (out / "references.json").write_text(json.dumps(
        {"generated": date.today().isoformat(), "sources": refs,
         "outside_absorption_list": {iso: len(keys) for iso, keys in missing.items()}}, separators=(",", ":")))

    manifest = {
        "generated": date.today().isoformat(),
        "model": {"package": "i2spec", "version": __version__, "parameters": DEFAULT_PARAMETERS, "git": _git_revision(),
                  "note": "Regenerate with `uv run i2spec web export` after a model change."},
        "constants": {"c2": C2, "mhz_per_cm": MHZ_PER_CM},
        "cutoff_cm_at_300K": s_min,
        "references": "references.json",
        "isotopologues": {},
    }
    reasons: dict = {}
    for iso in isotopologues:
        manifest["isotopologues"][iso] = _export_isotopologue(
            out, iso, masters[iso], s_min, n_shards, n_hfs, measured.get(iso, {}), catalog, reasons,
            C2, nuclear_spin_weight, Line)
    # each line's `w` indexes this list: why its uncertainty is what it is, and the flags that go with it
    manifest["reasons"] = [{"why": why, "flags": [[f, FLAG_CLASS.get(f, "none")] for f in flags]}
                           for why, flags in reasons]
    (out / "manifest.json").write_text(json.dumps(manifest, indent=1))
    total = sum(f.stat().st_size for f in out.rglob("*.json"))
    print(f"wrote {len(list(out.rglob('*.json')))} files, {total / 1e6:.1f} MB, to {out}")
    for iso, keys in missing.items():
        if keys:
            print(f"  {iso}: {len(keys)} measured lines lie outside the absorption line list (listed with zero strength), "
                  f"e.g. {', '.join(keys[:4])}")


def _export_isotopologue(out, iso, master, s_min, n_shards, n_hfs, measured, catalog, reasons, C2, nuclear_spin_weight,
                         Line):
    from .lookup import uncertainty

    lines = master.at(300.0, S_min=0.0)
    keys = np.array([line_key("R" if b > 0 else "P", int(j), int(u), int(v))
                     for b, j, u, v in zip(lines.branch, lines.J_lower, lines.v_upper, lines.v_lower)])
    is_measured = np.array([k in measured for k in keys], dtype=bool)
    keep = (lines.S >= s_min) | is_measured       # every measured line, however weak at 300 K
    cols = dict(nu=lines.nu[keep], S=lines.S[keep], el=lines.E_lower[keep], vu=lines.v_upper[keep],
                vl=lines.v_lower[keep], jl=lines.J_lower[keep], br=lines.branch[keep])
    # Measured lines the absorption list cannot hold -- emission lines from v'' = 48-54, far below its
    # range and with no population at any cell temperature -- are added with zero strength, so the app
    # still lists them and marks where they are.
    listed = set(keys.tolist())
    extra = [k for k in measured if k not in listed]
    if extra:
        model = catalog.model(iso)
        e00 = model.energy("X", 0, 0)
        rows = []
        for k in extra:
            vu_, rest = k.split("-")
            vl_, br_, J_ = re.match(r"(\d+)([PR])(\d+)", rest).groups()
            vu_, vl_, J_ = int(vu_), int(vl_), int(J_)
            try:
                nu_ = model.transition(vu_, vl_, J_, br_)
                rows.append((nu_, 0.0, model.energy("X", vl_, J_) - e00, vu_, vl_, J_, 1 if br_ == "R" else -1))
            except Exception as err:          # a level beyond every validated grid
                print(f"    {iso} {k}: no model position ({err})")
        if rows:
            for name, col in zip(("nu", "S", "el", "vu", "vl", "jl", "br"), zip(*rows)):
                cols[name] = np.concatenate([cols[name], np.array(col, dtype=cols[name].dtype)])
            print(f"    {iso}: {len(rows)} measured lines outside the absorption list, added with zero strength")
    # Stored in the order the explorer lists them, increasing wavelength (decreasing wavenumber), within
    # each shard and from one shard to the next, so the browser never has to sort them.
    order = np.argsort(-cols["nu"], kind="stable")
    nu, S, el = cols["nu"][order], cols["S"][order], cols["el"][order]
    vu, vl, jl, br = cols["vu"][order], cols["vl"][order], cols["jl"][order], cols["br"][order]
    keys = np.array([line_key("R" if b > 0 else "P", int(j), int(u), int(v)) for b, j, u, v in zip(br, jl, vu, vl)])
    is_measured = np.array([k in measured for k in keys], dtype=bool)
    # S(T) = s0 exp(-c2 E''/T) / Q(T), so s0 follows exactly from the strengths at 300 K
    s0 = S * master.partition_function(300.0) * np.exp(C2 * el / 300.0)
    # the position uncertainty and its reason, from the same rules the CLI and the TUI use (lookup.uncertainty)
    u_line, w_line = np.zeros(len(nu)), np.zeros(len(nu), dtype=int)
    for k in range(len(nu)):
        line = Line(iso, "R" if br[k] > 0 else "P", int(jl[k]), int(vu[k]), int(vl[k]), float(nu[k]), 0.0, 0.0, 300.0)
        u_line[k], why = uncertainty(line)
        w_line[k] = reasons.setdefault((why, line.flags), len(reasons))

    # ΔJ = 0 hyperfine patterns for every line, shown by the explorer once the view is narrow enough
    t0 = time.time()
    patterns = hyperfine_patterns(iso, [("R" if b > 0 else "P", int(j), int(u), int(v)) for b, j, u, v in zip(br, jl, vu, vl)])
    print(f"    {iso}: hyperfine patterns of all {len(nu):,} lines in {time.time() - t0:.0f} s", flush=True)

    # even slices of the absorption range; zero-strength rows below it (emission lines) join the first
    edges = np.linspace(nu[S > 0].min(), nu.max() + 1e-6, n_shards + 1)
    edges[0] = min(edges[0], nu.min())
    shards = []
    folder = out / "lines" / iso
    hfs_folder = out / "hfs" / iso
    for f in (folder, hfs_folder):
        if f.exists():
            shutil.rmtree(f)
        f.mkdir(parents=True)
    for k in reversed(range(n_shards)):             # shortest wavelength first, like the rows inside
        sel = np.flatnonzero((nu >= edges[k]) & (nu < edges[k + 1]))
        if not sel.size:
            continue
        rows = [i for i, g in enumerate(sel) if is_measured[g]]
        payload = {
            # 1e-6 cm-1 is 0.03 MHz, well under the 3 MHz the best-known lines are good to
            "nu": [round(float(v), 6) for v in nu[sel]],
            "s0": [float(f"{v:.4g}") for v in s0[sel]],
            "el": [round(float(v), 2) for v in el[sel]],
            "vu": [int(v) for v in vu[sel]],
            "vl": [int(v) for v in vl[sel]],
            "j": [int(j) * (1 if b > 0 else -1) for j, b in zip(jl[sel], br[sel])],  # sign carries the branch
            "u": [float(f"{v:.3g}") for v in u_line[sel]],                         # 1σ position, MHz
            "w": [int(v) for v in w_line[sel]],                                  # ... and why: manifest.reasons
            "m": rows,                                                            # rows with measurements
            "ms": [measured[keys[sel[i]]] for i in rows],                        # ... and what they are
        }
        name = f"{len(shards):03d}.json"
        (folder / name).write_text(json.dumps(payload, separators=(",", ":"), ensure_ascii=False))
        # the same rows' hyperfine patterns, in a file of their own: only a narrow view loads it
        hfs_rows = [patterns[g] or [[], [], None] for g in sel]
        (hfs_folder / name).write_text(json.dumps({"o": [p[0] for p in hfs_rows], "s": [p[1] for p in hfs_rows],
                                                   "x": [p[2] for p in hfs_rows]},
                                                  separators=(",", ":")))
        shards.append({"file": f"lines/{iso}/{name}", "hfs": f"hfs/{iso}/{name}", "nu0": round(float(edges[k]), 4),
                       "nu1": round(float(edges[k + 1]), 4), "n": int(sel.size), "n_measured": len(rows)})

    # the partition function needs the X term values and their nuclear-spin weights
    J = np.arange(master.x_levels.shape[0])
    partition = {"g": [int(nuclear_spin_weight(int(j), iso)) for j in J],
                 "levels": [[round(float(e), 4) for e in row] for row in master.x_levels]}
    (out / f"partition_{iso}.json").write_text(json.dumps(partition, separators=(",", ":")))

    # Hyperfine components for a selection of lines. Strength alone would put every pattern in the
    # green, where the band heads are, and leave 633 nm bare, so take the strongest lines overall
    # *and* the strongest few of each shard. Lines measured with hyperfine resolution are always included.
    strongest = np.argsort(S)[::-1][:n_hfs]
    per_shard = []
    for k in range(n_shards):
        sel = np.flatnonzero((nu >= edges[k]) & (nu < edges[k + 1]))
        per_shard.extend(sel[np.argsort(S[sel])[::-1][:PER_SHARD_HFS]])
    resolved = [i for i in np.flatnonzero(is_measured)
                if any(src["i"] or any(f[0] for f in src["f"]) for src in measured[keys[i]])]
    wanted = list(strongest) + per_shard + resolved
    hfs, t0 = {}, time.time()
    for n, k in enumerate(dict.fromkeys(int(i) for i in wanted)):
        line = Line(iso, "R" if br[k] > 0 else "P", int(jl[k]), int(vu[k]), int(vl[k]), float(nu[k]),
                    float(S[k]), float(el[k]), 300.0)
        try:
            comps = catalog.components(line, main_only=False)
        except Exception:
            continue
        # Keep every labelled (main, ΔF = ΔJ) component and any other that reaches 2% of the strongest.
        # No ΔF ≠ ΔJ component does in practice, so this drops two thirds of the file for nothing visible.
        peak = max(c.strength for c in comps)
        kept = [c for c in comps if c.label or c.strength >= 0.02 * peak]
        total = sum(c.strength for c in kept) or 1.0
        hfs[keys[k]] = {"o": [round(c.offset_MHz, 3) for c in kept],      # 1 kHz, finer than the model is right to
                        "s": [round(c.strength, 5) for c in kept],
                        "l": [c.label or "" for c in kept]}
        # the weaker ones that share a level with a kept component, for the crossovers (strengths on the
        # same scale as "s")
        weak = weak_links(comps, kept, total)
        if weak:
            hfs[keys[k]]["x"] = [[o, round(f * total, 7), iu, il] for o, f, iu, il in weak]
        if n % 200 == 0:
            print(f"    {iso} hyperfine {n}/{len(wanted)} ({time.time() - t0:.0f} s)", flush=True)
    (out / f"hfs_{iso}.json").write_text(json.dumps(hfs, separators=(",", ":")))

    print(f"  {iso}: {len(nu):,} lines ({int(is_measured.sum()):,} measured) in {len(shards)} shards, "
          f"{len(hfs)} hyperfine patterns")
    return {"label": LABELS[iso], "n_lines": int(len(nu)), "n_measured": int(is_measured.sum()),
            "nu_min": round(float(nu.min()), 4), "nu_max": round(float(nu.max()), 4), "shards": shards,
            "order": "wavelength",           # rows and shards run from short to long wavelength
            "partition": f"partition_{iso}.json", "hfs": f"hfs_{iso}.json"}


# --- site and server ---------------------------------------------------------------------------------

def has_data(data=DATA):
    return (Path(data) / "manifest.json").exists()


def build_docs(dest):
    """The documentation (MkDocs Material, mkdocs.yml) into ``dest``."""
    from mkdocs.commands.build import build as mkdocs_build
    from mkdocs.config import load_config

    config = load_config(str(MKDOCS), site_dir=str(Path(dest).resolve()))
    mkdocs_build(config)


def build(dest, data=DATA, docs=True):
    """The explorer and its data at ``dest``, the documentation at ``dest/docs``: the static site."""
    dest = Path(dest)
    if not has_data(data):
        raise FileNotFoundError(f"no exported data in {data}; run `i2spec web export` first")
    dest.mkdir(parents=True, exist_ok=True)
    for name in SITE_FILES:
        shutil.copy2(WEB / name, dest / name)
    if docs:
        build_docs(dest / "docs")
    if (dest / "data").exists():
        shutil.rmtree(dest / "data")
    shutil.copytree(data, dest / "data")
    (dest / ".nojekyll").write_text("")          # serve the files as they are
    return dest


def _lan_address():
    """This machine's address on the local network, as other devices would reach it."""
    import socket
    try:
        with socket.socket(socket.AF_INET, socket.SOCK_DGRAM) as s:
            s.connect(("192.0.2.1", 9))          # TEST-NET: nothing is sent, it only picks the route
            return s.getsockname()[0]
    except OSError:
        return None


def serve(port=8777, open_browser=True, directory=ROOT / "site", host="127.0.0.1"):
    """Serve the site until interrupted: on this machine only by default, on the network with host="0.0.0.0"."""
    handler = functools.partial(http.server.SimpleHTTPRequestHandler, directory=str(directory))
    handler.log_message = lambda *args: None
    for p in range(port, port + 20):
        try:
            server = http.server.ThreadingHTTPServer((host, p), handler)
            break
        except OSError:
            continue
    else:
        raise OSError(f"no free port in {port}-{port + 19}")
    url = f"http://localhost:{server.server_address[1]}/"
    print(f"Iodine Line Explorer at {url}  (Ctrl-C to stop)", flush=True)
    if host not in ("127.0.0.1", "localhost"):
        lan = _lan_address()
        print(f"on the local network: http://{lan or '<this machine>'}:{server.server_address[1]}/", flush=True)
    if open_browser:
        webbrowser.open(url)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        server.server_close()
