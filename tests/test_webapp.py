"""The web app's export pieces that need no model: sources, measurements, the site build, the CLI."""
import json

import pytest

from i2spec import webapp
from i2spec.cli import build_parser
from i2spec.observations import load_all


def test_every_source_is_cited_and_linked():
    sources = webapp.sources()
    ids = [s["id"] for s in sources]
    assert len(ids) == len(set(ids))
    assert len({s["short"] for s in sources}) == len(sources), "short names must tell sources apart"
    assert {s["kind"] for s in sources} == {"precision", "atlas"}
    for s in sources:
        assert len(s["citation"]) > 30 and s["n_lines"] > 0
    # the theses and two printed atlas volumes (the 1978 one photographed, Orsay Partie IV) are the only sources
    # with no online copy
    unlinked = [s["id"] for s in sources if not s["link"]]
    assert unlinked == ["bodermann1998c", "chen2005thesis_a", "chen2005thesis_b", "eickhoff1994thesis",
                        "gerstenkorn1978b", "orsay1983_part4"]
    assert all(s["link"].startswith("https://doi.org/") for s in sources if s["doi"])


def test_measurements_carry_every_row():
    calls = []

    def model(iso, branch, J, vu, vl, rank):
        calls.append(rank)
        return 0.0

    measured = webapp.measurements(model)
    sources = webapp.sources()
    by_source = {}
    for lines in measured.values():
        for entries in lines.values():
            for e in entries:
                n = by_source.setdefault(e["s"], [0, 0])
                n[0] += len(e["f"])
                n[1] += e["i"]
    for k, s in enumerate(sources):
        if not s["n_rows"]:                 # prints only what others reported first: cited on their rows
            assert k not in by_source and s["n_also"] > 0, s["id"]
            continue
        absolute, intervals = by_source[k]
        assert absolute == s["n_absolute"], s["id"]
        assert intervals >= s["n_intervals"], s["id"]     # an interval to another line counts for both
    # the residual is the published value minus the model's, here minus zero
    r56 = measured["127I2"]["32-0R56"]
    # one record per measurement: Jones 2002's a10 is cited to Jones first, then to the CIPM 532 nm value that
    # averages it, which is a record of its own; Jones 2002's 532 nm row carries no correction
    jones = next(e for e in r56 if sources[e["s"]]["id"] == "jones2002a")
    assert [sources[c]["id"] for c in jones["c"]] == ["bipm2012a"]
    cipm = next(e for e in r56 if sources[e["s"]]["id"] == "bipm2012a")
    # the 2005 CIPM list (shown since the third blind test) prints the same value: cited on the CIPM record
    assert [f[1] for f in cipm["f"]] == [563260223.513] and [sources[c]["id"] for c in cipm["c"]] == ["felder2005a"]
    a10 = next(f for f in jones["f"] if f[0] == "a10")
    assert a10[1] == a10[3] == 563260223.5144 and len(a10) == 4
    assert any(r is not None for r in calls)             # components are asked for by rank


def test_every_dataset_has_a_doi_or_a_document():
    from i2spec.catalog import entry_of
    for ds in load_all():
        e = entry_of(ds.meta)
        assert e.get("doi") or "http" in str(e.get("url") or "") + ds.meta["source"] or ds.id == "bodermann1998c", ds.id


def test_build_copies_the_site(tmp_path):
    data = tmp_path / "data"
    (data / "lines").mkdir(parents=True)
    (data / "manifest.json").write_text(json.dumps({"isotopologues": {}}))
    site = webapp.build(tmp_path / "site", data=data)
    for name in webapp.SITE_FILES:
        assert (site / name).exists(), name
    assert (site / "data" / "manifest.json").exists() and (site / ".nojekyll").exists()
    # the documentation sits under docs/, flat, with its figures and without its build scripts
    for page in ("index", "calculation", "levels", "hyperfine", "references"):
        assert (site / "docs" / f"{page}.html").exists(), page
    assert (site / "docs" / "figures" / "potentials.svg").exists()
    assert not (site / "docs" / "hooks.py").exists() and not (site / "docs" / "figures" / "make_figures.py").exists()
    # the references page is written from the data at build time
    refs = (site / "docs" / "references.html").read_text(encoding="utf-8")
    assert "10.1088/0026-1394/44/5/003" in refs and "Salami &amp; Ross 2005" in refs


def test_explorer_links_to_the_docs():
    html = (webapp.WEB / "index.html").read_text(encoding="utf-8")
    assert 'href="docs/index.html"' in html and 'href="docs/references.html"' in html
    # and the documentation leads back: its logo and first navigation entry
    config = (webapp.ROOT / "mkdocs.yml").read_text(encoding="utf-8")
    assert "homepage: ../index.html" in config and ": ../index.html" in config.split("\nnav:")[1].splitlines()[1]


def test_docs_figures_exist_and_follow_the_theme():
    import re
    docs = webapp.ROOT / "documentation"
    names = set()
    for page in docs.glob("*.md"):
        names |= set(re.findall(r'data-svg="figures/([a-z_]+)\.svg"', page.read_text(encoding="utf-8")))
    assert len(names) >= 6
    for name in names:
        svg = (docs / "figures" / f"{name}.svg").read_text()
        assert "var(--" in svg and not re.search(r"#0[1-8]0[1-8]0[1-8]", svg), name


def test_web_command_defaults():
    args = build_parser().parse_args(["web"])
    assert (args.action, args.port, args.no_browser) == ("serve", 8777, False)
    args = build_parser().parse_args(["web", "build", "out"])
    assert (args.action, args.dest) == ("build", "out")


def test_reference_generator_is_safe_to_import_in_spawned_workers(monkeypatch, capsys):
    import runpy
    from i2spec import intensity, lookup

    def unexpected_calculation(*args, **kwargs):
        raise AssertionError("Importing the reference generator must not calculate a spectrum")

    monkeypatch.setattr(intensity, "intensity_model", unexpected_calculation)
    monkeypatch.setattr(lookup, "Catalog", unexpected_calculation)
    module = runpy.run_path(str(webapp.WEB / "make_reference.py"), run_name="__mp_main__")
    assert callable(module["make_reference"]) and callable(module["main"])
    assert capsys.readouterr().out == ""


def test_web_commands_need_a_checkout_and_say_so_at_once(monkeypatch, tmp_path, capsys):
    import pytest

    from i2spec import cli

    webapp.check_checkout()                          # this is a checkout
    monkeypatch.setattr(webapp, "OBSERVATIONS_DIR", tmp_path / "data" / "observations")
    monkeypatch.setattr(webapp, "export", lambda **kw: (_ for _ in ()).throw(AssertionError("exported")))
    for action in ("build", "export", "serve"):      # as from an installed wheel: no data/ beside the package
        assert cli.main(["web", action, "--no-browser"]) == 2, action
        assert "needs a source checkout" in capsys.readouterr().err
    monkeypatch.undo()
    monkeypatch.setattr(webapp, "MKDOCS", tmp_path / "mkdocs.yml")
    with pytest.raises(FileNotFoundError, match="mkdocs.yml"):
        webapp.build(tmp_path / "site", data=tmp_path)


def test_web_build_warns_when_s_min_differs_from_the_export(monkeypatch, capsys):
    from i2spec import cli

    monkeypatch.setattr(webapp, "has_data", lambda data=None: True)
    monkeypatch.setattr(webapp, "data_cutoff", lambda data=None: 1e-24)
    monkeypatch.setattr(webapp, "data_inputs", lambda data=None: "same")
    monkeypatch.setattr(webapp, "export_inputs", lambda: "same")
    monkeypatch.setattr(webapp, "export", lambda **kw: pytest.fail("exported"))
    monkeypatch.setattr(webapp, "build", lambda dest, **kw: dest)
    assert cli.main(["web", "build", "out"]) == 0
    assert "--s-min" not in capsys.readouterr().err
    assert cli.main(["web", "build", "out", "--s-min", "1e-26"]) == 0
    assert "--s-min 1e-24, not 1e-26" in capsys.readouterr().err


def test_web_build_exports_again_only_when_its_inputs_changed(monkeypatch, capsys):
    from i2spec import cli

    exported = []
    monkeypatch.setattr(webapp, "has_data", lambda data=None: True)
    monkeypatch.setattr(webapp, "data_cutoff", lambda data=None: 1e-24)
    monkeypatch.setattr(webapp, "export_inputs", lambda: "now")
    monkeypatch.setattr(webapp, "export", lambda **kw: exported.append(kw))
    monkeypatch.setattr(webapp, "build", lambda dest, **kw: dest)
    monkeypatch.setattr(webapp, "data_inputs", lambda data=None: "now")
    assert cli.main(["web", "build", "out"]) == 0 and not exported
    monkeypatch.setattr(webapp, "data_inputs", lambda data=None: "before")
    assert cli.main(["web", "build", "out"]) == 0 and len(exported) == 1
    assert "changed since the last export" in capsys.readouterr().out
    # serving locally keeps the older data, and says so
    monkeypatch.setattr(webapp, "serve", lambda **kw: None)
    assert cli.main(["web", "--no-browser"]) == 0 and len(exported) == 1
    assert "changed since the last export" in capsys.readouterr().err
    assert cli.main(["web", "inputs"]) == 0 and capsys.readouterr().out.strip() == "now"


def test_manifest_and_shard_bounds(tmp_path):
    import http.server

    assert not webapp.has_data(tmp_path) and webapp.data_cutoff(tmp_path) is None
    webapp._write_atomic(tmp_path / "manifest.json", json.dumps({"cutoff_cm_at_300K": 1e-25}))
    assert webapp.has_data(tmp_path) and webapp.data_cutoff(tmp_path) == 1e-25
    assert [p.name for p in tmp_path.iterdir()] == ["manifest.json"]   # no temporary file left behind
    # shard bounds are rounded outward, so a line stored to 1e-6 cm-1 near an edge stays inside its shard
    assert webapp._floor4(18788.123456) == 18788.1234 <= round(18788.123456, 6)
    assert webapp._ceil4(18788.123449) == 18788.1235 >= round(18788.123449, 6)
    assert webapp._floor4(18788.1234) <= 18788.1234 <= webapp._ceil4(18788.1234)
    # the server logs nothing per request
    assert webapp._QuietHandler.log_message is not http.server.SimpleHTTPRequestHandler.log_message


def test_overview_sums_each_kind_per_bin_and_temperature():
    import math
    import numpy as np
    from i2spec.intensity import C2
    nu = np.array([100.5, 100.6, 101.5, 109.9])
    s0, el, kind = np.array([1e-18, 2e-18, 3e-18, 0.0]), np.array([0.0, 100.0, 200.0, 0.0]), np.array([0, 0, 2, 1])
    q = [2.0] * len(webapp.OVERVIEW_T)
    ov = webapp._overview(nu, s0, el, kind, 100.0, 110.0, C2, q)
    bins = webapp.OVERVIEW_BINS
    assert len(ov["p"]) == 3 and all(len(rows) == len(webapp.OVERVIEW_T) and len(rows[0]) == bins for rows in ov["p"])
    T, k0 = webapp.OVERVIEW_T[2], int(0.5 / 10 * bins)              # both model lines fall in one bin
    S = (1e-18 + 2e-18 * math.exp(-C2 * 100.0 / T)) / 2.0
    assert ov["p"][0][2][k0] == round(1000 * math.log10(S))
    assert sum(1 for v in ov["p"][0][2] if v) == 1 and sum(1 for v in ov["p"][2][2] if v) == 1
    assert not any(ov["p"][1][2])                       # the atlas line has no strength: an empty bin
    # the atlas line still marks its bin, as the precision line marks its own
    assert ov["m"][int(9.9 / 10 * bins)] == "1" and ov["m"][int(1.5 / 10 * bins)] == "2" and len(ov["m"]) == bins


def test_export_inputs_follow_what_the_export_reads(tmp_path, monkeypatch):
    # two data sets: one loaded, one pending; only the loaded one is part of the digest
    for name, meta in (("loaded1999a", 'source = "x"\n'), ("pending1999a", 'source = "x"\nexclude = "pending"\n')):
        (tmp_path / name).mkdir()
        (tmp_path / name / "meta.toml").write_text(meta)
        (tmp_path / name / "data.csv").write_text("a\n1\n")
    monkeypatch.setattr(webapp, "OBSERVATIONS_DIR", tmp_path)
    monkeypatch.setattr(webapp, "blind_scored", lambda: frozenset())
    first = webapp.export_inputs()
    assert webapp.export_inputs() == first
    (tmp_path / "pending1999a" / "data.csv").write_text("a\n2\n")
    assert webapp.export_inputs() == first
    (tmp_path / "loaded1999a" / "data.csv").write_text("a\n2\n")
    assert webapp.export_inputs() != first
    # a scored blind set is loaded, so it counts
    second = webapp.export_inputs()
    monkeypatch.setattr(webapp, "blind_scored", lambda: frozenset({"pending1999a"}))
    assert webapp.export_inputs() != second
    # the exported data carry the digest they were made from
    webapp._write_atomic(tmp_path / "manifest.json", json.dumps({"inputs": "abc"}))
    assert webapp.data_inputs(tmp_path) == "abc" and webapp.data_inputs(tmp_path / "none") is None


def test_hyperfine_patterns_are_cached(tmp_path, monkeypatch):
    monkeypatch.setenv("I2SPEC_CACHE", str(tmp_path))
    calls = []

    def compute(iso, lines, workers=None):
        calls.append(list(lines))
        return [("[1.0]", "[1000]", "null") for _ in lines]

    monkeypatch.setattr(webapp, "_compute_hyperfine_patterns", compute)
    lines = [("R", 56, 32, 0), ("P", 53, 32, 0)]
    first = webapp.hyperfine_patterns("127I2", lines)
    again = webapp.hyperfine_patterns("127I2", lines)
    assert len(calls) == 1 and [list(p) for p in first] == again        # read back from the cache
    webapp.hyperfine_patterns("127I2", lines[:1])                        # other lines: another key
    webapp.hyperfine_patterns("127I2", lines, cache=False)
    assert len(calls) == 3 and len(list(tmp_path.glob("hfs_patterns_127I2_*.json.gz"))) == 2
