"""The web app's export pieces that need no model: sources, measurements, the site build, the CLI."""
import json

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
    # the thesis and the printed Orsay Partie IV volume are the only sources with no online copy
    unlinked = [s["id"] for s in sources if not s["link"]]
    assert unlinked == ["bodermann1998c", "orsay1983_part4"]
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
        absolute, intervals = by_source[k]
        assert absolute == s["n_absolute"], s["id"]
        assert intervals >= s["n_intervals"], s["id"]     # an interval to another line counts for both
    # the residual is the published value minus the model's, here minus zero
    r56 = measured["127I2"]["32-0R56"]
    # bipm2012a's a10 defers to its sources (observations.load_all); Jones 2002's 532 nm row carries no correction
    jones = next(e for e in r56 if sources[e["s"]]["id"] == "jones2002a")
    a10 = next(f for f in jones["f"] if f[0] == "a10")
    assert a10[1] == a10[3] == 563260223.5144
    assert any(r is not None for r in calls)             # components are asked for by rank


def test_every_dataset_has_a_doi_or_a_document():
    for ds in load_all():
        assert ds.meta.get("doi") or "http" in ds.meta["source"] or ds.id == "bodermann1998c", ds.id


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
