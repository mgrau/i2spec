"""MkDocs hooks: the references page is written from the data files at build time, so it always lists
exactly the sources the model and the explorer use (i2spec.webapp.sources); and every number of the model a page
shows, written {{ key }} in its source, is computed from the released model at build time (model_numbers.py)."""
import html
import importlib.util
import re
import sys
from pathlib import Path

MARKER = "<!-- sources -->"


def _nm(r):
    if not r:
        return "–"
    return f"{1e7 / r[0]:.1f}" if r[0] == r[1] else f"{1e7 / r[1]:.1f}–{1e7 / r[0]:.1f}"


def _u(v):
    return f"{v:.1f} MHz" if v >= 1 else (f"{v * 1e3:.0f} kHz" if v >= 1e-3 else f"{v * 1e3:.2f} kHz")


def _link(s):
    if s["doi"]:
        return f'[{s["doi"]}](https://doi.org/{s["doi"]})'
    if s["link"]:
        host = s["link"].split("/")[2].removeprefix("www.")
        return f"[{host}]({s['link']})"
    return "no online copy"


def _source(s, extra=""):
    return (f'**{s["short"]}** · {_link(s)}<br><span class="cite">{html.escape(s["citation"])}</span>'
            + (f'<br><span class="cite">{extra}</span>' if extra else ""))


def _rows(s):
    iso = model_numbers().ISO
    n_int = s["n_intervals"]
    parts = [f'{s["n_absolute"]} absolute' if s["n_absolute"] else "",
             f'{n_int} interval{"" if n_int == 1 else "s"}' if n_int else ""]
    also = s.get("n_also", 0)
    if not s["n_rows"]:              # every row printed first elsewhere: cited under the paper that reported it
        return (" ".join(iso.get(i, i) for i in s["isotopologues"])
                + f' · prints again {also} measurement{"" if also == 1 else "s"} listed under the paper'
                  f'{"" if also == 1 else "s"} that first reported {"it" if also == 1 else "them"}')
    return (f'{s["n_rows"]} row{"" if s["n_rows"] == 1 else "s"} ({", ".join(p for p in parts if p)}) · '
            + " ".join(iso.get(i, i) for i in s["isotopologues"])
            + (f' · also prints {also} measurement{"" if also == 1 else "s"} listed under the papers that '
               f'report{"s" if also == 1 else ""} {"it" if also == 1 else "them"}' if also else ""))


def _tables():
    from i2spec.webapp import sources

    src = sources()
    year = lambda s: int("".join(c for c in s["id"] if c.isdigit())[:4] or 0)  # noqa: E731
    precision = sorted((s for s in src if s["kind"] == "precision" and s.get("fit", True)),
                       key=lambda s: (year(s), s["short"]))
    blind = sorted((s for s in src if s["kind"] == "precision" and not s.get("fit", True)),
                   key=lambda s: (year(s), s["short"]))
    atlas = [s for s in src if s["kind"] == "atlas"]
    n = lambda kind: sum(s["n_lines"] for s in src if s["kind"] == kind)  # noqa: E731
    out = [f'{len(precision)} precision data sets in the fit, {len(blind)} more in no fit (held out for the '
           f'blind tests, or set aside), covering {n("precision"):,} lines together, and {len(atlas)} atlases with '
           f'{n("atlas"):,} assigned lines.', "",
           '## <span class="kind precision"></span>Sub-Doppler and frequency-comb measurements', "",
           "Absolute frequencies of hyperfine components or line centers, and hyperfine intervals, with "
           "uncertainties between 0.25 kHz and a few MHz. They determine the potentials and the hyperfine "
           "functions of the model. Under each citation: the "
           "number of absolute frequencies and hyperfine intervals, and the isotopologues.", "",
           "| source | lines | range (nm) | median u |",
           "|---|--:|--:|--:|"]
    for s in precision:
        out.append(f'| {_source(s, _rows(s))} | {s["n_lines"]:,} | {_nm(s["nu"])} | {_u(s["u_MHz"])} |')
    if blind:
        out += ["", '## <span class="kind precision"></span>In no fit', "",
                "Most were transcribed after the fit and held out of it for the blind tests: the model's predictions "
                "of them, and their uncertainties, were frozen before their values were read, then scored once "
                "(results on the [uncertainty page](uncertainty.md)). The others were set aside, each for the "
                "reason given under it. The explorer shows them all, marked “not fitted”.", "",
                "| source | lines | range (nm) | median u |", "|---|--:|--:|--:|"]
        for s in blind:
            why = "" if "blind tests" in s.get("why", "blind tests") else \
                f'<br><span class="cite">Set aside: {html.escape(s["why"])}.</span>'
            out.append(f'| {_source(s, _rows(s))}{why} | {s["n_lines"]:,} | {_nm(s["nu"])} | {_u(s["u_MHz"])} |')
    out += ["", '## <span class="kind atlas"></span>Fourier-transform atlases', "",
            "Doppler-limited line centers, assigned in this work by comparison with the model. The Salami & Ross "
            "and APO atlases are in the core fit; Orsay part I and Partie IV determine the levels beyond it "
            "(X v″ = 18–25, B v′ = 51–80). The others are shown for comparison and used in no fit (the 1977 and "
            "1978 atlases calibrate the uncertainties of atlas-class lines). The "
            "calibration of each atlas is given with it.", "",
            "| source | lines | range (nm) | median u |", "|---|--:|--:|--:|"]
    for s in atlas:
        fitted = "" if s.get("fit", True) else " Not used in the fit."
        out.append(f'| {_source(s)}<br><span class="cite">Calibration: {html.escape(s["uncertainty"])}{fitted}</span> | '
                   f'{s["n_lines"]:,} | {_nm(s["nu"])} | {_u(s["u_MHz"])} |')
    return "\n".join(out)


PLACEHOLDER = re.compile(r"\{\{\s*([A-Za-z0-9_.]+)\s*\}\}")
_NUMBERS = "i2spec_docs_model_numbers"


def model_numbers():
    """documentation/model_numbers.py, loaded by path (the documentation folder is not a package) and kept in
    sys.modules until the file changes. MkDocs loads this file afresh for every build (each rebuild of
    ``mkdocs serve``; the site build of tests/test_webapp.py after tests/test_docs.py), and the numbers take tens of
    seconds: this way they are computed once per process, not once per build."""
    path = Path(__file__).with_name("model_numbers.py")
    mtime = path.stat().st_mtime_ns
    module = sys.modules.get(_NUMBERS)
    if module is None or getattr(module, "_mtime", None) != mtime:
        spec = importlib.util.spec_from_file_location(_NUMBERS, path)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        module._mtime = mtime
        sys.modules[_NUMBERS] = module
    return module


def fill(markdown, where="page"):
    """Every {{ key }} of a page replaced by its number; an unknown key stops the build. A value of several lines
    (a table's rows) that stands alone on its line keeps that line's indentation on each of its lines."""
    if not PLACEHOLDER.search(markdown):
        return markdown
    values = model_numbers().numbers()

    def one(match):
        if match[1] not in values:
            raise KeyError(f"{where}: no generated number {match[1]!r} (documentation/model_numbers.py)")
        prefix = markdown[markdown.rfind("\n", 0, match.start()) + 1:match.start()]
        value = values[match[1]]
        return value.replace("\n", "\n" + prefix) if not prefix.strip() else value
    return PLACEHOLDER.sub(one, markdown)


def on_page_markdown(markdown, page, config, files):
    if MARKER in markdown:
        markdown = markdown.replace(MARKER, f'<div class="refs" markdown>\n\n{_tables()}\n\n</div>')
    return fill(markdown, page.file.src_path)
