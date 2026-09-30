"""MkDocs hooks: the references page is written from the data files at build time, so it always lists
exactly the sources the model and the explorer use (i2spec.webapp.sources)."""
import html

MARKER = "<!-- sources -->"
ISO = {"127I2": "¹²⁷I₂", "129I2": "¹²⁹I₂", "127I129I": "¹²⁷I¹²⁹I"}


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
    n_int = s["n_intervals"]
    parts = [f'{s["n_absolute"]} absolute' if s["n_absolute"] else "",
             f'{n_int} interval{"" if n_int == 1 else "s"}' if n_int else ""]
    return (f'{s["n_rows"]} row{"" if s["n_rows"] == 1 else "s"} ({", ".join(p for p in parts if p)}) · '
            + " ".join(ISO.get(i, i) for i in s["isotopologues"]))


def _tables():
    from i2spec.webapp import sources

    src = sources()
    year = lambda s: int("".join(c for c in s["id"] if c.isdigit())[:4] or 0)  # noqa: E731
    precision = sorted((s for s in src if s["kind"] == "precision"), key=lambda s: (year(s), s["short"]))
    atlas = [s for s in src if s["kind"] == "atlas"]
    n = lambda kind: sum(s["n_lines"] for s in src if s["kind"] == kind)  # noqa: E731
    out = [f'{len(precision)} precision data sets covering {n("precision"):,} lines, and {len(atlas)} atlases '
           f'with {n("atlas"):,} assigned lines.', "",
           '## <span class="kind precision"></span>Sub-Doppler and frequency-comb measurements', "",
           "Absolute frequencies of hyperfine components or line centres, and hyperfine intervals, with "
           "uncertainties between 0.25 kHz and a few MHz. They determine the level energies, the level corrections "
           "and the hyperfine parameters. Under each citation: the number of absolute frequencies and hyperfine intervals, and the isotopologues.", "",
           "| source | lines | range (nm) | median u |",
           "|---|--:|--:|--:|"]
    for s in precision:
        out.append(f'| {_source(s, _rows(s))} | {s["n_lines"]:,} | {_nm(s["nu"])} | {_u(s["u_MHz"])} |')
    out += ["", '## <span class="kind atlas"></span>Fourier-transform atlases', "",
            "Doppler-limited line centres, assigned in this work by comparison with the model. They determine "
            "levels outside the range of the precision data (X v″ = 18–25, from the Orsay atlas) and are used to "
            "test the model elsewhere. The calibration of each atlas is given with it.", "",
            "| source | lines | range (nm) | median u |", "|---|--:|--:|--:|"]
    for s in atlas:
        out.append(f'| {_source(s)}<br><span class="cite">Calibration: {html.escape(s["uncertainty"])}</span> | '
                   f'{s["n_lines"]:,} | {_nm(s["nu"])} | {_u(s["u_MHz"])} |')
    return "\n".join(out)


def on_page_markdown(markdown, page, config, files):
    if MARKER in markdown:
        return markdown.replace(MARKER, f'<div class="refs" markdown>\n\n{_tables()}\n\n</div>')
    return markdown
