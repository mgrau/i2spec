"""What make_figures.py and pipeline.py share: the placeholder colors the figures are drawn with, and save(), which
writes a figure's SVG to this folder with each placeholder replaced by a CSS variable of the site (web/site.css), so
the pages that inline it follow the light and dark themes."""
import re
from pathlib import Path

import matplotlib.pyplot as plt

OUT = Path(__file__).resolve().parent
# placeholder colors -> the site's tokens
C = dict(ink="#010101", ink2="#020202", ink3="#030303", rule="#040404", trace="#050505", pick="#060606",
         mark="#070707", surface="#080808", u1="#0a0a01", u2="#0a0a02", u3="#0a0a03", u4="#0a0a04", u5="#0a0a05")
TOKENS = {**{C[f"u{k}"]: f"var(--unc-{k})" for k in range(1, 6)}, C["ink"]: "var(--ink)", C["ink2"]: "var(--ink-2)",
          C["ink3"]: "var(--ink-3)", C["rule"]: "var(--rule)", C["trace"]: "var(--trace)", C["pick"]: "var(--pick)",
          C["mark"]: "var(--mark)", C["surface"]: "var(--surface)"}


def save(fig, name, **savefig):
    """fig as OUT/<name>.svg (savefig's keyword arguments as given), its colors and fonts the site's."""
    path = OUT / f"{name}.svg"
    fig.savefig(path, **savefig)
    plt.close(fig)
    svg = path.read_text()
    for hexcode, var in TOKENS.items():
        svg = svg.replace(hexcode, var)
    svg = re.sub(r"font-family:[^;\"]*", "font-family:var(--body)", svg)
    svg = re.sub(r"font: ([\d.]+)px '[^']*'", r"font-size:\1px;font-family:var(--body)", svg)
    svg = re.sub(r'<svg([^>]*?) width="[^"]*" height="[^"]*"', r'<svg\1', svg, count=1)
    svg = re.sub(r"<metadata>.*?</metadata>", "", svg, flags=re.S)
    path.write_text(svg)
    print(f"{path.name}: {len(svg) / 1e3:.0f} kB")
