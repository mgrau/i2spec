"""Full-resolution reading crops of the Orsay I2 atlas, part IV (19700-20035 cm-1, Gerstenkorn & Luc 1983).

The table region of a photo (a box in fractions of the rotated photo) is cut into column groups, found
at vertical white gaps wider than --gap px, and each group into --bands horizontal bands that overlap
by --overlap px, then flattened and enlarged by --scale. Each tile is kept below ~1500 px on its long
side, so that a viewer shows it without downsampling; aim for glyphs of 25 px or more. Digits are read
by eye from the tiles (orsay_ocr.py: OCR misreads about a third of them).

usage: orsay4_crops.py PHOTO OUTDIR [--rot=R] [--box=x0,y0,x1,y1] [--gap=35] [--bands=2] [--overlap=80] [--scale=2]
  plate page:           --box=0.05,0.58,0.98,0.95 --bands=2 --scale=2   (4 column groups of N sigma eps I)
  classification page:  --box=0.05,0.08,0.98,0.95 --gap=150 --bands=3 --scale=1
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from PIL import Image

from orsay_ocr import binarize, flatten, glyphs_only, load, runs


def tiles(path, outdir, rot=0, box=(0, 0, 1, 1), gap=35, bands=2, overlap=80, scale=2.0, pad=20):
    im = flatten(load(path, rot))
    W, H = im.size
    x0, y0, x1, y1 = int(box[0] * W), int(box[1] * H), int(box[2] * W), int(box[3] * H)
    region = im.crop((x0, y0, x1, y1))
    ink = glyphs_only(binarize(region), min_size=10)
    rows = np.flatnonzero(ink.sum(1) > 3)
    ya, yb = (rows.min(), rows.max()) if len(rows) else (0, region.height)
    cols = runs(ink[ya:yb].sum(0) > 2, gap, min_len=60)
    stem = Path(path).stem.replace("PXL_20260929_", "")
    out = []
    for ci, (a, b) in enumerate(cols):
        step = (yb - ya) / bands
        for bi in range(bands):
            t = max(int(ya + bi * step) - overlap - pad, 0)
            u = min(int(ya + (bi + 1) * step) + overlap + pad, region.height)
            tile = region.crop((max(a - pad, 0), t, min(b + pad, region.width), u))
            if scale != 1:
                tile = tile.resize((int(tile.width * scale), int(tile.height * scale)), Image.LANCZOS)
            f = Path(outdir) / f"{stem}_c{ci}_b{bi}.png"
            tile.save(f)
            out.append((f, tile.size, (x0 + a, y0 + t, x0 + b, y0 + u)))
    return out


def main(argv):
    opts = dict(a[2:].split("=", 1) for a in argv if a.startswith("--"))
    path, outdir = [a for a in argv if not a.startswith("--")]
    Path(outdir).mkdir(parents=True, exist_ok=True)
    for f, size, bbox in tiles(path, outdir, rot=int(opts.get("rot", 0)),
                               box=tuple(float(v) for v in opts.get("box", "0,0,1,1").split(",")),
                               gap=int(opts.get("gap", 35)), bands=int(opts.get("bands", 2)),
                               overlap=int(opts.get("overlap", 80)), scale=float(opts.get("scale", 2))):
        warn = "\tTOO WIDE: columns merged, lower --gap" if size[0] > 1600 else ""
        print(f"{f}\t{size[0]}x{size[1]}\tphoto px {tuple(int(v) for v in bbox)}{warn}")


if __name__ == "__main__":
    main(sys.argv[1:])
