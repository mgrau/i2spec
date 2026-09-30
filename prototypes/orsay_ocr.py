"""Transcribe the line tables of the Orsay I2 atlas (Gerstenkorn, Vergès & Chevillard 1982) from photographs.

Each table row is  N  sigma  eps  I  : tick number, vacuum wavenumber (cm-1, 4 decimals), pointing
uncertainty (mk, one decimal), absorption depth (integer). Nothing is trusted from OCR alone. A row is
accepted only if two OCR passes with different preprocessing read it identically and it fits the
table's own structure: N consecutive, sigma increasing, eps in 0.1-9.9, I a small integer. Everything
else is listed for reading by eye. docs/research/orsay-atlas-11000-14000.md.

usage: orsay_ocr.py orient PHOTO_DIR           -> PHOTO_DIR/orient.tsv (rotation, sigma range per photo)
       orsay_ocr.py page PHOTO [--rot=R]       -> rows of one photo, with the verdict for each
"""
from __future__ import annotations

import io
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
from PIL import Image, ImageFilter, ImageOps

SIGMA = re.compile(r"1[1-4]\d{3}\.\d{4}")


def tesseract(img: Image.Image, psm: int = 6) -> str:
    buf = io.BytesIO(); img.save(buf, "PNG")
    p = subprocess.run(["tesseract", "stdin", "stdout", "--psm", str(psm), "-c",
                        "tessedit_char_whitelist=0123456789. "], input=buf.getvalue(), capture_output=True)
    return p.stdout.decode("utf8", "replace")


def load(path, rot=0, max_side=None):
    im = Image.open(path)
    if max_side:
        im.draft("RGB", (max_side, max_side))
    im = ImageOps.exif_transpose(im).convert("L")
    return im.rotate(rot, expand=True) if rot else im


def flatten(im: Image.Image, radius=40) -> Image.Image:
    """Divide out the illumination gradient (a heavy blur is the background), then stretch."""
    a = np.asarray(im, dtype=float)
    bg = np.asarray(im.filter(ImageFilter.GaussianBlur(radius)), dtype=float) + 1
    b = np.clip(a / bg * 235, 0, 255).astype(np.uint8)
    return ImageOps.autocontrast(Image.fromarray(b), cutoff=1)


def orient(path):
    """The rotation (0/90/180/270) under which OCR finds the most wavenumbers, and those wavenumbers."""
    best = (-1, 0, [])
    for rot in (0, 90, 180, 270):
        im = flatten(load(path, rot, max_side=1600))
        s = [float(x) for x in SIGMA.findall(tesseract(im).replace(" ", ""))]
        s = [x for x in s if 11000 <= x <= 14200]
        if len(s) > best[0]:
            best = (len(s), rot, s)
    return best[1], best[2]


ROW = re.compile(r"^\s*(\d{1,4})\s+(1[1-4]\d{3})\s*\.\s*(\d{4})\s+(\d)\s*\.\s*(\d)\s+(\d{1,2})\s*$")


def binarize(im: Image.Image) -> np.ndarray:
    """True where there is ink, by Otsu's threshold on the flattened page."""
    a = np.asarray(im, dtype=np.uint8)
    hist = np.bincount(a.ravel(), minlength=256).astype(float)
    w = np.cumsum(hist); m = np.cumsum(hist * np.arange(256)); mt = m[-1]; wt = w[-1]
    between = (mt * w - m * wt) ** 2 / np.maximum(w * (wt - w), 1)
    return a < int(np.argmax(between))


def runs(mask, min_gap, min_len=1):
    """(start, stop) of runs of True in a 1-D mask, merging runs closer than min_gap."""
    idx = np.flatnonzero(mask)
    if not len(idx):
        return []
    out, lo, prev = [], idx[0], idx[0]
    for i in idx[1:]:
        if i - prev > min_gap:
            out.append((lo, prev + 1)); lo = i
        prev = i
    out.append((lo, prev + 1))
    return [r for r in out if r[1] - r[0] >= min_len]


def glyphs_only(ink: np.ndarray, max_size=110, min_size=4) -> np.ndarray:
    """Keep only connected ink blobs of character size: drops the desk, fingers, page edges and specks."""
    from scipy import ndimage
    lab, n = ndimage.label(ink)
    sl = ndimage.find_objects(lab)
    keep = np.zeros(n + 1, bool)
    for k, s_ in enumerate(sl, start=1):
        h, w = s_[0].stop - s_[0].start, s_[1].stop - s_[1].start
        keep[k] = min_size <= max(h, w) <= max_size
    return keep[lab]


def segment(ink: np.ndarray, col_gap=90, row_gap=1):
    """Column groups (by vertical white gaps wider than col_gap px), then text rows within each."""
    ink = glyphs_only(ink, min_size=18)          # characters only: no speckle, no decimal points (OCR sees the original)
    cols = runs(ink.sum(0) > 2, col_gap, min_len=120)
    groups = []
    for x0, x1 in cols:
        sub = ink[:, x0:x1]
        rows = runs(sub.sum(1) > 0, row_gap, min_len=15)
        rows = [(y0, y1) for y0, y1 in rows if y1 - y0 < 80]       # one text line is ~35 px; drop anything taller
        groups.append(((x0, x1), rows))
    return groups


def parse(text):
    m = ROW.match(text.strip())
    if not m:
        return None
    n, a, b, e1, e2, i = m.groups()
    return int(n), float(f"{a}.{b}"), float(f"{e1}.{e2}"), int(i)


def read_page(path, rot):
    """Every text row of the photo: (group, y, bbox, reading A, reading B, parsed or None)."""
    im = flatten(load(path, rot))
    ink = binarize(im)
    out = []
    for gi, ((x0, x1), rows) in enumerate(segment(ink)):
        for y0, y1 in rows:
            pad = 6
            box = (max(x0 - pad, 0), max(y0 - pad, 0), min(x1 + pad, im.width), min(y1 + pad, im.height))
            crop = im.crop(box)
            a = tesseract(crop.resize((crop.width * 2, crop.height * 2), Image.LANCZOS), psm=7)
            bw = Image.fromarray(np.where(binarize(crop), 0, 255).astype(np.uint8))
            b = tesseract(bw.resize((bw.width * 3, bw.height * 3), Image.NEAREST).filter(ImageFilter.MedianFilter(3)), psm=7)
            pa, pb = parse(a), parse(b)
            out.append(dict(group=gi, box=box, a=a.strip(), b=b.strip(), row=pa if pa is not None and pa == pb else None,
                            pa=pa, pb=pb))
    return out


def check(rows):
    """Structural verdicts on the agreed rows, in N order: consecutive N, increasing sigma, sane eps and I."""
    good = sorted([r for r in rows if r["row"] is not None], key=lambda r: r["row"][0])
    for k, r in enumerate(good):
        n, s, e, i = r["row"]
        why = []
        if not 0.1 <= e <= 9.9: why.append("eps")
        if not 1 <= i <= 99: why.append("I")
        if k and n - good[k - 1]["row"][0] != 1: why.append(f"gap after {good[k-1]['row'][0]}")
        if k and not 0 < s - good[k - 1]["row"][1] < 3: why.append("sigma step")
        r["why"] = why
    return good


def words(img: Image.Image, scale=1.0, psm=6):
    """Tesseract word boxes on the whole page: (x, y, w, h, text, conf), in unscaled page pixels."""
    if scale != 1.0:
        img = img.resize((int(img.width * scale), int(img.height * scale)), Image.LANCZOS)
    buf = io.BytesIO(); img.save(buf, "PNG")
    p = subprocess.run(["tesseract", "stdin", "stdout", "--psm", str(psm), "-c", "tessedit_char_whitelist=0123456789. ",
                        "tsv"], input=buf.getvalue(), capture_output=True)
    out = []
    for line in p.stdout.decode("utf8", "replace").splitlines()[1:]:
        f = line.split("\t")
        if len(f) == 12 and f[11].strip():
            x, y, w, h = (int(v) / scale for v in f[6:10])
            out.append((x, y, w, h, f[11].strip(), float(f[10])))
    return out


TOKEN_SIGMA = re.compile(r"^1[1-4]\d{3}\.\d{4}$")
TOKEN_EPS = re.compile(r"^\d\.\d$")
TOKEN_INT = re.compile(r"^\d{1,4}$")


def records(ws, row_tol=18):
    """Assemble (N, sigma, eps, I) records from word boxes: a sigma token anchors each record; N is the
    integer just left of it on the same line, eps and I the next two tokens to its right."""
    ws = sorted(ws, key=lambda t: (t[1], t[0]))
    recs = []
    for k, (x, y, w, h, t, c) in enumerate(ws):
        if not TOKEN_SIGMA.match(t):
            continue
        yc = y + h / 2
        same = sorted([u for u in ws if abs(u[1] + u[3] / 2 - yc) < row_tol], key=lambda u: u[0])
        i = next(j for j, u in enumerate(same) if u[0] == x and u[4] == t)
        left = same[i - 1] if i > 0 and x - (same[i - 1][0] + same[i - 1][2]) < 3 * h else None
        right = [u for u in same[i + 1:i + 3]]
        n = int(left[4]) if left and TOKEN_INT.match(left[4]) else None
        e = float(right[0][4]) if len(right) > 0 and TOKEN_EPS.match(right[0][4]) else None
        i_ = int(right[1][4]) if len(right) > 1 and TOKEN_INT.match(right[1][4]) and len(right[1][4]) <= 2 else None
        recs.append(dict(N=n, sigma=float(t), eps=e, I=i_, box=(x, y, w, h)))
    return recs


def crops(path, rot, outdir, pad=40):
    """Write the table region of a photo as two full-resolution halves (left and right column groups),
    split at the widest vertical gap near the middle, for reading by eye."""
    im = load(path, rot)
    ink = glyphs_only(binarize(flatten(im)), min_size=18)
    ys, xs = np.nonzero(ink)
    # the table: rows and columns with a real amount of glyph ink (ignores stray marks and page numbers)
    rows = np.flatnonzero(ink.sum(1) > 8); cols = np.flatnonzero(ink.sum(0) > 3)
    y0, y1 = max(rows.min() - pad, 0), min(rows.max() + pad, im.height)
    x0, x1 = max(cols.min() - pad, 0), min(cols.max() + pad, im.width)
    prof = ink[y0:y1, x0:x1].sum(0)
    w = x1 - x0
    lo = int(0.3 * w)
    gaps = runs(prof[lo:int(0.7 * w)] == 0, 1)          # blank column runs; cut in the middle of the widest
    cut = x0 + (lo + (max(gaps, key=lambda g: g[1] - g[0])[0] + max(gaps, key=lambda g: g[1] - g[0])[1]) // 2
                if gaps else w // 2)
    stem = Path(path).stem.replace("PXL_20260925_", "")
    out = []
    for tag, (a, b) in (("L", (x0, cut)), ("R", (cut, x1))):
        f = Path(outdir) / f"{stem}_{tag}.png"
        im.crop((a, y0, b, y1)).save(f)
        out.append(f)
    return out


def main(argv):
    if argv[0] == "crops":
        for f in crops(argv[1], int(argv[2]), argv[3]):
            print(f)
        return
    if argv[0] == "words":
        path = argv[1]
        rot = int(next((a.split("=")[1] for a in argv if a.startswith("--rot=")), 0))
        im = flatten(load(path, rot))
        for scale in (1.0, 1.5):
            recs = records(words(im, scale))
            print(f"scale {scale}: {len(recs)} records")
            for r in recs:
                print("  ", r["N"], f"{r['sigma']:.4f}", r["eps"], r["I"])
        return
    if argv[0] == "page":
        path = argv[1]
        rot = int(next((a.split("=")[1] for a in argv if a.startswith("--rot=")), 0))
        rows = read_page(path, rot)
        good = check(rows)
        bad = [r for r in rows if r["row"] is None]
        for r in good:
            print("OK " if not r["why"] else "?  ", *r["row"], " ".join(r["why"]))
        for r in bad:
            print("XX ", r["box"], "|", r["a"], "|", r["b"])
        print(f"{len(rows)} text rows: {len(good)} agreed ({sum(1 for r in good if not r['why'])} clean), {len(bad)} to read")
        return
    if argv[0] == "orient":
        d = Path(argv[1])
        out = []
        for p in sorted(d.glob("PXL_*.jpg")):
            rot, s = orient(p)
            med = f"{np.median(s):.1f}" if s else "-"
            lo = f"{min(s):.4f}" if s else "-"
            hi = f"{max(s):.4f}" if s else "-"
            line = f"{p.name}\t{rot}\t{len(s)}\t{lo}\t{med}\t{hi}"
            print(line, flush=True)
            out.append(line)
        (d / "orient.tsv").write_text("photo\trot\tn_sigma\tmin\tmedian\tmax\n" + "\n".join(out) + "\n")


if __name__ == "__main__":
    main(sys.argv[1:])
