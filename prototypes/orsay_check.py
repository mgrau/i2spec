"""Check the hand-transcribed Orsay atlas pages (one TSV per table page: N, sigma, eps, I).

Every page is checked on its own and against its neighbours: N consecutive within and across pages,
sigma strictly increasing, eps in 0.1-9.9 or blank (the atlas prints *** where it has none), I a small
integer. Steps in sigma above 1.5 cm-1 or below 0.002 cm-1 are listed for a second look: they are legal
but are where a misread digit in the integer part or the first decimal would show.

usage: orsay_check.py PAGES_DIR
"""
import sys
from pathlib import Path


def load(p):
    rows = []
    for k, line in enumerate(Path(p).read_text().splitlines(), 1):
        f = line.split("\t")
        if len(f) != 4:
            raise SystemExit(f"{p}:{k}: expected 4 fields, got {len(f)}: {line!r}")
        rows.append((int(f[0]), float(f[1]), float(f[2]) if f[2] else None, int(f[3]), f[1]))
    return rows


def main(d):
    pages = sorted(Path(d).glob("p*.tsv"))
    prev, total, problems = None, 0, 0
    for p in pages:
        rows = load(p)
        total += len(rows)
        for a, b in zip([prev] + rows[:-1] if prev else rows[:-1], rows if prev else rows[1:]):
            if b[0] != a[0] + 1:
                print(f"{p.name}: N jumps {a[0]} -> {b[0]}"); problems += 1
            step = b[1] - a[1]
            if not step > 0:
                print(f"{p.name}: sigma not increasing at N={b[0]}: {a[1]:.4f} -> {b[1]:.4f}"); problems += 1
            elif step > 1.5 or step < 0.002:
                print(f"{p.name}: look again at N={b[0]}: step {step:.4f} cm-1 ({a[1]:.4f} -> {b[1]:.4f})")
        for n, s, e, i, raw in rows:
            if len(raw.split(".")[1]) != 4:
                print(f"{p.name}: N={n} sigma not 4 decimals: {raw}"); problems += 1
            if e is not None and not 0.1 <= e <= 9.9:
                print(f"{p.name}: N={n} eps {e}"); problems += 1
            if not 1 <= i <= 99:
                print(f"{p.name}: N={n} I {i}"); problems += 1
        print(f"{p.name}: N {rows[0][0]}-{rows[-1][0]} ({len(rows)}), sigma {rows[0][1]:.4f}-{rows[-1][1]:.4f}")
        prev = rows[-1]
    print(f"{len(pages)} pages, {total} lines, {problems} problems")


if __name__ == "__main__":
    main(sys.argv[1])
