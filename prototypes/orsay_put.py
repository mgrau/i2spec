"""Turn a typed page ("N sigma eps I" per line, *** for a missing eps) into the checked TSV format.
usage: orsay_put.py TYPED.txt OUT.tsv"""
import sys
out = []
for line in open(sys.argv[1]):
    f = line.split()
    if not f:
        continue
    if len(f) != 4:
        raise SystemExit(f"bad line: {line!r}")
    out.append("\t".join([f[0], f[1], "" if f[2] == "***" else f[2], f[3]]))
open(sys.argv[2], "w").write("\n".join(out) + "\n")
