"""A minimal reader for Bruker OPUS files (the .0, .1, ... files of the IFS 125HR and its kin).

An OPUS file starts with a directory: from byte 24, 12-byte entries (block type, channel type, text type, spare;
int32 length in 4-byte words; int32 byte offset), ended by an offset <= 0. Parameter blocks are lists of
(3-letter name + NUL, uint16 type, uint16 size in 2-byte words, value): type 0 int32, 1 float64, 2-4 a string.
A data block of block type 7 holds float32 values; its x axis and point count are in the "data status"
parameter block with the same channel type (FXV, LXV, NPT, CSF), of block type 23 in the IAG files (31 in
some other versions). Channel types: 4 the
single-channel sample spectrum (ScSm), 8 the sample interferogram (IgSm), 12 the sample phase (PhSm).

usage: opus.py FILE      list the blocks and the main parameters
"""
from __future__ import annotations

import struct
import sys
from pathlib import Path

import numpy as np

CHANNELS = {4: "ScSm", 8: "IgSm", 12: "PhSm", 132: "ScRf", 136: "IgRf", 140: "PhRf"}
STATUS = (23, 31)                     # block types of a data-status parameter block


def directory(raw):
    """[(block type, channel type, text type, length in words, offset)] of an OPUS file. The header gives the
    directory's offset and its number of entries (bytes 12-24); failing that, entries are read from byte 24
    until an offset <= 0."""
    first, most, n = struct.unpack("<iii", raw[12:24])
    if not (first == 24 and 0 < n <= most < 1000):
        first, n = 24, (len(raw) - 24) // 12
    out = []
    for k in range(n):
        c = first + 12 * k
        btype, ctype, ttype, _ = raw[c:c + 4]
        length, offset = struct.unpack("<ii", raw[c + 4:c + 12])
        if offset <= 0:
            break
        out.append((btype, ctype, ttype, length, offset))
    return out


def parameters(raw, offset):
    """{name: value} of the parameter block at offset."""
    out, i = {}, offset
    while i + 8 <= len(raw):
        name = raw[i:i + 3].decode("ascii", "replace")
        ptype, size = struct.unpack("<HH", raw[i + 4:i + 8])
        if name == "END":
            break
        value = raw[i + 8:i + 8 + 2 * size]
        if ptype == 0:
            out[name] = struct.unpack("<i", value[:4])[0]
        elif ptype == 1:
            out[name] = struct.unpack("<d", value[:8])[0]
        else:
            out[name] = value.split(b"\0", 1)[0].decode("latin-1")
        i += 8 + 2 * size
    return out


def read(path, channel="ScSm"):
    """(x, y, status, all parameter blocks) of one data block: x from FXV to LXV in NPT points (ascending),
    y scaled by CSF. channel: ScSm, IgSm or PhSm."""
    raw = Path(path).read_bytes()
    ch = {v: k for k, v in CHANNELS.items()}[channel]
    blocks = directory(raw)
    status = data = None
    params = {}
    for btype, ctype, ttype, length, offset in blocks:
        if btype in STATUS and ctype == ch:
            status = parameters(raw, offset)
        elif btype == 7 and ctype == ch:
            data = (offset, length)
        elif btype not in (7, 11, 15) and btype not in STATUS:
            try:
                params[(btype, ctype, ttype)] = parameters(raw, offset)
            except (struct.error, UnicodeDecodeError):
                pass
    if status is None or data is None:
        raise ValueError(f"{path}: no {channel} block (blocks: {[(b, c, t) for b, c, t, _, _ in blocks]})")
    n = int(status["NPT"])
    y = np.frombuffer(raw, dtype="<f4", count=n, offset=data[0]).astype(float) * float(status.get("CSF", 1.0))
    x = np.linspace(float(status["FXV"]), float(status["LXV"]), n)
    if x[0] > x[-1]:
        x, y = x[::-1], y[::-1]
    return x, y, status, params


def main(argv):
    raw = Path(argv[0]).read_bytes()
    for btype, ctype, ttype, length, offset in directory(raw):
        print(f"block type {btype:3d} channel {ctype:3d} ({CHANNELS.get(ctype, '?')}) text {ttype:3d}: {length:9d} words at {offset}")
    x, y, status, params = read(argv[0])
    print("ScSm status:", {k: status[k] for k in ("DPF", "NPT", "FXV", "LXV", "CSF", "DXU", "DAT", "TIM") if k in status})
    for key, p in params.items():
        keep = {k: p[k] for k in ("RES", "APF", "PHZ", "ZFF", "HFQ", "LFQ", "LWN", "NSS", "SRC", "BMS", "DTC", "APT", "VEL", "HFW", "LFW",
                                  "TSC", "INS", "DAT", "TIM", "SNM", "EXP") if k in p}
        if keep:
            print(key, keep)
    print(f"x {x[0]:.4f}-{x[-1]:.4f} cm-1, step {(x[-1] - x[0]) / (len(x) - 1):.6f}, y {np.nanmin(y):.4g}-{np.nanmax(y):.4g}")


if __name__ == "__main__":
    main(sys.argv[1:])
