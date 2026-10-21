"""Assign parts II-IV of the Orsay 11 000-14 000 cm-1 atlas (12 900-14 100 cm-1, 500 C) to model transitions.

Part I (790 C) is assigned by orsay_assign.py. Parts II-IV were recorded together, in the 30 cm cell at
500 C with the cold point at 60 C (5 torr), and their lower levels (v'' <= 17) are known to about 1 MHz
from comb data, so the same method works with less ambiguity:

1. the model line list at 773.15 K over each part's range, thinned per 30 cm-1 to 1.5 times the atlas's
   own line count there (orsay_assign.thin_range);
2. each atlas line takes the strongest model line within the window if it is DOMINANCE times stronger than
   any other there and no other atlas line claims it; three passes, 360, 180 and 120 MHz, the later ones
   after a smooth per-v'' offset in J(J+1) fitted to the previous pass;
3. the chance rate from a null that keeps the number of atlas lines per cm-1 but scrambles their positions
   within each cm-1;
4. each part's scale (ppb) fitted to the assigned lines whose X and B levels both lie inside the J ranges
   of the comb-referenced level corrections, as in orsay_calibration.py.

The overlaps II/III (13 350-13 392 cm-1) and III/IV (13 800-13 918 cm-1) are the same recording printed
twice (orsay-atlas-11000-14000.md); each printed line is kept once, from the lower part.

usage: orsay_assign_parts.py [--parameters=SET]
       -> data/atlas_lines/orsay1982_part{2,3,4}_assigned.csv and .toml
"""
from __future__ import annotations

import csv
import sys
from collections import Counter
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import orsay_assign as A                                                  # noqa: E402
from i2spec.constants import DEFAULT_PARAMETERS                           # noqa: E402
from i2spec.intensity import SHARED_GRID, X_RMIN, line_list               # noqa: E402
from i2spec.level_corrections import J_MARGIN, load_level_corrections     # noqa: E402
from i2spec.model import RovibronicModel                                  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
ATLAS = ROOT / "data" / "atlas_lines"
MHZ, T = 29979.2458, 773.15
PARTS = (2, 3, 4)
NAMES = {2: "II", 3: "III", 4: "IV"}
DUPLICATE_CM = 2e-4          # a line printed twice reads the same to the last digit; this allows one digit


def model_lines(pset, lo, hi, part):
    f = A.CACHE / f"orsay_ll_part{part}_{pset}.npz"
    if f.exists():
        return dict(np.load(f))
    g = dict(SHARED_GRID)
    n_inner = max(int(round((g["rmin"] - X_RMIN) / g["step"])), 0)
    gx = dict(g, rmin=g["rmin"] - n_inner * g["step"], nlev=48)
    m = RovibronicModel("127I2", pset, grids={"X": gx, "B": dict(g, nlev=70)}, solver="dvr")
    ll = line_list(m, T, lo, hi, S_min=1e-28)
    d = dict(nu=np.asarray(ll.nu), S=np.asarray(ll.S), vu=np.asarray(ll.v_upper), vl=np.asarray(ll.v_lower),
             J=np.asarray(ll.J_lower), br=np.asarray(ll.branch))
    o = np.argsort(d["nu"])
    d = {k: v[o] for k, v in d.items()}
    np.savez(f, **d)
    return d


def smooth_shift(M, obs, a):
    """A per-v'' offset quadratic (linear for few lines) in J(J+1), clipped at 3 robust sigma: the next pass's
    expected obs - model."""
    ok = np.array([x is not None for x in a])
    idx = np.array([x for x in a if x is not None], dtype=int)
    r = obs[ok] - M["nu"][idx]
    vl, y = M["vl"][idx], M["J"][idx] * (M["J"][idx] + 1) / 1e4
    coef = {}
    for v in np.unique(vl):
        m = vl == v
        if m.sum() >= 8:
            c = np.polyfit(y[m], r[m], 2 if m.sum() > 30 else 1)
            res = r[m] - np.polyval(c, y[m])
            good = np.abs(res) < 3 * 1.4826 * np.median(np.abs(res))
            coef[v] = np.polyfit(y[m][good], r[m][good], 2 if good.sum() > 30 else 1)

    def shift(ii, coef=coef):
        out = np.zeros(len(ii))
        for v, c in coef.items():
            m = M["vl"][ii] == v
            out[m] = np.polyval(c, M["J"][ii][m] * (M["J"][ii][m] + 1) / 1e4)
        return out
    return shift


def scramble(obs, rng):
    """The null: the same number of lines in each cm-1, at random positions within it."""
    base = np.floor(obs)
    return np.sort(base + rng.random(len(obs)))


def calibrate(obs, eps, M, a, C):
    """The part's scale (ppb) from the assigned lines whose levels the comb-referenced corrections cover."""
    def covered(st, v, J):
        if (st == "X" and v <= 10) or (st == "B" and v == 0):
            return True
        r = C.coverage[st].get(v)
        return r is not None and r[0] - J_MARGIN <= J <= r[1] + J_MARGIN
    sel = [(k, ai) for k, ai in enumerate(a) if ai is not None and M["vl"][ai] <= 17
           and covered("X", int(M["vl"][ai]), int(M["J"][ai]))
           and covered("B", int(M["vu"][ai]), int(M["J"][ai]) + (1 if M["br"][ai] > 0 else -1))]
    k = np.array([s[0] for s in sel]); ai = np.array([s[1] for s in sel])
    r = (obs[k] - M["nu"][ai]) * MHZ
    x = obs[k] * MHZ * 1e-9
    w = 1 / np.hypot(eps[k], 20.0) ** 2
    good = np.ones(len(r), bool)
    for _ in range(6):
        ppb = np.sum(w[good] * x[good] * r[good]) / np.sum(w[good] * x[good] ** 2)
        res = r - ppb * x
        s = 1.4826 * np.median(np.abs(res[good]))
        good = np.abs(res) < 3 * s
    err = 1 / np.sqrt(np.sum(w[good] * x[good] ** 2)) * np.sqrt(np.sum(w[good] * res[good] ** 2) / (good.sum() - 1))
    return ppb, err, int(good.sum()), s


def main(argv):
    opts = {a.split("=")[0]: a.split("=", 1)[1] for a in argv if "=" in a}
    pset = opts.get("--parameters", DEFAULT_PARAMETERS)
    C = load_level_corrections(__import__("json").loads(
        (ROOT / "src/i2spec/data" / f"{pset}.json").read_text())["level_corrections"])
    rng = np.random.default_rng(1)
    kept_sigmas = []                      # (part, sigma) of the lines kept so far, for the printed-twice overlaps
    summary = {}
    for part in PARTS:
        rows = list(csv.DictReader(open(ATLAS / f"orsay1982_part{part}.csv")))
        obs = np.array([float(r["sigma_cm1"]) for r in rows])
        eps = np.array([float(r["eps_mk"]) if r["eps_mk"] else 9.0 for r in rows]) * 1e-3 * MHZ
        M = A.thin_range(model_lines(pset, obs.min() - 2, obs.max() + 2, part), obs, A.STRONG_FACTOR)
        shift = None
        for w in (0.012, 0.006, 0.004):
            a = A.assign(obs, M, w, shift)
            shift = smooth_shift(M, obs, a)
        a = A.assign(obs, M, 0.004, shift)
        n_ok = sum(x is not None for x in a)
        nulls = [sum(x is not None for x in A.assign(scramble(obs, rng), M, 0.004, shift)) for _ in range(5)]
        ppb, err, n_anchor, mad = calibrate(obs, eps, M, a, C)
        # printed twice: a line of this part that reads the same as one already kept from a lower part
        prev = np.array([s for p, s in kept_sigmas if p < part])
        dup = np.zeros(len(obs), bool)
        if prev.size:
            j = np.searchsorted(prev, obs)
            for k, x in enumerate(obs):
                near = [prev[i] for i in (j[k] - 1, j[k]) if 0 <= i < len(prev)]
                dup[k] = any(abs(x - y) <= DUPLICATE_CM for y in near)
        kept_sigmas += [(part, x) for x in obs[~dup]]
        kept_sigmas.sort(key=lambda t: t[1])
        out_rows = []
        for rr, x, ai, d in zip(rows, obs, a, dup):
            if ai is None or d:
                continue
            out_rows.append([rr["N"], rr["sigma_cm1"], rr["eps_mk"], int(M["vu"][ai]), int(M["vl"][ai]), int(M["J"][ai]),
                             "R" if M["br"][ai] > 0 else "P", f"{M['nu'][ai]:.5f}", f"{(x - M['nu'][ai]) * MHZ:.1f}",
                             f"{M['S'][ai]:.3e}"])
        name = f"orsay1982_part{part}"
        with open(ATLAS / f"{name}_assigned.csv", "w", newline="") as fh:
            wr = csv.writer(fh)
            wr.writerow(["N", "sigma_cm1", "eps_mk", "v_upper", "v_lower", "J_lower", "branch", "model_cm1",
                         "obs_minus_model_MHz", "S_model"])
            wr.writerows(out_rows)
        r = np.array([float(x[8]) for x in out_rows])
        resid = r - ppb * 1e-9 * np.array([float(x[1]) for x in out_rows]) * MHZ
        vls = Counter(int(x[4]) for x in out_rows)
        vus = Counter(int(x[3]) for x in out_rows)
        summary[part] = dict(n=len(rows), assigned=n_ok, null=float(np.mean(nulls)), dup=int(dup.sum()),
                             kept=len(out_rows), ppb=ppb, err=err, anchors=n_anchor, mad=mad,
                             resid_med=float(np.median(resid)), resid_mad=float(1.4826 * np.median(np.abs(resid - np.median(resid)))),
                             vl=(min(vls), max(vls)), vu=(min(vus), max(vus)), lo=obs.min(), hi=obs.max())
        s = summary[part]
        print(f"part {NAMES[part]}: {s['n']} lines, {s['assigned']} assigned ({s['null']:.0f} by chance in the scrambled "
              f"null, {100 * s['null'] / s['assigned']:.1f} %), {s['dup']} printed twice, {s['kept']} kept; v' {s['vu'][0]}-{s['vu'][1]}, "
              f"v'' {s['vl'][0]}-{s['vl'][1]}; scale {ppb:+.1f} +- {err:.1f} ppb on {n_anchor} anchors (MAD {mad:.0f} MHz); "
              f"after it, obs - model median {s['resid_med']:+.1f} MHz, robust spread {s['resid_mad']:.0f} MHz", flush=True)
        write_toml(part, s, pset)
    return summary


def write_toml(part, s, pset):
    name = f"orsay1982_part{part}"
    (ATLAS / f"{name}.toml").write_text(f'''id = "{name}"
unit = "cm-1"
isotopologue = "127I2"
kind = "atlas line positions"
citation = """S. Gerstenkorn, J. Vergès, J. Chevillard, Atlas du spectre d'absorption de la molécule d'iode, 11 000-14 000 cm-1 (Laboratoire Aimé Cotton, CNRS II, Orsay, 1982), part {NAMES[part]}"""
url = "https://books.google.com/books/about/Atlas_du_spectre_d_absorption_de_la_mole.html?id=2Qy8QwAACAAJ"
instrument = """Fourier transform spectrometer, 30 cm cell at 500 °C with the cold point at 60 °C (5 torr); line table of wavenumber, precision (mK) and depth, {s['lo']:.1f}-{s['hi']:.1f} cm-1"""
calibration = """The atlas's own wavenumber scale, which reads {s['ppb']:+.1f} +/- {s['err']:.1f} ppb against the model (scale_ppb: observed = true (1 + scale_ppb 1e-9)), fitted to {s['anchors']} assigned lines whose X and B levels lie inside the J ranges of the comb-referenced level corrections (prototypes/orsay_assign_parts.py). The web app divides it out; after that the lines read {s['resid_med']:+.1f} MHz (median) from the model, robust spread {s['resid_mad']:.0f} MHz"""
scale_ppb = {s['ppb']:.1f}
extraction = """Transcribed by eye from photographs of the printed table, with structural checks (consecutive line numbers, monotonic wavenumbers) and the overlaps between parts; assigned against the model line list at 773 K thinned to the atlas's line density, 120 MHz window, a strongest candidate 3x any other (prototypes/orsay_assign_parts.py, with the production model of the time). {s['assigned']} of {s['n']} lines assigned, about {s['null']:.0f} ({100 * s['null'] / s['assigned']:.1f} %) by chance in a null that scrambles positions within each cm-1. Lines printed twice in the overlap with a lower part are kept once, from the lower part ({s['dup']} here). docs/research/orsay-atlas-11000-14000.md"""
n_lines = {s['kept']}
''')


if __name__ == "__main__":
    main(sys.argv[1:])
