"""Assign the Orsay part I lines (data/atlas_lines/orsay1982_part1.csv) to model transitions.

1. Model line list at the cell temperature, 1063 K (cached in the scratchpad).
2. The model list is thinned to lines strong enough to be in the atlas: per 30 cm-1 window, the
   STRONG_FACTOR x N strongest, N the number of atlas lines there. Unthinned, the list is so dense
   (200 lines per cm-1) that shifted copies of the atlas get three quarters as many "assignments" as
   the real one.
   Each atlas line takes the strongest model line within WINDOW; the assignment is kept only if that
   line is DOMINANCE times stronger than any other model line in the window (a blend is not an
   assignment) and no other atlas line claims the same model line.
3. Independent check, not using the B state: two lines from the same upper level (v', J') to the same
   v'' - R(J'-1) and P(J'+1) - differ by X(v'', J'+1) - X(v'', J'-1). Their observed difference against
   the model's tests both assignments at once.

usage: orsay_assign.py [--window=MHz] [--parameters=SET]  -> data/atlas_lines/orsay1982_part1_assigned.csv
"""
import csv, sys
from pathlib import Path
import numpy as np
from i2spec.model import RovibronicModel
from i2spec.intensity import line_list, SHARED_GRID, X_RMIN

MHZ, T = 29979.2458, 1063.15
DOMINANCE = 3.0
STRONG_FACTOR = 1.5
CACHE = Path(__file__).resolve().parent / "out"          # model line lists, cached between runs (not committed)
CACHE.mkdir(exist_ok=True)


def model_lines(pset):
    f = CACHE / f"orsay_ll_{pset}.npz"
    if f.exists():
        return dict(np.load(f))
    g = dict(SHARED_GRID)
    n_inner = max(int(round((g["rmin"] - X_RMIN) / g["step"])), 0)
    gx = dict(g, rmin=g["rmin"] - n_inner * g["step"], nlev=48)
    m = RovibronicModel("127I2", pset, grids={"X": gx, "B": dict(g, nlev=70)}, solver="dvr")
    ll = line_list(m, T, 10995, 13015, S_min=1e-28)
    d = dict(nu=np.asarray(ll.nu), S=np.asarray(ll.S), vu=np.asarray(ll.v_upper), vl=np.asarray(ll.v_lower),
             J=np.asarray(ll.J_lower), br=np.asarray(ll.branch))
    o = np.argsort(d["nu"]); d = {k: v[o] for k, v in d.items()}
    np.savez(f, **d)
    return d


def thin(M, obs, factor=STRONG_FACTOR, W=30.0):
    keep = np.zeros(len(M["nu"]), bool)
    for lo in np.arange(10995.0, 13015.0, W):
        n = int(np.sum((obs >= lo) & (obs < lo + W)))
        sel = np.flatnonzero((M["nu"] >= lo) & (M["nu"] < lo + W))
        keep[sel[np.argsort(M["S"][sel])[::-1][: int(factor * n)]]] = True
    return {k: v[keep] for k, v in M.items()}


def thin_range(M, obs, factor=STRONG_FACTOR, W=30.0):
    """thin() over the observed range, whatever it is (the parts beyond part I)."""
    keep = np.zeros(len(M["nu"]), bool)
    for lo in np.arange(obs.min() - 1, obs.max() + 1, W):
        n = int(np.sum((obs >= lo) & (obs < lo + W)))
        sel = np.flatnonzero((M["nu"] >= lo) & (M["nu"] < lo + W))
        keep[sel[np.argsort(M["S"][sel])[::-1][: int(factor * n)]]] = True
    return {k: v[keep] for k, v in M.items()}


def assign(obs, M, window_cm, shift=None):
    """shift: optional callable(model index array) -> expected obs-model offset in cm-1 (iteration)."""
    nu = M["nu"] + (shift(np.arange(len(M["nu"]))) if shift else 0.0)
    order = np.argsort(nu); nus = nu[order]
    out = []
    for k, x in enumerate(obs):
        lo, hi = np.searchsorted(nus, [x - window_cm, x + window_cm])
        cand = order[lo:hi]
        if len(cand) == 0:
            out.append(None); continue
        s = M["S"][cand]; best = cand[np.argmax(s)]
        others = np.sort(s)[::-1]
        if len(others) > 1 and others[0] < DOMINANCE * others[1]:
            out.append(None); continue
        out.append(int(best))
    # a model line claimed by two atlas lines is assigned to neither
    ids = [a for a in out if a is not None]
    dup = {a for a in ids if ids.count(a) > 1}
    return [None if a in dup else a for a in out]


def main(argv):
    opts = {a.split("=")[0]: a.split("=", 1)[1] for a in argv if "=" in a}
    pset = opts.get("--parameters", "i2spec2026f")
    rows = list(csv.DictReader(open("data/atlas_lines/orsay1982_part1.csv")))
    obs = np.array([float(r["sigma_cm1"]) for r in rows])
    M = thin(model_lines(pset), obs, float(opts.get("--strong", STRONG_FACTOR)))
    print(f"{len(M['nu'])} model lines kept; {len(obs)} atlas lines")
    shift = None
    for it, w in enumerate((0.012, 0.006, 0.004)):
        a = assign(obs, M, w, shift)
        ok = np.array([x is not None for x in a])
        idx = np.array([x for x in a if x is not None])
        r = obs[ok] - M["nu"][idx]
        # per-v'' smooth offset in J(J+1) for the next pass
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
        print(f"pass {it}: window {w*MHZ:.0f} MHz, assigned {ok.sum()}/{len(obs)}, "
              f"median |obs-model| {np.median(np.abs(r))*MHZ:.0f} MHz")
    # final residuals after the per-level smooth shift
    a = assign(obs, M, 0.004, shift)
    with open("data/atlas_lines/orsay1982_part1_assigned.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["N", "sigma_cm1", "eps_mk", "v_upper", "v_lower", "J_lower", "branch", "model_cm1", "obs_minus_model_MHz", "S_model"])
        for rr, x, ai in zip(rows, obs, a):
            if ai is None:
                continue
            w.writerow([rr["N"], rr["sigma_cm1"], rr["eps_mk"], int(M["vu"][ai]), int(M["vl"][ai]), int(M["J"][ai]),
                        "R" if M["br"][ai] > 0 else "P", f"{M['nu'][ai]:.5f}", f"{(x - M['nu'][ai])*MHZ:.1f}", f"{M['S'][ai]:.3e}"])
    A = list(csv.DictReader(open("data/atlas_lines/orsay1982_part1_assigned.csv")))
    print(f"final: {len(A)} assigned")
    from collections import Counter
    c = Counter(int(x["v_lower"]) for x in A)
    for v in sorted(c):
        rv = np.array([float(x["obs_minus_model_MHz"]) for x in A if int(x["v_lower"]) == v])
        Js = [int(x["J_lower"]) for x in A if int(x["v_lower"]) == v]
        print(f"   v''={v:2d}: {c[v]:4d} lines, J'' {min(Js):3d}-{max(Js):3d}, obs-model median {np.median(rv):+7.0f} MHz, spread {1.4826*np.median(np.abs(rv-np.median(rv))):5.0f}")
    # combination differences: same (v', J') upper level, same v'', R(J'-1) and P(J'+1)
    key = {}
    for x in A:
        Jl, br = int(x["J_lower"]), x["branch"]
        Ju = Jl + 1 if br == "R" else Jl - 1
        key.setdefault((int(x["v_upper"]), Ju, int(x["v_lower"])), {})[br] = x
    cd = []
    for (vu, Ju, vl), d in key.items():
        if "R" in d and "P" in d:
            o = float(d["R"]["sigma_cm1"]) - float(d["P"]["sigma_cm1"])
            m = float(d["R"]["model_cm1"]) - float(d["P"]["model_cm1"])
            cd.append((vl, Ju, (o - m) * MHZ))
    cd = np.array(cd)
    print(f"combination differences (R-P from one upper level): {len(cd)} pairs; "
          f"|obs-model| median {np.median(np.abs(cd[:,2])):.0f} MHz, "
          f"{np.mean(np.abs(cd[:,2]) < 150)*100:.0f}% within 150 MHz")
    np.save(CACHE / "orsay_cd.npy", cd)


if __name__ == "__main__":
    main(sys.argv[1:])
