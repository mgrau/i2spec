"""Held-out bake-off of three level-correction schemes on the comb-referenced 127I2 rows.

usage: uv run --group research python prototypes/bakeoff.py [rows] [run] [figure]
       (no argument: all three stages; `rows` builds the cache, `run` the folds, `figure` the plot)

Rows: those of prototypes/level_corrections_fit.py (127I2 absolute frequencies and inter-line intervals,
r = observed - bare potentials of i2spec2026n, sigma_eff = hypot(sigma, 0.3 MHz)), restricted to
B v' <= 50 and X v'' <= 17. Each row is +B(v', J') - X(v'', J'') (absolute) or the difference of two such
(interval).

Schemes
  A  per-level polynomials in y = J(J+1)/1e4, exactly the fit of level_corrections_fit.py (degree rule,
     ridge prior, robust reweighting, delta X = 0 below v'' = 11, delta B(0) = 0); a level with no
     parameters in the training fit gets zero. Held-out sigma: the lookup rule (corrected levels within
     15 in J: 0.3 MHz floor + per-level covariance/discrepancy or held-out rms in quadrature; else the
     region value of lookup.uncertainty).
  B  Gaussian-process correction surfaces f_B(v, y), f_X(v, y): Matern-5/2 in v x squared-exponential in y,
     plus a per-v independent SE-in-y term, plus a white row term tau; the posterior is exact linear algebra
     on the rows (linear functionals), hyperparameters by type-II maximum likelihood on the training rows.
  C  linear potential correction dV_s(R) = sum_k c_k B_k(R) (cubic B-splines over the data's R range),
     plus a centrifugal term dq_s(R) J(J+1) hbar^2/(2 mu R^2) (fewer splines); the MLR-built B levels
     (v' >= 44) get their own dV = dV_B + dV_ext. First-order (Hellmann-Feynman) level shifts
     <psi_vJ| dV |psi_vJ>, ridge + second-difference penalty with the strengths chosen by set-blocked
     cross-validation inside each training set.

Splits: leave-one-data-set-out; leave-one-level-out (every level scheme A corrects that >= 2 rows touch);
the prediction test of paper/analysis/predtest_2026k.py (train on the sets available before 2026-09-29,
predict the absolute rows of the sets transcribed later).

Outputs: prototypes/out/bakeoff_rows.pkl (cache: rows and level densities), prototypes/out/bakeoff.json,
docs/figures/bakeoff.png. Notes in docs/research/model-bakeoff.md.
"""
from __future__ import annotations

import json
import os
import pickle
import sys
import time
from collections import defaultdict
from pathlib import Path

for _k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_k, "2")

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "prototypes/out"
CACHE = OUT / "bakeoff_rows.pkl"
BASE = "i2spec2026n"
FLOOR = 0.3
V_MAX_B, V_MAX_X = 50, 17          # the bake-off region
V_EXT_B = 44                       # from here the B levels come from the MLR (i2spec2026n "extended")
#: density bins for <psi|f(R)|psi>: f is integrated as sum_i m_i f(R_i) with m_i the psi^2 mass of bin i
R_LO, R_HI, R_BIN = 2.10, 8.00, 0.002
#: the unmeasured levels whose predictions are compared (B: the line R(J'') v'-0; X: X(v'', J) - X(0, J))
SAMPLE_B, SAMPLE_X = (7, 12, 19, 27, 38), (5, 8)
SAMPLE_J = tuple(range(0, 201, 10))
#: paper/analysis/predtest_2026k.py: the sets transcribed after i2spec2026k was fixed
NEW = ["huang2013a", "hong2009a", "ye1999a", "holzwarth2001a", "jones2002a", "goncharov2004a", "simonsen2000a",
       "huang2018a", "manzoor2024a", "fan2014a", "kobayashi2015a", "hauden2024a", "grieser1994a", "badr2006a",
       "hong2001b", "sakagami2020a", "yoshii2019a", "tanabe2022a", "arie1993a", "arie1994a", "cheng2001a",
       "zhang2001a", "hong2001a", "hong2002a", "hong2000a", "hong2004a", "hong2004b", "hong2002b", "fang2006a"]


# --------------------------------------------------------------------------------------------- rows

def build_rows():
    """The rows of level_corrections_fit.py (its loop, verbatim in substance) and the psi^2 of every level."""
    from level_corrections_fit import BarePredictor
    from i2spec.constants import MHZ_PER_CM
    from i2spec.observations import load_all

    t0 = time.time()
    pred = BarePredictor(BASE)
    rows = []
    for ds, raw in zip(load_all(), load_all(raw=True)):
        assert ds.id == raw.id and len(ds) == len(raw), ds.id
        scale = MHZ_PER_CM if ds.unit == "cm-1" else 1.0
        for o, o_raw in zip(ds.observations, raw.observations):
            ref = o.ref_line or o.line
            if o.line.isotopologue != "127I2" or ref.isotopologue != "127I2":
                continue
            if o.kind == "interval" and ref == o.line:
                continue
            if max(o.line.v_upper, ref.v_upper) > V_MAX_B or max(o.line.v_lower, ref.v_lower) > V_MAX_X:
                continue
            r = (o.value - pred(o, ds.unit)) * scale
            levels = defaultdict(float)
            lines = [(o.line, +1.0)] + ([(ref, -1.0)] if o.kind == "interval" else [])
            for line, sign in lines:
                Ju = line.J_lower + (1 if line.branch == "R" else -1)
                levels[("B", line.v_upper, Ju)] += sign
                levels[("X", line.v_lower, line.J_lower)] -= sign
            levels = {k: s for k, s in levels.items() if s != 0}
            nus = [pred.position(ln) / MHZ_PER_CM for ln, _ in lines]
            rows.append(dict(set=ds.id, kind=o.kind, line=str(o.line), ref=str(ref) if o.kind == "interval" else None,
                             comp=o.component, r=float(r), sig_eff=float(np.hypot(o.uncertainty * scale, FLOOR)),
                             sig=float(o_raw.uncertainty * scale), sig_full=float(o.uncertainty * scale),
                             v_upper=o.line.v_upper, v_lower=o.line.v_lower, J_lower=o.line.J_lower,
                             branch=o.line.branch,
                             lines=[(ln.branch, ln.J_lower, ln.v_upper, ln.v_lower, float(nu)) for (ln, _), nu in zip(lines, nus)],
                             levels=levels))
    print(f"{len(rows)} rows in B v' <= {V_MAX_B}, X v'' <= {V_MAX_X} ({time.time() - t0:.0f} s)", flush=True)

    # psi^2 of every level the rows touch, plus the sample levels, binned in R
    need = {k for r in rows for k in r["levels"]}
    need |= {("B", v, J + 1) for v in SAMPLE_B for J in SAMPLE_J} | {("X", v, J) for v in (0,) + SAMPLE_X for J in SAMPLE_J}
    need |= {("B", v, J + 1) for v in range(0, V_MAX_B + 1) for J in (0, 100, 200)}
    need |= {("X", v, J) for v in range(0, V_MAX_X + 1) for J in (0, 100, 200)}
    m = pred.model("127I2")
    edges = np.arange(R_LO, R_HI + R_BIN / 2, R_BIN)
    dens = {}
    by_sJ = defaultdict(list)
    for st, v, J in need:
        by_sJ[(st, J)].append(v)
    for (st, J), vs in sorted(by_sJ.items()):
        for src in ("pub", "ext"):
            vv = [v for v in vs if (st == "B" and v >= V_EXT_B) == (src == "ext")]
            if not vv:
                continue
            solver = m.states[st] if src == "pub" else m.extended[st]
            E, psi = solver.wavefunctions(J)
            idx = np.clip(np.searchsorted(edges, solver.R) - 1, 0, len(edges) - 2)
            for v in vv:
                w = solver.W * psi[:, v] ** 2
                dens[(st, v, J)] = (np.bincount(idx, w, len(edges) - 1).astype(np.float32),
                                    np.bincount(idx, w / solver.R**2, len(edges) - 1).astype(np.float32))
    print(f"{len(dens)} level densities ({time.time() - t0:.0f} s)", flush=True)
    OUT.mkdir(exist_ok=True)
    with open(CACHE, "wb") as f:
        pickle.dump(dict(rows=rows, dens=dens, edges=edges, mu=m._mu), f)


# --------------------------------------------------------------------------------------------- data

class Data:
    """Rows, the level incidence matrix, and the sample queries, from the cache."""

    def __init__(self):
        import scipy.sparse as sp
        d = pickle.load(open(CACHE, "rb"))
        self.rows, self.dens, self.edges, self.mu = d["rows"], d["dens"], d["edges"], d["mu"]
        R = self.rows
        self.N = len(R)
        self.r = np.array([x["r"] for x in R])
        self.sig_eff = np.array([x["sig_eff"] for x in R])
        self.sig = np.array([x["sig"] for x in R])
        self.sig_full = np.array([x["sig_full"] for x in R])
        self.sets = np.array([x["set"] for x in R])
        self.kind = np.array([x["kind"] for x in R])
        self.line = np.array([x["line"] for x in R])
        self.precise = (self.sig <= 0.3) & (self.sig_full <= 0.5)
        # the sample queries: name -> {level: sign}
        self.queries, self.qmeta = {}, []
        for v in SAMPLE_B:
            for J in SAMPLE_J:
                self.queries[f"B{v} R({J})"] = {("B", v, J + 1): 1.0, ("X", 0, J): -1.0}
                self.qmeta.append(("B", v, J))
        for v in SAMPLE_X:
            for J in SAMPLE_J:
                self.queries[f"X{v} J={J}"] = {("X", v, J): 1.0, ("X", 0, J): -1.0}
                self.qmeta.append(("X", v, J))
        keys = sorted({k for x in R for k in x["levels"]} | {k for q in self.queries.values() for k in q})
        self.levels = keys
        self.lidx = {k: i for i, k in enumerate(keys)}
        self.state = np.array([k[0] for k in keys])
        self.v = np.array([k[1] for k in keys])
        self.J = np.array([k[2] for k in keys])
        self.y = self.J * (self.J + 1) / 1e4

        def mat(dicts):
            ii, jj, vv = [], [], []
            for i, lv in enumerate(dicts):
                for k, s in lv.items():
                    ii.append(i); jj.append(self.lidx[k]); vv.append(s)
            return sp.csr_matrix((vv, (ii, jj)), shape=(len(dicts), len(keys)))
        self.S = mat([x["levels"] for x in R])
        self.qnames = list(self.queries)
        self.Q = mat([self.queries[n] for n in self.qnames])
        # which (state, v) each row touches
        self.touch = [{(k[0], k[1]) for k in x["levels"]} for x in R]


# --------------------------------------------------------------------------------------------- scheme A

def region_sigma(branch, J_lower, v_upper, v_lower, nu):
    """lookup.uncertainty's fallback for a 127I2 line with no (covering) level correction, v' <= 50, v'' <= 17."""
    if v_upper > 43:
        return 15.0
    if nu >= 15000.0:
        return 5.0 if v_upper >= 31 else 3.0
    if nu >= 13250.0:
        return 5.0
    if nu >= 12270.0:
        return 3.0 if v_upper == 0 else 8.0
    return 50.0


class SchemeA:
    """level_corrections_fit.py's per-level polynomial fit, on a training mask."""

    name = "A"

    def __init__(self, D: Data):
        import level_corrections_fit as lcf
        self.D, self.lcf = D, lcf

    def fit(self, train):
        lcf, D = self.lcf, self.D
        rows = D.rows
        precise = defaultdict(list)
        for i in np.flatnonzero(train & D.precise):
            for (st, v, J), s in rows[i]["levels"].items():
                if (st == "X" and v >= lcf.X_MIN_V) or (st == "B" and v >= 1):
                    precise[(st, v)].append(J)
        deg = {k: lcf.degree(Js) for k, Js in precise.items()}
        names = sorted({(st, v, k) for (st, v), d in deg.items() for k in range(d + 1)}, key=str)
        idx = {n: i for i, n in enumerate(names)}
        A = np.zeros((D.N, len(names)))
        for i, x in enumerate(rows):
            for (st, v, J), s in x["levels"].items():
                if (st, v) in deg:
                    y = J * (J + 1) / 1e4
                    for k in range(deg[(st, v)] + 1):
                        A[i, idx[(st, v, k)]] += s * y**k
        r, w = D.r, 1 / D.sig_eff

        def prior_of(n):
            st, v, _ = n
            if st == "B" and v > lcf.V_PERTURBED:
                return 1e5
            if (st == "B" and v > lcf.V_MAX_B) or (st == "X" and v > lcf.V_MAX_X):
                return 100.0
            return 5.0
        ridge = np.diag([1 / prior_of(n) ** 2 for n in names])

        def fit(mask, cov=False):
            ww = w.copy()
            for _ in range(8):
                Aw, rw = A[mask] * ww[mask, None], r[mask] * ww[mask]
                c = np.linalg.solve(Aw.T @ Aw + ridge, Aw.T @ rw)
                z = np.abs((r - A @ c) * w)
                ww = w * np.minimum(1.0, 3.0 / np.maximum(z, 1e-9))
            if cov:
                Aw = A[mask] * ww[mask, None]
                C = np.linalg.inv(Aw.T @ Aw + ridge)
                dof = max(int(mask.sum()) - A.shape[1], 1)
                chi2 = float(np.sum(((r[mask] - A[mask] @ c) * ww[mask]) ** 2)) / dof
                return c, C * max(chi2, 1.0)
            return c
        c_all, C_all = fit(train, cov=True)
        used = A.any(axis=1) & train
        loo = r - A @ c_all
        for ln in set(D.line[used]):
            m = (D.line == ln) & train
            loo[m] = r[m] - A[m] @ fit(train & ~m)
        sig_full = D.sig_full

        def partner_ok(lv, st):
            return all((k0 != "X" or k1 <= lcf.V_MAX_X) if st == "B" else (k0 != "B" or k1 <= lcf.V_MAX_B)
                       for (k0, k1, _), s_ in lv.items() if s_ != 0)
        lev = {}
        for (st, v), d in deg.items():
            touch = np.array([(st, v) in D.touch[i] and partner_ok(rows[i]["levels"], st) for i in range(D.N)]) \
                & train & D.precise
            ids = [idx[(st, v, k)] for k in range(d + 1)]
            cov_v = C_all[np.ix_(ids, ids)]
            n_lines = len(set(D.line[touch]))
            md = None
            if n_lines >= 2:
                pv = np.array([sum(sg * (J * (J + 1) / 1e4) ** k * (J * (J + 1) / 1e4) ** l * cov_v[k, l]
                                   for (k0, k1, J), sg in rows[i]["levels"].items() if (k0, k1) == (st, v)
                                   for k in range(d + 1) for l in range(d + 1)) for i in np.flatnonzero(touch)])
                rr = float(np.sqrt(np.mean(loo[touch] ** 2)))
                md = float(np.sqrt(max(rr ** 2 - np.mean(sig_full[touch] ** 2 + np.abs(pv)), 0.0)))
            Js = precise[(st, v)]
            lev[(st, v)] = dict(c=np.array([c_all[i] for i in ids]), cov=cov_v, md=md or 0.0,
                                cover=(min(Js), max(Js)))
        self.lev = lev
        return self

    # the model's shift (LevelCorrections.shift with clamp_J) and lookup's per-level sigma
    def shift(self, st, v, J):
        L = self.lev.get((st, v))
        if L is None:
            return 0.0
        return float(self._phi(L, J) @ L["c"])

    @staticmethod
    def _phi(L, J):
        n = len(L["c"])
        y = J * (J + 1) / 1e4
        lo, hi = L["cover"]
        Jc = min(max(J, max(lo - 15, 0)), hi + 15)
        yc = Jc * (Jc + 1) / 1e4 if n > 2 else y
        return np.array([y ** k if k < 2 else yc ** k for k in range(n)])

    def covers(self, st, v, J):
        L = self.lev.get((st, v))
        return L is not None and L["cover"][0] - 15 <= J <= L["cover"][1] + 15

    def line_sigma(self, branch, J, vu, vl, nu):
        Ju = J + (1 if branch == "R" else -1)
        x_ok = vl <= 10 or self.covers("X", vl, J)
        b_ok = vu == 0 or self.covers("B", vu, Ju)
        if not (x_ok and b_ok):
            return region_sigma(branch, J, vu, vl, nu)
        parts = [0.3]
        for st, v, JJ in (("X", vl, J), ("B", vu, Ju)):
            if (st == "X" and v <= 10) or (st == "B" and v == 0):
                continue
            L = self.lev[(st, v)]
            ph = self._phi(L, JJ)
            parts.append(np.sqrt(max(ph @ L["cov"] @ ph, 0.0) + L["md"] ** 2))
        return float(np.sqrt(np.sum(np.square(parts))))

    def predict(self, idx):
        D = self.D
        mu = np.array([sum(s * self.shift(*k) for k, s in D.rows[i]["levels"].items()) for i in idx])
        sm = np.array([np.sqrt(sum(self.line_sigma(*ln) ** 2 for ln in D.rows[i]["lines"])) for i in idx])
        return mu, np.hypot(sm, D.sig_full[idx])

    def predict_queries(self):
        D = self.D
        mu = np.array([sum(s * self.shift(*k) for k, s in D.queries[n].items()) for n in D.qnames])
        sd = []
        for st, v, J in D.qmeta:            # B: the line R(J) v'-0 (visible); X: the level against X(0), no line
            sd.append(self.line_sigma("R", J, v, 0, 16000.0) if st == "B" else
                      3.0 if not self.covers("X", v, J) else self.line_sigma("R", J, 0, v, 16000.0))
        return mu, np.array(sd)


# --------------------------------------------------------------------------------------------- scheme B

class SchemeB:
    """Gaussian-process correction surfaces for B and X; rows are linear functionals of them."""

    name = "B"
    #: log-hyperparameters: per state (a1, l_v, l_y, a2, l_y2), then tau
    NAMES = [f"{s}_{p}" for s in ("B", "X") for p in ("amp", "len_v", "len_y", "amp_v", "len_y_v")] + ["tau"]
    LO = np.log([1e-3, 0.5, 0.05, 1e-3, 0.05] * 2 + [1e-4])
    HI = np.log([300, 100, 30, 300, 30] * 2 + [10])
    TH0 = np.log([3.0, 5.0, 1.0, 1.0, 1.0] * 2 + [0.1])

    def __init__(self, D: Data):
        self.D = D
        self.st = {}
        for s in ("B", "X"):
            ii = np.flatnonzero(D.state == s)
            v, y = D.v[ii].astype(float), D.y[ii]
            dv = np.abs(v[:, None] - v[None, :])
            self.st[s] = dict(ii=ii, dv=dv, same=(dv == 0).astype(float), dy2=(y[:, None] - y[None, :]) ** 2,
                              S=D.S[:, ii].tocsr(), Q=D.Q[:, ii].tocsr())
        self.theta = self.TH0.copy()

    def kern(self, th, s, grad=False):
        P = self.st[s]
        o = 0 if s == "B" else 5
        a1, lv, ly, a2, ly2 = np.exp(th[o:o + 5])
        q = np.sqrt(5.0) * P["dv"] / lv
        m52 = (1 + q + q * q / 3) * np.exp(-q)
        se = np.exp(-P["dy2"] / (2 * ly * ly))
        se2 = np.exp(-P["dy2"] / (2 * ly2 * ly2))
        k1, k2 = a1 * a1 * m52 * se, a2 * a2 * P["same"] * se2
        K = k1 + k2
        if not grad:
            return K
        dm52 = (q * q / 3) * (1 + q) * np.exp(-q)
        dK = [2 * k1, a1 * a1 * dm52 * se, k1 * P["dy2"] / (ly * ly), 2 * k2, k2 * P["dy2"] / (ly2 * ly2)]
        return K, dK

    def _rowK(self, tr, th, grad=False):
        K = np.zeros((tr.size, tr.size))
        parts = {}
        for s in ("B", "X"):
            S = self.st[s]["S"][tr]
            out = self.kern(th, s, grad)
            Ks = out[0] if grad else out
            K += np.asarray(S @ (S @ Ks).T)          # S K_s S^T (K_s symmetric)
            parts[s] = (S, out)
        return K, parts

    def nlml(self, th, tr, noise2):
        from scipy.linalg import cho_factor, cho_solve
        K, parts = self._rowK(tr, th, grad=True)
        tau2 = np.exp(2 * th[-1])
        K[np.diag_indices_from(K)] += noise2 + tau2
        try:
            cf = cho_factor(K, lower=True)
        except np.linalg.LinAlgError:
            return 1e10, np.zeros_like(th)
        r = self.D.r[tr]
        al = cho_solve(cf, r)
        val = 0.5 * r @ al + np.sum(np.log(np.diag(cf[0]))) + 0.5 * tr.size * np.log(2 * np.pi)
        M = np.outer(al, al) - cho_solve(cf, np.eye(tr.size))
        g = np.zeros_like(th)
        for s, o in (("B", 0), ("X", 5)):
            S, (_, dK) = parts[s]
            G = np.asarray((S.T @ (S.T @ M).T).T)         # S^T M S
            for k in range(5):
                g[o + k] = -0.5 * np.sum(G * dK[k])
        g[-1] = -0.5 * np.trace(M) * 2 * tau2
        return float(val), g

    def fit(self, train, theta0=None):
        from scipy.optimize import minimize
        from scipy.linalg import cho_factor, cho_solve
        tr = np.flatnonzero(train)
        th = self.theta.copy() if theta0 is None else theta0.copy()
        sig = self.D.sig_eff[tr].copy()
        w = sig.copy()
        for it in range(4):
            if it in (0, 3):
                res = minimize(self.nlml, th, args=(tr, w ** 2), jac=True, method="L-BFGS-B",
                               bounds=list(zip(self.LO, self.HI)), options=dict(maxiter=200 if it == 0 else 60))
                th = res.x
            K, _ = self._rowK(tr, th)
            Kn = K.copy()
            Kn[np.diag_indices_from(Kn)] += w ** 2 + np.exp(2 * th[-1])
            cf = cho_factor(Kn, lower=True)
            al = cho_solve(cf, self.D.r[tr])
            z = np.abs((self.D.r[tr] - K @ al) / sig)
            w = sig * np.maximum(1.0, z / 3.0)            # A's reweighting: weight x min(1, 3/|z|)
        self.th, self.tr, self.cf, self.al = th, tr, cf, al
        self.nll = float(res.fun)
        return self

    def _post(self, Sq_of):
        from scipy.linalg import cho_solve
        Kqt, kqq = 0.0, 0.0
        for s in ("B", "X"):
            P = self.st[s]
            Ks = self.kern(self.th, s)
            Sq, St = Sq_of(P), P["S"][self.tr]
            A = np.asarray((Sq @ Ks))                     # nq x L
            Kqt = Kqt + np.asarray((St @ A.T).T)
            kqq = kqq + np.asarray((Sq.multiply(A)).sum(axis=1)).ravel()
        mu = Kqt @ self.al
        var = kqq - np.sum(Kqt * cho_solve(self.cf, Kqt.T).T, axis=1)
        return mu, np.sqrt(np.maximum(var, 0.0))

    def predict(self, idx):
        mu, sd = self._post(lambda P: P["S"][idx])
        return mu, np.sqrt(sd ** 2 + self.D.sig_eff[idx] ** 2 + np.exp(2 * self.th[-1]))

    def predict_queries(self):
        return self._post(lambda P: P["Q"])


# --------------------------------------------------------------------------------------------- scheme C

def bspline_basis(x, lo, hi, n, k=3):
    """n clamped B-splines of degree k on [lo, hi] at x (zero outside)."""
    from scipy.interpolate import BSpline
    inner = np.linspace(lo, hi, n - k + 1)
    t = np.r_[[lo] * k, inner, [hi] * k]
    out = np.zeros((x.size, n))
    for i in range(n):
        c = np.zeros(n); c[i] = 1.0
        out[:, i] = np.nan_to_num(BSpline(t, c, k, extrapolate=False)(x))
    return out


class SchemeC:
    """First-order potential (and centrifugal) corrections: level shift = <psi| dV + dq J(J+1) h^2/2muR^2 |psi>."""

    name = "C"
    NV = {"B": 30, "X": 16}
    NQ = {"B": 8, "X": 6}
    SV_GRID = (0.1, 0.3, 1.0, 3.0, 10.0, 30.0, 100.0, 300.0)   # MHz: scale of dV's coefficients' second differences
    SQ_GRID = (0.3, 3.0, 30.0, 300.0, 3000.0)                 # 1e-7: the same for dq
    PRIOR = 30.0                                # MHz (dV) / 1e-7 x 10 (dq): the weak ridge on the coefficients
    Q_UNIT = 1e-7

    def __init__(self, D: Data):
        from i2spec.constants import HBAR2_2U, MHZ_PER_CM
        self.D = D
        centres = 0.5 * (D.edges[1:] + D.edges[:-1])
        # the R range of the data: 1e-3 to 1 - 1e-3 of each data level's psi^2
        data_levels = {k for x in D.rows for k in x["levels"]}
        self.range = {}
        for s in ("B", "X"):
            lo, hi = np.inf, -np.inf
            for k in data_levels:
                if k[0] != s:
                    continue
                cm = np.cumsum(D.dens[k][0]) / np.sum(D.dens[k][0])
                lo, hi = min(lo, centres[np.searchsorted(cm, 1e-3)]), max(hi, centres[np.searchsorted(cm, 1 - 1e-3)])
            self.range[s] = (float(lo), float(hi))
        Bv = {s: bspline_basis(centres, *self.range[s], self.NV[s]) for s in ("B", "X")}
        Bq = {s: bspline_basis(centres, *self.range[s], self.NQ[s]) for s in ("B", "X")}
        kap = HBAR2_2U / D.mu * MHZ_PER_CM * self.Q_UNIT
        # parameter blocks: V_B, V_X, V_ext, q_B, q_X
        nb = [self.NV["B"], self.NV["X"], self.NV["B"], self.NQ["B"], self.NQ["X"]]
        off = np.r_[0, np.cumsum(nb)]
        self.blocks = dict(zip(("VB", "VX", "Vext", "qB", "qX"), zip(off[:-1], off[1:])))
        self.p = int(off[-1])
        G = np.zeros((len(D.levels), self.p))
        for i, (s, v, J) in enumerate(D.levels):
            m, m2 = (np.asarray(a, float) for a in D.dens[(s, v, J)])
            a, b = self.blocks["V" + s]
            G[i, a:b] = m @ Bv[s]
            if s == "B" and v >= V_EXT_B:
                a, b = self.blocks["Vext"]
                G[i, a:b] = m @ Bv[s]
            a, b = self.blocks["q" + s]
            G[i, a:b] = J * (J + 1) * kap * (m2 @ Bq[s])
        self.G = G
        self.X = np.asarray(D.S @ G)
        self.XQ = np.asarray(D.Q @ G)
        self.pen = {}
        for name, (a, b) in self.blocks.items():
            n = b - a
            D2 = np.diff(np.eye(n), 2, axis=0)
            self.pen[name] = (a, b, D2.T @ D2)

    def penalty(self, sv, sq):
        P = np.eye(self.p) / self.PRIOR ** 2
        for name, (a, b, DD) in self.pen.items():
            P[a:b, a:b] += DD / (sv if name[0] == "V" else sq) ** 2
        return P

    def _solve(self, mask, P, cov=False):
        X, r, w = self.X, self.D.r, 1 / self.D.sig_eff
        ww = w.copy()
        for _ in range(8):
            Xw = X[mask] * ww[mask, None]
            c = np.linalg.solve(Xw.T @ Xw + P, Xw.T @ (r[mask] * ww[mask]))
            z = np.abs((r - X @ c) * w)
            ww = w * np.minimum(1.0, 3.0 / np.maximum(z, 1e-9))
        if not cov:
            return c
        Xw = X[mask] * ww[mask, None]
        C = np.linalg.inv(Xw.T @ Xw + P)
        edf = float(np.trace(C @ (Xw.T @ Xw)))
        chi2 = float(np.sum(((r[mask] - X[mask] @ c) * ww[mask]) ** 2)) / max(mask.sum() - edf, 1.0)
        return c, C * max(chi2, 1.0), chi2

    CV = "set"

    def fit(self, train, k=6):
        D = self.D
        if self.CV == "set":        # inner folds of whole data sets
            sets = sorted(set(D.sets[train]), key=lambda s: (-np.sum(D.sets[train] == s), s))
            fold = {s: i % k for i, s in enumerate(sets)}
            f_of = np.array([fold.get(s, -1) for s in D.sets])
        else:                       # inner folds of whole levels: by v' (v'' for the v' = 0 bands)
            f_of = np.array([(x["v_upper"] if x["v_upper"] else 100 + x["v_lower"]) % k for x in D.rows])
        best = None
        for sv in self.SV_GRID:
            for sq in self.SQ_GRID:
                P = self.penalty(sv, sq)
                loss, cvr, cvs = 0.0, np.zeros(D.N), np.zeros(D.N)
                for f in range(k):
                    te = train & (f_of == f)
                    c, C, _ = self._solve(train & (f_of != f), P, cov=True)
                    cvr[te] = D.r[te] - self.X[te] @ c
                    cvs[te] = np.sqrt(np.einsum("ij,jk,ik->i", self.X[te], C, self.X[te]) + D.sig_eff[te] ** 2)
                    z = np.abs(cvr[te] / D.sig_eff[te])
                    loss += np.sum(np.where(z < 3, 0.5 * z * z, 3 * z - 4.5))
                if best is None or loss < best[0]:
                    best = (loss, sv, sq, cvr.copy(), cvs.copy())
        _, self.sv, self.sq, cvr, cvs = best
        # the discrepancy term: the white per-row sigma that makes the inner-CV residuals' median |z| = 0.674
        rr, ss = cvr[train], cvs[train]
        lo, hi = 0.0, 100.0
        for _ in range(60):
            tau = 0.5 * (lo + hi)
            lo, hi = (tau, hi) if np.median(np.abs(rr) / np.hypot(ss, tau)) > 0.6745 else (lo, tau)
        self.tau = 0.5 * (lo + hi)
        self.c, self.C, self.chi2 = self._solve(train, self.penalty(self.sv, self.sq), cov=True)
        return self

    def predict(self, idx, with_tau=True):
        X = self.X[idx]
        return X @ self.c, np.sqrt(np.einsum("ij,jk,ik->i", X, self.C, X) + self.D.sig_eff[idx] ** 2
                                   + (self.tau ** 2 if with_tau else 0.0))

    def predict_queries(self, with_tau=True):
        X = self.XQ
        return X @ self.c, np.sqrt(np.einsum("ij,jk,ik->i", X, self.C, X) + (self.tau ** 2 if with_tau else 0.0))

    def dV(self, s, R):
        a, b = self.blocks["V" + s]
        return bspline_basis(np.atleast_1d(R), *self.range[s], self.NV[s]) @ self.c[a:b]


# --------------------------------------------------------------------------------------------- splits

def splits(D: Data):
    """name -> list of (fold label, train mask, test index array)."""
    out = {}
    out["loso"] = [(s, D.sets != s, np.flatnonzero(D.sets == s)) for s in sorted(set(D.sets))]
    A = SchemeA(D).fit(np.ones(D.N, bool))
    lolo = []
    for key in sorted(A.lev, key=str):
        m = np.array([key in t for t in D.touch])
        if m.sum() >= 2:
            lolo.append((f"{key[0]}{key[1]}", ~m, np.flatnonzero(m)))
    out["lolo"] = lolo
    new = np.isin(D.sets, NEW)
    out["predtest"] = [("new sets", ~new, np.flatnonzero(new & (D.kind == "frequency")))]
    return out


def run_scheme(name):
    D = Data()
    S = splits(D)
    cls = dict(A=SchemeA, B=SchemeB, C=SchemeC)[name]
    model = cls(D)
    t0 = time.time()
    full = model.fit(np.ones(D.N, bool))
    qmu, qsd = full.predict_queries()
    res = dict(queries=dict(names=D.qnames, mu=qmu.tolist(), sd=qsd.tolist()))
    if name == "B":
        res["theta_full"] = dict(zip(SchemeB.NAMES, np.exp(full.th).tolist()))
        th_full = full.th.copy()
    if name == "C":
        res["full"] = dict(sv=full.sv, sq=full.sq, chi2=full.chi2, tau=full.tau, range=full.range,
                           queries_sd_cov=full.predict_queries(with_tau=False)[1].tolist(),
                           dV={s: dict(R=np.linspace(*full.range[s], 60).tolist(),
                                       dV=full.dV(s, np.linspace(*full.range[s], 60)).tolist()) for s in ("B", "X")})
    mu_in, _ = full.predict(np.arange(D.N))
    res["in_sample_resid"] = (D.r - mu_in).tolist()
    print(f"[{name}] full fit {time.time() - t0:.0f} s", flush=True)
    for split, folds in S.items():
        out = []
        for label, train, test in folds:
            m = model.fit(train, th_full) if name == "B" else model.fit(train)
            mu, sd = m.predict(test)
            extra = {}
            if name == "B":
                extra = dict(theta=np.exp(m.th).tolist())
            if name == "C":
                extra = dict(sv=m.sv, sq=m.sq, tau=m.tau, sd_cov=m.predict(test, with_tau=False)[1].tolist())
            out.append(dict(label=label, test=test.tolist(), mu=mu.tolist(), sd=sd.tolist(), **extra))
        res[split] = out
        print(f"[{name}] {split}: {len(folds)} folds, {time.time() - t0:.0f} s", flush=True)
    with open(OUT / f"bakeoff_{name}.pkl", "wb") as f:
        pickle.dump(res, f)


# --------------------------------------------------------------------------------------------- report

COLORS = {"A": "#2a78d6", "B": "#eb6834", "C": "#1baf7a"}
LABELS = {"A": "A per-level polynomials", "B": "B Gaussian process", "C": "C potential correction"}


def metrics(res, z):
    res, z = np.asarray(res), np.asarray(z)
    if not res.size:
        return None
    a = np.abs(z)
    return dict(n=int(res.size), rms=float(np.sqrt(np.mean(res ** 2))), median_abs=float(np.median(np.abs(res))),
                within1=float(np.mean(a <= 1)), within2=float(np.mean(a <= 2)), within3=float(np.mean(a <= 3)),
                rms_z=float(np.sqrt(np.mean(z ** 2))), median_abs_z=float(np.median(a)))


def report():
    D = Data()
    R = {s: pickle.load(open(OUT / f"bakeoff_{s}.pkl", "rb")) for s in "ABC"}
    out = dict(rows=dict(n=D.N, precise=int(D.precise.sum()), absolute=int(np.sum(D.kind == "frequency")),
                         region=f"B v' <= {V_MAX_B}, X v'' <= {V_MAX_X}", base=BASE),
               splits={}, per_fold={}, queries={}, scheme_B=R["B"].get("theta_full"), scheme_C=R["C"].get("full"),
               in_sample={})
    subsets = dict(all=lambda t: np.ones(t.size, bool), precise=lambda t: D.precise[t],
                   precise_absolute=lambda t: D.precise[t] & (D.kind[t] == "frequency"),
                   absolute=lambda t: D.kind[t] == "frequency", interval=lambda t: D.kind[t] == "interval")
    for s in "ABC":
        ri = np.array(R[s]["in_sample_resid"])
        out["in_sample"][s] = dict(rms_precise=float(np.sqrt(np.mean(ri[D.precise] ** 2))),
                                   median_precise=float(np.median(np.abs(ri[D.precise]))))
    for split in ("loso", "lolo", "predtest"):
        out["splits"][split] = {}
        out["per_fold"][split] = {}
        for s in "ABC":
            t = np.concatenate([f["test"] for f in R[s][split]]).astype(int)
            mu = np.concatenate([f["mu"] for f in R[s][split]])
            sd = np.concatenate([f["sd"] for f in R[s][split]])
            res = D.r[t] - mu
            z = res / sd
            out["splits"][split][s] = {k: metrics(res[f(t)], z[f(t)]) for k, f in subsets.items()}
            pf = {}
            for f in R[s][split]:
                tt = np.array(f["test"], int)
                rr = D.r[tt] - np.array(f["mu"])
                pm = D.precise[tt]
                pf[f["label"]] = dict(n=int(tt.size), n_precise=int(pm.sum()),
                                      rms_precise=float(np.sqrt(np.mean(rr[pm] ** 2))) if pm.any() else None,
                                      median_precise=float(np.median(np.abs(rr[pm]))) if pm.any() else None,
                                      rms_z_precise=float(np.sqrt(np.mean((rr[pm] / np.array(f["sd"])[pm]) ** 2))) if pm.any() else None,
                                      bare_rms_precise=float(np.sqrt(np.mean(D.r[tt][pm] ** 2))) if pm.any() else None)
            out["per_fold"][split][s] = pf
        t = np.concatenate([f["test"] for f in R["C"][split]]).astype(int)
        res = D.r[t] - np.concatenate([f["mu"] for f in R["C"][split]])
        z = res / np.concatenate([f["sd_cov"] for f in R["C"][split]])
        out["splits"][split]["C_cov"] = {k: metrics(res[f(t)], z[f(t)]) for k, f in subsets.items()}
        out["splits"][split]["C_tau"] = [f["tau"] for f in R["C"][split]]
        # bare: no correction at all, for reference
        t = np.concatenate([f["test"] for f in R["A"][split]]).astype(int)
        out["splits"][split]["bare"] = {k: metrics(D.r[t][f(t)], np.zeros(int(f(t).sum()))) for k, f in subsets.items()}
    # leave-one-level-out by region: the visible B levels, the MLR-built B levels, the NIR levels
    groups = {"B v'=3-43": lambda l: l[0] == "B" and 3 <= int(l[1:]) <= 43,
              "B v'=44-50": lambda l: l[0] == "B" and int(l[1:]) >= 44,
              "NIR: B v'=1-2, X v''=11-17": lambda l: l in ("B1", "B2") or l[0] == "X"}
    out["lolo_by_region"] = {}
    for g, fn in groups.items():
        out["lolo_by_region"][g] = {}
        for s in "ABC":
            res, z, bare = [], [], []
            for f in R[s]["lolo"]:
                if fn(f["label"]):
                    t = np.array(f["test"], int); p = D.precise[t]
                    rr = D.r[t] - np.array(f["mu"])
                    res += list(rr[p]); z += list((rr / np.array(f["sd"]))[p]); bare += list(D.r[t][p])
            out["lolo_by_region"][g][s] = metrics(res, z)
        out["lolo_by_region"][g]["bare"] = metrics(bare, np.zeros(len(bare)))
    for s in "ABC":
        q = R[s]["queries"]
        out["queries"][s] = {n: dict(mu=m, sd=d) for n, m, d in zip(q["names"], q["mu"], q["sd"])}
    (OUT / "bakeoff.json").write_text(json.dumps(out, indent=1))
    print_tables(out)
    figure(out, D)


def print_tables(out):
    for split in ("loso", "lolo", "predtest"):
        for sub in ("precise", "all"):
            print(f"\n{split} [{sub}]  scheme   n    rms   median |r|  <=1s  <=2s  <=3s  rms z")
            for s in ("bare", "A", "B", "C", "C_cov"):
                m = out["splits"][split][s][sub]
                if m:
                    print(f"{'':22s}{s:6s}{m['n']:4d} {m['rms']:7.2f} {m['median_abs']:7.2f}   "
                          + ("" if s == "bare" else f"{m['within1']:.2f}  {m['within2']:.2f}  {m['within3']:.2f}  {m['rms_z']:5.2f}"))
    for g, d in out["lolo_by_region"].items():
        print(f"\nlolo {g} [precise]: " + "; ".join(f"{s} rms {m['rms']:.2f} med {m['median_abs']:.2f}" + (
            "" if s == "bare" else f" 1/2/3s {m['within1']:.2f} {m['within2']:.2f} {m['within3']:.2f} rms z {m['rms_z']:.2f}")
            for s, m in d.items()) + f" (n {d['A']['n']})")
    print("\nunmeasured levels (MHz): scheme mean +- sd at J = 20, 60, 100, 150")
    for v in SAMPLE_B:
        for s in "ABC":
            q = out["queries"][s]
            print(f"B v'={v:2d} R(J) v'-0  {s}: " + "  ".join(f"{q[f'B{v} R({J})']['mu']:+7.2f}+-{q[f'B{v} R({J})']['sd']:5.2f}" for J in (20, 60, 100, 150)))
    for v in SAMPLE_X:
        for s in "ABC":
            q = out["queries"][s]
            print(f"X v''={v} - X 0     {s}: " + "  ".join(f"{q[f'X{v} J={J}']['mu']:+7.2f}+-{q[f'X{v} J={J}']['sd']:5.2f}" for J in (20, 60, 100, 150)))


def figure(out, D):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                         "axes.edgecolor": "#52514e", "axes.labelcolor": "#0b0b0b", "xtick.color": "#52514e",
                         "ytick.color": "#52514e", "axes.grid": True, "grid.color": "#e6e5e0", "grid.linewidth": 0.6})
    fig = plt.figure(figsize=(11, 8.2), facecolor="#fcfcfb")
    gs = fig.add_gridspec(2, 3, height_ratios=[1, 1.05], hspace=0.42, wspace=0.28)
    # (a) leave-one-level-out rms per level
    ax = fig.add_subplot(gs[0, :2])
    pf = out["per_fold"]["lolo"]
    labs = [l for l in pf["A"] if pf["A"][l]["n_precise"]]
    labs.sort(key=lambda l: (l[0] != "B", int(l[1:])))
    x = np.arange(len(labs))
    for k, s in enumerate("ABC"):
        yv = [max(pf[s][l]["rms_precise"], 1e-3) for l in labs]
        ax.scatter(x + (k - 1) * 0.22, yv, s=16, color=COLORS[s], label=LABELS[s], zorder=3, edgecolor="#fcfcfb", linewidth=0.5)
    ax.scatter(x, [max(pf["A"][l]["bare_rms_precise"], 1e-3) for l in labs], marker="_", s=60, color="#52514e",
               label="bare potentials", zorder=2)
    ax.set_yscale("log")
    ax.set_xticks(x, [l.replace("B", "B ").replace("X", "X ") for l in labs], rotation=90, fontsize=7)
    ax.set_ylabel("held-out rms, precise rows (MHz)")
    ax.set_title("(a) leave one level out: every row touching the level predicted from the rest", loc="left", fontsize=9)
    fig.legend(*ax.get_legend_handles_labels(), frameon=False, fontsize=8.5, ncol=4, loc="lower center",
               bbox_to_anchor=(0.5, 0.92))
    # (b) calibration: fraction within k sigma
    ax = fig.add_subplot(gs[0, 2])
    ks = np.array([1, 2, 3])
    ideal = [0.683, 0.954, 0.997]
    for split, ls in (("loso", "-"), ("lolo", "--"), ("predtest", ":")):
        for s in "ABC":
            m = out["splits"][split][s]["precise"] if split != "predtest" else out["splits"][split][s]["all"]
            ax.plot(ks, [m["within1"], m["within2"], m["within3"]], ls, color=COLORS[s], lw=1.6, marker="o", ms=4)
    ax.plot(ks, ideal, color="#0b0b0b", lw=0.8, marker="x", ms=4)
    ax.text(1.05, 0.66, "Gaussian", fontsize=7, va="top", color="#0b0b0b")
    ax.set_xticks(ks, ["1σ", "2σ", "3σ"])
    ax.set_xlim(0.8, 3.6)
    ax.set_ylabel("fraction of held-out rows within")
    ax.set_title("(b) calibration\n— by set, -- by level, ··· new sets", loc="left", fontsize=9)
    # (c) unmeasured levels
    Js = np.array(SAMPLE_J)
    for k, v in enumerate((12, 27, 38)):
        ax = fig.add_subplot(gs[1, k])
        for s in "ABC":
            q = out["queries"][s]
            mu = np.array([q[f"B{v} R({J})"]["mu"] for J in Js])
            sd = np.array([q[f"B{v} R({J})"]["sd"] for J in Js])
            ax.plot(Js, mu, color=COLORS[s], lw=2, label=LABELS[s])
            ax.fill_between(Js, mu - sd, mu + sd, color=COLORS[s], alpha=0.12, lw=0)
        ax.axhline(0, color="#52514e", lw=0.6)
        ax.set_xlabel("J″ of R(J″) " + f"{v}-0")
        if k == 0:
            ax.set_ylabel("predicted correction (MHz), ±1σ")
        ax.set_title(f"({'cde'[k]}) B v′ = {v}, never measured", loc="left", fontsize=9)
    fig.savefig(ROOT / "docs/figures/bakeoff.png", dpi=150, bbox_inches="tight", facecolor="#fcfcfb")
    print("wrote docs/figures/bakeoff.png")


if __name__ == "__main__":
    args = sys.argv[1:] or ["rows", "run", "figure"]
    if "rows" in args or not CACHE.exists():
        build_rows()
    for s in "ABC":
        if "run" in args or f"run{s}" in args:
            run_scheme(s)
    if "figure" in args or "report" in args:
        report()
