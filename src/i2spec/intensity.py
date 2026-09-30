"""Line intensities of the I₂ B-X system. Inputs and sources: docs/research/intensity-inputs.md.

Integrated absorption cross section of one rovibronic line (Tellinghuisen, JCP 134, 084301 (2011), eq. 1):

    S = 2π² ν / (3 ε₀ h c) · s_J'J'' / (2J''+1) · |<v'J'|μ_e(R)|v''J''>|² · f(v'', J''),

with s_J'J'' = J'' for P and J''+1 for R lines, and f the fractional population of the lower level.
In practical units S [cm] = 4.1625e-19 · ν[cm⁻¹] · s/(2J''+1) · |μ|²[D²] · f. The matrix elements
use the J-dependent wavefunctions of both states, so Herman-Wallis effects are included.

Because f = g_ns (2J''+1) exp(-c2 E''/T) / Q(T), a line list is stored in the temperature-independent
form S₀ = S_UNIT · ν · s · |μ|² · g_ns, with S(T) = S₀ exp(-c2 E''/T) / Q(T) (MasterLineList).
"""

from __future__ import annotations

import hashlib
import os
from dataclasses import dataclass, fields
from math import pi
from pathlib import Path

import warnings

import numpy as np

from .constants import MHZ_PER_CM, ISOTOPOLOGUES, NUCLEAR_SPIN
from .model import RovibronicModel

_H, _C, _EPS0 = 6.62607015e-34, 299792458.0, 8.8541878188e-12
DEBYE = 1e-21 / _C  # C m
#: S [cm] per ν [cm⁻¹] per μ² [D²].
S_UNIT = 2 * pi**2 * 100 * DEBYE**2 / (3 * _EPS0 * _H * _C) * 100
#: Second radiation constant hc/k, cm K.
C2 = 1.438776877

#: Grid of the B state, whose points X shares, so that transition matrix elements are plain sums over it.
#: B cannot start lower: its nonadiabatic function alpha(R) is a polynomial fitted where the levels
#: are, and it diverges below ~2.35 Å (+19 at 2.33 Å, +1171 at 2.30, -19 000 at 2.25). On a grid
#: reaching 2.20 Å the B levels collapse by 10²² cm⁻¹ at J = 150.
SHARED_GRID = dict(rmin=2.33, rmax=7.0, step=0.005)
#: X is solved on the same points extended inward to about X_RMIN. Stopping at 2.33 Å put a hard wall
#: against the inner turning points of high v'' (2.39 Å at v'' = 25, 2.35 Å at v'' = 39) and raised
#: those levels by 51 MHz at v'' = 25, 1.6 GHz at 30 and up to 86 GHz. B wavefunctions vanish below
#: 2.33 Å, so the overlaps only need the shared points.
X_RMIN = 2.10
#: Line positions and lower-level energies come from the B-spline solver, which converges to under
#: 1 kHz, not from the sinc-DVR, which is off by 0.24 MHz at R(56) 32-0 and by up to ~3 MHz elsewhere.
#: The DVR eigenvectors are still used for the intensities. The boxes match the DVR's, so that
#: both solvers see the same box states above the dissociation asymptote.
POSITION_GRIDS = {"X": dict(rmin=2.10, rmax=7.0, h=0.01, order=10, nlev=70),
                  "B": dict(rmin=2.35, rmax=7.0, h=0.01, order=10, nlev=110)}
#: A DVR level takes the energy of the B-spline level nearest to it, if that is within this (cm⁻¹).
MATCH_TOLERANCE = 0.01
#: An eigenvector with more than EDGE_WEIGHT of its probability within EDGE (Å) of the outer grid
#: boundary is a box state (continuum above the dissociation limit or the centrifugal barrier).
EDGE, EDGE_WEIGHT = 0.5, 1e-6
_CACHE_VERSION = 7   # 4: level corrections of the parameter set; 5: extended-range potentials (i2spec2026d);
                     # 6: the NIR band correction (i2spec2026e); 7: DVR levels beyond the position model's
                     # physical ones are continuum on every grid


def mu_tellinghuisen2011(R):
    """B-X transition moment |μ_e(R)| in debye: μ_e² = c1 exp[-c2 (R - c3)²]/R², Tellinghuisen 2011 eq. (10)."""
    return np.sqrt(21.5 * np.exp(-1.18 * (R - 3.65) ** 2)) / R


def nuclear_spin_weight(J, isotopologue):
    """Nuclear-spin statistical weight of the X(0g+) level J."""
    a, b = ISOTOPOLOGUES[isotopologue]
    if a != b:
        return round((2 * NUCLEAR_SPIN[a] + 1) * (2 * NUCLEAR_SPIN[b] + 1))
    return sum(2 * I + 1 for I in range(round(2 * NUCLEAR_SPIN[a]) + 1) if (I - J) % 2 == 0)


def intensity_model(isotopologue="127I2", nlev_x=40, nlev_b=70, grid=None, x_rmin=X_RMIN):
    """A sinc-DVR RovibronicModel for intensities: B on ``grid``, X on the same points extended inward.

    Its eigenvectors give the transition moments. Its energies are only good to a few MHz, so
    master_line_list takes positions from the B-spline solver instead.
    """
    g = {**SHARED_GRID, **(grid or {})}
    n_inner = max(int(round((g["rmin"] - x_rmin) / g["step"])), 0)
    gx = dict(g, rmin=g["rmin"] - n_inner * g["step"], nlev=nlev_x)
    return RovibronicModel(isotopologue, grids={"X": gx, "B": dict(g, nlev=nlev_b)}, solver="dvr")


#: For an intensity model whose B box reaches the dissociation region (intensity_model(grid=NEAR_DISSOCIATION)):
#: the position grid follows it out, with the levels a 12 A box holds (docs/design/near-dissociation.md).
NEAR_DISSOCIATION = dict(rmax=12.0)


def position_model(model) -> RovibronicModel:
    """The B-spline model whose level energies master_line_list uses for ``model``'s lines."""
    rmax = model.grids["B"]["rmax"]
    if rmax <= POSITION_GRIDS["B"]["rmax"]:
        return RovibronicModel(model.isotopologue, model.parameters, grids=POSITION_GRIDS)
    far = {"X": POSITION_GRIDS["X"], "B": dict(POSITION_GRIDS["B"], rmax=rmax, nlev=150)}
    return RovibronicModel(model.isotopologue, model.parameters, grids=far)


def match_levels(dvr, exact, tol=MATCH_TOLERANCE, n_index=0):
    """Each DVR energy replaced by the matching one of ``exact`` (sorted), and a mask of those matched.

    The first ``n_index`` levels match by index: the same v in both solvers, however far apart the
    energies are, which is what an extended-range potential set needs (its X levels at v'' >= 48 sit
    cm⁻¹ from the published curve's, whose eigenvectors give the intensities). ``n_index`` is the length
    of the physical prefix of ``exact`` (``physical_prefix``): beyond it, box states interleave with the
    quasi-bound levels behind the centrifugal barrier and the two solvers need not order them alike,
    so there the nearest energy within ``tol`` is taken. Unmatched levels keep their DVR energy.
    """
    dvr = np.asarray(dvr, dtype=float)
    idx = np.searchsorted(exact, dvr).clip(1, len(exact) - 1)
    lo, hi = exact[idx - 1], exact[idx]
    nearest = np.where(np.abs(dvr - lo) <= np.abs(hi - dvr), lo, hi)
    ok = np.abs(nearest - dvr) < tol
    n = min(n_index, len(dvr), len(exact))
    if n:
        nearest[:n] = exact[:n]
        ok[:n] = True
    return np.where(ok, nearest, dvr), ok


def physical_prefix(solver, J, edge=EDGE, weight=EDGE_WEIGHT):
    """How many of a B-spline solver's levels at J, counted from v = 0, are physical (not box states):
    the same criterion as ``bound``, on the quadrature grid with its weights."""
    _, psi = solver.wavefunctions(J)
    near_edge = (solver.W[:, None] * psi ** 2)[solver.R > solver.R[-1] - edge].sum(axis=0)
    ok = near_edge < weight
    return int(ok.size if ok.all() else np.argmin(ok))


def shared_offset(X, B):
    """Index of B's first grid point in X's grid; X must extend B's points inward."""
    k = int(round((B.R[0] - X.R[0]) / (B.R[1] - B.R[0])))
    if k < 0 or X.R.size != k + B.R.size or not np.allclose(X.R[k:], B.R, atol=1e-9):
        raise ValueError("X must be on B's grid points, extended inward; build the model with intensity_model()")
    return k


def partition_function(model, T, rel_tol=1e-10):
    """Rovibrational partition function of X, with nuclear-spin weights, energies relative to X(0, 0)."""
    e00 = model.energy("X", 0, 0)
    Q, J = 0.0, 0
    while True:
        e = model.states["X"].levels(J) - e00
        term = nuclear_spin_weight(J, model.isotopologue) * (2 * J + 1) * np.exp(-C2 * e / T).sum()
        Q += term
        if J > 10 and term < rel_tol * Q:
            return Q
        J += 1


def bound(vectors, R, edge=EDGE, weight=EDGE_WEIGHT):
    """Mask of eigenvectors (columns) that are bound, i.e. carry little probability near the outer grid edge."""
    return (vectors[R > R[-1] - edge] ** 2).sum(axis=0) < weight


def bound_prefix(energies, vectors, R, e00):
    """Mask of the levels below the first box state, and the energy (above X(0, 0)) half a level
    spacing above the last of them, where the continuum takes over (nan if fewer than two levels)."""
    ok = bound(vectors, R)
    n = ok.size if ok.all() else int(np.argmin(ok))
    cut = energies[n - 1] + (energies[n - 1] - energies[n - 2]) / 2 - e00 if n >= 2 else np.nan
    return np.arange(ok.size) < n, cut


@dataclass
class LineList:
    nu: np.ndarray  # cm⁻¹
    S: np.ndarray  # integrated cross section per molecule, cm
    v_upper: np.ndarray
    v_lower: np.ndarray
    J_lower: np.ndarray
    branch: np.ndarray  # +1 for R, -1 for P
    E_lower: np.ndarray  # cm⁻¹ above X(0, 0)
    T: float

    def __len__(self):
        return len(self.nu)

    def label(self, k):
        b = "R" if self.branch[k] > 0 else "P"
        return f"{b}({self.J_lower[k]}) {self.v_upper[k]}-{self.v_lower[k]}"


@dataclass
class MasterLineList:
    """Temperature-independent line list; S(T) = strength0 · exp(-c2 E''/T) / Q(T)."""

    nu: np.ndarray
    strength0: np.ndarray  # S_UNIT · ν · s_J'J'' · |<v'J'|μ|v''J''>|² · g_ns(J''), cm
    v_upper: np.ndarray
    v_lower: np.ndarray
    J_lower: np.ndarray
    branch: np.ndarray
    E_lower: np.ndarray
    x_levels: np.ndarray  # X term values above X(0, 0), shape (J, v), for the partition function
    upper_cut: np.ndarray  # per J': B term value (above X(0, 0)) above which the list has no levels
    isotopologue: str

    def __len__(self):
        return len(self.nu)

    def partition_function(self, T):
        J = np.arange(self.x_levels.shape[0])
        weight = np.array([nuclear_spin_weight(j, self.isotopologue) for j in J]) * (2 * J + 1)
        return float((weight[:, None] * np.exp(-C2 * self.x_levels / T)).sum())

    def at(self, T, nu_min=-np.inf, nu_max=np.inf, S_min=0.0) -> LineList:
        S = self.strength0 * np.exp(-C2 * self.E_lower / T) / self.partition_function(T)
        keep = (self.nu >= nu_min) & (self.nu <= nu_max) & (S >= S_min)
        return LineList(self.nu[keep], S[keep], self.v_upper[keep], self.v_lower[keep], self.J_lower[keep],
                        self.branch[keep], self.E_lower[keep], T=T)

    def save(self, path):
        np.savez(path, **{f.name: getattr(self, f.name) for f in fields(self)})

    @classmethod
    def load(cls, path):
        z = np.load(path)
        return cls(**{f.name: z[f.name] for f in fields(cls) if f.name != "isotopologue"}, isotopologue=str(z["isotopologue"]))


def cache_dir() -> Path:
    """Directory for cached line lists: $I2SPEC_CACHE or ~/.cache/i2spec."""
    return Path(os.environ.get("I2SPEC_CACHE", Path.home() / ".cache" / "i2spec"))


def _cache_key(model, *args):
    text = repr((_CACHE_VERSION, model.parameters, model.isotopologue, model.solver, sorted(model.grids["X"].items()),
                 sorted(model.grids["B"].items()), sorted((k, sorted(v.items())) for k, v in POSITION_GRIDS.items()),
                 MATCH_TOLERANCE, args))
    return hashlib.sha1(text.encode()).hexdigest()[:16]


class _Build:
    """What one process needs to turn a range of J into lines: the DVR model, its position model, and the
    constants of the build. Built once per process (the serial path uses one in-process)."""

    def __init__(self, isotopologue, parameters, grids, solver, nu_min, nu_max, temps, S_min, mu):
        self.model = RovibronicModel(isotopologue, parameters, grids=grids, solver=solver)
        self.X, self.B = self.model.states["X"], self.model.states["B"]
        self.k0 = shared_offset(self.X, self.B)
        self.iso = isotopologue
        self.exact = position_model(self.model)
        self.e00 = self.exact.energy("X", 0, 0)
        self.top = {st: (self.exact.extended[st] if self.exact.extended is not None else self.exact.states[st])
                    for st in ("X", "B")}
        self.nu_min, self.nu_max, self.temps, self.S_min = nu_min, nu_max, np.asarray(temps), S_min
        self.mu_R = mu(self.B.R)

    def x_levels(self, J):
        """X term values of the physical levels at J (above X(0,0)), and which matched a B-spline level."""
        e_matched, ok = match_levels(self.X.levels(J), self.exact.levels("X", J),
                                     n_index=physical_prefix(self.top["X"], J))
        return e_matched - self.e00, ok

    def lines(self, Js, levels, x_matched, Q):
        """Every line from the X levels at Js; B levels are solved here (J' = J +- 1) and reused within Js."""
        cols = {k: [] for k in ("nu", "strength0", "v_upper", "v_lower", "J_lower", "branch", "E_lower")}
        upper, cuts, unmatched = {}, {}, 0
        for J in Js:
            _, c_x = self.X.states(J)
            E = levels[J]
            e_x = E + self.e00
            boltz = (np.exp(-C2 * E[None, :] / self.temps[:, None]) / Q[:, None]).max(axis=0)
            weighted = self.mu_R[:, None] * c_x[self.k0:]
            ok_x = bound(c_x, self.X.R)
            g = nuclear_spin_weight(J, self.iso)
            for branch, J_up, s in ((+1, J + 1, J + 1), (-1, J - 1, J)):
                if J_up < 0:
                    continue
                if J_up not in upper:
                    e_b, c_b = self.B.states(J_up)
                    ok, _ = bound_prefix(e_b, c_b, self.B.R, self.e00)
                    n_phys = physical_prefix(self.top["B"], J_up)
                    # a DVR level beyond the position model's last physical level has no position: treat it as
                    # continuum. The published curve holds a few more levels in the box than the extended MLR
                    # does, near the limit on the 12 A grid and, since mlr_b_2026d (wider at long range), at
                    # v' ~ 63-66 on the default 7 A one, where the MLR's levels already touch the wall
                    ok[n_phys:] = False
                    matched, good = match_levels(e_b[ok], self.exact.levels("B", J_up), n_index=n_phys)
                    e_b = np.concatenate([matched, e_b[~ok]])
                    _, cuts[J_up] = bound_prefix(e_b, c_b, self.B.R, self.e00)
                    upper[J_up] = (e_b, c_b, ok, np.concatenate([good, np.ones((~ok).sum(), bool)]))
                e_b, c_b, ok_b, b_matched = upper[J_up]
                nu = e_b[:, None] - e_x[None, :]
                corr = self.exact.corrections
                if corr is not None and corr.band is not None and len(e_b):
                    label = "R" if branch > 0 else "P"
                    for vx in range(corr.band["v_lower"][0], min(corr.band["v_lower"][1] + 1, nu.shape[1])):
                        nu[0, vx] += corr.band_shift(0, vx, J, label) / MHZ_PER_CM
                strength0 = S_UNIT * nu * s * (c_b.T @ weighted) ** 2 * g
                keep = ((nu >= self.nu_min) & (nu <= self.nu_max) & (strength0 * boltz[None, :] >= self.S_min)
                        & ok_b[:, None] & ok_x[None, :])
                # a line whose DVR level has no counterpart in the position model (a quasi-bound level the two
                # solvers order differently, at the B asymptote) is dropped rather than carried at DVR accuracy
                matched_pair = b_matched[:, None] & x_matched[J][None, :]
                unmatched += int((keep & ~matched_pair).sum())
                keep &= matched_pair
                vb, vx = np.nonzero(keep)
                cols["nu"].append(nu[vb, vx])
                cols["strength0"].append(strength0[vb, vx])
                cols["v_upper"].append(vb)
                cols["v_lower"].append(vx)
                cols["J_lower"].append(np.full(vb.size, J))
                cols["branch"].append(np.full(vb.size, branch))
                cols["E_lower"].append(E[vx])
            upper.pop(J - 1, None)
        return {k: (np.concatenate(v) if v else np.array([])) for k, v in cols.items()}, cuts, unmatched


_WORKER: "_Build | None" = None


_THREAD_VARS = ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "MKL_NUM_THREADS")


def _worker_init(args):
    global _WORKER
    _WORKER = _Build(*args)


def _worker_x(Js):
    return [(J, *_WORKER.x_levels(J)) for J in Js]


def _worker_lines(job):
    Js, levels, x_matched, Q = job
    return _WORKER.lines(Js, levels, x_matched, Q)


#: Worker processes for master_line_list; I2SPEC_WORKERS overrides, 1 runs in-process.
def _n_workers():
    return max(1, int(os.environ.get("I2SPEC_WORKERS", min(8, os.cpu_count() or 1))))


def master_line_list(model, nu_min=11000.0, nu_max=20100.0, T_range=(200.0, 600.0), S_min=1e-28,
                     mu=mu_tellinghuisen2011, cache=True, workers=None) -> MasterLineList:
    """All bound-bound B-X lines of ``model`` with nu_min <= ν <= nu_max and S >= S_min at some T in T_range.

    ``model`` must come from intensity_model(): its eigenvectors give the strengths. Positions, lower-level
    energies and the partition function use the B-spline levels of position_model(model), matched to the
    DVR levels by energy. Upper levels stop at the first box state of each J' (see bound_prefix);
    absorption above that is continuum.Continuum's. The result is cached in cache_dir() unless
    cache=False.

    The work is per J and runs in ``workers`` processes (default: I2SPEC_WORKERS, else up to 8): first the
    X levels and the partition function, then the lines, in blocks of consecutive J so each B level is
    solved once per block. workers=1 is the same computation in-process.
    """
    key = _cache_key(model, float(nu_min), float(nu_max), tuple(map(float, T_range)), float(S_min),
                     f"{mu.__module__}.{mu.__qualname__}")
    path = cache_dir() / f"linelist_{model.isotopologue}_{key}.npz"
    if cache and path.exists():
        return MasterLineList.load(path)

    temps = np.unique([float(T_range[0]), float(np.sqrt(T_range[0] * T_range[1])), float(T_range[1])])
    args = (model.isotopologue, model.parameters, model.grids, model.solver, nu_min, nu_max, temps, S_min, mu)
    n = _n_workers() if workers is None else max(1, int(workers))
    pool = None
    main = __import__("__main__")
    if n > 1 and not getattr(main, "__file__", None):
        n = 1          # spawned workers re-import __main__; from stdin or a REPL there is none, so stay in-process
    if n > 1:
        import multiprocessing as mp
        saved = {k: os.environ.get(k) for k in _THREAD_VARS}
        os.environ.update({k: "1" for k in _THREAD_VARS})     # each worker one BLAS thread; set before spawning
        try:
            pool = mp.get_context("spawn").Pool(n, initializer=_worker_init, initargs=(args,))
        finally:
            for k, v in saved.items():
                if v is None:
                    os.environ.pop(k, None)
                else:
                    os.environ[k] = v
        run_x = lambda blocks: [r for part in pool.map(_worker_x, blocks) for r in part]      # noqa: E731
        run_lines = lambda jobs: pool.map(_worker_lines, jobs)                                # noqa: E731
    else:
        _worker_init(args)
        run_x = lambda blocks: [r for b in blocks for r in _worker_x(b)]                     # noqa: E731
        run_lines = lambda jobs: [_worker_lines(j) for j in jobs]                             # noqa: E731
    try:
        # X term values and partition functions, until a J contributes < 1e-9 of Q at the highest temperature
        levels, x_matched, Q, J, done = [], [], np.zeros_like(temps), 0, False
        block = 8
        while not done:
            for Jr, e, ok in run_x([list(range(J + k * block, J + (k + 1) * block)) for k in range(max(n, 1))]):
                if done:
                    break
                term = nuclear_spin_weight(Jr, model.isotopologue) * (2 * Jr + 1) * np.exp(-C2 * e[None, :] / temps[:, None]).sum(axis=1)
                levels.append(e)
                x_matched.append(ok)
                Q += term
                if Jr > 10 and term[-1] < 1e-9 * Q[-1]:
                    done = True
            J += max(n, 1) * block
        nJ = len(levels)
        size = max(4, -(-nJ // (4 * n)))
        blocks = [list(range(a, min(a + size, nJ))) for a in range(0, nJ, size)]
        jobs = [(b, {j: levels[j] for j in b}, {j: x_matched[j] for j in b}, Q) for b in blocks]
        results = run_lines(jobs)
    finally:
        if pool is not None:
            pool.close()
            pool.join()

    upper_cut = np.full(nJ + 1, np.nan)
    cols = {k: [] for k in ("nu", "strength0", "v_upper", "v_lower", "J_lower", "branch", "E_lower")}
    unmatched = 0
    for c, cuts, u in results:
        for k in cols:
            cols[k].append(c[k])
        for Ju, cut in cuts.items():
            if Ju <= nJ:
                upper_cut[Ju] = cut
        unmatched += u
    if unmatched > 100:
        raise RuntimeError(f"{unmatched} lines use a DVR level with no B-spline level within {MATCH_TOLERANCE} cm⁻¹; "
                           "widen POSITION_GRIDS['nlev'] or check the grids")
    if unmatched:
        warnings.warn(f"{unmatched} lines at the B asymptote dropped: no position-model level within {MATCH_TOLERANCE} cm⁻¹")
    arrays = {k: np.concatenate(v) for k, v in cols.items()}
    arrays["v_upper"] = arrays["v_upper"].astype(int); arrays["v_lower"] = arrays["v_lower"].astype(int)
    arrays["J_lower"] = arrays["J_lower"].astype(int); arrays["branch"] = arrays["branch"].astype(int)
    order = np.argsort(arrays["nu"], kind="stable")
    master = MasterLineList(**{k: v[order] for k, v in arrays.items()}, x_levels=np.array(levels),
                            upper_cut=upper_cut, isotopologue=model.isotopologue)
    if cache:
        path.parent.mkdir(parents=True, exist_ok=True)
        master.save(path)
    return master


def line_list(model, T, nu_min, nu_max, S_min=1e-26, mu=mu_tellinghuisen2011) -> LineList:
    """Lines with nu_min <= ν <= nu_max and S >= S_min at temperature T (not cached; see master_line_list)."""
    return master_line_list(model, nu_min, nu_max, (T, T), S_min, mu, cache=False).at(T, nu_min, nu_max, S_min)
