"""Continuous absorption of I₂: A←X, C←X, and B←X above the discrete line list.

Model of Tellinghuisen, J. Chem. Phys. 135, 054301 (2011) (docs/research/continuum-model.md), eq. (2):

    σ(ν) = S_UNIT · ν · Σ_{v'',J''} f(v'', J'') · |<E'|μ_e(R)|v'' J''>|²,    E' = E'' + ν,

with |E'> normalized per unit energy (cm⁻¹), J' = J'' (a Q branch stands in for P and R) and f the
Boltzmann population of X. Nuclear-spin weights average out.

* A←X and C←X: the paper's Table I potentials and transition moments, with G_ab = 1 and the
  tabulated μ. A is a pseudocontinuum: its potential is flat (U = T_e) beyond R_e, which folds the
  weak A←X bands into the continuum. Beyond the fitted exponential wall (R > 2.806 Å) the A curve
  continues as a Morse inner branch matched in value and slope, standing in for the RKR curve the
  paper used. The C curve continues as an R⁻ⁿ tail below 200 cm⁻¹, which only matters beyond 800 nm.
* B←X: the Hannover B potential with the paper's exponential inner wall below their crossing
  (2.715 Å; the Hannover extension is far too soft there), and the line-list μ_B(R). Absorption is
  counted only above the last B level of the line list for each J (MasterLineList.upper_cut), or
  above the dissociation limit if no cut is given. Continuum functions are normalized with the WKB
  amplitude inside the B well, so below the dissociation limit this is the level-averaged
  absorption of the bands missing from the line list, and above it the bound-free continuum.
"""

from __future__ import annotations

import hashlib
from dataclasses import dataclass
from math import log

import numpy as np
from numpy.polynomial import polynomial as P
from scipy.optimize import brentq, minimize_scalar

from .constants import HBAR2_2U, ISOTOPOLOGUES, REFERENCE_ISOTOPOLOGUE, reduced_mass, DEFAULT_PARAMETERS
from .intensity import C2, S_UNIT, cache_dir, mu_tellinghuisen2011
from .potentials import load_potentials
from .solver import RadialSolver

#: Tellinghuisen 2011, Table I: U = U_ref + A0 + B0 exp(Σ a_i z^i), z = R - R0 (cm⁻¹, Å);
#: μ = Σ μ_i (R - 2.7 Å)^i in debye. U_ref is T_e(A), T_e(B), or the X asymptote for C.
T11C = {
    "A": dict(R0=2.7, A0=-226.70, B0=3401.305, a=(-6.36, -4.0), mu=(0.2845, -0.048)),
    "B": dict(R0=2.7, A0=-1783.226, B0=4408.223, a=(-4.63,)),
    "C": dict(R0=2.83, A0=-211.0, B0=4599.3, a=(-3.297, 0.505, -0.46, 0.0, -0.33), mu=(0.4714, -0.16)),
}
A_WELL = 1639.9  # D_e(A), cm⁻¹
A_JOIN = 2.805916  # Å, end of the fitted A wall (paper eq. 5)
C_TAIL = 200.0  # cm⁻¹ above the asymptote, where the C curve switches to an R⁻ⁿ tail
RMIN, RMAX, STEP = 2.15, 4.5, 0.004  # X grid (Å); the Numerov step is STEP/4
WINDOWS = {"A": (4.5, 5.0), "B": (3.3, 3.8), "C": (4.5, 5.0)}  # Å, where the WKB amplitude is fitted
E_STEPS = {"A": 20.0, "B": 10.0, "C": 20.0}  # cm⁻¹, upper-state energy grids
POP_MIN = 1e-6  # X levels with exp(-c2 E''/T_max) below this are left out
GRID_STEP = 2.0  # cm⁻¹, grid on which σ is evaluated before interpolation
_CACHE_VERSION = 1


def _wall(p, R):
    """A0 + B0 exp(a1 z + a2 z² + ...), z = R - R0, and its derivative."""
    z = np.asarray(R, dtype=float) - p["R0"]
    coef = (0.0, *p["a"])
    e = p["B0"] * np.exp(P.polyval(z, coef))
    return p["A0"] + e, e * P.polyval(z, P.polyder(coef))


def a_potential(De_x):
    """A state (pseudocontinuum), cm⁻¹ above the X minimum; De_x is the X dissociation energy."""
    Te = De_x - A_WELL
    U1, dU1 = _wall(T11C["A"], A_JOIN)
    y = 1 + np.sqrt(U1 / A_WELL)  # Morse A_WELL (exp(-β(R - Re)) - 1)² matched in value and slope
    beta = -dU1 / (2 * A_WELL * y * (y - 1))
    Re = A_JOIN + log(y) / beta

    def U(R):
        R = np.asarray(R, dtype=float)
        morse = A_WELL * (np.exp(-beta * (R - Re)) - 1) ** 2
        return Te + np.where(R < A_JOIN, _wall(T11C["A"], R)[0], np.where(R < Re, morse, 0.0))

    return U


def c_potential(De_x):
    """C(1u) state, cm⁻¹ above the X minimum."""
    p = T11C["C"]
    Rt = brentq(lambda R: _wall(p, R)[0] - C_TAIL, 3.0, 4.5)
    n = -Rt * _wall(p, Rt)[1] / C_TAIL

    def U(R):
        R = np.asarray(R, dtype=float)
        return De_x + np.where(R < Rt, _wall(p, np.minimum(R, Rt))[0], C_TAIL * (Rt / R) ** n)

    return U


def b_potential(hannover_b):
    """B state for the continuum: Hannover, with the Tellinghuisen wall inside their crossing."""
    V = lambda R: float(hannover_b(np.array([R]))[0])  # noqa: E731
    Te = minimize_scalar(V, bounds=(2.8, 3.3), method="bounded").fun
    wall = lambda R: Te + _wall(T11C["B"], R)[0]  # noqa: E731
    join = brentq(lambda R: wall(R) - V(R), 2.69, 2.74)
    return lambda R: np.where(np.asarray(R) < join, wall(R), hannover_b(np.asarray(R, dtype=float)))


def _moment(coef):
    return lambda R: P.polyval(np.asarray(R) - 2.7, coef)


def _overlaps(V, h, E, c, stride, n_x, window, weighted):
    """|<E|f_k>|² (D²/cm⁻¹) for each column f_k of ``weighted``, sampled every ``stride`` points.

    u(R; E) is the outward Numerov solution of -c u'' + V u = E u on a grid of step h that starts
    deep inside the repulsive wall. It is normalized per unit energy by fitting
    u = k^(-1/2) (a sin φ + b cos φ), φ = ∫k dR, over the grid indices ``window`` and setting
    a² + b² = 1/(π c). Energies classically forbidden anywhere in the window get zero.
    """
    lo, hi = window
    last = max(hi - 1, (n_x - 1) * stride)
    q = h * h / (12 * c)
    samples, win = np.zeros((n_x, E.size)), np.zeros((hi - lo, E.size))
    u_prev, u_cur = np.zeros_like(E), np.full_like(E, 1e-30)
    f_prev, f_cur = 1 + q * (E - V[0]), 1 + q * (E - V[1])
    for i in range(1, last + 1):
        f_next = 1 + q * (E - V[i + 1])
        u_next = ((12 - 10 * f_cur) * u_cur - f_prev * u_prev) / f_next
        if i % stride == 0 and i // stride < n_x:
            samples[i // stride] = u_cur
        if lo <= i < hi:
            win[i - lo] = u_cur
        u_prev, u_cur, f_prev, f_cur = u_cur, u_next, f_cur, f_next
        if i % 64 == 0:
            big = np.abs(u_cur) > 1e100
            if big.any():
                s = np.abs(u_cur[big])
                u_prev[big] /= s
                u_cur[big] /= s
                samples[:, big] /= s
                win[:, big] /= s
    k2 = (E[None, :] - V[lo:hi, None]) / c
    ok = (k2 > 0).all(axis=0)
    k = np.sqrt(np.where(k2 > 0, k2, 1.0))
    phi = h * (np.cumsum(k, axis=0) - 0.5 * (k + k[0]))
    s, co = np.sin(phi) / np.sqrt(k), np.cos(phi) / np.sqrt(k)
    scale = np.abs(win).max(axis=0)
    scale[scale == 0] = 1.0
    win /= scale
    Sss, Scc, Ssc = (s * s).sum(0), (co * co).sum(0), (s * co).sum(0)
    Su, Cu = (s * win).sum(0), (co * win).sum(0)
    det = Sss * Scc - Ssc**2
    a, b = (Scc * Su - Ssc * Cu) / det, (Sss * Cu - Ssc * Su) / det
    norm = np.sqrt(1 / (np.pi * c) / (a * a + b * b)) / scale
    m2 = ((samples * np.where(ok, norm, 0.0)).T @ weighted) ** 2
    return m2


@dataclass
class Continuum:
    """Temperature-independent continuum data; σ(ν, T) from cross_section()."""

    J: np.ndarray  # rotational quantum numbers of the samples (bin centres, spacing dJ)
    dJ: float
    e_lower: np.ndarray  # (nJ, nv) X term values above X(0, 0)
    energies: dict  # state -> upper-state energy grid, cm⁻¹ above X(0, 0)
    m2: dict  # state -> (nJ, nE, nv) |<E'|μ|v''J>|², D²/cm⁻¹
    cut: np.ndarray  # (nJ,) B energy (above X(0, 0)) below which the line list holds the absorption
    isotopologue: str

    def populations(self, T):
        """Fraction of molecules in each (J sample, v'') at temperature T."""
        w = self.dJ * (2 * self.J[:, None] + 1) * np.exp(-C2 * self.e_lower / T)
        return w / w.sum()

    def cross_section(self, nu, T, states=("A", "B", "C")):
        """Continuum cross section (cm²) at wavenumbers nu (cm⁻¹) and temperature T."""
        nu = np.asarray(nu, dtype=float)
        grid = np.arange(nu.min() - GRID_STEP, nu.max() + 2 * GRID_STEP, GRID_STEP)
        pop = self.populations(T)
        total = np.zeros_like(grid)
        for state in states:
            E, m2 = self.energies[state], self.m2[state]
            for j, v in zip(*np.nonzero(pop > 1e-12)):
                Ep = grid + self.e_lower[j, v]
                y = np.interp(Ep, E, m2[j, :, v], left=0.0, right=0.0)
                if state == "B":
                    y[Ep < self.cut[j]] = 0.0
                total += pop[j, v] * y
        return np.interp(nu, grid, S_UNIT * grid * total)

    def at(self, T):
        """Callable ν -> σ(ν) at temperature T, e.g. for spectrum.apparent_cross_section(continuum=...)."""
        return lambda nu: self.cross_section(nu, T)

    def save(self, path):
        arrays = dict(J=self.J, dJ=self.dJ, e_lower=self.e_lower, cut=self.cut, isotopologue=self.isotopologue)
        for s in self.m2:
            arrays[f"E_{s}"], arrays[f"m2_{s}"] = self.energies[s], self.m2[s]
        np.savez(path, **arrays)

    @classmethod
    def load(cls, path):
        z = np.load(path)
        states = [k[3:] for k in z.files if k.startswith("m2_")]
        return cls(z["J"], float(z["dJ"]), z["e_lower"], {s: z[f"E_{s}"] for s in states},
                   {s: z[f"m2_{s}"] for s in states}, z["cut"], str(z["isotopologue"]))


def continuum_model(isotopologue="127I2", upper_cut=None, *, T_max=600.0, dJ=5, nu_max=26500.0,
                    parameters=DEFAULT_PARAMETERS, cache=True) -> Continuum:
    """Continuum for T <= T_max and ν <= nu_max (cm⁻¹).

    upper_cut: MasterLineList.upper_cut of the line list the continuum complements, or None for B←X
    above the dissociation limit only. dJ: spacing of the J'' samples. Cached in cache_dir().
    """
    if isotopologue not in ISOTOPOLOGUES:
        raise ValueError(f"unknown isotopologue {isotopologue!r}")
    cut_hash = hashlib.sha1(np.ascontiguousarray(upper_cut, dtype=float).tobytes()).hexdigest() if upper_cut is not None else None
    key = hashlib.sha1(repr((_CACHE_VERSION, isotopologue, parameters, float(T_max), float(dJ), float(nu_max),
                             cut_hash)).encode()).hexdigest()[:16]
    path = cache_dir() / f"continuum_{isotopologue}_{key}.npz"
    if cache and path.exists():
        return Continuum.load(path)

    pots = load_potentials(parameters)
    mu = reduced_mass(isotopologue)
    ratio = reduced_mass(REFERENCE_ISOTOPOLOGUE) / mu
    c = HBAR2_2U / mu
    E_pop = T_max * log(1 / POP_MIN) / C2
    X = RadialSolver(pots["X"], mu, rmin=RMIN, rmax=RMAX, step=STEP, nlev=80, mass_ratio=ratio)
    e00 = X.energy(0, 0)
    nv = int((X.levels(0) - e00 < E_pop).sum())
    X = RadialSolver(pots["X"], mu, rmin=RMIN, rmax=RMAX, step=STEP, nlev=nv, mass_ratio=ratio)

    J = []
    while X.energy(0, (dJ - 1) / 2 + len(J) * dJ) - e00 < E_pop:
        J.append((dJ - 1) / 2 + len(J) * dJ)
    J = np.array(J)
    if upper_cut is None:
        cut = np.full(J.size, pots["B"].De - e00)
    else:
        known = np.flatnonzero(np.isfinite(upper_cut))
        cut = np.interp(J, known, np.asarray(upper_cut)[known], right=-np.inf)

    De_x = pots["X"].De
    upper = {"A": (a_potential(De_x), _moment(T11C["A"]["mu"]), De_x - A_WELL - e00),
             "B": (b_potential(pots["B"]), mu_tellinghuisen2011, np.nanmin(np.where(np.isfinite(cut), cut, np.nan)) - 10),
             "C": (c_potential(De_x), _moment(T11C["C"]["mu"]), De_x - e00)}
    h, stride = STEP / 4, 4
    n_fine = int(round((max(RMAX, max(w[1] for w in WINDOWS.values())) - RMIN) / h)) + 3
    R = RMIN + np.arange(n_fine) * h
    energies = {s: np.arange(e_lo, nu_max + E_pop, E_STEPS[s]) for s, (_, _, e_lo) in upper.items()}
    m2 = {s: np.zeros((J.size, energies[s].size, nv), dtype=np.float32) for s in upper}
    e_lower = np.zeros((J.size, nv))
    for j, Jj in enumerate(J):
        e_x, c_x = X.states(Jj)
        e_lower[j] = e_x - e00
        for s, (U, moment, _) in upper.items():
            window = tuple(int(round((w - RMIN) / h)) for w in WINDOWS[s])
            V = U(R) - e00 + c * Jj * (Jj + 1) / R**2
            weighted = (moment(X.R) * np.sqrt(STEP))[:, None] * c_x
            m2[s][j] = _overlaps(V, h, energies[s], c, stride, X.R.size, window, weighted)
    result = Continuum(J, float(dJ), e_lower, energies, m2, cut, isotopologue)
    if cache:
        path.parent.mkdir(parents=True, exist_ok=True)
        result.save(path)
    return result
