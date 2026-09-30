"""Radial Schrödinger equation in a B-spline basis, with knots at the potential's joins.

The Hannover potentials are only C¹ where their pieces join (V'' jumps at R_I and R_O). That limits
the convergence of grid methods such as the sinc-DVR in solver.py. Here the basis has a knot at
every join with continuity C³, the regularity of the exact eigenfunctions there, and matrix
elements are integrated piecewise. High-order convergence is restored.
"""

from __future__ import annotations

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.interpolate import BSpline
from scipy.linalg import eigh

from .constants import HBAR2_2U


def breakpoints(rmin, rmax, h, joins=(), mesh=()):
    """Nearly uniform breakpoints with spacing <= h that include every join inside (rmin, rmax).

    ``mesh`` coarsens the long range: ((r1, h1), (r2, h2), ...) uses spacing <= h1 beyond r1, h2 beyond
    r2, and so on (r increasing). Near dissociation the local wavelength grows from 0.05 A in the well
    to ~2 A at 10 A, so a box that holds the last levels need not be uniformly fine.
    """
    steps = [(rmin, h), *[(r, hr) for r, hr in mesh if rmin < r < rmax]]
    edges = sorted({rmin, rmax, *(r for r, _ in steps), *(j for j in joins if rmin < j < rmax)})
    points = [rmin]
    for a, b in zip(edges[:-1], edges[1:]):
        step = [hr for r, hr in steps if r <= a][-1]
        points.extend(np.linspace(a, b, max(1, int(np.ceil((b - a) / step))) + 1)[1:])
    return np.array(points)


class BSplineSolver:
    """Bound levels E(v, J), in cm⁻¹, of one electronic state; same interface as solver.RadialSolver.

    Basis: B-splines of the given order (degree order-1) on breakpoints spaced <= h. Knot multiplicity
    is order-4 at the joins (continuity C³) and 1 elsewhere. ψ(rmin) = ψ(rmax) = 0. Matrix elements use
    Gauss-Legendre quadrature on every knot interval. The effective radial Hamiltonian is that of
    solver.RadialSolver.
    """

    def __init__(self, potential, mu, *, rmin, rmax, h=0.02, order=10, nlev=60, mass_ratio=1.0, joins=None,
                 mesh=()):
        if joins is None:      # a potential without joins (an MLR, say) is smooth: no repeated knots
            joins = (potential.RI, potential.RO) if hasattr(potential, "RI") else ()
        joins = [j for j in joins if rmin < j < rmax]
        bp = breakpoints(rmin, rmax, h, joins, mesh)
        interior = []
        for x in bp[1:-1]:
            interior += [x] * (order - 4 if any(abs(x - j) < 1e-12 for j in joins) else 1)
        t = np.r_[[rmin] * order, interior, [rmax] * order]
        n_basis = len(t) - order

        xg, wg = leggauss(order + 4)
        a, b = bp[:-1, None], bp[1:, None]
        self.R = (0.5 * (b - a) * xg + 0.5 * (a + b)).ravel()
        self.W = (0.5 * (b - a) * wg).ravel()
        values, slopes = [], []
        for i in range(1, n_basis - 1):  # dropping the end functions enforces ψ = 0 at rmin and rmax
            c = np.zeros(n_basis)
            c[i] = 1.0
            spline = BSpline(t, c, order - 1, extrapolate=False)
            values.append(np.nan_to_num(spline(self.R)))
            slopes.append(np.nan_to_num(spline.derivative()(self.R)))
        self._B = np.array(values).T
        dB = np.array(slopes).T
        W = self.W[:, None]
        self.S = self._B.T @ (W * self._B)
        self._T = HBAR2_2U / mu * (dB.T @ (W * dB))
        V = potential(self.R) + potential.adiabatic(self.R, mass_ratio)
        rot = HBAR2_2U / mu * (1 + potential.nonadiabatic(self.R, mass_ratio)) / self.R**2
        self._V = self._B.T @ (W * V[:, None] * self._B)
        self._rot = self._B.T @ (W * rot[:, None] * self._B)
        self.nlev = nlev
        self._cache: dict[int, np.ndarray] = {}

    def levels(self, J: int) -> np.ndarray:
        """Lowest ``nlev`` eigenvalues (cm⁻¹) for rotational quantum number J.

        LAPACK's simple driver, which returns every eigenvalue, is twice as fast here as the expert
        one asked for a subset (9.0 against 17.6 ms at 397 basis functions), and agrees with it to
        4e-6 MHz. Taking all of them and slicing is the cheaper way to get the lowest few.
        """
        if J not in self._cache:
            H = self._T + self._V + J * (J + 1) * self._rot
            self._cache[J] = eigh(H, self.S, eigvals_only=True, driver="gv")[:self.nlev]
        return self._cache[J]

    def energy(self, v: int, J: int) -> float:
        return float(self.levels(J)[v])

    def wavefunctions(self, J: int):
        """(E, psi) for rotational level J: energies in cm⁻¹ and psi on the quadrature grid.

        psi has shape (len(self.R), nlev) and is normalised so that sum(W psi²) = 1, which is what
        the Hellmann–Feynman derivative ∂E/∂p = <psi|∂V/∂p|psi> expects (see i2spec.fitting).
        """
        H = self._T + self._V + J * (J + 1) * self._rot
        # divide and conquer, again for every eigenpair: 13.3 ms against 22.6 for the subset driver,
        # agreeing to 2e-6 MHz in energy and 1e-12 in psi^2
        energies, c = eigh(H, self.S, driver="gvd")
        energies, c = energies[:self.nlev], c[:, :self.nlev]
        self._cache.setdefault(J, energies)
        return energies, self._B @ c
