"""Radial Schrödinger equation on a sinc-DVR grid (Colbert & Miller, J. Chem. Phys. 96, 1982 (1992))."""

from __future__ import annotations

import numpy as np
from scipy.linalg import eigh

from .constants import HBAR2_2U


class RadialSolver:
    """Bound rovibrational levels E(v, J), in cm⁻¹, of one electronic state.

    Effective radial Hamiltonian (Salumbides et al. 2008, eq. 5):
    -hbar²/(2 mu) d²/dR² + V(R) + V_corr(R) + hbar² [1 + alpha(R)] J(J+1) / (2 mu R²).
    """

    def __init__(self, potential, mu, *, rmin, rmax, step, nlev, mass_ratio=1.0):
        self.R = np.arange(rmin, rmax + step / 2, step)
        self.nlev = nlev
        k = np.arange(self.R.size)
        d = (k[:, None] - k[None, :]).astype(float)
        with np.errstate(divide="ignore"):
            T = 2.0 * (-1.0) ** d / d**2
        T[k, k] = np.pi**2 / 3
        self._T = HBAR2_2U / mu / step**2 * T
        self._V = potential(self.R) + potential.adiabatic(self.R, mass_ratio)
        self._rot = HBAR2_2U / mu * (1 + potential.nonadiabatic(self.R, mass_ratio)) / self.R**2
        self._cache: dict[int, np.ndarray] = {}

    def levels(self, J: int) -> np.ndarray:
        """Lowest ``nlev`` eigenvalues (cm⁻¹) for rotational quantum number J."""
        if J not in self._cache:
            H = self._T + np.diag(self._V + J * (J + 1) * self._rot)
            self._cache[J] = eigh(H, eigvals_only=True, subset_by_index=[0, self.nlev - 1])
        return self._cache[J]

    def states(self, J: int, points: int | None = None):
        """Lowest ``nlev`` eigenvalues (cm⁻¹) and orthonormal eigenvectors (columns) for J.

        On a DVR grid a vector holds ψ(R_i)·√step, so <a|f(R)|b> = Σ_i a_i f(R_i) b_i.

        ``points`` solves on the first that many grid points only, a box ending there, and pads the
        vectors with zeros: the same as a solver built on the shorter grid. For levels that vanish well
        inside it, that is the full-grid answer at a fraction of the cost (the solve is O(n³)).
        """
        if points is None or points >= self.R.size:
            H = self._T + np.diag(self._V + J * (J + 1) * self._rot)
            energies, vectors = eigh(H, subset_by_index=[0, self.nlev - 1])
            self._cache.setdefault(J, energies)
            return energies, vectors
        H = self._T[:points, :points] + np.diag(self._V[:points] + J * (J + 1) * self._rot[:points])
        energies, inner = eigh(H, subset_by_index=[0, self.nlev - 1])
        vectors = np.zeros((self.R.size, inner.shape[1]))
        vectors[:points] = inner
        return energies, vectors

    def energy(self, v: int, J: int) -> float:
        return float(self.levels(J)[v])
