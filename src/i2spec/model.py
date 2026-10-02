"""Rovibronic level energies and B-X transition wavenumbers."""

from __future__ import annotations

import numpy as np

from . import hfs_params
from .constants import (DEFAULT_PARAMETERS, HBAR2_2U, ISOTOPOLOGUES, MHZ_PER_CM, NUCLEAR_SPIN, REFERENCE_ISOTOPOLOGUE,
                        reduced_mass)
from .hyperfine import level_structure, line_components
from .bspline import BSplineSolver
from .potentials import load_extended, load_potentials, parameter_set
from .level_corrections import corrections_for
from .solver import RadialSolver

# Grids must start above R = -b Rm (2.27 Å for B), where the BOC functions have a pole.
DEFAULT_GRIDS = {
    # The X grid reaches v'' = 69: nesterenko2019 and matyugin2012 measure emission transitions down
    # to v'' = 48, 53 and 54, which need both the larger box and the extra levels. Extending it
    # changes v'' <= 24 by under 0.01 kHz, at about 5x the cost per solve.
    "bspline": {"X": dict(rmin=2.10, rmax=6.0, h=0.01, order=10, nlev=70),
                "B": dict(rmin=2.35, rmax=8.0, h=0.01, order=10, nlev=60)},
    "dvr": {"X": dict(rmin=2.10, rmax=4.0, step=0.004, nlev=25),
            "B": dict(rmin=2.35, rmax=8.0, step=0.005, nlev=60)},
}
SOLVERS = {"bspline": BSplineSolver, "dvr": RadialSolver}

#: Smallest B-spline box that reproduces the full-box levels to under 1 kHz, per highest v needed.
#: Solving every J on a box large enough for the most extreme level we ever use is wasteful: the
#: cost is O(n^3), and almost all data need only low v. Verified against rmax = 6.0 (X) and 8.0 (B).
ADAPTIVE_GRIDS = {
    "X": [(17, dict(rmin=2.10, rmax=3.6, h=0.01, order=10, nlev=20)),
          (24, dict(rmin=2.10, rmax=4.0, h=0.01, order=10, nlev=30)),
          (33, dict(rmin=2.10, rmax=4.4, h=0.01, order=10, nlev=40)),
          (54, dict(rmin=2.10, rmax=5.0, h=0.01, order=10, nlev=60)),
          (69, dict(rmin=2.10, rmax=6.0, h=0.01, order=10, nlev=70)),
          # to the X dissociation limit (Martin 1986 reaches v'' = 108): a 40 A box graded beyond 5 A, as for B
          (115, dict(rmin=2.10, rmax=40.0, h=0.01, order=10, nlev=118, mesh=((5.0, 0.02), (8.0, 0.05), (15.0, 0.1))))],
    "B": [(17, dict(rmin=2.35, rmax=4.0, h=0.01, order=10, nlev=26)),
          (33, dict(rmin=2.35, rmax=4.5, h=0.01, order=10, nlev=40)),
          (43, dict(rmin=2.35, rmax=5.5, h=0.01, order=10, nlev=55)),
          (59, dict(rmin=2.35, rmax=8.0, h=0.01, order=10, nlev=60)),
          # to the dissociation limit: every level bound by more than 0.3 cm-1 (v' <= 73 at J = 0, 72 at
          # J = 26, 61 at J = 100) agrees with a 30 A box to < 0.1 MHz (docs/design/near-dissociation.md)
          (73, dict(rmin=2.35, rmax=12.0, h=0.01, order=10, nlev=80)),
          # the last bound levels (the Orsay atlas Partie IV reaches v' = 79): a 40 A box, graded beyond 5 A
          # where the local wavelength grows; equal to a uniform 0.01 A, 30 A box to < 1e-4 MHz for every
          # level bound by more than 0.3 cm-1 (docs/research/orsay-atlas-19700-20035.md)
          # Above the asymptote the box also holds continuum states outside the centrifugal barrier; nlev is
          # large enough that the resonances behind it (v' <= 77 at J' <= 96, docs/research/quasibound-b.md)
          # survive the cut, and _far_levels keeps only the states localised inside the barrier.
          (90, dict(rmin=2.35, rmax=40.0, h=0.01, order=10, nlev=260, mesh=((5.0, 0.02), (8.0, 0.05), (15.0, 0.1))))],
}


def resonance_ladder(solver, potential, J, mu, mass_ratio, inside=0.9):
    """The levels of ``solver`` at J with the box states of the continuum removed: every state below the
    asymptote, and above it those with at least ``inside`` of their probability inside the centrifugal
    barrier (quasi-bound resonances; their tunnelling widths are < 1e-8 MHz for the levels the Orsay atlas
    Partie IV measured). v then counts physical levels, whatever the box puts between them."""
    E, psi = solver.wavefunctions(J)
    E = np.asarray(E, dtype=float)
    asymptote = float(potential(np.array([400.0]))[0])
    R = solver.R
    veff = potential(R) + potential.adiabatic(R, mass_ratio) + \
        HBAR2_2U / mu * (1 + potential.nonadiabatic(R, mass_ratio)) * J * (J + 1) / R**2
    outer = R > 4.5
    if not outer.any() or veff[outer].max() <= asymptote:
        return E[E < asymptote]                    # no barrier above the asymptote: all such states are box states
    Rb = R[outer][int(np.argmax(veff[outer]))]
    w = solver.W[:, None] * psi**2
    frac = (w * (R[:, None] < Rb)).sum(axis=0) / w.sum(axis=0)
    return E[(E < asymptote) | (frac >= inside)]


def grid_for(state: str, v_max: int) -> dict:
    """The cheapest validated B-spline grid that is exact for levels up to ``v_max``."""
    for limit, grid in ADAPTIVE_GRIDS[state]:
        if v_max <= limit:
            return grid
    raise ValueError(f"no validated {state} grid reaches v = {v_max}")


class RovibronicModel:
    """Hyperfine-free B-X model of one isotopologue from a published parameter set.

    solver: "bspline" (default) puts knots at the potential joins and converges to ~0.1 kHz.
    "dvr" (sinc-DVR) gives eigenvectors on a grid, which the intensity code needs, but is only good
    to ~0.3 MHz for levels near the joins.
    """

    def __init__(self, isotopologue="127I2", parameters=DEFAULT_PARAMETERS, grids=None, solver="bspline",
                 potentials=None):
        if isotopologue not in ISOTOPOLOGUES:
            raise ValueError(f"unknown isotopologue {isotopologue!r}; choose from {sorted(ISOTOPOLOGUES)}")
        self.isotopologue = isotopologue
        self.parameters = parameters
        self.solver = solver
        #: ``potentials`` replaces one or both published curves, for trying a new form (Phase C).
        self.potentials = {**load_potentials(parameters), **(potentials or {})}
        mu = reduced_mass(isotopologue)
        ratio = reduced_mass(REFERENCE_ISOTOPOLOGUE) / mu
        self._mu, self._ratio, self._far = mu, ratio, {}
        self.grids = {**DEFAULT_GRIDS[solver], **(grids or {})}
        self.states = {s: SOLVERS[solver](self.potentials[s], mu, mass_ratio=ratio, **self.grids[s]) for s in ("X", "B")}
        self._reference = None
        #: Measured level corrections (level_corrections.py), named by the parameter set; ¹²⁷I₂ only, and
        #: not when the potentials are replaced (a new potential form must earn its own).
        self.corrections = (corrections_for(parameters)
                            if isotopologue == REFERENCE_ISOTOPOLOGUE and not potentials else None)
        #: The set's extended-range potentials (an MLR pair): levels from ``from_v`` up come from them.
        #: Not when a potential is replaced by the caller, which is a study of that potential alone.
        ext, self.from_v = load_extended(parameters) if not potentials else (None, None)
        self.extended_potentials = ext
        self.extended = ({s: SOLVERS[solver](ext[s], mu, mass_ratio=ratio, **self.grids[s]) for s in ("X", "B")}
                         if ext else None)
        self._levels: dict[tuple[str, int], np.ndarray] = {}
        self._hfs_levels: dict[tuple, list] = {}
        self._hfs_table = None

    def levels(self, state: str, J: int) -> np.ndarray:
        """Term values (cm⁻¹) of v = 0, 1, ... at J, with the parameter set's level corrections applied."""
        key = (state, J)
        if key not in self._levels:
            e = np.array(self.states[state].levels(J), dtype=float)
            if self.extended is not None:
                v0 = self.from_v[state]
                e[v0:] = np.asarray(self.extended[state].levels(J), dtype=float)[v0:len(e)]
            if self.corrections is not None:
                e = e + self.corrections.shift_levels(state, J, len(e)) / MHZ_PER_CM
            self._levels[key] = e
        return self._levels[key]

    def energy(self, state: str, v: int, J: int) -> float:
        """Term value (cm⁻¹) above the X-state potential minimum, level corrections included.

        A level above the default grid's reach (B v' >= 60) is taken from a larger validated grid,
        built on first use (ADAPTIVE_GRIDS), so the default solve stays cheap.
        """
        if v >= len(self.levels(state, J)):
            return float(self._far_levels(state, J, v)[v])
        if self.corrections is None and self.extended is None:
            return self.states[state].energy(v, J)
        return float(self.levels(state, J)[v])

    def _far_levels(self, state: str, J: int, v: int) -> np.ndarray:
        """levels() on the cheapest validated grid that reaches v (B-spline only)."""
        if self.solver != "bspline":
            raise IndexError(f"{state} v = {v} is beyond this {self.solver} grid ({len(self.levels(state, J))} levels)")
        grid = grid_for(state, v)
        key = (state, grid["rmax"], grid["nlev"])
        if key not in self._far:
            mk = lambda pot: SOLVERS["bspline"](pot, self._mu, mass_ratio=self._ratio, **grid)   # noqa: E731
            self._far[key] = (mk(self.potentials[state]),
                              mk(self.extended_potentials[state]) if self.extended is not None else None, {})
        solver, ext, cache = self._far[key]
        if J not in cache:
            e = np.array(solver.levels(J), dtype=float)
            if ext is not None:
                v0 = self.from_v[state]
                e_ext = np.asarray(ext.levels(J), dtype=float)
                if "mesh" in grid:                       # the dissociation grid: resonances, not box states
                    e_ext = resonance_ladder(ext, self.extended_potentials[state], J, self._mu, self._ratio)
                n = min(len(e), len(e_ext))
                e = np.r_[e[:v0], e_ext[v0:n]]
            if self.corrections is not None:
                e = e + self.corrections.shift_levels(state, J, len(e)) / MHZ_PER_CM
            cache[J] = e
        return cache[J]

    def _reference_term(self, state: str, v: int) -> float:
        """¹²⁷I₂ term value of (v, J = 0) above X(0, 0), cm⁻¹: the energy in the hyperfine formulae.

        From the published curves alone, without level corrections or the extended-range potentials:
        the BKT02/S06 formulae are empirical functions of Hannover's own term values, and the
        hyperfine table was fitted on top of them with those energies.
        """
        if self._reference is None:
            self._reference = (self if self.isotopologue == REFERENCE_ISOTOPOLOGUE else
                               RovibronicModel(REFERENCE_ISOTOPOLOGUE, self.parameters, self.grids, self.solver))
        r = self._reference
        if v < r.states[state].nlev:
            return r.states[state].energy(v, 0) - r.states["X"].energy(0, 0)
        key = ("raw", state)          # beyond the default grid: the published curve on a larger validated grid
        if key not in r._far:
            r._far[key] = SOLVERS["bspline"](r.potentials[state], r._mu, mass_ratio=r._ratio, **grid_for(state, v))
        return r._far[key].energy(v, 0) - r.states["X"].energy(0, 0)

    def transition(self, v_upper: int, v_lower: int, J_lower: int, branch: str) -> float:
        """Wavenumber (cm⁻¹) of the P or R line (v'-v'') with lower-state J''."""
        if branch == "R":
            J_upper = J_lower + 1
        elif branch == "P":
            J_upper = J_lower - 1
        else:
            raise ValueError("branch must be 'P' or 'R' (Q lines are forbidden for 0u+ - 0g+)")
        nu = self.energy("B", v_upper, J_upper) - self.energy("X", v_lower, J_lower)
        if self.corrections is not None and self.corrections.band is not None:
            nu += self.corrections.band_shift(v_upper, v_lower, J_lower, branch) / MHZ_PER_CM
        return nu

    def hyperfine_components(self, v_upper: int, v_lower: int, J_lower: int, branch: str, dJ: int = 2,
                             corrections=None, table="default", position=True):
        """Hyperfine components of one B-X line.

        Returns (nu0, components): the hyperfine-free line frequency in MHz, and the components,
        whose offsets (MHz) are relative to nu0. dJ = 2 includes the J ± 2 couplings needed for
        kHz accuracy; dJ = 0 is faster and good to ~1 MHz, enough for Doppler-limited spectra.
        Parameters (hfs_params.line_states): Bodermann et al. (2002) for ¹²⁷I₂ with v' <= 43;
        Salumbides et al. (2006) otherwise, scaled by nuclear moments for ¹²⁹I₂ and ¹²⁷I¹²⁹I, with
        one C per nucleus. Above v' = 53 the B-state parameters are frozen at that energy.

        On top of those formulae, ``table`` adds the measured B-state corrections of hfs_table for
        ¹²⁷I₂: "default" uses the table the parameter set names, None the bare formulae (what a fit of corrections
        to the formulae must see), or pass an hfs_table.HyperfineTable.

        position=False skips the line frequency (nu0 is None). With dJ = 0 the components then need no
        rotational level energies, only the J = 0 term values the parameters depend on, which is what makes a
        pattern for every line of a list cheap.
        """
        if table == "default":
            table = self._default_hfs_table()
        nu0 = self.transition(v_upper, v_lower, J_lower, branch) * MHZ_PER_CM if position else None
        J_upper = J_lower + 1 if branch == "R" else J_lower - 1
        a, b = ISOTOPOLOGUES[self.isotopologue]
        x_params, b_params = hfs_params.line_states(self.isotopologue, v_upper, v_lower,
                                                    self._reference_term("B", v_upper),
                                                    self._reference_term("X", v_lower), corrections, table)
        eqQ_ratio, C_ratio = hfs_params.nucleus_ratios(self.isotopologue)
        spins = dict(i1=NUCLEAR_SPIN[a], i2=NUCLEAR_SPIN[b], eqQ_ratio=eqQ_ratio, C_ratio=C_ratio, dJ=dJ)
        lower = self._hyperfine_levels("X", v_lower, J_lower, x_params, "g" if a == b else None, spins,
                                       memo=corrections is None)
        upper = self._hyperfine_levels("B", v_upper, J_upper, b_params, "u" if a == b else None, spins,
                                       memo=corrections is None)
        return nu0, line_components(upper, lower, J_upper, J_lower)

    def _hyperfine_levels(self, state, v, J, params, symmetry, spins, memo=True):
        """hyperfine.level_structure of (state, v, J), kept for the other lines that share the level (the P
        and R lines of a band, and every band from or to it). The key holds the parameters themselves, so
        a level is reused only where it would come out the same; not at all when ``memo`` is False (a fit,
        whose corrections change every call)."""
        energy = lambda Jn: self.energy(state, v, Jn) * MHZ_PER_CM          # noqa: E731
        if not memo:
            return level_structure(J, params, energy, symmetry=symmetry, **spins)
        key = (state, v, J, symmetry, tuple(sorted(spins.items())),
               tuple(params(Jn) for Jn in range(J - spins["dJ"], J + spins["dJ"] + 1, 2) if Jn >= 0))
        if key not in self._hfs_levels:
            self._hfs_levels[key] = level_structure(J, params, energy, symmetry=symmetry, **spins)
        return self._hfs_levels[key]

    def _default_hfs_table(self):
        """The hyperfine table the parameter set names (hfs_table.default_table), looked up once."""
        if self._hfs_table is None:
            from .hfs_table import default_table
            from .potentials import parameter_set
            self._hfs_table = default_table(parameter_set(self.parameters).get("hyperfine_table", "b_state_lines"))
        return self._hfs_table
