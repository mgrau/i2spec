"""Global least-squares refit of the potentials to the observation sets (M5, Phase B).

What is fitted: the power-series coefficients a_i of the X and B "X-representation" potentials. The
inner wall and the outer branch are not free — after every change to a_i they are re-derived from C¹
continuity at R_I and R_O (``XRepPotential.with_continuous_extensions``), so the potential stays
continuous by construction rather than by constraint.

What is held fixed, and why:

* **The hyperfine parameters.** ``bodermann1998b`` measures both an absolute frequency and a
  hyperfine splitting of the same lines at v'' = 16 and 17: the splitting is reproduced to 36 kHz
  while the positions are out by 2.5-5.8 MHz. The error is in the term values, not the hyperfine
  Hamiltonian, so letting the hyperfine constants vary in this stage would only let them absorb
  potential error. Each observation's hyperfine offset is therefore computed once from the starting
  model and carried as a constant; ``refresh_hyperfine`` recomputes it between outer iterations,
  which is enough because the offsets depend on the potentials only through rotational spacings.
* **The Born-Oppenheimer corrections**, for the same reason one stage earlier: the isotopologue
  data is a handful of lines and would not constrain them here.

The Jacobian is analytic, by Hellmann-Feynman:

    dE(v,J)/dp = <psi_vJ| dV/dp |psi_vJ>

so one eigensolve per (isotopologue, state, J) yields the derivative with respect to *every*
parameter, instead of one eigensolve per parameter. Measured on the X state: 412x faster than a
one-sided numerical Jacobian, and the two agree to about 1e-5 relative. dV/dp itself is taken by
central differences of the potential on the quadrature grid, which costs two function evaluations
per parameter and keeps the re-derived wall and tail correctly inside the derivative.
"""

from __future__ import annotations

import time
from dataclasses import dataclass, replace

import numpy as np
from scipy.optimize import least_squares

from .constants import MHZ_PER_CM, REFERENCE_ISOTOPOLOGUE, reduced_mass, DEFAULT_PARAMETERS
from .bspline import BSplineSolver
from .model import DEFAULT_GRIDS, RovibronicModel, grid_for
from .observations import component_rank, load_all
from .potentials import load_potentials


#: Relative step used for dV/dp on the quadrature grid.
DERIVATIVE_STEP = 1e-6

#: Long-range parameters of the X state, with the 1-sigma uncertainties of the measurement they come
#: from: Bacis, Cerny & Martin, J. Mol. Spectrosc. 118, 434 (1986), which is Knoeckel 2004's ref. 47.
#: Its abstract gives De = 12 547.335(128) cm-1, C6 = 1.48(12)e6, C8 = 3.86(1.20)e7, C10 = 1.0(0.5)e8.
#: C8 is known to 31% and C10 to 50%, so these are genuinely loose priors, not decoration.
LONG_RANGE_PRIORS = {"X": {"De": 0.128, "C6": 0.12e6, "C8": 1.20e7, "C10": 0.5e8}}


def read_knob(pot, name):
    """Value of a fitted potential parameter: 'a3', 'De' or 'C6'."""
    if name.startswith("a"):
        return pot.a[int(name[1:])]
    if name == "De":
        return pot.De
    if name.startswith("C"):
        return pot.C[int(name[1:])]
    raise KeyError(name)


def set_knob(pot, name, value):
    """Copy of ``pot`` with one fitted parameter changed; extensions are re-derived by the caller."""
    if name.startswith("a"):
        a = list(pot.a)
        a[int(name[1:])] = value
        return replace(pot, a=tuple(a))
    if name == "De":
        return replace(pot, De=float(value))
    if name.startswith("C"):
        C = dict(pot.C)
        C[int(name[1:])] = float(value)
        return replace(pot, C=C)
    raise KeyError(name)


@dataclass(frozen=True)
class Term:
    """One signed level energy entering a residual: coefficient * E(isotopologue, state, v, J)."""

    isotopologue: str
    state: str
    v: int
    J: int
    coefficient: float


@dataclass
class Datum:
    """One observation reduced to a linear combination of level energies plus a constant."""

    key: str
    group: str
    dataset: str
    value: float          # MHz
    uncertainty: float    # MHz
    terms: tuple[Term, ...]
    constant: float = 0.0  # hyperfine offsets, MHz

    def model(self, energies) -> float:
        total = self.constant
        for t in self.terms:
            total += t.coefficient * energies[(t.isotopologue, t.state, t.J)][t.v]
        return total


@dataclass
class FitResult:
    parameters: dict
    x: np.ndarray
    x0: np.ndarray
    cost: float
    chi2_per_datum: float
    rms_before: dict
    rms_after: dict
    covariance: np.ndarray | None
    success: bool
    message: str
    iterations: int
    seconds: float


class GlobalFit:
    """Least-squares refit of the X and B series coefficients to a set of observations."""

    def __init__(self, datasets=None, parameters=DEFAULT_PARAMETERS, fit_states=("X", "B"),
                 grids=None, max_uncertainty=None, isotopologues=None,
                 prior_tau=30.0, fixed=(("X", 0),), shape_sigma=0.01, shape_points=120,
                 calibration_sigma=0.3, min_calibration_group=5, calibration_bound=3.0,
                 long_range=("X",)):
        """prior_tau sets the Gaussian prior on the coefficients, centred on the published values.

        It is expressed in the natural metric of the problem: each coefficient is allowed to move the
        data by about prior_tau standard deviations on its own, i.e. sigma_prior_i = prior_tau/||J_i||.
        In the scaled units the optimiser uses this makes the prior block exactly (1/prior_tau) times
        the identity, which is perfectly conditioned - a plain relative prior is not, because the
        coefficients differ by ten orders of magnitude in how strongly they move the levels.

        ``fixed`` names coefficients held at their published values. The default holds X.a[0], which
        is exactly 0 in the published set because it defines the energy origin: shifting it moves
        every X level, so it is perfectly degenerate with B.a[0] (the electronic term energy) and
        only the difference is observable. Fitting both wastes a direction and lets a meaningless
        gauge mode drift.

        ``shape_sigma``/``shape_points`` add a prior in *function* space: V(R) is pinned to the
        published curve at a grid of radii, with a width that opens up where the data actually have
        probability density and closes to ``shape_sigma`` (cm-1) where they do not. Without it the
        coefficient prior alone does not stop the fit moving V by more than a wavenumber in regions
        carrying 1e-22 of the density, because the series is global: adjusting the well necessarily
        wags the extrapolated wall and tail.

        ``calibration_sigma`` gives every data group with at least ``min_calibration_group`` points a
        free additive offset in MHz, with that Gaussian prior. Different laboratories' frequency
        chains disagree at this level -- dube2004a documents its own scale systematic -- and without
        the offsets a constant inter-laboratory difference is indistinguishable from a defect of the
        model. ``calibration_bound`` (MHz) is a HARD bound on them, and it is not optional: the
        starting residuals at v'' = 53 and 54 are about 5e5 sigma, and against those a soft Gaussian
        prior is no constraint at all -- a free linear parameter per group is simply cheaper than
        bending the potential, so the offsets absorb the physics and every data set gets worse. An
        offset of more than a few MHz is not a credible frequency-chain error in any case.

        ``long_range`` names the states whose dissociation energy and dispersion coefficients are
        fitted alongside the series. This is what lets the fit reach high v'' at all: with only the
        series coefficients free, the response of v'' = 53 to the parameters is 94% collinear with
        that of v'' <= 17, so raising one drags the other and no improving direction exists. De and
        C6, C8, C10 act where the high-v'' wavefunctions have amplitude and the low-v'' ones have
        none. They carry the uncertainties of the measurement they come from (LONG_RANGE_PRIORS).

        The prior is not a numerical convenience. The data constrain only 29 of the 46 directions
        (v' = 19-31 and 38-57 and v'' = 7-10 have no data at all), and the published values carry the
        information from the Knoeckel 2004 fit set, which no longer exists (see
        docs/research/data-availability.md). The prior is how that lost data still constrains the
        result. Set prior_tau=None for an unregularised fit."""
        self.datasets = list(datasets) if datasets is not None else load_all()
        self.parameters = parameters
        self.fit_states = tuple(fit_states)
        self.grids = {**DEFAULT_GRIDS["bspline"], **(grids or {})}
        self.potentials = load_potentials(parameters)
        self._models: dict[str, RovibronicModel] = {}
        self._solver_cache: dict = {}

        # the fitted parameters, as an explicit list of (state, knob) pairs
        self.fixed = tuple(fixed or ())
        self.knobs: list[tuple[str, str]] = []
        for s in self.fit_states:
            for i in range(len(self.potentials[s].a)):
                if (s, i) not in self.fixed:
                    self.knobs.append((s, f"a{i}"))
        self.long_range = tuple(long_range or ())
        self.knob_prior = {}                     # index -> sigma, for knobs with their own prior
        for s in self.long_range:
            for name, sigma in LONG_RANGE_PRIORS.get(s, {}).items():
                if name.startswith("C") and int(name[1:]) not in self.potentials[s].C:
                    continue
                self.knob_prior[len(self.knobs)] = sigma
                self.knobs.append((s, name))
        self.n_knobs = len(self.knobs)

        self.prior_tau = prior_tau
        self.shape_sigma = shape_sigma
        self.shape_points = shape_points
        self.calibration_sigma = calibration_sigma
        self.min_calibration_group = min_calibration_group
        self.calibration_bound = calibration_bound
        # The optimiser works in scaled units u = (a - a0)/scale. The scale is set from the
        # Jacobian's column norms at the starting point (see _set_scale), so that a unit step in u
        # moves the residuals by about one sigma. A relative scale does not work here: the
        # coefficients differ enormously in how strongly they move the levels, and a step that is
        # negligible for one is catastrophic for another.
        self.scale = None
        self._shape = None
        self.prior_sigma = None                 # filled once the scale is known

        self.data: list[Datum] = []
        self.skipped: list[tuple[str, str]] = []
        self._build(max_uncertainty, isotopologues)

        # the parameter vector is [potential knobs, per-group calibration offsets]
        self.n_parameters = self.n_knobs + self.n_offsets
        self.x0 = np.concatenate([
            np.array([read_knob(self.potentials[s], k) for s, k in self.knobs], float),
            np.zeros(self.n_offsets)])
        self.x0_full = self.x0

    # ---------------------------------------------------------------- setting up the residuals

    def _model(self, isotopologue) -> RovibronicModel:
        if isotopologue not in self._models:
            self._models[isotopologue] = RovibronicModel(isotopologue, self.parameters)
        return self._models[isotopologue]

    def _line_terms(self, line, sign):
        """E_B(v', J') - E_X(v'', J''), in MHz, as signed terms."""
        J_upper = line.J_lower + 1 if line.branch == "R" else line.J_lower - 1
        return (Term(line.isotopologue, "B", line.v_upper, J_upper, sign * MHZ_PER_CM),
                Term(line.isotopologue, "X", line.v_lower, line.J_lower, -sign * MHZ_PER_CM))

    def _offset(self, line, component):
        """Hyperfine offset of a component from the line centre, MHz; 0 for the centre itself."""
        if not component:
            return 0.0
        rank = component_rank(component)
        model = self._model(line.isotopologue)
        _, comps = model.hyperfine_components(line.v_upper, line.v_lower, line.J_lower, line.branch)
        main = sorted([c for c in comps if c.label], key=lambda c: c.offset)
        if rank > len(main):
            raise LookupError(f"component {component} beyond the {len(main)} main components")
        return float(main[rank - 1].offset)

    def _build(self, max_uncertainty, isotopologues):
        """Reduce every observation to (linear combination of level energies) + constant."""
        for ds in self.datasets:
            scale = MHZ_PER_CM if ds.meta["unit"] == "cm-1" else 1.0
            for o in ds.observations:
                iso = o.line.isotopologue
                if isotopologues and iso not in isotopologues:
                    continue
                if max_uncertainty and o.uncertainty * scale > max_uncertainty:
                    continue
                try:
                    if o.kind == "frequency":
                        terms = list(self._line_terms(o.line, +1))
                        constant = self._offset(o.line, o.component)
                    else:
                        ref = o.ref_line or o.line
                        constant = self._offset(o.line, o.component) - self._offset(ref, o.ref_component)
                        if ref == o.line:
                            # a splitting inside one line: the rovibronic centres cancel exactly, so
                            # this datum constrains the hyperfine Hamiltonian, not the potentials
                            terms = []
                        else:
                            terms = list(self._line_terms(o.line, +1)) + list(self._line_terms(ref, -1))
                except Exception as exc:  # noqa: BLE001 - record and skip rather than abort the build
                    self.skipped.append((f"{o.line} {o.component or ''}".strip(), str(exc)))
                    continue
                self.data.append(Datum(
                    key=f"{o.line} {o.component or 'centre'}", group=o.group or ds.id, dataset=ds.id,
                    value=o.value * scale, uncertainty=o.uncertainty * scale,
                    terms=tuple(terms), constant=constant))

        # Data with no terms cannot move the potentials while the hyperfine parameters are fixed.
        # They are kept for reporting but excluded from the fit itself.
        self.hyperfine_only = [d for d in self.data if not d.terms]
        self.data = [d for d in self.data if d.terms]

        # every (isotopologue, state, J) the residuals need, and the highest v required at each
        self.levels_needed: dict[tuple[str, str, int], int] = {}
        for d in self.data:
            for term in d.terms:
                k = (term.isotopologue, term.state, term.J)
                self.levels_needed[k] = max(self.levels_needed.get(k, 0), term.v)

        # a free additive offset per data group, for groups large enough to separate it from noise
        from collections import Counter
        counts = Counter(d.group for d in self.data)
        self.offset_groups = ([g for g, n in sorted(counts.items()) if n >= self.min_calibration_group]
                              if self.calibration_sigma else [])
        self.offset_index = {g: i for i, g in enumerate(self.offset_groups)}
        self.n_offsets = len(self.offset_groups)

    # ---------------------------------------------------------------- the forward model

    def _offsets(self, x):
        """The per-group calibration offsets carried at the end of the parameter vector, in MHz."""
        x = np.asarray(x, float)
        return x[self.n_knobs:] if self.n_offsets else np.zeros(0)

    def _potentials(self, x, bump=None):
        """Potentials at parameters x; ``bump`` = (knob index, delta) perturbs one knob."""
        x = np.asarray(x, float)
        pots = dict(self.potentials)
        for i, (s, name) in enumerate(self.knobs):
            v = x[i] + (bump[1] if bump and bump[0] == i else 0.0)
            pots[s] = set_knob(pots[s], name, v)
        for s in set(s for s, _ in self.knobs):
            pots[s] = pots[s].with_continuous_extensions()
        return pots

    def _boxes(self):
        """(isotopologue, state, box) -> the J values solved on that box, and J -> its box.

        Each J is solved on the cheapest grid that is exact for the highest v it needs, rather than
        on one box sized for the most extreme level in the whole data set. Since the eigensolve is
        O(n^3) and almost every J needs only low v, this is where the time goes.
        """
        if getattr(self, "_box_map", None) is not None:
            return self._box_map
        groups, which = {}, {}
        for (iso, state, J), vmax in self.levels_needed.items():
            g = grid_for(state, vmax)
            key = (iso, state, g["rmax"], g["nlev"])
            groups.setdefault(key, (g, []))[1].append(J)
            which[(iso, state, J)] = key
        self._box_map = (groups, which)
        return self._box_map

    def _solvers(self, x):
        pots = self._potentials(x)
        groups, _ = self._boxes()
        out = {}
        for key, (grid, _Js) in groups.items():
            iso, state = key[0], key[1]
            mu = reduced_mass(iso)
            ratio = reduced_mass(REFERENCE_ISOTOPOLOGUE) / mu
            out[key] = BSplineSolver(pots[state], mu, mass_ratio=ratio, **grid)
        return out, pots

    def energies(self, x):
        """{(isotopologue, state, J): array of level energies in cm-1}."""
        solvers, _ = self._solvers(x)
        _, which = self._boxes()
        return {k: solvers[which[k]].levels(k[2]) for k in self.levels_needed}

    def residuals(self, x):
        """(model - observed)/uncertainty for every datum, at physical parameters x.

        Neither prior is included here: this is the data residual used for reporting.
        ``scaled_residuals`` is what the optimiser minimises, and it appends both priors.
        """
        E = self.energies(x)
        off = self._offsets(x)
        out = np.empty(len(self.data))
        for i, d in enumerate(self.data):
            model = d.model(E)
            if self.n_offsets and d.group in self.offset_index:
                model += off[self.offset_index[d.group]]
            out[i] = (model - d.value) / d.uncertainty
        return out

    def _level_derivatives(self, x):
        """dE/d(knob) for every level the residuals need, by Hellmann-Feynman."""
        solvers, pots = self._solvers(x)
        x = np.asarray(x, float)
        dV = {}
        for key, solver in solvers.items():
            state = key[1]
            cols = np.zeros((len(solver.R), self.n_knobs))
            for i, (s, name) in enumerate(self.knobs):
                if s != state:
                    continue
                step = DERIVATIVE_STEP * max(abs(x[i]), 1.0)
                up = self._potentials(x, bump=(i, +step))[state]
                dn = self._potentials(x, bump=(i, -step))[state]
                cols[:, i] = (up(solver.R) - dn(solver.R)) / (2 * step)
            dV[key] = cols
        _, which = self._boxes()
        dE = {}
        for (iso, state, J), vmax in self.levels_needed.items():
            key = which[(iso, state, J)]
            solver = solvers[key]
            _, psi = solver.wavefunctions(J)
            if state in self.fit_states:
                dE[(iso, state, J)] = np.einsum("g,gk,gv->vk", solver.W, dV[key],
                                                psi[:, :vmax + 1] ** 2)
            else:
                dE[(iso, state, J)] = np.zeros((vmax + 1, self.n_knobs))
        return dE, dV, solvers, pots

    def jacobian(self, x):
        """d(data residual)/dx. One eigensolve per (isotopologue, state, J) covers every parameter."""
        dE, _, _, _ = self._level_derivatives(x)
        out = np.zeros((len(self.data), self.n_parameters))
        for i, d in enumerate(self.data):
            for term in d.terms:
                out[i, :self.n_knobs] += term.coefficient * dE[(term.isotopologue, term.state, term.J)][term.v]
            if self.n_offsets and d.group in self.offset_index:
                out[i, self.n_knobs + self.offset_index[d.group]] = 1.0
            out[i] /= d.uncertainty
        return out

    def shape_residuals(self, x):
        """(V_fitted - V_published)/sigma at the pinned radii, per fitted state."""
        shape = self._shape_grid()
        if not shape:
            return np.zeros(0)
        pots = self._potentials(x)
        return np.concatenate([(pots[s](R) - V0) / sig for s, (R, sig, V0) in shape.items()])

    def shape_jacobian(self, x):
        shape = self._shape_grid()
        if not shape:
            return np.zeros((0, self.n_parameters))
        x = np.asarray(x, float)
        blocks = []
        for s, (R, sig, _V0) in shape.items():
            block = np.zeros((len(R), self.n_parameters))
            for i, (state, name) in enumerate(self.knobs):
                if state != s:
                    continue
                step = DERIVATIVE_STEP * max(abs(x[i]), 1.0)
                up = self._potentials(x, bump=(i, +step))[s]
                dn = self._potentials(x, bump=(i, -step))[s]
                block[:, i] = (up(R) - dn(R)) / (2 * step) / sig
            blocks.append(block)
        return np.vstack(blocks)

    def _set_scale(self):
        """scale_i = 1/||J_i||: one unit of u changes the residuals by roughly one sigma."""
        if self.scale is not None:
            return self.scale
        j = np.vstack([self.jacobian(self.x0), self.shape_jacobian(self.x0)])
        norms = np.linalg.norm(j, axis=0)
        floor = max(norms.max() * 1e-12, 1e-300)
        self.scale = 1.0 / np.maximum(norms, floor)
        if self.prior_tau is not None:
            self.prior_sigma = self.prior_tau * self.scale
        return self.scale

    def to_x(self, u):
        """Physical parameters from the scaled vector the optimiser sees."""
        return self.x0_full + self._set_scale() * np.asarray(u, float)

    def scaled_residuals(self, u):
        x = self.to_x(u)
        parts = [self.residuals(x), self.shape_residuals(x)]
        if self.prior_tau is not None:
            pri = np.asarray(u, float) / self.prior_tau
            for i, sigma in self.knob_prior.items():       # measured parameters keep their own prior
                pri[i] = (x[i] - self.x0[i]) / sigma
            parts.append(pri)
        if self.n_offsets:
            parts.append(self._offsets(x) / self.calibration_sigma)
        return np.concatenate(parts)

    def scaled_jacobian(self, u):
        x = self.to_x(u)
        parts = [self.jacobian(x) * self.scale, self.shape_jacobian(x) * self.scale]
        if self.prior_tau is not None:
            blk = np.eye(self.n_parameters) / self.prior_tau
            for i, sigma in self.knob_prior.items():
                blk[i, :] = 0.0
                blk[i, i] = self.scale[i] / sigma
            parts.append(blk)
        if self.n_offsets:
            blk = np.zeros((self.n_offsets, self.n_parameters))
            for i in range(self.n_offsets):
                blk[i, self.n_knobs + i] = self.scale[self.n_knobs + i] / self.calibration_sigma
            parts.append(blk)
        return np.vstack(parts)

    # ---------------------------------------------------------------- reporting and running

    def rms_by_dataset(self, x):
        """Weighted-free rms of (model - observed) in MHz, per data set."""
        E = self.energies(x)
        out = {}
        for name in sorted({d.dataset for d in self.data}):
            d_of = [d for d in self.data if d.dataset == name]
            resid = np.array([d.model(E) - d.value for d in d_of])
            out[name] = float(np.sqrt(np.mean(resid ** 2)))
        return out

    def report(self, x, label=""):
        E = self.energies(x)
        resid = np.array([d.model(E) - d.value for d in self.data])
        sig = np.array([d.uncertainty for d in self.data])
        chi2 = float(np.mean((resid / sig) ** 2))
        print(f"{label:<10} n={len(resid)}  rms {resid.std():9.4f} MHz  "
              f"max |{np.abs(resid).max():9.3f}|  chi2/n {chi2:11.2f}")
        return resid

    def _shape_grid(self):
        """Radii at which V is pinned, and the width of the pin at each.

        The width opens where the fitted levels actually have probability density and closes to
        shape_sigma where they do not, so the prior constrains the extrapolation without fighting
        the data.
        """
        if getattr(self, "_shape", None) is not None:
            return self._shape
        if not self.shape_sigma:
            self._shape = {}
            return self._shape
        solvers, pots = self._solvers(self.x0)
        _, which = self._boxes()
        weight = {}
        for d in self.data:
            for term in d.terms:
                k = (term.isotopologue, term.state, term.J, term.v)
                weight[k] = weight.get(k, 0.0) + 1.0 / d.uncertainty ** 2
        self._shape = {}
        for state in self.fit_states:
            iso = REFERENCE_ISOTOPOLOGUE
            # the widest box in use for this state, so every level has a home on the grid
            key = max([k for k in solvers if k[0] == iso and k[1] == state], key=lambda k: k[2])
            s = solvers[key]
            dens = np.zeros(len(s.R))
            for (i2, st, J, v), w in weight.items():
                if st != state or i2 != iso or v >= s.nlev:
                    continue
                _, psi = s.wavefunctions(J)
                dens += w * psi[:, v] ** 2
            dens *= s.W
            if dens.sum() > 0:
                dens /= dens.max()
            R = np.linspace(s.R.min(), s.R.max(), self.shape_points)
            rho = np.interp(R, np.sort(s.R), dens[np.argsort(s.R)])
            sigma = self.shape_sigma * np.sqrt(1.0 + rho / 1e-6)   # wide where the data look
            self._shape[state] = (R, sigma, pots[state](R).copy())
        return self._shape

    def run(self, stages=((1e4, 25), (100.0, 25), (5.0, 50)), verbose=2):
        """Fit, walking the robust-loss scale down from nearly quadratic to robust.

        The starting point is 14-21 cm-1 (about 5e5 sigma) away from the v'' = 53 and 54 data of
        nesterenko2019, far outside the linear regime. A robust loss applied from the start would
        simply ignore those points; a quadratic loss applied at the end would let the single v' = 58
        datum, which is 41 GHz out because the published B potential fails above v' = 44, dominate
        everything. So f_scale is graduated: each stage starts from the previous stage's solution.

        The first scale is deliberately not enormous. Under soft_l1 the gradient from a residual is
        bounded by f_scale, so a very large first scale would hand the whole fit to the eighteen
        nesterenko2019 points, which carry 2-235 kHz uncertainties and start 5e5 sigma out.
        """
        t0 = time.perf_counter()
        before = self.rms_by_dataset(self.x0)
        scale = self._set_scale()
        lo = np.full(self.n_parameters, -np.inf)
        hi = np.full(self.n_parameters, np.inf)
        if self.n_offsets and self.calibration_bound:
            k = slice(self.n_knobs, self.n_parameters)
            lo[k] = -self.calibration_bound / scale[k]
            hi[k] = +self.calibration_bound / scale[k]
        u = np.zeros(self.n_parameters)
        res = None
        for f_scale, max_nfev in stages:
            res = least_squares(self.scaled_residuals, u, jac=self.scaled_jacobian, method="trf",
                                bounds=(lo, hi), loss="soft_l1", f_scale=f_scale,
                                max_nfev=max_nfev, verbose=verbose)
            u = res.x
            if verbose:
                x = self.to_x(u)
                r = self.residuals(x)
                mhz = r * np.array([d.uncertainty for d in self.data])
                ok = np.abs(mhz) < 100
                print(f"  [stage f_scale={f_scale:g}] bulk rms {mhz[ok].std():.4f} MHz over {ok.sum()} "
                      f"data; worst {np.abs(mhz).max():.0f} MHz", flush=True)
        x = self.to_x(u)
        after = self.rms_by_dataset(x)
        cov = None
        try:
            j = self.scaled_jacobian(u)
            _, sv, VT = np.linalg.svd(j, full_matrices=False)
            keep = sv > sv[0] * 1e-12
            cov = (VT[keep].T / sv[keep] ** 2) @ VT[keep]
        except Exception:  # noqa: BLE001
            pass
        pots = self._potentials(x)
        return FitResult(
            parameters={s: tuple(map(float, pots[s].a)) for s in self.fit_states},
            x=x, x0=self.x0, cost=float(res.cost),
            chi2_per_datum=float(np.mean(self.residuals(x) ** 2)),
            rms_before=before, rms_after=after, covariance=cov,
            success=bool(res.success), message=str(res.message),
            iterations=int(res.nfev), seconds=time.perf_counter() - t0)
