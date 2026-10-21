"""balanced-expanded: research-balanced continued to X v''<=108 and B v'<=80 (release 20261020).

Every coefficient of research-20261018-balanced is kept (read from that release, SHA256-checked). Beyond its
data the MLR exponent and rotational function continue as their own Taylor polynomials plus Bernstein tail
terms, X gains a C6 change, a local cubic-spline correction on 4.2-6.6 A and, shared with B, a dissociation
energy change. Levels are solved on the fit's own grids (one grid per range of v, long boxes for the highest
levels, quasi-bound levels kept inside their centrifugal barrier).

Hyperfine: balanced's own constants in its domain (B v'<=43, X v''<=17). Beyond, tabulated constants of the
second-order model (B: anchored background + second-order sums over the nine 2P3/2 + 2P1/2 states; X: eqQ
continued to the separated-atom value), and at B v' = 77, 78 the 1g(1Pig) level that Pique et al. 1986 found
near resonance, diagonalized with B. Tables hold every level the potentials bind (or hold behind their
barrier) up to J = 203.
"""
from __future__ import annotations

from math import comb, factorial, sqrt
import numpy as np
from numpy.polynomial.chebyshev import chebder, chebval
from scipy.interpolate import BSpline
from scipy.linalg import eigh

from .bspline import BSplineSolver
from .constants import HBAR2_2U, MHZ_PER_CM, reduced_mass
from .hyperfine import (HyperfineLevel, HyperfineLevels, HyperfineParameters, _coefficient_terms, build_blocks,
                        level_structure)
from .wigner import wigner_3j as w3j, wigner_6j as w6j

IA = 2.5


class _Tail:
    """A Chebyshev series in z continued past z_k by its Taylor polynomial (order k) plus a_j t^(k+1) B_j(t)."""

    def __init__(self, c, shift, scale, zk, side, k, a, end=None, vanish=0):
        c = np.asarray(c, float)
        self.shift, self.scale, self.zk, self.side, self.k = shift, scale, zk, side, k
        self.a, self.end, self.vanish = np.asarray(a, float), end, vanish
        x = (zk-shift)/scale
        self.taylor = [chebval(x, chebder(c, i) if i else c)/scale**i/factorial(i) for i in range(k+1)]
        self.span = abs((end if end is not None else side)-zk)

    def region(self, z):
        return (z-self.zk)*self.side > 0

    def __call__(self, z):
        t = np.where(self.region(z), (z-self.zk)*self.side/self.span, 0.)
        if self.end is not None:
            t = np.minimum(t, 1.)
        m = len(self.a)
        dz = z-self.zk
        out = sum(d*dz**i for i, d in enumerate(self.taylor))
        for j, aj in enumerate(self.a):
            out = out+aj*t**(self.k+1)*(1-t)**self.vanish*comb(m-1, j)*t**j*(1-t)**(m-1-j)
        return out


class _ExpandedPotential:
    """V(R) and the rotational function of one state: balanced's MLR with its tails (BSplineSolver protocol)."""

    def __init__(self, state, base, structure, values):
        cfg, frozen = base["config"], base["values"][state]
        self.state, self.cfg, self.frozen = state, cfg, frozen
        self.values = values[state]
        self.C = {int(n): c for n, c in base["inherited_potential_coefficients"][state].items()}
        knots, orders = structure["knots"], structure["orders"][state]
        re = np.exp(frozen["logRe"])
        q = cfg["qx" if state == "X" else "qb"]
        self.ref = re*cfg["ref_x" if state == "X" else "ref_b"]
        self.yq = lambda R: (R**q-self.ref**q)/(R**q+self.ref**q)
        n = cfg["nx" if state == "X" else "nb"]
        c = [frozen[f"c{i}"] for i in range(n)]
        shift, scale = (.2, .6) if state == "X" and cfg["mapped_x"] else (0., 1.)
        k, v = structure["continuity"][state], self.values
        self.tails = [_Tail(c, shift, scale, self.yq(knots[f"{state}_out"]), 1, k,
                            [v[f"to{j}"] for j in range(1, orders["outer"]+1)], vanish=structure["vanish"][state])]
        if f"{state}_in" in knots:
            self.tails.append(_Tail(c, shift, scale, self.yq(knots[f"{state}_in"]), -1, k,
                                    [v[f"ti{j}"] for j in range(1, orders["inner"]+1)],
                                    end=self.yq(structure["inner_end"])))
        self.c, self.shift, self.scale = np.asarray(c), shift, scale
        na = cfg["alpha_x" if state == "X" else "alpha_b"]
        self.alpha = [frozen[f"alpha{i}"] for i in range(na)]
        self.ashift, self.ascale = (.3, .4) if cfg["mapped_alpha"] else (0., 1.)
        kr = knots[f"{state}_rot"]
        self.rotation = _Tail(self.alpha, self.ashift, self.ascale, (kr-re)/(kr+re), 1, structure["rotation_continuity"],
                              [v[f"ro{j}"] for j in range(1, orders["rotation"]+1)])
        self.dispersion = [int(n) for n in structure["dispersion"][state]]
        self.new_cn_unit, self.de_unit = structure["new_cn_unit"], structure["de_unit"]
        self.free_de = structure["free_de"][state]
        self.window = structure["window_knots"][state]
        self.window_unit = structure["window_unit"]
        if self.window:
            k = np.asarray(self.window, float)
            self._splines = [BSpline.basis_element(k[i:i+5], extrapolate=False) for i in range(len(k)-4)]
            self._w = np.array([v[f"w{j}"] for j in range(1, len(self._splines)+1)])

    def __call__(self, R):
        R = np.asarray(R, float)
        cfg, frozen = self.cfg, self.frozen
        re, de = np.exp(frozen["logRe"]), np.exp(frozen["logDe"])
        p = cfg["px" if self.state == "X" else "pb"]
        u = sum(c/R**n for n, c in self.C.items())
        ue = float(sum(c/re**n for n, c in self.C.items()))
        y = (R**p-self.ref**p)/(R**p+self.ref**p)
        yq = self.yq(R)
        yp = (R**p-re**p)/(R**p+re**p)
        outer = self.tails[0]
        if self.dispersion:
            S = np.where(outer.region(yq), ((yq-outer.zk)/outer.span)**(outer.k+1), 0.)
            u = u+S*sum(self.values[f"dC{n}" if n in self.C else f"C{n}"]*(self.C.get(n, self.new_cn_unit))/R**n
                        for n in self.dispersion)
        exponent = chebval((yq-self.shift)/self.scale, self.c)
        for tail in self.tails:
            exponent = np.where(tail.region(yq), tail(yq), exponent)
        beta = y*np.log(2*de/ue)+(1-y)*exponent
        shape = u/ue*np.exp(np.clip(-beta*yp, -100, 100))
        V = frozen.get("Te", 0.)+de*(1-shape)**2
        if self.free_de:
            switch = np.where(outer.region(yq), np.minimum((yq-outer.zk)/outer.span, 1.)**(outer.k+1), 0.)
            V = V+self.values["dDe"]*self.de_unit*switch
        if self.window:
            basis = np.nan_to_num(np.stack([s(R) for s in self._splines], axis=-1))
            V = V+self.window_unit*basis@self._w
        return V

    def adiabatic(self, R, mass_ratio=1.):
        if mass_ratio != 1.:
            raise ValueError("research models are restricted to reference-isotope 127I2")
        return np.zeros_like(np.asarray(R, float))

    def nonadiabatic(self, R, mass_ratio=1.):
        if mass_ratio != 1.:
            raise ValueError("research models are restricted to reference-isotope 127I2")
        R = np.asarray(R, float)
        re = np.exp(self.frozen["logRe"])
        ya = (R-re)/(R+re)
        A = chebval((ya-self.ashift)/self.ascale, self.alpha)
        A = np.where(self.rotation.region(ya), self.rotation(ya), A)
        return 2*re/(R+re)*A


def _sign(x):
    return -1. if round(x) % 2 else 1.


def _element(k, Jp, Ip, J, I, F, omega, parity, eps, value):
    """<1g eps' J' I' F|H^k(a) + H^k(b)|B 0u+ J I F> per unit vibrational overlap, MHz (Chen et al. 2005,
    eqs 4, 5, 9, 10)."""
    dJ, dI = Jp-J, Ip-I
    s_b = _sign(dI) if parity == "u" else -_sign(dI)
    if 1+s_b == 0.:
        return 0.
    f = _sign(omega)*value/w3j(k, IA, IA, 0, IA, -IA)
    if omega:
        f *= (1+eps*_sign(dJ))/sqrt(2.)
    angular = (_sign(F+dI+2*IA+dJ+k)*sqrt((2*Jp+1)*(2*J+1)*(2*Ip+1)*(2*I+1))*w3j(Jp, k, J, -omega, omega, 0)
               * w6j(Jp, k, J, I, F, Ip)*w6j(k, IA, IA, IA, Ip, I))
    return angular*f*(1+s_b)


class ExpandedLevels:
    """E(state, v, J) of the expanded potentials on the fit's grids (one grid per range of v)."""

    def __init__(self, base, structure, values):
        self.potentials = {s: _ExpandedPotential(s, base, structure, values) for s in ("X", "B")}
        self.grids = structure["grids"]
        self.joins = structure["joins"]
        self._solvers, self._levels = {}, {}
        self.asymptote = {s: float(p(np.array([400.]))[0]) for s, p in self.potentials.items()}

    def tier(self, state, v):
        for i, (limit, _) in enumerate(self.grids[state]):
            if v <= limit:
                return i
        raise ValueError(f"no grid holds {state} v={v}")

    def _kept(self, state, tier, J):
        key = state, tier, J
        if key not in self._levels:
            grid = dict(self.grids[state][tier][1])
            if (state, tier) not in self._solvers:
                g = {k: (tuple(map(tuple, x)) if k == "mesh" else x) for k, x in grid.items() if k != "nlev"}
                self._solvers[state, tier] = BSplineSolver(self.potentials[state], reduced_mass("127I2"),
                                                           joins=self.joins[state], nlev=grid["nlev"], **g)
            solver, pot = self._solvers[state, tier], self.potentials[state]
            E, psi = solver.wavefunctions(J)
            if "mesh" in grid:
                # the fit engine's rule: bound levels and quasi-bound resonances (>= 90% inside the barrier)
                R, asym = solver.R, self.asymptote[state]
                w = solver.W[:, None]*psi**2
                veff = pot(R)+J*(J+1)*HBAR2_2U/reduced_mass("127I2")/R**2*(1+pot.nonadiabatic(R))
                outer = R > 4.5
                if veff[outer].max() > asym:
                    Rb = R[outer][int(np.argmax(veff[outer]))]
                    frac = (w*(R[:, None] < Rb)).sum(0)/w.sum(0)
                    E = E[(E < asym) | (frac >= .9)]
                else:
                    E = E[E < asym]
            self._levels[key] = E
        return self._levels[key]

    def energy(self, state, v, J):
        E = self._kept(state, self.tier(state, v), J)
        if v >= len(E):
            raise ValueError(f"balanced-expanded has no {state} v={v} level at J={J} (dissociated)")
        return float(E[v])

    def densities(self, state, v, J):
        """(R, W psi^2 of the kept levels, solver) on the grid that holds v: what _kept keeps, in its order."""
        tier = self.tier(state, v)
        E = self._kept(state, tier, J)
        solver, pot = self._solvers[state, tier], self.potentials[state]
        energies, psi = solver.wavefunctions(J)
        density = solver.W[:, None]*psi**2
        keep = np.isin(energies, E)
        return solver.R, density[:, keep], tier


class ExpandedHyperfine:
    """Constants (eqQ MHz; C, d, delta kHz) of the out-of-core levels, and the B v' 77/78 resonance data."""

    def __init__(self, tables):
        self.B = {(int(r[0]), int(r[1])): r[2:] for r in tables["B"]["rows"]}
        self.X = {(int(r[0]), int(r[1])): r[2:] for r in tables["X"]["rows"]}
        self.x_spin_spin = tables["X"]["spin_spin"]
        self.resonance = {int(v): dict(r, constants={int(J): c for J, c in r["constants"].items()},
                                       overlaps={tuple(map(int, k.split(","))): s for k, s in r["overlaps"].items()})
                          for v, r in tables["resonance"]["levels"].items()}
        self.resonance_state = tables["resonance"]

    def constants(self, state, v, J):
        table = self.B if state == "B" else self.X
        if (v, J) not in table:
            raise ValueError(f"balanced-expanded has no hyperfine constants for {state} v={v} J={J}")
        c = table[(v, J)]
        if state == "B":
            return HyperfineParameters(eqQ=c[0], C=c[1], d=c[2], delta=c[3])
        return HyperfineParameters(eqQ=c[0], C=c[1], **self.x_spin_spin)


def resonance_structure(hf, levels, v, J):
    """HyperfineLevels of B (v, J) with the near-resonant 1g(1Pig) level diagonalized in (P86 II, sect. 3.3):
    B J-2..J+2 blocks with their constants less that level's own second-order term, plus the 1g level's
    rotational levels J-2..J+2 (both Omega-doubling components) with its own hyperfine, coupled through the
    1g(1Pig) electronic couplings times the vibrational overlaps. The vectors keep B's components only."""
    r = hf.resonance[v]
    omega, parity, _, v1, v2 = r["coupling"]
    values = {1: v1, 2: v2}
    Jb = [j for j in (J-2, J, J+2) if j >= 0]
    for j in Jb:
        if j not in r["constants"]:
            raise ValueError(f"balanced-expanded's resonance tables lack B v={v} J={j}")
    E0, Ezero = levels.energy("B", v, J), levels.energy("B", v, 0)

    def params(j):
        c = r["constants"][j]
        return HyperfineParameters(eqQ=c[0], C=c[1], d=c[2], delta=c[3])
    rot = lambda j: (levels.energy("B", v, j)-E0)*MHZ_PER_CM         # noqa: E731
    blocks = build_blocks(Jb, params, rot, IA, IA, "u")
    out = []
    for F, (basis, HB) in blocks.items():
        gbasis = []
        for Jg in range(max(J-2, 1), J+3):
            for eps in (1, -1):
                partners = [(jb, ib) for jb, ib in basis if eps*(-1)**(Jg-jb) == 1]
                if not partners:
                    continue
                parity_I = (partners[0][1]+1) % 2
                for Ig in range(int(2*IA)+1):
                    if Ig % 2 == parity_I and abs(Jg-Ig) <= F <= Jg+Ig:
                        gbasis.append((Jg, eps, Ig))
        nb, ng = len(basis), len(gbasis)
        H = np.zeros((nb+ng, nb+ng))
        H[:nb, :nb] = HB
        for a, (Jg, eps, Ig) in enumerate(gbasis):
            Eg = (Ezero-E0)*MHZ_PER_CM+r["dE_MHz"]+r["B_MHz"]*Jg*(Jg+1)
            for b, (Jg2, eps2, Ig2) in enumerate(gbasis):
                if Jg2 != Jg or eps2 != eps:
                    continue
                c = _coefficient_terms(Jg, Ig, Jg, Ig2, F, IA, IA, 1., 1.)
                eqq = hf.resonance_state["eqQ_1g_MHz"]*(1-3/(Jg*(Jg+1)))
                H[nb+a, nb+b] += c[2]*eqq+c[0]*hf.resonance_state["C_1g_MHz_JJ1"]/(Jg*(Jg+1))*1e3
            H[nb+a, nb+a] += Eg
            for i, (jb, ib) in enumerate(basis):
                if (Jg, jb) not in r["overlaps"]:
                    raise ValueError(f"balanced-expanded's resonance tables lack the overlap {Jg}-{jb} at v={v}")
                x = sum(_element(k, Jg, Ig, jb, ib, F, omega, parity, eps, values[k]) for k in (1, 2))
                H[i, nb+a] = H[nb+a, i] = x*r["overlaps"][(Jg, jb)]
        w, U = np.linalg.eigh(H)
        for k in range(nb+ng):
            b_part = U[:nb, k]
            weight = float(b_part@b_part)
            if weight < 1e-6:
                continue
            central = sum(b_part[i]**2 for i, (jb, _) in enumerate(basis) if jb == J)
            if central < .5*weight and weight > .5:
                continue                    # a B level of J +- 2: it belongs to that rotational level
            dom = int(np.argmax(np.abs(U[:, k])))
            I = basis[dom][1] if dom < nb else gbasis[dom-nb][2]
            out.append(HyperfineLevel(F, I, float(w[k]), basis, b_part))
    return HyperfineLevels(sorted(out, key=lambda lv: lv.energy))


class ExpandedUncertainty:
    """1-sigma of balanced-expanded's predictions (release 20261020 uncertainty.json; recipe in its `recipe`).

    sigma_stat: the linearized covariance of the fit, Cov = F F^T over the free parameters (new terms, frozen core
    potential coefficients, inherited dispersion log factors, calibrations), propagated through each level's
    Hellmann-Feynman derivatives, dE/dp = <psi| dV/dp + J(J+1) hbar^2/(2 mu R^2) dq/dp |psi>, the potential
    derivatives by central differences. A line's statistical sigma is |a_B - a_X| with a = F^T dE/dp.
    sigma_model, sigma_HF and the data flags are the frozen smooth functions and tables of the same file."""

    def __init__(self, table, base, structure, values):
        stat = table["statistical"]
        self.names = list(stat["parameters"])
        self.F = np.asarray(stat["factor"], float)
        self.step = float(stat["fd_step"])
        self.base, self.structure, self.values = base, structure, values
        self.table = table
        self._derivatives, self._factors = {}, {}

    def _perturbed(self, state, name, delta):
        """The state's potential with one free parameter moved by delta (dDe moves both states' shift)."""
        import copy
        base, values = copy.deepcopy(self.base), copy.deepcopy(self.values)
        s, n = name.split(":", 1)
        if n.startswith("logC"):
            power = n[4:]
            base["inherited_potential_coefficients"][s][power] *= np.exp(delta)
        elif n in base["values"].get(s, {}):
            base["values"][s][n] += delta
        else:
            for t in (("X", "B") if n == "dDe" else (s,)):
                values[t][n] += delta
        return _ExpandedPotential(state, base, self.structure, values)

    def _stepsize(self, name):
        s, n = name.split(":", 1)
        if n.startswith("logC"):
            value = 0.
        elif n in self.base["values"].get(s, {}):
            value = self.base["values"][s][n]
        else:
            value = self.values[s][n]
        return self.step*max(1., abs(value))

    def _potential_derivatives(self, state, tier, R):
        """dV/dp and dq/dp (q the rotational function) on a grid, one column per free parameter."""
        key = state, tier
        if key not in self._derivatives:
            dV, dq = np.zeros((len(R), len(self.names))), np.zeros((len(R), len(self.names)))
            for j, name in enumerate(self.names):
                s, n = name.split(":", 1)
                if s == "calibration" or (s != state and n != "dDe"):
                    continue
                h = self._stepsize(name)
                plus, minus = self._perturbed(state, name, h), self._perturbed(state, name, -h)
                dV[:, j] = (plus(R)-minus(R))/(2*h)
                dq[:, j] = (plus.nonadiabatic(R)-minus.nonadiabatic(R))/(2*h)
            self._derivatives[key] = dV, dq
        return self._derivatives[key]

    def level_factor(self, levels, state, v, J):
        """a = F^T dE/dp (MHz) of one level."""
        tier = levels.tier(state, v)
        key = state, tier, J
        if key not in self._factors:
            R, density, _ = levels.densities(state, v, J)
            dV, dq = self._potential_derivatives(state, tier, R)
            rot = J*(J+1)*HBAR2_2U/reduced_mass("127I2")/R**2
            self._factors[key] = (density.T@(dV+rot[:, None]*dq))*MHZ_PER_CM@self.F
        factors = self._factors[key]
        if v >= len(factors):
            raise ValueError(f"balanced-expanded has no {state} v={v} level at J={J} (dissociated)")
        return factors[v]

    # the frozen smooth functions and tables

    def _near(self, kind, state, v, J):
        runs = self.table["coverage"][kind].get(f"{state}:{v}", ())
        margin = self.table["coverage"]["margin_J"]
        return any(lo-margin <= J <= hi+margin for lo, hi in runs)

    def _span(self, state, v):
        runs = self.table["coverage"]["all"].get(f"{state}:{v}")
        return (runs[0][0], runs[-1][1]) if runs else None

    def status(self, state, v, J):
        """measured: training rows within 15 in J at this v; interpolated: inside this v's J range, or between
        neighbours within 2 in v whose ranges (+-15) hold J; else extrapolated."""
        if self._near("all", state, v, J):
            return "measured"
        margin, reach = self.table["coverage"]["margin_J"], self.table["coverage"]["neighbor_v"]
        span = self._span(state, v)
        if span and span[0] <= J <= span[1]:
            return "interpolated"
        holds = lambda w: (s := self._span(state, w)) is not None and s[0]-margin <= J <= s[1]+margin   # noqa: E731
        if any(holds(v-k) for k in range(1, reach+1)) and any(holds(v+k) for k in range(1, reach+1)):
            return "interpolated"
        return "extrapolated"

    def line_class(self, core, upper, lower):
        """'precise' (both core levels within 15 in J of precision rows at their v), 'atlas', or 'beyond'."""
        (vu, Ju), (vl, Jl) = upper, lower
        if vu > core["B"] or vl > core["X"]:
            return "beyond"
        return "precise" if self._near("precise", "B", vu, Ju) and self._near("precise", "X", vl, Jl) else "atlas"

    def model_sigma(self, core, state, v, J, line_class):
        """One level's share of sigma_model (MHz)."""
        m = self.table["model_error"]
        if v <= core[state]:
            return m[f"core {'atlas' if line_class == 'beyond' else line_class}"]["sigma_MHz"]/sqrt(2.)
        f = m[f"{state} beyond"]
        J = min(J, f["J_clip"])
        return float(np.exp(f["log_sigma0"]+f["slope_per_v"]*(v-f["v0"])+f["J2_coefficient"]*(J/100)**2))

    def kappa(self, core, state, v):
        """Fractional error of the hyperfine pattern's scale from one state's levels (held beyond the validated v)."""
        hf = self.table["hyperfine"]
        if v <= core[state]:
            return hf["core"]["kappa"]
        if state == "B":
            b = hf["B beyond"]
            return float(np.exp(b["log_kappa0"]+b["power"]*np.log1p((min(v, b["v_range"][1])-b["v_ref"])/b["v_scale"])))
        x = hf["X beyond"]
        return x["kappa_17"]*np.exp(x["slope_per_v"]*(min(v, x["validated_v"][1])-17))

    def hyperfine_terms(self, core, vu, vl):
        """(a MHz, kappa) of a line: the beyond-core values of each state, the core's for a core line."""
        hf = self.table["hyperfine"]
        kc = hf["core"]["kappa"]
        k2 = self.kappa(core, "B", vu)**2+self.kappa(core, "X", vl)**2-kc**2
        a = hf["B beyond"]["a_MHz"] if vu > core["B"] else hf["core"]["a_MHz"]
        return a, sqrt(max(k2, kc**2))

    def hyperfine_status(self, state, v):
        for top, text in self.table["hyperfine"]["status"][state]:
            if v <= top:
                return text
        return "extrapolated"


def install(model, spec, tables, base):
    """Attach the expanded levels and hyperfine to a ResearchModel instance (see research._ExpandedModel)."""
    model._levels_x = ExpandedLevels(base.coefficients, spec["structure"], spec["values"])
    model._hf = ExpandedHyperfine(tables)
    model._base = base


def structure(model, state, v, J):
    core = model._spec["core"]
    if state == "B" and v in model._hf.resonance:
        return resonance_structure(model._hf, model._levels_x, v, J)
    energy = model._energy(state, v, J)
    if v <= core[state]:
        couplings = lambda j: model._base._couplings(state, v, j)      # noqa: E731
    else:
        couplings = lambda j: model._hf.constants(state, v, j)        # noqa: E731
    return level_structure(J, couplings, lambda j: (model._energy(state, v, j)-energy)*MHZ_PER_CM,
                           symmetry="g" if state == "X" else "u", dJ=2)

