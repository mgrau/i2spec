"""B 0u+ ~ 1g(1Pig) hyperfine resonance at v' = 77 and 78, by diagonalization (Pique et al., J. Physique 47, 1917
(1986) [P86 II], sect. 3.3).

Near the 2P3/2 + 2P1/2 limit a level of the 1g(1Pig) ("1'g") state lies a few GHz from B v' = 77 and 78, close
enough that the second-order sums (second_order) fail (P86 II: patterns no effective constants can
make; the forbidden 1g - X lines at the 78-0 head). Here that one 1g vibrational level is taken out of the
1g(1Pig) resolvent and put into the hyperfine Hamiltonian of B (v', J' - 2 ... J' + 2), with its own rotational
levels J = J' - 2 ... J' + 2 and both Omega-doubling components eps = +-1, F by F (P86 II's 45 x 45 matrices):

  B block: build_blocks over J' and J' +- 2 with B's constants from hyperfine_beyond minus that
           level's own second-order term (projected like the rest);
  1g block: E_1g(J) - E_B(v', J') with E_1g(J) = E_B(v', 0) + dE + B_1g J(J + 1) (P86 II's fitted dE and B_1g),
           and the 1g level's own hyperfine, P86 II's separated-atom values: eqQ = -224 MHz (Omega = 1 form,
           eqQ (1 - 3/J(J + 1)) in the J-diagonal eqQ operator), C = 3886/J(J + 1) MHz;
  B-1g: C05 eq 4 (second_order.element, k = 1, 2) with the 1g(1Pig) couplings of the hyperfine_beyond
        model (separated-atom functions mixed at theta_g) times the vibrational overlap of the model curve's level
        nearest in energy to the resonance.

Inherited (P86 II): dE(v+1, 77) = 4.5 +- 0.5 GHz, B_v+1 = 120 +- 5 MHz; dE(v+2, 78) = 1.8 +- 0.3 GHz,
B_v+2 = 100 +- 10 MHz; the 1g level's eqQ and C. v' = 76 (perturbed by level v, B_v = 155 MHz) has no printed dE
and stays second order. A line's components carry the B character of their upper levels, so the weak extra lines
of 1g-dominated levels come out with their strengths.

"""
from __future__ import annotations

import numpy as np
from scipy.linalg import eigh

from ..constants import MHZ_PER_CM
from ..hyperfine import HyperfineLevel, HyperfineLevels, HyperfineParameters, _coefficient_terms, build_blocks
from . import hyperfine_beyond as hb
from . import second_order as so
from .continuation import SolverPotential

STATE = '1g(1Pig)'
PIQUE = {77: dict(dE=4500., dE_sigma=500., B=120., B_sigma=5.),          # MHz, P86 II sect. 3.3.2-3.3.3
         78: dict(dE=1800., dE_sigma=300., B=100., B_sigma=10.)}
G_EQQ, G_C = -224., 3886.                                               # MHz; C = G_C/J(J+1)


def _phase(psi):
    """+1 or -1: the sign that makes the innermost lobe of a radial function (on its grid, R increasing) positive."""
    first = np.flatnonzero(abs(psi) > .01*np.max(abs(psi)))[0]
    return 1. if psi[first] > 0 else -1.


class Resonance:
    def __init__(self, constants, b, theta, levels=None, cache=None, orders=None):
        """constants: the perturbers' curve constants; b, theta: the hyperfine_beyond fit (orders: its background
        orders where they are not the core fit's)."""
        self.constants, self.b, self.theta, self.cache = constants, np.asarray(b, float), dict(theta), cache
        self.orders = orders
        curves, values = hb.potential()
        bpot = SolverPotential(curves['B'], values['B'])
        limit = float(bpot(np.array([300.]))[0])
        self.radial = so.Radial(bpot, curves['B'].joins, {STATE: hb.perturber_curves(limit, self.constants)[STATE]})
        self.levels = levels or hb.Levels()
        self.coup = so.molecular_couplings(self.theta['u'], self.theta['g'])
        self.model = None
        self._g = {}

    def prepare(self, vJ):
        self.model = hb.Model(self.levels, hb.fill_terms(sorted(vJ), self.constants, cache=self.cache),
                              orders=self.orders)

    def g_levels(self, Jp):
        """(energies cm-1, wavefunction coefficients) of the 1g(1Pig) model curve's levels at rotation J'."""
        if Jp not in self._g:
            s = self.radial.p[STATE]
            self._g[Jp] = eigh(s._T+s._V+(Jp*(Jp+1)-1)*s._rot, s.S)
        return self._g[Jp]

    def overlap(self, n, Jp, v, J):
        """<1g level n at J'|B v J>, each radial function in the phase convention of the effective Hamiltonian
        (positive innermost lobe): its J, J +- 2 mixing assumes <v J|v J +- 2> > 0, which an eigensolver's arbitrary
        signs do not give."""
        s = self.radial.p[STATE]
        _, psi = self.radial.level(J)
        e, c = self.g_levels(Jp)
        return float(c[:, n]@(s._B.T@(s.W*psi[:, v])))*_phase(psi[:, v])*_phase(s._B@c[:, n])

    def resonant_index(self, v, dE):
        """The model level nearest in energy to E_B(v, J=1) + dE (rotation J = 1)."""
        E, _ = self.radial.level(1)
        e, _ = self.g_levels(1)
        return int(np.argmin(np.abs(e-(E[v]+dE/MHZ_PER_CM))))

    def b_constants(self, v, J, n):
        """B constants of (v, J) without the resonant level's second-order term."""
        c = self.model.constants(self.b, self.theta, v, J)
        E, _ = self.radial.level(J)
        G = {}
        for Jp in range(max(J-2, 1), J+3):
            e, _ = self.g_levels(Jp)
            s2 = self.overlap(n, Jp, v, J)**2/((E[v]-e[n])*MHZ_PER_CM)
            G[Jp] = s2
        blocks = {}
        for F in range(max(0, J-5), J+6):
            bl = so.second_order_block(J, F, STATE, lambda k, kk, Jp: G[Jp], couplings=self.coup)
            if len(bl[0]):
                blocks[F] = bl
        fit, _, _ = so.project(J, blocks)
        return c-np.array([fit[p] for p in hb.PARAMS])

    def structure(self, v, J):
        """(HyperfineLevels of B (v, J) with the resonant 1g level, B-only vectors; info)."""
        p = PIQUE[v]
        n = self.resonant_index(v, p['dE'])
        Jb = [j for j in (J-2, J, J+2) if j >= 0]
        consts = {j: self.b_constants(v, j, n) for j in Jb}
        E0 = self.levels.energy(v, J)
        Ezero = self.levels.energy(v, 0)

        def params(j):
            c = consts[j]
            return HyperfineParameters(eqQ=c[0], C=c[1], d=c[2], delta=c[3])
        rot = lambda j: (self.levels.energy(v, j)-E0)*MHZ_PER_CM                 # noqa: E731
        blocks = build_blocks(Jb, params, rot, so.IA, so.IA, 'u')
        omega, parity, v0, v1, v2 = self.coup[STATE]
        values = {1: v1, 2: v2}
        out, info = [], dict(model_level=n, overlaps={})
        for F, (basis, HB) in blocks.items():
            gbasis = []
            for Jg in range(max(J-2, 1), J+3):
                for eps in (1, -1):
                    # partner B levels: J_b with eps (-1)^(Jg - J_b) = +1; their I parity fixes the 1g level's
                    partners = [(jb, ib) for jb, ib in basis if eps*(-1)**(Jg-jb) == 1]
                    if not partners:
                        continue
                    parity_I = (partners[0][1]+1) % 2
                    for Ig in range(int(2*so.IA)+1):
                        if Ig % 2 == parity_I and abs(Jg-Ig) <= F <= Jg+Ig:
                            gbasis.append((Jg, eps, Ig))
            nb, ng = len(basis), len(gbasis)
            H = np.zeros((nb+ng, nb+ng))
            H[:nb, :nb] = HB
            for a, (Jg, eps, Ig) in enumerate(gbasis):
                Eg = (Ezero-E0)*MHZ_PER_CM+p['dE']+p['B']*Jg*(Jg+1)
                for bb, (Jg2, eps2, Ig2) in enumerate(gbasis):
                    if Jg2 != Jg or eps2 != eps:
                        continue
                    c = _coefficient_terms(Jg, Ig, Jg, Ig2, F, so.IA, so.IA, 1., 1.)
                    eqq = G_EQQ*(1-3/(Jg*(Jg+1)))
                    H[nb+a, nb+bb] += c[2]*eqq+c[0]*G_C/(Jg*(Jg+1))*1e3
                H[nb+a, nb+a] += Eg
                for i, (jb, ib) in enumerate(basis):
                    S = info['overlaps'].setdefault((Jg, jb), self.overlap(n, Jg, v, jb))
                    x = sum(so.element(k, Jg, Ig, jb, ib, F, omega, parity, eps, values[k]) for k in (1, 2))*S
                    H[i, nb+a] = H[nb+a, i] = x
            w, U = np.linalg.eigh(H)
            for k in range(nb+ng):
                b_part = U[:nb, k]
                weight = float(b_part@b_part)
                if weight < 1e-6:
                    continue
                central = sum(b_part[i]**2 for i, (jb, _) in enumerate(basis) if jb == J)
                if central < .5*weight and weight > .5:
                    continue                    # a B level of J +- 2: belongs to that rotational level
                dom = int(np.argmax(np.abs(U[:, k])))
                I = basis[dom][1] if dom < nb else gbasis[dom-nb][2]
                out.append(HyperfineLevel(F, I, float(w[k]), basis, b_part))
        info['overlaps'] = {f'{a}-{b}': s for (a, b), s in info['overlaps'].items()}
        return HyperfineLevels(sorted(out, key=lambda lv: lv.energy)), info
