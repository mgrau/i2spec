"""Second-order hyperfine couplings of B 0u+ through the states of its 2P3/2 + 2P1/2 dissociation limit.

Chen, de Jong & Ye, JOSA B 22, 951 (2005) [C05], after Broyer, Vigue & Lehmann, J. Physique 39, 591 (1978) and
Pique et al., J. Physique 47, 1909 (1986) [P86]. A state p of the limit couples to B through the off-diagonal
rotational operator V0 (k = 0, 1u only), the magnetic dipole hyperfine operator (k = 1) and the electric
quadrupole one (k = 2). In the separated-atom basis their electronic parts are constants, C_k<V^k>, from atomic
hyperfine constants (P86 Table II); V0 carries the radial factor hbar^2/(2 mu R^2).

Matrix element (C05 eqs 4, 5, 9, 10), B = |0u+ v J I F>, p = |Omega' eps' v' J' I' F>:
    <p|H^k(a)+H^k(b)|B> = (-1)^(F+dI+2Ia+dJ+k) [(2J'+1)(2J+1)(2I'+1)(2I+1)]^1/2 (J' k J; -W' W' 0)
                          {J' k J; I F I'} {k Ia Ia; Ia I' I} f_k (1 + s_b),
    f_k = (-1)^W' C_k<V^k> / (k Ia Ia; 0 Ia -Ia) <v'|v> (1/sqrt2)^|eps'| [1 + eps' (-1)^dJ]   (eps' = 0 for W' = 0)
    s_b = (-1)^dI for u states, -(-1)^dI for g states.
The second-order effective Hamiltonian inside B (v, J) is
    H2 = sum_p sum_{J', I'} <B|H|p> <p|H|B> G_p(J'),   G_p(J') = <v J| O (E_vJ - H_p(J'))^-1 O' |v J>,
the resolvent of p's radial Hamiltonian on B's grid: every bound level and the discretized continuum at once.
Projected onto the effective operators of i2spec.hyperfine it gives the second-order C, delta, d, eqQ; C05 eqs
13-16 give the same constants in closed form (chen_constants; tests/test_fit_hyperfine.py checks both).

Perturbers (perturber_curves()), as C05 built them: the nine states of the limit from de Jong, Visscher &
Nieuwpoort, JCP 107, 9046 (1997) inside 7 A (data/theory/dejong1997; their Table III constants, as Morse
curves, since the paper prints no point-wise curves) joined to Saute & Aubert-Frecon, JCP 77, 5639 (1982)
outside (data/theory/saute1982, limit + C5/R^5 + C6/R^6). The Omega = 1 pairs are mixtures of the separated-atom
functions; the first-order C5 matrix (long_range_quadrupole()) reproduces Saute's C5 at this limit and fixes
the long-range mixing at 30 degrees for both pairs.
"""
from __future__ import annotations

import csv
from math import factorial, sqrt
from pathlib import Path

import numpy as np
from scipy.linalg import solve_banded

from ..bspline import BSplineSolver
from ..constants import HBAR2_2U, MHZ_PER_CM, reduced_mass
from ..hyperfine import HyperfineLevel, HyperfineLevels, _coefficient_terms, allowed_spins, build_blocks
from ..wigner import wigner_3j as w3j, wigner_6j as w6j

ROOT = Path(__file__).resolve().parents[3]

IA = 2.5
# B and every perturbing state share one B-spline basis and quadrature (the resolvent needs B's wavefunction in
# the perturber's basis); 25 A holds B's last levels and the decay of the resolvent below the limit
GRID = dict(rmin=2.35, rmax=25., h=.02, order=10)
SQRT2 = sqrt(2.)
# 127I atomic hyperfine constants, MHz: A and B of 2P3/2, and A of 2P1/2 as implied by P86 Table II
# ((A12 +- A32)/2 = 3709, 2882 give 6590.7, 6591.3)
ATOMIC = dict(A32=827.265, A12=6591., B32=1146.356)
# P86 Table II: (Omega, g/u, V0 / sqrt(2J(J+1)) at R = 10 A, beta g_I <V1>, 1/2 eQ <V2>), MHz, separated atoms.
# For the tests: in the spherical-tensor normalization of C05 eq 5 these are, per nucleus, 2 (-sqrt2)^|dW| times
# A<j_q> for k = 1 and 4, -4 sqrt2, 4 sqrt2 (dW = 0, 1, 2) times C_2<V^2_q> for k = 2 (separated_atom()).
TABLE2 = {
    "0g+": (0, "g", 0., 0., -573.),
    "0u-": (0, "u", 0., -2882., 0.),
    "0g-": (0, "g", 0., 3709., 0.),
    "B": (0, "u", 0., 0., -573.),
    "1'g": (1, "g", 0., -3491., 0.),
    "1''g": (1, "g", 0., 1013., 810.),
    "1'u": (1, "u", 170., 5830., 0.),
    "1''u": (1, "u", 100., 1013., 810.),
    "2g": (2, "g", 0., 0., 810.),
    "2u": (2, "u", 0., 0., 810.),
}
R_TABLE = 10.         # A, where P86 evaluate V0


def _sign(x):
    return -1. if round(x) % 2 else 1.


def _separated_atom_states():
    """{name: (Omega', g/u, [(product state, coefficient)])}, the separated-atom functions of 2P3/2 + 2P1/2.

    Product states (j_A, m_A; j_B, m_B), nucleus a on atom A. With c1..c4 = (3/2 1/2; 1/2 -1/2), (3/2 -1/2; 1/2 1/2),
    (1/2 -1/2; 3/2 1/2), (1/2 1/2; 3/2 -1/2): B = (c1 + c2 + c3 + c4)/2 (sigma_v pairs c1-c2 and c3-c4, inversion
    c1-c3 and c2-c4; the phases are fixed, up to labels, by P86's V0 for 1'u, 3/sqrt2 rather than 1/sqrt2).
    Each Omega' > 0 state is given by its +Omega' component, (x + s y)/sqrt2 with y the atom-exchanged x; s = +1
    for u, -1 for g: |1'> from x = (3/2 1/2; 1/2 1/2), |1''> from (3/2 3/2; 1/2 -1/2), |2> from (3/2 3/2; 1/2 1/2).
    """
    def exchange(s):
        return (s[2], s[3], s[0], s[1])

    c = [(1.5, .5, .5, -.5), (1.5, -.5, .5, .5), (.5, -.5, 1.5, .5), (.5, .5, 1.5, -.5)]

    def state(x, parity):
        s = 1. if parity == "u" else -1.
        return [(x, 1/SQRT2), (exchange(x), s/SQRT2)]

    def zero(sign_sigma, parity):
        s = 1. if parity == "u" else -1.
        return [(c[0], .5), (c[1], .5*sign_sigma), (c[2], .5*s), (c[3], .5*s*sign_sigma)]

    states = {"0g+": (0, "g", zero(1, "g")), "0u-": (0, "u", zero(-1, "u")), "0g-": (0, "g", zero(-1, "g")),
              "B": (0, "u", zero(1, "u"))}
    for name, x in (("1'", (1.5, .5, .5, .5)), ("1''", (1.5, 1.5, .5, -.5)), ("2", (1.5, 1.5, .5, .5))):
        for parity in "gu":
            states[name+parity] = (int(x[1]+x[3]), parity, state(x, parity))
    return states


def separated_atom(A32=ATOMIC["A32"], A12=ATOMIC["A12"], B32=ATOMIC["B32"]):
    """Electronic couplings of B 0u+ to each state of the 2P3/2 + 2P1/2 limit, pure separated-atom functions
    (_separated_atom_states()). Returns {name: (Omega', g/u, <J+^el>, C_1<V^1_W'(a)>, C_2<V^2_W'(a)>)}, the last
    two in MHz (spherical components, nucleus a only), <J+^el> = <j+(A) + j+(B)> (the operator of V0, C05 eq 11).
    """
    def ladder(j, m):          # <j m+1|j_+|j m>
        return sqrt(j*(j+1)-m*(m+1))

    def spherical(k, j, mp, q, m):     # <j m'|T^k_q|j m> / <j j|T^k_0|j j>
        return _sign(j-mp)*w3j(j, k, j, -mp, q, m)/w3j(j, k, j, -j, 0, j)

    A = {1.5: A32, .5: A12}
    B = {1.5: B32, .5: 0.}

    def atom_a(k, bra, ket):
        """C_k<bra|V^k_q(a)|ket> between product states: atom A changes, atom B must not."""
        (jA, mA, jB, mB), (jA2, mA2, jB2, mB2) = bra, ket
        if jA != jA2 or jB != jB2 or mB != mB2:
            return 0.
        q = mA-mA2
        if k == 1:          # Ia A <j m'|j_q|j m>, j_q spherical
            return IA*A[jA]*_sign(jA-mA)*w3j(jA, 1, jA, -mA, q, mA2)*sqrt(jA*(jA+1)*(2*jA+1)) if abs(q) <= 1 else 0.
        return B[jA]/4*spherical(2, jA, mA, q, mA2) if abs(q) <= 2 and B[jA] else 0.

    def j_plus(bra, ket):
        (jA, mA, jB, mB), (jA2, mA2, jB2, mB2) = bra, ket
        if jA != jA2 or jB != jB2:
            return 0.
        out = 0.
        if mB == mB2 and mA == mA2+1:
            out += ladder(jA, mA2)
        if mA == mA2 and mB == mB2+1:
            out += ladder(jB, mB2)
        return out

    states = _separated_atom_states()
    ket = states["B"][2]
    out = {}
    for name, (omega, parity, bra) in states.items():
        element = lambda f: sum(a*b*f(x, y) for x, a in bra for y, b in ket)      # noqa: E731
        out[name] = (omega, parity, element(j_plus) if omega == 1 else 0.,
                     element(lambda x, y: atom_a(1, x, y)), element(lambda x, y: atom_a(2, x, y)))
    return out


def couplings():
    """{name: (Omega', g/u, V0 at R_TABLE (MHz, per sqrt(2J(J+1))), C_1<V^1>, C_2<V^2>)} from separated_atom()."""
    rot = HBAR2_2U/reduced_mass("127I2")/R_TABLE**2*MHZ_PER_CM
    return {name: (o, p, j*rot, v1, v2) for name, (o, p, j, v1, v2) in separated_atom().items()}


COUPLINGS = couplings()

# The nine perturbing states of 2P3/2 + 2P1/2 (C05 Fig. 2 labels; de Jong et al. 1997 Table III names). Each is
# a separated-atom function of COUPLINGS, or for the Omega = 1 pairs a rotation of |1'> and |1''> by theta:
# the 1'-like state is cos(theta)|1'> + sin(theta)|1''>, its partner -sin(theta)|1'> + cos(theta)|1''>.
# saute: the row of data/theory/saute1982 (state, chen2005 label); dejong: (state, parentage) of dejong1997.
PERTURBERS = {
    "(3)0g+": dict(basis=("0g+", None), saute="(3)0g+", dejong=("(3)0g+", "")),
    "0g-": dict(basis=("0g-", None), saute="0g-(3Pig)", dejong=("0g-", "3Pig")),
    "(3)0u-": dict(basis=("0u-", None), saute="(3)0u-", dejong=("(3)0u-", "")),
    "1g(1Pig)": dict(basis=("g", 0), saute="1'g = 1g(1Pig)", dejong=("1g", "1Pig")),
    "1g(3Sigmag-)": dict(basis=("g", 1), saute="1''g = 1g(3Sigmag-)", dejong=("1g", "3Sigmag-")),
    "1u(3Sigmau+)": dict(basis=("u", 0), saute="1'u = 1u(3Sigmau+)", dejong=("1u", "3Sigmau+")),
    "C1u": dict(basis=("u", 1), saute="1''u = C:1u(3Sigmau+)", dejong=("C1u", "3Sigmau+")),
    "b'2u": dict(basis=("2u", None), saute="2u", dejong=("b'2u", "3Deltau")),
    "2g(1Deltag)": dict(basis=("2g", None), saute="2g", dejong=("2g", "1Deltag")),
}
# Long-range mixing of each Omega = 1 pair: the eigenvectors of the first-order C5 interaction, 30 degrees
# (|coefficients| sqrt3/2, 1/2) for both pairs (long_range_quadrupole()); P86 use 0.84 (32.9 deg) for 1g, C05
# fit 0.99 (8.1 deg) and 0.92 (23.1 deg) for 1u, in their own phase convention.
THETA_SA = dict(u=30., g=30.)
DATA = ROOT/'data'


def long_range_quadrupole():
    """First-order quadrupole-quadrupole (C5) matrix of the 2P3/2 + 2P1/2 states in separated_atom()'s basis.

    The 2P1/2 atom has no quadrupole moment, so the interaction is pure exchange (A: 3/2 -> 1/2, B: 1/2 -> 3/2),
    one reduced matrix element squared times angular factors: V = sum_m 24/((2-m)!(2+m)!) Q_m(A) Q_-m(B) / R^5
    (Saute & Aubert-Frecon 1982 eq 4). Scaled so that B 0u+ has Saute's C5 = -3.59e5 cm-1 A^5, it reproduces
    their Table II at this limit; the Omega = 1 eigenvectors give the long-range mixing angles (THETA_SA).
    Returns ({name: C5}, {g/u: (eigenvalues, eigenvectors over (1', 1''))}), C5 in 1e5 cm-1 A^5.
    """
    reduced = {(1.5, .5): 1., (.5, 1.5): -1.}         # <j||Q||j'> up to a common factor (hermitian tensor)

    def q_element(j, m, jp, mp, q):
        r = reduced.get((j, jp), 0.)
        return _sign(j-m)*w3j(j, 2, jp, -m, q, mp)*r if r and abs(q) <= 2 else 0.

    def V(x, y):
        return sum(24/(factorial(2-m)*factorial(2+m))*q_element(x[0], x[1], y[0], y[1], m)
                   * q_element(x[2], x[3], y[2], y[3], -m) for m in range(-2, 3))

    states = {n: s[2] for n, s in _separated_atom_states().items()}
    el = lambda a, b: sum(ca*cb*V(x, y) for x, ca in states[a] for y, cb in states[b])      # noqa: E731
    scale = -3.59/el("B", "B")
    c5 = {n: scale*el(n, n) for n in states}
    pairs = {}
    for p in "gu":
        M = scale*np.array([[el(a+p, b+p) for b in ("1'", "1''")] for a in ("1'", "1''")])
        pairs[p] = np.linalg.eigh(M)
    return c5, pairs


def molecular_couplings(theta_u=THETA_SA["u"], theta_g=THETA_SA["g"]):
    """COUPLINGS of the nine PERTURBERS, the Omega = 1 pairs rotated by theta (degrees). Every coupling is a
    matrix element <p|O|B>, linear in p, so a mixture's couplings are the same mixture of COUPLINGS."""
    out = {}
    for name, spec in PERTURBERS.items():
        first, which = spec["basis"]
        if which is None:
            out[name] = COUPLINGS[first]
            continue
        t = np.radians(theta_u if first == "u" else theta_g)
        c = (np.cos(t), np.sin(t)) if which == 0 else (-np.sin(t), np.cos(t))
        a, b = COUPLINGS["1'"+first], COUPLINGS["1''"+first]
        out[name] = (a[0], a[1], *(c[0]*x+c[1]*y for x, y in zip(a[2:], b[2:])))
    out["B"] = COUPLINGS["B"]
    return out


OMEGA = {name: c[0] for name, c in {**COUPLINGS, **molecular_couplings()}.items()}


def _rows(path):
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def perturber_constants():
    """{name: dict(C5, C6 (cm-1 A^5, A^6), De (cm-1), re (A), we (cm-1), and the experimental De, re, we where
    de Jong et al. print them)} from data/theory/saute1982 and dejong1997."""
    saute = {r["chen2005"]: r for r in _rows(DATA/"theory/saute1982/data.csv") if r["chen2005"]}
    dejong = {(r["state"], r["parentage"]): r for r in _rows(DATA/"theory/dejong1997/data.csv")
              if r["limit"] == "P3/2+P1/2"}
    ev = 8065.54393734921       # cm-1 per eV (CODATA 2018)
    out = {}
    for name, spec in PERTURBERS.items():
        s, d = saute[spec["saute"]], dejong[spec["dejong"]]
        out[name] = dict(C5=float(s["C5_1e5_cm_A5"])*1e5, C6=float(s["C6_1e6_cm_A6"])*1e6,
                         De=float(d["De_eV"])*ev, re=float(d["re_A"]), we=float(d["we_cm"]),
                         exp=dict(De=float(d["De_exp_eV"])*ev if d["De_exp_eV"] else None,
                                  re=float(d["re_exp_A"]) if d["re_exp_A"] else None,
                                  we=float(d["we_exp_cm"]) if d["we_exp_cm"] else None))
    return out


def perturber_curves(limit, join=7., width=.5, experimental=False, overrides=None):
    """{name: V(R)} of the PERTURBERS, cm-1 on B's energy scale (limit = B's dissociation energy there).

    C05's construction: inside `join` the de Jong et al. 1997 curve, outside Saute & Aubert-Frecon's
    limit + C5/R^5 + C6/R^6, the inner curve shifted vertically to meet the outer one at `join`. De Jong et al.
    print only spectroscopic constants, so the inner curve here is the Morse function of their De, re, we
    (experimental where they print them, if `experimental`), and the two are blended over join +- width with a
    quintic smoothstep instead of meeting with a kink. overrides: {name: {key: value}} replaces constants.
    """
    mu = reduced_mass('127I2')
    curves = {}
    for name, c in perturber_constants().items():
        c = dict(c)
        if experimental:
            c.update({k: v for k, v in c["exp"].items() if v is not None})
        c.update((overrides or {}).get(name, {}))
        a = c["we"]/(2*sqrt(c["De"]*HBAR2_2U/mu))

        def inner(R, c=c, a=a):
            return c["De"]*((1-np.exp(-a*(R-c["re"])))**2-1)

        def outer(R, c=c):
            return c["C5"]/R**5+c["C6"]/R**6

        shift = float(outer(np.array(join))-inner(np.array(join)))

        def V(R, inner=inner, outer=outer, shift=shift):
            R = np.asarray(R, dtype=float)
            s = np.clip((R-join+width)/(2*width), 0., 1.)
            s = s**3*(10-15*s+6*s*s)
            Rs = np.maximum(R, .5)
            return limit+(1-s)*(inner(R)+shift)+s*outer(Rs)
        curves[name] = V
    return curves


class PerturberPotential:
    """A perturbing state's potential (cm-1 above the X minimum, R in A), as BSplineSolver wants it."""

    def __init__(self, function):
        self.function = function

    def __call__(self, R):
        return self.function(R)

    def adiabatic(self, R, mass_ratio=1.):
        return np.zeros_like(R)

    def nonadiabatic(self, R, mass_ratio=1.):
        return np.zeros_like(R)


class Radial:
    """B levels and the resolvents of the perturbing states, on one B-spline basis.

    green(name, v, J, J') = {(k, k'): <v J|O_k (E_vJ - H_p(J'))^-1 O_k'|v J>} in MHz^-1, O_0 = (R_TABLE/R)^2 (the
    radial factor of V0), O_1 = O_2 = 1, H_p(J') with centrifugal term J'(J'+1) - Omega'^2. Solved as one banded
    linear system: every bound level of p and its continuum discretized in the box, with no sum over states.
    """

    def __init__(self, b_potential, joins, perturbers, grid=GRID, nlev=90):
        mu = reduced_mass('127I2')
        self.order = grid['order']
        self.b = BSplineSolver(b_potential, mu, joins=joins, nlev=nlev, **grid)
        self.p = {name: BSplineSolver(PerturberPotential(f), mu, joins=joins, nlev=1, **grid) for name, f in perturbers.items()}
        R = self.b.R
        self.ops = np.column_stack([(R_TABLE/R)**2, np.ones_like(R)])
        self.levels = {}

    def level(self, J):
        if J not in self.levels:
            self.levels[J] = self.b.wavefunctions(J)
        return self.levels[J]

    def _banded(self, A):
        w = self.order-1
        n = len(A)
        ab = np.zeros((2*w+1, n))
        for d in range(-w, w+1):
            diag = np.diagonal(A, d)
            if d >= 0:
                ab[w-d, d:] = diag
            else:
                ab[w-d, :n+d] = diag
        return ab

    def green(self, name, v, J, Jp, energy=None):
        E, psi = self.level(J)
        E = E[v] if energy is None else energy
        s = self.p[name]
        omega = OMEGA[name]
        A = E*s.S-(s._T+s._V+(Jp*(Jp+1)-omega**2)*s._rot)
        b = s._B.T@(s.W[:, None]*self.ops*psi[:, v:v+1])          # (basis, 2): O_0 chi, O_1 chi
        x = solve_banded((self.order-1, self.order-1), self._banded(A), b)
        g = b.T@x/MHZ_PER_CM
        return {(k, kk): g[min(k, 1), min(kk, 1)] for k in (0, 1, 2) for kk in (0, 1, 2)}


def f_factor(k, omega, eps, dJ, value):
    """C05 eqs 5 and 9 without the vibrational overlap: (-1)^W' C_k<V^k> / 3j, times the Omega-doubling factor."""
    f = _sign(omega)*value/w3j(k, IA, IA, 0, IA, -IA)
    if omega:
        f *= (1 + eps*_sign(dJ))/SQRT2
    return f


def element(k, Jp, Ip, J, I, F, omega, parity, eps, value):
    """<p W' eps' J' I' F|H^k(a) + H^k(b)|B 0u+ J I F> per unit vibrational overlap, MHz (C05 eqs 4, 9, 10)."""
    dJ, dI = Jp-J, Ip-I
    s_b = _sign(dI) if parity == "u" else -_sign(dI)
    if 1+s_b == 0.:
        return 0.
    angular = (_sign(F+dI+2*IA+dJ+k)*sqrt((2*Jp+1)*(2*J+1)*(2*Ip+1)*(2*I+1))*w3j(Jp, k, J, -omega, omega, 0)
               * w6j(Jp, k, J, I, F, Ip)*w6j(k, IA, IA, IA, Ip, I))
    return angular*f_factor(k, omega, eps, dJ, value)*(1+s_b)


def coriolis(J, omega, eps, value):
    """<1u eps'=+1 J I F|V0|B J I F> per unit radial integral <v'|(R_TABLE/R)^2|v>, MHz (C05 eq 12, P86)."""
    if omega != 1 or eps != 1 or not value:
        return 0.
    return -sqrt(2*J*(J+1))*value


def perturber_channels(name, J):
    """The (J', I'-independent) channels of state `name` reachable from B J: (J', eps') pairs."""
    omega = OMEGA[name]
    eps_values = (0,) if omega == 0 else (1, -1)
    return [(Jp, eps) for Jp in range(max(J-2, omega), J+3) for eps in eps_values]


def second_order_block(J, F, name, G, i1=IA, i2=IA, couplings=None):
    """H2 inside B (J, F) from state `name`: matrix over B's allowed I, MHz.

    G(k, k', J') -> the radial resolvent element (MHz^-1 per unit electronic coupling squared), with k = 0 the
    Coriolis operator (radial factor (R_TABLE/R)^2) and k = 1, 2 a unit operator. couplings: COUPLINGS by default,
    or e.g. molecular_couplings().
    """
    omega, parity, v0, v1, v2 = (couplings or COUPLINGS)[name]
    spins = [I for I in allowed_spins(J, i1, i2, "u") if abs(J-I) <= F <= J+I]
    H = np.zeros((len(spins), len(spins)))
    values = {0: v0, 1: v1, 2: v2}
    for Jp, eps in perturber_channels(name, J):
        for Ip in range(0, int(2*IA)+1):
            if not abs(Jp-Ip) <= F <= Jp+Ip:
                continue
            row = {}
            for k in (0, 1, 2):
                if not values[k]:
                    continue
                if k == 0:
                    e = [coriolis(J, omega, eps, v0) if (Jp == J and Ip == I) else 0. for I in spins]
                else:
                    e = [element(k, Jp, Ip, J, I, F, omega, parity, eps, values[k]) for I in spins]
                row[k] = np.array(e)
            for k, a in row.items():
                for kk, b in row.items():
                    if k == 0 and kk == 0:
                        continue           # the V0 x V0 term only shifts the rotational level
                    H += np.outer(a, b)*G(k, kk, Jp)
    return spins, H


def second_order_matrix(basis, F, v, name, radial, couplings=None):
    """H2 of state `name` over the B basis [(J_b, I)] of one F block, MHz.

    Van Vleck's symmetrized form for the rotational neighbours J_b = J, J +- 2 that share the block:
    <a|H2|b> = sum_p <a|H|p><p|H|b> (G_a + G_b)/2, G_a the resolvent at B (v, J_a)'s energy.
    """
    omega, parity, v0, v1, v2 = (couplings or COUPLINGS)[name]
    values = {0: v0, 1: v1, 2: v2}
    n = len(basis)
    H = np.zeros((n, n))
    Jb = sorted({J for J, _ in basis})
    for Jp in range(max(min(Jb)-2, omega), max(Jb)+3):
        G = {J: radial.green(name, v, J, Jp) for J in Jb if abs(J-Jp) <= 2}
        if not G:
            continue
        for eps in ((0,) if omega == 0 else (1, -1)):
            for Ip in range(0, int(2*IA)+1):
                if not abs(Jp-Ip) <= F <= Jp+Ip:
                    continue
                e = {}
                for k in (0, 1, 2):
                    if not values[k]:
                        continue
                    if k == 0:
                        e[k] = np.array([coriolis(J, omega, eps, v0) if (J == Jp and I == Ip) else 0.
                                         for J, I in basis])
                    else:
                        e[k] = np.array([element(k, Jp, Ip, J, I, F, omega, parity, eps, values[k])
                                         if J in G else 0. for J, I in basis])
                if not any(np.any(x) for x in e.values()):
                    continue
                g = {(k, kk): np.array([G[J][k, kk] if J in G else 0. for J, _ in basis]) for k in e for kk in e}
                for k in e:
                    for kk in e:
                        if k == 0 and kk == 0:
                            continue          # V0 x V0: a rotational shift, part of B's own rotational function
                        a = e[k][:, None]*e[kk][None, :]
                        H += a*(g[k, kk][:, None]+g[k, kk][None, :])/2
    return H


def structure(J, v, params, rot_energy, radial, names, dJ=2, couplings=None):
    """Hyperfine levels of B (v, J) with H2 of the states `names` added to the effective Hamiltonian of params
    (i2spec.hyperfine.level_structure, symmetry u, J +- dJ neighbours)."""
    J_values = [Jn for Jn in range(J-dJ, J+dJ+1, 2) if Jn >= 0]
    E0 = rot_energy(J)
    blocks = build_blocks(J_values, params, lambda Jn: rot_energy(Jn)-E0, IA, IA, "u")
    levels = []
    for F, (basis, H) in blocks.items():
        H = H.copy()
        for name in names:
            H += second_order_matrix(basis, F, v, name, radial, couplings)
        energies, vectors = np.linalg.eigh(H)
        central = np.array([Jb == J for Jb, _ in basis])
        weights = np.where(central[:, None], vectors**2, 0.)
        for k in np.flatnonzero(weights.sum(axis=0) > .5):
            levels.append(HyperfineLevel(F, basis[int(np.argmax(weights[:, k]))][1], float(energies[k]), basis,
                                         vectors[:, k]))
    return HyperfineLevels(sorted(levels, key=lambda level: level.energy))


def project(J, blocks, i1=IA, i2=IA):
    """Least-squares C (kHz), delta (kHz), eqQ (MHz), d (kHz) and a scalar (MHz) of the Delta J = 0 blocks.

    blocks: {F: (spins, H)}. Returns (constants, rms misfit MHz, rms of H MHz).
    """
    rows, rhs = [], []
    for F, (spins, H) in blocks.items():
        for a, Ia in enumerate(spins):
            for b, Ib in enumerate(spins):
                c = _coefficient_terms(J, Ia, J, Ib, F, i1, i2, 1., 1.)
                rows.append([*c, 1. if a == b else 0.])
                rhs.append(H[a, b])
    A, y = np.array(rows), np.array(rhs)
    x, *_ = np.linalg.lstsq(A, y, rcond=None)
    return dict(C=x[0], delta=x[1], eqQ=x[2], d=x[3], scalar=x[4]), float(np.sqrt(np.mean((A@x-y)**2))), \
        float(np.sqrt(np.mean(y**2)))


def second_order_constants(radial, v, J, names, couplings=None):
    """{name: (constants, misfit)} of B (v, J): each state's H2 inside J (Delta J = 0 blocks) projected onto the
    four effective constants (project(): C, delta, d in kHz, eqQ in MHz), with its J'-dependent resolvent."""
    out = {}
    for name in names:
        G = {Jp: radial.green(name, v, J, Jp) for Jp in range(max(J-2, OMEGA[name]), J+3)}
        blocks = {F: second_order_block(J, F, name, lambda k, kk, Jp: G[Jp][k, kk], couplings=couplings)
                  for F in range(max(0, J-int(2*IA)), J+int(2*IA)+1)}
        fit, misfit, _ = project(J, {F: b for F, b in blocks.items() if len(b[0])})
        out[name] = (fit, misfit)
    return out


def chen_constants(name, G, couplings=None):
    """C05 eqs 13-16, closed forms, for one state with J'-independent resolvent elements G(k, k')."""
    omega, parity, v0, v1, v2 = (couplings or COUPLINGS)[name]
    chi = 1. if parity == "g" else -1.
    dbl = 2. if omega else 1.
    norm = IA*(IA+1)*(2*IA+1)
    f = {k: _sign(omega)*v/w3j(k, IA, IA, 0, IA, -IA) for k, v in ((1, v1), (2, v2))}
    # f_k(a, 0u+, W') = f_k(a, W', 0u+) for real elements (hermiticity of the spherical components)
    f0 = v0/SQRT2
    out = dict(C=0., delta=0., d=0., eqQ=0.)
    if omega == 1 and parity == "u" and v0:
        out["C"] = -4/sqrt(norm)*f0*f[1]*G(0, 1)*1e3
        out["eqQ"] += -4*SQRT2*4*f0*f[2]*G(0, 2)
    out["delta"] = -2/(3*norm)*chi*dbl*f[1]**2*G(1, 1)*1e3
    out["d"] = -1/(3*norm)*chi*(3*omega**2-2)*dbl*f[1]**2*G(1, 1)*1e3
    for k in (1, 2):
        for kk in (1, 2):
            if not f[k] or not f[kk]:
                continue
            out["eqQ"] += (10*_sign(2*IA+omega)*dbl*2*f[k]*f[kk]*G(k, kk)
                           * w3j(k, kk, 2, omega, -omega, 0)*w6j(IA, 2, IA, k, IA, kk)*w3j(2, IA, IA, 0, IA, -IA))
    return out
