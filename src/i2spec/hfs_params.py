"""I₂ hyperfine parameters from published interpolation formulae.

BKT02: Bodermann, Knöckel & Tiemann, Eur. Phys. J. D 19, 31 (2002), eqs. (10)-(15), for ¹²⁷I₂.
S06: Salumbides et al., Mol. Phys. 104, 2641 (2006), Tables 2-4: a refit to all three isotopologues
that extends the B-state formulae to v' <= 53. The other isotopologues use the ¹²⁷I₂ formulae times
nuclear-moment ratios (S06 eq. 8, Table 1), evaluated with ¹²⁷I₂ energies at the same v.

E_x is the energy of X(v'', J''=0) above X(0, 0); E_b is the energy of B(v', J'=0) above X(0, 0).
Both are in cm⁻¹ and are taken from the ¹²⁷I₂ rovibronic model. eqQ is returned in MHz; C, d and
delta in kHz. BKT02 stated validity: v' <= 43 (E_b < 19 500 cm⁻¹), v'' <= 17, J'' <= 200. Stated 2σ
uncertainties (§7): ΔeqQ ±50 kHz, ΔC ±2 kHz, Δδ ±3 kHz, Δd ±4 kHz.

Deviations from the printed texts (docs/research/isotopologue-hyperfine.md):

BKT02 eqs. (12) and (13) are used as the B-state values δ_B and d_B. The paper says to add δ_X and
d_X to them. Evaluated at the levels of the paper's own Table 4 (17 lines, E from our model), eq.
(12) - Table 4 Δδ has mean +3.8 kHz as printed but +0.1 kHz (rms 1.3 kHz) after subtracting δ_X.
Likewise, eq. (13) - Table 4 Δd has mean +1.7 kHz but +0.15 kHz (rms 2.3 kHz) after subtracting
d_X. Only the second reading agrees with the stated 2 kHz fit scatter; S06 (p. 2649) confirms it.

S06 eqs. (6)-(7) are used with the sign of the Ω = 0 perturber term flipped. As printed, δ_B(v' = 32)
is -78 kHz against -6 kHz from BKT02, and the whole S06 set misses the BIPM R(56)32-0 intervals by
245 kHz rms; flipped, by 57 kHz. BKT02 gives 24 kHz there (S06 C_B accounts for the difference), so
line_states keeps BKT02 for ¹²⁷I₂ at v' <= 43.
"""

from dataclasses import dataclass
from math import exp

from .constants import ISOTOPOLOGUES
from .hyperfine import HyperfineParameters


@dataclass(frozen=True)
class Corrections:
    """Additive corrections to the published hyperfine formulae (M5 stage 2).

    Zero reproduces BKT02 and S06 exactly, so nothing changes unless a fit asks for it. Units follow
    HyperfineParameters: eqQ in MHz; C, d and delta in kHz. Each correction is a low-order polynomial
    in s = (v + 1/2)/10 and y = J(J + 1)/10^4, which are both of order one over the data, so every
    coefficient carries the units of its parameter. Suffixes name the monomial: ``C_b_sy`` multiplies
    s*y. The published formulae carry their own v and J dependence in the same variables.

    BKT02 section 7 gives 2 sigma uncertainties of 50 kHz on eqQ, 2 kHz on C, 3 kHz on delta and
    4 kHz on d; halved, those are the natural priors on the constant terms.
    """

    eqQ_b: float = 0.0
    eqQ_b_s: float = 0.0
    eqQ_b_y: float = 0.0
    C_b: float = 0.0
    C_b_s: float = 0.0
    C_b_s2: float = 0.0
    C_b_y: float = 0.0
    C_b_sy: float = 0.0
    C_b_y2: float = 0.0
    d_b: float = 0.0
    d_b_s: float = 0.0
    d_b_y: float = 0.0
    delta_b: float = 0.0
    delta_b_s: float = 0.0
    delta_b_y: float = 0.0
    eqQ_x: float = 0.0
    eqQ_x_s: float = 0.0
    C_x: float = 0.0
    C_x_s: float = 0.0

    def apply_b(self, p: HyperfineParameters, v, J) -> HyperfineParameters:
        s, y = (v + 0.5) / 10, J * (J + 1) / 1e4
        return HyperfineParameters(
            p.eqQ + self.eqQ_b + self.eqQ_b_s * s + self.eqQ_b_y * y,
            p.C + self.C_b + self.C_b_s * s + self.C_b_s2 * s * s + self.C_b_y * y + self.C_b_sy * s * y
            + self.C_b_y2 * y * y,
            p.d + self.d_b + self.d_b_s * s + self.d_b_y * y,
            p.delta + self.delta_b + self.delta_b_s * s + self.delta_b_y * y)

    def apply_x(self, p: HyperfineParameters, v, J) -> HyperfineParameters:
        s = (v + 0.5) / 10
        return HyperfineParameters(p.eqQ + self.eqQ_x + self.eqQ_x_s * s, p.C + self.C_x + self.C_x_s * s,
                                   p.d, p.delta)


#: The published formulae, unmodified.
NO_CORRECTIONS = Corrections()

#: Upper end of the validated range of the BKT02 B-state formulae (v' = 43 is at 19 434 cm⁻¹).
E_B_MAX = 19450.0

#: X-state spin-spin parameters, kHz, used for all v'' <= 17 (paper §6.2, from Wallerand et al. 1999).
DELTA_X = 3.705
D_X = 1.524


def eqq_x(v, J):
    u, x = v + 0.5, J * (J + 1)
    return -2452.2916 - 0.542 * u + 0.4534e-1 * u**2 - 0.1927e-3 * x + 0.694e-5 * u * x


def eqq_b(v, J):
    u, x = v + 0.5, J * (J + 1)
    return -487.879 - 1.8621 * u + 0.12511e-3 * u**3 - 0.1281e-3 * x - 0.225e-5 * u * x - 0.308e-9 * x**2


def c_x(v, E_x):
    return 1.9245 + 0.01356 * (v + 0.5) - 15098 / (E_x - 12340)


def c_b(v, J, E_b):
    u, x = v + 0.5, J * (J + 1)
    return -4.016 - 0.1501 * u - 3.957e-4 * x - 1.767e-5 * u * x - (110704 + 1.862 * x) / (E_b - 19986)


def _bump(E_b):
    return 22.41 * exp(-((E_b - 16787) ** 2) / 260267)


def delta_b(E_b):
    """δ_B, eq. (12) (see the module docstring)."""
    return 31.57 - 32452 / (E_b - 19896) + 126257 / (E_b - 20687) - _bump(E_b)


def d_b(E_b):
    """d_B, eq. (13) (see the module docstring)."""
    return 19.56 + 32452 / (E_b - 19896) + 0.5 * (126257 / (E_b - 20687) - _bump(E_b))


def x_state(v, J, E_x):
    return HyperfineParameters(eqQ=eqq_x(v, J), C=c_x(v, E_x), d=D_X, delta=DELTA_X)


def b_state(v, J, E_b):
    return HyperfineParameters(eqQ=eqq_b(v, J), C=c_b(v, J, E_b), d=d_b(E_b), delta=delta_b(E_b))


#: Upper end of the validated range of the S06 B-state formulae (v' = 53 is at 19 772 cm⁻¹).
#: Above it the formulae are held at this energy. That keeps them finite, but it is wrong: Chen 2004 measures
#: C_B rising to 2.2 MHz at v' = 70 while the held value stays at 500-650 kHz, and the hyperfine intervals
#: of v' = 55-70 lines come out 63 MHz rms (median) away from those Chen's parameters give, up to 470 MHz
#: (docs/design/hyperfine-fit.md, "Chen 2004").
E_B_MAX_S06 = 19785.0

#: Nuclear quadrupole-moment and magnetic-moment (μ/I) ratios to ¹²⁷I (S06 Table 1).
GAMMA_Q = {127: 1.0, 129: 0.701213}
GAMMA_MU = {127: 1.0, 129: 0.6655}


def eqq_x_s06(v, J):
    u, x = v + 0.5, J * (J + 1)
    return -2452.285 - 0.5474 * u + 0.4487e-1 * u**2 - 0.2089e-3 * x + 0.6965e-5 * u * x


def eqq_b_s06(v, J):
    u, x = v + 0.5, J * (J + 1)
    return (-488.086 - 1.83777 * u + 0.99774e-4 * u**3 + 0.80818e-8 * u**5
            - 0.14020e-3 * x - 0.32217e-5 * u * x + 0.2716e-7 * u**2 * x - 0.3377e-9 * x**2)


def c_b_s06(v, J, E_b):
    u, x = v + 0.5, J * (J + 1)
    return (28.89 + 0.9234e-1 * u + 0.2241e-2 * x + 0.1027e-4 * u * x
            + (36672 - 3388 * u + 9.7985 * x - 0.2201 * u * x) / (E_b - 20140))


def _perturbers_s06(v, J, E_b):
    """Ω = 0 and Ω = 1 perturber terms of S06 eqs. (6)-(7); the Ω = 0 sign is flipped (module docstring)."""
    u2x = (v + 0.5) ** 2 * J * (J + 1)
    return -(39698 + 0.2475e-3 * u2x) / (E_b - 19973), (106923 + 0.1712e-3 * u2x) / (E_b - 20534)


def delta_b_s06(v, J, E_b):
    omega0, omega1 = _perturbers_s06(v, J, E_b)
    return 20.88 + omega0 + omega1 - _bump(E_b)


def d_b_s06(v, J, E_b):
    omega0, omega1 = _perturbers_s06(v, J, E_b)
    return 25.93 - omega0 + 0.5 * (omega1 - _bump(E_b))


def _scales(isotopologue):
    """S06 eq. (8): factors for eqQ and C of nucleus 1, and for the spin-spin terms (omitted for 127I129I, as in S06)."""
    a, b = ISOTOPOLOGUES[isotopologue]
    return GAMMA_Q[a], GAMMA_MU[a], GAMMA_MU[a] * GAMMA_MU[b] if a == b else 0.0


def nucleus_ratios(isotopologue):
    """(eqQ_ratio, C_ratio): eqQ and C of nucleus 2 over those of nucleus 1."""
    a, b = ISOTOPOLOGUES[isotopologue]
    return GAMMA_Q[b] / GAMMA_Q[a], GAMMA_MU[b] / GAMMA_MU[a]


def x_state_s06(v, J, E_x, isotopologue="127I2"):
    q, m, s = _scales(isotopologue)
    return HyperfineParameters(eqQ=q * eqq_x_s06(v, J), C=m * c_x(v, E_x), d=s * D_X, delta=s * DELTA_X)


def b_state_s06(v, J, E_b, isotopologue="127I2"):
    q, m, s = _scales(isotopologue)
    return HyperfineParameters(eqQ=q * eqq_b_s06(v, J), C=m * c_b_s06(v, J, E_b), d=s * d_b_s06(v, J, E_b),
                               delta=s * delta_b_s06(v, J, E_b))


def line_states(isotopologue, v_upper, v_lower, E_b, E_x, corrections=None, table=None):
    """Parameter functions J -> HyperfineParameters of the (lower, upper) levels of one B-X line.

    BKT02 for ¹²⁷I₂ up to E_B_MAX (v' <= 43); S06 otherwise, with E_b frozen at E_B_MAX_S06 above
    v' = 53, because the formulae have poles near the dissociation limit. ``corrections`` adds the
    fitted offsets of ``Corrections``; ``table`` (an hfs_table.HyperfineTable) adds the measured
    corrections of the B state for ¹²⁷I₂. None for both leaves the published formulae untouched.
    """
    c = corrections or NO_CORRECTIONS
    tab = (lambda p, J: table.apply(p, v_upper, J)) if (table is not None and isotopologue == "127I2") \
        else (lambda p, J: p)                                                       # noqa: E731
    if isotopologue == "127I2" and E_b <= E_B_MAX:
        lower = lambda J: c.apply_x(x_state(v_lower, J, E_x), v_lower, J)          # noqa: E731
        upper = lambda J: tab(c.apply_b(b_state(v_upper, J, E_b), v_upper, J), J)  # noqa: E731
        return lower, upper
    E_b = min(E_b, E_B_MAX_S06)
    lower = lambda J: c.apply_x(x_state_s06(v_lower, J, E_x, isotopologue), v_lower, J)   # noqa: E731
    upper = lambda J: tab(c.apply_b(b_state_s06(v_upper, J, E_b, isotopologue), v_upper, J), J)   # noqa: E731
    return lower, upper
