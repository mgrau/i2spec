"""Hyperfine structure of I₂ rovibrational levels and of B-X line components.

Effective Hamiltonian (Broyer, Vigué & Lehmann, J. Phys. 39, 591 (1978)):

    H = eqQ·H_EQ + C·H_SR + d·H_TSS + δ·H_SSS,

set up in the coupled basis |J, (i1 i2) I; F⟩ and diagonalized for each F. Couplings with ΔJ = ±2
and ΔI = ±2 are included, as in the fits of Bodermann et al., EPJD 19, 31 (2002), §3. When the two
nuclei differ (¹²⁷I¹²⁹I), eqQ and C differ between them and ΔI = ±1 couplings appear as well.
eqQ_ratio and C_ratio are the nucleus-2 constants relative to nucleus 1 (Salumbides et al., Mol.
Phys. 104, 2641 (2006)).

Within a ΔJ = 0 block the operators reduce to:

    H_EQ  = −Σ_n [3(I_n·J)² + (3/2)(I_n·J) − I_n²J²] / [2 i_n (2 i_n − 1)(2J − 1)(2J + 3)]

(the Townes-Schawlow sign for an Ω = 0 state; the opposite sign mirrors the R(56)32-0 pattern and
misses the BIPM intervals by ~100 MHz)
    H_SR  = I1·J + C_ratio·I2·J
    H_TSS = [3(I1·J)(I2·J) + 3(I2·J)(I1·J) − 2 (I1·I2) J(J+1)] / [(2J − 1)(2J + 3)]

(the ΔJ = 0 projection of the dipolar form I1·I2 − 3(I1·n̂)(I2·n̂); a factor of 5 in front, as
sometimes quoted, misses the BIPM R(56)32-0 intervals by ~0.7 MHz instead of ~25 kHz)
    H_SSS = I1·I2

The off-diagonal blocks come from the same tensor operators:
K_Q · Σ_n √6/(4 i_n(2 i_n − 1)) [I_n⊗I_n]²·C²(n̂) and K_T · [I1⊗I2]²·C²(n̂).
tests/test_hyperfine.py checks both against brute-force constructions.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import sqrt

import numpy as np

from .wigner import wigner_3j as w3j
from .wigner import wigner_6j as w6j
from .wigner import wigner_9j as w9j

#: Factors that make the ΔJ = 0 blocks equal the forms in the module docstring.
K_Q = 1.0
K_T = -sqrt(6.0)


@dataclass(frozen=True)
class HyperfineParameters:
    """Effective hyperfine parameters of one rovibrational level: eqQ in MHz; C, d, delta in kHz."""

    eqQ: float
    C: float
    d: float
    delta: float


def _sign(x):
    """(-1)**x for integer-valued x."""
    return -1.0 if round(x) % 2 else 1.0


def rot_reduced(Jp, J, k):
    """<J'||C^k(n)||J> between Ω = 0 rotational states."""
    return _sign(Jp) * sqrt((2 * J + 1) * (2 * Jp + 1)) * w3j(Jp, k, J, 0, 0, 0)


def _vec_reduced(j):
    return sqrt(j * (j + 1) * (2 * j + 1))


def _quad_reduced(i):
    """<i||[I⊗I]^(2)||i> (Edmonds 7.1.1)."""
    return _sign(2 * i + 2) * sqrt(5) * w6j(1, 1, 2, i, i, i) * _vec_reduced(i) ** 2


def _one_nucleus(n, Ip, I, i1, i2, k, red):
    """<(i1 i2) I'||T^k(n)||(i1 i2) I> for an operator acting on nucleus n (Edmonds 7.1.7, 7.1.8)."""
    norm = sqrt((2 * I + 1) * (2 * Ip + 1)) * red
    if n == 1:
        return _sign(i1 + i2 + I + k) * norm * w6j(i1, Ip, i2, I, i1, k)
    return _sign(i1 + i2 + Ip + k) * norm * w6j(i2, Ip, i1, I, i2, k)


def _spin_spin_reduced(Ip, I, i1, i2):
    """<(i1 i2) I'||[I1⊗I2]^(2)||(i1 i2) I> (Edmonds 7.1.5)."""
    return sqrt(5 * (2 * I + 1) * (2 * Ip + 1)) * w9j(i1, i1, 1, i2, i2, 1, Ip, I, 2) * _vec_reduced(i1) * _vec_reduced(i2)


def matrix_element(Jp, Ip, J, I, F, p: HyperfineParameters, i1, i2, eqQ_ratio=1.0, C_ratio=1.0):
    """<J' I' F|H|J I F> in MHz. eqQ_ratio, C_ratio = eqQ, C of nucleus 2 over those of nucleus 1."""
    value = 0.0
    if Jp == J and Ip == I:
        value += 1e-3 * p.C * (F * (F + 1) - I * (I + 1) - J * (J + 1)) / 2
        value += 1e-3 * p.delta * (I * (I + 1) - i1 * (i1 + 1) - i2 * (i2 + 1)) / 2
    if Jp == J and C_ratio != 1.0:  # C I·J + (C_ratio - 1) C I2·J, Edmonds 7.1.6
        value += (1e-3 * p.C * (C_ratio - 1.0) * _sign(J + Ip + F) * w6j(F, Ip, J, 1, J, I) * _vec_reduced(J)
                  * _one_nucleus(2, Ip, I, i1, i2, 1, _vec_reduced(i2)))
    rot = rot_reduced(Jp, J, 2)
    if rot == 0.0:
        return value
    geometry = _sign(J + Ip + F) * w6j(F, Ip, Jp, 2, J, I) * rot
    quad = sum(ratio * sqrt(6) / (4 * i * (2 * i - 1)) * _one_nucleus(n, Ip, I, i1, i2, 2, _quad_reduced(i))
               for n, (i, ratio) in enumerate(((i1, 1.0), (i2, eqQ_ratio)), start=1))
    return value + geometry * (K_Q * p.eqQ * quad + K_T * 1e-3 * p.d * _spin_spin_reduced(Ip, I, i1, i2))


def allowed_spins(J, i1, i2, symmetry):
    """Total nuclear spins I that combine with rotational level J.

    symmetry: "g" or "u" for identical fermionic nuclei in a 0g+ or 0u+ state (I ≡ J or I ≢ J mod 2);
    None for distinguishable nuclei (no restriction).
    """
    spins = range(round(abs(i1 - i2)), round(i1 + i2) + 1)
    if symmetry is None:
        return list(spins)
    if i1 != i2 or round(2 * i1) % 2 == 0:
        raise ValueError("exchange symmetry is implemented for identical fermionic nuclei only")
    return [I for I in spins if ((I - J) % 2 == 0) == (symmetry == "g")]


def _mean(a: HyperfineParameters, b: HyperfineParameters) -> HyperfineParameters:
    return HyperfineParameters((a.eqQ + b.eqQ) / 2, (a.C + b.C) / 2, (a.d + b.d) / 2, (a.delta + b.delta) / 2)


def build_blocks(J_values, params, rot_energy, i1=2.5, i2=2.5, symmetry="g", eqQ_ratio=1.0, C_ratio=1.0):
    """Hamiltonian blocks {F: (basis, H)} spanning the rotational levels J_values.

    params(J) -> HyperfineParameters; rot_energy(J) -> rotational energy in MHz.
    """
    states = [(J, I) for J in J_values for I in allowed_spins(J, i1, i2, symmetry)]
    blocks = {}
    for F in sorted({F for J, I in states for F in range(abs(J - I), J + I + 1)}):
        basis = [(J, I) for J, I in states if abs(J - I) <= F <= J + I]
        H = np.array([[matrix_element(Ja, Ia, Jb, Ib, F, params(Ja) if Ja == Jb else _mean(params(Ja), params(Jb)),
                                      i1, i2, eqQ_ratio, C_ratio) for Jb, Ib in basis] for Ja, Ia in basis])
        if not np.allclose(H, H.T, atol=1e-9):
            raise RuntimeError(f"hyperfine block F={F} is not symmetric")
        H += np.diag([rot_energy(J) for J, _ in basis])
        blocks[F] = (basis, H)
    return blocks


@dataclass
class HyperfineLevel:
    F: int
    I: int  # dominant total nuclear spin
    energy: float  # MHz, relative to the hyperfine-free rovibrational level
    basis: list
    vector: np.ndarray


def level_structure(J, params, rot_energy, *, i1=2.5, i2=2.5, symmetry="g", eqQ_ratio=1.0, C_ratio=1.0, dJ=2):
    """Hyperfine levels of rotational level J, with rotational neighbours J ± 2, ..., J ± dJ in the basis.

    rot_energy(J) is in MHz; only differences to rot_energy(J) matter.
    """
    E0 = rot_energy(J)
    J_values = [Jn for Jn in range(J - dJ, J + dJ + 1, 2) if Jn >= 0]
    levels = []
    blocks = build_blocks(J_values, params, lambda Jn: rot_energy(Jn) - E0, i1, i2, symmetry, eqQ_ratio, C_ratio)
    for F, (basis, H) in blocks.items():
        energies, vectors = np.linalg.eigh(H)
        central = np.array([Jb == J for Jb, _ in basis])
        for k, energy in enumerate(energies):
            weights = np.where(central, vectors[:, k] ** 2, 0.0)
            if weights.sum() > 0.5:
                levels.append(HyperfineLevel(F, basis[int(np.argmax(weights))][1], float(energy), basis, vectors[:, k]))
    return sorted(levels, key=lambda level: level.energy)


@dataclass
class Component:
    offset: float  # MHz, relative to the hyperfine-free rovibronic line
    strength: float  # relative line strength; sums to 1 over the line
    F_upper: int
    F_lower: int
    I_upper: int  # dominant total nuclear spin of the upper level
    I_lower: int
    label: str | None = None  # a1, a2, ... for main components (ΔF = ΔJ, same I), by increasing frequency


def dipole_reduced(Jp, Fp, J, F, I):
    """<J' I F'||C^1(n)||J I F> for an Ω = 0 - 0 transition; the nuclear spins are spectators."""
    return _sign(Jp + I + F + 1) * sqrt((2 * F + 1) * (2 * Fp + 1)) * w6j(Jp, Fp, I, F, J, 1) * rot_reduced(Jp, J, 1)


def line_components(upper, lower, J_upper, J_lower, threshold=1e-10):
    """Components of one rovibronic line from the hyperfine levels of its upper and lower levels.

    Main components (labels a1, a2, ... by increasing frequency) have ΔF = ΔJ. They are chosen by
    strongest-first matching, one per upper and lower level, because I is not a good quantum
    number at high J: the quadrupole coupling mixes the degenerate I states.
    """
    found = []
    for iu, u in enumerate(upper):
        for il, low in enumerate(lower):
            if abs(u.F - low.F) > 1:
                continue
            amplitude = 0.0
            for a, (Ja, Ia) in enumerate(u.basis):
                for b, (Jb, Ib) in enumerate(low.basis):
                    if Ia == Ib and abs(Ja - Jb) == 1:
                        amplitude += u.vector[a] * low.vector[b] * dipole_reduced(Ja, u.F, Jb, low.F, Ia)
            if amplitude**2 > threshold:
                found.append((iu, il, Component(u.energy - low.energy, amplitude**2, u.F, low.F, u.I, low.I)))
    total = sum(c.strength for _, _, c in found)
    for _, _, c in found:
        c.strength /= total
    used_upper, used_lower, main = set(), set(), []
    for iu, il, c in sorted(found, key=lambda t: -t[2].strength):
        if c.F_upper - c.F_lower == J_upper - J_lower and iu not in used_upper and il not in used_lower:
            used_upper.add(iu)
            used_lower.add(il)
            main.append(c)
    for n, c in enumerate(sorted(main, key=lambda c: c.offset), start=1):
        c.label = f"a{n}"
    return sorted((c for _, _, c in found), key=lambda c: c.offset)
