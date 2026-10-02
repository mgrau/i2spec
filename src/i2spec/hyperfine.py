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
from functools import lru_cache
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


def _coefficient_terms(Jp, Ip, J, I, F, i1, i2, eqQ_ratio, C_ratio):
    """matrix_element split by parameter: the factors multiplying C, delta, eqQ and d (MHz per unit)."""
    c = dl = q = t = 0.0
    if Jp == J and Ip == I:
        c += 1e-3 * (F * (F + 1) - I * (I + 1) - J * (J + 1)) / 2
        dl += 1e-3 * (I * (I + 1) - i1 * (i1 + 1) - i2 * (i2 + 1)) / 2
    if Jp == J and C_ratio != 1.0:
        c += (1e-3 * (C_ratio - 1.0) * _sign(J + Ip + F) * w6j(F, Ip, J, 1, J, I) * _vec_reduced(J)
              * _one_nucleus(2, Ip, I, i1, i2, 1, _vec_reduced(i2)))
    rot = rot_reduced(Jp, J, 2)
    if rot != 0.0:
        geometry = _sign(J + Ip + F) * w6j(F, Ip, Jp, 2, J, I) * rot
        quad = sum(ratio * sqrt(6) / (4 * i * (2 * i - 1)) * _one_nucleus(n, Ip, I, i1, i2, 2, _quad_reduced(i))
                   for n, (i, ratio) in enumerate(((i1, 1.0), (i2, eqQ_ratio)), start=1))
        q = geometry * K_Q * quad
        t = geometry * K_T * 1e-3 * _spin_spin_reduced(Ip, I, i1, i2)
    return c, dl, q, t


@lru_cache(maxsize=8192)
def _block_coefficients(J_values, i1, i2, symmetry, eqQ_ratio, C_ratio):
    """{F: (basis, J of each basis state, coefficient matrices of C, delta, eqQ, d)} for build_blocks.

    The Hamiltonian is linear in the four parameters, and the factors depend on the angular momenta
    alone, so they are computed once per set of rotational levels and reused for every vibrational level.
    """
    states = [(J, I) for J in J_values for I in allowed_spins(J, i1, i2, symmetry)]
    out = {}
    for F in sorted({F for J, I in states for F in range(abs(J - I), J + I + 1)}):
        basis = [(J, I) for J, I in states if abs(J - I) <= F <= J + I]
        terms = np.array([[_coefficient_terms(Ja, Ia, Jb, Ib, F, i1, i2, eqQ_ratio, C_ratio) for Jb, Ib in basis]
                          for Ja, Ia in basis])                     # (n, n, 4)
        coef = np.ascontiguousarray(np.moveaxis(terms, -1, 0))
        if not np.allclose(coef, np.swapaxes(coef, 1, 2), atol=1e-12):
            raise RuntimeError(f"hyperfine block F={F} is not symmetric")
        for a in coef:
            a.flags.writeable = False
        out[F] = (basis, np.array([J for J, _ in basis]), coef)
    return out


def build_blocks(J_values, params, rot_energy, i1=2.5, i2=2.5, symmetry="g", eqQ_ratio=1.0, C_ratio=1.0):
    """Hamiltonian blocks {F: (basis, H)} spanning the rotational levels J_values.

    params(J) -> HyperfineParameters; rot_energy(J) -> rotational energy in MHz. Between two rotational
    levels the parameters are the mean of the two levels' ones.
    """
    J_values = tuple(int(J) for J in J_values)
    coefficients = _block_coefficients(J_values, i1, i2, symmetry, float(eqQ_ratio), float(C_ratio))
    index = {J: k for k, J in enumerate(J_values)}
    p = np.array([[pj.C, pj.delta, pj.eqQ, pj.d] for pj in map(params, J_values)])        # (nJ, 4)
    pairs = (p[:, None, :] + p[None, :, :]) / 2                                           # the mean off the diagonal
    pairs[np.arange(len(J_values)), np.arange(len(J_values))] = p
    rot = np.array([rot_energy(J) for J in J_values], dtype=float)
    blocks = {}
    for F, (basis, Js, coef) in coefficients.items():
        k = np.array([index[J] for J in Js])
        P = pairs[k[:, None], k[None, :]]                                                  # (n, n, 4)
        H = np.einsum("knm,nmk->nm", coef, P)
        H += np.diag(rot[k])
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
    J_values = [Jn for Jn in range(J - dJ, J + dJ + 1, 2) if Jn >= 0]
    levels = []
    if len(J_values) == 1:          # J alone: its rotational energy is the zero, and need not be known
        relative = lambda Jn: 0.0   # noqa: E731
    else:
        E0 = rot_energy(J)
        relative = lambda Jn: rot_energy(Jn) - E0   # noqa: E731
    blocks = build_blocks(J_values, params, relative, i1, i2, symmetry, eqQ_ratio, C_ratio)
    for F, (basis, H) in blocks.items():
        energies, vectors = np.linalg.eigh(H)
        central = np.array([Jb == J for Jb, _ in basis])
        weights = np.where(central[:, None], vectors ** 2, 0.0)
        for k in np.flatnonzero(weights.sum(axis=0) > 0.5):
            levels.append(HyperfineLevel(F, basis[int(np.argmax(weights[:, k]))][1], float(energies[k]), basis,
                                         vectors[:, k]))
    return HyperfineLevels(sorted(levels, key=lambda level: level.energy))


class HyperfineLevels(list):
    """The hyperfine levels of one rotational level, by energy (a list), keeping the matrix of their
    eigenvectors that line_components uses, so that the lines sharing a level build it once."""

    _embedding = None


def _embedding(levels):
    """(blocks, U): the F blocks of ``levels`` as ((F, basis), ...) and the eigenvectors as the columns of U,
    each in the rows of its own block, so that every amplitude of a line is one matrix product."""
    cached = getattr(levels, "_embedding", None)
    if cached is not None:
        return cached
    order, start = {}, 0
    for level in levels:
        if level.F not in order:
            order[level.F] = (start, tuple(level.basis))
            start += len(level.basis)
    U = np.zeros((start, len(levels)))
    for k, level in enumerate(levels):
        a = order[level.F][0]
        U[a:a + len(level.basis), k] = level.vector
    result = (tuple((F, basis) for F, (_, basis) in order.items()), U)
    if isinstance(levels, HyperfineLevels):
        levels._embedding = result
    return result


@lru_cache(maxsize=16384)
def _dipole_operator(blocks_upper, blocks_lower):
    """dipole_reduced between every basis state of the upper and of the lower blocks (as from _embedding)."""
    sizes = lambda blocks: np.cumsum([0] + [len(b) for _, b in blocks])   # noqa: E731
    su, sl = sizes(blocks_upper), sizes(blocks_lower)
    D = np.zeros((su[-1], sl[-1]))
    for m, (Fu, bu) in enumerate(blocks_upper):
        for n, (Fl, bl) in enumerate(blocks_lower):
            if abs(Fu - Fl) <= 1:
                D[su[m]:su[m + 1], sl[n]:sl[n + 1]] = _dipole_matrix(bu, Fu, bl, Fl)
    D.flags.writeable = False
    return D


@dataclass
class Component:
    offset: float  # MHz, relative to the hyperfine-free rovibronic line
    strength: float  # relative line strength; sums to 1 over the line
    F_upper: int
    F_lower: int
    I_upper: int  # dominant total nuclear spin of the upper level
    I_lower: int
    label: str | None = None  # a1, a2, ... for main components (ΔF = ΔJ, same I), by increasing frequency
    upper_level: int | None = None  # which hyperfine eigenstate of the upper (lower) level: unlike (I, F),
    lower_level: int | None = None  # unique, since I is only approximately good at high J


def dipole_reduced(Jp, Fp, J, F, I):
    """<J' I F'||C^1(n)||J I F> for an Ω = 0 - 0 transition; the nuclear spins are spectators."""
    return _sign(Jp + I + F + 1) * sqrt((2 * F + 1) * (2 * Fp + 1)) * w6j(Jp, Fp, I, F, J, 1) * rot_reduced(Jp, J, 1)


@lru_cache(maxsize=65536)
def _dipole_matrix(basis_upper, F_upper, basis_lower, F_lower):
    """dipole_reduced between the coupled basis states of an upper and a lower F block (rows: upper)."""
    D = np.zeros((len(basis_upper), len(basis_lower)))
    for a, (Ja, Ia) in enumerate(basis_upper):
        for b, (Jb, Ib) in enumerate(basis_lower):
            if Ia == Ib and abs(Ja - Jb) == 1:
                D[a, b] = dipole_reduced(Ja, F_upper, Jb, F_lower, Ia)
    D.flags.writeable = False
    return D


def line_components(upper, lower, J_upper, J_lower, threshold=1e-10):
    """Components of one rovibronic line from the hyperfine levels of its upper and lower levels.

    Main components (labels a1, a2, ... by increasing frequency) have ΔF = ΔJ. They are chosen by
    strongest-first matching, one per upper and lower level, because I is not a good quantum
    number at high J: the quadrupole coupling mixes the degenerate I states.
    """
    blocks_u, U = _embedding(upper)
    blocks_l, L = _embedding(lower)
    amplitudes = U.T @ _dipole_operator(blocks_u, blocks_l) @ L     # zero unless |F' - F| <= 1
    squared = amplitudes ** 2
    iu, il = np.nonzero(squared > threshold)                          # the components, upper level by upper level
    if not iu.size:
        return []
    F_u, F_l = np.array([u.F for u in upper])[iu], np.array([low.F for low in lower])[il]
    offset = np.array([u.energy for u in upper])[iu] - np.array([low.energy for low in lower])[il]
    strength = squared[iu, il]
    strength = strength / np.cumsum(strength)[-1]                     # summed in order, as sum() would
    # main components: strongest first (a stable sort, so ties keep the order above), one per level
    used_upper, used_lower, main = set(), set(), []
    candidates = np.argsort(-strength, kind="stable")
    up, low = iu.tolist(), il.tolist()
    for k in candidates[(F_u - F_l)[candidates] == J_upper - J_lower].tolist():
        if up[k] not in used_upper and low[k] not in used_lower:
            used_upper.add(up[k])
            used_lower.add(low[k])
            main.append(k)
    labels = [None] * iu.size
    for n, k in enumerate(sorted(main, key=lambda k: offset[k]), start=1):
        labels[k] = f"a{n}"
    I_u, I_l = [u.I for u in upper], [low.I for low in lower]
    rows = zip(offset.tolist(), strength.tolist(), F_u.tolist(), F_l.tolist(), iu.tolist(), il.tolist(), labels)
    found = [Component(o, s, fu, fl, I_u[u], I_l[low], lab, upper_level=u, lower_level=low)
             for o, s, fu, fl, u, low, lab in rows]
    return [found[k] for k in np.argsort(offset, kind="stable").tolist()]
