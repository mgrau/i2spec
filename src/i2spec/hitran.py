"""Line lists in the HITRAN 2004 fixed-width format (160 characters per transition, ``.par``), with a
tabulated partition function, for HAPI, RADIS and other line-by-line codes.

    i2spec hitran 500 520 --out i2_500-520     # writes i2_500-520.par and i2_500-520_q_127I2.txt

I₂ is not a HITRAN molecule, so the molecule number is a user choice (``molecule``, default 0) and the
isotopologue numbers are local (ISOTOPOLOGUE_IDS). Conventions, as in HITRAN:

* S at T_ref = 296 K, in cm⁻¹/(molecule cm⁻²), and E″ above X(v=0, J=0).
* The partition function includes the nuclear-spin degeneracy, so g″ = (2J″+1) g_ns(J″) and
  g′ = (2J′+1) g_ns(J″) (the nuclear spin is not changed by the transition); with hyperfine components,
  g = 2F+1 for each hyperfine level, which sums to the same partition function.
* The isotopic abundance is 1 for each isotopologue: S is per molecule of that isotopologue. ¹²⁷I is all
  of natural iodine; ¹²⁹I is not natural, so its lines must be scaled by the user's own abundance.
* Pressure broadening and shifts are not part of the model; γ_air, γ_self, n_air and δ_air are written as
  given (default zero: Doppler profiles, right for a pure-iodine cell at its vapour pressure).
* The position uncertainty code (Ierr, first digit) is HITRAN's decade code of this model's 1σ; the
  intensity code is 2 ("average or estimate").
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from .constants import MHZ_PER_CM

#: Local isotopologue numbers, in order of natural abundance as HITRAN numbers them
ISOTOPOLOGUE_IDS = {"127I2": 1, "127I129I": 2, "129I2": 3}
#: Isotopologue masses (u), for the Doppler width a line-by-line code computes
MASSES = {"127I2": 2 * 126.904473, "127I129I": 126.904473 + 128.904984, "129I2": 2 * 128.904984}
T_REF = 296.0
C_CM = 2.99792458e10     # cm/s
C2 = 1.4387768775        # cm K


@dataclass
class Record:
    """One transition: a rovibronic line, or one hyperfine component of it."""

    iso: str
    nu: float            # cm-1
    S: float             # cm-1/(molecule cm-2) at T_REF
    A: float             # s-1
    E_lower: float       # cm-1
    v_upper: int
    v_lower: int
    J_lower: int
    branch: str          # "P" or "R"
    g_upper: float
    g_lower: float
    sigma_MHz: float
    F_upper: int | None = None
    F_lower: int | None = None


def einstein_a(S, nu, E_lower, g_upper, Q, T=T_REF):
    """A from the HITRAN relation S = A g′ exp(-c2 E″/T) (1 - exp(-c2 ν/T)) / (8π c ν² Q), abundance 1."""
    return S * 8 * math.pi * C_CM * nu ** 2 * Q / (g_upper * math.exp(-C2 * E_lower / T) * -math.expm1(-C2 * nu / T))


def _err_code(sigma_cm):
    """HITRAN's decade code for a position uncertainty in cm-1: 1 for [0.1, 1) down to 9 for < 1e-8; 0 if unknown."""
    if not (sigma_cm > 0) or sigma_cm >= 1.0:
        return 0
    return min(9, max(1, int(math.floor(-math.log10(sigma_cm))) + 1))


def _fortran_e(x, width=10, digits=3):
    """Fortran E10.3: '1.234E-19'; HITRAN writes the mantissa below 10 with one digit before the point."""
    return f"{x:{width}.{digits}E}"


def _fortran_f(x, width, digits):
    """Fortran Fw.d: drop the leading zero when it does not fit (HITRAN writes γ = .0710 in F5.4)."""
    s = f"{x:.{digits}f}"
    if len(s) > width:
        s = s.replace("0.", ".", 1)
    if len(s) > width:
        raise ValueError(f"{x} does not fit F{width}.{digits}")
    return f"{s:>{width}}"


def format_record(r: Record, molecule=0, gamma_air=0.0, gamma_self=0.0, n_air=0.0, delta_air=0.0) -> str:
    """The 160-character HITRAN 2004 line for one record."""
    iso_id = ISOTOPOLOGUE_IDS[r.iso]
    J_upper = r.J_lower + (1 if r.branch == "R" else -1)
    # global quanta, A15 each: the electronic state and v
    gq_up, gq_low = f"{'B':>8}{r.v_upper:>7d}", f"{'X':>8}{r.v_lower:>7d}"
    # local quanta, A15 each, as HITRAN's diatomic classes: upper 10X A5 (F'); lower 5X A1 I3 A1 A5 (Br J'' Sym F'')
    lq_up = " " * 10 + (f"{r.F_upper:>5.1f}" if r.F_upper is not None else " " * 5)
    lq_low = " " * 5 + r.branch + f"{r.J_lower:>3d}" + " " + (f"{r.F_lower:>5.1f}" if r.F_lower is not None else " " * 5)
    ierr = f"{_err_code(r.sigma_MHz / MHZ_PER_CM)}2" + "0000"
    iref = "0" * 12
    line = (f"{molecule:>2d}{iso_id:1d}{r.nu:12.6f}{_fortran_e(r.S)}{_fortran_e(r.A)}"
            f"{_fortran_f(gamma_air, 5, 4)}{_fortran_f(gamma_self, 5, 3)}{r.E_lower:10.4f}"
            f"{_fortran_f(n_air, 4, 2)}{_fortran_f(delta_air, 8, 6)}"
            f"{gq_up}{gq_low}{lq_up}{lq_low}{ierr}{iref} {r.g_upper:7.1f}{r.g_lower:7.1f}")
    assert len(line) == 160, (len(line), line)
    assert J_upper >= 0
    return line


def records(catalog, isotopologue, nu_min, nu_max, S_min=0.0, hyperfine=True, dJ=0):
    """The transitions of one isotopologue between nu_min and nu_max (cm-1) with S(296 K) >= S_min.

    With ``hyperfine`` each line is replaced by its main (ΔF = ΔJ) components, positions from the model's
    hyperfine calculation (``dJ`` 0, the fast ΔJ = 0 matrix, or 2 with the J ± 2 couplings) and strengths
    the line's, shared in the components' proportions (normalised over the main components)."""
    from .intensity import nuclear_spin_weight
    from .lookup import NEAR_FROM, Line, uncertainty

    master = catalog.master(isotopologue, near_dissociation=nu_max > NEAR_FROM)
    Q = master.partition_function(T_REF)
    lines = master.at(T_REF, nu_min, nu_max, S_min=S_min)
    model = catalog.model(isotopologue) if hyperfine else None
    out = []
    for k in range(len(lines.nu)):
        nu, S, El = float(lines.nu[k]), float(lines.S[k]), float(lines.E_lower[k])
        vu, vl, J = int(lines.v_upper[k]), int(lines.v_lower[k]), int(lines.J_lower[k])
        br = "R" if lines.branch[k] > 0 else "P"
        Ju = J + (1 if br == "R" else -1)
        gns = nuclear_spin_weight(J, isotopologue)
        sigma = uncertainty(Line(isotopologue, br, J, vu, vl, nu, S, El, T_REF))[0]
        comps = None
        if hyperfine:
            try:
                _, comps = model.hyperfine_components(vu, vl, J, br, dJ=dJ, position=False)
                comps = [c for c in comps if c.label]
            except Exception:              # a level beyond every validated hyperfine grid: the line alone
                comps = None
        if not comps:
            gu, gl = (2 * Ju + 1) * gns, (2 * J + 1) * gns
            out.append(Record(isotopologue, nu, S, einstein_a(S, nu, El, gu, Q), El, vu, vl, J, br, gu, gl, sigma))
            continue
        total = sum(c.strength for c in comps)
        for c in comps:
            nu_c, S_c = nu + c.offset / MHZ_PER_CM, S * c.strength / total
            gu, gl = 2 * c.F_upper + 1, 2 * c.F_lower + 1
            out.append(Record(isotopologue, nu_c, S_c, einstein_a(S_c, nu_c, El, gu, Q), El, vu, vl, J, br, gu, gl,
                              sigma, c.F_upper, c.F_lower))
    out.sort(key=lambda r: r.nu)
    return out


def partition_table(catalog, isotopologue, T=None):
    """(T, Q(T)) in HITRAN's tabulated form, 1 K steps; the same Q the line strengths use."""
    master = catalog.master(isotopologue)
    T = np.arange(100.0, 1001.0) if T is None else np.asarray(T, dtype=float)
    return T, np.array([master.partition_function(t) for t in T])


def write(path_stem, catalog, isotopologues, nu_min, nu_max, S_min=0.0, hyperfine=True, dJ=0, molecule=0,
          **broadening):
    """Write <stem>.par (every isotopologue, by wavenumber) and <stem>_q_<iso>.txt; returns the record count."""
    from pathlib import Path

    stem = Path(path_stem)
    recs = []
    for iso in isotopologues:
        recs += records(catalog, iso, nu_min, nu_max, S_min=S_min, hyperfine=hyperfine, dJ=dJ)
        T, Q = partition_table(catalog, iso)
        with open(f"{stem}_q_{iso}.txt", "w") as f:
            f.write(f"# partition function of {iso} (molecule {molecule}, isotopologue {ISOTOPOLOGUE_IDS[iso]}, "
                    f"mass {MASSES[iso]:.6f} u), with nuclear-spin degeneracy; T (K), Q\n")
            for t, q in zip(T, Q):
                f.write(f"{t:7.1f} {q:.6e}\n")
    recs.sort(key=lambda r: r.nu)
    with open(f"{stem}.par", "w") as f:
        for r in recs:
            f.write(format_record(r, molecule=molecule, **broadening) + "\n")
    return len(recs)
