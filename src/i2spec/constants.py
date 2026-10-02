"""Physical constants and isotope data. Energies in cm⁻¹, lengths in Å, masses in u.

Constants are pinned here rather than taken from scipy.constants. Switching SciPy's CODATA
version (2018 -> 2022) changes ħ²/(2u) by 1.4e-9, which moves B-X lines by tens of kHz.
"""

from math import pi

_H = 6.62607015e-34  # Planck constant, J s (exact)
_C = 299792458.0  # speed of light, m/s (exact)
_U = 1.66053906892e-27  # atomic mass constant, kg (CODATA 2022)

#: MHz per cm⁻¹.
MHZ_PER_CM = _C / 1e4

#: ħ²/(2 u Å²) in cm⁻¹.
HBAR2_2U = _H / (8 * pi**2 * _C * _U) * 1e18

#: Relative atomic masses (NIST Atomic Weights and Isotopic Compositions).
ATOMIC_MASS = {127: 126.9044719, 129: 128.9049837}

#: Nuclear spins.
NUCLEAR_SPIN = {127: 2.5, 129: 3.5}

#: Isotopologue label -> mass numbers of the two nuclei.
ISOTOPOLOGUES = {"127I2": (127, 127), "129I2": (129, 129), "127I129I": (127, 129)}

#: Isotopologue to which the Hannover Born-Oppenheimer-correction functions refer.
REFERENCE_ISOTOPOLOGUE = "127I2"


def reduced_mass(isotopologue: str) -> float:
    """Reduced mass in u, from atomic masses."""
    a, b = ISOTOPOLOGUES[isotopologue]
    ma, mb = ATOMIC_MASS[a], ATOMIC_MASS[b]
    return ma * mb / (ma + mb)

#: The parameter set the model uses unless told otherwise: the published Hannover potentials plus the
#: changes i2spec has validated (see the file's "changes" and docs/design/parameter-sets.md).
DEFAULT_PARAMETERS = "i2spec2026o"
#: The published model (Salumbides et al. 2008, unchanged), kept for comparison
PUBLISHED_PARAMETERS = "hannover2008"

#: I(2P3/2) + I(2P1/2) above 2 I(2P3/2), cm⁻¹: the B-state asymptote sits this far above the X-state one.
ATOMIC_SPLITTING = 7602.9762
