"""Sub-Doppler (saturated-absorption) spectra of single B-X lines.

Saturation spectroscopy punches narrow resonances through the Doppler profile of a line:

* a **Lamb dip** at every hyperfine component, from the v_z = 0 velocity class;
* a **crossover** halfway between two components that share a level, from the two velocity classes
  v_z = ±(ν_j − ν_i)/2, in which one component is resonant with the forward beam and the other with
  the backward beam. Two components that share no level produce no crossover.

Amplitudes come from the weak-saturation rate-equation limit, in which the signal is bilinear in
the pump and probe transition strengths (Letokhov & Chebotayev 1977; Demtröder, *Laser
Spectroscopy*, §2.3 and §7.2):

    dip(i)        = S_i**p
    crossover(ij) = 2 · w_ij · (S_i S_j)**(p/2) · exp(−ln2 · Δ_ij² / Δν_D²)

* p = ``exponent`` is 2 in that limit, because pump and probe each contribute one factor of the
  transition strength. Strong saturation flattens it towards 1; the calculated traces of
  Salumbides et al. (2006) correspond to p = 1.
* The exponential is the fraction of molecules in the velocity class a crossover uses, and the
  factor 2 counts the two classes. It is what makes crossovers fade out beyond the Doppler width.
* w_ij is 1 for a shared lower level (V-type) and ``lambda_weight`` for a shared upper level
  (Λ-type). The two mechanisms differ, so the default of 1 is a deliberate assumption, not a
  result: see docs/research/sub-doppler.md.

The amplitudes are relative and sum to 1, so this module gives the *shape* of a Doppler-free
spectrum, not its depth. docs/research/sub-doppler.md derives all of it and lists the assumptions.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import log, sqrt

import numpy as np

from .constants import ATOMIC_MASS, ISOTOPOLOGUES, MHZ_PER_CM
from .spectrum import doppler_fwhm

#: Crossovers are dropped beyond this multiple of the Doppler FWHM, where the velocity class the
#: crossover needs holds exp(-ln2 · reach²) of the molecules: 3 leaves 2e-4 of the peak density.
REACH = 3.0


@dataclass(frozen=True)
class Resonance:
    """One Doppler-free feature: a Lamb dip, or a crossover between two components."""

    offset: float  # MHz, relative to the hyperfine-free rovibronic line
    amplitude: float  # relative; the resonances of one line sum to 1
    kind: str  # "dip" or "crossover"
    sharing: str | None  # crossovers: "lower" (V-type) or "upper" (Λ-type); None for a dip
    components: tuple[int, int]  # indices into the component list; (i, i) for a dip
    labels: tuple[str | None, str | None]  # component labels, None for components without one


def shared_level(a, b):
    """"lower", "upper" or None: which level two hyperfine components have in common.

    Levels are identified by (I, F). I is the dominant total nuclear spin, which the quadrupole
    coupling mixes at high J, so this labelling can fail where two levels of one F have nearly
    equal I content (docs/research/sub-doppler.md).
    """
    lower = (a.I_lower, a.F_lower) == (b.I_lower, b.F_lower)
    upper = (a.I_upper, a.F_upper) == (b.I_upper, b.F_upper)
    if lower and upper:
        return None  # the same transition twice
    return "lower" if lower else "upper" if upper else None


def resonances(components, doppler_width, *, exponent=2.0, lambda_weight=1.0, reach=REACH,
               threshold=1e-5) -> list[Resonance]:
    """Lamb dips and crossovers of one line, strongest first, with amplitudes summing to 1.

    components: from RovibronicModel.hyperfine_components. doppler_width: Doppler FWHM in MHz, the
    same units as the component offsets. threshold: drop resonances weaker than this fraction of
    the strongest one.
    """
    found = [Resonance(c.offset, c.strength**exponent, "dip", None, (i, i), (c.label, c.label))
             for i, c in enumerate(components)]
    for i, a in enumerate(components):
        for j, b in enumerate(components[i + 1:], start=i + 1):
            sharing = shared_level(a, b)
            gap = b.offset - a.offset
            if sharing is None or abs(gap) > reach * doppler_width:
                continue
            weight = 1.0 if sharing == "lower" else lambda_weight
            amplitude = (2 * weight * (a.strength * b.strength) ** (exponent / 2)
                         * np.exp(-log(2) * (gap / doppler_width) ** 2))
            found.append(Resonance((a.offset + b.offset) / 2, float(amplitude), "crossover", sharing,
                                   (i, j), (a.label, b.label)))
    strongest = max((r.amplitude for r in found), default=0.0)
    kept = [r for r in found if r.amplitude > threshold * strongest]
    total = sum(r.amplitude for r in kept)
    return sorted((Resonance(r.offset, r.amplitude / total, r.kind, r.sharing, r.components, r.labels)
                   for r in kept), key=lambda r: -r.amplitude)


def doppler_width(isotopologue, nu0, T):
    """Doppler FWHM in MHz at line frequency nu0 (MHz) and temperature T (K)."""
    a, b = ISOTOPOLOGUES[isotopologue]
    return float(doppler_fwhm(nu0 / MHZ_PER_CM, T, ATOMIC_MASS[a] + ATOMIC_MASS[b]) * MHZ_PER_CM)


def line_resonances(model, v_upper, v_lower, J_lower, branch, T, *, dJ=2, **options):
    """(nu0, resonances) of one B-X line at temperature T (K); nu0 is the line centre in MHz."""
    nu0, components = model.hyperfine_components(v_upper, v_lower, J_lower, branch, dJ=dJ)
    return nu0, resonances(components, doppler_width(model.isotopologue, nu0, T), **options)


def _lorentzian(x, fwhm):
    return 1.0 / (1.0 + (2 * np.asarray(x, dtype=float) / fwhm) ** 2)


def lineshape(x, fwhm, *, harmonic=0, modulation=0.0, points=256):
    """Saturation line shape at detuning x (MHz), peak 1 for an unmodulated resonance.

    A Lorentzian of the given FWHM, or, when the laser frequency is modulated sinusoidally with
    peak-to-peak width ``modulation``, the in-phase ``harmonic``-th Fourier component of it:

        S_n(x) = (2/π) ∫₀^π L(x + (m/2) cos θ) cos(nθ) dθ.

    For m ≪ FWHM this tends to the n-th derivative of the Lorentzian. The BIPM mise en pratique
    defines the recommended iodine frequencies for third-harmonic detection (n = 3) with about
    1 MHz peak-to-peak modulation; the 3f signal crosses zero at the centre of each resonance,
    which is what a laser locks to.
    """
    x = np.asarray(x, dtype=float)
    if harmonic == 0 and modulation == 0.0:
        return _lorentzian(x, fwhm)
    theta = (np.arange(points) + 0.5) * np.pi / points  # midpoint rule, exact for smooth periodic integrands
    weight = np.cos(harmonic * theta) * (2.0 / points)
    return (weight * _lorentzian(x[..., None] + 0.5 * modulation * np.cos(theta), fwhm)).sum(-1)


def signal(nu, resonances, fwhm, *, harmonic=0, modulation=0.0, points=256, wings=100.0):
    """Doppler-free signal on ``nu`` (MHz, same origin as the resonance offsets).

    Every resonance contributes ``amplitude × lineshape(ν − offset)``, so the result is positive at
    each resonance for harmonic = 0: it is the increase in probe transmission, not the absorption.
    ``fwhm`` (MHz) is the homogeneous width the caller wants: natural width, transit-time and
    pressure broadening and saturation broadening together. Modulated shapes are interpolated from
    an internal kernel of step FWHM/256, truncated at ``wings`` half-widths, which holds the error
    below about 1e-5 of the peak even across the steep zero crossing a laser locks to.
    """
    nu = np.asarray(nu, dtype=float)
    offsets = np.array([r.offset for r in resonances])
    amplitudes = np.array([r.amplitude for r in resonances])
    if offsets.size == 0:
        return np.zeros_like(nu)
    if harmonic == 0 and modulation == 0.0:
        return (amplitudes * _lorentzian(nu[..., None] - offsets, fwhm)).sum(-1)
    step = fwhm / 256
    grid = np.arange(-wings * (fwhm + modulation), wings * (fwhm + modulation) + step, step)
    kernel = lineshape(grid, fwhm, harmonic=harmonic, modulation=modulation, points=points)
    out = np.zeros_like(nu)
    for offset, amplitude in zip(offsets, amplitudes):
        out += amplitude * np.interp(nu - offset, grid, kernel, left=0.0, right=0.0)
    return out


def transit_time_width(diameter_mm, T, isotopologue="127I2"):
    """Transit-time broadening FWHM in MHz for a beam of the given diameter (Demtröder §3.4).

    Γ_tt = 2 √(2 ln2) · v̄ / (π d) for a Gaussian beam of 1/e² diameter d; it is the floor on the
    width of a Doppler-free resonance in a cell, a few tens of kHz for a millimetre beam.
    """
    a, b = ISOTOPOLOGUES[isotopologue]
    mass = (ATOMIC_MASS[a] + ATOMIC_MASS[b]) * 1.66053906892e-27
    speed = sqrt(8 * 1.380649e-23 * T / (np.pi * mass))
    return 2 * sqrt(2 * log(2)) * speed / (np.pi * diameter_mm * 1e-3) * 1e-6
