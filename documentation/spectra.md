# Line intensities and spectra

## Line strength

The integrated absorption cross section of a line is (Tellinghuisen 2011)

$$
S = \frac{2\pi^2\,\nu}{3\varepsilon_0 h c}\;\frac{s_{J'J''}}{2J''+1}\;\left|\langle v'J'|\,\mu_e(R)\,|v''J''\rangle\right|^2\; f(v'',J''), \tag{1}
$$

where s = J″ for P lines and s = J″ + 1 for R lines is the Hönl–London factor, μₑ(R) is the electronic
transition moment and f is the fractional population of the lower level,

$$
f = \frac{g_\text{ns}(J'')\,(2J''+1)\,e^{-c_2 E''/T}}{Q(T)},\qquad Q(T) = \sum_{v'',J''} g_\text{ns}(J'')\,(2J''+1)\,e^{-c_2E''/T}, \tag{2}
$$

with g<sub>ns</sub> the nuclear-spin statistical weight and c₂ = hc/k<sub>B</sub> = 1.4388 cm K. The
transition moment is that of Tellinghuisen (2011),
|μₑ(R)|² = 21.5 exp[−1.18 (R − 3.65)²] / R² D², with R in Å.

The matrix element in Eq. (1) is computed with the J-dependent vibrational wavefunctions of both
states, so the variation of the intensity with J within a band (Herman–Wallis effect) is included.
Only f depends on temperature. The line list therefore stores the temperature-independent quantity
S₀ = S Q(T) exp(c₂E″/T), from which S is obtained for any T.

<figure class="fig" markdown="span">
<div class="svg" data-svg="figures/wavefunctions.svg" role="img" aria-label="Vibrational wavefunctions in the X and B potentials and the transition moment"></div>
<figcaption markdown="span">**Figure 1.** Vibrational wavefunctions of X v″ = 0 and of B v′ = 5 and v′ = 32, drawn at their energies, and the transition moment μₑ(R) (lower panel). The band strength is proportional to the square of the overlap ⟨v′|μₑ|v″⟩. The v″ = 0 wavefunction is localised near R = 2.67 Å (shaded). The overlap is largest for upper levels whose inner turning point lies in this region, where the upper wavefunction has its largest amplitude; for v″ = 0 this is the case near v′ = 30, at 500–540 nm.</figcaption>
</figure>

Upper levels are included up to the last bound level of each J′. Absorption to the continuum above the
B-state dissociation limit is computed separately by the Python package (`i2spec.continuum`).

## Absorption spectrum

At temperature T each line is Doppler broadened to a Gaussian with full width at half maximum

$$
\Delta\nu_D = \nu\,\sqrt{\frac{8k_BT\ln 2}{mc^2}} \approx 4.5\times10^{-8}\,\nu\,\sqrt{T/\mathrm{K}} \quad (^{127}\mathrm{I}_2), \tag{3}
$$

about 430 MHz at 532 nm and 293 K. The absorption cross section is

$$
\sigma(\nu) = \sum_i S_i\,g(\nu-\nu_i), \tag{4}
$$

with g a normalised Gaussian (a Voigt profile, including pressure broadening, in the Python
package). A cell of length L transmits

$$
\mathcal T(\nu) = \exp\left[-\sigma(\nu)\,N\,L\right], \tag{5}
$$

where the number density N is determined by the iodine vapour pressure at the cold-finger temperature
(Tellinghuisen 2011, Eq. 8).

<figure class="fig" markdown="span">
<div class="svg" data-svg="figures/spectrum.svg" role="img" aria-label="Cross section and cell transmission near 532.25 nm"></div>
<figcaption markdown="span">**Figure 2.** Absorption cross section near 532 nm at 20 °C (upper panel) and transmission of a 10 cm cell with the cold finger at 20 °C (lower panel). Each peak is a single rovibronic line; the strongest are labelled.</figcaption>
</figure>

The Python package also computes sub-Doppler spectra: Lamb dips, crossover resonances between
components that share a level, and the third-harmonic line shape to which the BIPM frequencies refer.
