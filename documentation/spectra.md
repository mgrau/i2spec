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
states. This includes the variation of the intensity with J within a band (Herman–Wallis effect). The
wavefunctions are the model's, solved on a uniform grid (a sinc discrete-variable representation, 0.005 Å) on
which the transition moment is also evaluated. The line positions and E″ come from the model's levels.

The strengths of bands whose upper level lies near a node of the overlap with the lower level depend sharply on
where each potential lies in R. Line positions alone leave this partly open, because the shape of a potential
trades against its rotational function; measured intensities test it. Three Fourier-transform cell spectra were
compared: Debus 2023 and the two pyodine cells, with 28 000–33 000 unblended lines each. The line strength implied
by each line's integrated absorbance, relative to its spectral window's column density, scatters about the model by
1.7 %, 3.3 % and 4.5 % (robust standard deviations of the logarithm about its median). With the wavefunctions of
the 2008 potentials the scatter is 2.3 %, 4.1 % and 4.7 %. The v″ = 1 bands near the overlap's node show the
difference best. 41-1 lies 7–9 % below the 2008 potentials' strengths in all three spectra and within 2 % of the
model's. 40-1 shows the same in the SONG spectrum (−10 % against +1 %) but not in the Fischer spectrum (+4 %
against +16 %).

Of the factors in Eq. (1), only f depends on temperature. The line list therefore stores
S₀ = S Q(T) exp(c₂E″/T), from which S follows for any T.

<figure class="fig" markdown="span">
<div class="svg" data-svg="figures/wavefunctions.svg" role="img" aria-label="Vibrational wavefunctions in the X and B potentials and the transition moment"></div>
<figcaption markdown="span">**Figure 1.** Vibrational wavefunctions of X v″ = 0 and of B v′ = 5 and v′ = 32, drawn at their energies, and the transition moment μₑ(R) (lower panel). The band strength is proportional to the square of the overlap ⟨v′|μₑ|v″⟩. The v″ = 0 wavefunction is localized near R = 2.67 Å (shaded). The overlap is largest for upper levels whose inner turning point lies in this region, where the upper wavefunction has its largest amplitude; for v″ = 0 this is the case near v′ = 30, at 500–540 nm.</figcaption>
</figure>

Upper levels are included up to the last bound level of each J′. The Python package computes the absorption to the
continuum above the B-state dissociation limit separately (`i2spec.continuum`).

## Absorption spectrum

At temperature T each line is Doppler broadened to a Gaussian with full width at half maximum

$$
\Delta\nu_D = \nu\,\sqrt{\frac{8k_BT\ln 2}{mc^2}} \approx 4.5\times10^{-8}\,\nu\,\sqrt{T/\mathrm{K}} \quad (^{127}\mathrm{I}_2), \tag{3}
$$

about 430 MHz at 532 nm and 293 K. The absorption cross section is

$$
\sigma(\nu) = \sum_i S_i\,g(\nu-\nu_i), \tag{4}
$$

with g a normalized Gaussian; in the Python package g is a Voigt profile that includes pressure broadening.
A cell of length L transmits

$$
\mathcal T(\nu) = \exp\left[-\sigma(\nu)\,N\,L\right], \tag{5}
$$

where the number density N is determined by the iodine vapor pressure at the cold-finger temperature
(Tellinghuisen 2011, Eq. 8).

<figure class="fig" markdown="span">
<div class="svg" data-svg="figures/spectrum.svg" role="img" aria-label="Cross section and cell transmission near 532.25 nm"></div>
<figcaption markdown="span">**Figure 2.** Absorption cross section near 532 nm at 20 °C (upper panel) and transmission of a 10 cm cell with the cold finger at 20 °C (lower panel). Each peak is a single rovibronic line; the strongest are labeled.</figcaption>
</figure>

The Python package also computes sub-Doppler spectra: Lamb dips, crossover resonances between
components that share a level, and the third-harmonic line shape to which the BIPM frequencies refer.

## Worked example

!!! example "R(56) 32–0: strength and Doppler width at 20 °C"
    At T = {{ we.T }} K, Eqs. (1)–(3) give for the line:

    | quantity | value |
    |---|--:|
    | E″, the [lower level](levels.md#worked-example) above X(0, 0) | {{ we.e_lower }} cm⁻¹ |
    | nuclear-spin weight g<sub>ns</sub> (even J″) | {{ we.g_ns }} |
    | Hönl–London factor s = J″ + 1 | {{ we.honl_london }} |
    | partition function Q(T) | {{ we.Q }} |
    | fractional population f of X(0, 56) | {{ we.population }} |
    | line strength S | {{ we.S }} cm per molecule |
    | Doppler width (FWHM), Eq. (3) | {{ we.doppler }} MHz |

    The strength is shared among the [hyperfine components](hyperfine.md#worked-example): a10 carries a fraction
    {{ we.comp_strength }} of it. The pattern spans {{ we.span }} MHz, about twice the Doppler width, so the
    components overlap in a Doppler-limited spectrum and are resolved only by sub-Doppler methods.

    ```python
    from i2spec.lookup import Catalog

    line = Catalog().line("R(56) 32-0", temperature=293.15)  # the first call builds the line list (minutes)
    line.strength          # {{ we.S_raw }} cm per molecule
    line.E_lower           # {{ we.E_strength_raw }} cm⁻¹
    line.doppler_fwhm_MHz  # {{ we.doppler_raw }} MHz
    ```
