# Overview of the calculation

## Steps

A line list is computed in five steps, once per isotopologue. The result does not depend on
temperature and is stored. Spectra for a given temperature and absorption cell are computed from it
in a sixth step.

<figure class="fig" markdown="span">
<div class="svg" data-svg="figures/pipeline.svg" role="img" aria-label="The six steps of the calculation, each with a small plot from the model, and the inputs each step takes"></div>
<figcaption markdown="span">**Figure 1.** The calculation in six steps, each drawn from the model; the inputs of each step are listed beneath it. (1) Levels of X and B, every sixth v at J = 0. (2) Corrections to the B levels at J′ = 40: per-level fits to comb data (points, 1σ) and the Gaussian process that fills levels without their own data (line and band, v′ = 3–35). (3) The 32–0 band as a Fortrat diagram. (4) The 15 strong hyperfine components of R(56) 32–0 under a 300 K Doppler profile. (5) Band strengths at 300 K for v″ = 0, 1 and 2. (6) The cross section near 532 nm and the transmission of a 10 cm cell at 20 °C.</figcaption>
</figure>

| step | quantity | section |
|---|---|---|
| 1 | level energies E(v, J) of X and B | [Rovibrational levels](levels.md) |
| 2 | empirical corrections to individual levels | [Rovibrational levels](levels.md#empirical-level-corrections) |
| 3 | positions of all allowed P and R lines | [Line positions](lines.md) |
| 4 | hyperfine components of each line | [Hyperfine structure](hyperfine.md) |
| 5 | temperature-independent line strengths | [Line intensities and spectra](spectra.md) |
| 6 | cross section and transmission at a given temperature | [Line intensities and spectra](spectra.md#absorption-spectrum) |

<figure class="fig" markdown="span">
<div class="svg" data-svg="figures/potentials.svg" role="img" aria-label="Potential-energy curves of the X and B states with vibrational levels and the 532 nm transition"></div>
<figcaption markdown="span">**Figure 2.** Potential-energy curves of X¹Σg⁺ and B³Π(0u⁺) with their J = 0 vibrational levels, drawn between the classical turning points (all levels of X up to v″ = 17, every second level of B up to v′ = 43, then every third level). Levels outside the range of the 2008 fit (Salumbides et al.) are drawn faint. The arrow marks R(56) 32–0 at 532 nm. X dissociates to two ground-state atoms, I(²P3/2) + I(²P3/2); B dissociates to I(²P3/2) + I(²P1/2).</figcaption>
</figure>

## Notation and units

| symbol | meaning |
|---|---|
| v″, J″ | vibrational and rotational quantum numbers of the lower level (X state) |
| v′, J′ | the same for the upper level (B state) |
| R(J″) v′–v″, P(J″) v′–v″ | a line with J′ = J″ + 1 (R branch) or J′ = J″ − 1 (P branch) |
| ν | wavenumber in cm⁻¹; 1 cm⁻¹ = 29 979.2458 MHz |
| E″ | energy of the lower level above X(v″ = 0, J″ = 0), in cm⁻¹ |
| a1, a2, … | main hyperfine components in order of increasing frequency |

Energies of both states are measured from the minimum of the X potential. Wavelengths are vacuum
wavelengths unless stated otherwise.
