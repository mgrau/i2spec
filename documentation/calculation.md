# Overview of the calculation

## Steps

The calculation has six steps. Steps 1–4 compute the line list: the levels, the positions of the lines,
their hyperfine components and their strengths. Step 5 computes a spectrum from it for a given temperature
and absorption cell. Step 6 gives the uncertainty of every position and hyperfine component in the list.
The line list and its uncertainties do not depend on temperature; they are computed once per isotopologue
and stored.

<figure class="fig" markdown="span">
<div class="svg" data-svg="figures/pipeline.svg" role="img" aria-label="The six steps of the calculation, each with a small plot from the model, and the inputs each step takes"></div>
<figcaption markdown="span">**Figure 1.** The calculation in six steps, each drawn from the model (¹²⁷I₂); the inputs of each step are listed beneath it. (1) Levels of X and B, every sixth v at J = 0. (2) The 32–0 band as a Fortrat diagram. (3) The 15 strong hyperfine components of R(56) 32–0 under a 300 K Doppler profile. (4) Band strengths at 300 K for v″ = 0, 1 and 2. (5) The cross section near 532 nm and the transmission of a 10 cm cell at 20 °C. (6) The 1σ uncertainty of the R(20) lines of the progressions v′–0 and 5–v″: below 1 MHz in the core, growing beyond it as the data thin out.</figcaption>
</figure>

| step | quantity | section |
|---|---|---|
| 1 | level energies E(v, J) of X and B | [Rovibrational levels](levels.md) |
| 2 | positions of all allowed P and R lines | [Line positions](lines.md) |
| 3 | hyperfine components of each line | [Hyperfine structure](hyperfine.md) |
| 4 | temperature-independent line strengths | [Line intensities and spectra](spectra.md) |
| 5 | cross section and transmission at a given temperature | [Line intensities and spectra](spectra.md#absorption-spectrum) |
| 6 | the uncertainty of every position and hyperfine component | [Uncertainty estimates](uncertainty.md) |

## Notation and units

| symbol | meaning |
|---|---|
| v″, J″ | vibrational and rotational quantum numbers of the lower level (X state) |
| v′, J′ | the same for the upper level (B state) |
| R(J″) v′–v″, P(J″) v′–v″ | a line with J′ = J″ + 1 (R branch) or J′ = J″ − 1 (P branch) |
| ν | wavenumber in cm⁻¹; 1 cm⁻¹ = 29 979.2458 MHz |
| E″ | energy of the lower level above X(v″ = 0, J″ = 0), in cm⁻¹ |
| a1, a2, …, a15 or a21 | main hyperfine components of a line, in order of increasing frequency: a1–a15 for even J″ and a1–a21 for odd J″ in ¹²⁷I₂ (a1–a28 and a1–a36 in ¹²⁹I₂, a1–a48 in ¹²⁷I¹²⁹I) |

Energies of both states are measured from the minimum of the X potential. Wavelengths are vacuum
wavelengths unless stated otherwise.
