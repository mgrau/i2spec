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
| – | which data fix which parameters | [The fit](fit.md) |

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

## Worked example

One line is followed through every step: R(56) 32–0 of ¹²⁷I₂, the 532 nm frequency standard, and its hyperfine
component a10. Each page of the model ends with a box giving the numbers of its step for this line, computed by the
package when this site was built.

| step | for R(56) 32–0 | page |
|---|---|---|
| 1 | E<sub>X</sub>(0, 56) = {{ we.ex }} cm⁻¹, E<sub>B</sub>(32, 57) = {{ we.eb }} cm⁻¹ | [levels](levels.md#worked-example) |
| 2 | ν = {{ we.nu }} cm⁻¹ ({{ we.nm }} nm in vacuum) | [line positions](lines.md#worked-example) |
| 3 | a10 at {{ we.a10 }} MHz from ν: {{ we.comp_MHz }} MHz | [hyperfine](hyperfine.md#worked-example) |
| 4–5 | S = {{ we.S }} cm per molecule at 20 °C, Doppler width {{ we.doppler }} MHz | [spectra](spectra.md#worked-example) |
| 6 | 1σ: {{ we.sigma }} MHz for ν, {{ we.sigma_comp }} MHz for a10 | [uncertainty](uncertainty.md#worked-example) |
| fit | R(56) 32–0's own measurements were held out of the fit | [the fit](fit.md#worked-example) |

The same numbers from Python, with the public interface:

```python
from i2spec import MHZ_PER_CM, ResearchModel, RovibronicModel
from i2spec.lookup import Catalog

m = ResearchModel("i2spec")                 # the released model, ¹²⁷I₂
levels = RovibronicModel("127I2")           # its levels, cm⁻¹ above the X minimum

levels.energy("X", 0, 56)                   # {{ we.ex_raw }}
levels.energy("B", 32, 57)                  # {{ we.eb_raw }}
nu = m.transition(32, 0, 56, "R")           # {{ we.nu_raw }} cm⁻¹
nu * MHZ_PER_CM                             # {{ we.nu_MHz_raw }} MHz
m.frequency("R(56) 32-0", "a10")            # {{ we.comp_MHz_raw }} MHz
m.uncertainty("R(56) 32-0").sigma_MHz       # {{ we.sigma_raw }} MHz
m.uncertainty("R(56) 32-0", "a10")          # {{ we.sigma_comp }} MHz, with its parts
Catalog().line("R(56) 32-0").strength       # {{ we.S_raw }} cm per molecule at 20 °C
```
