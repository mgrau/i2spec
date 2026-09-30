# Roadmap

**Goal:** an open, installable tool that computes the B–X absorption spectrum of molecular iodine, with
hyperfine structure, from a physical model fitted to all available data, and gives an uncertainty for
every line.

**Reference point:** the model of the Hannover group as implemented in IodineSpec (Knöckel, Bodermann &
Tiemann 2004; Salumbides *et al.* 2008; Bodermann, Knöckel & Tiemann 2002). `docs/research/iodinespec5.md`
compares the two.

## Current state

| part | status | record |
|---|---|---|
| Rovibrational levels | published potentials (Salumbides *et al.* 2008) inside their fitted range; Morse/long-range potentials fitted in this work beyond it; measured level corrections from all frequency-comb-referenced data | `docs/design/parameter-sets.md`, `docs/design/mlr-x.md` |
| X v″ = 18–25 | measured from the Orsay atlas (1982), 15–52 MHz per level | `docs/research/orsay-atlas-11000-14000.md` |
| B v′ = 51–79 | refitted to 2 009 lines of the Orsay atlas Partie IV (1983) and the comb lines there, with level corrections: 20–250 MHz per level | `docs/research/orsay-atlas-19700-20035.md` |
| Hyperfine structure | effective Hamiltonian with ΔJ = ±2 couplings; published parameter formulae plus a table of measured B-state corrections | `docs/design/hyperfine-fit.md`, `docs/research/isotopologue-hyperfine.md` |
| Intensities and spectra | transition-moment matrix elements, partition function, Doppler and Voigt profiles, cell transmission, bound–free continuum, sub-Doppler spectra | `docs/research/intensity-inputs.md`, `docs/research/continuum-model.md`, `docs/research/sub-doppler.md` |
| Validation | 35 precision data sets and three atlases; better than IodineSpec5 on 19 of 22 comparable sets | `docs/research/spectra-validation.md`, `docs/research/iodinespec5.md` |
| Uncertainties | per line, from held-out validation | `docs/design/uncertainty.md` |
| Interfaces | Python package, command line, terminal browser, static web explorer, documentation site | `README.md`, `docs/design/web-app.md` |
| Line list | ¹²⁷I₂, ¹²⁹I₂, ¹²⁷I¹²⁹I, 11 000–20 100 cm⁻¹, levels to the B dissociation limit | `docs/design/near-dissociation.md` |

## Planned

1. **B levels near dissociation.** Done for the levels the Orsay atlas Partie IV measured (i2spec2026k).
   Left: the quasi-bound levels above the asymptote (about 70 atlas lines), a coupled-channel treatment
   of the 1g perturbation, and a comb measurement anywhere at v′ = 63–79.
2. **First release (v0.1).** Package data shipped with the package, a PyPI release, a citable version
   (Zenodo), and continuous testing.
3. **Uncertainty from the fit.** Parameter covariance and a model-discrepancy term, validated by blocked
   cross-validation, in place of the validation-based estimates.
4. **Isotope shifts.** A refit of the Born–Oppenheimer corrections if the isotopologue difference data of
   Salumbides *et al.* (2008) become available (`docs/research/data-availability.md`).
5. **A compiled forward model** shared by the Python package and the web explorer (Rust, compiled to
   WebAssembly), once the model form is stable.
6. **Explorer features:** sub-Doppler spectra, export of line lists, and a search by laser wavelength
   (including harmonics of a fundamental). Done (`documentation/explorer.md`). Left: crossover
   resonances in the explorer, which need the (I, F) level labels of the weak components in the export.

## Decisions

| decision | choice |
|---|---|
| package name | `i2spec` |
| licence | MIT, copyright "i2spec contributors" |
| language | Python for the model and the fits; a compiled forward model later |
| potential form | the published power-series potentials inside their range; Morse/long-range potentials outside it |
| data sources | published measurements and their supplements only |
