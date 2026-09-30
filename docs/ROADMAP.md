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
| Rovibrational levels | published potentials (Salumbides *et al.* 2008) inside their fitted range; Morse/long-range potentials fitted in this work beyond it; measured level corrections from all frequency-comb-referenced data; a Gaussian-process correction for the visible B levels without their own (i2spec2026n) | `docs/design/parameter-sets.md`, `docs/design/mlr-x.md` |
| X v″ = 18–25 | measured from the Orsay atlas (1982), 15–52 MHz per level | `docs/research/orsay-atlas-11000-14000.md` |
| B v′ = 51–79 | refitted to all 2 079 unblended lines of the Orsay atlas Partie IV (1983), including 70 from quasi-bound levels above the asymptote, and the comb lines there, with level corrections: 20–250 MHz per level | `docs/research/orsay-atlas-19700-20035.md`, `docs/research/quasibound-b.md` |
| Hyperfine structure | effective Hamiltonian with ΔJ = ±2 couplings; published parameter formulae plus a table of measured B-state corrections from every set in use (i2spec2026n) | `docs/design/hyperfine-fit.md`, `docs/research/isotopologue-hyperfine.md` |
| Intensities and spectra | transition-moment matrix elements, partition function, Doppler and Voigt profiles, cell transmission, bound–free continuum, sub-Doppler spectra | `docs/research/intensity-inputs.md`, `docs/research/continuum-model.md`, `docs/research/sub-doppler.md` |
| Validation | 63 precision data sets, four atlases and one set of laser scans; better than IodineSpec5 on 19 of 22 comparable sets; a held-out comparison of correction schemes (`docs/research/model-bakeoff.md`) | `docs/research/spectra-validation.md`, `docs/research/iodinespec5.md` |
| Uncertainties | per line, from held-out validation | `docs/design/uncertainty.md` |
| Interfaces | Python package, command line, terminal browser, static web explorer, documentation site | `README.md`, `docs/design/web-app.md` |
| Line list | ¹²⁷I₂, ¹²⁹I₂, ¹²⁷I¹²⁹I, 11 000–20 100 cm⁻¹, levels to the B dissociation limit | `docs/design/near-dissociation.md` |

## Planned

1. **B levels near dissociation.** Done for every level the Orsay atlas Partie IV measured, including the
   quasi-bound ones (i2spec2026n). Left: their lines in the exported list (the 12 Å intensity grid stops at
   the first box state), a coupled-channel treatment of the 1g perturbation, and a comb measurement
   anywhere at v′ = 63–79.
2. **First release (v0.1).** Package data shipped with the package, a PyPI release, a citable version
   (Zenodo), and continuous testing.
3. **Uncertainty from the fit.** Done for the level corrections (covariance and a discrepancy term,
   i2spec2026l) and the Gaussian-process levels (posterior σ, i2spec2026n). Left: the potentials' own
   covariance, and a level-sharing correction for the X state and the near infrared, where the bake-off
   found no scheme that predicts an unmeasured level.
4. **Isotope shifts.** A refit of the Born–Oppenheimer corrections if the isotopologue difference data of
   Salumbides *et al.* (2008) become available (`docs/research/data-availability.md`).
5. **The Kato Doppler-free atlas** (526–667 nm, about 3 MHz), requested from its project leader; the
   biggest data set not yet in the fit. The B–X band strength at 540–630 nm, if an absolute measurement with
   a known column appears (`docs/research/bx-band-strength.md`).
6. **A compiled forward model** shared by the Python package and the web explorer (Rust, compiled to
   WebAssembly), once the model form is stable.
7. **Explorer features:** sub-Doppler spectra, export of line lists, and a search by laser wavelength
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
