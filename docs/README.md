# Notes

These notes record how the model was built and tested: the literature and data it rests on, the design
of each part, and the validation behind every number the package reports. They are written as working
records, in the order the work was done, and some describe intermediate states that later notes
supersede; `design/parameter-sets.md` gives the lineage of the current parameter set. For a description
of the model as it stands, see the documentation site (`documentation/`).

## Research

| note | subject |
|---|---|
| [`01-software-survey.md`](research/01-software-survey.md) | existing software for I₂ spectra, and general diatomic codes |
| [`02a-broadband-data.md`](research/02a-broadband-data.md) | Doppler-limited data: atlases, cross sections, line shapes |
| [`02b-precision-data.md`](research/02b-precision-data.md) | sub-Doppler and frequency-comb measurements |
| [`03a-theory-physics.md`](research/03a-theory-physics.md) | review of the physics: potentials, long range, BO breakdown, hyperfine structure, intensities |
| [`03b-methods.md`](research/03b-methods.md) | review of fitting and uncertainty-quantification methods |
| [`data-availability.md`](research/data-availability.md) | what the repository contains, and data that exist but are not available |
| [`hannover-model-reproduction.md`](research/hannover-model-reproduction.md) | reproducing the published potentials and hyperfine formulae |
| [`isotopologue-hyperfine.md`](research/isotopologue-hyperfine.md) | hyperfine parameters of ¹²⁹I₂ and ¹²⁷I¹²⁹I |
| [`bipm-hyperfine-tables.md`](research/bipm-hyperfine-tables.md) | the BIPM *mise en pratique* tables |
| [`intensity-inputs.md`](research/intensity-inputs.md) | transition moment, partition function, line strengths |
| [`continuum-model.md`](research/continuum-model.md) | bound–free absorption |
| [`cross-section-data.md`](research/cross-section-data.md) | absolute cross sections |
| [`spectra-validation.md`](research/spectra-validation.md) | the model against the Salami & Ross atlas and absolute cross sections |
| [`sub-doppler.md`](research/sub-doppler.md) | Lamb dips, crossover resonances, third-harmonic line shapes |
| [`nir-anchor-data-2000.md`](research/nir-anchor-data-2000.md), [`nir-model-2010.md`](research/nir-model-2010.md) | the near-infrared data and local models |
| [`atlas-line-positions.md`](research/atlas-line-positions.md) | line positions fitted to atlas spectra |
| [`orsay-atlas-11000-14000.md`](research/orsay-atlas-11000-14000.md) | the Orsay atlas and the X levels v″ = 18–25 |
| [`beyond-hannover.md`](research/beyond-hannover.md) | where the published model can be improved, and how |
| [`iodinespec5.md`](research/iodinespec5.md) | comparison with IodineSpec5 |

## Design

| note | subject |
|---|---|
| [`observations.md`](design/observations.md) | the format of the measurement data sets |
| [`parameter-sets.md`](design/parameter-sets.md) | every parameter set, what changed and why |
| [`mlr-x.md`](design/mlr-x.md) | the Morse/long-range potentials |
| [`fitting.md`](design/fitting.md), [`global-fit.md`](design/global-fit.md) | fitting the potentials |
| [`hyperfine-fit.md`](design/hyperfine-fit.md) | fitting the hyperfine parameters |
| [`near-dissociation.md`](design/near-dissociation.md) | levels up to the B dissociation limit |
| [`uncertainty.md`](design/uncertainty.md) | per-line uncertainty estimates |
| [`web-app.md`](design/web-app.md) | the web explorer |

## Terms used in the notes

The notes refer to the plan the work followed:

- **Stages 1–3:** the software survey, the data catalogue and the theory review.
- **Phase A:** an open implementation of the published model. **Phase B:** a refit of the same model form
  to newer data. **Phase C:** improved model forms (Morse/long-range potentials, measured level
  corrections).
- **M1–M8:** milestones: M1 reproduce the published potentials; M2 hyperfine Hamiltonian; M3 intensities
  and spectra; M4 measurement data sets; M5 fitting; M6 uncertainty quantification; M7 first release;
  M8 web explorer.
- **Levers 1–7:** the ways the published model could be improved, set out in
  `research/beyond-hannover.md`.
