"""The fit that produces the model's coefficients from the measurements.

Stages, in order (each reads the previous stage's output; release.rebuild runs them all):

  rows, reduction     the training table: every measured line center, interval between lines and level energy, as a
                      linear functional of term values; built from the data files with the model's own hyperfine
  core, core_hyperfine
                      the potentials of X v''<=17 and B v'<=43, J<=200, and the B hyperfine functions, alternated
  continuation        tail terms past the core's edge: X to v''=108, B to v'=80, the core frozen
  hyperfine_beyond, second_order, x_hyperfine, b1g_resonance
                      the hyperfine beyond the core: B background + second-order sums, X eqQ to the separated atoms,
                      the B~1g(1Pig) resonance at v'77, 78
  uncertainty         the covariance of the fit, sigma_model and sigma_HF from cross-validation
  isotopologues       129I2 and 127I129I: the isotope shift
  package             levels, hyperfine tables and reference predictions of a release
  registry            blind tests: forecasts frozen before the data are read, scored once

The release's own inputs (training rows, hyperfine records, cross-validation outputs) are frozen in
data/training/<release>/; from them every stage returns the released numbers (tests/test_fit*.py; the full refits
run with pytest --slow).
"""
