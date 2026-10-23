# Training inputs of release 20261021

The inputs the fits of the released model read, frozen so that `i2spec.fit` can rebuild the release from the
repository. They are derived from the data sets in `data/` (see `i2spec.fit.reduction`), not new measurements.

| file | what | read by |
|---|---|---|
| `rows.csv.gz` | 47 045 rows: every hyperfine-free line center, interval between lines and X level energy of the fits, in MHz, with its sigma, category and term-value functional (`levels`). Columns `core` and `continuation` give each row's role in that fit (`train:<n>` or `benchmark:<n>`, n its position). | `i2spec.fit.rows.release` |
| `hyperfine.csv.gz` | 1 426 intervals between two components of one line (1 151 in the core), Chen et al. 2004's effective B constants and Pique et al. 1986's | `i2spec.fit.hyperfine_rows.release` |
| `core_shape_nodes.json` | extra radial nodes of the core fit's single-well constraint | `i2spec.fit.core.shape_nodes` |
| `hyperfine_potential.json` | the continuation's tail values on which the hyperfine beyond the core was fitted (before the out-of-core rows were re-reduced with that hyperfine) | `i2spec.fit.hyperfine_beyond.potential` |
| `cross_validation.json.gz` | outputs of the cross-validation refits the uncertainty model is built from: the continuation's four v folds, the core fit's folds, the nested hyperfine folds, the benchmark predictions | `i2spec.fit.uncertainty.cross_validation` |

The rows were built when the core fit was made. Rebuilt from today's data files with the model's own hyperfine
(`i2spec.fit.reduction.release_rows`), every row returns, with no value or sigma moved by more than 45 Hz. The core
rows carry the centers of the previous round of the alternation between potentials and hyperfine (50 Hz away), and
the benchmark rows' measurement sigmas in `cross_validation.json.gz` are those the uncertainty build used, from the
reduction before the model's own hyperfine.

The exact line 127I2 P(13) 43-0 is left out of every fit and calibration.
