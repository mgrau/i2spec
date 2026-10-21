# Prospective forecast numerical and assignment feasibility audit

Completed: 2026-10-18T12:27:44.376208+00:00

**All 384 component comparisons passed the unchanged 1 kHz numerical threshold.** The frozen registry, molecular coefficients, source assets and evaluation implementation remained hash-identical. No measurement, fit, correction or forecast revision was performed.

## Scope and controls

The audit read `src/i2spec/data/research/20261018/validation/prospective_v1` and its four frozen model records: 192 point forecasts, 48 component targets and 16 rovibrational lines. Each fresh `ResearchModel` received a transient deep copy of its specification with only these numerical grids changed:

| Check | X grid | B grid |
|---|---|---|
| Mesh refinement | Original 2.1–4.2 A box, h=0.005 A | Original 2.35–6.1 A box, h=0.005 A |
| Box expansion | 2.0–4.4 A, original h=0.01 A | 2.25–6.3 A, original h=0.01 A |

The B-spline order, number of computed levels, molecular coefficients, hyperfine functions, component-label algorithm and constants remained fixed. These are two convergence perturbations of the B-spline representation, rather than a different solver family or an analytic error bound. Two single-thread subprocesses completed eight finite jobs, totaling 128 line evaluations and 384 component comparisons. The script never writes into the frozen registry or packaged assets.

## Results

| Requirement | Result |
|---|---:|
| Maximum finer-mesh frequency change | 18.715858 Hz |
| Maximum expanded-box frequency change | 0.834465 Hz |
| Largest radial wavefunction normalization error | 1.055e-14 |
| Component comparisons within 1 kHz | 384 / 384 |
| Consistent angular assignments across variants | 192 / 192 |
| Required a1/a3/a10 frequencies and strengths finite and positive | All |
| Target and J +/- 2 neighboring levels below asymptote | All; smallest margin 669.302 cm^-1 |
| Frozen input files unchanged | Yes |
| Numerical/assignment failures | None |

The largest mesh difference occurs for absolute_accuracy, 127I2 R(51) 39-0, a1. The largest box difference occurs for balanced, 127I2 P(163) 25-0, a1. Full per-role, per-variant and per-target values, bound-state margins, angular labels and failures are retained in `jobs/*.json`, `summary.json` and `component_checks.csv`.

## Wavelength and conditional hyperfine-strength feasibility

Vacuum component wavelengths span approximately **520.810–631.423 nm** across all four models. The chosen components contain **4.006–5.494%** of their line’s total normalized hyperfine strength, or **72.9–100%** of its strongest modeled component. These are conditional fractions within a single rovibrational line. They do not assess the electronic/vibrational transition moment, thermal populations, absolute absorption, neighboring rovibrational blends, available lasers or laboratory sensitivity.

The following table uses the primary balanced model. Wavelengths are derived from the original frozen frequencies, while strengths come from the finer-mesh calculation. Each range covers the planned a1/a3/a10 components. `balanced_lab_feasibility.csv` preserves the individual component values and their nearest other modeled hyperfine component separation. A separation alone does not establish experimental resolvability.

| Target line | Vacuum wavelength range, nm | Relative strength vs strongest component |
|---|---:|---:|
| P(53) 32-0 | 532.242130–532.242580 | 0.827–1.000 |
| R(51) 32-0 | 532.020409–532.020859 | 0.824–1.000 |
| P(61) 32-0 | 532.646314–532.646767 | 0.848–1.000 |
| R(59) 32-0 | 532.390492–532.390945 | 0.846–1.000 |
| P(53) 39-0 | 521.022319–521.022766 | 0.827–1.000 |
| R(51) 39-0 | 520.809845–520.810293 | 0.824–1.000 |
| P(79) 39-1 | 528.557822–528.558304 | 0.881–1.000 |
| R(77) 39-1 | 528.232214–528.232696 | 0.879–1.000 |
| P(159) 25-0 | 556.398421–556.398944 | 0.939–1.000 |
| R(157) 25-0 | 555.671553–555.672074 | 0.939–1.000 |
| P(163) 25-0 | 556.949139–556.949666 | 0.940–1.000 |
| R(161) 25-0 | 556.202714–556.203239 | 0.940–1.000 |
| P(63) 33-4 | 555.858531–555.859028 | 0.853–1.000 |
| R(61) 33-4 | 555.574314–555.574811 | 0.851–1.000 |
| P(33) 1-0 | 631.422299–631.422925 | 0.736–1.000 |
| R(31) 1-0 | 631.229022–631.229650 | 0.729–1.000 |

All required components retain the same F/I/upper-state-index/lower-state-index assignments between the two numerical variants. This establishes formal consistency of the modeled labels; independent experimental line/component assignment is still required. P/R partners share an upper rovibrational J, but the same a-rank does not guarantee the same upper hyperfine eigenstate or exact cancellation of its hyperfine energy.

## Interpretation and reproducibility

The finite prospective target plan is numerically feasible at the requested accuracy under both checks. This result does not measure predictive error, certify laboratory access or calibrate probability coverage. The new-measurement validation gate remains pending external measurements. The registered envelopes retain their existing nonprobabilistic meaning.

Reproduce with `OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python prototypes/explore/physical_search/prospective_numerical_audit.py --workers 2`. Finished jobs are retained; use a distinct output directory for a new audit version. Five targeted tests verify transient-only grid changes, frozen-point differences, normalized relative strengths, missing-label failures and nonpositive-strength failures. Input SHA256 hashes and implementation/environment provenance are in `manifest.json`.
