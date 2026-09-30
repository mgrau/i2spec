# Research scripts

The scripts behind each fit, parameter set and figure in `docs/`. They are records of how a number
was produced, not library code: the package in `src/i2spec/` never imports them. Run them from the
repository root, e.g. `uv run python prototypes/orsay_fit.py`; the ones that draw figures or read the
APO atlas's HDF5 file need `uv run --group research ...`. Outputs go to `prototypes/out/` (not
committed); third-party inputs live in `data/external/` (not committed; see its READMEs).

## Level corrections and the NIR (parameter sets 2026b–f)
| script | what |
|---|---|
| `level_corrections_fit.py` | level corrections δE(v, J) from every comb-referenced ¹²⁷I₂ line |
| `nir_corrections.py` | the first, NIR-only version of those corrections (2026a) |
| `nir_band_fit.py` | the line-level correction of the v′ = 0 → v″ = 12–17 bands |
| `bodermann1998c_make.py` | transcribes Bodermann's thesis tables into `data/observations/bodermann1998c` |
| `bo_refit.py` | the B-state Born–Oppenheimer constant, from the isotopologue lines |
| `error_map.py` | what the default model achieves, by region |

## Potentials: MLR curves (Phase C, 2026d and 2026h)
| script | what |
|---|---|
| `mlr_x.py`, `mlr_x_fit.py`, `mlr_x_refit.py` | MLR X trials: against the published curve, to levels, to data |
| `mlr_x_levels.py`, `mlr_b_levels.py` | the alternating X/B level fits that produced `mlr_*_2026c` |
| `mlr_b_fit.py` | the first MLR B fit |
| `mlr_joint_fit.py`, `mlr_joint_uncertainty.py` | X and B fitted together, and its covariance |
| `mlr_from_rkr.py`, `global_fit.py`, `global_fit_evaluate.py` | a from-scratch direct-potential fit |
| `mlr_evaluate.py` | any MLR pair against every observation |
| `x_bridge.py` | bridging the unmeasured X gap v″ = 18–47 |
| `mlr_x_atlas.py`, `mlr_x_atlas_figure.py` | the X refit to the Orsay-atlas levels (`mlr_x_2026d`) and its figure |
| `near_dissociation_lines.py` | the line list to the B dissociation limit |
| `hannover_baseline.py` | rebuilding the Hannover 2004 model (feasibility) |

## Hyperfine
| script | what |
|---|---|
| `hyperfine_stage2.py` | the Stage 2 refit (`docs/design/hyperfine-fit.md`) |
| `hyperfine_chen2004.py` | Chen 2004's B-state parameters against the formulae |
| `hfs_measured_table.py` | builds the measured B-state table `i2spec.hfs_table` ships |
| `hfs_table_validate.py` | leave-one-out test of that table |

## The Orsay atlas, 11 000–14 000 cm⁻¹ (`docs/research/orsay-atlas-11000-14000.md`)
| script | what |
|---|---|
| `orsay_ocr.py`, `orsay_put.py`, `orsay_check.py` | photographs → typed pages → checked TSVs |
| `orsay4_crops.py`, `orsay4_reconcile.py`, `orsay4_dataset.py` | Partie IV (19 700–20 035 cm⁻¹): crops for reading, plate against classification against the Dunham constants, the line lists |
| `orsay_assign.py` | assignment against the model list thinned to the atlas density |
| `orsay_calibration.py` | the wavenumber scale of each part |
| `orsay_fit.py` | X v″ = 18–25 levels from the assigned lines (2026f) |
| `orsay_part1_test.py` | registration of part I against the model, window by window |
| `orsay_other_volumes.py` | which levels the other volumes of the series would measure |
| `orsay_levels_figure.py` | `docs/figures/orsay_x_levels.png` |

## FTS atlases, spectra and cross sections (`docs/research/spectra-validation.md`, `atlas-line-positions.md`)
| script | what |
|---|---|
| `atlas_lines.py`, `atlas_dataset.py` | line positions from the Salami & Ross and APO/NIST spectra |
| `compare_salami_ross.py`, `scan_salami_ross.py`, `band_shifts_salami_ross.py` | the model against the Salami & Ross atlas |
| `plot_salami_ross.py`, `plot_scan_salami_ross.py` | their figures |
| `compare_cross_sections.py`, `check_continuum.py` | absolute cross sections and the continuum |
| `nolleke_assign.py`, `nolleke_combdiff.py`, `nolleke_transmission.py` | the Nölleke 2018 915–985 nm spectra (no registration; documented negative) |

## IodineSpec5 comparison (`docs/research/iodinespec5.md`)
IodineSpec5 and its output are not part of the repository. These read the output files of a local copy
of the program (`$IODINESPEC5_OUTPUT`, or `data/external/iodinespec5/output/`).

| script | what |
|---|---|
| `iodinespec5_compare.py` | its output against i2spec and the data |
| `iodinespec5_gap.py`, `iodinespec5_s06.py` | its Dunham levels above v″ = 17; its hyperfine formulae |
