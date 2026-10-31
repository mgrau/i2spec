# Frozen release: October 30, 2026 — the final form with its uncertainty calibrated

Release 20261029 (the final form) with its uncertainty model refitted. The coefficients, frequencies and reference
predictions are those of release 20261029; only `expanded/uncertainty.json` differs. Reached by id:

```python
i2spec.ResearchModel("research-20261030-i2spec", isotopologue)
```

Release 20261021 stays the default.

## What changed from release 20261029

The uncertainty model keeps its form (σ_stat, σ_model and σ_HF). Release 20261029 fitted σ_model and σ_HF to the
cross-validation residuals of the training rows. This release fits them to those and to the residuals of the
calibrate data sets together (`data/training/roles_20261029`): sets the model never trained on, predicted as they
were in the third blind test.

- **Calibrate rows.** 2,432 rows. Left out: compilation rows whose source is also held, lines the training data hold
  inside a CIPM mean (`held_out_in_training` of `outcome.json`), record rows, diagnostic exclusions, and values
  printed again in a later set (each measurement once).
- **Correlations.** In the calibrate rows, each level is a random effect with the model's own per-level σ_model, so
  that lines sharing a level are correlated as the model assumes. Each laboratory family has a random offset of its
  absolute frequencies: measurement error, which is not part of a prediction's σ. In the intraline intervals, the
  pattern-scale error is one per line and the floor one per component, so a line counts once however many of its
  components were measured.
- **Two terms beyond the core.** Each is ≥ 0, so it can only add uncertainty. The B hyperfine floor grows as
  exp(g (v′ − 44)) up to the last validated v′ (`a_growth_per_v`). B σ_model beyond the core gains a (v − v0)² term
  (`quadratic_per_v2`), and does not fall below the core's per-level share at v′ 44. Earlier releases read both as 0.
- **Isotopologue floors.** Kept from release 20261029. The calibrate rows of the other isotopologues are covered
  within 2σ with or without them.

| term | this release | 20261029 |
|---|---|---|
| core σ_HF floor a | 63 kHz | 32 kHz |
| core σ_HF pattern-scale error κ | 2.1 × 10⁻⁴ | 7 × 10⁻¹² |
| core σ_model, lines near precision training rows | 0.162 MHz | 0.139 MHz |
| core σ_model, other core lines | 0.288 MHz | 0.213 MHz |
| B σ_model at v′ 44 / 60 / 70 / 77 (J = 50) | 0.8 / 4.6 / 39 / 285 MHz | 0.9 / 5.7 / 40 / 28 MHz |
| B σ_HF floor at v′ 44 / 60 / 68 | 0.063 / 1.8 / 9.4 MHz | 0.073 MHz at every v′ |

The laboratory offset of absolute frequencies in the core is 0.17 MHz.

## Coverage

| rows | within 1σ / 2σ, this release | 20261029 |
|---|---|---|
| calibrate rows, all | 80.2 % / 96.2 % | 60.1 % / 83.5 % |
| calibrate ¹²⁷I₂ rows in the core | 79.8 % / 95.2 % | 60.3 % / 86.2 % |
| calibrate hyperfine intervals in the core | 82 % / 96 % | 56 % / 84 % |
| calibrate hyperfine intervals, B v′ 50–79 | 46–85 % / 93–97 % | 12–60 % / 35–80 % |

Every region of the cross-validation is covered as well as in release 20261029 or better. Two groups fall below
90 % within 2σ: the held-out levels at B v′ 50–59 (88 %, as in release 20261029) and three intervals at B v′ 70–79,
where the model is at its limit. The 1978 atlas volumes (B v′ ≤ 62), on their own scale, fall within 2σ for 94–98 %
of their lines beyond the core. With a scatter of 35–51 MHz, they cannot set σ_model there.

## Verification

The references of release 20261029 are replayed by `ResearchModel(id).verify()` and reproduce as they did there.

## Locked test

The locked-test sets of `data/training/roles_20261029` (17 sets, kept out of every fit and of the uncertainty
refit), scored once with this release in `validation/locked_test`. The registry, with the criteria, strata and
rows left out in its protocol, was frozen and committed before any value was read; `results.json` is the
registry's score and `validation/locked_test_results.json` the summary below. Every set was also in the third
blind test of release 20261029 (same forecasts, the earlier 1-sigma), and its protocol says so.

z = r / √(σ_model² + σ²). Target: 95 % within 2σ; floor: 90 %; a stratum of fewer than 10 rows is reported,
not judged. Left out before scoring: 127I2 P(13) 43-0 (31 rows) and the lines a training CIPM mean contains (61 rows of
borde1981a, knox1971a and landsberg1976a); knox1971a's ten intervals that tesic1975a prints again with the same
values count once.

| stratum | rows | within 1σ / 2σ | rms z | result |
|---|---|---|---|---|
| all | 164 | 66 % / 88 % | 1.70 | below the floor |
| ¹²⁷I₂ absolute frequencies | 44 | 68 % / 98 % | 1.10 | target met |
| ¹²⁷I₂ intervals between lines | 2 | 100 % / 100 % | 0.67 | not judged |
| ¹²⁷I₂ hyperfine intervals | 86 | 57 % / 78 % | 2.18 | below the floor |
| other isotopologues | 32 | 84 % / 100 % | 0.63 | target met |

Under the 1-sigma of release 20261029 the same rows give 59 % / 85 % (all) and 44 % / 72 % (hyperfine intervals).

One set dominates the shortfall. cruz1994a (20 intervals of R(115) 20-1, a scanned Fabry–Perot scale) is within
2σ for 5 of them: its intervals are 0.88 % smaller than the model's, in proportion to their size (one scale factor
leaves an rms of 0.61 MHz, the 0.54–0.62 MHz of the paper's own fit). Its printed ΔeqQ, 1911.8(0.9) MHz, is
0.84 % below the model's 1928.1 MHz. A scale error of the interferometer or an error of the model's eqQ at B v′ 20,
J′ 116 would each give this; the set cannot tell them apart. The other 144 rows are within 2σ for 96.5 %;
the hyperfine intervals without cruz1994a for 94 % (62 of 66). These are stated beside the score, not instead of it.

The hyperfine constants of wakasugi1988a and the X eqQ of wakasugi1988 have no 1-sigma in the model: rms residual
3.5 MHz (B − X eqQ, rms r/σ 2.3), 6 kHz (C), 16 kHz (d), 5 kHz (δ) and 7.9 MHz (X eqQ, rms r/σ 1.7).
