# Frozen research model validation

The coefficient search and numerical acceptance checks are finished. The four
[research models](research-models.md) are explicitly selectable. Their shared
functions contain no independent line or level corrections. A genuinely
prospective experimental test is **awaiting new laboratory measurements**;
none were available when this release was prepared.

## What the current evidence supports

The balanced representative has 61 fitted molecular coefficients, 18 inherited
coefficients and four additional fitted atlas calibration terms. Its reused
absolute-component benchmark RMS is **156.36 kHz**, with **576.20 kHz** worst
repeat-level error and a signed mean error of **−90.45 kHz**. Its intraline
benchmark RMS is **58.56 kHz**. These are descriptive errors, not prediction
standard deviations.

All four models passed the complete seven-fold, two-start exact-J protocol,
physical/numerical acceptance and frozen-reference replay. The installed API
reproduces all 2,608 reference values within 1 Hz. Starts were selected by
training objectives and acceptance checks, never by their held-out errors.
The finite search is not a proof of a global optimum.

The seven folds overlap and were reused during architecture exploration.
Balanced B(v=39) whole-band absolute RMS is **3.44 MHz**, substantially worse
than its component benchmark. The 60 distinct benchmark components span only
34 lines and ten connected line/source groups. Counting components, repeats
or overlapping folds as independent calibration trials would overstate the
evidence.

The uncertainty audit freezes two **empirical discrepancy envelopes** around
the balanced point forecasts:

| Envelope | Half-width | Frozen recipe and interpretation |
|---|---:|---|
| Routine | 0.6 MHz | Round up the future cohort's worst repeat error plus its maximum accepted-start prediction spread to 0.1 MHz. |
| Stress | 13 MHz | Round up the largest repeat error plus the largest accepted-start spread across the future cohort and seven folds to 1 MHz. |

The full-precision inputs are 0.5761994123458862 + 0.008391499519348145 MHz
and 3.509594202041626 + 9.068304181098938 MHz respectively. Neither envelope
has an assigned coverage probability. The routine envelope covers all 60
balanced benchmark components but none of the four B(39) fold components.
The wide stress envelope incorporates observed optimizer sensitivity and
band/source extrapolation error; it does not establish MHz accuracy everywhere.
Fixed accuracy thresholds of **0.1, 0.3 and 1 MHz** are scored alongside both
envelopes so that a broad interval cannot disguise poor frequency accuracy.
Other models use the same declared comparison thresholds; these widths are
calibrated from balanced results only.

Saved hyperfine derivatives condition on the fitted potential. They omit the
potential/hyperfine/calibration cross-covariance and inherited-constant
uncertainty. They cannot supply a full prediction covariance. A nominal 95%
prediction interval remains unestablished; neither these empirical envelopes
nor collecting the finite target list automatically establishes that claim.

## Prospective measurement plan

The immutable registry is installed under
`src/i2spec/data/research/20261018/validation/prospective_v1/` in a checkout,
or `i2spec/data/research/20261018/validation/prospective_v1/` in the installed
package. Its manifest SHA256 is
`e0cdfbaba8fc4758a1700edfeceff9b92b0af76a56b737b081a9a56dfa9a32f4`.
It fixes all four models' 192 point forecasts and envelope endpoints, the
balanced primary model, target identities, exposure metadata and scoring
rules before any new laboratory outcomes. Existing registries cannot be
overwritten by the preparation command.

The plan requests a1, a3 and a10 on **16 lines in eight P/R pairs**:

| Region | P/R lines | Band |
|---|---|---|
| Routine | P(53)/R(51), P(61)/R(59) | 32–0 |
| B(39) stress | P(53)/R(51) | 39–0 |
| B(39) hot-band stress | P(79)/R(77) | 39–1 |
| High-J stress | P(159)/R(157), P(163)/R(161) | 25–0 |
| X(4) stress | P(63)/R(61) | 33–4 |
| Low-B stress | P(33)/R(31) | 1–0 |

The wavelengths span approximately **521–631 nm in vacuum**. The accompanying
`lab_targets.csv` lists each primary prediction and wavelength. Laser access,
signal strength, assignments, cell conditions and traceability still require
experimental review. Record infeasible targets and the reason before seeing
outcomes; revise a separately versioned plan if necessary. Do not silently
drop difficult or inaccurate outcomes.

All 192 prospective frequencies were also checked with a halved radial mesh
and, separately, expanded radial boxes. The 384 component comparisons passed
the 1 kHz tolerance: the largest mesh change was **18.72 Hz** and the largest
box change **0.84 Hz**. Wavefunction normalization and angular assignments
passed, and frozen coefficients/forecasts remained byte-identical. This checks
numerical stability; it does not establish experimental accuracy or observability.

The two lines in each pair share an upper rovibrational J. Equal component
rank does not guarantee a common upper hyperfine eigenstate; their component
frequency difference does not exactly cancel hyperfine structure. Fourteen
of these lines have no resolved measurements in the checked precision CSVs;
six have no assigned entries in either those CSVs or the checked atlas CSVs.
That audit does not prove the identities were never exposed elsewhere.
Fresh measurement independence and unseen-transition prediction are distinct.

The exact `127I2 P(13) 43-0` exclusion applies to fitting, calibration,
prospective targets and headline scoring. The source archive remains intact.

## Record and score new measurements

From a checkout:

```sh
I2SPEC_VALIDATION=src/i2spec/data/research/20261018/validation
uv run python -m i2spec.validation verify "$I2SPEC_VALIDATION/prospective_v1"
cp "$I2SPEC_VALIDATION/measurements_template.json" /tmp/iodine-measurements.json
# Fill in the copied measurement rows when new lab results exist.
uv run python -m i2spec.validation score "$I2SPEC_VALIDATION/prospective_v1" \
  /tmp/iodine-measurements.json --out /tmp/iodine-validation-results.json
```

For an installed package, locate the same directory with
`importlib.resources.files("i2spec").joinpath("data/research/20261018/validation")`.
The original templates and frozen registries should remain unchanged.

Each comparable row needs a source ID/link, conditions, a source-supported
`frequency_basis: "zero_pressure_unperturbed"`, a timezone-aware
`acquired_utc` after freezing, and explicit
`independent_of_fitting_archive: true` and
`outcome_unseen_before_freeze: true` declarations. These are provenance
declarations to audit against laboratory records, not facts software can
verify independently. Keep finite-pressure readings and all correction
evidence in the source records; do not estimate a correction from the model
residual. A missing value remains pending, and a value whose convention is
unknown remains noncomparable.

Record `standard_uncertainty_MHz` if the source gives a standard uncertainty.
For an expanded uncertainty use `expanded_uncertainty_MHz` and the explicit
`coverage_factor`. Without that factor, the program does not infer a standard
uncertainty. Residuals divided by measurement uncertainty describe disagreement
with that measurement alone; they are not model-normalized confidence scores.

The report preserves each measured residual, scores unique physical components
against median repeated outcomes, retains the worst repeat disagreement and
reports line/source counts, connected groups and each declared stratum.
Derived intraline differences use the lowest measured component rank on each
line as reference; their errors are reported separately without assuming
independent measurements or propagating an unknown covariance. New outcomes
never change model coefficients, point forecasts, envelope widths or the
primary model. The output is a separate file and cannot overwrite an existing
report. Completion of the measurement list requires scientific review of
errors, conditions and intended tolerances before a final adequacy decision.

To prepare a different experiment, copy `plan_v1.json`, declare the targets
and primary model before any outcomes, and use:

```sh
uv run python -m i2spec.validation freeze /tmp/new-iodine-plan.json \
  /tmp/new-iodine-forecast-registry
```

## Published external candidates

The byte-preserved `literature_v1/` registry predates this package integration.
It contains predictions for the 2025
[R(158) 25–0 study](https://doi.org/10.35848/1347-4065/ae0225) and the
[PGMF R(56) 32–0 study](https://doi.org/10.1088/2631-8695/ada487).
The exact first outcome was unavailable. The second appeared in a search
snippet before the point-registry freeze, after coefficients were frozen,
and is explicitly retrospective. Primary methods establishing correction and
uncertainty conventions were unavailable for both, so **zero literature
outcomes are currently scored as comparable**. The old R(158) unresolved
center was already represented in fitted atlases. These two cases cannot
replace prospective lab validation or establish probability coverage.

The portable `uncertainty_summary.json`, exact
`envelope_calibration_inputs.json`, pending reports and readiness provenance
preserve the evidence without requiring the ignored fitting workspace.

## Blind comparison on the pending data sets

On October 20, 2026 Matt decided that the 36 data sets transcribed that day,
which carry `exclude = "pending: ..."`, stay out of every fit and serve as a
blind test. The procedure:

1. The balanced-expanded release was frozen.
2. `prototypes/explore/physical_search/blind_pending.py freeze` wrote the
   registry `src/i2spec/data/research/20261020/validation/pending_blind_v1/`
   at 2026-10-20T20:24:59Z. Its targets come from the identity columns of
   each `data.csv` only: line, component, kind, reference line and component.
   The value and uncertainty columns were never parsed.
3. Four models then predicted every target:
   - balanced-expanded (primary);
   - research-balanced, within its domain;
   - i2spec2026o, with its measured B hyperfine table;
   - the Hannover 2008 potentials with the bare BKT02/S06 formulae.
4. `score` read the values only afterwards. It first checks the registry and
   data-file hashes, and nothing was refitted or reselected.

The exact P(13) 43-0 line is a diagnostic target only.

The 580 targets fall into three groups, classified from the repository's
`[overlaps]` declarations before scoring:

- **365 blind.**
- **162 seen.** The same measurement, or a BIPM mean containing it, is already
  in training.
- **53 superseded.** Training holds the same quantity re-measured or
  re-evaluated by the same authors: Sansonetti 1990 against 1997, and Ikeda
  2020 against 2022.

The research models predict every 127I2 target; the 23 129I2 rows are outside
their domain.

Median |observed − predicted| on blind 127I2 targets:

| Blind targets | n | Balanced-expanded | i2spec2026o | Hannover 2008 |
|---|---:|---:|---:|---:|
| Absolute, line measured before by others | 78 | 124 kHz | 17 kHz | 1.06 MHz |
| Absolute, line never measured before | 25 | 621 kHz | 633 kHz | 1.20 MHz |
| Intraline, line measured before | 182 | 61 kHz | 9 kHz | 33 kHz |
| Intraline, new line (Weller v′51, Picard 1996) | 40 | 248 kHz | 20 kHz | 196 kHz |
| Inter-line, line measured before | 5 | 281 kHz | 144 kHz | 1.34 MHz |
| Inter-line, new lines (Bertinetto 1987, v′9) | 17 | 950 kHz | 4.51 MHz | 4.70 MHz |
| Beyond balanced's domain (Weller v′51, Bordé v′58/62 intervals) | 48 | 728 kHz | 8.5 kHz | 2.47 MHz |

Inside balanced's domain, balanced-expanded and balanced give identical
results. On lines nobody had measured before, the shared-function model ties
i2spec2026o on absolute frequencies and is about five times better on
inter-line intervals. Other examples, rms:

| Set | Balanced | i2spec2026o | Hannover |
|---|---:|---:|---:|
| de Melo 2024, R(158) 25-0 against its reference | 0.59 MHz | 18.7 MHz | 18.3 MHz |
| Powel 2021, NIR | 0.08 MHz | 2.2 MHz | 5.7 MHz |
| Shiner 1993 | 0.20 MHz | 0.50 MHz | 0.88 MHz |
| Ikeda 2020, R(73) 46-0 (superseded) | 0.30 MHz | 0.96 MHz | 1.3 GHz |

Where other measurements of the same line are in training, i2spec2026o is
5–10 times better. Its per-level corrections and measured B hyperfine table
interpolate those measurements; balanced's core benchmark is 156 kHz.

Beyond the core, the only blind data are B hyperfine intervals, and there the
six-parameter physical model loses clearly:

- Weller 2026, R(41) 51-0: 0.77 MHz rms. The errors grow with the interval
  size, consistent with eqQ about 0.5 MHz too small.
- Bordé 1981, R(26) 62-0 and R(98) 58-1: up to 2.5 MHz.

The nested cross-validation had suggested about 250 kHz at v′44–53. Its
interval-by-interval residuals, extracted on October 21, already showed 1.6–2.0 MHz
at v′58 (BIPM 2005, R(98) 58-1) and v′62 (Goncharov 2007). Its held-out Chen 2004
eqQ constants were off by 0.4–1.5 MHz at v′42–53, in the same direction as Weller's
residuals. The cross-validation had therefore seen this failure. i2spec2026o
interpolates Chen 2004's constants measured at the same v′: five J′ values at
v′51, and v′61 and 63 around v′62. Its table also holds a per-line fit of
BIPM 2005's R(98) 58-1. So its 10–40 kHz there is near-interpolation rather
than an independent prediction. The Hannover formulae freeze above v′53 and
miss by 10–80 MHz.

Every model misses five absolute rows by 12–113 MHz, but these carry 30–60 MHz
uncertainties: Baird 1979 P(62) 17-1, and Zhang 1999 R(60) 8-4, P(119) 9-4,
P(32) 6-3 and R(38) 6-3. They dominate any MHz rms, so medians and r/σ are
reported. As the protocol requires, they stay in the scores.

These outcomes have now been seen. Any model change made in response is no
longer blind for these sets. `results.json` beside the registry holds every
residual and the pre-declared summaries.

### Coverage of the uncertainty model (post hoc)

The uncertainty model of balanced-expanded (October 21; see
[research models](research-models.md#uncertainties-of-balanced-expanded-october-21-2026))
was built from the training fit and its cross-validation only. It was then frozen
(`uncertainty.json`, SHA256 `b5618c43…`) before it was scored on these rows. The
scores are **post hoc** all the same: the accuracy results above had been read the
day before, so this is a calibration check, not a blind test.

`python prototypes/explore/physical_search/expanded_uncertainty.py coverage` first
checks that the model still reproduces the frozen forecasts; it does, exactly. It
then writes `coverage_posthoc.json` beside the registry. With z = residual /
√(σ_model² + σ_meas²):

| Rows | n | within 1σ | within 2σ | median σ_model |
|---|---:|---:|---:|---:|
| blind | 347 | 66% | 91% | 0.058 MHz |
| blind, absolute | 103 | 66% | 88% | 0.34 MHz |
| blind, intraline, core | 174 | 66% | 91% | 0.053 MHz |
| blind, intraline, beyond (Weller, Bordé) | 48 | 52% | 94% | 0.91 MHz |
| blind, between lines (Bertinetto) | 22 | 95% | 100% | 1.35 MHz |
| seen | 126 | 70% | 95% | 0.053 MHz |
| superseded | 53 | 57% | 100% | 0.18 MHz |

Nominal values are 68% and 95%.

- **Beyond the core.** The B hyperfine scale error, fitted to the held-out
  intervals and constants, gives Weller's R(41) 51-0 a 1.1 MHz σ on its widest
  interval. Weller is within 2σ on 19 of 20 rows.
- **Outlier set.** Nishiyama 2013 is the one clear outlier: a median miss of
  0.64 MHz against a median σ of 0.16 MHz, with 4 of 15 rows within 2σ.
- **Over-coverage.** The between-line intervals are over-covered.

Nothing was tuned on these rows. A blind test of the uncertainties needs
measurements made after the model was frozen, scored through a prospective
registry (above).

### Second blind test: 35 sets of October 21, with frozen uncertainties

On October 21 Matt decided that all new data stay a blind test for balanced-expanded ("let's keep all new data
as a blind test for now for this model"). The 35 sets transcribed that day went into a second registry,
`pending_blind_v2`, frozen with `blind_pending.py freeze … --after pending_blind_v1` before any of their values
was read. This time the primary forecasts also carried the model's frozen 1σ (`uncertainty.json`, SHA256
`b5618c43…`), so the coverage below is a genuinely blind test of the uncertainties.

The registry has 698 targets: 474 blind, 172 seen through training, 51 tied by `[overlaps]` to v1 sets whose
values had been read (`scored_before`), and 1 superseded. Balanced-expanded predicts the 502 ¹²⁷I₂ targets;
the other 196 are other isotopologues. The snapshot checks `data.csv` only and records `meta.toml` for
information. It was scored once, and no meta had changed.

Median |observed − predicted| on blind targets that all three models predict:

| Blind targets | n | Balanced-expanded | i2spec2026o | Hannover 2008 |
|---|---:|---:|---:|---:|
| Absolute frequencies | 70 | 0.50 MHz | 0.12 MHz | 17 MHz |
| Intraline intervals | 199 | 76 kHz | 26 kHz | 113 kHz |
| Inter-line intervals | 76 | 110 kHz | 138 kHz | 847 kHz |

Coverage of the frozen σ (z = residual / √(σ_model² + σ_meas²)):

| Blind rows | n | within 1σ | within 2σ |
|---|---:|---:|---:|
| all | 345 | 61% | 90% |
| absolute | 70 | 61% | 100% |
| intraline | 199 | 49% | 83% |
| inter-line | 76 | 93% | 99% |

Nominal values are 68% and 95%. The absolute and inter-line uncertainties hold; the intraline (hyperfine)
uncertainty is too small. Kremser 1994 is the one set far outside: 10 rows at rms z = 11.
`results.json` in the registry holds every row.
