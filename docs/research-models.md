# Selectable frozen research models

The October 18, 2026 release preserves four exact-J research representatives.
The October 20, 2026 release adds a fifth model,
[balanced, expanded](#balanced-expanded-october-20-2026), which continues the
balanced model to near dissociation. Select a model explicitly through
`i2spec.research.ResearchModel` or the `research` CLI.

Since October 21, 2026 balanced-expanded is the preferred model (Matt: "Since it
is better on lines it wasn't fit to, we will prefer it"). For ¹²⁷I₂ lines inside
its domain, the following give its position, hyperfine components and
[uncertainty](#uncertainties-of-balanced-expanded-october-21-2026):

- the line lookup;
- the terminal browser;
- the command line;
- the explorer.

i2spec2026o remains the parameter set for line strengths, the other
isotopologues and the lines beyond that domain; each line says which model gave
it. See [preferred in the apps](#preferred-in-the-apps).

| Canonical model ID | Role alias | Fitted molecular | Inherited molecular | Fitted calibration | Absolute RMS kHz | Intraline RMS kHz |
|---|---|---:|---:|---:|---:|---:|
| `research-20261018-compact` | `compact_value` | 56 | 18 | 4 | 190.737 | 85.879 |
| `research-20261018-absolute` | `absolute_accuracy` | 61 | 11 | 4 | 155.537 | 108.003 |
| `research-20261018-balanced` | `balanced` | 61 | 18 | 4 | 156.363 | 58.557 |
| `research-20261018-precision` | `interval_precision` | 68 | 11 | 4 | 196.106 | 29.752 |
| `research-20261020-balanced-expanded` | `balanced_expanded` | 105 | 73 | 5 | 156.347 | 58.557 |

The October 18 scores were frozen with that release. The expanded model's
scores were computed on October 20 against the current data files, which give
156.347 kHz for balanced itself, so the core is unchanged.

Each model uses shared MLR potentials, shared rotational functions and fitted
radial hyperfine functions. There are no independent line shifts or
vibrational-level offsets. The numerical grids and all coefficients remain
frozen. Fixed-X models inherit the eleven published X hyperfine coefficients;
anchored-X models inherit the two ground-state anchors and two spin-spin
constants. The dispersion constants are counted separately in the inventory.

## Command line

```sh
uv run i2spec research list
uv run i2spec research list --json
uv run i2spec research predict 'P(53) 32-0' --model research-20261018-balanced --component a1
uv run i2spec research predict 'P(53) 32-0' --model absolute_accuracy --json
uv run i2spec research verify --model research-20261018-precision
```

`predict` and `verify` require `--model`. Omitting `--component` predicts the
hyperfine-free center. JSON includes the selected model, inventory and domain.
The exact `127I2 P(13) 43-0` line remains excluded from fitting and headline
scores; an explicitly requested prediction is labeled diagnostic.

## Python

```python
from i2spec import ResearchModel, research_models

for info in research_models():
    print(info.id, info.inventory, info.domain)

m = ResearchModel("research-20261018-balanced")
prediction = m.predict("P(53) 32-0", "a1")
print(prediction.frequency_MHz, prediction.wavenumber_cm1)

nu = m.transition(32, 0, 53, "P")  # hyperfine-free center, cm^-1
nu0, components = m.hyperfine_components(32, 0, 53, "P")  # center and offsets, MHz
print(m.verify())  # 60 absolute and 592 intraline frozen reference predictions
```

`frequency(line, component=None)` returns MHz. A line can be a label or the
tuple `(branch, J_lower, v_upper, v_lower)`. `transition` and
`hyperfine_components` follow `RovibronicModel`'s argument order and units.
Hyperfine predictions retain the fitted J±2 mixing (`dJ=2`); they apply no
measured correction table. Component objects returned to callers are copies.
The `coefficients` property returns a copy for inspection. Model metadata and
coefficient inventories are immutable dataclasses.

The domain is **127I2, X v≤17, B v≤43, J_lower≤200**, with P J_lower≥1.
Out-of-domain quantum numbers, unsupported isotopes and unknown model IDs
raise `ValueError`; no production or extrapolation fallback is substituted.
This API supplies frequencies and hyperfine patterns, not a line-intensity
catalog.

## Validation and limits

The absolute benchmark contains 60 distinct components from 34 transitions;
the intraline benchmark contains 592 intervals. The final architectures have
complete seven-fold, two-start exact-J validation. Coefficient starts are
selected by training objectives and physical/numerical acceptance. These
benchmarks and overlapping folds were reused during architecture exploration.
They are exploratory evidence, not fresh blind or nested validation, and the
finite search does not certify a global optimum. Whole-band prediction,
particularly B(v=39), remains weaker than the headline component benchmark.

Benchmark RMS values are not per-prediction standard deviations. The four
October 18 models carry no per-prediction uncertainty. Balanced-expanded has its
own, described [below](#uncertainties-of-balanced-expanded-october-21-2026). Calibration, prospective prediction and uncertainty studies are
separate evidence; these frozen coefficients are not refitted during loading
or prediction.

The [frozen validation workflow](research-model-validation.md) includes the
uncertainty audit, 48 prospective component targets, immutable predictions and
measurement-scoring commands. New laboratory measurements are still pending;
computational completion does not establish blind experimental validation.

All seven source JSON artifacts per role are installed under
`i2spec/data/research/20261018/`: coefficients/configuration, reference
predictions, selection, audits, source provenance, verification and verification
provenance. The release manifest preserves SHA256 hashes and the original
inventories. Loading checks coefficient hashes; `verify()` also checks reference
hashes and reproduces all 652 saved predictions to a 100 Hz threshold. Full
integrated replay differs by under 1 Hz for every role. Replay establishes
numerical reproduction, not independent predictive accuracy. Runtime evaluation
requires only the installed `i2spec` package and its NumPy/SciPy dependencies;
the ignored research workspace and observation archive are unnecessary.

## Balanced, expanded (October 20, 2026)

`research-20261020-balanced-expanded` (role `balanced_expanded`) continues
research-balanced to **X v≤108 and B v≤80**, with J_lower≤200. It keeps every
balanced coefficient. They are read at load time from the October 18 release,
which must match the SHA256 recorded in the new release. It has no fallback to
i2spec2026o and no line or level terms.

```python
from i2spec import ResearchModel

m = ResearchModel("balanced_expanded")
m.frequency("R(41) 51-0", "a1")       # MHz
m.frequency("R(20) 0-100")            # hyperfine-free center, X v''=100
m.verify()                            # 842 centers, 2,586 components, 1,434 intervals
```

**Potentials.** Past balanced's data, each frozen MLR exponent and rotational
function is continued by its own Taylor polynomial plus Bernstein tail terms:

- X: 12 outer exponent terms, 1 inner and 1 rotational;
- B: 10 outer exponent terms and 3 rotational, with both knots at 4.45 Å.

X also gains a change of C6 and nine cubic B-splines on 4.2–6.6 Å. One change of
the dissociation energy is shared by X and B, so their limits keep the atomic
splitting. An Orsay IV offset is the only new calibration term.

Levels are solved on the final fit's own grids, one per range of v. The highest
levels sit in a 40.6 Å box. Quasi-bound levels are kept when at least 90% of
their density lies inside the centrifugal barrier. The packaged levels
reproduce every training row of the fit to within 1.4 kHz.

**Hyperfine.** Inside balanced's domain (B v′≤43, X v″≤17) the model uses
balanced's own constants. Beyond it, constants come from tables:

- **B**: a background anchored to balanced's core predictions, plus
  second-order sums over the nine 2P3/2+2P1/2 states (the framework of Chen et
  al. 2005, with de Jong 1997 and Saute 1982 curves). The 1u mixing angle and
  the eqQ continuation were fitted. The 1g(1Πg) long-range scale was chosen by
  training cost.
- **X**: BKT02 eqQ continued to the separated-atom value, using Bacis et al.
  1980. The bounded C_X rests on one measured value and is unvalidated.
- **B v′77–78**: the near-resonant 1g(1Πg) level (Pique et al. 1986) is
  diagonalized with B. Its weak extra components appear among the line's
  components.

The tables cover levels bound below the dissociation asymptote: up to J′≈180
at v′44, falling to J′≈25 at v′80. A quasi-bound level has a center but no
hyperfine components. Asking for one raises `ValueError`.

**Inventory.**

- 105 fitted molecular coefficients:
  - 61 from balanced;
  - 38 new potential terms;
  - 6 new hyperfine parameters: the 1u angle, the eqQ curvature, the 1g(1Πg)
    scale, two X eqQ shape parameters and the X C slope.
- 73 inherited: balanced's 18, plus 55 atomic, curve and Pique constants.
- 5 calibration terms.
- 24 hyperfine coefficients that are derived rather than fitted: the B
  background (14), anchored to balanced's core predictions, and the X eqQ and
  C polynomials (5 each), anchored to BKT02 in the core. `model.json` lists
  them under `inventory_detail`.

**Evidence.** The release's `validation.json` holds:

- the four-fold cross-validation by v;
- the nested hyperfine cross-validation;
- the held-out comb benchmark;
- a same-rows comparison with i2spec2026o and the Hannover 2008 potentials;
- the isotopologue check.

On those rows, the expanded model scores as follows (rms):

| Rows | Balanced-expanded | i2spec2026o | Hannover 2008 |
|---|---:|---:|---:|
| Core comb, held out (79) | 0.26 MHz | 0.13 MHz (fitted to them) | 3.4 MHz |
| X v″26–47, Martin | 43 MHz | 143 MHz | ~400 GHz |
| X v″61–89, Martin | 104 MHz | 156 MHz | ~5 THz |
| X v″90–108, Martin | 332 MHz | 2.2 THz | 3.5 THz |
| B v′51–62, APO atlas | 37 MHz | 611 MHz | 64 GHz |
| B v′51–77, Orsay IV | 62–124 MHz | 62–121 MHz | 20–30 GHz |
| B v′78+, Orsay IV | 193 MHz | 137 MHz | 31 GHz |

The 36 pending data sets were then scored blind: see
[the blind comparison](research-model-validation.md#blind-comparison-on-the-pending-data-sets).
The headline results:

- Inside balanced's domain the expanded model is balanced.
- Outside it, the only blind measurements are B hyperfine intervals at v′51,
  58 and 62. There the physical hyperfine model misses by 0.1–2.5 MHz, where
  i2spec2026o's measured-constant table is within about 40 kHz.

**Limits.**

- 127I2 only, J_lower≤200.
- D0 and C6 are not separately determined: their fold-to-fold spread is
  0.012–0.016 cm⁻¹ and 2.5–4.2%.
- X v″≥90 is held out at 0.66 GHz rms.
- Out-of-core hyperfine is good to 0.1–2.5 MHz, not kHz.
- The v′77–78 resonance inherits Pique's ±5–10 MHz rotational constants. That
  is a GHz-scale detuning uncertainty at J′≈27, and it moves those components
  by up to 14 MHz with the source of the B rotational energies.


## Uncertainties of balanced-expanded (October 21, 2026)

`ResearchModel("balanced_expanded").uncertainty(line, component=None, ref_component=None, ref_line=None)`
returns the 1σ (MHz) of any prediction the model makes:

- a line center;
- a hyperfine component;
- an intraline interval (`component` minus `ref_component` of the same line);
- an interval between two lines (with `ref_line`).

It returns the parts and the data status with it:

```python
m = ResearchModel("balanced_expanded")
u = m.uncertainty("R(41) 51-0", "a21")
u.sigma_MHz, u.statistical_MHz, u.model_MHz, u.hyperfine_MHz   # 1.50 = sqrt(0.63² + 1.28² + 0.46²)
u.status, u.line_class, u.hyperfine_status                    # measured, beyond, validated: held-out intervals
m.uncertainty("R(41) 51-0", "a21", "a1").sigma_MHz           # 1.10: an interval has no center error
```

`uncertainty` also works from the command line:
`i2spec research predict 'R(41) 51-0' --model balanced_expanded --component a21`.

The model has three parts. Each was built from the final fit and its
cross-validation only, then frozen as `uncertainty.json` in the release, with its
recipe, checks and input hashes. No pending-set outcome entered it.

**Statistical.** This is the linearized covariance of the final fit. Its
parameters are:

- the 38 free new terms;
- the 47 frozen core potential coefficients, refitted jointly in the
  linearization, so their covariance and cross terms come from the same rows;
- the 7 inherited dispersion coefficients, with a 10% prior;
- the 5 calibrations, with their own fits as priors.

It is the sandwich covariance of the weighted estimator actually used, with each
data category scaled by its Birge ratio (at least 1). It is computed from an SVD
of the weighted Jacobian: the normal equations have condition number 9×10¹³ and
leak spurious variance into well-determined lines. It reaches each line through
the Hellmann–Feynman derivatives of its two levels, so lines that share levels
or parameters are correlated.

**Model.** This is the error of predicting a level from the data around it. It
was fitted by marginal likelihood to held-out residuals, with each held-out
level a random effect. Each row's measurement σ was removed in quadrature; atlases
and level constants were scaled by their in-sample Birge ratios. Martin's stated
0.005 cm⁻¹ is 2.6 times too large.

| Region | σ_model | Held-out evidence |
|---|---|---|
| Core line, both levels within 15 in J of comb or Doppler-free data | 0.12 MHz | 47 comb lines (research-balanced's folds) |
| Core line, otherwise | 0.67 MHz | 40 lines, including a whole level held out (B v′39: 3.4 MHz) |
| B v′≥44 level | 0.34 MHz · e^{0.136(v′−44)} · e^{1.67(J/100)²} | 170 held-out levels (vmod folds), J clipped at 186 |
| X v″≥18 level | 16 MHz · e^{0.10(J/100)²} | 574 held-out levels; slopes constrained ≥ 0 |

**Hyperfine.** A component of a line carries an extra σ_HF = √(a² + (κ·o)²), where
o is its offset from the center. An interval Δ within a line carries
√(2a² + (κΔ)²) alone. Both were fitted to held-out intervals.

| Region | a | κ (fraction of an interval) |
|---|---|---|
| Core | 0.037 MHz | 2.9×10⁻⁵ |
| B v′44–75 | 0.064 MHz | κ = e^{−10.38}(1 + (v′−42)/2)^{2.11}: 1.3×10⁻⁴ at v′44, 0.009 at v′70, 0.013 at v′75; held constant beyond |
| X v″18–54 | core a | log-linear from the core to 0.025 at v″51 (within-line scatter of the emission components); held constant beyond |

The B scale term comes from two held-out sources: the intervals of the nested
cross-validation (v′44–62), and the eqQ constants of Chen 2004 and Pique 1986
(v′42–75), converted through |eqQ′ − eqQ″|. An exponential in v′ would fit the
intervals, but it overshoots the held-out constants at v′65–75 by a factor of 5–9.

**Data status.** A level is:

- *measured* when training rows reach its v within 15 in J;
- *interpolated* when it lies inside its v's J range, or between neighbors
  within 2 in v;
- *extrapolated* otherwise.

A line takes the worse of its two levels.

Hyperfine status reads as follows:

| Region | Hyperfine status |
|---|---|
| B v′≤43, X v″≤17 | core |
| B v′44–62 | validated: held-out intervals |
| B v′63–75 | validated: held-out constants |
| B v′77–78 | resonance (Pique 1986), not validated |
| X v″48–54 | validated: emission components |
| X v″18–47 | interpolated |
| B v′76, 79–80 and X v″>54 | extrapolated |

**Typical values** (1σ of the center; of a15 − a1):

| Line | Center | Interval | Status |
|---|---:|---:|---|
| R(56) 32-0 | 0.42 MHz | 0.058 MHz | measured, precise |
| R(150) 20-2 | 0.88 MHz | 0.058 MHz | measured, atlas |
| R(41) 51-0 | 1.4 MHz | 1.1 MHz | measured |
| R(98) 58-1 | 13 MHz | 2.9 MHz | measured |
| R(20) 70-0 | 14 MHz | 7.7 MHz | measured |
| R(10) 77-0 | 35 MHz | 12 MHz | measured; hyperfine at the 1g resonance |
| R(50) 10-30 | 28 MHz | 0.34 MHz | measured |
| R(10) 2-90 | 91 MHz | 18 MHz | extrapolated |
| R(3) 1-108 | 182 MHz | 18 MHz | extrapolated |

**Checks.**

- Calibration on the held-out rows the model was fitted to: 94–97% within 2σ in
  every band beyond the core. Per held-out level, 79% are within 1σ and 95%
  within 2σ. Core comb lines are over-covered: 95% within 2σ, rms z 0.71. The
  statistical part dominates there, through the fit's 0.5 MHz comb floor.
- The runtime reproduces the prototype's level factors within 0.6% (grid tiers
  differ) and 69 frozen reference values exactly.
- Two restarts that reach the same cost differ by up to 38 MHz at X v″78–90.
  The statistical σ there is 35–90 MHz.
- The fold-refit jackknife (g = 4) exceeds the new terms' statistical σ by a
  factor of 2–4 at B v′56–76. The linearization understates the parameters'
  sensitivity there; the model term, fitted to the same folds' held-out
  residuals, carries it.

**Post-hoc coverage on the pending sets.** The accuracy results of the blind
comparison had been read before this model existed, so the following is not a
blind test. The recipe, however, came from the training fit and its
cross-validation and was not tuned on these rows. With z = residual /
√(σ_model² + σ_meas²):

| Pending rows | n | within 1σ | within 2σ |
|---|---:|---:|---:|
| blind, all | 347 | 66% | 91% |
| blind, absolute components | 103 | 66% | 88% |
| blind, intraline, core | 174 | 66% | 91% |
| blind, intraline, beyond the core (Weller, Bordé) | 48 | 52% | 94% |
| blind, between lines (Bertinetto) | 22 | 95% | 100% |
| seen (in training through compilations) | 126 | 70% | 95% |

Nominal values are 68% and 95%. The one set clearly outside is Nishiyama 2013:
its 15 rows miss by a median 0.64 MHz, against a median σ of 0.16 MHz, so only 27%
are within 2σ. The only blind coverage test left is a prospective one, with new
measurements scored against frozen predictions (see the
[validation workflow](research-model-validation.md)).

**Limits.**

- The core model term has no J dependence: its held-out evidence is comb lines
  at J ≤ 150. About 32,000 exported core lines with J″ > 150 lie beyond their
  levels' data. They are flagged extrapolated, and their σ (1–3 MHz) is
  probably optimistic.
- Beyond the core, the model term is the error of predicting a level from its
  neighbors. For levels with their own data it is an upper estimate.
- X hyperfine at v″18–47 has no data.
- The resonance lines' hyperfine σ is extrapolated.

## Preferred in the apps

`i2spec.preferred` makes balanced-expanded the apps' model. It solves every level
of the domain once, together with its statistical factor and data status (about
80 s on 12 cores). The result is cached in `~/.cache/i2spec`, keyed by
the release files' SHA256.

- **Line list.** `lookup.Catalog` replaces each ¹²⁷I₂ line's position with
  balanced-expanded's wherever it holds both levels. The line strengths stay
  i2spec2026o's.
- **Uncertainties and flags.** `lookup.uncertainty` and `Line.flags` report its
  σ and data status. `Line.source` names the model that gave the line.
- **Hyperfine.** `Catalog.components` gives its hyperfine components, with the
  J±2 couplings, each with its own σ.
- **Explorer export.** For every such line the export carries balanced-expanded's
  position, uncertainty and hyperfine pattern. It also carries i2spec2026o's
  position minus it (`d26`, shown as "i2spec2026o" in the line's details).
- **Measured data.** The explorer shows each measurement with its own
  uncertainty and its residual. Since October 21 that includes the 36 data sets
  of the scored blind test, marked "blind test, not fitted"; they still carry
  their `exclude` line, so no fit uses them. Sets transcribed later stay out of
  the export until a registry scores them.

`I2SPEC_MODEL=i2spec2026o` switches the apps back to i2spec2026o throughout.
