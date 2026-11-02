# The fit

This page lists which parameters of the model are fitted, to which data, and which are taken from the literature.
The model has {{ fit.fitted }} fitted molecular coefficients, {{ fit.inherited }} inherited ones and
{{ fit.calibration }} calibration terms of the atlases. Every coefficient is shared by many levels. The
[references page](references.md) lists the data sets with their sources. The paper describes the methods: the
weights, the solvers and the choice of each functional form.

| stage | fitted | data | levels |
|---|---|---|---|
| [core fit](#the-core-fit) | {{ fit.n_potential }} potential, {{ fit.n_hf_core }} B hyperfine and {{ fit.n_shift }} isotope-shift coefficients; {{ fit.core_calibration }} atlas calibrations | {{ fit.core_rows }} rows: precision lines, two atlases, level constants; {{ fit.hf_core_trained }} hyperfine intervals; {{ fit.n_iso_centers }} isotopologue lines | X v″ ≤ 17, B v′ ≤ 43, J ≤ 200 |
| [continuation](#beyond-the-core) | {{ fit.n_tails }} terms of the potentials; 1 atlas offset | {{ fit.continuation_rows }} rows: the core's, and every row beyond it | X v″ ≤ 108, B v′ ≤ 80 |
| [hyperfine beyond the core](#hyperfine-beyond-the-core) | {{ fit.n_hf_beyond }}: seven for B (1u and 1g mixing, curvatures of eqQ, C, d and δ, 1g scale), three for X | {{ fit.hf_beyond_rows }} intervals, {{ fit.hf_chen }} + {{ fit.hf_pique }} B constants, X eqQ near the limit | B v′ ≥ 44, X v″ ≥ 18 |
| [uncertainty model](#the-uncertainty-model) | σ<sub>model</sub> and σ<sub>HF</sub> | cross-validation; {{ fit.calibrate_rows }} rows of {{ fit.calibrate_sets }} calibration sets | all |

*Table 1. The stages of the fit, in order. Each stage keeps the coefficients of the stages before it.*

## The core fit

The core fit determines, in one weighted least-squares problem, the potentials of X and B
(Eqs. (2) and (3) of the [levels page](levels.md)), the radial functions of the B hyperfine parameters
(Eq. (3) of the [hyperfine page](hyperfine.md)) and the isotope shift (Eq. (4) of the levels page):

- X: R<sub>e</sub>, D<sub>e</sub>, {{ fit.x_c }} exponent and {{ fit.x_alpha }} rotational coefficients;
  B: R<sub>e</sub>, T<sub>e</sub>, {{ fit.b_c }} exponent and {{ fit.b_alpha }} rotational coefficients
  ({{ fit.n_potential }} in all; D<sub>e</sub> of B follows from the X limit);
- the B hyperfine functions: {{ fit.hf_eqq }}, {{ fit.hf_c }}, {{ fit.hf_d }} and {{ fit.hf_delta }} coefficients for
  eqQ, C, d and δ;
- the isotope shift a₀ and a₁.

Measured hyperfine components are reduced to hyperfine-free line centers with the model's own hyperfine structure.
The fit alternates between the potentials and the hyperfine functions until both are consistent. Its data are of
¹²⁷I₂ with X v″ ≤ 17 and B v′ ≤ 43, apart from the isotopologue lines:

| data | rows | sets |
|---|--:|--:|
| absolute frequencies of components or line centers (comb and other), and intervals between lines | {{ fit.core_comb }} ({{ fit.core_intervals }} intervals) | {{ fit.core_comb_sets }} |
| Doppler-free absolute frequencies | {{ fit.core_doppler }} | {{ fit.core_doppler_sets }} |
| hyperfine intervals within a line | {{ fit.hf_core_trained }} | {{ fit.hf_core_sources }} |
| Salami & Ross atlas (2005); scale offset {{ fit.sr_offset }} MHz, slope {{ fit.sr_slope }} MHz per 1000 cm⁻¹ | {{ fit.core_sr }} | 1 |
| NIST atlas of the APO cell; scale offset {{ fit.apo_offset }} MHz, slope {{ fit.apo_slope }} MHz per 1000 cm⁻¹ | {{ fit.core_apo }} | 1 |
| Orsay atlas part I, X v″ = {{ fit.orsay1_core_v }} | {{ fit.core_orsay1 }} | 1 |
| X level constants of Martin *et al.* (1986), v″ = {{ fit.martin_core_v }} | {{ fit.core_martin }} | 1 |
| isotopologue line centers of the CIPM compilation (2003) | {{ fit.n_iso_centers }} | – |

*Table 2. The data of the core fit. Each atlas's scale is calibrated against the precision lines it shares with an
offset and a slope ({{ fit.sr_shared }} and {{ fit.apo_shared }} lines).*

Each measurement enters once: where a compilation (the CIPM and BIPM tables) reprints a primary paper, the primary
paper is used. The six isotopologue lines ({{ fit.iso_centers }}) separate the shape of each potential from its
rotational function, which the ¹²⁷I₂ levels alone leave open.

**The benchmark.** The {{ fit.benchmark_sets }} precision data sets transcribed last were held out of the core fit
({{ fit.benchmark_rows }} line rows). The fit predicts their {{ fit.benchmark_absolute }} absolute component
frequencies to {{ fit.benchmark_absolute_rms }} MHz rms and their {{ fit.benchmark_intraline }} hyperfine intervals
to {{ fit.benchmark_intraline_rms }} MHz. The 515 nm line {{ fit.exclusion }} is kept out of every fit and score as a
diagnostic.

## Beyond the core

The **continuation** keeps every coefficient of the core fit and adds {{ fit.n_tails }} terms that act only beyond
it (see [Beyond the core](levels.md#beyond-the-core)): for X {{ fit.x_outer }} outer, {{ fit.x_inner }} inner and
{{ fit.x_rot }} rotational term, a change of C₆ and {{ fit.x_splines }} B-splines; for B {{ fit.b_outer }} outer and
{{ fit.b_rot }} rotational terms; and one change of the dissociation energy shared by both. Together with one
offset of the Orsay Partie IV scale ({{ fit.orsay4_offset }} MHz), they were fitted to {{ fit.continuation_rows }}
rows: the data beyond the core, plus every training row of the core fit again to hold the core in place.

| levels | training rows and the rms of their residuals |
|---|---|
{{ fit.region_table }}

*Table 3. The data beyond the core: {{ fit.cont_orsay1 }} lines of the Orsay atlas part I (Gerstenkorn, Vergès &
Chevillard 1982), {{ fit.cont_orsay4 }} unblended lines of Partie IV (Gerstenkorn & Luc 1983), {{ fit.cont_martin }}
level energies from the constants of Martin *et al.* (1986) and {{ fit.cont_martin_dissoc }} from Martin *et al.*
(1983), {{ fit.cont_emission }} comb-referenced emission lines, and the comb-referenced and atlas lines of
B v′ = 44–62.*

Beyond the core, the model is tested by removing every row of a level and predicting the level from its neighbors
({{ fit.continuation_folds }} folds by v). Held out this way:

- the comb lines of B v′ = 44–50 are predicted to 1.5 MHz;
- the atlas lines to 22–200 MHz, growing towards the B limit;
- the X level constants to 46–122 MHz;
- the last X levels (v″ ≥ 90) to about 0.7 GHz.

These held-out errors set the model part of the uncertainty there.

## Hyperfine beyond the core

- **B v′ = 44–80.** Each parameter is a background continued from the core fit's functions plus second-order sums
  over the nine molecular states of the I(²P3/2) + I(²P1/2) limit. The curves and couplings of those states are
  taken from the literature. Fitted here are the mixing angles of the two 1u and the two 1g states, the curvature of
  the eqQ continuation and those of the C, d and δ backgrounds. They are fitted to the {{ fit.hf_beyond_rows }}
  measured hyperfine intervals beyond the core and to the effective B constants of Chen *et al.* (2004,
  {{ fit.hf_chen }} values) and Pique *et al.* (1986, {{ fit.hf_pique }}). The long-range scale of the 1g(¹Πg)
  curve was chosen by training cost. The backgrounds of d and δ have orders 5 and 4, fitted to the core fit's own
  values. This form was selected after the blind tests had been scored, so no blind test of it remains; in
  cross-validation it improves eight regions and worsens none.
- **B v′ = 77–78.** The near-resonant 1g level is taken from Pique *et al.* and diagonalized with B. At these two
  levels the second-order sums use a separate fit: the 1u angle fitted alone, the 1g angle at its separated-atom
  value of 30° and a 1g(¹Πg) scale of 0.65. No measurement at v′ = 77–78 tests the coupling that the fitted 1g angle
  would give.
- **X v″ ≥ 18.** eqQ continues from the formula of the core to its separated-atom value along a curve with two
  coefficients fitted to the eqQ of Bacis *et al.* (1980). A slope of C rests on one value and is not validated.
- **X v″ ≤ 17.** The formulae of Bodermann, Knöckel & Tiemann (2002), inherited unchanged (eleven constants).

That makes {{ fit.n_hf_beyond }} fitted hyperfine parameters beyond the core. The other isotopologues take every
hyperfine parameter from the same functions, scaled by the moments of each nucleus.

## Isotopologues

The parameters specific to ¹²⁹I₂ and ¹²⁷I¹²⁹I are the two coefficients a₀ and a₁ of the isotope shift. The core fit
determines them from the {{ fit.n_iso_centers }} CIPM line centers above, each tied through measured intervals to a
measured ¹²⁷I₂ frequency.
The power p = 5 of Eq. (4) of the levels page was chosen by leave-one-line-out error among p = 0, 3, 5 and 8.

## The uncertainty model

The three parts of the [uncertainty](uncertainty.md) come from different data:

- **σ<sub>stat</sub>**: the covariance of the continuation fit, over the same rows.
- **σ<sub>model</sub>**: the residuals of predictions made without the measurement: {{ fit.core_folds }}
  cross-validation folds of the core fit (the benchmark, {{ fit.core_level_folds }} single levels and
  {{ fit.core_source_folds }} single sources), the {{ fit.continuation_folds }} folds by v beyond the core, and the
  calibration sets.
- **σ<sub>HF</sub>**: hyperfine intervals held out in {{ fit.hyperfine_folds }} nested folds by v′, the B constants of
  Chen *et al.* and Pique *et al.*, and the intervals of the calibration sets.

Data sets transcribed after the fit were each given one role, assigned by set and laboratory rather than by their
residuals:

- {{ fit.calibrate_sets }} sets calibrate the uncertainty model ({{ fit.calibrate_rows }} rows);
- {{ fit.locked_sets }} sets were locked away and scored once ([the locked test](uncertainty.md#the-locked-test));
- {{ fit.record_sets }} sets (level and constants tables, widths, lifetimes, duplicates) are kept as a record.

A set was eligible to enter the fit if it improved the held-out error of some region by more than two standard
errors without worsening any. No set met that test, and those sets went to calibration. The explorer marks every set
outside the fit as “not fitted”.

## Worked example

!!! example "R(56) 32–0 in the fit"
    The 532 nm line enters the fit only through its {{ we.train_atlas }} atlas line centers. Its
    {{ we.held_out_precision }} comb frequencies and {{ we.hf_held_out }} hyperfine intervals belong to the benchmark
    held out of the core fit, and so do all {{ we.b_precision_held_out }} precision rows of its upper level,
    B v′ = 32. In the fit, that level is reached by {{ we.b_all }} training rows, of which {{ we.b_precision }}
    are precision rows, and through the coefficients it shares with every other level of the core; the lower
    level, X v″ = 0, by {{ we.x_all }} rows, {{ we.x_precision }} of them precision rows.

    The a10 frequency of the [hyperfine page](hyperfine.md#worked-example), {{ we.comp_MHz }} MHz, is therefore a
    prediction: it differs from the measurement of {{ we.measured_short }} by {{ we.residual }} MHz, or
    {{ we.z }} σ. The benchmark residuals were used afterwards, with the other cross-validation folds, to set the
    model part of the uncertainty.
