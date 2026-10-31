# Uncertainty estimates

Every line position, hyperfine component and hyperfine interval comes with a standard (1σ) uncertainty from the
model's frozen uncertainty model, which has three parts:

$$
\begin{aligned}
\sigma_\text{center}^2 &= \sigma_\text{stat}^2 + \sigma_\text{model}^2,\\
\sigma_\text{component}^2 &= \sigma_\text{center}^2 + \sigma_\text{HF}^2,\qquad \sigma_\text{HF}^2 = a^2 + (\kappa\,o)^2,\\
\sigma_\text{interval}^2 &= 2a^2 + (\kappa\,\Delta)^2.
\end{aligned} \tag{1}
$$

o is a component's offset from the line center; an interval Δ within one line carries the hyperfine part alone,
since its two components share the center. The statistical part comes from the fit; the model and hyperfine parts
from predictions of measurements the fit did not see (its cross-validation and the calibration sets), and a further
set, kept out of both, tested the result once ([the locked test](#the-locked-test)). The data behind each part are
listed on [The fit](fit.md#the-uncertainty-model); the derivations are in the paper.

## Statistical part

σ<sub>stat</sub> is the uncertainty the data leave in {{ unc.n_parameters }} parameters (the {{ fit.n_tails }} terms
beyond the core, the {{ fit.n_potential }} core potential coefficients, the seven dispersion coefficients with a
10 % prior and the five atlas calibrations), from the covariance of the weighted fit with each data category scaled
by its Birge ratio where that exceeds 1. It reaches a line through the derivatives of its two level energies
(Hellmann–Feynman theorem), so lines that share levels or parameters are correlated, and the uncertainty of a
difference of two lines includes that correlation (`ResearchModel.uncertainty(line, component, ref_component,
ref_line)`).

## Model part

σ<sub>model</sub> is the error of predicting a level from the data around it, which σ<sub>stat</sub> cannot see.
It was fitted by maximum likelihood to held-out residuals, each held-out level sharing one error among its rows and
each laboratory's absolute frequencies an offset ({{ unc.lab_offset }} MHz in the core: measurement error, not part
of a prediction's σ).

| levels | σ<sub>model</sub> | evidence |
|---|---|---|
| core line, both levels within 15 in J of comb or Doppler-free data | {{ unc.core_precise }} MHz | {{ unc.core_precise_groups }} lines held out; {{ unc.core_precise_cal }} calibration rows |
| core line, otherwise | {{ unc.core_atlas }} MHz | {{ unc.core_atlas_groups }} lines held out; {{ unc.core_atlas_cal }} calibration rows |
| B v′ ≥ 44, per level | {{ unc.b0 }} MHz × e<sup>{{ unc.b_slope }} (v′ − 44) + {{ unc.b_quad }} (v′ − 44)²</sup> × e<sup>{{ unc.b_j2 }} (J/100)²</sup> | {{ unc.b_groups }} held-out levels (J used up to {{ unc.b_clip }}); {{ unc.b_cal }} calibration rows |
| X v″ ≥ 18, per level | {{ unc.x0 }} MHz × e<sup>{{ unc.x_j2 }} (J/100)²</sup> | {{ unc.x_groups }} held-out levels |

*Table 1. The model part. A core line's value is shared equally by its two levels.*

The B part does not fall below a core level's share. **Above the fitted J**, beyond the J<sub>r</sub> the training
rows reach at a level's v (or within 2 in v) plus 15, the model part grows as J(J + 1)/J<sub>r</sub>(J<sub>r</sub> + 1),
the scale of the centrifugal term. This is a rule, not a fit, and such levels are marked *extrapolated*.

## Hyperfine part

The floor a and the fractional error κ of the scale of the pattern were fitted to hyperfine intervals the fit did
not use, the pattern-scale error shared by the components of a line:

| levels | a | κ |
|---|---|---|
| core | {{ unc.hf_a }} kHz | {{ unc.hf_kappa }} |
| B v′ = 44–75 | {{ unc.hf_a }} kHz × e<sup>{{ unc.hf_growth }} (v′ − 44)</sup>: {{ unc.hf_a_beyond }} | {{ unc.hf_kappa_beyond }} |
| X v″ = 18–54 | {{ unc.hf_a }} kHz | from the core value to {{ unc.x_kappa51 }} at v″ = 51 |

*Table 2. The hyperfine part.*

Beyond the last validated level κ is held constant and the component is marked extrapolated; at B v′ = 77–78,
perturbed by a near-resonant 1g level, it is marked as not validated.

## Data status

| status | meaning |
|---|---|
| *measured* | training rows reach both levels' v within 15 in J |
| *interpolated* | a level lies inside its v's J range, or between neighbors within 2 in v |
| *extrapolated* | a level lies beyond the data in J or v: treat the value with care |

Lines whose levels the model does not hold (B v′ > 80 for ¹²⁷I₂, or a level beyond its barrier rule) carry the flag
*outside the fitted domain* and a rough bound of 1 GHz.

<figure class="fig" markdown="span">
<div class="svg" data-svg="figures/uncertainty_map.svg" role="img" aria-label="Map of position uncertainty over upper and lower vibrational quantum numbers"></div>
<figcaption markdown="span">**Figure 1.** Uncertainty of the R(20) line of each band v′–v″ (blank: a level is dissociated at this J). Inside the core (v″ ≤ 17, v′ ≤ 43) almost every line is known to better than 1 MHz; beyond it the uncertainty grows as the data thin out, to 0.1–0.5 GHz for the last levels.</figcaption>
</figure>

| line | center | a15 − a1 | status |
|---|---:|---:|---|
{{ unc.typical_table }}

*Table 3. Uncertainties of a few ¹²⁷I₂ lines.*

## Isotopologues

¹²⁹I₂ and ¹²⁷I¹²⁹I lines have the same three parts, read at the equivalent ¹²⁷I₂ level (the same J = 0 energy and
ρ J(J + 1)). The statistical part is conditioned on the {{ fit.n_iso_centers }} isotopologue lines of the fit; each B
level adds the error of the isotope shift and the spread of the equally good forms of Eq. (4) of the
[levels page](levels.md#isotopologues); each component adds a floor of 0.10 MHz. A B level more than 2 in v from the
lines that fix the shift (v′ = 6–14) is marked *extrapolated*.

| line | center | a1 | a2 − a1 | status |
|---|---:|---:|---:|---|
{{ unc.typical_iso_table }}

*Table 4. Uncertainties of a few isotopologue lines.*

## How well the uncertainties work

On the calibration sets, to which the uncertainty model was fitted together with the cross-validation:

| rows | n | within 1σ | within 2σ |
|---|---:|---:|---:|
| all | {{ unc.cal_rows }} | {{ unc.cal_1 }} | {{ unc.cal_2 }} |
| ¹²⁷I₂ in the core | {{ unc.cal_core_rows }} | {{ unc.cal_core_1 }} | {{ unc.cal_core_2 }} |
| hyperfine intervals within a line, core | {{ unc.cal_hf_core_rows }} | {{ unc.cal_hf_core_1 }} | {{ unc.cal_hf_core_2 }} |
| hyperfine intervals within a line, B v′ = 50–79 | {{ unc.cal_hf_b_rows }} | {{ unc.cal_hf_b_1 }} | {{ unc.cal_hf_b_2 }} |

*Table 5. Coverage on the calibration sets (68 % and 95 % for a normal distribution).*

### The locked test

{{ fit.locked_sets }} data sets ({{ unc.locked_rows }} measurements) were kept out of every fit and of the
uncertainty model, and scored once. The criteria were recorded before any value was read: 95 % of the rows within
2σ as the target, 90 % as the floor, for all rows and for each kind of row (a kind with fewer than 10 rows reported
but not judged). Nothing was changed after scoring.

| rows | n | within 1σ | within 2σ | rms z | result |
|---|---:|---:|---:|---:|---|
{{ unc.locked_table }}

*Table 6. The locked test. z is the residual over √(σ² + the measurement's own uncertainty²).*

!!! warning "Known limitation"
    Hyperfine intervals within a line are covered within 2σ for {{ unc.locked_hf_2 }} of the locked-test
    measurements, below the floor of 90 %. Absolute frequencies ({{ unc.locked_abs_2 }}) and the other
    isotopologues meet the target.

One data set causes most of the shortfall: the 20 intervals of R(115) 20-1 of Cruz *et al.* (1994), measured on a
scanned Fabry–Perot scale, are 0.88 % smaller than the model's, which a scale error of the interferometer or an
error of the model's quadrupole constant at B v′ = 20 would each explain. Without that set, 94 % of the hyperfine
intervals fall within 2σ (given beside the score, not in place of it).

## Worked example

!!! example "R(56) 32–0: the uncertainty of the position and of a10"
    The line is a core line near comb data, with data status *{{ we.status }}*. Eq. (1) gives:

    | part | value |
    |---|--:|
    | σ<sub>stat</sub>, from the fit's covariance through its two levels | {{ we.sigma_stat }} MHz |
    | σ<sub>model</sub>, {{ we.sigma_level }} MHz per level (Table 1) | {{ we.sigma_model }} MHz |
    | **σ<sub>center</sub>** = √(σ<sub>stat</sub>² + σ<sub>model</sub>²) | **{{ we.sigma }} MHz** |
    | σ<sub>HF</sub> of a10: a = {{ we.hf_a }} kHz, κ = {{ we.hf_kappa }}, o = {{ we.a10 }} MHz | {{ we.sigma_hf }} MHz |
    | **σ of a10** | **{{ we.sigma_comp }} MHz** |
    | σ of the interval a15 − a1 = {{ we.span }} MHz | {{ we.sigma_interval }} MHz |

    The model's a10, {{ we.comp_MHz }} MHz, differs from the value of {{ we.measured_short }} by
    {{ we.residual }} MHz, z = {{ we.z }}. That measurement was not in the fit ([The fit](fit.md#worked-example)).

    ```python
    from i2spec import ResearchModel

    m = ResearchModel("i2spec")
    u = m.uncertainty("R(56) 32-0", "a10")
    u.sigma_MHz, u.statistical_MHz, u.model_MHz, u.hyperfine_MHz, u.status
    m.uncertainty("R(56) 32-0", "a15", "a1").sigma_MHz   # {{ we.sigma_interval }} MHz
    ```
