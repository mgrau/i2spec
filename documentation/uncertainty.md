# Uncertainty estimates

Every line position, every hyperfine component and every hyperfine interval comes with a standard (1σ)
uncertainty, computed by the model's frozen uncertainty model, which has three parts:

$$
\begin{aligned}
\sigma_\text{center}^2 &= \sigma_\text{stat}^2 + \sigma_\text{model}^2,\\
\sigma_\text{component}^2 &= \sigma_\text{center}^2 + \sigma_\text{HF}^2,\\
\sigma_\text{interval} &= \sqrt{2a^2 + (\kappa\,\Delta)^2}.
\end{aligned} \tag{1}
$$

The statistical part comes from the fit, the model part from the errors of predictions the fit did not see,
and the hyperfine part from the same kind of test for hyperfine intervals. Two components of one line share
the line's center, so a hyperfine interval Δ within a line carries the hyperfine part alone. The model and
hyperfine parts were fitted to residuals of predictions that did not use the measurement predicted: the
cross-validation of the training fit, and data sets kept out of the fit for this purpose (the calibration sets).
A further set of measurements, kept out of both, tested the result once.

Measured values have their own uncertainties. The explorer lists every measurement of a line with its
stated uncertainty beside the model's position and σ, so the two can be compared directly.

## Statistical part

σ<sub>stat</sub> is the uncertainty the data leave in the fitted coefficients, propagated to the line. Its
parameters are:

- the 37 coefficients of the continuation beyond the core;
- the 46 potential coefficients of the core fit, refitted jointly in the linearization so that their
  covariance, and their correlation with the new terms, come from the same data;
- the seven inherited dispersion coefficients, each with a 10 % prior;
- the five calibration terms (atlas scales and the Orsay offset), with their own fits as priors.

The covariance is that of the weighted least-squares estimator actually used. Each data category is scaled
by its own Birge ratio, taken as at least 1:

| category | Birge ratio |
|---|---:|
| NIST APO atlas | 1.52 |
| Salami & Ross atlas | 1.25 |
| Orsay part I | 1.15 |
| Orsay Partie IV | 1.20 |
| comb | 1.17 |
| Doppler-free, emission | 0.62, 0.38 (used as 1) |
| Martin level constants | 0.38 (used as 1); near the X limit 2.22 |

It is computed from a singular-value decomposition of the weighted Jacobian rather than the normal
equations, which are badly conditioned. Eleven directions with singular values below 10⁻⁸ of the
largest are combinations of core coefficients that no level depends on; they are dropped.

A line's σ<sub>stat</sub> follows from the derivatives of its two level energies with respect to every
parameter. The derivatives are computed from the wavefunctions (Hellmann–Feynman theorem). Lines that share
levels or parameters are therefore correlated, and the uncertainty of a difference of two lines is computed
with that correlation.

## Model part

σ<sub>model</sub> is the error of predicting a level from the data around it: a property of the functional
form, which σ<sub>stat</sub> cannot see. It was fitted to residuals of predictions that did not use the
measurement being predicted:

- in the core, the seven cross-validation folds of the core fit. These hold out the newer data sets, the
  precision data of one level at a time, and two single sources;
- beyond the core, four folds that each hold out every row of one vibrational level in four, which is then
  predicted from its neighbors;
- the calibration sets: 2,432 measurements from data sets that no fit used, predicted by the final fit. A line
  that the training data hold inside a CIPM mean, and a value printed again in a later set, count once.

Each residual's measurement uncertainty and its fold's σ<sub>stat</sub> were removed in quadrature. Atlas and
level-constant rows were scaled by their in-sample Birge ratios; the stated precision of Martin *et al.*
(1986) is 2.6 times larger than their scatter. Every held-out level shares one error among its rows (a
random effect). Among the calibration sets, each laboratory's absolute frequencies share an offset as well:
this is measurement error, 0.17 MHz between laboratories in the core, and is not part of a prediction's σ.
The variances were fitted by maximum likelihood:

| levels | σ<sub>model</sub> | evidence |
|---|---|---|
| core line, both levels within 15 in J of comb or Doppler-free data | 0.16 MHz | 47 lines held out; 264 calibration rows |
| core line, otherwise | 0.29 MHz | 40 lines held out, among them a whole level (B v′ = 39); 42 calibration rows |
| B v′ ≥ 44, per level | 0.52 MHz × e<sup>0.046 (v′ − 44) + 0.004 (v′ − 44)²</sup> × e<sup>1.65 (J/100)²</sup> | 170 held-out levels (J used up to 186); 68 calibration rows |
| X v″ ≥ 18, per level | 17 MHz × e<sup>0.11 (J/100)²</sup> | 574 held-out levels (slopes constrained to be non-negative) |

*Table 1. The model part of the uncertainty. A core line's value is shared equally by its two levels.*

The B part grows from 0.5 MHz at v′ = 44 (J = 0) to 3 MHz at v′ = 60, 26 MHz at v′ = 70 and 0.5 GHz at
v′ = 80, and does not fall below a core level's share. The likelihood barely constrains the quadratic term;
its value was chosen from an explicit profile.

## Hyperfine part

A component at offset o from the line center carries σ<sub>HF</sub> = √(a² + (κo)²). The constant a is the
floor of a component; κ is a fractional error of the scale of the pattern. Both were fitted to hyperfine
intervals that the fit did not use:

| levels | a | κ |
|---|---|---|
| core | 63 kHz | 2.1 × 10⁻⁴ |
| B v′ = 44–75 | 63 kHz × e<sup>0.21 (v′ − 44)</sup>: 0.22 MHz at v′ = 50, 1.8 MHz at v′ = 60, 14 MHz at v′ = 70 | 5 × 10⁻⁵ at v′ = 44, 1.5 × 10⁻³ at v′ = 55, 9 × 10⁻³ at v′ = 70, 0.013 at v′ = 75 |
| X v″ = 18–54 | 63 kHz | from the core value to 0.025 at v″ = 51 |

*Table 2. The hyperfine part of the uncertainty.*

In the calibration sets the pattern-scale error is shared by the components of a line, and the floor is
counted once per component, so a line counts once however many of its components were measured. The B values
rest on three kinds of held-out data:

- the intervals of the nested cross-validation (v′ = 44–62);
- the B-state quadrupole constants measured by Chen *et al.* (2004) and Pique *et al.* (1986)
  (v′ = 42–75), converted through the difference of the upper and lower quadrupole constants;
- 1,102 hyperfine intervals of the calibration sets.

The X value at v″ = 48–54 is the scatter of the emission components about their line centers. Beyond the
last validated level κ is held constant, and the component is marked extrapolated. At B v′ = 77–78 the
hyperfine structure is perturbed by a near-resonant 1g level (Pique *et al.* 1986); there it is marked as not
validated.

## Data status

Each line is marked with how directly the data constrain its levels:

| status | meaning |
|---|---|
| *measured* | training rows reach both levels' v within 15 in J |
| *interpolated* | a level lies inside its v's J range, or between neighbors within 2 in v |
| *extrapolated* | a level lies beyond the data in J or v: treat the value with care |

About 23 % of the exported ¹²⁷I₂ lines are extrapolated. Most of them (47 000) are core lines with J″ > 150,
beyond the J their levels' data reach. Their σ (median 1.8 MHz) grows with J through the rule
[below](#above-the-fitted-j), which is not fitted: treat it as an estimate.

<figure class="fig" markdown="span">
<div class="svg" data-svg="figures/uncertainty_map.svg" role="img" aria-label="Map of position uncertainty over upper and lower vibrational quantum numbers"></div>
<figcaption markdown="span">**Figure 1.** Uncertainty of the R(20) line of each band v′–v″ (blank: no such line, a level is dissociated at this J). Inside the core fit (v″ ≤ 17, v′ ≤ 43) 95 % of the lines are known to better than 1 MHz, and to 0.2–0.5 MHz where comb data reach both levels. Beyond it the uncertainty is 1–15 MHz for B v′ = 44–67 (comb data and the APO atlas), 20–600 MHz for B v′ ≥ 68, 17–30 MHz for X v″ = 18–60 (the Orsay atlas and the level constants of Martin *et al.*), and 0.1–0.5 GHz for the last X levels.</figcaption>
</figure>

## Typical values

| line | center | a15 − a1 | status |
|---|---:|---:|---|
| R(56) 32-0 (532 nm) | 0.39 MHz | 0.20 MHz | measured |
| R(150) 20-2 | 0.58 MHz | 0.20 MHz | measured |
| R(41) 51-0 | 1.3 MHz | 0.64 MHz | measured |
| R(98) 58-1 | 12 MHz | 2.8 MHz | measured |
| R(20) 70-0 | 28 MHz | 22 MHz | measured |
| R(50) 10-30 | 28 MHz | 1.1 MHz | measured |
| R(10) 2-90 | 91 MHz | 15 MHz | extrapolated |

*Table 3. Uncertainties of a few lines.*

## How well the uncertainties work

On the cross-validation residuals, the uncertainties hold 93–96 % of the rows within 2σ in every band of B beyond
the core and 95–97 % in most bands of X. Taking each held-out level's mean residual instead, 79 % are within 1σ
and 96 % within 2σ, against 68 % and 95 % for a normal distribution. Core comb lines are over-covered (rms
z = 0.5–0.7), because the statistical part carries the fit's 0.3 MHz floor on the comb data. On the calibration
sets, which the uncertainty model was fitted to as well:

| rows | n | within 1σ | within 2σ |
|---|---:|---:|---:|
| all | 2,432 | 80 % | 96 % |
| ¹²⁷I₂ in the core | 1,150 | 80 % | 95 % |
| hyperfine intervals within a line, core | 844 | 82 % | 96 % |
| hyperfine intervals within a line, B v′ = 50–79 | 858 | 46–85 % | 93–97 % |

*Table 4. Coverage on the calibration sets.*

Earlier versions of the uncertainty model, fitted to the cross-validation alone, were tested blind three times
on data sets transcribed later. Their hyperfine part was too small: in the largest test, 82 % of 2,022 hyperfine
intervals within a line fell within 2σ. Most of those data sets then became the calibration sets, and seventeen
the locked test below.

### The locked test

Seventeen data sets (164 measurements) were kept out of every fit and of the uncertainty model, and scored once
against it. The criteria were frozen and recorded before any value was read: 95 % of the rows within 2σ as the
target, 90 % as the floor, for all rows and for each kind of row, a kind with fewer than 10 rows reported but not
judged. Nothing was changed after scoring.

| rows | n | within 1σ | within 2σ | rms z | result |
|---|---:|---:|---:|---:|---|
| all | 164 | 66 % | 88 % | 1.70 | below the floor |
| ¹²⁷I₂ absolute frequencies | 44 | 68 % | 98 % | 1.10 | target met |
| ¹²⁷I₂ intervals between lines | 2 | 100 % | 100 % | 0.67 | not judged |
| ¹²⁷I₂ hyperfine intervals within a line | 86 | 57 % | 78 % | 2.18 | below the floor |
| ¹²⁹I₂ and ¹²⁷I¹²⁹I | 32 | 84 % | 100 % | 0.63 | target met |

*Table 5. The locked test. z is the residual over √(σ² + the measurement's own uncertainty²).*

!!! warning "Known limitation"
    Hyperfine intervals within a line are covered within 2σ for 78 % of the locked-test measurements, below
    the floor of 90 %. Absolute frequencies and the other isotopologues meet the target.

One data set causes most of the shortfall. Cruz *et al.* (1994) measured 20 intervals of R(115) 20-1 on a
scanned Fabry–Perot scale; 5 of them fall within 2σ. Their intervals are 0.88 % smaller than the model's, in
proportion to their size, and one scale factor leaves an rms of 0.61 MHz, as in the paper's own fit. The
paper's difference of the upper and lower quadrupole constants, 1911.8(0.9) MHz, is 0.84 % below the model's
1928.1 MHz. A scale error of the interferometer and an error of the model's quadrupole constant at B v′ = 20,
J′ = 116 would each give this, and the set cannot tell them apart. The other 144 rows are within 2σ for 96.5 %,
and the hyperfine intervals without this set for 94 % (62 of 66); these numbers are given beside the score, not
in place of it.

## Above the fitted J

The fit's data reach J = 200, and at most levels much less. Beyond the J the training rows reach at a level's v
(or within 2 in v), plus 15, the level's model part grows as J(J + 1), the size of the centrifugal term through
which the rotational function and the shape of the potential act: σ<sub>model</sub> × J(J + 1)/J<sub>r</sub>(J<sub>r</sub> + 1). This
is a rule, not fitted. Such levels are marked *extrapolated*. The 332 lines of the Orsay atlas part I with
J″ = 200–271 scatter about the model by 50 MHz up to J″ = 225, the atlas's own scatter at lower J, and by
65–100 MHz above.

## Isotopologues

The lines of ¹²⁹I₂ and ¹²⁷I¹²⁹I have the same three parts, read from the ¹²⁷I₂ tables at the equivalent
level (the same J = 0 energy and ρ J(J + 1)), with three additions:

- **Statistical part, conditioned.** The ¹²⁷I₂ levels do not separate the shape of a potential from its
  rotational function, which an isotopologue sees scaled by ρ: combinations of coefficients that leave every
  ¹²⁷I₂ level unchanged move the isotopologue levels as J(J + 1). The covariance of the fit and of the two
  isotope-shift coefficients is therefore conditioned on the six CIPM isotopologue lines. It leaves
  0.09–0.32 MHz at those lines, and grows with the distance from them in v and J.
- **Isotope shift.** Each B level adds the leave-one-line-out rms of the shift fit, 0.081 MHz for ¹²⁹I₂ and
  half that for ¹²⁷I¹²⁹I, and the spread of the equally good forms of Eq. (5) of the
  [levels page](levels.md#isotopologues) (Table 4 there).
- **Hyperfine.** Each component adds a floor of 0.10 MHz, the excess of the CIPM intervals over their stated
  uncertainties.

A B level whose ¹²⁷I₂-equivalent v lies more than 2 from the lines that fix the shift (v′ = 6–14) is marked
*extrapolated*.

| line | center | a1 | a2 − a1 | status |
|---|---:|---:|---:|---|
| ¹²⁹I₂ P(54) 8-4 (633 nm) | 0.21 MHz | 0.25 MHz | 0.17 MHz | interpolated |
| ¹²⁹I₂ R(20) 15-1 | 0.72 MHz | 0.73 MHz | 0.17 MHz | interpolated |
| ¹²⁹I₂ R(56) 32-0 (532 nm) | 3.4 MHz | 3.4 MHz | 0.17 MHz | extrapolated |
| ¹²⁹I₂ R(41) 51-0 | 12 MHz | 12 MHz | 0.41 MHz | extrapolated |
| ¹²⁷I¹²⁹I P(33) 6-3 (633 nm) | 0.19 MHz | 0.24 MHz | 0.17 MHz | interpolated |
| ¹²⁷I¹²⁹I R(56) 32-0 | 1.7 MHz | 1.7 MHz | 0.17 MHz | extrapolated |

*Table 6. Uncertainties of a few isotopologue lines.*

Lines whose levels the model does not hold (B v′ > 80 for ¹²⁷I₂, or a level beyond its barrier rule) carry
the flag *outside the fitted domain* and an estimated uncertainty of at least 1 GHz.
