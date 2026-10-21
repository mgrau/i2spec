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
the line's center, so a hyperfine interval Δ within a line carries the hyperfine part alone. The model was
built from the training fit and its cross-validation only, and frozen before it was compared with any new
measurement.

Measured values have their own uncertainties. The explorer lists every measurement of a line with its
stated uncertainty beside the model's position and σ, so the two can be compared directly.

## Statistical part

σ<sub>stat</sub> is the uncertainty the data leave in the fitted coefficients, propagated to the line. Its
parameters are:

- the 38 coefficients of the continuation beyond the core;
- the 47 potential coefficients of the core fit, refitted jointly in the linearization so that their
  covariance, and their correlation with the new terms, come from the same data;
- the seven inherited dispersion coefficients, each with a 10 % prior;
- the five calibration terms (atlas scales and the Orsay offset), with their own fits as priors.

The covariance is that of the weighted least-squares estimator actually used. Each data category is scaled
by its own Birge ratio, taken as at least 1:

| category | Birge ratio |
|---|---:|
| NIST APO atlas | 1.52 |
| Salami & Ross atlas | 1.25 |
| Orsay part I | 1.14 |
| Orsay Partie IV | 1.20 |
| comb, Doppler-free, emission | 0.77, 0.61, 0.36 (used as 1) |
| Martin level constants | 0.38 (used as 1) |

It is computed from a singular-value decomposition of the weighted Jacobian rather than the normal
equations, whose condition number is 9 × 10¹³. Eleven directions with singular values below 10⁻⁸ of the
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
  predicted from its neighbors.

Each residual's measurement uncertainty and its fold's σ<sub>stat</sub> were removed in quadrature. Atlas and
level-constant rows were scaled by their in-sample Birge ratios; the stated precision of Martin *et al.*
(1986) is 2.6 times larger than their scatter. Every held-out level shares one error among its rows (a
random effect), and the variance was fitted by maximum likelihood:

| levels | σ<sub>model</sub> | evidence |
|---|---|---|
| core line, both levels within 15 in J of comb or Doppler-free data | 0.12 MHz | 47 lines |
| core line, otherwise | 0.67 MHz | 40 lines, among them a whole level held out (B v′ = 39: 3.4 MHz) |
| B v′ ≥ 44, per level | 0.34 MHz × e<sup>0.136 (v′ − 44)</sup> × e<sup>1.67 (J/100)²</sup> | 170 held-out levels (J used up to 186) |
| X v″ ≥ 18, per level | 16 MHz × e<sup>0.10 (J/100)²</sup> | 574 held-out levels (slopes constrained to be non-negative) |

*Table 1. The model part of the uncertainty.*

## Hyperfine part

A component at offset o from the line center carries σ<sub>HF</sub> = √(a² + (κo)²). The constant a is the
floor of a component; κ is a fractional error of the scale of the pattern. Both were fitted to hyperfine
intervals that the fit did not use:

| levels | a | κ |
|---|---|---|
| core | 37 kHz | 2.9 × 10⁻⁵ |
| B v′ = 44–75 | 64 kHz | 1.3 × 10⁻⁴ at v′ = 44, 2 × 10⁻³ at v′ = 55, 9 × 10⁻³ at v′ = 70, 0.013 at v′ = 75 |
| X v″ = 18–54 | 37 kHz | from the core value to 0.025 at v″ = 51 |

*Table 2. The hyperfine part of the uncertainty.*

The B scale error rests on two kinds of held-out data:

- the intervals of the nested cross-validation (v′ = 44–62);
- the B-state quadrupole constants measured by Chen *et al.* (2004) and Pique *et al.* (1986)
  (v′ = 42–75), converted through the difference of the upper and lower quadrupole constants.

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

About 15 % of the exported lines are extrapolated. Most of them (32 000) are core lines with J″ > 150,
beyond the J their levels' data reach. The core's model part has no J dependence, because its evidence is comb
lines at J ≤ 150, so the σ of these lines (typically 1–3 MHz) grows only through the statistical part and is
probably optimistic.

<figure class="fig" markdown="span">
<div class="svg" data-svg="figures/uncertainty_map.svg" role="img" aria-label="Map of position uncertainty over upper and lower vibrational quantum numbers"></div>
<figcaption markdown="span">**Figure 1.** Uncertainty of the R(20) line of each band v′–v″ (blank: no such line, a level is dissociated at this J). Inside the core fit (v″ ≤ 17, v′ ≤ 43) the lines are known to better than 1 MHz, and to 0.2–0.5 MHz where comb data reach both levels. Beyond it the uncertainty is 1–10 MHz for B v′ = 44–67 (comb data and the APO atlas), 10–30 MHz for X v″ = 18–60 (the Orsay atlas and the level constants of Martin *et al.*), and 0.1–0.2 GHz for the last X levels.</figcaption>
</figure>

## Typical values

| line | center | a15 − a1 | status |
|---|---:|---:|---|
| R(56) 32-0 (532 nm) | 0.42 MHz | 0.058 MHz | measured |
| R(150) 20-2 | 0.88 MHz | 0.058 MHz | measured |
| R(41) 51-0 | 1.4 MHz | 1.1 MHz | measured |
| R(98) 58-1 | 13 MHz | 2.9 MHz | measured |
| R(20) 70-0 | 14 MHz | 7.7 MHz | measured |
| R(50) 10-30 | 28 MHz | 0.34 MHz | measured |
| R(10) 2-90 | 91 MHz | 18 MHz | extrapolated |

*Table 3. Uncertainties of a few lines.*

## How well the uncertainties work

On the held-out residuals used to fit them, the uncertainties hold 94–97 % of the rows within 2σ in every
band beyond the core. Taking each held-out level's mean residual instead, 79 % are within 1σ and 95 % within
2σ, against 68 % and 95 % for a normal distribution. Core comb lines are over-covered (95 % within 2σ,
rms z = 0.71), because the statistical part carries the fit's 0.5 MHz floor on the comb data.

Thirty-six data sets transcribed on October 20, 2026 were kept out of every fit and compared with the
model's frozen predictions that day. The uncertainty model was built the day after, from the training fit and
its cross-validation, and not tuned on them. Its coverage of those rows is therefore a calibration
check, not a blind test:

| rows | n | within 1σ | within 2σ |
|---|---:|---:|---:|
| independent of the fit | 347 | 66 % | 91 % |
| absolute frequencies | 103 | 66 % | 88 % |
| hyperfine intervals beyond the core (B v′ = 51–62) | 48 | 52 % | 94 % |
| in the fit through compilations | 126 | 70 % | 95 % |

*Table 4. Coverage of the uncertainties on data the fit never used.*

One set, Nishiyama *et al.* (2013), lies well outside: a median deviation of 0.64 MHz against a median σ of
0.16 MHz.

A further 35 data sets were held out as a blind test of the uncertainties: the model's predictions and their
σ were frozen for every target first, and the measurements were then scored once:

| rows | n | within 1σ | within 2σ |
|---|---:|---:|---:|
| independent of the fit | 345 | 61 % | 90 % |
| absolute frequencies | 70 | 61 % | 100 % |
| intervals between lines | 76 | 93 % | 99 % |
| hyperfine intervals within a line | 199 | 49 % | 83 % |

*Table 5. The blind test of the uncertainties.*

The center uncertainties hold. The hyperfine part is too small: intervals within a line fall outside 2σ about
three times as often as they should.

## Above the fitted J

The fit's data reach J = 200, and at most levels much less. Beyond the J the training rows reach at a level's v
(or within 2 in v), plus 15, the level's model part grows as J(J + 1), the size of the centrifugal term through
which the rotational function and the shape of the potential act: σ<sub>model</sub> × J(J + 1)/J<sub>r</sub>(J<sub>r</sub> + 1). This
is a rule, not fitted. Such levels are marked *extrapolated*. The 321 lines of the Orsay atlas part I with
J″ = 200–271, which no fit used, scatter about the model by 50 MHz up to J″ = 225, the atlas's own scatter at lower J,
and by 65–100 MHz above.

## Isotopologues

The lines of ¹²⁹I₂ and ¹²⁷I¹²⁹I have the same three parts, read from the ¹²⁷I₂ tables at the equivalent
level (the same J = 0 energy and ρ J(J + 1)), with three additions:

- **Statistical part, conditioned.** The ¹²⁷I₂ levels do not separate the shape of a potential from its
  rotational function, which an isotopologue sees scaled by ρ: combinations of coefficients that leave every
  ¹²⁷I₂ level unchanged move the isotopologue levels as J(J + 1). The covariance of the fit and of the two
  isotope-shift coefficients is therefore conditioned on the seven CIPM isotopologue lines. It leaves
  0.09–0.31 MHz at those lines, and grows with the distance from them in v and J.
- **Isotope shift.** Each B level adds the leave-one-line-out rms of the shift fit, 0.078 MHz for ¹²⁹I₂ and
  half that for ¹²⁷I¹²⁹I, and the spread of the equally good forms of Eq. (5) of the
  [levels page](levels.md#isotopologues) (Table 4 there).
- **Hyperfine.** Each component adds a floor of 0.10 MHz, the excess of the CIPM intervals over their stated
  uncertainties.

A B level whose ¹²⁷I₂-equivalent v lies more than 2 from the lines that fix the shift (v′ = 6–14) is marked
*extrapolated*.

| line | center | a1 | a2 − a1 | status |
|---|---:|---:|---:|---|
| ¹²⁹I₂ P(54) 8-4 (633 nm) | 0.18 MHz | 0.21 MHz | 0.15 MHz | interpolated |
| ¹²⁹I₂ R(20) 15-1 | 0.74 MHz | 0.75 MHz | 0.15 MHz | interpolated |
| ¹²⁹I₂ R(56) 32-0 (532 nm) | 3.5 MHz | 3.5 MHz | 0.15 MHz | extrapolated |
| ¹²⁹I₂ R(41) 51-0 | 12 MHz | 12 MHz | 0.20 MHz | extrapolated |
| ¹²⁷I¹²⁹I P(33) 6-3 (633 nm) | 0.15 MHz | 0.19 MHz | 0.15 MHz | interpolated |
| ¹²⁷I¹²⁹I R(56) 32-0 | 1.8 MHz | 1.8 MHz | 0.15 MHz | extrapolated |

*Table 6. Uncertainties of a few isotopologue lines.*

Lines whose levels the model does not hold (B v′ > 80 for ¹²⁷I₂, or a level beyond its barrier rule) carry
the flag *outside the fitted domain* and an estimated uncertainty of at least 1 GHz.
