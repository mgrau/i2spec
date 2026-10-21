# Rovibrational levels

i2spec computes the levels of ¹²⁷I₂, ¹²⁹I₂ and ¹²⁷I¹²⁹I from one potential-energy curve per electronic state,
fitted to every measurement used, with no correction to any individual level. The other isotopologues see the
same curves at their own reduced mass, with a Born–Oppenheimer correction fitted to their measured isotope shifts
(see [Isotopologues](#isotopologues)).

## Radial Hamiltonian

The X and B states both have Ω = 0, so J is the rotational quantum number of the nuclear motion. For each
state the level energies E<sub>vJ</sub> and wavefunctions ψ<sub>vJ</sub>(R) are solutions of

$$
\left[-\frac{\hbar^2}{2\mu}\frac{d^2}{dR^2} + V(R) + \frac{\hbar^2\,[1+q(R)]}{2\mu R^2}\,J(J+1)\right]\psi_{vJ}(R) = E_{vJ}\,\psi_{vJ}(R), \tag{1}
$$

where R is the internuclear distance and μ = m(¹²⁷I)/2 = 63.452 236 u is the reduced mass of ¹²⁷I₂, from the
atomic mass of ¹²⁷I, 126.904 472 u. V(R) is the potential of the state. The rotational function q(R) is the
nonadiabatic correction to the centrifugal term:

$$
q(R) = \frac{2R_e}{R+R_e}\sum_{i} \alpha_i\,T_i\!\left(\frac{y-0.3}{0.4}\right),\qquad y=\frac{R-R_e}{R+R_e}, \tag{2}
$$

with T<sub>i</sub> the Chebyshev polynomials, which are defined on −1 … 1. The fixed shift 0.3 and scale 0.4 map
y = −0.1 … 0.7, that is R = 0.82 R<sub>e</sub> … 5.7 R<sub>e</sub>, onto that interval. They are part of the
definition of q(R), not fitted. The core levels sample y = −0.05 … 0.07 (X) and −0.07 … 0.21 (B), the lower end
of the interval. q(R) has 2 coefficients for X and 6 for B (Table 1). Beyond 3.08 Å (X) and 4.45 Å (B) it
continues as described [below](#continuation-to-both-dissociation-limits). ¹²⁷I₂ is the reference isotopologue:
for it there is no adiabatic correction term, and the other isotopologues' corrections are relative to it.

| i | α<sub>i</sub>, X ¹Σg⁺ | α<sub>i</sub>, B ³Π(0u⁺) |
|--:|--:|--:|
| 0 | 1.311 810 × 10⁻⁵ | −9.546 549 × 10⁻³ |
| 1 | 2.058 615 × 10⁻⁵ | −1.575 361 × 10⁻² |
| 2 | | −8.614 338 × 10⁻³ |
| 3 | | −3.017 776 × 10⁻³ |
| 4 | | −6.025 725 × 10⁻⁴ |
| 5 | | −6.359 487 × 10⁻⁵ |

*Table 1. Coefficients of the rotational function q(R), Eq. (2), from the core fit.*

## Potential-energy curves

<figure class="fig" markdown="span">
<div class="svg" data-svg="figures/potentials.svg" role="img" aria-label="Potential-energy curves of the X and B states with vibrational levels and the 532 nm transition"></div>
<figcaption markdown="span">**Figure 1.** Potential-energy curves of X¹Σg⁺ and B³Π(0u⁺) with their J = 0 vibrational levels, drawn between the classical turning points: every level of X up to v″ = 17 and every second level of B up to v′ = 43, then every third level of each up to the dissociation limits. Levels beyond the core fit are drawn faint. The arrow marks R(56) 32–0 at 532 nm. X dissociates to two ground-state atoms, I(²P3/2) + I(²P3/2); B dissociates to I(²P3/2) + I(²P1/2).</figcaption>
</figure>

### Morse/long-range form

Both potentials have the Morse/long-range (MLR) form of Le Roy & Henderson (2007):

$$
\begin{aligned}
V(R) &= T_e + D_e\left[1-\frac{u(R)}{u(R_e)}\,e^{-\beta(R)\,y_p^{e}(R)}\right]^2,\\
u(R) &= \sum_n\frac{C_n}{R^n},\qquad y_n^{r}(R)=\frac{R^n-r^n}{R^n+r^n},
\end{aligned} \tag{3}
$$

$$
\begin{aligned}
\beta(R) &= y_p^\text{ref}(R)\,\ln\frac{2D_e}{u(R_e)} + \left[1-y_p^\text{ref}(R)\right]\sum_i c_i\,T_i(z),\\
z &= \begin{cases} [y_q^\text{ref}(R)-0.2]/0.6 & \text{X,} \\ y_q^\text{ref}(R) & \text{B,}\end{cases}
\end{aligned} \tag{4}
$$

where y<sup>e</sup> uses r = R<sub>e</sub> and y<sup>ref</sup> uses r = R<sub>ref</sub>. For X the shift 0.2 and
scale 0.6 map y<sub>q</sub><sup>ref</sup> = −0.4 … 0.8 onto −1 … 1, as in Eq. (2). For any values of the fitted coefficients
c<sub>i</sub>, Eq. (3) gives V(R<sub>e</sub>) = T<sub>e</sub> and V(R) → T<sub>e</sub> + D<sub>e</sub> − u(R) at long range. The dispersion
coefficients C<sub>n</sub> are taken from the literature, so the long-range behavior does not depend on the
fit. Table 2 gives the settings and the main constants.

| | X ¹Σg⁺ | B ³Π(0u⁺) |
|---|---|---|
| p, q, R<sub>ref</sub> | 5, 3, R<sub>e</sub> | 6, 4, 1.5 R<sub>e</sub> |
| exponent coefficients c<sub>i</sub> | 14 | 20 |
| rotational coefficients α<sub>i</sub> | 2 | 6 |
| R<sub>e</sub> | 2.66638 Å | 3.02672 Å |
| D<sub>e</sub> | 12 547.34 cm⁻¹ | 4 381.25 cm⁻¹ |
| T<sub>e</sub> | 0 | 15 769.066 cm⁻¹ |
| C<sub>n</sub> (inherited), cm⁻¹ Åⁿ | C₆ 1.48 × 10⁶, C₈ 3.86 × 10⁷, C₁₀ 1.0 × 10⁸ | C₅ 3.161 × 10⁵, C₆ 1.506 × 10⁶, C₈ 2.48 × 10⁷, C₁₀ 4.2 × 10⁸ |

*Table 2. The MLR potentials of the core fit.*

These potentials, with the B hyperfine functions of the [hyperfine page](hyperfine.md), make up the **core
fit**. Its 61 coefficients (47 for the potentials, 14 for the B hyperfine functions)
were fitted in one least-squares problem to two kinds of data:

- every precision measurement of ¹²⁷I₂ with X v″ ≤ 17 and B v′ ≤ 43: frequency-comb and
  Doppler-free absolute frequencies of individual hyperfine components, and hyperfine intervals;
- the lines of two Fourier-transform atlases (Salami & Ross 2005; the NIST scan of the APO cell). Each
  atlas's scale was calibrated against the precision data with an offset and a slope.

The 29 precision data sets transcribed last were held out of the fit as a benchmark. The 515 nm reference line P(13) 43-0 is left out of every fit and score and kept as a
diagnostic. The core fit predicts the 60 held-out absolute component frequencies to 0.156 MHz rms and the
592 held-out hyperfine intervals to 0.059 MHz.

### Continuation to both dissociation limits

The core levels sample each potential only between their classical turning points: X 2.43–3.06 Å and
B 2.65–4.66 Å. The **continuation** keeps every coefficient of the core fit and continues each curve beyond
fixed knots:

- **Exponent and rotational tails.** Past a knot z<sub>k</sub>, the exponent sum of Eq. (4) and the rotational
  function of Eq. (2) are replaced by their own Taylor polynomial at the knot, plus new terms
  t<sup>k+1</sup>(1 − t)<sup>2</sup> B<sub>j</sub>(t). Here t ∈ [0, 1] is the distance past the knot in the bounded coordinate,
  and B<sub>j</sub> are Bernstein polynomials. The joins are smooth to second order for the X exponent and to
  first order for the B exponent and both rotational functions. The factor (1 − t)<sup>2</sup> makes every new
  term vanish at R → ∞, so the asymptote stays the inherited dispersion series.
    - X: 12 outer terms beyond 3.08 Å, 1 inner term inside 2.42 Å (to the solver's inner end, 2.1 Å) and
      1 rotational term beyond 3.08 Å.
    - B: 10 outer terms and 3 rotational terms beyond 4.45 Å.
- **A change of C₆ for X**, switched on beyond its knot as t³: +3.3 %.
- **A local correction to the X potential**: nine cubic B-splines on knots 4.2, 4.4, …, 6.6 Å, added to
  V(R). Beyond about 5 Å the MLR exponent barely reshapes V, and the levels near the X limit turn there.
- **One change of the dissociation energy**, shared by X and B so that their limits keep the atomic
  fine-structure splitting: +0.0135 cm⁻¹, switched on beyond the outer knots.

The new terms (38 free, with an offset of the Orsay Partie IV atlas) were fitted to 46 964 rows: every
training row of the core fit, which protects the core, together with the measurements beyond it:

| levels | data | held out, rms |
|---|---|---|
| X v″ = 18–25 | 3 700 lines of the Orsay atlas part I (Gerstenkorn, Vergès & Chevillard 1982) | 55 MHz |
| X v″ = 26–89 | level constants of Martin *et al.* (1986) | 44–122 MHz |
| X v″ = 48–54 | comb-referenced emission lines (Matyugin 2012, Nesterenko 2019) | 37 MHz |
| X v″ = 90–108 | levels of Martin *et al.* (1983) at v″ = 93, 97, 102, 103, 106 | 0.66 GHz |
| B v′ = 44–50 | comb-referenced lines (Matsunaga 2024, Yoshiki 2023, Ikeda 2022, …), APO and Salami & Ross atlases | 1.4 MHz (comb), 22–33 MHz (atlases) |
| B v′ = 51–62 | comb-referenced lines (Sakamoto 2024, Goncharov 2007, BIPM 2005), APO atlas | 10.5 MHz (comb, 3 rows), 38 MHz (atlas) |
| B v′ = 51–80 | 2 079 unblended lines of the Orsay atlas Partie IV (Gerstenkorn & Luc 1983) | 62–207 MHz |

*Table 3. Data beyond the core, and how well a level is predicted when all its rows are held out (four
folds by v, the final settings).*

The held-out column is the hardest test: every row of a level is removed and the level is predicted from
its neighbors. Within the core, the continued model reproduces the core fit. Every training row of the final
fit is reproduced to 1.4 kHz.

<figure class="fig" markdown="span">
<div class="svg" data-svg="figures/published_vs_mlr.svg" role="img" aria-label="Difference between the 2008 potentials and the model's for X and B"></div>
<figcaption markdown="span">**Figure 2.** The potentials of Salumbides *et al.* (2008) minus the model's, for X (left) and B (right). The shaded interval is the core, where the model's curves are the frozen core fit. The 2008 X curve was fitted only to v″ ≤ 17; beyond 3.1 Å it differs by up to 340 cm⁻¹, which moves its high levels by up to 20 cm⁻¹. Inside the core the two models give the same B levels (within 0.0003 cm⁻¹ at J = 0, v′ ≤ 43), yet their B potentials differ by up to 2 cm⁻¹: a potential and its rotational function q(R) are not determined separately by the levels, and each model's q(R) takes up part of the shape.</figcaption>
</figure>

## Numerical solution

Equation (1) is solved variationally in a basis of B-splines of order 10. The knots are uniform, with a
spacing of 0.01 Å (0.005 Å for X v″ ≥ 25), and each join of the potential (the knots above and the ends of
the X correction) is a knot of reduced continuity. Each range of v has its own radial box:

| levels | box |
|---|---|
| X v″ ≤ 17, 18–24, 25–33, 34–54, 55–69 | 2.1–4.2, 4.6, 5.0, 5.6 and 6.6 Å |
| X v″ ≥ 70 | 2.1–40.6 Å, graded mesh |
| B v′ ≤ 17, 18–33, 34–43, 44–59 | 2.35–4.6, 5.1, 6.1 and 8.6 Å |
| B v′ ≥ 60 | 2.35–40.6 Å, graded mesh |

Near the limits the J > 0 potentials have a centrifugal barrier. A level above the asymptote is kept as a
quasi-bound level when at least 90 % of its density lies inside the barrier. Halving the spacing and
lengthening the boxes changes the X levels of the training rows by at most 1.3 kHz, and every other
training row by less than 10 kHz. The one exception is B v′ = 44 at J′ = 186, a quasi-bound level whose
index changes on the refined grid.

The fitted data reach J = 200. Above it every level is solved in the longest box, where the same barrier rule
keeps the quasi-bound levels; the line lists reach J″ = 441. The levels of the other isotopologues use the box
of the ¹²⁷I₂ level of about the same energy. The line strengths use the wavefunctions of the same potentials, on a
grid suited to the transition moment (see [Line intensities and spectra](spectra.md)).

## No level corrections

The model has no level- or line-specific terms. Every level comes from the two potentials, and a
measurement of one line constrains all the levels through the shared coefficients. The uncertainty of each
level is therefore computed, not looked up: see [Uncertainty estimates](uncertainty.md).

## Isotopologues

¹²⁹I₂ (μ = 64.452 492 u) and ¹²⁷I¹²⁹I (μ = 63.948 453 u) are computed from the same two potentials. The breakdown
of the Born–Oppenheimer approximation enters in the form of the 2008 model, with ρ = μ(¹²⁷I₂)/μ:

- the rotational function of each state becomes ρ q(R);
- the B potential gains an adiabatic correction relative to X, which vanishes towards the separated atoms:

$$
\Delta V_B(R) = (1-\rho)\,(1-y)^{p}\,(a_0 + a_1\,y),\qquad y = \frac{R-R_e}{R+R_e}, \tag{5}
$$

with p = 5, R<sub>e</sub> the B equilibrium distance, a₀ = −0.1517 cm⁻¹ and a₁ = −1.761 cm⁻¹. 1 − ρ is 0.015 52 for
¹²⁹I₂ and 0.007 76 for ¹²⁷I¹²⁹I. No ¹²⁷I₂ coefficient changes.

a₀ and a₁ were fitted to the isotope shifts of the CIPM compilation (BIPM 2003): six ¹²⁹I₂ lines at 633 and
532 nm (v′ = 6–14, J″ = 33–113) and one ¹²⁷I¹²⁹I line at 633 nm, each tied through measured intervals to a measured
¹²⁷I₂ frequency, with the line centers taken from the model's own hyperfine patterns. Mass scaling alone puts the
¹²⁹I₂ lines 73–76 MHz above their measured positions and the ¹²⁷I¹²⁹I line 36 MHz above, half as much, as 1 − ρ
predicts. A constant shift predicts each line from the other six to 1.2 MHz rms; Eq. (6) does so to 0.08 MHz. The
power p was chosen by that leave-one-line-out rms among p = 0, 3, 5 and 8 (0.18, 0.13, 0.08 and 0.17 MHz).

| B v′ | 0 | 14 | 32 | 43 | 60 | 78 |
|---|---:|---:|---:|---:|---:|---:|
| ¹²⁹I₂, MHz | −70.8 | −75.2 | −79.5 | −78.7 | −63.6 | −19.5 |
| spread of the forms p = 3, 8 | 0.1 | 0.3 | 2.7 | 6.8 | 17.5 | 19.4 |

*Table 4. The isotope shift of B levels (J = 0) of ¹²⁹I₂, Eq. (5), and how far the equally good forms with p = 3
and 8 differ from it; ¹²⁷I¹²⁹I has half these values. The spread is part of each line's uncertainty.*

The ground-state levels that Cerny, Bacis & Vergès (1986) derived from fluorescence, v″ = 4–108, test the X curve
at the other masses. For ¹²⁹I₂ the term values G(v) − G(8) agree with the model within their stated
uncertainties (median |z| 0.3–1.3 by range of v). For ¹²⁷I¹²⁹I they lie 0.1–0.4 GHz below the model (median
|z| 1.2–3.0); an error of the potential would show twice as strongly in ¹²⁹I₂.

The model holds every level with X v″ ≤ 108 (109 for ¹²⁹I₂), B v′ ≤ 80 (81 for ¹²⁹I₂) and J ≤ 446 that is bound
or held behind its centrifugal barrier.
