# Rovibrational levels

i2spec computes the levels of ¹²⁷I₂ with **balanced-expanded**. It has one potential-energy curve per
electronic state, fitted to every measurement used, with no correction to any individual level. The other
isotopologues, and the few ¹²⁷I₂ lines outside balanced-expanded's domain, use the earlier model,
i2spec2026o, described [at the end of this page](#the-earlier-model-i2spec2026o).

## Radial Hamiltonian

The X and B states both have Ω = 0, so J is the rotational quantum number of the nuclear motion. For each
state the level energies E<sub>vJ</sub> and wavefunctions ψ<sub>vJ</sub>(R) are solutions of

$$
\left[-\frac{\hbar^2}{2\mu}\frac{d^2}{dR^2} + V(R) + \frac{\hbar^2\,[1+q(R)]}{2\mu R^2}\,J(J+1)\right]\psi_{vJ}(R) = E_{vJ}\,\psi_{vJ}(R), \tag{1}
$$

where R is the internuclear distance and μ the reduced nuclear mass of ¹²⁷I₂. V(R) is the potential of the
state. The rotational function q(R) is the nonadiabatic correction to the centrifugal term:

$$
q(R) = \frac{2R_e}{R+R_e}\sum_{i} \alpha_i\,T_i\!\left(\frac{y-0.3}{0.4}\right),\qquad y=\frac{R-R_e}{R+R_e}, \tag{2}
$$

with T<sub>i</sub> the Chebyshev polynomials. It has 2 coefficients for X and 6 for B. ¹²⁷I₂ is the reference
isotopologue, so there is no adiabatic correction term.

## Potential-energy curves

### Morse/long-range form

Both potentials have the Morse/long-range (MLR) form of Le Roy & Henderson (2007):

$$
V(R) = T_e + D_e\left[1-\frac{u(R)}{u(R_e)}\,e^{-\beta(R)\,y_p^{e}(R)}\right]^2,\qquad u(R)=\sum_n\frac{C_n}{R^n},\qquad y_n^{r}(R)=\frac{R^n-r^n}{R^n+r^n}, \tag{3}
$$

$$
\beta(R) = y_p^\text{ref}(R)\,\ln\frac{2D_e}{u(R_e)} + \left[1-y_p^\text{ref}(R)\right]\sum_i c_i\,T_i(z),\qquad z = \begin{cases} [y_q^\text{ref}(R)-0.2]/0.6 & \text{X} \\ y_q^\text{ref}(R) & \text{B,}\end{cases} \tag{4}
$$

where y<sup>e</sup> uses r = R<sub>e</sub> and y<sup>ref</sup> uses r = R<sub>ref</sub>. For any values of the fitted coefficients
c<sub>i</sub>, Eq. (3) gives V(R<sub>e</sub>) = T<sub>e</sub> and V(R) → T<sub>e</sub> + D<sub>e</sub> − u(R) at long range. The dispersion
coefficients C<sub>n</sub> are taken from the literature, so the long-range behavior does not depend on the
fit. Table 1 gives the settings and the main constants.

| | X ¹Σg⁺ | B ³Π(0u⁺) |
|---|---|---|
| p, q, R<sub>ref</sub> | 5, 3, R<sub>e</sub> | 6, 4, 1.5 R<sub>e</sub> |
| exponent coefficients c<sub>i</sub> | 14 | 20 |
| rotational coefficients α<sub>i</sub> | 2 | 6 |
| R<sub>e</sub> | 2.66638 Å | 3.02672 Å |
| D<sub>e</sub> | 12 547.34 cm⁻¹ | 4 381.25 cm⁻¹ |
| T<sub>e</sub> | 0 | 15 769.066 cm⁻¹ |
| C<sub>n</sub> (inherited), cm⁻¹ Åⁿ | C₆ 1.48 × 10⁶, C₈ 3.86 × 10⁷, C₁₀ 1.0 × 10⁸ | C₅ 3.161 × 10⁵, C₆ 1.506 × 10⁶, C₈ 2.48 × 10⁷, C₁₀ 4.2 × 10⁸ |

*Table 1. The MLR potentials of the core fit, research-balanced (release 20261018).*

These potentials, with the B hyperfine functions of the [hyperfine page](hyperfine.md), make up the **core
fit**, research-balanced. Its 61 coefficients (47 for the potentials, 14 for the B hyperfine functions)
were fitted in one least-squares problem to two kinds of data:

- every precision measurement of ¹²⁷I₂ with X v″ ≤ 17 and B v′ ≤ 43: frequency-comb and
  Doppler-free absolute frequencies of individual hyperfine components, and hyperfine intervals;
- the lines of two Fourier-transform atlases (Salami & Ross 2005; the NIST scan of the APO cell). Each
  atlas's scale was calibrated against the precision data with an offset and a slope.

The 29 precision data sets transcribed after the parameter set i2spec2026k was fixed were held out as a
benchmark. The 515 nm reference line P(13) 43-0 is left out of every fit and score and kept as a
diagnostic. The core fit predicts the 60 held-out absolute component frequencies to 0.156 MHz rms and the
592 held-out hyperfine intervals to 0.059 MHz.

### Continuation to both dissociation limits

The core levels sample each potential only between their classical turning points: X 2.43–3.06 Å and
B 2.65–4.66 Å. **balanced-expanded** keeps every coefficient of the core fit and continues each curve beyond
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

*Table 2. Data beyond the core, and how well a level is predicted when all its rows are held out (four
folds by v, the final settings).*

The held-out column is the hardest test: every row of a level is removed and the level is predicted from
its neighbors. Within the core, balanced-expanded reproduces the core fit. Every training row of the final
fit is reproduced to 1.4 kHz.

<figure class="fig" markdown="span">
<div class="svg" data-svg="figures/published_vs_mlr.svg" role="img" aria-label="Difference between the 2008 potentials and balanced-expanded's for X and B"></div>
<figcaption markdown="span">**Figure 1.** The potentials of Salumbides *et al.* (2008), the basis of the earlier model, minus balanced-expanded's, for X (left) and B (right). The shaded interval is the core, where balanced-expanded's curves are the frozen core fit. The 2008 X curve was fitted only to v″ ≤ 17; beyond 3.1 Å it differs by up to 340 cm⁻¹, which moves its high levels by up to 20 cm⁻¹. Inside the core the two models give the same B levels (within 0.0003 cm⁻¹ at J = 0, v′ ≤ 43), yet their B potentials differ by up to 2 cm⁻¹: a potential and its rotational function q(R) are not determined separately by the levels, and each model's q(R) takes up part of the shape.</figcaption>
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

The wavefunctions for the line strengths are computed with the earlier model (see
[Line intensities and spectra](spectra.md)).

## No level corrections

balanced-expanded has no level- or line-specific terms. Every level comes from the two potentials, and a
measurement of one line constrains all the levels through the shared coefficients. The uncertainty of each
level is therefore computed, not looked up: see [Uncertainty estimates](uncertainty.md).

## The earlier model, i2spec2026o

i2spec2026o gives the levels of ¹²⁹I₂ and ¹²⁷I¹²⁹I. It also gives the ¹²⁷I₂ levels that balanced-expanded does
not hold: J″ > 200, B v′ > 80, and levels beyond the dissociation limit at their J. It combines three
pieces:

- **The 2008 potentials** of Salumbides *et al.* (2008), in the range they were fitted to (X v″ ≤ 17,
  B v′ ≤ 43). They are power series in ξ = (R − R<sub>m</sub>)/(R + bR<sub>m</sub>) with exponential inner walls and
  dispersion-plus-exchange outer branches, with Born–Oppenheimer correction functions. These give the
  isotope shifts. Their form goes back to Ogilvie (1981), Šurkus *et al.* (1984) and Samuelis *et al.*
  (2000); Knöckel, Bodermann & Tiemann (2004) fitted it to I₂.
- **MLR potentials** fitted beyond that range.
- **Empirical corrections to individual levels**, polynomials in J(J + 1), fitted to the comb lines, with
  a Gaussian process for the B levels that have no data of their own.
