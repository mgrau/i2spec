# Rovibrational levels

i2spec computes the levels of ¹²⁷I₂, ¹²⁹I₂ and ¹²⁷I¹²⁹I from one potential-energy curve and one rotational function
per electronic state, with no correction to any individual level. The other isotopologues see the same curves at
their own reduced mass (see [Isotopologues](#isotopologues)). Which data fix which coefficients is described on
[The fit](fit.md); the derivations are in the paper.

## Radial Hamiltonian

The X and B states both have Ω = 0, so J is the rotational quantum number of the nuclear motion. For each state
the level energies E<sub>vJ</sub> and wavefunctions ψ<sub>vJ</sub>(R) are solutions of

$$
\left[-\frac{\hbar^2}{2\mu}\frac{d^2}{dR^2} + V(R) + \frac{\hbar^2\,[1+q(R)]}{2\mu R^2}\,J(J+1)\right]\psi_{vJ}(R) = E_{vJ}\,\psi_{vJ}(R), \tag{1}
$$

where R is the internuclear distance, μ = m(¹²⁷I)/2 = 63.452 236 u the reduced mass of ¹²⁷I₂ and V(R) the
potential of the state. The rotational function q(R) is the nonadiabatic correction to the centrifugal term,

$$
q(R) = \frac{2R_e}{R+R_e}\sum_{i} \alpha_i\,T_i\!\left(\frac{y-0.3}{0.4}\right),\qquad y=\frac{R-R_e}{R+R_e}, \tag{2}
$$

with T<sub>i</sub> the Chebyshev polynomials; the fixed shift 0.3 and scale 0.4 are part of the definition. q(R) has
{{ fit.x_alpha }} coefficients for X and {{ fit.b_alpha }} for B. ¹²⁷I₂ is the reference isotopologue and has no
adiabatic correction. Energies are measured from the minimum of the X potential.

## Potential-energy curves

<figure class="fig" markdown="span">
<div class="svg" data-svg="figures/potentials.svg" role="img" aria-label="Potential-energy curves of the X and B states with vibrational levels and the 532 nm transition"></div>
<figcaption markdown="span">**Figure 1.** Potential-energy curves of X¹Σg⁺ and B³Π(0u⁺) with their J = 0 vibrational levels: every level of X up to v″ = 17 and every second level of B up to v′ = 43, then every third level of each up to the dissociation limits (drawn faint beyond the core). The arrow marks R(56) 32–0 at 532 nm. X dissociates to I(²P3/2) + I(²P3/2), B to I(²P3/2) + I(²P1/2).</figcaption>
</figure>

Both potentials have the Morse/long-range (MLR) form of Le Roy & Henderson (2007):

$$
\begin{aligned}
V(R) &= T_e + D_e\left[1-\frac{u(R)}{u(R_e)}\,e^{-\beta(R)\,y_p^{e}(R)}\right]^2,\qquad
u(R) = \sum_n\frac{C_n}{R^n},\qquad y_n^{r}(R)=\frac{R^n-r^n}{R^n+r^n},\\
\beta(R) &= y_p^\text{ref}(R)\,\ln\frac{2D_e}{u(R_e)} + \left[1-y_p^\text{ref}(R)\right]\sum_i c_i\,T_i(z),\qquad
z = \begin{cases} [y_q^\text{ref}(R)-0.2]/0.6 & \text{X,} \\ y_q^\text{ref}(R) & \text{B,}\end{cases}
\end{aligned} \tag{3}
$$

where y<sup>e</sup> uses r = R<sub>e</sub> and y<sup>ref</sup> uses r = R<sub>ref</sub>. Whatever the fitted c<sub>i</sub>,
V(R<sub>e</sub>) = T<sub>e</sub> and V(R) → T<sub>e</sub> + D<sub>e</sub> − u(R) at long range, so the long-range
behavior is the inherited dispersion series. The B potential has no free dissociation energy: its limit is that of
X plus the atomic fine-structure splitting, 7 602.9762 cm⁻¹, exactly.

| | X ¹Σg⁺ | B ³Π(0u⁺) |
|---|---|---|
| p, q, R<sub>ref</sub> | 5, 3, R<sub>e</sub> | 6, 4, 1.5 R<sub>e</sub> |
| exponent coefficients c<sub>i</sub> | {{ fit.x_c }} | {{ fit.b_c }} |
| rotational coefficients α<sub>i</sub> | {{ fit.x_alpha }} | {{ fit.b_alpha }} |
| R<sub>e</sub> | {{ we.x_re }} Å | {{ we.b_re }} Å |
| D<sub>e</sub> | {{ we.x_de }} cm⁻¹ | {{ we.b_de }} cm⁻¹ (from the limit of X) |
| T<sub>e</sub> | 0 | {{ we.b_te }} cm⁻¹ |
| C<sub>n</sub> (inherited), cm⁻¹ Åⁿ | C₆ 1.48 × 10⁶, C₈ 3.86 × 10⁷, C₁₀ 1.0 × 10⁸ | C₅ 3.161 × 10⁵, C₆ 1.506 × 10⁶, C₈ 2.48 × 10⁷, C₁₀ 4.2 × 10⁸ |

*Table 1. The MLR potentials inside the core (X v″ ≤ 17, B v′ ≤ 43). Every coefficient is returned by
`ResearchModel("i2spec").coefficients`.*

### Beyond the core

The core levels sample X only at 2.43–3.06 Å and B at 2.65–4.66 Å. Beyond fixed knots (3.08 Å for X, 4.45 Å for B)
each curve is continued to its dissociation limit, keeping every coefficient of the core:

- the exponent sum and the rotational function are replaced, past their knot, by their own Taylor polynomial at
  the knot plus new terms that vanish at R → ∞, so the asymptote stays the dispersion series
  (X: {{ fit.x_outer }} outer, {{ fit.x_inner }} inner and {{ fit.x_rot }} rotational term; B: {{ fit.b_outer }}
  outer and {{ fit.b_rot }} rotational terms);
- for X, a change of C₆ switched on beyond the knot, and {{ fit.x_splines }} cubic B-splines on 4.2–6.6 Å added to
  V(R), where the levels near the X limit turn;
- one change of the dissociation energy, shared by X and B so that their limits keep the atomic splitting.

<figure class="fig" markdown="span">
<div class="svg" data-svg="figures/published_vs_mlr.svg" role="img" aria-label="Difference between the 2008 potentials and the model's for X and B"></div>
<figcaption markdown="span">**Figure 2.** The potentials of Salumbides *et al.* (2008) minus the model's (shaded: the core). Inside the core both give the same B levels, yet their B curves differ by up to 2 cm⁻¹: the ¹²⁷I₂ levels do not separate a potential from its rotational function.</figcaption>
</figure>

## Numerical solution

Equation (1) is solved variationally in a basis of B-splines of order 10 on uniform knots (0.01 Å; 0.005 Å for
X v″ ≥ 25), with a knot of reduced continuity at every join of the potential, in a radial box for each range of v
(up to 40.6 Å, graded, near the limits). Above the asymptote a level is kept as quasi-bound when at least 90 % of
its density lies inside the centrifugal barrier. Halving the knot spacing moves no training line by more than
10 kHz, except one quasi-bound level whose index changes. Since every level comes from the two potentials, its
uncertainty is computed, not looked up (see [Uncertainty estimates](uncertainty.md)).

## Isotopologues

¹²⁹I₂ (μ = 64.452 492 u) and ¹²⁷I¹²⁹I (μ = 63.948 453 u) use the same two potentials. With ρ = μ(¹²⁷I₂)/μ, the
rotational function of each state becomes ρ q(R), and the B potential gains an adiabatic correction relative to X
that vanishes towards the separated atoms:

$$
\Delta V_B(R) = (1-\rho)\,(1-y)^{p}\,(a_0 + a_1\,y),\qquad y = \frac{R-R_e}{R+R_e}, \tag{4}
$$

with p = 5, R<sub>e</sub> that of B, a₀ = {{ we.iso_a0 }} cm⁻¹ and a₁ = {{ we.iso_a1 }} cm⁻¹, fitted to
{{ fit.n_iso_centers }} measured isotopologue line centers (see [The fit](fit.md#isotopologues)); 1 − ρ is 0.015 52
for ¹²⁹I₂ and 0.007 76 for ¹²⁷I¹²⁹I. Mass scaling alone would put the ¹²⁹I₂ lines 73–76 MHz above their measured
positions; Eq. (4) predicts each line from the others to 0.08 MHz rms. Equally good forms (p = 3 and 8) differ from
it by up to 25 MHz near the B limit, and that spread is part of each isotopologue line's uncertainty. The
fluorescence term values of Cerny, Bacis & Vergès (1986), v″ = 4–108, agree with the model within their
uncertainties for ¹²⁹I₂; for ¹²⁷I¹²⁹I they lie 0.1–0.4 GHz above it.

The model holds every level with X v″ ≤ 108 (109 for ¹²⁹I₂), B v′ ≤ 80 (81 for ¹²⁹I₂) and J ≤ 446 that is bound
or held behind its centrifugal barrier.

## Worked example

!!! example "R(56) 32–0: the two levels"
    The 532 nm line joins X(v″ = 0, J″ = 56) to B(v′ = 32, J′ = 57). Both lie inside the core.

    | | X v″ = 0 | B v′ = 32 |
    |---|--:|--:|
    | E(v, J = 0), cm⁻¹ | {{ we.ex0 }} | {{ we.eb0 }} |
    | classical turning points at J = 0 | {{ we.x_turn }} Å | {{ we.b_turn }} Å |
    | B<sub>v</sub> = [E(v, 2) − E(v, 0)]/6, cm⁻¹ | {{ we.x_bv }} | {{ we.b_bv }} |
    | E(v, J) at the line's J, cm⁻¹ | {{ we.ex }} (J″ = 56) | {{ we.eb }} (J′ = 57) |
    | rotational energy E(v, J) − E(v, 0), cm⁻¹ | {{ we.ex_rot }} | {{ we.eb_rot }} |

    The lower level lies E″ = {{ we.e_lower }} cm⁻¹ above X(0, 0), which sets its population on the
    [spectra page](spectra.md#worked-example). Its difference from the upper level is the
    [line position](lines.md#worked-example).

    ```python
    from i2spec import RovibronicModel

    levels = RovibronicModel("127I2")   # the released model's levels, cm⁻¹ above the X minimum
    levels.energy("X", 0, 56)           # {{ we.ex_raw }}
    levels.energy("B", 32, 57)          # {{ we.eb_raw }}
    levels.potentials["B"](3.0)         # V_B at 3.0 Å, cm⁻¹
    ```
