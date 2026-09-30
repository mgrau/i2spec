# Rovibrational levels

## Radial Hamiltonian

The X and B states both have Ω = 0, so J is the rotational quantum number of the nuclear motion. For
each state the level energies E<sub>vJ</sub> and wavefunctions ψ<sub>vJ</sub>(R) are solutions of

$$
\left[-\frac{\hbar^2}{2\mu}\frac{d^2}{dR^2} + V(R) + V_\text{ad}(R) + \frac{\hbar^2\,[1+\alpha(R)]}{2\mu R^2}\,J(J+1)\right]\psi_{vJ}(R) = E_{vJ}\,\psi_{vJ}(R), \tag{1}
$$

where R is the internuclear distance, μ the reduced nuclear mass and V(R) the Born–Oppenheimer
potential. The functions V<sub>ad</sub>(R) and α(R) describe the adiabatic and nonadiabatic corrections
to the Born–Oppenheimer approximation (Salumbides *et al.* 2008):

$$
V_\text{ad}(R) = \left(1-\frac{\mu_\text{ref}}{\mu}\right)\left(\frac{2R_m}{R+R_m}\right)^{5}\sum_i v_i\,\xi^i,
\qquad
\alpha(R) = \frac{\mu_\text{ref}}{\mu}\,\frac{2R_m}{R+R_m}\sum_i \alpha_i\,\xi^i , \tag{2}
$$

with μ<sub>ref</sub> the reduced mass of ¹²⁷I₂ and ξ defined in Eq. (3). V<sub>ad</sub> vanishes for ¹²⁷I₂
and shifts the levels of the other isotopologues; α(R) modifies the effective rotational constant.
Both corrections are used for the B state only. The constant term v₀ of V<sub>ad</sub> was refitted in
this work to the BIPM tables of ¹²⁹I₂ and ¹²⁷I¹²⁹I, which reduced their rms deviation from 7.1 MHz to
2.4 MHz.

## Potential-energy curves

Two analytical forms are used. Each is used only for the levels where it has been tested against
measurements.

### Power-series potentials (Salumbides *et al.* 2008)

Within the range R<sub>I</sub> ≤ R ≤ R<sub>O</sub> the potential is a power series in the reduced coordinate ξ:

$$
V(R) = \sum_{i} a_i\,\xi^i,\qquad \xi = \frac{R-R_m}{R+b\,R_m}. \tag{3}
$$

For R < R<sub>I</sub> it is continued as A<sub>I</sub> exp[−B<sub>I</sub>(R − R<sub>I</sub>)], and for R > R<sub>O</sub> as
D<sub>e</sub> − Σ<sub>n</sub> C<sub>n</sub>/R<sup>n</sup> − A<sub>O</sub> exp[−B<sub>O</sub>(R − R<sub>O</sub>)]. The coefficients are those
of Salumbides *et al.* (2008): 14 coefficients a<sub>i</sub> for X and 32 for B. They were fitted to
data for X v″ ≤ 17 and B v′ ≤ 43, and within this range they reproduce the precision measurements to
a few MHz. They are referred to below as the 2008 potentials.

!!! note "Origin of the functional form"
    Expansions of a diatomic potential in a reduced coordinate that stays finite as R → ∞ predate their use
    for iodine. Simons, Parr & Finlan (1973) expanded in (R − R<sub>e</sub>)/R. Ogilvie (1981) introduced
    2(R − R<sub>e</sub>)/(R + R<sub>e</sub>) and, more generally, the family
    W<sub>mn</sub> = (m + n)(R − R<sub>e</sub>)/(mR + nR<sub>e</sub>) with integer m and n, which is ξ of Eq. (3) up to a
    constant factor, with b = n/m. Šurkus, Rakauskas & Bolotin (1984) generalized it to
    (R<sup>p</sup> − R<sub>e</sub><sup>p</sup>)/(R<sup>p</sup> + bR<sub>e</sub><sup>p</sup>); p = 1 is ξ with b a free real parameter.
    The model used here combines this series with an exponential inner wall and a dispersion-plus-exchange
    outer branch, joined smoothly at R<sub>I</sub> and R<sub>O</sub>. That construction is due to Tiemann and
    co-workers (for example Samuelis *et al.* 2000, for Na₂, where b is chosen to minimize the number of
    coefficients). For I₂ it was fitted by Knöckel, Bodermann & Tiemann (2004) and refitted, with the
    Born–Oppenheimer correction functions of Eq. (2), by Salumbides *et al.* (2008).

### Morse/long-range potentials

Outside the range of the fit the 2008 potentials deviate strongly from measurements: by
+7.8 cm⁻¹ at X v″ = 48, −21 cm⁻¹ at X v″ = 54 and about 2 cm⁻¹ at B v′ = 58. These levels have their
outer turning points where the potential is determined by the exponential term that joins the power
series to the long-range expansion, and this term is not constrained by the data. For these levels
i2spec uses Morse/long-range (MLR) potentials (Le Roy & Henderson 2007) fitted in this work:

$$
V(R) = D_e\left[1-\frac{u(R)}{u(R_e)}\,e^{-\beta(R)\,y_p^\text{eq}(R)}\right]^2,\qquad u(R)=\sum_n\frac{C_n}{R^n},\qquad y_p^{r}(R)=\frac{R^p-r^p}{R^p+r^p}, \tag{4}
$$

$$
\beta(R) = y_p^\text{ref}\,\beta_\infty + \left(1-y_p^\text{ref}\right)\sum_i \beta_i\,\left(y_q^\text{ref}\right)^i,\qquad \beta_\infty=\ln\frac{2D_e}{u(R_e)} . \tag{5}
$$

For any values of the fitted coefficients β<sub>i</sub>, Eq. (4) satisfies V(R<sub>e</sub>) = 0 and
V(R) → D<sub>e</sub> − u(R) at large R. The dissociation energy D<sub>e</sub> and the dispersion coefficients
C<sub>n</sub> therefore keep their measured values, and the long-range behaviour does not depend on the
fit.

The X potential is fitted to the levels v″ ≤ 17 of the 2008 potential, to the levels v″ = 48, 53 and 54 derived
from emission measurements, and to the levels v″ = 18–25 derived from the Orsay atlas. The B potential
is fitted to all precision measurements and to 2 009 lines of the Orsay atlas Partie IV (Gerstenkorn & Luc
1983), which reaches v′ = 79, 4 cm⁻¹ below the dissociation limit. Its centrifugal term carries one
fitted constant beyond 5 Å in addition to the published Born–Oppenheimer function, which was fitted
only inside the well. The potential alone reproduces the atlas to 30–35 MHz (median) at v′ = 51–66 and
80–100 MHz at v′ = 67–79; level corrections from the atlas bring every measured level to its data.

<figure class="fig" markdown="span">
<div class="svg" data-svg="figures/published_vs_mlr.svg" role="img" aria-label="Difference between the 2008 and MLR potentials for X and B"></div>
<figcaption markdown="span">**Figure 1.** Difference between the 2008 potentials (Salumbides *et al.*) and the MLR potentials for X (left) and B (right). The shaded interval is the range of turning points of the levels included in the 2008 fit. Within it the two potentials agree; outside it they differ by up to several hundred cm⁻¹ for X and several tens of cm⁻¹ for B.</figcaption>
</figure>

Table 1 lists the origin of the levels used by the model.

| levels | potential | tested against |
|---|---|---|
| X v″ = 0–17 | Salumbides *et al.* (2008) | frequency-comb measurements, 0.01–1 MHz |
| X v″ = 18–25 | MLR, with level corrections from the Orsay atlas | 3 700 atlas lines, 15–50 MHz per level |
| X v″ = 26–47 | MLR | levels of Martin *et al.* (1986), 0.1–0.3 GHz |
| X v″ = 48–54 | MLR | emission lines, 7–19 MHz |
| B v′ = 0–43 | Salumbides *et al.* (2008) | frequency-comb measurements |
| B v′ = 44–50 | MLR, with level corrections | frequency-comb measurements, 1–15 MHz |
| B v′ = 51–79 | MLR, with level corrections from the Orsay atlas Partie IV | 2 009 atlas lines and 3 comb lines, 20–250 MHz per level |
| B v′ = 80–86 | MLR | nothing measured |

*Table 1. Source of the level energies.*

## Numerical solution

Equation (1) is solved variationally in a basis of B-splines of order 10 with a knot spacing of about
0.01 Å. The 2008 potentials have a discontinuous second derivative at R<sub>I</sub> and R<sub>O</sub>; with a
uniform grid this limits the accuracy near these points to about 0.3 MHz. The B-spline basis
therefore has a knot of reduced continuity (C³) at each join, and the matrix elements are integrated
separately on each interval by Gauss–Legendre quadrature. The generalised eigenvalue problem
**Hc** = E**Sc** is then converged to about 0.1 kHz.

The size of the radial box is chosen for each vibrational level as the smallest box that reproduces
the energies of a larger box to better than 1 kHz. It ranges from 3.6 Å for X v″ ≤ 17 to 12 Å for
the highest B levels, which lie within a few cm⁻¹ of the dissociation limit. The wavefunctions needed for the transition moments are
computed with a sinc discrete-variable representation (Colbert & Miller 1992) on a common grid; their
energies are not used for line positions.

## Empirical level corrections

After the potentials are applied, the frequency-comb measurements still differ from the model by
1–6 MHz. The differences vary smoothly with J within a vibrational level and change from one level
to the next. They are therefore described as a correction to each level,

$$
\delta E(\text{state}, v, J) = \sum_{k=0}^{K} c_k\,y^k,\qquad y = J(J+1)/10^4, \tag{6}
$$

fitted to all precision measurements by robust regression with a ridge penalty. A correction changes
all lines that share the level. It is applied to ¹²⁷I₂ only, and only within 15 units of J of the
measured lines; outside that range the terms of first order in y are extrapolated and the
higher-order terms are held at their values at the boundary.

With the corrections the frequency-comb data sets are reproduced to 0.02–0.3 MHz, and lines left out
of the fit are predicted to 0.1–0.8 MHz.
