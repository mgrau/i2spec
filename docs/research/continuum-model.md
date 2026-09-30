# Continuum absorption of I₂: A←X, C←X and B←X bound–free (M3)

*Read from the full texts, 2026-09-15. The numbers were taken from the `pdftotext` layer and checked against the page images. Equations and tables are cited by paper and page.*

**Sources.**
* **[T73]** J. Tellinghuisen, J. Chem. Phys. 58, 2821 (1973): "Resolution of the visible-infrared absorption spectrum of I₂ into three contributing transitions".
* **[T11C]** J. Tellinghuisen, J. Chem. Phys. 135, 054301 (2011): least-squares analysis of overlapped bound–free spectra and predissociation data; the C(1u, ¹Π) state.
* **[T11B]** J. Tellinghuisen, J. Chem. Phys. 134, 084301 (2011): B–X transition moment and ε_c at 520–635 nm.
* **[T82]** J. Tellinghuisen, J. Chem. Phys. 76, 4736 (1982). **Not available locally**; it is only quoted here as the 2011 papers cite it.

Page numbers are journal pages (for example p. 2826, or p. 054301-9). Note that "the C state" in the 2011 papers is the 1u(¹Π) state, which is labelled "1u(¹Π)" in [T73] and "B′" in older literature ([T11C] ref. 2).

---

## 0. Bottom line

1. **Use the [T11C] model.**
   * ε_ν = 108.861 · ν · G_ab · ⟨|⟨ε′J″|μ_e(R)|υ″J″⟩|²⟩_Boltzmann, with the continuum energy-normalized per cm⁻¹ (eq. 2).
   * The potentials are exponential polynomials (eq. 4), and the transition moments are linear or quadratic in R − 2.7 Å (Table I).
   * **Use G_ab = 1 with the tabulated μ.** Table I footnote (e): the μ are for "nondegenerate" transitions; if G = 2, divide by √2.
   * This model supersedes [T73] and [T82]. **[T82] is not needed.**
2. **Parameters** (§2.2; energies in cm⁻¹, R in Å, μ in D):

   | State | Energy zero | R₀ | A₀ | B₀ | aᵢ | μ(R) |
   |---|---|---|---|---|---|---|
   | A | T_e(A) = D_e(X) − 1639.9 | 2.7 | −226.70 | 3401.305 | a₁ = −6.36(17), a₂ = −4.0(1.1) | 0.2845(11) − 0.048(33) z |
   | C | D_e(X) | 2.83 | −211.0(2.5) | 4599.3(2.5) | a₁ = −3.297(23), a₂ = 0.505(49), a₃ = −0.46(11), a₄ = 0, a₅ = −0.33(7) | 0.4714(32) − 0.16(11) z |
   | B (inner wall) | T_e(B) | 2.7 | −1783.226 | 4408.223 | a₁ = −4.63(11) | 1.0055(12) + 0.72(3) z − 0.17(42) z² |

   In μ(R), z = R − 2.7 for all three states.
3. **Verified numerically** (§7). Implementing this with i2spec's Hannover X reproduces [T11C] Table II:
   * to **0.3–0.7% at 400–500 nm** (0 °C and 35 °C);
   * to ≤ 1.1% at 650–800 nm;
   * to 1–2% elsewhere.
4. **B←X bound–free: compute it, don't parameterize** (§9.2).
   * Use Hannover B for R ≥ 2.7154 Å and the [T11C] exponential inner wall below. **The Hannover inner extension (R < 2.647 Å) is about 1100 cm⁻¹ too soft at 2.50 Å and is unusable for the continuum.**
   * Use energy-normalized Numerov continuum functions; a box-discretized continuum is only a cross-check.
   * Close the line-list gap near the dissociation limit (§9.4).
5. **Validation targets** (§6): [T11C] Table II (390–900 nm at 0 and 35 °C; λ ≤ 495 nm is truly continuous), plus:
   * ε(500.2 nm, 35 °C) = 587.4(8);
   * the point where ε is independent of T: 486.1 nm, ε ≈ 408, from Table II. The text of [T11C] says 487 nm and 416 (§7.7);
   * ε(436 nm, room T) = 31.0 ± 0.4;
   * [T11B] ε_c at 530–635 nm, with the cell at about 309 K (not "room temperature").
6. **Missing** (§8):
   * the [T11C] supplementary component bands;
   * the A-state RKR wall between 2.806 Å and R_e(A), needed for λ ≳ 810 nm;
   * any validation data above 393 K. At 700 K the far wings are uncertain by about 10–20%, mostly from the poorly determined slopes of μ_C and μ_A.
7. **Typos in [T11C]:**
   * "3321.8 cm⁻¹ above the X-state limit" should be 3221.8 (from Hannover, T_e(B) − D_e(X) = 3221.73).
   * "1 D = 3.3456 C m" should be 3.33564×10⁻³⁰ C m.

---

## 1. Units and conversions

| Quantity | Relation | Source |
|---|---|---|
| Decadic molar absorptivity ε (L mol⁻¹ cm⁻¹) | A = ε l C, with l in cm and C in mol/L | [T11B] eq. (2) |
| Cross section | σ = 3.8235×10⁻²¹ ε cm² | from [T73] eq. (12): k_ν/N = 2.303×10³ ε_ν/N₀ (N₀ = Avogadro's number) |
| Oscillator strength | f = 4.32×10⁻⁹ ∫ε dν (ν in cm⁻¹) | [T73] footnote 34 |

---

## 2. [T11C]: the current model (use this one)

### 2.1 Cross section (p. 054301-4)

**Eq. (1)**, the cross section from level (υ″, J″) to the continuum level at energy ε′:

  σ_ν,ab(υ″, J″) = (2π²ν / 3ε₀hc) · G_ab · |⟨ε′J″| μ_e(R) |υ″J″⟩|²

* The sum over the R, P and (Q) branches is replaced by a single Q-branch transition, J′ = J″, which may be real or hypothetical.
* G_ab is the electronic degeneracy factor. **G_ab = 1 for B←X** (0←0 in Hund's case c) and **G_ab = 2 for A←X and C←X** (both 1←0).

**Eq. (2)**, the molar absorptivity:

  **ε_ν = 108.861 · ν · G_ab · [ |⟨ε′J″| μ_e(R) |υ″J″⟩|² ]_avg**

* ν is in cm⁻¹ and μ_e in debye.
* For this constant to give L mol⁻¹ cm⁻¹, the continuum function |ε′⟩ must be **energy-normalized per cm⁻¹**. Check: 2π²/(3ε₀hc) · (1 D)² · (1 cm⁻¹) · N_A/(1000 ln 10) = 108.86. This is the same constant as the bound–bound 108.862 in [T11B] eq. (3), with the unit-area line shape S_L(ν) replaced by the continuum density.
* The paper prints "1D = 3.3456 C m". That is a typo for 3.33564×10⁻³⁰ C m; [T11B] prints 3.3356.

**Averaging and numerics** (p. 054301-4):
* "[…] averaged over the Boltzmann distribution of absorbing υ″ and J″ levels."
* "[…] the rotational averaging can be accomplished with negligible error using 5 appropriately chosen J″ levels, and with very little error using just the average J″."
* "Nuclear-spin degeneracy effects average out and can be ignored."
* "I have included the first 10 υ″ levels, which covers > 99.9% of the population at 120 °C." At 700 K this is **not** enough; see §7.
* The bound and free wavefunctions come from "standard numerical methods that may be considered exact".
* "The potential functions used here are effective potentials including the centrifugal term proportional to J(J+1)/R²."
* For X, [T11C] used the RKR curve from the constants of Martin et al. (1986). It notes that the quantal X potential of Salumbides et al. (2008) differs from it by only about 10⁻⁵ Å in the low-υ region (p. 054301-5). So using our Hannover X potential is consistent.
* Footnote 31: "The computed bound-free absorption is insensitive to the actual dissociation limit, provided the potential curve is smooth and the continuum wavefunctions are properly energy normalized."
* Fig. 6 compares the average-J analysis with the 5-J analysis.

### 2.2 Potentials (eq. 4, p. 054301-4; Table I, p. 054301-9)

**Eq. (4)**, used for all three upper states at small R:

  U(z) = A₀ + B₀ exp(a₁z + a₂z² + …), with z ≡ R − R₀

R₀ is 2.7 Å for A and B, and 2.83 Å for C (Table I, footnote a).

**Table I** (standard errors in the last digits in parentheses; units cm⁻¹, Å, D):

| Parameter | A | B | C |
|---|---|---|---|
| A₀ | −226.70 | −1783.226 | −211.0 (2.5) |
| B₀ | 3401.305 | 4408.223 | 4599.3 (2.5) |
| a₁ | −6.36 (17) | −4.63 (11) | −3.297 (23) |
| a₂ | −4.0 (1.1) | | 0.505 (49) |
| a₃ | | | −0.46 (11) |
| a₅ | | | −0.33 (7) |
| μ₀ | 0.2845 (11) | 1.0055 (12) | 0.4714 (32) |
| μ₁ | −0.048 (33) | 0.72 (3) | −0.16 (11) |
| μ₂ | | −0.17 (42) | |
| χ²_ν (whole fit) | | 1.186 | |

**Table I footnotes, verbatim or nearly so:**
* **(a)** "Potential curve parameters as defined in Eq. (4), with R₀ = 2.7 for the A and B states, 2.83 for C; transition moment functions defined as **|μ_e| = μ₀ + μ₁z + μ₂z², with z = (R − 2.7)**."
  * For the moments, z = R − 2.7 Å for **all three** transitions, C included, even though the C potential uses R₀ = 2.83.
* **(b)** "Standard errors in parentheses, in terms of final digits. Units cm⁻¹, Å, and debye, with **energies defined relative to state T_e values (A and B) or dissociation asymptote (C)**. Other parameters: wavelength corrections (additive, nm) — −0.14(7) for present data and 0.80(9) for Ref. 26; intensity scale factors — 1.00 for present data at 35 °C and 0.99 for Ref. 26 (both fixed), and 0.9963(8) for 64 °C. All parameters have been rounded systematically to preserve accuracy."
* **(c)**, attached to the B column: the predissociative rate coefficient P in eq. (3) is 3.22(3)×10¹⁰ s⁻¹ T⁻², for B→C magnetic predissociation (T = tesla).
* **(d)** For C, "Note absence of fourth-order term, which was dropped for statistical insignificance; the original fitted value was −0.07(94)."
* **(e)** "Transition moments fitted for assumed nondegenerate transitions. With allowance for electronic degeneracy (G_ab in Eq. (1)), the values of μ₀ and μ₁ for the A−X and C−X systems should be divided by √2."
  * **So either** use G_ab = 1 with the tabulated μ, **or** G_ab = 2 with μ/√2. The two are identical. **Do not** combine G_ab = 2 with the tabulated μ, which would double A←X and C←X.

**How each state is built.**
* **A state.**
  * Eq. (4) is used only at small R. It is attached to an RKR curve based on [Tellinghuisen, J. Chem. Phys. 118, 3532 (2003)], used up to υ = 29 on the repulsive branch.
  * The attachment points are given in **eq. (5)**, with energies relative to the A minimum:
    * U(2.804259 Å) = 1451.281 cm⁻¹
    * U(2.805916 Å) = 1431.379 cm⁻¹
  * A₀ and B₀ are set to make a smooth join; only a₁ and a₂ were fitted.
* **B state.**
  * Eq. (4) (one exponent) is attached to an RKR curve from Gerstenkorn & Luc constants (1985), smoothed below 2.63 Å.
  * The recommended fit (Table I) attaches at **υ = 33** (§IV C and D; Fig. 9). **The attachment R and energy are not printed.**
  * The χ² falls by 4.5% as the attachment moves from υ = 46 down to about 35, then stabilizes.
  * The IPA and Salumbides 2008 quantal B potentials were judged less suitable for this purpose (p. 054301-5):
    * The quantal potential "behaves anomalously on its repulsive branch above the high-υ limit of the data (υ = 43) […] it is displaced to small R by >0.0005 Å over the region of importance for absorption".
    * "an R shift of 0.0005 Å corresponds to an energy shift of 10–15 cm⁻¹, which represents a spectral shift of 0.3 nm below 500 nm."
* **Pseudocontinuum model for A and B** (p. 054301-4):
  * "Beyond the attachment points, the required potential energy points are obtained by interpolating on the A and B RKR curves, with **energy set to zero for R > R_e** (the pseudocontinuum model)."
  * This gives artificial dissociation limits "respectively 1639.9 cm⁻¹ below and 3321.8 cm⁻¹ above the X-state limit".
  * The reflection-type Franck–Condon intensity "[is] determined entirely by the shape of that wall and [is] completely insensitive to the nature of the attractive branch" (p. 054301-2).
  * The whole of B←X, discrete plus continuous, therefore becomes one smooth band.
* **C state.**
  * A₀ and B₀ are fitted along with the aᵢ.
  * The fit is constrained at intermediate R by the B→C magnetic predissociation rates (Vigué et al. 1981, Table VII) and at large R by 8 left-branch RKR turning points of the van der Waals well, υ = 0–7 (Inard et al. 1999), each with weight 2.7 (σ ≈ 0.6 cm⁻¹).
  * "The C curve is similarly extrapolated to zero at large R, but in this case by attaching an R⁻ⁿ segment beyond the last fitted RKR point (since Eq. (4) does not go to zero at large z)." **Neither n nor the attachment point is given.**
  * The C potential is most precise, about 1 cm⁻¹, at 2.83 Å. That is 0.07 Å below the B/C crossing, and "the B/C crossing occurs just above υ_B = 2" (Fig. 7 caption).

**The transition moments** (§IV C; Fig. 11; conclusion):
* A←X and C←X were fitted as linear in R, and B←X as quadratic.
* For B←X, μ₂ is statistically significant only for the υ = 46 attachment. "[It] is adequately linear for the statistically best fit".
* "The slope of the C ← X moment function changes sign for the same change in potentials". The negative slope agrees with theory [Zaitsevskii et al. 2000] but "is only marginally determined […] and smaller in magnitude than estimated from theory."

### 2.3 Data, fitted range and stated reliability

**Data used in the fit** (§II):
* Own spectra of ε at **35 °C and 64 °C**. Only 400–500 nm and 600–850 nm were used; 500–600 nm was dropped because Beer's law fails there. Resolution was 1 nm, sampled every 2 nm.
* The [T11B] ε_c and |μ_e|² values, omitting the two ε_c values at 520 nm.
* One "time-course" value: ε(500.2 nm, 35 °C) = 587.4 (8).
* Tamres & Bhat (1971), 410–495 nm at 60, 90 and 120 °C, with ±2 L mol⁻¹ cm⁻¹ errors.
* The B→C predissociation rates.
* The C-state RKR points.

**So the directly constrained temperature range is about 308–393 K** (35–120 °C). The 0 °C column of Table II is computed from the model.

**Reliability** (§IV D and §V):
* "For purely continuous absorption at λ ≤ 500 nm, the fitted molar absorptivity at 35 °C and its calculated counterpart at 0 °C should now be reliable within ∼0.5%, but full realization of such precision requires wavelength accuracy better than 0.2 nm."
* The statistical error in Fig. 10 "should be augmented by ∼0.5% in recognition of remaining intensity scale uncertainties."
* **Model error:** "Changes in the attachment point for the B potential most significantly affect the transition moment functions for B ← X and C ← X and the estimated C ← X spectrum. […] the differences in the latter exceed the statistical error; they are compensated by changes in the B ← X spectrum that preserve the greater precision of the total spectrum."
  * **The total is robust, but the B/C split is model dependent.** An implementation that swaps in a different B potential or μ_B should be validated on the total.
* The fit to the predissociation data is not statistically adequate: its probability is below 0.001 (§IV B).

### 2.4 [T11C] Table II: computed ε (L mol⁻¹ cm⁻¹) at 0 °C and 35 °C (p. 054301-9)

"Individual component bands provided as supplementary material (Ref. 25)". **That supplement is not available locally** (see §8).

**Caveat (§IV D).** "These results are from the pseudocontinuum model for the B state, so they are not easily related to experimental observations in the 500–630-nm region, where discrete vibrational structure remains prominent". So:
* For λ ≤ 495 nm the table is truly continuous: A + C + B bound–free.
* For 500–650 nm it includes smeared discrete B←X.
* Above about 650 nm it is mostly A←X, plus the smeared hot-band B←X.

| λ (nm) | 0 °C | 35 °C | λ (nm) | 0 °C | 35 °C | λ (nm) | 0 °C | 35 °C |
|---|---|---|---|---|---|---|---|---|
| 390 | 0.18 | 0.28 | 565 | 405.54 | 415.96 | 740 | 24.04 | 24.06 |
| 395 | 0.37 | 0.56 | 570 | 345.89 | 361.42 | 745 | 22.13 | 22.34 |
| 400 | 0.73 | 1.05 | 575 | 291.04 | 310.35 | 750 | 20.28 | 20.65 |
| 405 | 1.38 | 1.89 | 580 | 242.13 | 263.84 | 755 | 18.50 | 19.02 |
| 410 | 2.47 | 3.29 | 585 | 199.56 | 222.44 | 760 | 16.80 | 17.45 |
| 415 | 4.27 | 5.49 | 590 | 163.45 | 186.41 | 765 | 15.19 | 15.95 |
| 420 | 7.08 | 8.83 | 595 | 133.52 | 155.69 | 770 | 13.69 | 14.53 |
| 425 | 11.34 | 13.73 | 600 | 109.30 | 130.01 | 775 | 12.29 | 13.20 |
| 430 | 17.52 | 20.69 | 605 | 90.14 | 108.95 | 780 | 10.99 | 11.95 |
| 435 | 26.22 | 30.24 | 610 | 75.36 | 92.02 | 785 | 9.79 | 10.78 |
| 440 | 38.07 | 43.02 | 615 | 64.28 | 78.68 | 790 | 8.70 | 9.70 |
| 445 | 53.78 | 59.66 | 620 | 56.23 | 68.39 | 795 | 7.70 | 8.70 |
| 450 | 74.08 | 80.85 | 625 | 50.60 | 60.62 | 800 | 6.79 | 7.79 |
| 455 | 99.74 | 107.28 | 630 | 46.84 | 54.87 | 805 | 5.98 | 6.95 |
| 460 | 131.49 | 139.57 | 635 | 44.51 | 50.74 | 810 | 5.25 | 6.19 |
| 465 | 170.00 | 178.24 | 640 | 43.21 | 47.86 | 815 | 4.59 | 5.50 |
| 470 | 215.71 | 223.55 | 645 | 42.62 | 45.89 | 820 | 4.01 | 4.87 |
| 475 | 268.73 | 275.37 | 650 | 42.47 | 44.59 | 825 | 3.49 | 4.31 |
| 480 | 328.57 | 332.99 | 655 | 42.57 | 43.71 | 830 | 3.03 | 3.80 |
| 485 | 394.05 | 395.08 | 660 | 42.75 | 43.08 | 835 | 2.63 | 3.35 |
| 490 | 463.20 | 459.63 | 665 | 42.88 | 42.58 | 840 | 2.27 | 2.94 |
| 495 | 533.16 | 523.97 | 670 | 42.92 | 42.13 | 845 | 1.96 | 2.58 |
| 500 | 600.45 | 584.98 | 675 | 42.76 | 41.63 | 850 | 1.69 | 2.26 |
| 505 | 661.21 | 639.36 | 680 | 42.37 | 41.00 | 855 | 1.45 | 1.98 |
| 510 | 711.99 | 684.29 | 685 | 41.75 | 40.23 | 860 | 1.25 | 1.73 |
| 515 | 749.65 | 717.30 | 690 | 40.91 | 39.33 | 865 | 1.07 | 1.51 |
| 520 | 772.28 | 737.03 | 695 | 39.83 | 38.26 | 870 | 0.92 | 1.32 |
| 525 | 778.42 | 742.43 | 700 | 38.55 | 37.07 | 875 | 0.79 | 1.15 |
| 530 | 768.32 | 733.85 | 705 | 37.09 | 35.73 | 880 | 0.67 | 1.00 |
| 535 | 742.51 | 711.76 | 710 | 35.45 | 34.25 | 885 | 0.57 | 0.87 |
| 540 | 703.36 | 678.15 | 715 | 33.66 | 32.64 | 890 | 0.49 | 0.75 |
| 545 | 653.11 | 634.77 | 720 | 31.81 | 30.98 | 895 | 0.42 | 0.65 |
| 550 | 595.24 | 584.47 | 725 | 29.89 | 29.27 | 900 | 0.35 | 0.57 |
| 555 | 532.61 | 529.50 | 730 | 27.94 | 27.53 | | | |
| 560 | 468.56 | 472.66 | 735 | 25.98 | 25.79 | | | |

**Spot values quoted in the text** (§II and §IV D):

| Quantity | Value (L mol⁻¹ cm⁻¹) | Remarks |
|---|---|---|
| ε(500.2 nm, 35 °C), "time-course" measurement | **587.4 (8)** | Fitted input |
| ε(480 nm), measured at 35 °C and 64 °C | 333 and 337 | Compare 326(10) [T73], 340 (Lang & Strong 1965), 350 (Tamres & Bhat, claimed isosbestic point) |
| Temperature-independent point near room T | **487 nm, ε = 416** (text) | Where the 0 °C and 35 °C spectra intersect. Table II itself gives 486.1 nm and about 408 (§7.7) |
| Room-T estimate at 500 nm | **590 ± 4** | Compare Spietz 2006: 572 ± 6; Saiz-Lopez 2004: 599; Bauer 1998: 588 |
| Room-T calculated value at 436 nm | **31.0 ± 0.4** | Compare Saiz-Lopez: 40.0; Bauer: 36.9 ± 1.3 |
| Slope of ε near 480 nm | 3.5 %/nm | Wavelength errors matter |

The Saiz-Lopez et al. supplement spectrum "lies uniformly about 10 L mol⁻¹ cm⁻¹ above my 35° spectrum for the entire region 400–485 nm; at long wavelengths it drops to zero near 740 nm, largely missing the A ← X shoulder".

---

## 3. [T73]: the original resolution (historical; superseded for C)

### 3.1 Equations (pp. 2822–2823)

**Bound–bound**, eq. (2):

  ∫k_ν dν = (8π³ν/3hc)(N/q_r q_v)(g_A/g_B) s_J |R_e|² |⟨v′|v″⟩|² exp(−E/kT)

**Bound–continuum, "δ approximation"** (reflection approximation):
* Eq. (10): k_ν δν = (hν/4π) B_BA δN_ν
* **Eq. (11):**

  k_ν = (8π³ν/3hc)(N/q_v)(g_A/g_B) |R_e|² (dr/dν) Σ_v″ Ψ_v″²(r_ν) exp(−E_v″/kT)

  * (dr/dν) is "the reciprocal of the derivative of the upper-state potential curve at internuclear distance r_ν, fixed by the frequency ν".
  * Ψ_v″²(r) is the vibrational probability density of level v″.
  * q_v is the vibrational partition function, taken relative to v = 0.
* **Eq. (12):** k_ν/N = 2.303×10³ ε_ν/N₀
* Footnote 34: "the determination of |R_e|² for the 1u states involves a factor of 2 for the electronic degeneracy".

**Computational choices.**
* All continuum spectra used **J = 66**, "the average rotational number in the ground state". The FCFs were computed "for effective potential curves having J = 66" (footnote 37).
* The sums ran over v″ = 0–5, and up to 9 for B←X at 6700–6900 Å.
* X came from Le Roy's constants.

### 3.2 1u(¹Π) ← X, i.e. C←X (§V, p. 2827; Fig. 6)

* **Potential:** V(r) = **5.24×10⁷/r⁹** cm⁻¹, with r in Å, "relative to dissociation limit" (Table II footnote a).
  * The alternatives in Fig. 6 are 1.974×10⁷/r⁸ and 1.395×10⁸/r¹⁰. The r⁻⁹ form is preferred.
* **Strength:**
  * |R_e|² = **0.153 D²**, "assumed constant with wavelength".
  * f = 2.9×10⁻³.
  * Peak near **4985 Å, 200 L mol⁻¹ cm⁻¹**.
* **Errors:** "The maximum error in the peak position is estimated at ±50 Å, with the other quantities reliable within 10%."
* The crossing of the B curve's left branch at low v′ "seems reasonable".

### 3.3 A ← X (§VI, pp. 2828–2829; Fig. 8)

* **Potential for r ≤ 2.73 Å:** V = −4000 + 3.15×10⁸/r¹¹. This is an inverse-power repulsive curve "having an artificial dissociation limit 4000 cm⁻¹ below the actual limit".
* **Potential for r > 2.73 Å:** a Linnett function with m = 30 and r_e = 3.104 Å, using Brown's D_e = 660 cm⁻¹ and ω_e = 44 cm⁻¹.
  * Eq. (13): V(r) = [D_e/(m − x_e)] [x_e (r_e/r)^m − m exp(x_e − x_e r/r_e)]
  * Eq. (14): k_e r_e²/D_e = m x_e (m + 1 − x_e)/(m − x_e)
* **Strength:**
  * f = 6.2×10⁻⁴.
  * |R_e|² = **0.044 D²**. "Because of uncertainties in the spectrum for λ < 6400 Å, they may be in error by as much as 15%."
  * The peak is at **6730 Å (±50 Å), 41 ± 1 L mol⁻¹ cm⁻¹**.
* The Fig. 8 uncertainties include ±40% of the calculated B←X contribution.
* Brown's shallow A well (660 cm⁻¹) is now known to be wrong; the true D_e is about 1640 cm⁻¹, from [T11C] and Appadoo et al. 1996. **So the long-λ tail of this curve is obsolete.**

### 3.4 B ← X in [T73] (§VII, pp. 2830–2831)

* Bound–bound B←X was computed line by line (Doppler lines, J ≤ 200, v″ = 0–5), convolved with the trapezoidal slit.
* The B←X continuum (λ < 5200 Å) was computed with the δ approximation, eq. (11), using the repulsive branch in **Table IV** and the |R_e|² values in Table II.
* The repulsive branch above dissociation was *derived* by requiring |R_e|² to be linear in λ (and r̄). The turning points near the dissociation limit are estimated "accurate within 0.002 Å".

**Table IV** (p. 2830), the upper region of the B curve. Energies are relative to the B minimum, "at 15 770.45 cm⁻¹ above minimum of X curve".

| v | E (cm⁻¹) | r_min (Å) | r_max (Å) |
|---|---|---|---|
| Continuum | 7297 | 2.55 | |
| Continuum | 6815 | 2.56 | |
| Continuum | 6368 | 2.57 | |
| Continuum | 5953 | 2.58 | |
| Continuum | 5567 | 2.59 | |
| Continuum | 5208 | 2.60 | |
| Continuum | 4873 | 2.61 | |
| Continuum | 4559 | 2.62 | |
| 70 | 4351.3 | 2.626 | 7.369 |
| 68 | 4339.2 | 2.627 | 6.923 |
| 66 | 4323.8 | 2.628 | 6.545 |
| 64 | 4304.5 | 2.629 | 6.223 |
| 62 | 4280.8 | 2.630 | 5.948 |
| 60 | 4252.4 | 2.631 | 5.710 |
| 58 | 4218.9 | 2.632 | 5.501 |
| 56 | 4180.0 | 2.633 | 5.315 |

The recommended B-state constants are T_e = 15 770.45 cm⁻¹, D_e = 4381 cm⁻¹ and r_e = 3.0309 Å (Table V). The B limit, from Le Roy & Bernstein, is 20 044.0 ± 1.2 cm⁻¹ above X v″ = 0.

**On the continuum below the convergence limit** (p. 2831):
* "Calculations using Eq. (11) show that the absorption in the 5000-Å region is predominantly from v″ = 0 of the ground state, so continuum absorption from excited v″ levels is relatively unimportant within the banded region. However, absorption from excited rotational levels of v″ = 0 extends the continuum a considerable distance into the banded region, so that at 5061 Å the B–X extinction is estimated to be ∼40% continuum and ∼60% discrete."
* The cited detailed treatment is ref. 50: J. Tellinghuisen, "Continuous absorption below the band convergence limit in the halogen B–X transitions", J. Chem. Phys. 59, 849 (1973). **Not available locally.**

### 3.5 [T73] Table II: measured ε and the three-way resolution (p. 2826; room temperature, 22–27 °C)

* "Estimated maximum errors are ±3% or ±1 liter mole⁻¹ cm⁻¹, whichever is larger."
* ε_meas is at about 26 Å resolution. For 5100–6200 Å the values are zero-pressure extrapolations **and depend on the slit**, so they are not a continuum.
* One text-layer error was corrected from the page image: ε₁ᵤ←X = 1 at 6400 Å.

| λ (Å) | ε_meas | ε₁ᵤ←X ᵃ | ε_A←X ᵇ | ε_B←X ᶜ | \|R_e\|²_B←X (D²) ᵈ |
|---|---|---|---|---|---|
| 4200 | 12.5 | 4 | | 8 | |
| 4300 | 22.2 | 11 | | 11 | |
| 4400 | 43.7 | 27 | | 17 | |
| 4500 | 79 | 52 | | 27 | 0.87 |
| 4600 | 138 | 89 | | 50 | 0.81 |
| 4700 | 221 | 130 | | 91 | 0.81 |
| 4800 | 326 | 169 | | 157 | 0.81 |
| 4900 | 454 | 193 | | 261 | 0.85 |
| 5000 | 574 | 199 | | 375 | 0.85 |
| 5100 | 710 | 187 | | 523 | 0.91 |
| 5200 | 753 | 162 | | 591 | 0.88 |
| 5300 | 810 | 131 | | 679 | 0.97 |
| 5400 | 720 | 100 | 1.5 | 619 | 0.97 |
| 5500 | 600 | 72 | 2.5 | 525 | 0.95 |
| 5600 | 505 | 50 | 4.5 | 450 | 1.02 |
| 5700 | 341 | 34 | 6.8 | 301 | 0.99 |
| 5800 | 224 | 22 | 9.7 | 192 | 1.00 |
| 5900 | 174 | 14 | 13.2 | 147 | 1.03 |
| 6000 | 155 | 9 | 17.4 | 129 | 1.33 |
| 6100 | 71 | 5 | 21.9 | 44 | 1.08 |
| 6200 | 78 | 3 | 26.5 | 48 | 1.33 |
| 6300 | 47.0 | 2 | 31.0 | 14.0 | 0.90 |
| 6400 | 45.8 | 1 | 34.9 | 10.1 | |
| 6500 | 46.0 | | 38.0 | 7.0 | |
| 6600 | 45.2 | | 40.0 | 4.2 | |
| 6700 | 42.0 | | 40.9 | 2.4 | |
| 6800 | 41.7 | | 40.6 | 1.2 | |
| 6900 | 40.0 | | 39.4 | 0.5 | |
| 7000 | 38.1 | | 37.3 | | |
| 7100 | 34.5 | | 34.5 | | |
| 7200 | 31.0 | | 31.0 | | |
| 7400 | 24.0 | | 24.0 | | |
| 7600 | 16.9 | | 16.9 | | |
| 7800 | 12.2 | | 12.2 | | |
| 8000 | 8.4 | | 8.4 | | |

Table II footnotes:
* ᵃ Calculated for J = 66, V(r) = 5.24×10⁷/r⁹.
* ᵇ Values for λ < 7000 Å calculated with V(r) = −4000 + 3.15×10⁸/r¹¹.
* ᶜ Values for λ > 6300 Å come from a linear extrapolation of |R_e|² from λ < 6000 Å.
* ᵈ Values for λ < 5200 Å come from the δ approximation. For λ ≥ 5200 Å the slit function was trapezoidal, with base from (λ−16) to (λ+13) Å and peak from (λ−13) to (λ+10) Å.

### 3.6 [T73] Table III: argon-line "continuum" absorptivities (p. 2827; 75-cm cell)

These are upper bounds on the continuum between B←X lines.

| λ (Å) | Appearance | ε_continuum |
|---|---|---|
| 5060.08 / 5062.07 | true continuum | 400 ± 20 ᵃ |
| 5187.75 | slight rotational-line absorption | 166 (+4, −8) |
| 5451.65 | some rotational-line absorption | 95 (+4, −8) |
| 5606.73 | true continuum | 52.5 ± 2 |
| 5912.08 | true continuum | 34 ± 3 |

* ᵃ "includes a contribution of roughly 200 from the 1u(¹Π)←X transition, with most of the rest being due to B←X continuum absorption from excited vibrational and rotational levels of the X state."
* [T11B] (p. 084301-7) re-evaluated the 560.673-nm line as **46 (3)** L mol⁻¹ cm⁻¹, about 4σ below the new ε_c.

---

## 4. The 1982 reassessment [T82], as quoted in 2011 (paper not available)

All that the 2011 papers say about [T82]:

* **Method** ([T11C] §I, p. 054301-1):
  * The C←X bound–free spectrum "was computed in trial-and-error fashion to match these estimates". Those estimates came from:
    * the continuum between B←X lines at λ > 520 nm (refs. 10–12);
    * I*/I photodissociation branching ratios at λ < 500 nm (Oldman et al. 1971; Wiesenfeld & Young 1981);
    * the I* → I photodissociation laser (Davis 1978).
  * "[…] and the remaining absorption was allocated to B ← X and used to estimate the R-dependence of this system's electronic transition strength |μ_e(R)|² in its continuum."
* **C potential** ([T11C] p. 054301-6): "Updated analysis [T82] of the C ← X absorption replaced the R⁻⁹ curve with one of form **R⁻⁹·⁵**, which lies very close to the CB 6–12 curve in the relevant R region; in fact both cross the B curve only 0.0001 Å apart, near **2.896 Å and 3535 cm⁻¹** above the first dissociation limit." **The 1982 coefficient is not quoted.**
* **C←X strength, 2011 compared with 1982:**
  * [T11C] abstract: "the C ← X spectrum is most altered from the previous analysis, being now ∼20% weaker".
  * [T11C] §V: "its **18% smaller integrated intensity**".
  * [T11B] abstract: "lower the C(1u) ← X transition strength by **25%**".
    * Reading: the 25% was inferred in [T11B] from ε_c at 520–590 nm alone. The global fit of [T11C] supersedes it, with 18% in integrated intensity. The paper does not reconcile the two numbers.
  * [T11B] p. 084301-7: the new ε_c "lie systematically low for shorter wavelengths [below 600 nm] […] This would mean a further reduction from its original estimated strength, in Ref. 17 [T73]."
    * So the order is T73 > T82 > T11 for C←X.
  * [T11B] §V: "The changes in the continuum absorption in the 520–590-nm region **cannot be accommodated by simply scaling down the previous computed C ← X spectrum**, meaning a refinement of the C potential curve will be needed."
    * **So do not model C←X as 0.75 × (the 1982 band).** Use the [T11C] potential.
* **A←X:**
  * [T11B]: "the A ← X system changed little in the Ref. 16 reassessment".
  * [T11C] p. 054301-9: "the A ← X spectrum and the small-R extension for the A potential curve have changed little since my 1982 report. The analysis was actually presented earlier, in conjunction with a first estimated RKR potential for the A state [Viswanathan, Sur & Tellinghuisen, J. Mol. Spectrosc. 86, 393 (1981)]. […] the 1981 curve lies 0.006 Å to smaller R, placing it ∼80 cm⁻¹ lower at that R. However, the current exponential extension and the previous **R⁻¹¹ curve** agree within the current statistical error for R = 2.43–2.73 Å, and differ by at most 30 cm⁻¹ in the range 2.53–2.76 Å."
* **B←X |μ_e|²** ([T11B]):
  * The new values "indicate a weaker rise in |μ_e(R)|² with increasing R — by 20% as compared with my 1982 assessment".
  * "the rate of increase with λ and R actually agrees better with my original analysis [T73] than my reassessment [T82]."
* **Ar-line ε_c** ([T11B] p. 084301-7):
  * The [T82] Table II reanalysis of the [T73] Ar-line data "for the 300 K measurements" was "largely confirmed".
  * The exception is 560.673 nm, now 46(3), which is insensitive to the assumed Ar kinetic temperature (800 K compared with the 500 K used before).

**Verdict: [T82] is not needed** to build the model, because [T11C] Table I is a complete and more precise parameter set. It would only allow a historical comparison of components. The [T11C] **supplementary material** (the individual component bands) is far more useful; see §8.

---

## 5. How the papers treat B←X bound–free

**Neither paper gives an analytic parameterization. Both compute B←X from a B-state potential and a transition moment**, with the inner wall above dissociation fitted to absorption data.

* **[T73]:**
  * Uses the δ (reflection) approximation, eq. (11), with the Table IV repulsive branch and |R_e|² ≈ 0.81–0.87 D² at 4500–5000 Å (Table II).
  * The inner wall was *derived* by requiring |R_e|² to vary linearly with λ.
  * It uses v″ = 0–5 and J = 66.
  * It points out that continuum absorption from excited *rotational* levels of v″ = 0 extends well into the banded region: at 5061 Å it is about 40% of the B←X absorption. The detailed treatment is in J. Chem. Phys. 59, 849 (1973), not available locally.
* **[T11B]** (p. 084301-7): "there is truly continuous B ← X absorption to the red of the spectroscopic dissociation limit, due to absorption from excited υ″ and J″ levels [JCP 59, 849]. This often-neglected behavior probably accounts for at least half of the apparent excess in ε_c at 520 nm relative to its value at 525 nm."
  * In the [T11B] simulation, the discrete lines were limited to υ″ = 0–9, υ′ = 0–55 and J″ = 0–210.
  * Anything else, including the true B←X continuum and any missing high-υ′ lines, was absorbed into the fitted ε_c. That is why the two 520-nm ε_c values are "questionable".
* **[T11C]:**
  * An exact quantum bound–free calculation for the B state.
  * The small-R wall is the one-term exponential in Table I, attached at υ = 33, with μ_B(R) quadratic in R.
  * The **pseudocontinuum** trick sets U_B = T_e(B) for R > R_e. It turns the *whole* B←X system, discrete bands included, into one smooth continuum with correct Franck–Condon envelopes at low resolution.
  * This is exact for λ < 498.9 nm (vacuum, υ″ = 0), which is truly bound–free. Above that it gives the local average of the discrete structure, which is what the 1-nm spectra measure.
  * The model error is B↔C compensation (§2.3).

---

## 6. Validation targets (summary)

| # | Target | T | What it tests | Source |
|---|---|---|---|---|
| V1 | ε at 390–495 nm, Table II (§2.4) | 0 °C, 35 °C | **A + C + B bound–free**, all truly continuous. Stated reliability about 0.5% for λ ≤ 500 nm, if the wavelength is good to 0.2 nm | [T11C] |
| V2 | Point where ε(0 °C) = ε(35 °C): Table II gives 486.1 nm, ε ≈ 408. The text says 487 nm and 416, which is inconsistent with the table (§7.7) | 0 → 35 °C | Temperature dependence | [T11C] |
| V3 | ε(500.2 nm) = 587.4 (8); 590 ± 4 at room T at 500 nm | 35 °C; room T | B-continuum edge plus C | [T11C] |
| V4 | ε(436 nm) = 31.0 ± 0.4 | room T | Blue wing: C + B | [T11C] |
| V5 | ε at 500–900 nm, Table II | 0 °C, 35 °C | Only if i2spec also runs the **pseudocontinuum** B←X, or averages its line-by-line B←X plus continuum over about 1 nm. At 650–900 nm the table is mostly A←X plus a smeared hot-band B←X | [T11C] |
| V6 | ε_c at 520–635 nm (Table II in `intensity-inputs.md`, or §7.3 below) | **cell 34–38 °C** (the side arm was at 21–27 °C) | **A + C (+ B bound–free from hot levels near 520 nm)**. These were also fitted by [T11C] | [T11B] |
| V7 | Ar-line values: 560.673 nm → 46 (3); 591.208 nm → 34 ± 3 | 22–27 °C | Old; loose check on A + C | [T73], [T11B] |
| V8 | A←X peak 6730 ± 50 Å, 41 ± 1; C←X peak 4985 ± 50 Å, 200 (±10%) | 22–27 °C | Historical. [T11C] lowered C by about 18–25%, so its 1973 peak should come out lower | [T73] |

`intensity-inputs.md` calls the [T11B] ε_c values "room temperature". [T11B] §III says "The cell body was typically at 34–38 °C", and "at the cell temperature 36 °C employed in the experiments, the most populous J″ level is 53". The t₂ column in that table is the *side-arm* temperature that sets the pressure. So compare at about 309 K.

Other landmarks, from Hannover:
* The B dissociation limit seen from X(0, 0) is at 20 043.22 cm⁻¹ = **498.92 nm (vacuum)**. Absorption from υ″ = 0, J″ = 0 is purely bound–free below that.
* The X limit is at 12 440.24 cm⁻¹ = 803.84 nm. C←X and A←X from υ″ = 0 can reach the continuum only at shorter λ.

---

## 7. Numerical check: reproducing [T11C] Table II from Table I

A scratch script, not in the repo, implements §2 as written. It uses no tuned parameters.

**Setup.**

| Item | What the script uses |
|---|---|
| Formula | Eq. (2) with **G_ab = 1 and the tabulated μ** |
| Continuum functions | Numerov, h = 0.001 Å, R = 2.15–10 Å. Energy normalization per cm⁻¹ from the WKB amplitude over the last 0.5 Å: a² = u²k + u′²/k is scaled to 1/(πħ²/2μ) |
| X state | i2spec Hannover X (sinc-DVR, pinned with `solver="dvr"`; `model.py` now defaults to a B-spline solver), with υ″ ≤ 29 |
| Rotational sum | J″ = 5, 15, …, 345, each with weight 10. Full Boltzmann sum, no nuclear-spin weights. Q branch with the J(J+1) centrifugal term |
| Energy references | T_e(A) = D_e(X) − 1639.9 = 10 907.44 cm⁻¹. T_e(B) = 15 769.068 cm⁻¹ (Hannover). C relative to D_e(X) = 12 547.34 cm⁻¹ |
| B wall | The [T11C] exponential for R < 2.67134 Å (the Hannover υ′ = 33 inner turning point), Hannover beyond. Pseudocontinuum: U = T_e for R > R_e |
| A beyond the join | For R > 2.805916 Å, a Morse inner branch with D = 1639.9, matched in value and slope (β = 2.0163 Å⁻¹, R_e = 3.1331 Å). **This stands in for the A RKR curve, which I don't have** |
| C tail | For R > 3.568 Å, where U_C = 200 cm⁻¹: 200·(3.568/R)^27.8 |

It runs in about 11 s for all states and all J.

### 7.1 Transcription checks

* **A:** eq. (4) with Table I gives U(2.804259) = 1451.283 and U(2.805916) = 1431.381 cm⁻¹. The paper's eq. (5) gives 1451.281 and 1431.379. ✓
* **C:** U_C(2.666 Å) = 7811.6, U_C(2.83) = 4388.3 and U_C(2.896) = 3496.5 cm⁻¹ above the asymptote. Compare the old R⁻⁹·⁵ and CB crossing of B at 2.896 Å and 3535 cm⁻¹.
* **B:** T_e(B) − D_e(X) = **3221.73** cm⁻¹ from Hannover. **The "3321.8 cm⁻¹" printed in [T11C] p. 054301-4 is therefore almost certainly a typo for 3221.8.** It does not matter for the model (footnote 31).

### 7.2 Result: computed total (A + C + pseudocontinuum B) ÷ [T11C] Table II

| λ range | 0 °C | 35 °C | Comment |
|---|---|---|---|
| 390 nm | 1.021 | 1.017 | ε is only 0.18–0.28 |
| **400–500 nm** | **0.994–1.003** | **0.993–1.002** | Truly continuous region, V1 |
| 505–640 nm | 0.987–1.022 | 0.987–1.019 | Pseudocontinuum B depends on the B wall at 2.63–2.8 Å (Hannover here, G&L RKR in [T11C]) |
| 650–800 nm | 1.002–1.011 | 1.003–1.009 | A←X dominates |
| 850–900 nm | 1.017–1.030 | 1.009–1.015 | Morse stand-in for the A wall beyond 2.806 Å |

This confirms the G_ab/√2 convention, the energy references, the 108.861 constant with per-cm⁻¹ normalization, and that the Hannover X potential is adequate. Using G = 2 with the tabulated μ would double A and C and would fail at once.

### 7.3 Component bands from the check (L mol⁻¹ cm⁻¹)

These stand in for the unavailable supplement.

| λ (nm) | A, 0 °C | C, 0 °C | B (pseudo), 0 °C | A, 35 °C | C, 35 °C | B (pseudo), 35 °C | B bound–free only (true B potential), 35 °C |
|---|---|---|---|---|---|---|---|
| 420 | 0.00 | 6.58 | 0.46 | 0.00 | 8.06 | 0.71 | 0.71 |
| 440 | 0.00 | 31.50 | 6.39 | 0.00 | 34.32 | 8.43 | 8.41 |
| 460 | 0.00 | 83.89 | 46.91 | 0.00 | 84.20 | 54.51 | 54.36 |
| 480 | 0.00 | 134.54 | 192.93 | 0.00 | 129.04 | 202.75 | 202.19 |
| 490 | 0.01 | 144.66 | 318.52 | 0.01 | 137.52 | 321.88 | 321.00 |
| 500 | 0.02 | 141.47 | 459.58 | 0.04 | 134.49 | 451.07 | 405.29 |
| 520 | 0.17 | 106.12 | 655.75 | 0.24 | 103.41 | 624.00 | 18.72 |
| 560 | 2.73 | 29.20 | 444.32 | 3.25 | 32.39 | 444.15 | 0.02 |
| 600 | 14.50 | 4.53 | 92.65 | 15.19 | 6.14 | 111.15 | 0.00 |
| 640 | 33.76 | 0.54 | 9.15 | 32.81 | 0.92 | 14.46 | 0.00 |
| 680 | 41.82 | 0.06 | 0.58 | 39.71 | 0.13 | 1.28 | 0.00 |
| 750 | 20.42 | 0.00 | 0.00 | 20.79 | 0.00 | 0.01 | 0.00 |
| 800 | 6.87 | 0.00 | 0.00 | 7.86 | 0.00 | 0.00 | 0.00 |

**Band properties.**
* **C←X:**
  * The peak is at about 490–495 nm, with ε = 144.7 at 0 °C and 139.1 at 300 K.
  * f(390–900 nm) = 2.16×10⁻³. That is 25% below [T73] (2.9×10⁻³, with a 200 peak at 498.5 nm).
  * If the "18% smaller" in [T11C] is relative to [T82], the implied 1982 f is about 2.6×10⁻³. **This is inference made here.**
* **A←X:** the peak is at about 675 nm, with ε = 41.9 at 0 °C and 40.3 at 300 K. f = 5.67×10⁻⁴. [T73] gave 6730 Å, 41 ± 1, and f = 6.2×10⁻⁴.

**Above the B dissociation limit, true and pseudocontinuum B agree.** The pseudocontinuum B and the true-potential B bound–free differ by ≤ 0.3% at 480–495 nm. That is the reflection principle behind [T11C]. At 500 nm and beyond, the pseudocontinuum also contains the smeared discrete bands.

### 7.4 [T11B] ε_c compared with the computed A + C at 308 K (V6)

| λ (nm) | A + C | + B bound–free | [T11B] Table II |
|---|---|---|---|
| 520 | 103.65 | 122.36 | 120.3 / 132.5 (flagged questionable) |
| 525 | 93.51 | 102.01 | 91.4 / 90.5 |
| 530 | 83.31 | 86.56 | 84.7 / 85.3 |
| 534 | 75.35 | 77.26 | 76.9 / 78.0 |
| 540 | 64.09 | 64.75 | 63.8 / 58.6 |
| 545 | 55.54 | 55.84 | 57.4 / 53.8 |
| 555 | 41.27 | 41.33 | 38.1 / 41.2 |
| 570 | 27.32 | 27.33 | 28.1 / 25.8 |
| 580 | 22.58 | 22.58 | 21.2 / 21.9 |
| 590 | 20.80 | 20.80 | 19.5 / 23.1 |
| 600 | 21.33 | 21.33 | 21.1 / 21.9 |
| 610 | 23.51 | 23.51 | 26.3 / 27.0 / 23.0 / 23.9 |
| 625 | 28.46 | 28.46 | 31.5 / 31.1 |
| 635 | 32.04 | 32.04 | 35.1 / 32.8 |

* [T11C] fitted A + C to these points; its pseudocontinuum B cannot hold them. The agreement is within about 1–2σ, where 1σ ≈ 2–4 including the background term.
* At 520–525 nm the real B←X bound–free from hot (υ″, J″), 8–19 L mol⁻¹ cm⁻¹, is not negligible. So **validate i2spec's A + C + B bound–free against [T11B] only at λ ≥ 530 nm**, and against A + C at 525 nm.

### 7.5 The B inner wall: Hannover compared with [T11C]

Energies are E − T_e(B) in cm⁻¹.

| R (Å) | Hannover | [T11C] eq. (4) | Difference |
|---|---|---|---|
| 2.45 | 9928.5 | 12243.6 | +2315 |
| 2.50 | 8222.9 | 9344.9 | +1122 |
| 2.55 | 6630.4 | 7045.2 | +415 |
| 2.60 | 5143.7 | 5220.7 | +77 |
| 2.62 | 4577.0 | 4601.3 | +24 |
| 2.64 | 4025.6 | 4036.6 | +11 |
| 2.66 | 3506.7 | 3521.9 | +15 |
| 2.68 | 3040.5 | 3052.7 | +12 |
| 2.70 | 2618.2 | 2625.0 | +7 |

* Below its R_I = 2.647 Å, the Hannover curve continues as A_I exp[−B_I(R − R_I)] with B_I = 1.37 Å⁻¹. That is far too soft.
* With the raw Hannover wall, the B bound–free, relative to the result with the [T11C] wall, falls to 0.000× at 390–410 nm, 0.19× at 440 nm and 0.44× at 450 nm, and rises to 1.03–1.04× at 480–485 nm. **The Hannover inner wall must not be used for the continuum.**
* At 2.64–2.70 Å, inside the discrete data range, [T11C] lies 7–15 cm⁻¹ above Hannover, i.e. about 0.0006 Å to larger R. [T11C] notes the same offset for the Salumbides 2008 quantal potential.
  * At the Hannover υ′ = 33 inner turning point (2.67134 Å; E − T_e = 3236.92 cm⁻¹), [T11C] is **+13.6 cm⁻¹** higher. The slopes are −23 306 ([T11C]) and −23 166 (Hannover) cm⁻¹/Å.
  * The two curves cross at **R = 2.7154 Å** (E − T_e = 2322 cm⁻¹).
* [T73] Table IV, rebased to T_e = 15 769.07, gives E = 7298 at 2.55 Å and 4874 at 2.61 Å. [T11C] gives 7045 and 4904.

### 7.6 Temperature dependence (extrapolated beyond the validated range)

**Fraction of X population outside the [T11C] υ″ range:**

| T | υ″ ≥ 10 | υ″ ≥ 20 |
|---|---|---|
| 273 K | 1.9×10⁻⁵ | |
| 308 K | 6.6×10⁻⁵ | |
| 500 K | 2.8×10⁻³ | |
| 700 K | 1.5×10⁻² | 3.0×10⁻⁴ |

So at 700 K, sum υ″ ≤ 25 or more.

**B←X bound–free from the true B potential (T11C wall), ε in L mol⁻¹ cm⁻¹.** The υ″ = 0, J = 0 limit is at 498.9 nm (vacuum):

| λ (nm) | 273 K | 308 K | 500 K | 700 K |
|---|---|---|---|---|
| 495 | 388.8 | 386.7 | 359.9 | 327.8 |
| 500 | 405.0 | 405.3 | 379.1 | 343.9 |
| 505 | 219.0 | 239.5 | 287.6 | 289.6 |
| 510 | 78.0 | 97.3 | 173.3 | 208.1 |
| 520 | 11.9 | 18.7 | 67.1 | 110.4 |
| 530 | 1.6 | 3.3 | 24.3 | 55.7 |
| 540 | 0.26 | 0.66 | 9.4 | 29.0 |
| 560 | 0.01 | 0.02 | 1.3 | 7.1 |
| 600 | 0 | 0 | 0.04 | 0.57 |

**Total A + C + pseudocontinuum B:**

| λ (nm) | 273 K | 308 K | 500 K | 700 K |
|---|---|---|---|---|
| 400 | 0.73 | 1.05 | 4.07 | 9.00 |
| 450 | 73.8 | 80.4 | 111.6 | 133.4 |
| 500 | 601.1 | 585.6 | 508.7 | 448.2 |
| 600 | 111.7 | 132.5 | 217.2 | 261.5 |
| 700 | 38.7 | 37.2 | 36.8 | 49.0 |
| 800 | 6.87 | 7.86 | 11.5 | 14.2 |
| 900 | 0.36 | 0.58 | 2.26 | 4.14 |

* C←X peak ε: 144.7 at 273 K, 110.0 at 500 K, 93.0 at 700 K. Its f is essentially constant.
* A←X peak ε: 41.9 at 273 K, 31.8 at 500 K, 26.8 at 700 K.

### 7.7 Spot values, wavelength convention and the B join

[T11C] does not say whether its λ values are in air or vacuum. I evaluated both, using n − 1 = 2.78×10⁻⁴ for air.

| Target | [T11C] | Computed, λ as vacuum | Computed, λ as air |
|---|---|---|---|
| ε(500.2 nm, 35 °C), the time-course input | 587.4 (8) | **587.8** | 589.4 |
| ε(500 nm, "room T") | 590 ± 4 | 590.9 / 590.0 / 589.2 at 296 / 298 / 300 K | 592.5 / 591.6 / 590.8 |
| ε(436 nm, "room T") | 31.0 ± 0.4 | 30.88 / 31.12 / 31.34 at 296 / 298 / 300 K. At 298 K: C 26.16, B 4.96 | 31.16 / 31.40 / 31.62 |
| ε(480 nm, 35 °C) | 333 (measured) | 331.8 | 333.4 |
| Point where ε(0 °C) = ε(35 °C) | "487 nm, ε = 416" (text) | **486.0 nm, 406.7** | 486.0 nm, 408.6 |
| Table II ratio at 400–500 nm | | 0.993–1.003 | 1.000–1.018 (1.001–1.006 at 430–500) |

**Findings.**
* **Air or vacuum:** the two readings cannot be told apart at the 0.5% level. Vacuum reproduces the heavily weighted 500.2-nm time-course input better: +0.4 against +2.0, with σ = 0.8. ~~Treat Table II λ as vacuum.~~ **Superseded by §10:** this check's normalization was biased high by 0.3–0.7%. Without the bias, air fits better, so treat Table II λ as air.
* **The temperature-independent point:** the text's "487 nm, ε = 416" does not match [T11C]'s own Table II. Linear interpolation of Table II puts the crossing at **486.1 nm, with ε ≈ 408–409**. The computed value (486.0 nm, 406.7) agrees with the table. **Use the table, not the text value, as target V2.**
* **The B join:** joining the [T11C] wall at its crossing with Hannover (2.7154 Å, continuous) rather than at υ′ = 33 (2.6713 Å, with a 13.6 cm⁻¹ step) changes ε_B,bf. At 308 K:

  | λ (nm) | Change |
  |---|---|
  | 450 | −0.1% |
  | 480 | −0.15% |
  | 495 | −0.6% |
  | 500 | −0.3% |
  | 505 | +0.4% |
  | 520 | +0.9% |

  This is a small model sensitivity, well inside the [T11C] model error. §7.2 used the υ′ = 33 join. **Prefer the continuous join** (§9.1).

---

## 8. Missing pieces and papers to obtain

| Priority | Item | Why | Needed for |
|---|---|---|---|
| 1 | **[T11C] supplementary material** (doi:10.1063/1.3616039, ref. 25): "the spectral data and other information", including the **individual component bands** at 0 and 35 °C and a residual plot | Lets A, B and C be validated separately, instead of relying on §7.3 | Validation |
| 2 | **The A-state potential from R ≈ 2.806 Å to R_e(A).** [T11C] used the RKR curve of J. Tellinghuisen, J. Chem. Phys. 118, 3532 (2003), up to υ = 29 on the repulsive branch. An alternative is Appadoo et al., J. Chem. Phys. 104, 903 (1996) | Needed for A←X at λ ≳ 810–850 nm. The Morse stand-in in §7 is 1–3% off there | Model |
| 3 | J. Tellinghuisen, J. Chem. Phys. 59, 849 (1973), "Continuous absorption below the band convergence limit in the halogen B–X transitions" | Theory of the B←X continuum from excited υ″ and J″; background for §9.4 | Background |
| 4 | **Absorption data above 393 K.** Sulzer & Wieland, Helv. Phys. Acta 25, 653 (1952): pure I₂ at 600–1050 °C, cited by [T73]. Tamres & Bhat, J. Phys. Chem. 75, 1057 (1971): 60–120 °C | These papers contain no validation above 393 K | Validation at 400–700 K |
| 5 | [T82], J. Chem. Phys. 76, 4736 (1982) | Historical only; superseded by [T11C] (§4) | Nothing |
| 6 | Gerstenkorn & Luc (1985) B-state RKR constants | Only for exact reproduction of the [T11C] pseudocontinuum B at 505–640 nm. With Hannover instead, the difference is 1–2% | Optional |

**Not needed:** the C-state R⁻ⁿ tail and the van der Waals RKR points (Inard et al. 1999).
* C←X from υ″ = 0 cannot reach the continuum beyond 803.8 nm.
* At λ ≤ 900 nm the absorption samples U_C only above about 150 cm⁻¹, i.e. at R ≲ 3.6 Å.

---

## 9. Recommendation for i2spec

### 9.1 Model

**Adopt the [T11C] model**, for all three transitions and the whole 250–700 K range:
* the quantum bound–free eq. (2), with the Table I potentials and moments;
* evaluated with i2spec's own X states;
* summed over the Boltzmann distribution at the requested T;
* σ_c = 3.8235×10⁻²¹ (ε_A + ε_C + ε_B,bf) cm².

**Why not an empirical ε(λ, T):**
* None is published in these papers, and [T11C] Table II covers only 0 and 35 °C.
* The B←X bound–free part is a truncated edge whose hot-band tail grows rapidly with T. At 520 nm it goes from 12 to 110 L mol⁻¹ cm⁻¹ between 273 and 700 K (§7.6). A smooth band shape in ν cannot represent that.
* A physical model extrapolates in T through the X populations alone.

**The [T73] δ approximation and its r⁻⁹ C curve are superseded.** C←X is now 25% weaker in f and has a different shape ([T11B]: "cannot be accommodated by simply scaling").

| Transition | Upper potential. Energy zero = X minimum; z = R − R₀ | μ_e(R) in D, with z = R − 2.7 Å | G_ab |
|---|---|---|---|
| **A←X** | For R < 2.805916 Å: U = T_e(A) + [−226.70 + 3401.305 exp(−6.36z − 4.0z²)], R₀ = 2.7, with T_e(A) = D_e(X) − 1639.9 = 10 907.44 cm⁻¹. For R ≥ 2.805916 Å: the A RKR repulsive branch (§8, item 2; interim: the Morse stand-in of §7) up to R_e(A). Beyond that, U = T_e(A) (pseudocontinuum) | 0.2845 − 0.048z | 1 |
| **C←X** | U = D_e(X) + [−211.0 + 4599.3 exp(−3.297z + 0.505z² − 0.46z³ − 0.33z⁵)], R₀ = 2.83. Beyond U_C ≈ 200 cm⁻¹ (R ≈ 3.57 Å), any smooth repulsive tail to 0 | 0.4714 − 0.16z | 1 |
| **B←X bound–free** | "Continuum B potential": for R < 2.7154 Å (its crossing with Hannover), U = T_e(B) + [−1783.226 + 4408.223 exp(−4.63z)], R₀ = 2.7. Beyond that, Hannover B with its true asymptote D_e(B) = 20 150.317 cm⁻¹ | [T11B] eq. (10), already in i2spec, so the lines and the continuum share one μ_B | 1 |

**Notes on the table.**
* **B join point:** the choice matters at the 0.1–0.9% level (§7.7). Prefer the continuous join at 2.7154 Å to the υ′ = 33 join, which has a 13.6 cm⁻¹ step.
* **μ_B:** [T11B] eq. (10) and the [T11C] Table I quadratic agree to within 1.5% over 2.45–2.7 Å. With eq. (10), the B bound–free is 0.985–1.000× at 390–500 nm. That is ≤ 0.2% on the total everywhere (§7.2).
* **Wavelengths:** treat the Table II λ values as air (§10; §7.7 said vacuum).

### 9.2 B←X: parameterize or compute?

**Compute it from the potentials**, for four reasons:
1. **No published parameterization of the B←X bound–free part alone exists.**
   * [T73] gives it only at room temperature, at 420–510 nm, and with the old C.
   * [T11C] folds it into a pseudocontinuum total. The component bands are in the unavailable supplement, and only at 0 and 35 °C.
2. **Its T-dependence is the strongest of the continuum pieces** (§7.6). It comes entirely from the X populations, which i2spec already has.
3. **It must meet the discrete line list exactly at the dissociation limit** (§9.4). Only a calculation with the same X states and μ_B can guarantee that.
4. **It works.** With the [T11C] inner wall, the calculation reproduces [T11C] to ≤ 0.7% at 400–500 nm (§7.2).

**But not with the Hannover B potential alone.** Below 2.647 Å its inner extension is about 1100 cm⁻¹ too soft at 2.50 Å. That would suppress B←X at 400–450 nm by factors of 2 to ∞ (§7.5). Hence the "continuum B potential" in §9.1.

**Energy-normalized or discretized continuum?**
* **Energy-normalized continuum (recommended).** Numerov outward integration plus the WKB amplitude, as in §7.
  * It is exact on the grid, and every energy is computed independently on any E grid.
  * All states and 35 J values take about 11 s.
* **Discretized continuum: a box of sinc-DVR or B-spline states.** Box states above threshold give the continuum density |⟨n|μ|υ″⟩|²/ΔE_n, where ΔE_n is the local level spacing. It reuses the existing solvers and matrix-element code, but:
  * The box must be long: L ≈ 15–20 Å for ΔE ≈ 10–15 cm⁻¹ at 5000 cm⁻¹ above threshold.
  * It must also be fine: step ≲ 0.005 Å to resolve k up to about 200 Å⁻¹ at 400 nm. That means dense eigenproblems of N ≈ 3000–4000 for every J and state.
  * The result needs smoothing over several box states.
  * Near threshold, box states mix with quasi-bound (centrifugal) resonances.
  * **Use it as a cross-check only.**

### 9.3 Numerics

| Item | Choice |
|---|---|
| R grid | From about 2.15 Å, deep in the classically forbidden region for every energy used, to ≥ 10 Å. h = 0.001 Å |
| Normalization | WKB amplitude over the last 0.5 Å, where a² = u²k + u′²/k is set to 1/(πħ²/2μ). Drop energies that are classically forbidden there |
| E′ grid | 10–20 cm⁻¹ steps, from each asymptote (or the line-list cutoff, §9.4) to ν_max + E″_max |
| Caching | Tabulate M²(E′; υ″, J″), which is independent of T. Cache it per isotopologue, like `MasterLineList`. ε(ν, T) is then a Boltzmann sum |
| Boltzmann sum at 700 K | υ″ ≤ 25–30. J″ on a 10-step grid up to about 350. The check used J″ = 5…345, a full sum; [T11C] says 5 J values, or even the average J, suffice |
| Rotation | Effective potentials with the J(J+1)/R² term, and J′ = J″ (Q branch). No nuclear-spin factor ([T11C]: it averages out) |
| Output grid | σ_c is smooth. Evaluate it on a grid of about 5 cm⁻¹ and interpolate onto the spectrum grid. No line shape |
| Degeneracy | G_ab = 1 with the Table I μ (or G = 2 with μ/√2), **never both** |
| Isotopologues | Same potentials and moments, with each isotopologue's own reduced mass |

### 9.4 The boundary between B←X lines and continuum

* **The gap.** As of the code read on 2026-09-15, `intensity_model` uses nlev_b = 60 and `SHARED_GRID` rmax = 7.0 Å, plus a box-state filter. So the line list stops near υ′ = 59.
  * The bound levels from υ′ ≈ 60 up to the last bound level lie within about 130 cm⁻¹ of D_e(B) ([T73] Table IV: υ′ = 60 at 4252.4, D_e = 4381 cm⁻¹).
  * Their bands sit at 498.9–502.1 nm for υ″ = 0, and redder for hot bands.
  * **They would be in neither the line list nor a strict E′ > D_e(B) continuum.**
* **Preferred fix:** extend the discrete list to all bound levels, and add the true bound–free absorption for E′ > D_e(B).
  * This needs a long grid, because r_max(υ′ = 70) = 7.37 Å ([T73] Table IV) and grows quickly above that. A B-spline basis with non-uniform knots should make this cheap.
* **Fallback:** add the B contribution from the *pseudocontinuum* B potential for final energies E′ > E_cut(J′). E_cut(J′) is the highest B level in the line list for that J′.
  * Above D_e(B) this equals the true bound–free to ≤ 0.3% (§7.3).
  * Below D_e(B) it gives the reflection envelope of the missing near-threshold bands, smooth rather than resolved.
* **Wall mismatch.** Inside the Hannover data range (2.647–2.70 Å), the [T11C] wall sits 7–15 cm⁻¹, about 0.0006 Å, above Hannover (§7.5). At the dissociation-limit turning point (about 2.62 Å) the offset is 24 cm⁻¹.
  * I estimate a step of ≲ 1–2% in the Franck–Condon envelope across D_e(B). **This is an estimate; I have not computed it.**
  * By default, leave the wall unshifted, because that reproduces the fitted absorption.
  * The alternative is to lower A₀ by about 11 cm⁻¹ so the wall joins Hannover at 2.647 Å. That moves the 480–500-nm B←X by about 0.25 nm, i.e. about +1% at 480 nm.

### 9.5 A←X caveats

* **Long wavelengths.** λ ≳ 810–850 nm needs the real A wall beyond 2.806 Å (§8, item 2).
* **Bound–bound A←X is smeared.** The pseudocontinuum A potential smooths the bound–bound A←X absorption into a continuum. That absorption lies at λ > 803.8 nm for υ″ = 0 and further red for hot bands.
  * This is adequate for a background, because the A←X lines are weak and dense.
  * It will not reproduce individual A←X lines.

### 9.6 Accuracy to expect at 250–700 K

* **273–393 K, λ ≤ 500 nm:**
  * [T11C] claims about 0.5% for its total.
  * The wavelength scale matters: ε changes by 3.5%/nm at 480 nm. [T11C] does not say whether λ is air or vacuum. Reading it as air reproduces Table II and the 500.2-nm input better (§10). The two readings differ by 0.14 nm, i.e. about 0.5%.
  * The i2spec check matches [T11C] Table II to 0.3–0.7%.
* **A + C at 530–635 nm:** agreement with [T11B] ε_c is within 1–2σ, about ±2–4 L mol⁻¹ cm⁻¹.
* **Above 393 K the model is an extrapolation.** As υ″ rises, the potentials are used outside the fitted absorption window (2.6–2.8 Å). The linear transition moments are the weak point:
  * μ₁(C) = −0.16(11) D/Å gives ±6% in μ_C at |R − 2.7| = 0.25 Å, i.e. ±12% in ε. The sign of μ₁(C) even flips with the choice of B attachment point.
  * μ₁(A) = −0.048(33) gives ±3% in μ_A, i.e. ±6% in ε.
  * Expect the band centres to be good to a few percent. At 700 K the far wings (λ ≲ 420 nm and λ ≳ 800 nm) may be off by about 10–20%, until high-T data are available (§8, item 4).

### 9.7 Tests to add

| Test | Tolerance |
|---|---|
| V1: [T11C] Table II at 400–495 nm, 0 °C and 35 °C | 1% |
| A←X region: [T11C] Table II at 650–800 nm | 1.5% |
| V6: A + C (+ B bound–free) against [T11B] ε_c at 530–635 nm | 2σ |
| Pseudocontinuum-B mode (if implemented) against [T11C] Table II at 500–900 nm | 2–3% |
| G_ab convention (as a unit test) | ε_A and ε_C are unchanged under (G = 2, μ/√2) |

---

## 10. Implementation in i2spec (2026-09-15)

**Code.** `src/i2spec/continuum.py` (`continuum_model`, `Continuum`), `tests/test_continuum.py`, `prototypes/check_continuum.py`.

**Model.** As in §9.1, with these choices:
* **A and C:** Table I, with G_ab = 1. As in §7, A uses the Morse stand-in beyond 2.806 Å, and C a slope-matched R⁻ⁿ tail below 200 cm⁻¹.
* **B:** Hannover B, with the [T11C] wall inside their crossing at 2.7154 Å. μ_B is [T11B] eq. (10), the same function the line list uses.
* **The boundary with the line list (§9.4)** combines the preferred fix with the fallback:
  * The line list now runs to v′ ≤ 62. `intensity_model` uses nlev_b = 70, and the 7 Å box sets the limit.
  * For each J′ the list stops at the last bound level below the first box state. Quasi-bound levels above a box state are no longer listed as lines. `MasterLineList.upper_cut[J′]` records the energy half a level spacing above that last level.
  * The continuum counts B←X only above the cut. It uses continuum functions of the true B potential, not the pseudocontinuum, normalized with the WKB amplitude inside the B well (3.3–3.8 Å).
  * Below D_e(B) this gives the level-averaged absorption of the bands missing from the list. Above D_e(B) it gives the bound–free continuum, with centrifugal-barrier resonances averaged over.

**Numerics.**
* Numerov at h = 0.001 Å from 2.15 Å, vectorized over energy. The energy grids are 20 cm⁻¹ for A and C and 10 cm⁻¹ for B.
* X states come from the sinc-DVR on 2.15–4.5 Å (step 0.004 Å), and the overlaps use DVR quadrature on that grid.
* **Energy normalization uses no derivatives.** u = k^(−1/2)(a sin φ + b cos φ), with φ = ∫k dR, is fitted over a window, and then a² + b² = 1/(πc) is imposed. The windows are 4.5–5.0 Å for A and C and 3.3–3.8 Å for B.
* J″ is sampled every 5, with the midpoint rule for Σ(2J+1).
* **Convergence is better than 0.01%.** Checked against: J step 1; windows moved to 7.5–8 Å (A, C) and 3.6–4.1 Å (B); h halved.
* **Exact check:** energy-normalized Airy functions for a linear potential agree to 5×10⁻⁴.
* The continuum for T ≤ 600 K (79 J samples, v″ ≤ 29) takes about 10 s. It is cached together with the line list it complements.

**The §7 check was biased high.**
* It computed a² = u²k + u′²/k with u′ from a second-order finite difference. That underestimates u′² by (kh)²/6, so M² comes out too large by about (kh)²/12.
* For C←X this is 0.7% at 420 nm and 0.3% at 480 nm.
* Part of the 0.3–0.7% agreement in §7.2, and the "vacuum" reading in §7.7, came from this bias.

**Results** (`prototypes/check_continuum.py`):

| Check | Table λ read as vacuum | Table λ read as air |
|---|---|---|
| Table II, 405–495 nm, 0 °C | 0.981–0.995 | **0.995–0.999** |
| Table II, 390–495 nm, 35 °C (mean) | 0.988 | **0.996** |
| Table II, 700–900 nm, 0 °C / 35 °C (mean) | 1.010 / 1.006 | **1.005 / 1.002** |
| ε(500.2 nm, 35 °C), lines + continuum at 0.5–2 nm resolution; [T11C] 587.4(8) | 585.3 | **587.0** |
| ε(436 nm, 298 K); [T11C] 31.0(4) | 30.84 | **31.12** |

**Treat the [T11C] wavelengths as air.** Every check favours that reading, which is typical of a spectrophotometer scale. This is inference made here; the paper does not say.

Also:
* **A + C against the [T11B] ε_c at 309 K:** 0.97–1.05 at 530–600 nm, 0.91–0.94 at 610–635 nm. The §7 check gives the same, so this is the [T11C] fit's own residual.
* **Lines + continuum against [T11C] Table II (pseudocontinuum) at 505–645 nm**, both smoothed over 4 nm: 0.994–1.033.
  * So the discrete B←X is up to 3% stronger than [T11C]'s pseudocontinuum B at 550–605 nm. The discrete B←X uses the Hannover potentials and the [T11B] μ(R).
  * That is within the difference between Tellinghuisen's two B representations, which use different B walls and μ_B forms.
* **σ(500.0 nm air), 0.58 nm, 298 K:** 2.252×10⁻¹⁸ cm², i.e. ε = 589.
  * [T11C]'s room-temperature value is 590 ± 4.
  * Spietz 2006 gives 2.186(21), i.e. ε = 572. That 3% discrepancy is between the two papers themselves; see `spectra-validation.md`.
