# Physics and theory for modelling the I₂ B–X spectrum: a review

*Compiled 2026-09-15.*

**Structure.** §§0–3 are a synthesis, drawn from the topic reviews in appendices A–E and from the full
texts of Knöckel *et al.* (2004) and Salumbides *et al.* (2008). The appendices carry provenance tags:
[FT] full text, [ABS] abstract only, [SNIP] search-engine snippet only (provisional), [EST] an estimate
made here, with the working shown. Where the synthesis and an appendix disagree, the synthesis is more
recent; for example, the appendices list the Salumbides 2008 parameters as unavailable, but they have
since been read in the copy in the VU Amsterdam repository.

---

## 0. Bottom line

* **What the Hannover model is.** The model inside IodineSpec is:
  * a single-channel direct potential fit in Tiemann's "X-representation";
  * plus effective Born–Oppenheimer-correction (BOC) functions for the B state;
  * plus the Bodermann 2002 hyperfine interpolation formulae.

  Its chain of versions is Knöckel 2004, then Salumbides 2008, plus the Liao 2010 local NIR model.
* **Its weak points are where post-2004 physics helps:**
  * long-range tails spliced on by hand, with constants fixed from 1985–86;
  * effective BOC functions that cannot be split between upper and lower state, and whose α(R) disagrees with the measured B-state g_J factors;
  * hyperfine formulae that fail near 514 nm (up to 0.8 MHz; Yoshiki 2023) and above v′ ≈ 44;
  * a ~3 MHz long-period systematic in the energy scale, seen by an independent comb-FTS test (Reiners 2024).
* **No competing potential exists.** Our searches found no other I₂ X or B potential published since 2004, and no MLR, DELR or EMO fit of I₂ at all (App. A). Modern relativistic ab initio curves are nowhere near spectroscopic accuracy (Rₑ to 0.03–0.06 Å, ωₑ to 3–8 %). They are useful only as shape priors (App. C).
* **Recommendation (§3).**
  * A global direct potential fit with damped MLR potentials, whose long-range tails take their functional form from theory.
  * A common-asymptote constraint tying the B and X dissociation limits together.
  * Le Roy/Watson-style BO-breakdown functions.
  * A hyperfine Hamiltonian diagonalized level by level. Its parameters are vibrational expectation values of R-dependent functions, plus explicit second-order perturber terms.
  * Keep the Hannover ξ-form as an import path for validation.
  * A single-channel treatment is valid up to v′ ≈ 77. Coupled channels beyond that are a later extension.

---

## 1. What the Hannover model is (from the full texts)

### 1.1 Rovibronic part

| Element | Knöckel, Bodermann & Tiemann, EPJD 28, 199 (2004) | Salumbides et al., EPJD 47, 171 (2008) |
|---|---|---|
| Potential | V = Σ aᵢXⁱ with X = (R − Rm)/(R + b·Rm) | same form, refit |
| B state | Rm = 3.02669183 Å, b = −0.75, a₀…a₃₁ | same Rm and b, new a₀…a₃₁ |
| X state | Rm = 2.66638233 Å, b = −0.60, a₀…a₁₃ | same Rm and b, new a₀…a₁₃ |
| Inner wall (R < R_I) | A_I·exp(−B_I(R − R_I)); R_I = 2.647 Å (B), 2.40 Å (X); A_I and B_I re-derived from continuity after every change | same |
| Outer branch (R > R_O) | De − Σ Cₙ/Rⁿ plus an exponential; R_O = 4.9 Å (B), 3.30 Å (X) | De − Σ Cₙ/Rⁿ **− A_O·exp(−B_O(R − R_O))**. The exponential is only for smoothness: "does not have a physical meaning of an exchange term" |
| Long-range constants (fixed, not fitted) | De from GLL 1991. B: C₅, C₆, C₈, C₁₀ (Gerstenkorn, Luc & Amiot 1985). X: C₆, C₈, C₁₀ (Bacis, Cerny & Martin 1986), **no C₅** | same |
| BO corrections (B state only) | α(R) = (μ_ref/μ)(2Rm/(R+Rm)) Σ αᵢXⁱ with α₀–α₄. V_corr = (1 − μ_ref/μ)·V_ad, with V_ad = (2Rm/(R+Rm)) Σ vᵢXⁱ and v₀–v₂. β = 0 | α₀–α₅. V_ad = (2Rm/(R+Rm))^l Σ vᵢXⁱ with l = 5 and v₀–v₃. X-state BOCs were tested but cannot be told apart from B-state ones (correlation) |
| Solver and fit | Numerov (Blatt 1967); MINUIT χ²; data uncertainty floor 3 MHz | same; 2.1 MHz for frequency differences |
| Data | ≈1500 rovibronic frequencies (hyperfine removed): v″ ≤ 17, v′ ≤ 43, J ≤ 238 | > 1900, including ≈290 ¹²⁹I₂ and ≈90 ¹²⁷I¹²⁹I (v′ 9–21, v″ 0–5). Adds new comb lines at 661, 564, 585 and 735 nm, and high-J 0–15 lines |
| Turning points covered by data | B 2.650–4.589 Å; X 2.427–3.079 Å | same |
| Stated accuracy | 2σ < 3 MHz at 526–667 and 776–815 nm; ±30 MHz at 514–526 nm; ≤ 60 MHz at 667–776 nm | ≈1.5 MHz (1σ) within the quantum-number range of the data. Tₑ(¹²⁷I₂) − Tₑ(¹²⁹I₂) = 94(11) MHz |

**Parameter tables.** The 2004 Table 4 has misprints. The IQO web page says so, but its erratum link is dead. The misprints are in the X-state column and in the B-state A_O/B_O, as shown in `hannover-model-reproduction.md`. **Use the 2008 Table 1.** Our prototype reproduces that set: all four joins are continuous to 10⁻⁴ cm⁻¹, T₀ = 15724.5872 cm⁻¹, and the NIR bands agree with the 2004 local Dunham model at the few-MHz level.

### 1.2 Hyperfine part

The interpolation formulae come from Bodermann's 1998 thesis and were published in EPJD 19, 31 (2002). Details are in App. D, §4c.

**Formulae**
* **eqQ, both states:** Dunham-like, y_kl (v+½)ᵏ[J(J+1)]ˡ. 11 parameters, SD 12 kHz. Valid for v′ < 44, v″ < 18, J″ < 240.
* **C_B:** a single pole in the vibrational energy plus J-dependent terms. Valid for v′ = 0–36.
* **d_B and δ_B:** two poles plus a Gaussian bump. The bump mimics the B″1u curve crossing near 16 803 cm⁻¹. Valid for E(v′) < 19 500 cm⁻¹.
* **X-state C, d, δ:** fixed at 3.154, 1.524 and 3.705 kHz.
* **Isotopologues:** Salumbides et al. 2006.
* **Update, 2026-09-15:** the published 2002 paper has now been read, and its coefficients differ from the thesis values in App. D. For example, the quadratic term of eqQ_X is +0.4534×10⁻¹ MHz, not 4.5×10⁻³. The published formulae are implemented in `src/i2spec/hfs_params.py`. Eqs. (12)–(13) turn out to give δ_B and d_B directly. The operator conventions (EQ sign, TSS normalization) and the validation against BIPM are in `hannover-model-reproduction.md`.

**Documented failures**
* Near 514 nm (v′ = 44–45), Yoshiki 2023 found residual SDs of 78–387 kHz, with single components off by up to 0.8 MHz.
* The formulae do not cover v′ ≥ 45 except by extrapolation.
* The perturbed levels v′ = 57–60 and 76–78 are outside their scope.

### 1.3 Known limitations

From the papers:
* Bands 33–0, 34–0 and 35–0 have systematic residuals of up to 5 MHz (2004).
* α(R) is inconsistent with the measured B-state g_J factors (2004, Fig. 7). Upper- and lower-state BOCs cannot be separated (2008).
* The long-range tails are fixed, not fitted. The X tail has no C₅. B-state C₅ is fixed at 3.161×10⁵ cm⁻¹ Å⁵, although the literature spans 2.78–3.16×10⁵ (App. B).
* Data gaps (2004):
  * 667–776 nm and 514–526 nm;
  * v″ = 8–11 constrained only by the Gerstenkorn–Luc atlas;
  * no precise link between v″ = 0–7 and v″ = 12–17.

From later tests:
* **Reiners et al. 2024** (comb-FTS, 515–630 nm) found a ~2 m s⁻¹ (~3 MHz) long-period wave, plus an oscillation that follows the band structure. They conclude improvements should "mainly focus on the energy scale (potential functions and hyperfine interaction)".
* **Our own comparison:** the 2008 global potentials and the 2004 local NIR model differ by tens of MHz at the edges of the NIR set (0–12 at low J, 0–17 at high J). This is unresolved; test it against measured data.

---

## 2. Advances since 2004 and what they buy us

| Area | Advance | What it buys | Detail |
|---|---|---|---|
| Potential forms | Damped MLR (Le Roy et al. 2011), DELR, CPE (Brakmane 2026). MLR has predicted levels correctly across a 5000 cm⁻¹ data gap (Li₂ c state). The Tiemann ξ-form shows oscillation artifacts (K₂; Brakmane 2026) and poor near-threshold extrapolation (Falke 2008). | Better interpolation and extrapolation, fewer correlated parameters, correct tails | App. A |
| Asymptotes and long range | B-state C₅ varies by ~10 % with the fitting model, so fit it with a prior. D₀(X) values spread by 0.06 cm⁻¹, so float it under the common-asymptote constraint with ΔE_fs(I) = 7602.9762(65) cm⁻¹ (NIST). An X-state quadrupole–quadrupole C₅ exists physically but has no known value. | Extrapolation to high v′ and v″ | App. B |
| Hyperfine data and theory | B-state constants up to v′ = 70–82 (Chen, Cheng & Ye 2004; Pique 1986). An eqQ_B(R) function fitting v′ < 57 to 0.21 MHz SD (Chen 2004). Second-order perturber calculations (Chen, de Jong & Ye 2005). Comb-era constants (NMIJ/Yokohama 2016–2024). Q(¹²⁹I)/Q(¹²⁷I) = 0.701213(15) (Pyykkö 2018). Duo can do homonuclear hyperfine as of 2026. | Physically based hyperfine functions with honest uncertainties, replacing interpolation formulae | App. D, C |
| Near dissociation | Single-channel model plus perturbative hyperfine is valid up to v′ ≈ 77. Few-channel coupling is needed for v′ ≈ 78–81, and an atomic-basis coupled-channel model for v′ ≳ 81. | Sets the model's validity domain | App. B, D |
| BO breakdown and isotopes | Effective per-atom ũ(R) and q̃(R) (Le Roy 1999; Watson 2004). With only two isotopes, the adiabatic and field-shift terms are degenerate. Estimated field shift: 1–30 MHz in Tₑ. Atomic isotope shift at the asymptote < 30 MHz. | Isotopologue predictions with honest correlations | App. C |
| Intensities | μ(R) known to < 2 % over 520–640 nm (Tellinghuisen 2011). Strong J dependence of line strengths (Herman–Wallis). Continuum decomposition. Modern cross sections: σ(500 nm) = 2.186(21)×10⁻¹⁸ cm² (Spietz 2006). | Absorption spectra, not just line positions | App. E |
| Line shapes | Predissociation widths Γ(v, J, F) (Vigué 1981; Pique 1983). Pressure-broadening and shift data are sparse. Recoil and transit-time effects in sub-Doppler spectroscopy. | Realistic spectra and metrology systematics | App. E, D |
| Ab initio | SO-CI curves (Asano & Yabushita 2003). EFG gives Q(¹²⁷I) = −696(12) mb. No modern hyperfine, μ(R) or coupling functions. | Shape priors only | App. C |

---

## 3. Recommended physics model

### 3.1 Rovibronic levels

* **Fit scope:** one global direct potential fit of X and B for all three isotopologues, with mass scaling.
* **Primary potential form:** damped MLR.
  * The long-range function u_LR takes its form from theory:
    * B state: C₅, C₆, C₈ (and C₁₀), with priors from App. B.
    * X state: C₆, C₈, C₁₀, plus C₅ if a theoretical value can be found.
  * Dₑ and rₑ are explicit parameters.
* **Common asymptote:** Tₑ(B) + Dₑ(B) = Dₑ(X) + 7602.9762 cm⁻¹, plus the atomic isotope shift. That shift is below 30 MHz and is held at 0 at first. Dₑ(X) floats.
* **BO breakdown:**
  * Effective ũ_X(R) and ũ_B(R): adiabatic plus field shift, per atom.
  * Nonadiabatic centrifugal term q̃_B(R).
  * Use priors or regularization, because the upper- and lower-state functions are correlated.
  * Use the measured B-state g_J factors (Broyer et al. 1975) as an external check on q̃_B.
* **Hannover ξ-form as a second representation:**
  * to import the 2008 parameters and reproduce IodineSpec;
  * to cross-check bias from the choice of functional form. An optional spline or IPA potential can serve as a further bias diagnostic.
* **Validity:** single-channel for v′ ≤ 77. Flag the known local perturbations (v′ = 57–60, 76–78) and either exclude or model them. Near-dissociation coupled channels are deferred.
* **X-state range:** a single MLR also allows the high-v″ fluorescence data (v″ up to 107) to be included. Those were outside the Hannover X window.

### 3.2 Hyperfine structure

* **Hamiltonian:** H = eqQ·H_EQ + C·H_SR + d·H_TSS + δ·H_SSS (Broyer, Vigué & Lehmann 1978).
  * Diagonalize it in the |v J I F⟩ basis with ΔJ = 0, ±2 and I-mixing, using the model's own rovibrational energies.
  * Do not use first-order formulae.
* **Parameter functions**
  * **eqQ:** eqQ(v, J) = ⟨vJ|eqQ(R)|vJ⟩, with eqQ(R) fitted for each state, plus a small J-dependent second-order term. Chen 2004 showed this works for B to 0.2 MHz for v′ < 57.
  * **C, d, δ:**
    * First step: refit Bodermann-type energy-denominator (pole) forms to all comb-era data. This gives a quick improvement plus an uncertainty model.
    * Target: explicit second-order sums Σ_p |⟨vJ|f_p|v_pJ_p⟩|²/(E − E_p) over a few model perturber potentials (1u, 0g±, 1g), with separated-atom matrix elements scaled by fitted factors.
  * **X state:** Yokozeki & Muenter 1980 and Hong 2001 constants, plus eqQ_X(v, J).
  * **Isotopologues:** scale eqQ by the Q ratio (0.701213) and the magnetic terms by the ratio of nuclear moments. Treat ¹²⁷I¹²⁹I as heteronuclear, with no exchange symmetry.
* **Output for each component:** frequency with uncertainty; assignment (J, I, F, plus the BIPM aₙ label); relative intensity from the eigenvectors; natural width.

### 3.3 Intensities and line shapes

* **Line strengths**
  * J-dependent matrix elements ⟨v′J′|μ(R)|v″J″⟩, with μ(R) from Tellinghuisen 2011 and Lamrini 1994 at larger R.
  * Hönl–London factors for Ω = 0–0.
  * Nuclear-spin weights: 15:21 for ¹²⁷I₂, 28:36 for ¹²⁹I₂, a uniform 48 for ¹²⁷I¹²⁹I.
  * Boltzmann populations, with a partition function computed using the same weights.
* **Continuum**
  * Bound–free ¹Πu(1u)←X and A←X (Tellinghuisen 1973/1982/2011), plus the B←X continuum from hot bands.
  * Validate against Spietz 2006 and Saiz-Lopez 2004.
* **Line shape**
  * A Voigt profile for each hyperfine component, with the Doppler width from T.
  * Natural width Γ(v, J, F) from the Vigué predissociation model (C_v², a_v²) plus the radiative rate.
  * Pressure broadening and shift from measured coefficients where they exist (flagged as such), otherwise user input.
  * At low resolution, convolve the transmittance, not the cross section.

### 3.4 Open questions

1. The published Bodermann 2002 coefficients. The paper is paywalled, and the thesis values in App. D may differ.
2. Which ¹²⁷I mass convention the Hannover group used. It matters at about 1 MHz in the NIR.
3. The X-state C₅ (Saute & Aubert-Frécon 1982; value not retrieved), the atomic quadrupole moment Θ(I), and the primary value of the ¹²⁷I ²P₁/₂ hyperfine A constant.
4. How to include the high-v″ X-state fluorescence data (Martin et al. 1986, v″ ≤ 107). This requires a potential valid well beyond 3.08 Å.
5. The report of an ≈ 4 kHz linewidth near 508 nm (Cheng 2002) against the predissociation model.
6. Whether Tellinghuisen 2016 ("direct potential fitting RKR") includes I₂.

---

## Appendices: topic reviews

* A. Potential representations and catalog of published I₂ potentials
* B. Long-range physics and atomic asymptote data
* C. Relativistic ab initio electronic structure; BO breakdown, isotope and field-shift effects
* D. Hyperfine theory; B-state perturbations and predissociation
* E. Intensities and line shapes

---

## Appendix A — Potential representations and catalog of published I₂ potentials


No published MLR, DELR or EMO direct-potential fit (DPF) exists for I₂ X or B. The only post-2004 I₂ potentials found are the Hannover ones (2004, 2008) plus Tellinghuisen's 2011 inner-wall refinement. The recommendation is a damped MLR direct-potential fit for both states, with the Tiemann ξ-form kept as a compatibility path. Details and sources follow.

#### 0. Corrections to the brief
- **The Knöckel 2004 DOI in the brief is wrong.** `10.1140/epjd/e2003-00303-9` does not resolve in Crossref. The correct DOI is **10.1140/epjd/e2003-00313-4** (Eur. Phys. J. D 28, 199–209, 2004), confirmed in Crossref and OpenAlex. OpenAlex lists 87 papers citing it.
- **Searches for other post-2004 I₂ X/B potentials came up empty.** Searched:
  - all forward citations of Knöckel 2004 (87) and Salumbides 2008 (31);
  - post-2004 citers of Barrow–Yee 1973, Luc 1980, Gerstenkorn–Luc 1985 and Tromp–Le Roy 1985 (via Semantic Scholar);
  - Crossref and Semantic Scholar keyword searches ("direct potential fit iodine", "MLR iodine", "I2 B state potential", and similar).

  The only I₂ X/B potentials from after 2004 are the Hannover ones (Knöckel 2004; Salumbides 2008, with ¹²⁹I₂ and ¹²⁷I¹²⁹I) and Tellinghuisen's 2011 refinement of the B-state inner wall from bound–free absorption. One item is unresolved: Tellinghuisen 2016, "A direct potential fitting RKR method" (JMS 330, 20). It cites Gerstenkorn–Luc 1985, but whether it fits I₂ was not checked.

#### 1. Potential forms

**(a) Tiemann/Hannover analytic form.** The equations as given in Falke et al., PRA 78, 012503 (2008; arXiv:0804.2949). The same form appears in Samuelis 2000 and Allard 2002.
- ξ(R) = (R − R_m)/(R + b·R_m)
- Intermediate region: U_IR = Σ_{i=0}^{n} a_i ξ^i
- Inner wall (R < R_inn): U_SR = A + B/R^{N_s}. A and B are set for a continuous join at R_inn; N_s = 12 or 6 in the K₂ paper.
- Long range (R > R_out): U_LR = U_∞ − C₆/R⁶ − C₈/R⁸ − C₁₀/R¹⁰ ± A_ex R^γ e^{−βR}
- Born–Oppenheimer (BO) corrections: the effective Hamiltonian is H_eff = −(ħ²/2μ) d/dR[1+β(R)]d/dR + U + U_ad + ħ²[1+α(R)]J(J+1)/(2μR²), with β set to zero. The two correction functions are:
  - α(R) = (μ_ref/μ)·(2R_m/(R+R_m))·Σ α_i ξ^i
  - U_ad(R) = (1 − μ_ref/μ)·(2R_m/(R+R_m))⁶·Σ v_i ξ^i
- The 1998 I₂ version (Bodermann's dissertation) used b = −0.6 (X) and −0.75 (B), 13 X plus 31 B coefficients, and 6 BO-breakdown parameters.
- Behaviour: ξ → 1 as R → ∞, so any finite polynomial tends to a constant. The correct tail therefore has to be spliced on by hand at a hand-chosen R_out. High-order a_i are unconstrained outside the data window, and the pieces join only by continuity.

**(b) MLR (Morse/long-range).** Defined in Le Roy & Henderson 2007, extended in Le Roy et al. 2009, with damping added in Le Roy et al. 2011.
- V(r) = D_e[1 − (u_LR(r)/u_LR(r_e))·e^{−β(r)·y_p^{eq}(r)}]²
- y_p^{x}(r) = (r^p − r_x^p)/(r^p + r_x^p)
- β(r) = y_p^{ref}·β_∞ + (1 − y_p^{ref})·Σ β_i (y_q^{ref})^i, with β_∞ = ln(2D_e/u_LR(r_e))
- u_LR = Σ D_m(r)·C_m/r^m; the damping functions D_m were added in 2011.
- Behaviour:
  - At long range, V → D_e − u_LR + u_LR²/(4D_e), so the C_m, D_e and r_e appear as explicit parameters.
  - The 2009 paper extends u_LR to an eigenvalue of a small long-range matrix, for states whose long-range character changes with R (Li₂ A state).
  - The 2011 abstract says damping gives "much more realistic short-range extrapolation behaviour".
  - The 2009 paper is understood to require p > m_last − m_1; this is not verified here.
- Variant: MLR3 (Coxon & Hajigeorgiou 2010).

**(c) EMO (expanded Morse oscillator).** Lee et al., JMS 194, 197 (1999).
- V = D_e[1 − e^{−β(r)(r−r_e)}]², with β(r) = Σ β_i (y_p^{ref})^i
- Its tail is exponential rather than inverse-power, so it is wrong near dissociation. It is fine for the well region or as a starting guess.

**(d) DELR (double-exponential/long-range).** Huang & Le Roy 2003, introduced for barrier states. The equations below are as used in Klincare et al., JCP 160, 064307 (2024), consulted in full text.
- U = [T_dis + (Az − B)z] − U_LR, with z = e^{−β(R)(R−R_e)}
- A = D_e − U_LR(R_e) − U_LR′(R_e)/β(R_e); B = D_e − U_LR(R_e) + A
- β(R) = Σ β_i y^i, with y = ((R/R_ref)^p − 1)/((R/R_ref)^p + 1)
- The inverse-power tail is exact by construction.

**(e) Pointwise/spline potentials (IPA).** The inverse-perturbation idea comes from Kosman & Hinze 1975. The spline implementation is Pashov, Jastrzębski & Kowalczyk, CPC 128, 622 (2000); the regularized version is Grochola et al. 2004. betaFIT (Le Roy & Pashov 2017) converts pointwise curves to analytic forms.
- There is no functional bias, but there is no built-in extrapolation either: analytic extensions are needed at both ends, and regularization is needed where data are thin.

**(f) Near-dissociation expansions (NDE).** Limiting LeRoy–Bernstein (1970) theory, for V ≈ D − C_n/r^n:
- D − G(v) ∝ (v_D − v)^{2n/(n−2)}
- ΔG ∝ (v_D − v)^{(n+2)/(n−2)}
- B_v ∝ (v_D − v)^{4/(n−2)}
- Refinements: Padé-corrected NDEs (Le Roy & Lam 1980), Comparat 2004, generalized LeRoy–Bernstein (Sovkov & Ivanov 2014), and mixed polynomial-plus-NDE with a switching function (Tellinghuisen 2003, I₂ A state).
- These are formulas for level energies, not a potential, so they give no wavefunctions for Franck–Condon factors or hyperfine expectation values.
- Consistency check (mine): for n = 5, ΔG ∝ (v_D−v)^{7/3} and B_v ∝ (v_D−v)^{4/3}. This matches Tromp 1983 Eqs. 1–2 and Barrow–Yee's observation that ΔG^{10/7} is linear in E.

**(g) Dunham expansions** (Luc 1980; Gerstenkorn & Luc 1985).
- They are polynomials in (v+½), which cannot follow the (v_D−v)^{10/3} behaviour near the limit.
- Barrow & Yee (1973 abstract): B-state levels and centrifugal distortion constants "do not follow simple polynomials in (v+½)".
- Gerstenkorn & Luc (1985) note that Barrow–Yee's G(v) departs from Dunham above v′ ≈ 50.
- Extrapolating the minimum from Dunham coefficients is itself uncertain (Ilieva, Iliev & Pashov 2016).

**Newer forms (2010–2026)**
- Extended Lennard-Jones (Hajigeorgiou 2010, 2016) and compact 6–7-parameter models (Hajigeorgiou 2022). The 2022 average accuracy is about 0.06 % of D_e (~8 cm⁻¹), far too coarse for MHz work.
- **CPE**, presumably a Chebyshev-polynomial expansion (the abstract does not expand the acronym), was used alongside MLR in Brakmane et al., JCP 165, 044303 (2026), a K₂ X DPF covering the bound, quasi-bound and continuum ranges. The inner limb was extrapolated with CCSD(T)/CBS ab initio results.
- "Morphing" of ab initio potentials (Meuwly & Hutson 1999). Ab-initio-seeded analytic DPFs: Klincare 2024 (DELR), and MLR fits to ab initio points (Dattani, arXiv:1509.07041).
- **No** Gaussian-process or machine-learning potential fitted directly to diatomic line lists; the searches returned only polyatomic ab initio work.
- Reference implementations: dPotFit (Le Roy 2017) supports EMO, MLR, DELR and Šurkus polynomials; LEVEL (Le Roy 2017) solves the radial equation.

#### 2. Interpolation vs. extrapolation: published evidence
- **The Tiemann ξ-form has a documented artifact.** Brakmane et al. 2026 (K₂ X) report "non-physical oscillation-like behavior" when LIF frequencies are reproduced with the Tiemann et al. 2020 analytic potential (Phys. Rev. Research 2, 013366) and with Dunham constants. Their CPE and MLR fits avoid it to within ±0.002 cm⁻¹.
- **The Tiemann form did not extrapolate to new near-threshold levels at the 10⁻³ cm⁻¹ level.** Falke 2008 (K₂) found the previously published potentials missed the new levels by about 20× the experimental uncertainty. Differences between pairs of asymptotic levels were off by 20–40×.
- **MLR has a verified large-gap prediction.** A 2011 MLR fit to Li₂ c(1³Σg⁺) bridged a data gap of more than 5000 cm⁻¹. Its predictions in the middle of the gap were later found correct to about 1 cm⁻¹ (stated in arXiv:1509.07041; the data are Semczuk et al., PRA 87, 052505, 2013, with ±600 kHz uncertainty). The same arXiv paper shows raw ab initio points dipping below the theoretical long-range curve while the MLR behaves correctly.
- **MLR is more compact than Dunham.** Coxon & Hajigeorgiou 2010 (Cs₂ X, 99.24 % of the well) needed "significantly less than half" as many parameters as a Dunham analysis. They also found hyperfine effects not negligible near the limit: the fitted D_e is an *effective* limit lying between hyperfine limits. The same will matter for I₂ B near its limit.
- **NDE mechanical consistency has caught bad data in I₂.** Tromp 1983 fitted an NDE to v′ = 1–62 alone. It reproduced Barrow–Yee's data for v′ = 63–77 and the vibrational spacings up to v′ = 78, and showed that Danyluk–King's B_v values for v′ > 78 were wrong.
- **For the I₂ A state, a mixed polynomial-plus-NDE form beat a pure NDE** (Tellinghuisen 2003): lower χ² and "more realistic extrapolation of B_v … to dissociation".

**Synthesis:** extrapolation is best for forms that build in the theoretical inverse-power tail and a physically reasonable inner wall, with few, weakly correlated parameters. That ranks damped MLR, DELR and CPE-with-tail above the spliced Tiemann polynomial, which ranks above spline/IPA, which ranks above Dunham. NDEs extrapolate well in energy near the limit, but they are not potentials.

#### 3. Catalog of published I₂ X and B representations

| Authors (year) | State | Form | Data | v / J range | Stated accuracy | DOI |
|---|---|---|---|---|---|---|
| Richards & Barrow (1964) | X | Klein numerical integration | unknown | unknown | Binding exceeds a Morse curve beyond ~6.5 Å (van der Waals) | 10.1039/tf9646000797 |
| LeRoy (1970a,b) | X | Reassignment; RKR potential | not retrieved | unknown | unknown | 10.1063/1.1673357; 10.1063/1.1673358 |
| Le Roy & Bernstein (1971) | halogens incl. I₂ B | LeRoy–Bernstein NDE (D, C₅) | literature G(v) | high v | unknown | 10.1016/0022-2852(71)90046-4 |
| Barrow & Yee (1973) | B, X | Term values; RKR (B); C₅/C₆/C₈ fit | ~7000 lines | B 4–77; X 0–5 | Limit ±0.015–0.033 cm⁻¹ | 10.1039/f29736900684 |
| Wei & Tellinghuisen (1974) | B–X | "Best" constants | unknown | unknown | unknown | 10.1016/0022-2852(74)90239-2 |
| Danyluk & King (1977) | B near limit | Two-photon term values | — | v′ 77–82 | B_v for v′ > 78 judged erroneous (Tromp 1983) | 10.1016/0301-0104(77)85144-6 |
| King et al. (1980) | ¹²⁹I₂ B near limit | Two-photon | unknown | unknown | unknown | 10.1016/0301-0104(80)85100-7 |
| Luc (1980) | B–X | Dunham | 14 000 FTS lines, 139 bands | v′ 1–62, v″ 0–9, J ≤ ~150 | σ = 0.001/0.0017/0.004 cm⁻¹ for J ≤ 50/100/150 | 10.1016/0022-2852(80)90269-6 |
| Tromp, Le Roy, Gerstenkorn & Luc (1983) | B near limit | NDE + RKR | Luc; Barrow–Yee; Danyluk–King; Gerstenkorn–Luc (v′ 78–80) | v′ ≤ 82 | See key-numbers table | 10.1016/0022-2852(83)90027-9 |
| Tromp & Le Roy (1985) | B–X | NDE for G_v, B_v and centrifugal constants | 14 712 lines | unknown | unknown | 10.1016/0022-2852(85)90318-2 |
| Gerstenkorn & Luc (1985) | B–X | Dunham: 45 coefficients + 1 scale factor | ~100 000 FTS lines, 11 000–20 040 cm⁻¹ | X v″ ≤ 19; B v′ ≤ 80 (1.6 cm⁻¹ below limit) | σ = 0.002 cm⁻¹ | 10.1051/jphys:01985004606086700 |
| Martin, Bacis, Churassy & Vergès (1986) | X | Piecewise Dunham; D_e, R_e | B→X LIF-FTS | v″ ≤ 107 measured (108–113 extrapolated) | unknown | 10.1016/0022-2852(86)90254-7 |
| Ashmore & Tellinghuisen (1986) | X | Polynomial + NDE, smoothness-constrained | Fluorescence, "virtually entire well" | unknown | unknown | 10.1016/0022-2852(86)90202-x |
| Gerstenkorn, Luc & Le Roy (1991) | B–X of ¹²⁷I¹²⁹I, ¹²⁹I₂ | Dunham isotope scaling; RKR for centrifugal constants | Derived from 1985 ¹²⁷I₂ constants | X 0–19; B 0–80 | BO breakdown neglected | 10.1139/p91-194 |
| Bodermann (1998 thesis) | B, X | Tiemann ξ-form + BO-breakdown terms | Hannover | B 2.65–4.47 Å; X 2.427–3.04 Å | ~10 MHz | — |
| Knöckel, Bodermann & Tiemann (2004) | B, X | Hannover (your scope) | — | — | IodineSpec < 3 MHz at 526–667 nm (Reiners 2024) | 10.1140/epjd/e2003-00313-4 |
| Salumbides et al. (2008) | B, X incl. ¹²⁹I₂, ¹²⁷I¹²⁹I | Hannover + BO corrections (your scope) | — | — | — | 10.1140/epjd/e2008-00045-y |
| Tellinghuisen (2011) | A, B (small R), C | Inner walls from quantitative bound–free absorption | Low-resolution absorption | continuum | — | 10.1063/1.3616039 |

**Key numbers as stated by the sources**

| Quantity | Value | Source |
|---|---|---|
| B limit relative to X(v=0, J=0) | 20 043.208 ± 0.033 (extrapolation); 20 043.220 ± 0.015 cm⁻¹ (long-range fit) | Barrow & Yee 1973 |
| B limit relative to X(v=0, J=0) | 20 043.16(2) cm⁻¹ | Tromp 1983 |
| v_D (B) | Highest bound v′ = 87; v_D = 87.32(4) | Barrow & Yee 1973; Tromp 1983 |
| C₅ (B) | 2.88(3)×10⁵ cm⁻¹ Å⁵ | Tromp 1983 |
| Limiting NDE (B) | ΔG ≈ 0.00711(v_D−v)^{7/3}; B_v ≈ 0.000274(v_D−v)^{4/3} cm⁻¹ | Tromp 1983, Eqs. 1–2 |
| D₀ (X) | 12 440.1; 12 440.18(2); 12 440.239 cm⁻¹ | Barrow–Yee 1973; Tromp 1983; Gerstenkorn–Luc–Le Roy 1991 |
| D_e(X); D_e(B); T_e(B); T₀ | 12 547.340(6); 4381.249(1); 15 769.068(2); 15 724.587 cm⁻¹ | Gerstenkorn–Luc–Le Roy 1991 |
| I(²P₃/₂–²P₁/₂) splitting | 7602.977(2) cm⁻¹ | Gerstenkorn–Luc–Le Roy 1991 (adopted value) |
| X vibrational range | Measured to v″ = 107; last bound ≈ 113 | Martin 1986 abstract |

Consistency check: D₀(X) + splitting = 12 440.18 + 7602.977 = 20 043.157 cm⁻¹, matching the B limit of 20 043.16(2).

**Data gaps**
- **B state:** Dunham fits at FTS accuracy cover v′ 0–80. The last ~1.6 cm⁻¹ below the limit (v′ 81–87) is essentially unobserved; the two-photon data reach v′ = 82 but are unreliable above 78.
- **X state:** absorption constants cover v″ ≤ 19, and LIF reaches v″ = 107.
- **The 1998 Hannover X window is narrow (an estimate).** A Morse model with ω_e = 214.5 cm⁻¹ (Y₁₀ scaled by 1/ρ, ρ = 0.99611, from the 1991 tables), D_e = 12 547 cm⁻¹, μ = 63.452 u and r_e ≈ 2.666 Å (r_e not re-verified). This gives a = 0.121778·ω_e·√(μ/D_e) = 1.86 Å⁻¹. The window edges then fall at about 3150 cm⁻¹ (outer wall, 3.04 Å) and 3900 cm⁻¹ (inner wall, 2.427 Å), i.e. v″ ≈ 15–19. Anything beyond that, such as hot bands and fluorescence to high v″, is extrapolation.
- **MHz-level data are concentrated at 526–667 nm.** This follows from where IodineSpec quotes < 3 MHz (Reiners 2024). Reiners also saw a ~2 m s⁻¹ (≈3.6 MHz at 550 nm, conversion) long-period residual plus an oscillation that follows the band structure.
- **Newer near-IR data exist whose v-coverage is not mapped here:** 750–780 nm (JOSA B 2010), 915–985 nm (Nölleke 2018), 1053–1068 nm (2019), and 14 400–14 710 cm⁻¹ (Lefrán Torres 2022/2023, where IodineSpec5 reproduced the lines "very accurately").

#### 4. Recommendation for pyodine

**Primary representation:** a damped MLR for both X and B, in one global DPF across all isotopologues.
- u_LR comes from theory: for B, C₅, C₆, C₈ and possibly C₁₀. For X, u_LR might need to be a long-range matrix eigenvalue if several 0g⁺ curves share the asymptote.
- Tie the two asymptotes: T_e(B) + D_e(B) = D_e(X) + 7602.977 cm⁻¹. Falke 2008 stress that a common-asymptote constraint like this matters.
- Add Le Roy-type or Tiemann α(R)/U_ad(R) BO-breakdown functions.
- Start from RKR or the Hannover potentials.
- Constrain the inner wall with ab initio results (Brakmane 2026 approach) and with Tellinghuisen 2011, which pyodine will need for continuum cross sections.

**Secondary forms:**
- The Tiemann ξ-form, to import the Hannover potentials and reproduce IodineSpec for validation.
- CPE as a cross-check.
- Regularized IPA/spline, to diagnose bias from the functional form.

**Beyond single-channel:** for B-state levels with v′ ≳ 80, coupled channels are recommended, with atomic hyperfine structure rather than a single-channel potential, following the Coxon–Hajigeorgiou "effective D_e" lesson.

| Form | Pros | Cons |
|---|---|---|
| MLR (damped) | Correct tail; smooth and analytic everywhere; explicit D_e, r_e, C_n; few parameters; proven large-gap prediction; dPotFit available for validation | Nonlinear fit needs good starting values; β_i correlations; choice of p, q, r_ref; no prior I₂ use |
| Tiemann ξ | Continuity with IodineSpec; nearly linear fit; flexible | Piecewise with a hand-spliced tail; many parameters (31 for B); oscillation artifacts documented in K₂; poor extrapolation |
| DELR / CPE | Correct tail; DELR suits barrier states; CPE validated in 2026 | Less tooling; no advantage over MLR for deep wells |
| EMO | Simple | Exponential tail makes it unusable near dissociation |
| IPA / spline | No functional bias | No extrapolation; regularization required; noisy derivatives for R-dependent expectation values |
| NDE / Dunham | Compact; NDE exact near the limit | Not a potential; Dunham diverges near the limit; use only as legacy cross-checks |

#### References (each checked via Crossref or OpenAlex unless noted)

**Potential forms and methods**
- Samuelis et al., PRA 63, 012710 (2000), 10.1103/physreva.63.012710 (Crossref + abstract)
- Allard, Pashov, Knöckel, Tiemann, PRA 66, 042503 (2002), 10.1103/physreva.66.042503 (Crossref + abstract)
- Falke et al., PRA 78, 012503 (2008), 10.1103/physreva.78.012503 (OpenAlex; equations read in arXiv:0804.2949)
- Tiemann et al., Phys. Rev. Research 2, 013366 (2020), 10.1103/physrevresearch.2.013366 (OpenAlex)
- Le Roy & Henderson, Mol. Phys. 105, 663 (2007), 10.1080/00268970701241656
- Le Roy et al., JCP 131, 204309 (2009), 10.1063/1.3264688 (Crossref + abstract)
- Le Roy, Haugen, Tao, Li, Mol. Phys. 109, 435 (2011), 10.1080/00268976.2010.527304 (Crossref + abstract)
- Coxon & Hajigeorgiou, JCP 132, 094105 (2010), 10.1063/1.3319739 (Crossref + abstract)
- Lee et al., JMS 194, 197 (1999), 10.1006/jmsp.1998.7789
- Huang & Le Roy, JCP 119, 7398 (2003), 10.1063/1.1607313
- Klincare et al., JCP 160, 064307 (2024), 10.1063/5.0188443 (Crossref; text read in arXiv:2312.17005)
- Brakmane et al., JCP 165, 044303 (2026), 10.1063/5.0343162 (Crossref abstract)
- Pashov et al., CPC 128, 622 (2000), 10.1016/s0010-4655(00)00010-2
- Kosman & Hinze, JMS 56, 93 (1975), 10.1016/0022-2852(75)90206-4
- Grochola et al., JCP 121, 5754 (2004), 10.1063/1.1785782
- Le Roy & Pashov, JQSRT 186, 210 (2017), 10.1016/j.jqsrt.2016.03.036 (betaFIT)
- Le Roy, JQSRT 186, 179 (2017), 10.1016/j.jqsrt.2016.06.002 (dPotFit)
- Le Roy, JQSRT 186, 167 (2017), 10.1016/j.jqsrt.2016.05.028 (LEVEL)
- LeRoy & Bernstein, JCP 52, 3869 (1970), 10.1063/1.1673585
- Le Roy & Bernstein, JMS 37, 109 (1971), 10.1016/0022-2852(71)90046-4
- Le Roy & Lam, CPL 71, 544 (1980), 10.1016/0009-2614(80)80221-1
- Le Roy, JCP 73, 6003 (1980), 10.1063/1.440134
- Comparat, JCP 120, 1318 (2004), 10.1063/1.1626539
- Sovkov & Ivanov, JCP (2014), 10.1063/1.4869981
- Ilieva, Iliev, Pashov, JMS 330, 28 (2016), 10.1016/j.jms.2016.09.009
- Hajigeorgiou: JMS 263, 101 (2010), 10.1016/j.jms.2010.07.003; JMS 330, 4 (2016), 10.1016/j.jms.2016.06.014; Mol. Phys. 120 (2022), 10.1080/00268976.2022.2133754
- Meuwly & Hutson, JCP 110, 8338 (1999), 10.1063/1.478744
- Araújo & Ballester, IJQC 121 (2021), 10.1002/qua.26808
- Dattani, arXiv:1509.07041 (2015) (full text read)
- Semczuk et al., PRA 87, 052505 (2013), 10.1103/physreva.87.052505
- Yukiya et al., JMS 283, 32 (2013), 10.1016/j.jms.2012.12.006 (Br₂ A–X DPF; potential form not verified); IBr DPF, JMS (2025), 10.1016/j.jms.2025.112006
- Tellinghuisen, JCP 118, 3532 (2003), 10.1063/1.1539849
- Tellinghuisen, JMS 330, 20 (2016), 10.1016/j.jms.2016.09.008 (content unverified)

**I₂ spectroscopy and potentials**
- Verma, JCP 32, 738 (1960), 10.1063/1.1730793
- Richards & Barrow (1964), 10.1039/tf9646000797
- LeRoy (1970), 10.1063/1.1673357 and 10.1063/1.1673358
- Barrow & Yee (1973), 10.1039/f29736900684 (abstract)
- Yee, CPL 21, 334 (1973), 10.1016/0009-2614(73)80149-6
- Wei & Tellinghuisen (1974), 10.1016/0022-2852(74)90239-2
- Danyluk & King (1977), 10.1016/0301-0104(77)85144-6
- King et al. (1980), 10.1016/0301-0104(80)85100-7
- Luc (1980), 10.1016/0022-2852(80)90269-6 (abstract via search)
- Gerstenkorn, Luc, Sinzelle (1980), 10.1051/jphys:0198000410120141900
- Tromp et al. (1983), 10.1016/0022-2852(83)90027-9 (full text read)
- Gerstenkorn & Luc, Laser Chem. 1, 83 (1983), 10.1155/lc.1.83
- Tromp & Le Roy (1985), 10.1016/0022-2852(85)90318-2
- Gerstenkorn & Luc (1985), 10.1051/jphys:01985004606086700 (HAL full text)
- Martin et al. (1986), 10.1016/0022-2852(86)90254-7 (abstract via search snippet)
- Ashmore & Tellinghuisen (1986), 10.1016/0022-2852(86)90202-x
- Appadoo et al., JCP 104, 903 (1996), 10.1063/1.470814
- Gerstenkorn, Luc & Le Roy (1991), 10.1139/p91-194 (full text read)
- Knöckel et al. (2004), 10.1140/epjd/e2003-00313-4
- Salumbides et al. (2008), 10.1140/epjd/e2008-00045-y
- Tellinghuisen, JCP (2011), 10.1063/1.3555623 and 10.1063/1.3616039 (abstracts)
- Reiners et al., A&A 690, A210 (2024), 10.1051/0004-6361/202451389 (Crossref + arXiv HTML)
- Wang et al., PASP 132, 014503 (2019), 10.1088/1538-3873/ab5021 (its "2σ < 30 MHz" claim is unverified; the PDF was unreadable)
- Lefrán Torres et al., JMS 387, 111668 (2022), 10.1016/j.jms.2022.111668
- Rodríguez Fernández et al., JMS (2023), 10.1016/j.jms.2023.111789
- Nölleke et al., JMS 346, 19 (2018), 10.1016/j.jms.2017.12.013

Sources consulted in full: Tromp 1983, Gerstenkorn–Luc 1985, Gerstenkorn–Luc–Le Roy 1991, and the Klincare and Dattani arXiv papers.

---

## Appendix B — Long-range physics and atomic asymptote data

**Scope:** long-range physics of I₂ X(0g⁺) and B(0u⁺) — the Cₙ coefficients, dissociation energies, atomic fine and hyperfine data, LeRoy–Bernstein/NDE theory, and where hyperfine structure at the asymptote breaks the single-channel model.

**Provenance tags:** FT = full text read; ABS = abstract read (OpenAlex, Semantic Scholar, INIS or publisher); TAB = value read in a table of a later full-text paper that cites the source; SNIP = abstract seen only through a search-engine snippet (provisional); EST = an estimate made here, with working shown.

**Conventions:** V(R) = D − Σ Cₙ/Rⁿ, with C > 0 attractive. 1 E_h a₀ⁿ = 219474.63 × (0.529177 Å)ⁿ cm⁻¹, so 1 a.u. equals 32 523 (n=3), 9107.3 (n=5), 4819.4 (n=6), 1349.6 (n=8) and 377.9 (n=10) cm⁻¹ Åⁿ.

---

#### 1. Long-range coefficients

##### Physics

**B 0u⁺, asymptote I(²P₃/₂)+I(²P₁/₂)**
- **No C₃ term.** Both atomic levels are 5p⁵, so they have the same parity and no E1 transition moment connects them. There is no resonant dipole–dipole term.
- **No direct quadrupole–quadrupole term.** A J=1/2 atom has no static quadrupole moment, so the product Θ_a·Θ_b vanishes.
- **What produces C₅.** The leading first-order term is a resonant (exchange) quadrupole–quadrupole interaction. It couples the degenerate |a:3/2, b:1/2⟩ and |a:1/2, b:3/2⟩ configurations through the E2 matrix element ⟨²P₃/₂‖Q‖²P₁/₂⟩. E2 is allowed here (ΔJ=1, same parity).
- **g/u pairing.** The g and u combinations pick up ± this matrix element. So g/u partners should have C₅ of opposite sign, scaled by Ω-dependent angular factors.
- **Support in the literature:**
  - Tromp et al. 1983 [FT] state that "both theory and experiment agree that n = 5" for the B state. They also use Chang's (1967) prediction that the repulsive 1u state at the same limit has C₅(1u) = −(2/3)·C₅(B).
  - Saute & Aubert-Frécon (1982) [ABS] computed C₅ and C₆ for all 23 states as "quadrupolar electrostatic + dispersion" terms, valid for R > 7 Å.
  - Chen, de Jong & Ye 2005 [FT] use these curves for all ten states at the 3/2+1/2 limit.
- **Caveat:** the Saute and Chang full texts were not consulted. The explicit ±g/u sign statement is an interpretation of the mechanism, not a quoted result.

**X 0g⁺, asymptote 2 I(²P₃/₂)**
- Both atoms have static quadrupole moments, so first-order Θ–Θ gives a C₅ term.
- Near 5 Å the X state is mainly |m_a| = |m_b| = 1/2 in character (Martin et al. 1983 [ABS]).
- Saute's X-state C₅ value was not retrieved.
- The atomic quadrupole moment Θ of I(²P₃/₂) was not found. Medved', Fowler & Hutson 2000 cover only O, F, S, Cl, Se and Br (from its title).
- The Hannover X-state tail contains no C₅ term at all (Table 1).

**Other long-range terms [EST]**
- *Magnetic-dipole resonance.* The 1315 nm line is an M1 transition, so a resonant M1–M1 C₃ exists. Taking ~½ × (4/3)μ_B² with μ_B = α/2 a.u. gives C₃ ≈ 9×10⁻⁶ a.u. ≈ 0.3 cm⁻¹ Å³. At 20 Å this is about 1 MHz, against 0.1 cm⁻¹ from the C₅ term. It equals C₅/R⁵ only near 1000 Å, beyond the outer turning point of the last bound level (~90 Å). It is irrelevant for bound levels.
- *Retardation.* Casimir–Polder retardation of C₆ matters only at R of hundreds of Å, so it is negligible for bound levels.

**Mixed-n regime [EST, using the Gerstenkorn–Luc–Amiot 1985 constants]**
- The ratio of the C₆ term to the C₅ term is C₆/(C₅R) = 0.68 at 7 Å, 0.48 at 10 Å, 0.40 at 12 Å and 0.24 at 20 Å.
- The C₈ term is 8 % of the C₅ term at 10 Å.
- Levels v′ = 73–82 have outer turning points at 9–14 Å, so they are not in the pure n=5 regime. Any C₅ from a LeRoy–Bernstein fit is therefore an effective value.
- Le Roy 1980 [ABS] showed that level energies still follow the limiting law, but B_v values are sensitive to the higher Cₙ.

##### Table 1 — Cₙ values (cm⁻¹ Åⁿ, with a.u. in brackets)

| State (limit) | Coefficients | Method | Source | Tag |
|---|---|---|---|---|
| B (3/2+1/2) | C₅ = 3.11(20)e5 [34.1] | LB fit, 1931 bandheads, v=55–72 | Le Roy & Bernstein 1971 | TAB (Tromp 1983, Table III) |
| B | C₅ = 2.886(6)e5 [31.69] | linearized LB, v=64–77 | Barrow & Yee 1973 | TAB |
| B | C₅ = 2.776(18)e5 [30.48] | LB, v=72–82 (two-photon; calibration errors per Tromp) | Danyluk & King 1977 | TAB |
| B | C₅ = 2.850(10)e5 [31.29] | linearized LB, v=73–80 | Gerstenkorn & Luc 1983 | TAB |
| **B** | **C₅ = 2.884(31)e5 [31.67]** | nonlinear LB fit, v=73–80 (recommended) | Tromp et al. 1983 | FT |
| **B** | **C₅ = 3.161e5, C₆ = 1.506e6, C₈ = 2.480e7, C₁₀ = 4.200e8** [34.71; 312.5; 1.838e4; 1.111e6] | multipole fit to the outer branch of an IPA potential (v=0–80), truncated at C₁₀ | Gerstenkorn, Luc & Amiot 1985 | TAB (Salumbides 2008, Table 1); method from SNIP. **Fixed in the Hannover model.** |
| B | C₆, C₈: values not retrieved | RKR outer turning points; results depend strongly on the assumed C₁₀/C₈ ratio | Le Roy 1974 | ABS |
| B | \|C₅\| = 3.39e5 [37.2], C₆ = 1.79e6 [371] | theory: multipole (quadrupole + dispersion) | Saute & Aubert-Frécon 1982 (CPL) | SNIP; sign convention not verified |
| B | improved n=5 constants from mass-reduced LB on ¹²⁷I₂ + ¹²⁹I₂ (v=71–79); numbers not retrieved | LB | King et al. 1980 | ABS |
| B | NDE fit to 14 712 lines; numbers not retrieved | NDE | Tromp & Le Roy 1985 | metadata only |
| 1u (3/2+1/2) | C₅ = −(2/3)·C₅(B) | theory (Chang 1967) | as used by Tromp 1983 | FT |
| X (3/2+3/2) | C₆ = 1.1(1)e6 [228]; last bound level v″ = 114 | LB on v″ = 83–96 | Koffend, Bacis & Field 1979 | SNIP |
| **X** | **C₆ = 1.48e6, C₈ = 3.86e7, C₁₀ = 1.0e8, no C₅** [307; 2.86e4; 2.65e5] | LIF-FTS long-range tests, three isotopologues | Bacis, Cerny & Martin 1986 | TAB (Salumbides 2008). **Hannover X tail.** |
| X | C₅: unknown (computed by Saute 1982, value not retrieved) | theory | — | — |
| B″1u (3/2+3/2) | D_e = 202.53(4) cm⁻¹; C₅ = −0.49(28)e5 [−5.4]; C₆ = 1.00(33)e6 [207]; C₈ = 0.40(6)e8 [2.96e4] | E→B″ LIF-FTS | Inard et al. 1999 | SNIP |
| A 1u (3/2+3/2) | C₅ = 59 500 [6.53]; C₆ = 2.01e6 [417] | theory (Saute 1982) | as quoted by Appadoo et al. 1996 | FT |

**Search engines misattribute the B″1u coefficients.** Some snippets assign the B″1u values to Cerny et al. 1986 (an X-state paper). They belong to Inard et al. 1999. Cerny's X-state long-range constants remain unknown.

**Post-2004 results.** No new experimental or ab initio Cₙ for I₂ X or B. Chen et al. 2005 still use Saute's 1982 curves, and Salumbides et al. 2008 fix the 1985/1986 values. Two targeted searches for 2008–2025 relativistic calculations returned nothing. This is a search result, not proof of absence.

**Model spread.** Experimental B-state C₅ values range from 2.78 to 3.16×10⁵ depending on the model (LB-only versus a multipole fit through C₁₀). That ~10 % model dependence exceeds every quoted uncertainty. Theory gives 3.39×10⁵.

#### 2. Dissociation energies

| Quantity | Value (cm⁻¹) | Source | Tag |
|---|---|---|---|
| D_e(X) | 12 547.340(6), including the Kaiser correction (Y₀₀″ = −0.010) | Gerstenkorn, Luc & Le Roy 1991 | FT |
| D₀(X) for ¹²⁷I₂ / ¹²⁷I¹²⁹I / ¹²⁹I₂ | 12 440.239 / 12 440.654 / 12 441.072 | GLL 1991, Table 3 | FT |
| D_e(B) | 4381.249(1) | Gerstenkorn, Luc & Amiot 1985, as quoted by GLL 1991 | FT |
| D₀(B) for the same three isotopologues | 4318.629 / 4318.872 / 4319.115 | GLL 1991 | FT |
| T_e(B); T₀ | 15 769.068(2); 15 724.587 | GLL 1991 | FT |
| B limit above X(v=0, J=0) | 20 044.0(1.2) [LB 1971]; 20 043.208(33) [Barrow & Yee]; 20 043.063(20) [Danyluk & King, rejected]; 20 043.176(16) [Gerstenkorn & Luc 1983]; **20 043.159(19)** [Tromp 1983]; 20 043.216 [sum of GLL's D₀(X) + 7602.977] | as listed | FT/TAB/ABS |
| B asymptote in the Hannover model | 20 150.317 above the X minimum (= 12 547.340 + 7602.977) | Salumbides 2008 | FT |
| D₀(X), other determinations | 12 440.9(1.1) [Le Roy 1970]; 12 440.083(145) [Martin 1983, LB on X, a′0g⁺ and a1g]; 12 440.200(20) [Gerstenkorn & Luc 1983]; 12 440.18(2) [Tromp 1983]; 12 440.243 held fixed [Appadoo 1996] | as listed | ABS/FT |

**The "best" D₀(X) values disagree.** They span 0.06 cm⁻¹ (about 1.8 GHz), which is roughly 3σ of Tromp's quoted uncertainty. The Hannover model fixes the Gerstenkorn/Luc values.

#### 3. Atomic iodine

- **Fine structure ²P₁/₂–²P₃/₂**
  - 7602.9762(65) cm⁻¹ from NIST ASD (source refs L22426 and L368) [FT output].
  - 7602.977(2) cm⁻¹ from Luc-Koenig et al. 1973 (FTS), as quoted by GLL 1991 [FT].
  - This corresponds to λ_vac = 1315.27 nm.
- **Isotope shift of the 1.3 µm line (¹²⁹I − ¹²⁷I):** < 1 mK (< 30 MHz), Engleman et al. 1980 [ABS]. Salumbides 2008 say it is experimentally unknown and force their adiabatic Born–Oppenheimer correction to zero at the asymptote [FT].
- **²P₁/₂ radiative lifetime:** 140 ms (M1 transition; integrated cross-section 1050 ± 250 fm²), Ha et al. 1995 [ABS].
- **¹²⁷I (nuclear spin I = 5/2)**
  - ²P₃/₂: A = 27.59 mK = 827.1 MHz, B = 38.18 mK = 1144.6 MHz, as quoted by Vossel et al. 2025 [SNIP]. The primary source, Jaccarino et al. 1954 (atomic-beam magnetic resonance; significant nuclear octupole term), was not accessed.
  - ²P₁/₂ A: no primary value found. **EST:** scaling Engleman's ¹²⁹I A(½) by the ²P₃/₂ A ratio (1.5035) gives ≈ 220 mK ≈ 6.6 GHz, neglecting hyperfine anomaly.
- **¹²⁹I (nuclear spin I = 7/2)**
  - Consistency check: GLL count 28 + 36 = 64 = (2I+1)² molecular hyperfine components.
  - ²P₃/₂: A = 18.35(1) mK (550.1 MHz), B = 26.55(7) mK (795.9 MHz).
  - ²P₁/₂: A = 146.32(2) mK (4386.6 MHz).
  - Source: Engleman et al. 1980 [ABS].
- **Nuclear quadrupole moments**
  - Q(¹²⁷I) = −696(12) mb from molecular EFG calculations, van Stralen & Visscher 2003 [ABS].
  - The value Pyykkö 2018 adopts is unverified; snippets give both −696(12) and −688.22 mb.
  - Q(¹²⁹I) is unverified; a snippet gives a ratio of 0.701213(15), i.e. −483 mb.
  - The atomic B ratio 26.55/38.18 is 0.695.
- **Nuclear magnetic moments μ:** not retrieved.
- **Static dipole polarizability α(I):** 32.9 ± 1.3 a.u. recommended (experimental). Spin–orbit CI gives 34.6 (²P₃/₂, M_J-resolved), 35.1 (²P₁/₂) and 33.0 ± 1.7 (²P₃/₂). Source: Schwerdtfeger & Nagle 2019 [FT].
- **Atomic Θ(²P₃/₂) and quadrupole polarizability:** unknown. Vigué, Saute & Aubert-Frécon 1983 derived dynamic-polarizability matrix elements from I₂ C₆ values [ABS]; numbers not retrieved.

#### 4. LeRoy–Bernstein theory and near-dissociation expansions

**The LB law.** Le Roy & Bernstein (1970, CPL and JCP) and, independently, Stwalley (1970 CPL) showed that for V = D − Cₙ/Rⁿ:
- D − G(v) = X₀(n)·(v_D − v)^{2n/(n−2)}.
- The fit needs at least four levels near the limit. It predicts all unobserved levels with an uncertainty no larger than the binding energy of the highest measured level [ABS].

**Applied to the B state.** For n = 5 the exponent is 10/3, and Tromp gives X₀(5) = 9170.912/[μ⁵C₅²]^{1/3}. With μ = 63.4522378 u this yields X₀ = 2.081×10⁻³ cm⁻¹. A first-principles WKB period integral reproduces this value.

**NDE extensions:**
- Le Roy & Lam 1980: rational-polynomial NDE corrections to the LB law.
- Le Roy 1980 [ABS]: theory of deviations from the limit, tested on B-state I₂.
- Tromp & Le Roy 1985: NDE fit of the B–X line set.
- Ashmore & Tellinghuisen 1986: combined polynomial/NDE for X.
- Appadoo et al. 1996 [FT]: NDE for the A state (n = 5, theoretical C₅, D₀(X) fixed at 12 440.243).

**Limitations of the NDE approach for I₂:**
1. **Mixed n.** Because the C₆ term is still large at the relevant distances (Section 1), the fitted "C₅" is an effective value and D, v_D and C₅ are correlated. Linearized fits gave unrealistically small errors (Tromp).
2. **B_v bias from the separated-atom rotational energy.** At the limit, J_atomic = 2 and Ω = 0, so E = G + B_v[𝒥(𝒥+1)+6] (Tromp eq. 5). This biases the effective B_v, though the effect is small for v = 75–82.
3. **Coriolis and spin–rotation coupling** to the repulsive 1u state at the same limit are too small to matter. They shift B_v by 0.15–0.26×10⁻⁴ cm⁻¹, against B_v ≈ 20–37×10⁻⁴ cm⁻¹ [FT].
4. **Missing channels.** A single potential ignores hyperfine and fine-structure channel coupling, and possible BO breakdown near the limit. GLL 1991 fit the ¹²⁹I₂ near-limit data of King et al. only to 0.032 cm⁻¹ rms [FT].

**Highest observed levels**
- **B state, ¹²⁷I₂:**
  - v′ = 80 at 1.6 cm⁻¹ below the limit (Gerstenkorn & Luc 1983 [ABS]).
  - Danyluk & King reached v′ = 82 by two-photon spectroscopy, but their B_v values for v′ ≥ 79 are wrong because of calibration errors.
  - JILA sub-Doppler hyperfine spectra reach the limit (Cheng 2002, 523–498 nm [ABS]), and Chen et al. 2005 analyse levels 70 ≤ v′ ≤ 82.
- **B state, ¹²⁹I₂:** v′ = 71–79, within 2 cm⁻¹ of the limit (King et al. 1980).
- **Hannover 2008 fit:** B-state data cover v′ ≤ 43 (turning points 2.650–4.589 Å). Beyond R_O = 4.90 Å the tail is fixed at literature values. X-state data cover v″ ≤ 17, with the tail fixed beyond 3.30 Å.

##### Table 4 — B-state levels above the observed range [EST, Tromp constants, J = 0]

| v′ | E_b (cm⁻¹) | E_b (GHz) | ΔG (cm⁻¹) | R_out (Å) |
|---|---|---|---|---|
| 78 | 3.55 | 106 | 1.27 | 9.6 |
| 79 | 2.43 | 72.9 | 0.97 | 10.3 |
| 80 | 1.59 | 47.6 | 0.72 | 11.3 |
| 81 | 0.973 | 29.2 | 0.51 | 12.4 |
| 82 | 0.548 | 16.4 | 0.34 | 13.9 |
| 83 | 0.274 | 8.2 | 0.21 | 16.0 |
| 84 | 0.114 | 3.4 | 0.11 | 19.1 |
| 85 | 0.035 | 1.04 | 0.050 | 24 |
| 86 | 0.0053 | 0.16 | 0.013 | 35 |
| 87 | 5×10⁻⁵ | 0.0015 | — | ~90 |

So seven bound levels (v′ = 81–87) lie above v′ = 80, all within 1 cm⁻¹ of the limit.

**X state**
- The highest observed level is v″ = 107 (Martin et al. 1986, via Appadoo 1996 [FT]). Koffend et al. 1979 observed v″ = 83–96.
- [EST] Using LB with n = 6, C₆ = 1.1–1.48×10⁶ and v_D = 114–115:
  - v″ = 107 is bound by about 4–8 cm⁻¹.
  - v″ = 96 is bound by about 75–105 cm⁻¹.
  - About seven levels (v″ = 108–114) remain unobserved.

#### 5. Hyperfine asymptotes and the limit of the single-channel model [EST]

**Working.** Atomic hyperfine energies were computed as
E_F = (A/2)K + B·[(3/2)K(K+1) − 2I(I+1)J(J+1)] / [4I(2I−1)J(2J−1)], with K = F(F+1) − I(I+1) − J(J+1).

**Atomic splittings**

| Isotope, level | Constants used | Level energies | Span |
|---|---|---|---|
| ¹²⁷I ²P₃/₂ | A = 827.1, B = 1144.6 MHz | F = 1: −3541; F = 2: −2803; F = 3: −836; F = 4: +3388 MHz | 6.93 GHz (0.231 cm⁻¹) |
| ¹²⁷I ²P₁/₂ | A ≈ 6.60 GHz (estimated) | F = 2: −11.5 GHz; F = 3: +8.2 GHz | 19.8 GHz (0.66 cm⁻¹) |
| ¹²⁹I ²P₃/₂ | Engleman 1980 | — | 6.37 GHz |
| ¹²⁹I ²P₁/₂ | Engleman 1980 | — | 17.5 GHz |

**Molecular asymptotes (molecular hyperfine span at each dissociation limit)**
- ¹²⁷I, 3/2+1/2 (B-state limit): 8 (F_a, F_b) asymptotes spanning about 26.7 GHz (0.89 cm⁻¹).
- ¹²⁷I, 3/2+3/2 (X-state limit): span 13.9 GHz (0.46 cm⁻¹).
- ¹²⁹I, 3/2+1/2: span 23.9 GHz (0.80 cm⁻¹).

**Where the single-channel picture fails**

(a) **Vibrational spacing becomes comparable to the hyperfine splitting.**
- ΔG equals the ²P₁/₂ splitting at E_b ≈ 1.4 cm⁻¹ (v′ ≈ 80.3).
- ΔG equals the full asymptotic span at E_b ≈ 2.1 cm⁻¹ (v′ ≈ 79.3).
- Beyond this, the second-order (Broyer–Vigué–Lehmann) effective hyperfine Hamiltonian stops being a valid perturbation expansion.

(b) **Electronic (Ω) splitting becomes comparable to the hyperfine splitting.**
- The splitting between Ω-states of the same limit is about |ΔC₅|/R⁵.
- Using C₅(1u) = −(2/3)C₅(B), this is about (5/3)·E_b at the outer turning point: 2.6 cm⁻¹ at v′ = 80 and 0.9 cm⁻¹ at v′ = 82, the latter equal to the hyperfine span.
- From v′ ≈ 82 onward (E_b ≲ 0.5 cm⁻¹, R ≳ 14 Å), Hund's case (c) recouples toward atomic F-coupling. The g/u label then fails; hyperfine-induced g/u mixing near this limit is what Pique et al. 1984 reported.

(c) **Consistency with the literature.** Chen et al. 2005 [FT] find the B–1g(¹Πg) hyperfine mixing λ < 0.1 for v′ ≤ 78, but up to 0.4 for v′ = 79–82, with the perturbing 1g levels "typically a few GHz" away. They therefore restrict their second-order treatment to v′ < 78. Local 1g perturbations also occur much lower, at v′ = 57–60 (Chen 2004 [ABS]).

(d) **X state.** The same ΔG criterion with n = 6 puts the crossover at E_b ≈ 0.5 cm⁻¹. That affects only the last one or two levels, far beyond any data.

**Validity ranges for the B state:**
- **v′ ≲ 77 (E_b ≳ 4 cm⁻¹):** a single potential plus the perturbative hyperfine Hamiltonian is adequate.
- **v′ ≈ 78–81:** explicit few-channel coupling to the neighbouring 1g, 1u, 0g⁻ and 0u⁻ levels is needed.
- **v′ ≳ 81 (E_b ≲ 1 cm⁻¹):** a coupled-channel treatment in an atomic (F_a, F_b) asymptotic basis is needed. Options are Tiemann-style coupled channels, or the ERCAR asymptotic-basis method (Vossel et al. 2025, demonstrated on HI).

**Consequences for the model** (inputs to your recommendation section):
- The Hannover long-range tail is pure extrapolation from constants fixed in 1985/86, and the X tail has no C₅.
- The B-state C₅ is uncertain at about the 10 % level from model dependence. It should be fitted with a prior, not fixed.
- D₀(X) should float in the fit, with the B limit tied to it through the NIST fine-structure splitting.

**Open or unknown:**
- the primary ¹²⁷I ²P₁/₂ A constant, and confirmation of the ¹²⁷I ²P₃/₂ A/B;
- Pyykkö's 2018 Q values;
- Θ(I) and Saute's X-state C₅;
- the Tromp & Le Roy 1985 and Le Roy 1974 numbers;
- the 0.06 cm⁻¹ spread in D₀(X);
- whether any post-2004 ab initio Cₙ exist.

#### References (DOI; how verified)

1. Tromp, Le Roy, Gerstenkorn, Luc, J. Mol. Spectrosc. 100, 82 (1983). 10.1016/0022-2852(83)90027-9 — FT (open PDF).
2. Tromp & Le Roy, JMS 109, 352 (1985). 10.1016/0022-2852(85)90318-2 — Crossref metadata.
3. Gerstenkorn & Luc, Laser Chem. 1, 83 (1983). 10.1155/LC.1.83 — S2 abstract.
4. Gerstenkorn, Luc, Amiot, J. Phys. 46, 355 (1985). 10.1051/jphys:01985004603035500 — Crossref; values via ref. 6.
5. Gerstenkorn, Luc, Le Roy, Can. J. Phys. 69, 1299 (1991). 10.1139/p91-194 — FT.
6. Salumbides, Eikema, Ubachs, Hollenstein, Knöckel, Tiemann, Eur. Phys. J. D 47, 171 (2008). 10.1140/epjd/e2008-00045-y — FT (VU repository).
7. Knöckel, Bodermann, Tiemann, EPJD 28, 199 (2004). 10.1140/epjd/e2003-00313-4 — Crossref.
8. Barrow & Yee, J. Chem. Soc. Faraday Trans. 2 69, 684 (1973). 10.1039/F29736900684 — Crossref.
9. Le Roy & Bernstein, JCP 52, 3869 (1970). 10.1063/1.1673585 — abstract; also CPL 5, 42 (1970). 10.1016/0009-2614(70)80125-7 — Crossref.
10. Le Roy & Bernstein, JMS 37, 109 (1971). 10.1016/0022-2852(71)90046-4 — Crossref.
11. Stwalley, CPL 6, 241 (1970). 10.1016/0009-2614(70)80230-5 — Crossref.
12. Le Roy & Lam, CPL 71, 544 (1980). 10.1016/0009-2614(80)80221-1 — Crossref.
13. Le Roy, JCP 73, 6003 (1980). 10.1063/1.440134 — abstract.
14. Le Roy, Can. J. Phys. 52, 246 (1974). 10.1139/p74-035 — abstract.
15. Le Roy, JCP 52, 2678 (1970). 10.1063/1.1673357 — abstract; and JCP 52, 2683 (1970). 10.1063/1.1673358 — Crossref.
16. Danyluk & King, Chem. Phys. 25, 343 (1977). 10.1016/0301-0104(77)85144-6 — Crossref.
17. King, Littlewood, Robins, Wijeratne, Chem. Phys. 50, 291 (1980). 10.1016/0301-0104(80)85100-7 — INIS abstract.
18. Saute & Aubert-Frécon, JCP 77, 5639 (1982). 10.1063/1.443770 — abstract; and CPL 86, 59 (1982). 10.1016/0009-2614(82)83117-5 — Crossref (values via snippet only).
19. Vigué, Saute, Aubert-Frécon, JCP 78, 4544 (1983). 10.1063/1.445293 — abstract.
20. Chang, Rev. Mod. Phys. 39, 911 (1967). 10.1103/RevModPhys.39.911 — Crossref.
21. Koffend, Bacis, Field, JMS 77, 202 (1979). 10.1016/0022-2852(79)90102-4 — Crossref (values via snippet).
22. Martin, Churassy, Bacis, Field, Vergès, JCP 79, 3725 (1983). 10.1063/1.446293 — abstract.
23. Martin, Bacis, Churassy, Vergès, JMS 116, 71 (1986). 10.1016/0022-2852(86)90254-7 — Crossref.
24. Cerny, Bacis, Vergès, JMS 116, 458 (1986). 10.1016/0022-2852(86)90140-2 — Crossref.
25. Bacis, Cerny, Martin, JMS 118, 434 (1986). 10.1016/0022-2852(86)90180-3 — Crossref; values via ref. 6.
26. Ashmore & Tellinghuisen, JMS 119, 68 (1986). 10.1016/0022-2852(86)90202-X — Crossref.
27. Appadoo et al., JCP 104, 903 (1996). 10.1063/1.470814 — FT.
28. Inard, Cerny, Nota, Bacis, Churassy, Skorokhodov, Chem. Phys. 243, 305 (1999). 10.1016/S0301-0104(99)00077-4 — Crossref (values via snippet).
29. Chen, de Jong, Ye, JOSA B 22, 951 (2005). 10.1364/JOSAB.22.000951 — FT (JILA PDF).
30. Chen, Cheng, Ye, JOSA B 21, 820 (2004). 10.1364/JOSAB.21.000820 — abstract.
31. Cheng, Chen, Yoon, Hall, Ye, Opt. Lett. 27, 571 (2002). 10.1364/OL.27.000571 — abstract.
32. Pique, Hartmann, Bacis, Churassy, Koffend, PRL 52, 267 (1984). 10.1103/PhysRevLett.52.267 — Crossref.
33. Broyer, Vigué, Lehmann, J. Phys. 39, 591 (1978). 10.1051/jphys:01978003906059100 — Crossref.
34. Luc-Koenig, Morillon, Vergès, Physica 70, 175 (1973). 10.1016/0031-8914(73)90287-5 — Crossref.
35. Luc-Koenig, Morillon, Vergès, Phys. Scr. 12, 199 (1975). 10.1088/0031-8949/12/4/004 — abstract.
36. Jaccarino, King, Satten, Stroke, Phys. Rev. 94, 1798 (1954). 10.1103/PhysRev.94.1798 — Crossref only.
37. Engleman, Keller, Palmer, Appl. Opt. 19, 2767 (1980). 10.1364/AO.19.002767 — abstract.
38. Ha et al., Ber. Bunsenges. Phys. Chem. 99, 384 (1995). 10.1002/bbpc.19950990322 — abstract.
39. van Stralen & Visscher, Mol. Phys. 101, 2115 (2003). 10.1080/0026897031000109428 — abstract.
40. Pyykkö, Mol. Phys. 116, 1328 (2018). 10.1080/00268976.2018.1426131 — abstract only; iodine values unverified.
41. Schwerdtfeger & Nagle, Mol. Phys. 117, 1200 (2019). 10.1080/00268976.2018.1535143 — FT (Massey table).
42. Katsoprinakis et al., PRA 87, 040101 (2013). 10.1103/PhysRevA.87.040101 — FT (arXiv:1301.6947).
43. Vossel, Tsakontsis, Weike, Eisfeld, PCCP 27, 5043 (2025). 10.1039/D4CP04170D — Europe PMC abstract.
44. Medved', Fowler, Hutson, Mol. Phys. 98, 453 (2000). 10.1080/002689700162450 — Crossref (title shows no iodine).
45. Stone, At. Data Nucl. Data Tables 90, 75 (2005). 10.1016/j.adt.2005.04.001 — Crossref (values not retrieved).
46. Lukashov, Petrov, Pravilov, *The Iodine Molecule* (Springer, 2018). 10.1007/978-3-319-70072-4 — Crossref.
47. NIST Atomic Spectra Database, I I levels, https://physics.nist.gov/asd — queried 2026-09-14.

---

## Appendix C — Relativistic ab initio; BO breakdown, isotope and field-shift effects

Scope: relativistic ab initio I₂ electronic structure, and Born–Oppenheimer breakdown, isotope effects and the field shift.

**Four corrections to the brief**
- **JCP 145, 074104 (2016)** is Almoukhalalati, Knecht, Jensen, Dyall & Saue, "Electron correlation within the relativistic no-pair approximation" (10.1063/1.4959452). It is not about nuclear size. The nuclear-size paper is Almoukhalalati, **Shee** & Saue, PCCP 18, 15406 (2016).
- **No I₂ paper by Kokh, Alekseyev & Buenker was found.** Their closest work is on HI (2000) and Cl₂ (2004).
- **Chen, de Jong & Ye (2005) is not an ab initio EFG calculation.** It is second-order hyperfine perturbation theory using separated-atom wavefunctions plus de Jong's ab initio potential curves (details under 3d).
- **Watson 2004 (JMS):** its abstract does not mention the field shift. Whether the paper covers it is not confirmed.

---

#### TOPIC 3 — Relativistic ab initio electronic structure

##### 3a. Spin–orbit-coupled potential curves

**Earlier baselines**
- **Teichteil & Pélissier 1994:** pseudopotential plus spin–orbit CI for the I₂ excited states.
- **Visscher & Dyall 1996:** four-component Dirac–Coulomb(–Gaunt) benchmarks for the X state of F₂ through At₂, including CCSD(T).
- **de Jong, Visscher & Nieuwpoort 1997** (four-component MOLFDIR code):
  - X state: well reproduced at CCSD(T). Relativity and core–valence correlation are required; the Gaunt term is negligible.
  - Excited and ionized states: from relativistic CI/CC, described as "generally in good agreement" with experiment and with Teichteil & Pélissier.
  - Their numerical tables were not available.

**2000–2026**

| Work | Method | States | Reported accuracy vs experiment |
|---|---|---|---|
| Zaitsevskii, Pazyuk, Stolyarov, Teichteil & Vallet 2000 | Pseudopotentials + many-body multipartitioning perturbation theory | B–X, A–X, B′(1u)–X transition moments | Used to simulate absorption and B-state radiative rates; numbers unknown |
| Vala, Kosloff & Harvey 2001 | MRCI + relativistic ECP; some SO matrix elements; 108-state DIM model with atomic SO | X and excited states (also I₂⁻, I₃) | I₂ "very good agreement with experimental data" (abstract) |
| Visscher, Eliav & Kaldor 2001 | Relativistic Fock-space CC (first molecular implementation) | Pilot calculations on I₂ and HgH | "good and balanced description" (abstract); numbers unknown |
| **Asano & Yabushita 2003** (full text read) | Contracted SOCI over 16 spin-free states from I(²P)+I(²P); RECP with (4s4p2d1f1g) basis | X, A, B and five 1u states, plus 1u–1u radial derivative couplings | See the table below |
| Kalemos, Valdés & Prosmiti 2012 | Multireference methods + ECP | E ³Πg only | Double minimum: ion-pair global minimum plus Rydberg local minimum |
| Alekseev 2014 | Ab initio | Ion-pair states (I⁺ ³P/¹D + I⁻) and transition moments between them | "rather good agreement". A 2020 follow-up on IBr/ICl/BrCl reports ΔTe ≈ 100 cm⁻¹ and Re within 0.02 Å |
| Tashiro et al. 2014 | SO-CASPT2 | A, A′ and other excited states | Not a precision benchmark |
| Chattopadhyay, Mahapatra & Chaudhuri 2014; Ghosh et al. 2016 and 2017 | Multireference CC; four-component SSMRCC and IVO-SSMRPT | X only (Re, ωe, De) | "Close agreement"; relativity weakens the I–I bond; numbers unknown |
| Maurice et al. 2015 | Two-step SO-CI, benchmarked against experiment and four-component results | X and low states | Numbers unknown (paywalled) |
| Weike & Eisfeld 2024 | ERCAR diabatic model: SO as a sum of atomic operators | I₂ fine-structure states (proof of principle) | "good agreement" with SO ab initio |

**Asano & Yabushita 2003, quoted (calc / exp):**

| State | Re (bohr) | De (eV) | ωe (cm⁻¹) | Deviations |
|---|---|---|---|---|
| X | 5.099 / 5.038 | 1.431 / 1.543 | 208.9 / 214.5 | ΔRe = +0.032 Å, De −7%, ωe −2.6% |
| A | 6.002 / 5.885 | 0.1415 / 0.2033 | 80.44 / 92.9 | De −30% |
| B | 5.815 / 5.715 | 0.4386 / 0.5304 | 115.5 / 125.7 | ΔRe = +0.053 Å, De −17%, ωe −8% |

**Not found (2000–2026):** I₂-specific EOM-CC with spin–orbit coupling (Wang, Tu & Wang 2014 is the method paper only), KRCI, DMRG, and a four-component 23-state valence study of neutral I₂. The same group has published Cl₂ (2008), Br₂ (2014), At₂ (2021) and I₂⁻ (2023).

**Assessment:** the best published curves reach Re to about 0.03–0.06 Å, ωe to 3–8%, and De to 7–30%. That is useful for shapes and topology, and many orders of magnitude short of spectroscopic precision.

##### 3b. B–X transition-moment function μ(R)
- **Ab initio:** Zaitsevskii et al. 2000 (quasirelativistic). No modern high-level μ(R) was found. No quantitative comparison with experiment was found.
- **Experiment:**
  - Koffend, Bacis & Field 1979: optically-pumped-laser gain, R-centroid 2.8–4.6 Å.
  - Tellinghuisen 1982: transition strengths.
  - Lamrini et al. 1994: 2.633 ≤ R_c ≤ 6.035 Å; |μ̄e|² falls by about 4 orders of magnitude.
  - **Tellinghuisen 2011:** |μe|² to <2% relative standard error over 520–640 nm, by least-squares spectral simulation.
  - Tellinghuisen 2011b: extends μe(R) and fixes the small-R parts of the A and B potentials from continuum absorption.
- **Use:** the ab initio shape is a candidate prior only outside the measured R-centroid window.

##### 3c. Couplings behind B-state predissociation
- **Tellinghuisen 1972:** Franck–Condon analysis of B→1u(¹Π) predissociation. Rate negligible near v′ = 14, sharp peak near v′ ≈ 5, weaker maximum near v′ ≈ 26.
- **Vigué, Broyer & Lehmann 1981 (three papers):** theory and experiment for natural, hyperfine and magnetic predissociation, using separated-atom models.
- **Pazyuk et al. 2001:** B-state predissociation from a semi-empirical inverse atoms-in-molecule model. The abstract was not available.
- **Asano & Yabushita 2003:** ab initio radial derivative couplings among the 1u(1)–1u(5) states. The 1u(3)/1u(4) pair is a crossing type in I₂ with strong spin–orbit mixing in the Franck–Condon region.
- **Tellinghuisen 2011b:** B→C predissociation data constrain the C(¹Πu) potential at intermediate R. The abstract notes that no derived C potential yields (text truncated).
- **Gap:** no modern ab initio R-dependent B–B″(1u), B–0g⁻ or B–1g spin–orbit or rotational matrix elements were found.

##### 3d. Ab initio hyperfine parameters
- **van Stralen & Visscher 2003** (EFG at I in 9 molecules including I₂; numbers from van Stralen's thesis Tables 4.3–4.4):
  - I₂ EFG at Re: DC-HF 15.924 a.u., CCSD(T) correlation −0.942 a.u., total ≈14.98 a.u.
  - Q(¹²⁷I) from I₂ = −697.0 mb. Recommended −696(12) mb; CCSD(T) mean absolute deviation 5.4 mb (~0.8%).
  - Back-calculation: the implied X-state νQ ≈ −2.45 GHz.
- **Pyykkö 2018:**
  - Current standard Q(¹²⁷I) = −688.22 mb, from the atomic ground state.
  - **Q(¹²⁹I)/Q(¹²⁷I) = 0.701213(15)**, giving Q(¹²⁹I) = −483 mb.
- **Chen, de Jong & Ye 2005:** second-order C_B, δ_B, d_B and eqQ_B for 3 ≤ v′ ≤ 82.
  - Inputs: de Jong 1997 curves for the six perturbing states at R < 7 Å (shifted vertically to match), and Saute & Aubert-Frécon (1982) C₅/C₆ curves beyond 7 Å.
  - C_B agrees with experiment for R_c ≳ 4 Å with 1u mixing α = 0.99. The separated-atom model fails below ~4 Å.
  - δ(0u⁻) goes wrong beyond R_c ≈ 8 Å; replacing the 0u⁻ long-range curve with the 2g curve fixes it. The results are very sensitive to C₅/C₆ and to the 1g(¹Πg) long-range curve.
  - The 1g coupling is non-perturbative for v′ ≥ 79 (λ up to 0.4).
- **Earlier and current models:** Vigué et al. 1979 and Pique et al. 1984 used separated-atom models. Bacis et al. 1980 used LCAO/separated-atom models to explain X-state eQq(v).
- **Not found:** modern ab initio eQq(R) curves for the B state, and magnetic hyperfine or spin-rotation calculations.
- **New tool:** Yin, Yurchenko & Tennyson 2026 add hyperfine to Duo for homonuclear diatomics (five magnetic terms plus quadrupole, driven by R-dependent curves). It was validated on H₂, D₂, N₂ and N₂⁺, not I₂.

##### 3e. Using ab initio results as priors
Yes, but only as shape constraints; their absolute accuracy (above) is far below what the data need. Good places to use them:
- inner walls above the data;
- unobserved perturber states for second-order hyperfine and predissociation, as Chen 2005 did;
- μ(R) outside the measured window;
- the shape of dq/dR for eQq extrapolation.

Precedents:
- "Morphing" (scaling ab initio curves in R and energy): Meuwly & Hutson 1999; Špirko 2016 (BeH).
- Duo supports refining ab initio curves (Yurchenko et al. 2016).
- Caution: ab initio adiabatic corrections match empirical BO-breakdown functions at Re for LiK and LiRb, but their shapes disagree for Li₂, LiK and LiRb (Lutz & Hutson 2016). Use them as soft priors, not hard constraints.

---

#### TOPIC 6 — BO breakdown, isotope effects, field shift

##### 6f. Formalism and I₂ data
- **Formalism:**
  - Bunker 1968, 1972; Bunker & Moss 1977; Watson 1973, 1980.
  - Watson 2004 (JMS): BO breakdown represented by three functions Qi, Ri, Si per atom.
  - Watson 2004 (Can. J. Chem.): RKR-type inversion; the correction functions are **not fully determinable**.
  - Le Roy 1999: per-atom adiabatic ũ(R) and centrifugal q̃(R) functions.
- **I₂ data:**
  - Salumbides et al. 2006: hyperfine of ¹²⁹I₂ and ¹²⁷I¹²⁹I (eQq and spin-rotation), with models predicting all three isotopologues.
  - **Salumbides et al. 2008:** more than 380 ¹²⁹I₂ and ¹²⁷I¹²⁹I frequency differences relative to ¹²⁷I₂; direct potential fit with **effective BO-correction functions valid for all three isotopologues**. Functional form and magnitudes unknown (paywalled).
  - Lutz & Hutson 2016 cite Knöckel 2004 and Salumbides 2008 as "indications" of BO breakdown in I₂.
  - Earlier ¹²⁹I₂ hyperfine: Pique et al. 1980; Gläser et al. 1981.

##### 6g. Field shift
- **Theory:**
  - Tiemann, Knöckel & Schlembach 1982: finite nuclear size treated as a perturbation, giving δ⟨r²⟩-dependent terms.
  - Schlembach & Tiemann 1982: Pb/Tl rotational anomalies are field shifts.
  - Knöckel & Tiemann 1982 (PbS) and 1984 (PbO).
  - Knecht & Saue 2012: TlI, PbTe and PbS. The field shift there sits at Tl/Pb, not iodine.
  - Almoukhalalati, Shee & Saue 2016: four-component CCSD field-shift parameters. They confirm Tiemann's findings, reject a later scaling factor, and show the contact-density approximation errs by ~10%. The abstract names only PbS; I₂ and HI do not appear.
  - Lutz & Hutson 2016: δE = (2π/3)·Z·(e²/4πε₀)·|Ψ(0)|²·δ⟨r²⟩. Field shifts are "a few MHz" for K and Rb (mass shifts are tens of MHz), and dominate for Yb₂ (tens of MHz).
  - **Nothing found computes the field shift for I₂ or HI.**
- **Nuclear radii:**
  - Angeli & Marinova 2013 (IAEA table): R(¹²⁷I) = 4.7500(81) fm. **There is no ¹²⁹I entry.**
  - Neighbours from the same table: Te 126/128/130 = 4.7266/4.7346/4.7423 fm; Xe 126/128/130 = 4.7722/4.7774/4.7818 fm.
  - Engleman, Keller & Palmer 1980: the ¹²⁹I–¹²⁷I isotope shift of the 1.315 µm (²P₁/₂–²P₃/₂) line is **< 1 mK (< 30 MHz)**.

**Estimates made here (working shown; not from sources)**
1. **δ⟨r²⟩(127→129):** using 2R̄·ΔR, Te(126→128) = 0.076 and Te(128→130) = 0.073 fm²; Xe(126→128) = 0.050 and Xe(128→130) = 0.042 fm². Interpolating gives **≈0.06 ± 0.02 fm²**.
2. **Field-shift coefficient:** (2π/3)·53 = 111.0; 0.06 fm² = 2.14×10⁻¹¹ a₀². So δE = 2.38×10⁻⁹ Eh per a₀⁻³ = **15.6 MHz per a₀⁻³ per nucleus**, or 31 MHz for ¹²⁹I₂ vs ¹²⁷I₂.
   - Check: the same formula for Rb (Z = 37, δ⟨r²⟩ = −0.0362 fm²) gives 6.6 MHz per a₀⁻³. That matches Lutz & Hutson Fig. 5c (Δρ ≈ 1–1.3 a₀⁻³ ↔ ~7–9 MHz).
3. **Te of B–X:** Δρ for I₂ is unknown. Assuming bonding-induced Δρ ≈ 0.1–1 a₀⁻³ (Rb₂-like or smaller for p-bonding) gives **~1–30 MHz**. This is consistent with the atomic bound (< 30 MHz) at the dissociation limit.
4. **Normal mass shift of the atomic line, for comparison:** 7603 cm⁻¹ × 6.70×10⁻⁸ = 5.1×10⁻⁴ cm⁻¹ ≈ 15 MHz.
5. **Rotational and vibrational constants:** assume ΔV_FS varies by 10 MHz over 0.5 Å. With k(X) = ωe²/(2BeRe²) = 8.66×10⁴ cm⁻¹/Å², this gives δBe/Be ≈ 6×10⁻⁹ (≈ 6 Hz) and δωe ≈ 50 kHz. Both are negligible; only Te and levels near dissociation feel MHz effects.

**Answer:** yes, the field shift plausibly reaches the MHz level, up to tens of MHz, mainly in Te.

| Effect, ¹²⁹I₂ vs ¹²⁷I₂ | Size | Basis |
|---|---|---|
| Mass-scaled isotope shift, e.g. band (10,0) | (ρ−1)·[ωe′·10.5 − ωe″·0.5] = −0.00778 × 1212.6 cm⁻¹ ≈ **−9.4 cm⁻¹ ≈ −283 GHz** | Estimate; NIST WebBook constants |
| δ(me/M) per atom | 6.70×10⁻⁸ | Arithmetic |
| Adiabatic (DBOC) at Re(X) | ~−30 MHz for ⁸⁵Rb⁸⁷Rb (source). Scaled by δ(me/M) (1.34×10⁻⁷ vs 1.48×10⁻⁷): **~tens of MHz** (≈27 MHz if the electronic factor is Rb₂-like) | Estimate anchored on Lutz & Hutson |
| B-state nonadiabatic 0u⁺–1u rotational coupling | Assume ΔB/B ≈ 2B⟨L₊⟩²/ΔE with ΔE ≈ 3000–5000 cm⁻¹: ~10⁻⁵, i.e. 100–300 MHz at J = 100. A mass-independent effective potential mis-scales this by 1.55%: **~2–5 MHz** | Estimate (ΔE assumed) |
| Field shift in Te | **~1–30 MHz** | Estimate, above |
| Field shift in B and ωe | ~Hz and ~50 kHz | Estimate |
| Atomic asymptote shift (1.315 µm) | < 30 MHz total | Engleman 1980 |
| Q(¹²⁹I)/Q(¹²⁷I) | 0.701213(15) | Pyykkö 2018 |
| ¹²⁹I atomic hyperfine constants | A(³/₂) = 18.35(1) mK, B(³/₂) = 26.55(7) mK, A(¹/₂) = 146.32(2) mK | Engleman 1980 |

##### 6h. Consequences for mass-scaled combined fits
- **Adiabatic and field-shift terms cannot be separated.** With only two iodine isotopes, the adiabatic term (∝ δ(1/M)) and the field-shift term (∝ δ⟨r²⟩) are both a single number per isotope. Their R-functions are therefore perfectly degenerate, and any fitted "effective BO function" (such as Salumbides 2008) is their sum.
- **¹²⁷I¹²⁹I only tests additivity:** it should show half the ¹²⁹I₂ correction.
- **Recommended parameterization:**
  - one effective ũ(R) per state;
  - a q̃_B(R) term, which is needed because of the 1u coupling;
  - the atomic asymptote isotope shift, bounded below 30 MHz;
  - hyperfine scaled by the nuclear-moment ratios, applied to R-functions so that each isotopologue gets its own vibrational averaging.
- **Heteronuclear symmetry:** ¹²⁷I¹²⁹I is not strictly g/u symmetric and has no spin-statistics restriction. I have not quantified the effect of this.

---

#### References (✓ = verified: Crossref, OpenAlex or Semantic Scholar record, or full text read)
1. Teichteil & Pélissier, Chem. Phys. 180, 1 (1994). 10.1016/0301-0104(93)E0395-C ✓
2. Visscher & Dyall, JCP 104, 9040 (1996). 10.1063/1.471636 ✓
3. de Jong, Visscher & Nieuwpoort, JCP 107, 9046 (1997). 10.1063/1.475194 ✓
4. Zaitsevskii et al., Mol. Phys. 98, 1973 (2000). 10.1080/00268970009483400 ✓
5. Vala, Kosloff & Harvey, JCP 114, 7413 (2001). 10.1063/1.1361248 ✓
6. Visscher, Eliav & Kaldor, JCP 115, 9720 (2001). 10.1063/1.1415746 ✓
7. Pazyuk et al., Mol. Phys. 99, 91 (2001). 10.1080/00268970109483856 ✓
8. Asano & Yabushita, Bull. Korean Chem. Soc. 24, 703 (2003). 10.5012/BKCS.2003.24.6.703 ✓ (full text)
9. van Stralen & Visscher, Mol. Phys. 101, 2115 (2003). 10.1080/0026897031000109428 ✓ (values from van Stralen's PhD thesis on diracprogram.org)
10. Chen, de Jong & Ye, JOSA B 22, 951 (2005). 10.1364/JOSAB.22.000951 ✓ (full text)
11. Kalemos, Valdés & Prosmiti, JPCA 116, 2366 (2012). 10.1021/jp3000202 ✓
12. Alekseev, Opt. Spectrosc. 116, 329 (2014). 10.1134/S0030400X14030023 ✓
13. Alekseeva & Alekseev, Russ. J. Phys. Chem. A 94, 1382 (2020). 10.1134/S0036024420070043 ✓
14. Tashiro et al., PTEP 2014, 013B02. 10.1093/ptep/ptt118 ✓
15. Chattopadhyay, Mahapatra & Chaudhuri, Mol. Phys. 112, 2720 (2014). 10.1080/00268976.2014.906675 ✓
16. Wang, Tu & Wang, JCTC 10, 5567 (2014). 10.1021/ct500854m ✓
17. Maurice et al., JCP 142, 094305 (2015). 10.1063/1.4913738 ✓
18. Ghosh, Chaudhuri & Chattopadhyay, JCP 145 (2016). 10.1063/1.4962911 ✓
19. Ghosh et al., JPCA (2017). 10.1021/acs.jpca.6b11348 ✓
20. Weike & Eisfeld, JCP (2024). 10.1063/5.0191529 ✓
21. Tellinghuisen, JCP 57, 2397 (1972). 10.1063/1.1678600 ✓
22. Brewer & Tellinghuisen, JCP 56, 3929 (1972). 10.1063/1.1677797 ✓
23. Koffend, Bacis & Field, JCP 70, 2366 (1979). 10.1063/1.437744 ✓
24. Tellinghuisen, JCP 76, 4736 (1982). 10.1063/1.442791 ✓
25. Lamrini et al., JCP 100, 8780 (1994). 10.1063/1.466732 ✓
26. Tellinghuisen, JCP 134, 084301 (2011). 10.1063/1.3555623 ✓
27. Tellinghuisen, JCP 135 (2011). 10.1063/1.3616039 ✓
28. Vigué, Broyer & Lehmann, J. Phys. 42, 937 / 949 / 961 (1981). 10.1051/jphys:01981004207093700; …094900; …096100 ✓
29. Vigué, Broyer & Lehmann, PRL 42, 883 (1979). 10.1103/PhysRevLett.42.883 ✓
30. Broyer, Vigué & Lehmann, J. Phys. 39, 591 (1978). 10.1051/jphys:01978003906059100 ✓
31. Bacis et al., JCP 73, 2641 (1980). 10.1063/1.440477 ✓
32. Pique et al., JCP 80, 1390 (1984). 10.1063/1.446888 ✓
33. Saute & Aubert-Frécon, JCP 77, 5639 (1982). 10.1063/1.443770 ✓
34. Pyykkö, Mol. Phys. 116, 1328 (2018). 10.1080/00268976.2018.1426131 ✓ (full text)
35. Meuwly & Hutson, JCP 110, 8338 (1999). 10.1063/1.478744 ✓
36. Špirko, JMS 330, 89 (2016). 10.1016/j.jms.2016.08.009 ✓
37. Yurchenko et al., CPC 202, 262 (2016). 10.1016/j.cpc.2015.12.021 ✓
38. Yin, Yurchenko & Tennyson, JCP 165 (2026). 10.1063/5.0346970 ✓
39. Watson, JMS 45, 99 (1973). 10.1016/0022-2852(73)90179-3 ✓
40. Watson, JMS 80, 411 (1980). 10.1016/0022-2852(80)90152-6 ✓
41. Watson, JMS 223, 39 (2004). 10.1016/j.jms.2003.09.007 ✓
42. Watson, Can. J. Chem. 82, 820 (2004). 10.1139/v04-049 ✓
43. Bunker, JMS 28, 422 (1968). 10.1016/0022-2852(68)90176-8 ✓
44. Bunker, JMS 42, 478 (1972). 10.1016/0022-2852(72)90224-X ✓
45. Bunker & Moss, Mol. Phys. 33, 417 (1977). 10.1080/00268977700100351 ✓
46. Le Roy, JMS 194, 189 (1999). 10.1006/jmsp.1998.7786 ✓
47. Lutz & Hutson, JMS 330, 43 (2016). 10.1016/j.jms.2016.08.007 ✓ (arXiv:1608.02141, full text)
48. Tiemann, Knöckel & Schlembach, Ber. Bunsenges. 86, 821 (1982). 10.1002/bbpc.19820860910 ✓
49. Schlembach & Tiemann, Chem. Phys. 68, 21 (1982). 10.1016/0301-0104(82)85077-5 ✓
50. Knöckel & Tiemann, Chem. Phys. 68, 13 (1982). 10.1016/0301-0104(82)85076-3 ✓
51. Knöckel & Tiemann, CPL 104, 83 (1984). 10.1016/0009-2614(84)85309-9 ✓
52. Tiemann et al., Chem. Phys. 67, 133 (1982). 10.1016/0301-0104(82)85027-1 ✓
53. Knecht & Saue, Chem. Phys. 401, 103 (2012). 10.1016/j.chemphys.2011.10.030 ✓
54. Almoukhalalati, Shee & Saue, PCCP 18, 15406 (2016). 10.1039/C6CP01913G ✓
55. Almoukhalalati, Knecht, Jensen, Dyall & Saue, JCP 145, 074104 (2016). 10.1063/1.4959452 ✓ (electron correlation, not nuclear size)
56. Angeli & Marinova, ADNDT 99, 69 (2013). 10.1016/j.adt.2011.12.006 ✓ (IAEA charge_radii.csv)
57. Engleman, Keller & Palmer, Appl. Opt. 19, 2767 (1980). 10.1364/AO.19.002767 ✓
58. Salumbides et al., Mol. Phys. 104, 2641 (2006). 10.1080/00268970600747696 ✓
59. Salumbides et al., EPJD 47, 171 (2008). 10.1140/epjd/e2008-00045-y ✓ (abstract only)
60. Knöckel, Bodermann & Tiemann, EPJD 28, 199 (2004). 10.1140/epjd/e2003-00313-4 ✓
61. Bodermann, Knöckel & Tiemann, EPJD 19, 31 (2002). 10.1140/epjd/e20020052 ✓
62. Pique, Stoeckel & Hartmann, Opt. Commun. 33, 23 (1980). 10.1016/0030-4018(80)90085-1 ✓
63. Gläser et al., Opt. Commun. 38, 119 (1981). 10.1016/0030-4018(81)90212-1 ✓
64. NIST Chemistry WebBook, I₂ constants (Huber–Herzberg compilation): https://webbook.nist.gov/cgi/cbook.cgi?ID=C7553562&Mask=1000 ✓

**Not verified:** Knöckel, Kröckertskothen & Tiemann, Chem. Phys. 93, 349 (1985). It is cited by Lutz & Hutson, but has no Crossref match.

---

## Appendix D — Hyperfine theory; B-state perturbations and predissociation

Scope: I₂ hyperfine theory, including what the Bodermann 2002 formulae actually are, and B-state predissociation and perturbations.

**Sources consulted.** Full texts: Broyer 1978, Vigué 1981 I–III, Pique 1983, Pique 1986 I–II (HAL), Chen/Cheng/Ye 2004 (JILA PDF), Yoshiki 2023 (accepted manuscript), and Bodermann's 1998 dissertation (your file). Abstracts only: Hong 2001/2002, Yokozeki & Muenter 1980, Chen/de Jong/Ye 2005, Salumbides 2006, Capelle & Broida 1973, Paisner & Wallenstein 1974. **The EPJD 19, 31 (2002) paper is paywalled and has no abstract in OpenAlex or Semantic Scholar.** So every Bodermann coefficient below comes from the 1998 thesis, which is the precursor, and the published coefficients may differ.

---

#### 4a. The effective hyperfine Hamiltonian (Broyer, Vigué & Lehmann, J. Physique 39, 591–609, 1978)

**Structure**
- The Hamiltonian is H_hf = H_hf(1) + H_hf(2) + H_hf(1,2).
  - Only even electric and odd magnetic nuclear multipoles are allowed, with rank k ≤ 2I.
  - The basis is |Ω v J (I₁I₂) I F M_F⟩, with I = I₁ + I₂ and F = I + J, in Hund's case (a) or (c).
- **Key result:** in Ω = 0 states, the odd-rank electron–nucleus terms vanish to first order, because the 3j symbol (J k J; 0 0 0) = 0 for odd k. What survives at first order:
  - electric quadrupole, EQ (k = 2);
  - hexadecapole (k = 4), estimated at ≲1 kHz and neglected;
  - direct nucleus–nucleus terms: C_D ≈ 0.15 kHz and d_D ≈ 0.15 kHz (Broyer p. 605).
  - Consequently every observed magnetic constant in X and B is second order.

**Second-order terms**

V is the off-diagonal rotational (gyroscopic) operator, −(ħ²/2μR²)[J₊(L₋+S₋) + J₋(L₊+S₊)]. It couples 0u⁺ only to 1u states (1g for the X state).

| Product | Result |
|---|---|
| V × H_MD | C_E·(I·J) — only Ω = 1 perturbers of the same g/u symmetry |
| H_MD(1) × H_MD(2) | Tensor d_E and scalar δ·(I₁·I₂) — any Ω_p = 0 or 1, either g or u |
| V × H_EQ | A J-dependent contamination of "eqQ" |
| Same-nucleus products | Absorbed into eqQ |
| H_EQ(1) × H_EQ(2) | New constants e, f, h (ranks 0, 2, 4), of the same order as d and δ |
| V × H_MO | An effective octupole, (oO)_E ≈ 2×10⁻⁶·C_E ≈ 0.5 Hz at v′ = 43 |

- Because of the last row, Broyer judged the 1–2 kHz octupole reported by Hackel et al. (PRL 35, 568, 1975) "unrealistic".
- Yokozeki & Muenter (1980) needed neither octupole nor hexadecapole terms and gave only upper limits.
- Chen 2004 added e, f, g and hH at v′ = 57–65. The fit standard deviation did not improve, and the extra constants came out poorly determined with alternating signs.

**Diagnostic relations** (Broyer Table IV; thesis eqs. 5.12–5.13)
- For a single perturber with Ω_p = 1, δ = 2d_E. For Ω_p = 0, δ = −d_E.
- The sign of δ tells whether the perturber is g or u.
- At B v′ = 43, C_E/d_E ≈ 1.8 and δ/d_E ≈ +0.03, so more than one perturber is involved; Broyer proposed 1u plus 0g⁻.

**The operator used ever since**
- H = eqQ·H_EQ + C·H_SR + d·H_TSS + δ·H_SSS, with the geometric factors g(I, I′, J, J′, F) tabulated by Broyer.
- EQ and tensor spin–spin are off-diagonal in both I (ΔI = 0, ±2) and J (ΔJ = 0, ±2).
- Modern codes (Bordé's FORTRAN, used by Chen 2004) diagonalize a basis spanning ΔJ = 0, ±2, ±4 built on E_v, B_v, D_v, H_v, L_v, M_v. The within-state ΔJ mixing is therefore exact, not folded into the constants.

**Nuclear-spin symmetry for ¹²⁷I₂** (I = 5/2, fermion; thesis §5.1.3)
- Total I = 0, 2, 4 are antisymmetric spin functions; I = 1, 3, 5 are symmetric.
- X(0g⁺): even J″ goes with even I, giving 1+5+9 = 15 sublevels; odd J″ goes with odd I, giving 3+7+11 = 21.
- B(0u⁺): the pairing is reversed.
- Electric-dipole transitions keep ΔI = 0, so main lines (ΔF = ΔJ) have 15 or 21 components. The ortho:para weight is 21:15.
- Pique 1986 I:
  - any g/u hyperfine coupling requires odd ΔI (±1, ±3) — a "flip-flop" of one nuclear spin relative to the other — and therefore destroys the ortho/para labels;
  - g–g and u–u couplings have even ΔI.
- Arithmetic: ¹²⁹I₂ (I = 7/2) has 28 or 36 sublevels with weights 9:7; ¹²⁷I¹²⁹I has no restriction. Salumbides 2006 models both.

#### 4b. Values and their v, J dependence

**Table H1.** eqQ in MHz; C, d, δ in kHz unless noted.

| Level | eqQ | C | d | δ | Source |
|---|---|---|---|---|---|
| X v″=0, J″=13 | −2452.5837(16) | 3.162(8) | 1.58(5) | 3.66(3) | Yokozeki & Muenter 1980 (molecular-beam magnetic resonance) |
| X v″=0, J dependence | −2452.556(2) − 1.64(5)×10⁻⁴ J(J+1) − 5(2)×10⁻⁹ J²(J+1)² | | | | Hong et al. 2001b (JOSA B 18, 379) |
| X values held fixed in modern fits | from formula | 3.154 | 1.524 | 3.705 | Bodermann 2002, as used by Yoshiki 2023 and Chen 2004 |
| B v′=11 | — | C_E 28.8(1.4) | d_E −25.6(2.8) | −11.7(1.6) | Landsberg 1976, via Broyer 1978 Table VI |
| B v′=32 | −544.049(14) − 2.110(43)×10⁻⁴ J′(J′+1) | | | | Hong et al. 2001a (JOSA B 18, 1416) |
| B v′=43, J′=12 | −558.669(8) | 190.13(12) | −100.2(7) | 0.2(4) | Yokozeki & Muenter 1980 (reanalysis of Hackel 1975) |
| B P(13) 43–0 | −558.613(18) | 190.361(78) | −98.99(62) | −0.83(56) | Chen 2004 Table 1 |
| B v′=44, J′-dependent | −559.680(10) − 2.03(7)×10⁻⁴ J′(J′+1) | 205.39(2) + 2.290(13)×10⁻³ J′(J′+1) | −108.7(2) − 1.69(15)×10⁻³ J′(J′+1) | 2.06(6) + 6.3(4)×10⁻⁴ J′(J′+1) | Yoshiki 2023 |
| B P(19) 49–0 | −564.6327(72) | 309.250(21) | −169.23(43) | 26.31(31) | Chen 2004 |
| B P(21) 60–0 | −569.937(20) | 783.323(37) | −377.3(1.3) | 346.12(86) | Chen 2004 |
| B P(35) 70–0 | −557.32(24) | 2087.35(37) | −809(12) | 2833.5(8.8) | Chen 2004 |
| B v′=71 | −555(2) | 2160(70) | −900(100) | 3560(40) | Pique 1986 II (listed in MHz there) |
| B v′=75 | −537(5) | 3650(180) | −2700(200) | 5680(80) | Pique 1986 II |
| B v′=82 | −290(20) | 16000(1000) | −50000(8000) | 98000(10000) | Pique 1986 II |
| B separated-atom first-order eqQ′₀ | −573 | | | | Pique 1986 II |

**Physical basis of the v, J dependence**
- **eqQ** is mostly the first-order term, −eQ⟨q(R)⟩. Its v and J dependence comes from the R-dependence of the electric field gradient, sampled through anharmonicity and centrifugal stretching.
  - The J(J+1) coefficients at v′ = 32 and 44 agree, at about −2.0×10⁻⁴ MHz.
  - Chen 2004 inverted the data into eqQ_B(R) = Σ aᵢ(R−R_e)ⁱ over 3–5 Å (69 levels with v′ < 57; SD 0.21 MHz). The constant term a₀ = −487.80(17) MHz agrees with Bodermann's Dunham-type y₀₀(B) = −487.806(59) MHz.
  - Second-order pieces:
    - V×H_EQ, which scales with J; the thesis bounds the crossing-1u contribution at ≲400 kHz for J′ = 200;
    - near-limit Ω = 1 perturbers, which drive eqQ′ from −555 to −290 MHz over v′ = 71–82 (Pique).
  - eqQ_B reaches an extremum near v′ ≈ 60 and then reverses (Chen 2004).
- **C** is purely second order. Only Ω = 1 perturbers of the same g/u symmetry contribute: two shallow 1u states for B, and the a(1g) state for X.
  - Because the overlap comes from near the outer turning point, C ≈ α + β/(E − E_asymptote). An early form quoted by Broyer is C(kHz) = 1.3×10⁵/(E_c − E_v) with E_c = 4400 cm⁻¹.
  - The J dependence enters through J-dependent Franck–Condon overlaps and energy denominators.
- **d and δ** get contributions from every perturber (0g±, 0u⁻, 1g′, 1g″, 1u′, 1u″).
  - Near the limit the net effect looks like a single "effective 0g" perturber (δ > 0, d < 0).
  - Around v′ ≈ 43, Ω = 1 contributions cancel δ_B to about 0 while d_B stays near −100 kHz.
  - Between E ≈ 16000 and 17000 cm⁻¹, d_B and δ_B are modulated by the level shift from the crossing B″1u state. Its form follows an Ai·Bi Airy-function product (thesis eqs. 5.18 and 8.17). C_B shows no such modulation, which means ⟨0u⁺|V|1u⟩ ≪ ⟨0u⁺|H_MD|1u⟩. This fits the predissociation ratio Y² = a_v²/C_v² ≈ 1120 (see 5f).
- **X state:** C_X rises steeply near the X dissociation limit because of a(1g). Thesis fits:
  - global: C_X = 1.699(20) + 17932(64)/(12440.18 − G(v″)) kHz;
  - local, for 6 < v″ < 20: C_X = 3.162 − 3.313(85)×10⁻⁵ G + 6.335(28)×10⁻⁸ G² kHz.

#### 4c. What the Bodermann interpolation formulae actually are (thesis 1998, §§5.1.2 and 8.3–8.4)

**Data**
- Bodermann's own measurements: v′ = 0–3 and v″ = 12–17 (NIR saturation spectroscopy, 3f technique).
- Literature values of ΔeqQ, ΔC, Δd and Δδ.
- The functional forms follow from the physics in 4b.

**eqQ, both states fitted jointly**
- Form: Σ y_kl (v+½)ᵏ[J(J+1)]ˡ. 73 data points, 11 parameters, SD 12.3 kHz. Valid for v′ < 44, v″ < 18, J″ < 240.
- B state (MHz): y₀₀ = −487.806(59), y₁₀ = −1.8644(28), y₃₀ = 1.2507(88)×10⁻⁴, y₀₁ = −1.66(14)×10⁻⁴, y₁₁ = −2.08(15)×10⁻⁶, y₀₂ = −2.47(33)×10⁻¹⁰.
- X state (MHz): y₀₀ = −2452.3281(73), y₁₀ = −0.533(13), y₂₀ = 4.519(65)×10⁻³, y₀₁ = −2.22(15)×10⁻⁴, y₁₁ = 6.58(35)×10⁻⁶.
- Predicted ΔeqQ uncertainty is typically 10–20 kHz, and below 50 kHz for strong bands at J = 30–110.
- An alternative, eqQ = eQ⟨v,J|Σ qₙfⁿ|v,J⟩ (eq. 5.22), was set aside because users would need wavefunctions. The thesis says the fitted rotational dependence "cannot be separated" in this model.
- **Check made here (not from a source):**
  - B state, v′ = 32, J = 0: −544.11 MHz, close to Hong's value.
  - B state, (43, 12): −558.65 MHz, close to Yokozeki & Muenter.
  - X state, v″ = 0, J″ = 13: −2452.633 MHz versus the measured −2452.5837 MHz. That 50 kHz offset is larger than the fit SD, so either I misread the OCR'd table or a convention differs. **Use the published 2002 coefficients.**

**C_B** (eq. 8.20, restricted to v′ = 0–36)
- C_B = 2.93(41) + 73900(900)/[E(v′) − 19684.6(9.1)] − 0.0508(93)·B_v′·J′(J′+1) + 3.30(59)×10⁻⁶·B_v′·J′(J′+1)·E(v′) kHz, with E and B in cm⁻¹ from Gerstenkorn & Luc.
- The sign of the pole term is garbled in the OCR; C_B > 0 requires 73900/(19684.6 − E).
- Why the cutoff: a single pole cannot describe both the R < 4.2 Å (v′ ≤ 36) and R > 4.2 Å regions. Bodermann built an "effective 1u potential" from the C_B data (thesis Fig. 8.10); it resembles 1u′ from the Saute & Aubert-Frécon C₅/C₆ long-range curves but lies lower.
- Performance: residuals 0.87 kHz; maximum deviation below 2.5 kHz for v′ ≤ 36, compared with 8.5 kHz for Razet's formula and 15 kHz for Arie's. v′ = 43 is described poorly.

**d_B and δ_B** (eqs. 8.21–8.22, for E(v′) < 19500 cm⁻¹)
- The fit uses:
  - two poles at E_D1 = 19873.529 and E_D2 = 20614.963 cm⁻¹, standing for effective Ω = 0 and Ω = 1 perturber groups, with pole strengths b = −30416(96) and c = 117100(202) kHz·cm⁻¹;
  - offsets a = 30.770(10) and a′ = 18.230(40) kHz;
  - a Gaussian bump of amplitude −23.090(72) kHz, centred at 16803.47 cm⁻¹ with width parameter f = 2.857×10⁵ cm⁻², standing in for the B″1u crossing.
- About 20 values of each were fitted; SD 0.8 kHz. Predicted 1σ is 1.5 kHz for δ and 2 kHz for d, with maximum deviations of 3 and 3.8 kHz.
- Above 19500 cm⁻¹ the thesis says the potentials of all six perturbing states would be needed.

**Claimed accuracy**
- Thesis abstract: about 30 kHz for splittings over 530–820 nm.
- EPJD 2002: a search-engine snippet (not the abstract itself) says ≤30 kHz for F − J = 0 components over 514–820 nm.

**Documented failures**
- Yoshiki 2023 compared the 2002 formulas with new 514 nm data:
  - P(34)44–0: splitting residuals SD 78 kHz, versus 7.6 kHz for their own local J-fit.
  - R(58)45–0 (J′ = 59): SD 387 kHz, maximum 803 kHz.
  - Formula versus fitted constants: eqQ′ −561.609(50) vs −561.487(4) MHz; C′ 235(2) vs 231.687(3) kHz; d′ −118(3) vs −124.50(19) kHz; δ′ 13(4) vs 7.2(2) kHz. Each is 1.5–2.5σ of the formula's own uncertainty.
  - Interpretation: the ΔC of 3.3 kHz is multiplied by ⟨I·J⟩ ~ 10², which yields the hundreds-of-kHz errors.
- By construction, the formulas leave out:
  - v′ ≥ 45 (below ~514 nm), except by extrapolation;
  - the perturbed levels v′ = 57–60 and 76–78;
  - band heads near dissociation.

**Formulas they replaced** (thesis Table 8.8)
- Arie & Byer (v′ = 26–62): ΔeqQ ±2 MHz, ΔC ±6 kHz, Δd = −½ΔC ±5 kHz.
- Razet & Picard 1997 (v′ = 6–43): ±200 kHz and ±3 kHz.

#### 4d. Second-order and strong hyperfine effects

**(i) Within the B state (ΔJ = ±2)**
- At ordinary J this mixing is handled by diagonalization.
- **Order-of-magnitude estimate (not from a source):**
  - Shift ~ |V_J,J±2|²/[B(4J+6)], assuming |V| ≈ 0.1–0.3|eqQ|.
  - X state, J = 50, B ≈ 1.12 GHz: denominator ≈ 231 GHz, so the shift is (245–735 MHz)²/231 GHz ≈ 0.3–2.3 MHz, and it depends on F.
  - Diagonalizing is therefore mandatory at the kHz level.
- Near the dissociation limit, B_v′ → 0. At v′ = 75, B = 154 MHz while the hyperfine structure spans about 2 GHz (Pique 1986 II). J′ = 0 and 2 mix, producing "superhyperfine" structure in band heads and shifting the J′ = 0 centre of gravity. Line-position shifts stay below 0.010 cm⁻¹ up to v′ = 82.

**(ii) Coupling to other electronic states and g/u mixing** (Pique, PRL 52, 267, 1984; J. Physique 47, 1909 and 1917, 1986)
- All nine other states at the ²P₃/₂ + ²P₁/₂ limit couple to B:
  - through V⁰ (gyroscopic), V¹ (magnetic dipole) and V² (quadrupole);
  - with matrix elements taken from separated-atom wavefunctions and atomic hyperfine constants.
  - Ω = 2 states contribute negligibly.
- Second-order sums over perturbers reproduce C′ quantitatively for v′ = 71–82 (theory 2.20 → 16.8 MHz, experiment 2.16 → 16 MHz), and eqQ′, d′ and δ′ semi-quantitatively (Tables I–IV).
- v′ = 76: every effective parameter except C′ is linear in J′(J′+1), and the δ′/d′ slope ratio is about 2, which points to a 1g perturber.
- v′ = 77 and 78 required a 45×45 matrix (B + 1g; ΔJ ≤ 2; three I values; c± components). This gave extra, formally forbidden 1g–X lines near the 78–0 band head: the first observation of g/u symmetry breaking in a homonuclear molecule.
- Size of the effect at v′ = 60: Chen 2004 finds ΔE = 39 MHz, H₁₂ = 54 MHz, and mixing coefficients α = 0.82, β = 0.57 for P(84)60–0. With a Franck–Condon overlap of 0.1, the electronic matrix element is about 540 MHz.
- Near the limit at low J′, δ′ dominates and I′ becomes nearly a good quantum number.
- Related work in other I₂ states: hyperfine coupling between the D0u⁺ and β1g ion-pair states (Baturo 2018), and among valence states at the ²P₁/₂ + ²P₁/₂ limit (Baturo 2016).

#### 4e. Coupled-channel models and R-dependent hyperfine functions

**What exists**
- **No coupled-channel radial treatment of I₂ B-state hyperfine structure near dissociation** was found in the Semantic Scholar forward citations of Pique 1986 I, Chen 2004 and Bodermann 2002 through 2026.
- The closest work:
  - Pique 1986: separated-atom couplings, long-range ab initio curves (Saute & Aubert-Frécon 1982) completed with Morse functions, Airy-function overlaps, and full-matrix diagonalization for B + 1g.
  - Chen/de Jong/Ye 2005: each of the six perturbing states' contributions to C_B, δ_B, d_B and eqQ_B computed from available potential curves and separated-atom wavefunctions.
  - Chen 2004: R-functions obtained by inverting data through LEVEL. C_B(R) holds to ±2% (relative) over v′ = 3–70. eqQ, d and δ as R-functions fail above v′ ≥ 56 because of the 1g perturbation.

**Assessment**
- Pure R-function expectation values are adequate for eqQ.
- For C, d and δ they are not adequate at the kHz level: ±2% of about 200 kHz is roughly 4 kHz, worse than Bodermann's ≤2.5 kHz. These constants depend on perturber energy denominators, not just on R.
- A physical replacement for the formulas that can realistically be built:
  - eqQ(R) expectation values (first order),
  - plus explicit second-order sums Σ_p⟨vJ|f_p(R)|v_p J_p⟩²/(E − E_p) over model perturber potentials (ab initio long-range part plus fitted parameters), with the separated-atom matrix elements scaled by a few fitted factors.
  - Bodermann's pole formulas are effectively a one-pole approximation of exactly this sum.
- A full Tiemann-style coupled-channel model (compare the NaRb X/a work, Pashov et al. PRA 72, 062505, 2005) is the correct description for v′ ≳ 70 and for the resonances at v′ = 57–60 and 76–78. But it faces:
  - **Channel-count estimate:** 16 Ω-components × 5 J values × 6 I values ÷ 2 ≈ 240 channels per (F, parity).
  - Nine perturber potentials that are known only roughly.
  - R-dependent electronic hyperfine matrix elements that are known only in the separated-atom limit.
  - This makes it a research project, not a drop-in replacement.
- Pique's model achieved only about MHz-level agreement.

#### 5f. B-state predissociation

**Rate formula** (Vigué 1981 II eq. 10; quadrupole predissociation neglected, so I is conserved)
- Γ_vJIF = [C_v² J(J+1) + (a_v²/3)(I(I+1) + [3(I·J)² + (3/2)I·J − I(I+1)J(J+1)] / [(2J−1)(2J+3)]) − √2 a_v C_v (I·J)] × [1 + p_v J(J+1) + q_v (J(J+1))²]
- Here I·J = [F(F+1) − I(I+1) − J(J+1)]/2, and the total decay rate is Γ = Γ_rad + Γ_coll + Γ_pred.
- In a magnetic field, α_v²B² terms and M_F-odd interference terms appear, so lifetimes depend on M_F.

**Mechanism**
- The B″1Πu state (1u, correlating with ²P₃/₂ + ²P₃/₂) crosses B at R ≈ 2.89 Å, near v′ ≈ 5.
- In the gyroscopic coupling, the L and S matrix elements (+0.378ε and −0.363ε) nearly cancel, so C_v is small.
- Hence Y² = a_v²/C_v² = 1120 ± 160 (paper II) and 1180 ± 120 (paper III).
- Scaling with Franck–Condon density (FCD, in cm): C_v² = 5.8×10⁵ FCD, a_v² = 6.7×10⁸ FCD, α_v² = 3.1×10¹⁰ FCD, in units of s⁻¹, s⁻¹ and s⁻¹T⁻².
- The FCD has maxima near v′ ≈ 5 and v′ ≈ 25 and a node near v′ ≈ 13–14.
- Tellinghuisen's 1972 estimates (C_v² ≈ 1700 s⁻¹ at v′ ≈ 5 and 450 s⁻¹ at v′ ≈ 25) were 5–10 times too large.
- Pique 1983, v′ = 43: Γ_rad = 0.314(18)×10⁶ s⁻¹, C_v = 5.85(65) s⁻½, a_v = −142(5) s⁻½, from lifetimes of individual hyperfine sublevels.

**Table P1.** From Vigué 1981 III Tables V–VIII (read from rendered pages).

| v′ | Γ_rad (10⁶ s⁻¹) | C_v² (s⁻¹) | a_v² (10³ s⁻¹) | α_v² (10⁵ s⁻¹T⁻²) | Γ_tot (10⁶ s⁻¹) |
|---|---|---|---|---|---|
| 5 | — | — | — | 174(10) | — |
| 7 | 1.05(9) | 287(21) | 324(26) | 141(2) | — |
| 8 | 0.95(8) | 187(13) | 232(18) | 96(3) | — |
| 10 | 0.78(9) | 95(7) | 110(17) | 41(1) | — |
| 11 | 0.75(10) | 51(7) | 40(15) | 20.8(2) | — |
| 13 | 0.78(6) | 3.0(5) | 3(3) | 1.10(12) | — |
| 14 | 0.76(3) | <2 | 4(+10/−4) | 0.54(1) | — |
| 18 | 0.78(6) | 33(4) | 51(8) | 28.5(1.2) | — |
| 21 | 0.70(25) | 64(9) | 88(50) | 43.3(8) | 1.82(65); Paisner 1.45(4) |
| 24–27 | 0.61 (interpolated) | 80(20) | — | 44.6(1.2) | — |
| 32 | — | 52(28) (v′ = 30–34) | — | 28.3(3.6) (v′ = 30–34) | 0.94(13); 0.92(3) |
| 40 | 0.38(10) | 26(10) | — | — | 0.62(10); 0.70(2) |
| 43 | 0.35(5) | — | 12(5) | — | 0.43(12); 0.44(2) |
| 62 | ≤0.08(1) | — | — | — | 0.083(14); 0.11(2) |
| 70 | ≤0.04(1) | — | — | — | 0.038(12) |

Capelle & Broida (1973) measured lifetimes from under 0.4 μs to over 7 μs across v′ ≈ 5–70.

**Effect on linewidths (estimates, Δν = Γ_tot/2π)**
- v′ = 21: 230–290 kHz. v′ = 32: about 150 kHz. v′ = 43: about 70 kHz. v′ = 62: 13–18 kHz. v′ = 70: about 6 kHz.
- At v′ = 7–10 and high J, predissociation dominates:
  - v′ = 7, J = 100: C_v²J(J+1) ≈ 2.9×10⁶ s⁻¹.
  - Hyperfine components of a single rovibrational level differ in width by up to ~10× (Vigué II Fig. 5b: γ spans roughly 1–15×10⁶ s⁻¹ across F at v′ = 7, J = 88).
- Magnetic predissociation is negligible in Earth's field (≈0.04 s⁻¹ at v′ ≈ 5, B = 50 μT). It becomes significant, about 10⁵ s⁻¹, only at 0.1 T.
- **Unresolved:** Cheng et al. 2002 report a narrowest linewidth of about 4 kHz near 508 nm, which does not fit the Γ_tot/2π estimates above. The paper was not available to check how they define linewidth.

#### 5g. Local perturbations of the B state

**Table G1**

| Region | Perturber | Effect | Source |
|---|---|---|---|
| v′ ≈ 3–17, strongest 5–13 | B″1u crossing near v′ ≈ 5 | Predissociation; F-dependent lifetimes; irregularities in d_B and δ_B; J-dependent rovibronic residuals of several MHz at v′ = 7–10 | Vigué 1981; Chen 2004; thesis ch. 9 |
| v′ = 57–60 (v′ = 59 J′ ≈ 22; P(69)58–0, R(18)59–0, P(77)60–0, P(84)60–0; also P(63)70–0) | c1g (1Πg) rotational coincidences | Irregular eqQ_B, d_B, δ_B (C_B unaffected, since gyroscopic coupling cannot reach 1g); g/u mixing; extra lines; fit SD 0.3–1 MHz or no fit possible | Chen 2004; Jewsbury 1993 |
| v′ = 76–78 | "1g" | Strong g/u mixing; quasi-resonances at J′ = 15 (c₋) and 36 (c₊) in v′ = 77; forbidden lines at the 78–0 band head | Pique 1984, 1986 |
| v′ ≥ 71 band heads | Within B, ΔJ = 2 | Superhyperfine structure; shifts up to 0.010 cm⁻¹ | Pique 1986 II |
| v′ ≳ 60–82 | Nine states at the ²P₃/₂ + ²P₁/₂ limit | Steep rise in C, d and δ; eqQ trend reverses | Chen 2004; Pique 1986 |

- The suggested anomalies at v′ ≈ 42–47 are **not** real. Chen 2004's fit SDs there are normal, 6–33 kHz.
- How Knöckel 2004 handled perturbations: not available to me; that is in your extraction of the paper.
- Salami & Ross, J. Mol. Spectrosc. 233, 157 (2005) is "A molecular iodine atlas in ascii format". It contains no perturbation analysis that could be found; its abstract was not consulted.

---

#### References (all verified)

- Broyer, Vigué, Lehmann, J. Physique 39, 591 (1978). 10.1051/jphys:01978003906059100 — Crossref and full text (HAL jpa-00208791)
- Hackel, Casleton, Kukolich, Ezekiel, PRL 35, 568 (1975). 10.1103/PhysRevLett.35.568 — Crossref
- Levenson & Schawlow, PRA 6, 10 (1972). 10.1103/PhysRevA.6.10 — Crossref
- Yokozeki & Muenter, JCP 72, 3796 (1980). 10.1063/1.439594 — Crossref and abstract
- Vigué, Broyer, Lehmann, PRL 42, 883 (1979). 10.1103/PhysRevLett.42.883 — Crossref
- Bacis et al., JCP 73, 2641 (1980). 10.1063/1.440477 — Crossref
- Bordé et al., J. Physique 42, 1393 (1981). 10.1051/jphys:0198100420100139300 — OpenAlex abstract
- Pique et al., PRL 52, 267 (1984). 10.1103/PhysRevLett.52.267 — Crossref
- Pique et al., J. Physique 47, 1909 (1986). 10.1051/jphys:0198600470110190900 — Crossref and full text (HAL)
- Pique et al., J. Physique 47, 1917 (1986). 10.1051/jphys:0198600470110191700 — OpenAlex and full text (HAL)
- Pique et al., J. Physique 44, 347 (1983). 10.1051/jphys:01983004403034700 — Crossref and abstract
- Vigué, Broyer, Lehmann, J. Physique 42, 937 / 949 / 961 (1981). 10.1051/jphys:01981004207093700, 10.1051/jphys:01981004207094900, 10.1051/jphys:01981004207096100 — full texts (HAL)
- Vigué, Broyer, Lehmann, JCP 62, 4941 (1975). 10.1063/1.430409; Broyer, Vigué, Lehmann, JCP 64, 4793 (1976). 10.1063/1.432067; Vigué, Broyer, Lehmann, J. Phys. B 10, L379 (1977). 10.1088/0022-3700/10/10/004 — Crossref
- Tellinghuisen, JCP 57, 2397 (1972). 10.1063/1.1678600; Capelle & Broida, JCP 58, 4212 (1973). 10.1063/1.1678977; Paisner & Wallenstein, JCP 61, 4317 (1974). 10.1063/1.1681737 — Crossref and abstracts
- Martínez, Martínez, Castaño, J. Mol. Spectrosc. 128, 554 (1988). 10.1016/0022-2852(88)90170-1; Špirko & Blabla, J. Mol. Spectrosc. 129, 59 (1988). 10.1016/0022-2852(88)90258-5 — Crossref
- Saute & Aubert-Frécon, JCP 77, 5639 (1982). 10.1063/1.443770; Gerstenkorn, Luc, Amiot, J. Physique 46, 355 (1985). 10.1051/jphys:01985004603035500; Tromp, Le Roy, Gerstenkorn, Luc, J. Mol. Spectrosc. 100, 82 (1983). 10.1016/0022-2852(83)90027-9 — Crossref
- Jewsbury, Ridley, Lawley, Donovan, J. Mol. Spectrosc. 157, 33 (1993). 10.1006/jmsp.1993.1003 — Crossref
- Arie & Byer, JOSA B 10, 1990 (1993). 10.1364/JOSAB.10.001990; Razet & Picard, Metrologia 34, 181 (1997). 10.1088/0026-1394/34/2/10 — Crossref
- B. Bodermann, Dissertation, Universität Hannover (1998) — repo.uni-hannover.de bitstream 4031c4ef-bc3e-4a2d-b716-b5fd25d82e53; full text read; no DOI
- Bodermann, Knöckel, Tiemann, EPJ D 19, 31 (2002). 10.1140/epjd/e20020052 — Crossref; not read (paywalled)
- Knöckel, Bodermann, Tiemann, EPJ D 28, 199 (2004). 10.1140/epjd/e2003-00313-4 — Crossref
- Hong et al., JOSA B 18, 379 (2001). 10.1364/JOSAB.18.000379; JOSA B 18, 1416 (2001). 10.1364/JOSAB.18.001416; JOSA B 19, 946 (2002). 10.1364/JOSAB.19.000946 — abstracts
- Cheng et al., Opt. Lett. 27, 571 (2002). 10.1364/OL.27.000571; Chen & Ye, Chem. Phys. Lett. 381, 777 (2003). 10.1016/j.cplett.2003.10.052 — Crossref
- Chen, Cheng, Ye, JOSA B 21, 820 (2004). 10.1364/JOSAB.21.000820 — full text (JILA PDF)
- Chen, de Jong, Ye, JOSA B 22, 951 (2005). 10.1364/JOSAB.22.000951 — abstract
- Salumbides et al., Mol. Phys. 104, 2641 (2006). 10.1080/00268970600747696 — abstract; Salumbides et al., EPJ D (2008). 10.1140/epjd/e2008-00045-y — Semantic Scholar record; volume and pages not retrieved
- Sakagami et al., JOSA B 37, 1027 (2020). 10.1364/JOSAB.385779; Ikeda et al., JOSA B 39, 2264 (2022). 10.1364/JOSAB.465499 — Crossref
- Yoshiki et al., EPJ D 77 (2023). 10.1140/epjd/s10053-023-00712-7 — accepted-manuscript full text
- Matsunaga et al., Photonics 11, 770 (2024). 10.3390/photonics11080770; Nishiyama et al., JOSA B (2024). 10.1364/JOSAB.531115 — abstracts
- Baturo et al., JCP 144, 184310 (2016). 10.1063/1.4948630 — Crossref; Baturo et al., J. Phys. B (2018). 10.1088/1361-6455/aab6e3 — Semantic Scholar record
- Pashov et al., PRA 72, 062505 (2005). 10.1103/PhysRevA.72.062505; Salami & Ross, J. Mol. Spectrosc. 233, 157 (2005). 10.1016/j.jms.2005.06.002 — Crossref

**Web sources consulted:** [ResearchGate: Bodermann 2002](https://www.researchgate.net/publication/226561566_Widely_usable_interpolation_formulae_for_hyperfine_splittings_in_the_127I2_spectrum) · [Chen 2004, JILA PDF](https://jila.colorado.edu/~junye/yelabsOLD/pubs/scienceArticles/2004/sArticle_2004_04_Chen_JOSAB.pdf) · [Optica: Chen 2004 abstract](https://opg.optica.org/josab/abstract.cfm?uri=josab-21-4-820) · [Hannover repository: Bodermann thesis](https://repo.uni-hannover.de/server/api/core/bitstreams/4031c4ef-bc3e-4a2d-b716-b5fd25d82e53/content) · [Springer: Yoshiki 2023](https://link.springer.com/article/10.1140/epjd/s10053-023-00712-7) · [MDPI: Matsunaga 2024](https://www.mdpi.com/2304-6732/11/8/770) · [Tiemann publications](https://www.iqo.uni-hannover.de/en/tiemann/publications)

---

## Appendix E — Intensities and line shapes


The main results:
- A spectrum model cannot use a single Franck–Condon factor per band. Line strengths depend on J through the R-dependent transition moment, by anything from about 1% to about 70% at J≈150 depending on the band.
- Low-resolution cross sections depend on column density, so the package has to compute transmittance line by line and convolve that, not convolve the cross section.
- Natural widths depend on v′, J′ and F through predissociation. A 10³–10⁴ range is observed across the band: about 4 kHz near 508 nm, a few hundred kHz at 532 nm.
- Several broadening numbers below come only from search-engine summaries and are marked "unverified".

#### (a) Hönl–London factors and the line-strength formula

**Branches.** For 0u⁺ ← 0g⁺ (a parallel transition with ΔΩ=0) only P and R branches exist; Q is forbidden. The standard factors are:
- H_R = J″+1 (J′=J″+1)
- H_P = J″ (J′=J″−1)
- ΣH = 2J″+1

The factor (2−δ₀,Λ′+Λ″) equals 1 here, so the normalization is unambiguous for our case. Hansson & Watson (2005) addressed confusion in the literature over diatomic HLFs, mainly for perpendicular transitions. They derive a singlet–singlet formula for levels of definite parity (their Eq. 23) that is consistent with the standard sum rules. I have this from the abstract via a search summary, not the full text. Watson (2008) covers multiplet transitions, which are not needed for B–X.

**Recommended working formulas** (standard; the constants are derived here):
- Line intensity:
  S [cm molecule⁻¹] = 4.1624×10⁻¹⁹ · ν̃ · [g_ns(J″) e^(−c₂E″/T)/Q(T)] · (1−e^(−c₂ν̃/T)) · H_J′J″ · |⟨v′J′|μ(R)|v″J″⟩|² [D²]
- Einstein coefficient: A = 3.1362×10⁻⁷ ν̃³ |μ|² H/(2J′+1) s⁻¹, with ν̃ in cm⁻¹ and μ in D.
- Partition function: Q = Σ_v″J″ g_ns(J″)(2J″+1)e^(−c₂E/T). The nuclear-spin weights must be the same ones used in S.
- Cross section: σ(ν̃) = Σᵢ Sᵢ f(ν̃−ν̃ᵢ) + σ_cont(ν̃,T).
- The matrix element should be computed with J-dependent radial wavefunctions. This builds in the Herman–Wallis effect (see c).

#### (b) Nuclear spin statistics and hyperfine component counts

¹²⁷I (I=5/2) and ¹²⁹I (I=7/2) are fermions. In X ¹Σg⁺, even-J″ levels pair with the exchange-antisymmetric total spins I_T = 0,2,4(,6), and odd-J″ levels with the symmetric ones.

**Table 1: spin weights and number of strong (ΔF=ΔJ) components per rovibronic line at high J**

| Species | Allowed I_T (even J″ / odd J″) | g_ns even : odd | Components (even J″ / odd J″) | Source |
|---|---|---|---|---|
| ¹²⁷I₂ | 0,2,4 / 1,3,5 | 15 : 21 (=5:7) | 15 / 21 | See note¹ |
| ¹²⁹I₂ | 0,2,4,6 / 1,3,5,7 | 28 : 36 (=7:9) | 28 / 36 | **Derived here** by angular-momentum counting; not found stated in a source |
| ¹²⁷I¹²⁹I | I_T = 1…6, no exchange symmetry | 48 for every J | 48 | **Derived here**: 6×8 = 48 |

¹ Sources for ¹²⁷I₂:
- The PGOPHER I₂ example sets SymWt=5, AsymWt=7.
- Hong et al. (2004) measured "all the 21 hyperfine components" of R(85)33–0, which has odd J″.
- The BIPM label sets run a1–a15 for even J″ (e.g. R(56)32–0) and up to a21 for odd J″ (e.g. P(83)33–0).

At low J the count is Σ over I_T of min(2I_T+1, 2J+1), as derived here. For example, J″=1 of ¹²⁷I₂ has 9 components. Salumbides et al. (2006) measured and modeled the hyperfine structure of ¹²⁹I₂ and ¹²⁷I¹²⁹I (quadrupole and spin–rotation terms) and gave predictions for all isotopomers.

**Relative intensities in linear absorption.** This is the pure coupling-case result, derived by me from standard angular-momentum algebra (F = J + I_T, with I_T a spectator):

s(F′←F″) = (2F′+1)(2F″+1) {J′ F′ I_T; F″ J″ 1}² / (2I_T+1)

This sums to 1 over F′ and F″ for each I_T.
- The ΔF=ΔJ components are nearly equal, each carrying about (2F+1)/[(2I_T+1)(2J+1)] of the line.
- The ΔF≠ΔJ components scale as O(J⁻²).

Real eigenstates are mixtures. The quadrupole term mixes J with J±2 and different I_T of the same parity, so intensities must come from diagonalizing H_hf and transforming the dipole matrix.
- Wakasugi et al. (1989), 14–1 band, J=0–30: the "forbidden" ΔF=0 and ΔF=−ΔJ components are as strong as the allowed ones for J≲5, fall off rapidly, and nearly vanish above J≈25. This agrees with theory.
- Cheng et al. (2002) found that hyperfine patterns change strongly toward the dissociation limit.

**Saturation intensities differ from linear ones.** Vigué et al. (1981, part II) cite Ducasse & Couillaud: in the low-saturation limit, the saturated-absorption signal of a component scales as (2F+1)/(Γ_rad + Γ_pred(v,J,F)). Bordé & Bordé (1979) give the density-matrix treatment of hyperfine intensities in saturation spectroscopy.

#### (c) The transition moment μ(R), Herman–Wallis effects, and Franck–Condon factors

**Table 2: B–X transition-moment determinations**

| Work | Method | R range | Result (as stated) |
|---|---|---|---|
| Brewer & Tellinghuisen 1972 | Absorption, quantum yield | — | Qualitative μ(R) prediction, cited by Koffend et al. |
| Tellinghuisen 1973a | Absorption | — | Supports a linear dependence of dipole strength on r-centroid |
| Koffend, Bacis & Field 1979 | Gain of an optically pumped I₂ laser | r̄ 2.8–4.6 Å | Relative μ(R); agrees with the 1972 prediction |
| Tellinghuisen 1982 | Absorption, reassessed | 2.6–2.8 Å | \|μe\|² stays close to the 1973 values |
| Lamrini et al. 1994 | LIF Fourier-transform intensities | r̄ 2.633–6.035 Å | \|μe\|² falls by about 4 orders of magnitude, "perhaps the largest variation yet observed" |
| Tellinghuisen 1997 | Inversion of radiative rates vs v via a sum rule | — | Peak 1.85±0.1 D² near R=3.3 Å; prior studies' quantitative consistency is "poor" |
| Tellinghuisen 2011 | Least-squares simulation of 520–640 nm absorption at 0.1 nm resolution | — | \|μe\|² with <2% relative standard error over most of the region |
| Pique et al. 1983 | Lifetimes of hyperfine sublevels, v′=43 | — | Γ_rad = (0.314±0.018)×10⁶ s⁻¹, i.e. τ_rad = 3.18 μs |

Functional forms and coefficients are not in the abstracts. Take them from Tellinghuisen 2011 and 1997, and Lamrini 1994.

**Herman–Wallis size: a Morse-potential illustration, order of magnitude only.** Morse X and B potentials with standard constants, sinc-DVR, R-branch J′=J″+1:

| Band | FCF at J″=0 | J″=56 | J″=100 | J″=150 | Change 0→150 |
|---|---|---|---|---|---|
| 32–0 (532 nm) | 3.15e-2 | 3.10e-2 | 3.01e-2 | 2.81e-2 | −11% |
| 43–0 (515 nm) | 2.20e-2 | 2.20e-2 | 2.20e-2 | 2.18e-2 | −1% |
| 11–5 (633 nm) | 7.7e-3 | 8.5e-3 | 1.01e-2 | 1.33e-2 | +73% |

- The r-centroids also shift by about 0.01 Å between J=0 and 150.
- The average slope of ln|μ|² implied by Lamrini (ln 10⁴ over 3.4 Å ≈ 2.7 Å⁻¹) suggests a further change of about 3% (an estimate).
- Conclusion: compute the full ⟨v′J′|μ(R)|v″J″⟩. J=0 FCF×HLF is not good enough.

**FCF sources:** Steinfeld et al. 1965; Tellinghuisen 1978 (JQSRT "Intensity factors for the I₂ B↔X band system").

#### (d) Continuum absorption and absolute cross sections

- **Decomposition (Tellinghuisen 1973a):**
  - The ¹Πu(1u) ← X continuum peaks at 498.5 nm with ε = 200 L mol⁻¹ cm⁻¹. That is 7.6×10⁻¹⁹ cm² by the conversion σ = 3.8235×10⁻²¹ ε. The ¹Πu curve behaves as r⁻⁹.
  - A(1u) ← X peaks at 673.0 nm with ε = 41, i.e. 1.6×10⁻¹⁹ cm².
- **Revisions of those strengths:**
  - 1982: both the ¹Πu and A strengths are about 10% lower.
  - 2011: the C(¹Πu) ← X strength is lowered by a further 25%, and the A estimates are supported.
- **Continuum below D₀ (Tellinghuisen 1973b):** absorption from thermally excited v″ carries the B←X continuum well into the discrete region below D₀ (the λ<499.5 nm limit for v″=0).

**Table 3: absolute cross-section data**

| Dataset | T | Resolution | Key value |
|---|---|---|---|
| Tellinghuisen 1973 | Room temperature | 2.3–2.9 nm slit | σ(500 nm) = (2.20±0.07)×10⁻¹⁸ cm² |
| Saiz-Lopez et al. 2004 | 295 K | 4 cm⁻¹ (0.1 nm), FTS, 182–750 nm | Peak (4.24±0.50)×10⁻¹⁸ cm² at 533.0 nm; 2.29(27)×10⁻¹⁸ cm² at 500 nm |
| Bauer et al. 2004 (tabulated by Spietz; not checked independently) | — | — | σ(500) = 2.25±0.09 |
| Spietz et al. 2006 | ~298 K | Grating 0.25 nm, plus FTS | σ(500) = (2.186±0.021)×10⁻¹⁸ cm², independent of vapour pressure; weighted mean 2.191±0.02 |
| Tellinghuisen 2011 | Room temperature | 0.1 nm, 520–640 nm | Discrete and continuous parts separated by simulation |

- **Continuum fraction (an estimate):** at 500 nm the ¹Πu share is roughly 24–35% (0.52–0.76 of 2.19×10⁻¹⁸ cm²). The rest is B←X continuum from hot v″ plus converging bands. The continuum fraction at the 533 nm peak is **unknown** from what I accessed; Tellinghuisen 2011 has it.
- **Saturation warning (Spietz 2006):** at ≥1 nm FWHM, apparent cross sections depend on the column density of the reference spectrum. The error is about 13% for atmospheric conditions and up to 45% at low pressure.
- **Bound–free tools:**
  - BCONT (Le Roy 1989; BCONT 2.2, report CP-650R2, 2004) includes fitting of μ(R) and V(R).
  - The RKR-like inversion of Child, Essén & Le Roy (1983).
  - Tellinghuisen's 1985 review of the Franck–Condon principle in bound–free transitions.

#### (e) Doppler width

FWHM Δν_D = (ν₀/c)·√(8kT ln2/m) = 7.162×10⁻⁷ ν₀ √(T/M), with M = 253.809 u. Computation:

| λ | 293.15 K | 300 K |
|---|---|---|
| 500 nm | 461.5 MHz | 466.9 MHz (0.0156 cm⁻¹) |
| 532 nm | 433.8 MHz | 438.8 MHz (0.0146 cm⁻¹) |
| 633 nm | 364.6 MHz | 368.8 MHz (0.0123 cm⁻¹) |
| 800 nm | 288.5 MHz | 291.8 MHz (0.0097 cm⁻¹) |

- The width scales as √T, i.e. +0.17% per K at 300 K.
- A line's hyperfine spread is several hundred MHz; Wolf used ΔeQq ≈ 1950 MHz at 675 nm. That is comparable to the Doppler width, so a Doppler-limited line is a sum of 15 or 21 Gaussians, which is how Wolf modeled it, not a single Voigt.

#### (f) Pressure broadening, shift and quenching

**Table 5: broadening and shift coefficients** (FWHM; 1 MHz/Torr = 7.50 kHz/Pa)

| Line | Perturber | Broadening | Shift | Source and status |
|---|---|---|---|---|
| R(56)32–0 a10, 532 nm | I₂ | 63 kHz/Pa | — | Fang, Wang & Shy 2006; number from a search summary, abstract not accessed |
| R(56)32–0 a1, 532 nm | I₂ | 74 kHz/Pa | −1.3 kHz/Pa | Attributed to Ye et al. 1999 by search summaries; **unverified** |
| R(56)32–0 a10 | I₂ | — | −4.2 kHz/Pa | Nevsky et al. 2001, as adopted in the BIPM mise en pratique (MeP) for 532 nm; verified |
| P(46)44–0 a10 (514.7 nm) and R(56)32–0 a10 | I₂ | 3f peak-to-peak width 280→370 and 330→380 kHz over 0.49→2.48 Pa, i.e. about 45 and 25 kHz/Pa (slopes derived here; not Lorentzian FWHM) | — | Hrabina et al. 2014, Table 3; verified |
| R(127)11–5, 633 nm | I₂ | — | −15 kHz/°C of cold-finger temperature at 15 °C, which I convert to about −9.6 kHz/Pa² | BIPM MeP 633; verified |
| (17,1) P(10), P(70) | Ar | 10.7±0.4 and 8.3±0.3 MHz/Torr (80, 62 kHz/Pa); no variation across the 15 components | — | Phillips & Perram 2008; **from a search summary** |
| Doppler-limited lines near 543 nm | Air | 9.98(17) MHz/Torr | −1.091(75) MHz/Torr | Fletcher & McDaniel 1995, as quoted by Wolf 2009 |
| 675 nm, J″=69 / 122 | Air | 9.05(6) / 8.69(12) MHz/Torr | −1.17(1) / −1.26(3) MHz/Torr | Wolf 2009 thesis (12 gases, 292–388 K); verified |
| 675 nm, J″=69 / 122 | H₂O | 11.6 / 12.0 MHz/Torr | −1.12 / −0.96 MHz/Torr | Wolf 2009; verified |

² Conversion: the BIPM pressures P(−15 °C) = 0.83 Pa and P(−5 °C) = 2.46 Pa give ΔH/R = 7522 K. That implies P(15 °C) ≈ 17.2 Pa and dP/dT ≈ 1.56 Pa/K, so −15 kHz/K ÷ 1.56 Pa/K ≈ −9.6 kHz/Pa.

**Quenching.**
- Capelle & Broida (1973): self-quenching cross sections of 47–90 Å² for v′≈5–70; He and H₂ ≲1 Å²; Xe up to 23 Å².
- Estimate: σ = 70 Å² at 300 K gives n·σ·v_rel ≈ 3.8×10⁴ s⁻¹ Pa⁻¹, i.e. about 6 kHz/Pa FWHM. Quenching is therefore only about 10% of the 63–74 kHz/Pa self-broadening; the rest is elastic or rotationally inelastic dephasing.
- The v′ dependence of self-broadening is **unknown**. Wolf's data show a weak J dependence.

#### (g) Natural linewidth

**Mechanism.** Γ = Γ_rad + Γ_pred(v,J,F) (Vigué, Broyer & Lehmann 1981). Predissociation goes through ¹Πu(1u) and has three terms:
- gyroscopic, Cv² J(J+1);
- hyperfine, ∝ av², which depends on the component;
- an interference term, ∝ av·Cv.

**Table 6: lifetimes and widths**

| Level or region | Data | Source |
|---|---|---|
| v′≈5–70 | τ ranges from <0.4 μs to >7 μs, strongly v′-dependent | Capelle & Broida 1973 |
| All v′ | 1/τ is linear in J′(J′+1) | Broyer, Vigué & Lehmann 1975 |
| All v′ | Cv² peaks at about 1700 s⁻¹ near v≈5 and 450 s⁻¹ near v≈25; Γ_rad falls with v; natural width including predissociation is "almost always less than 1 MHz"; sublevels F=J±5 decay very differently | Vigué et al. 1981, part II |
| v′=43 | Γ_rad = 0.314×10⁶ s⁻¹; Cv = 5.85±0.65 s^-½; av = −142±5 s^-½ | Pique et al. 1983 |
| 514.5 nm, v′=43, J′≈12 | Estimate: Γ ≈ 3.2–3.4×10⁵ s⁻¹, i.e. Γ/2π ≈ 51–54 kHz | Derived from Pique 1983 |
| 532 nm, v′=32 | Observed 3f width 330 kHz at 0.49 Pa (P(46)44–0 at 514.7 nm: 280 kHz). A search summary says natural widths go "down to ~100 kHz"; **unverified** | Hrabina 2014 |
| 523–498 nm | Narrowest observed width ≈4 kHz, near 508 nm. Inference: τ of order tens of μs | Cheng et al. 2002 |
| 501.7 nm, R(26)62–0 | Observed width 45 kHz (HWHM) | Goncharov et al. 2007 |
| 633 nm, v′=11, J′=128 | **Unknown**: Shotton & Chapman 1972 and Titov et al. 1997 were not accessed. If Cv² is 10²–10³ s⁻¹, the gyroscopic term alone would be 0.26–2.6 MHz (a scenario) | — |

Conversion used: Γ/2π = 1/(2πτ), so τ = 0.4, 1, 2 and 7 μs correspond to 398, 159, 80 and 23 kHz.

#### (h) Sub-Doppler spectroscopy: saturation, crossovers, recoil, transit time

**Operating conditions.**
- BIPM MeP 532: pump intensity (17±11) mW cm⁻², 1 MHz peak-to-peak FM with 3f detection, cold finger at −15 °C (0.83 Pa).
- BIPM MeP 633: 6 MHz peak-to-peak FM, cold finger at 15 °C, 10 mW intracavity power.

**Power broadening.** Hrabina 2014, Table 2: raising the pump from 2.3 to 23 mW cm⁻² widens the lines only from 280 to 320 kHz (514.7 nm) and 330 to 380 kHz (532 nm).

**Saturation intensity.**
- Titov et al. (1997) determined the 633 nm saturation parameter, but I have **not** seen the value.
- An order-of-magnitude estimate for strong 532 nm components is 10–100 mW cm⁻². Assumptions: A_line ≈ 10⁴ s⁻¹, Γ ≈ 10⁶ s⁻¹, σ₀ = (λ²/2π)(A_line/Γ), and lower-level relaxation limited by transit time. This fits the weak power broadening above.

**Detection schemes.**
- FM heterodyne: Hall et al. 1981.
- Modulation transfer spectroscopy (MTS): Shirley 1982; Camy, Bordé & Ducloy 1982. Used for I₂ at 532 nm by Eickhoff & Hall 1995 and Ye et al. 1999.
- Wave-front curvature shifts: Hall & Bordé 1976.

**Crossover resonances** appear midway between transitions that share a level, and the I₂ hyperfine manifolds are full of them. Bordé, Camy & Decomps (1979) measured the recoil shift at 514.5 nm to within 12% of theory. They used the fact that crossovers with a common upper level and those with a common lower level have opposite recoil shifts.

**Recoil doublet** (computed here, splitting h/(mλ²)):

| λ | Splitting | Half-splitting |
|---|---|---|
| 500 nm | 6.29 kHz | 3.14 kHz |
| 515 nm | 5.93 kHz | 2.96 kHz |
| 532 nm | 5.55 kHz | 2.78 kHz |
| 633 nm | 3.92 kHz | 1.96 kHz |

The doublet cannot be resolved at 532 nm (widths ~300 kHz). It is comparable to the ~4 kHz widths near 508 nm, and it is a kHz-level systematic shift everywhere.

**Transit-time broadening** (computed here): FWHM ≈ (√(2 ln2)/π)(v/w) = 0.375 v/w. With v_p = 140 m/s at 300 K this gives 105, 53, 26 and 10.5 kHz for beam radii w = 0.5, 1, 2 and 5 mm.

**What the package should output.**
1. A hyperfine-resolved line list with ν, its uncertainty, the assignment including F′/F″ and the BIPM component label, S(T), A, a Γ_nat model, and pressure coefficients with verified/unverified flags.
2. A Doppler-limited σ(ν; T, P): a sum of Voigt profiles for each component plus the continuum. For comparison with a spectrometer, convolve the *transmittance* with the instrument function.
3. Optionally, a saturation/MTS/3f simulator including crossovers and recoil, with the saturation parameters left user-adjustable.

#### References

| # | Reference | DOI | How verified |
|---|---|---|---|
| 1 | A. Hansson, J. K. G. Watson, J. Mol. Spectrosc. 233, 169 (2005) | 10.1016/j.jms.2005.06.009 | Crossref |
| 2 | J. K. G. Watson, J. Mol. Spectrosc. 252, 5 (2008) | 10.1016/j.jms.2008.04.014 | Crossref |
| 3 | J. Tellinghuisen, J. Chem. Phys. 58, 2821 (1973) | 10.1063/1.1679584 | OpenAlex abstract |
| 4 | J. Tellinghuisen, J. Chem. Phys. 59, 849 (1973) | 10.1063/1.1680103 | OpenAlex abstract |
| 5 | J. Tellinghuisen, J. Chem. Phys. 76, 4736 (1982) | 10.1063/1.442791 | OpenAlex abstract |
| 6 | J. Tellinghuisen, J. Chem. Phys. 106, 1305 (1997) | 10.1063/1.473971 | OpenAlex abstract |
| 7 | J. Tellinghuisen, J. Chem. Phys. 134 (2011) | 10.1063/1.3555623 | OpenAlex abstract |
| 8 | J. Tellinghuisen, JQSRT 19, 149 (1978) | 10.1016/0022-4073(78)90074-2 | Crossref |
| 9 | J. Tellinghuisen, Adv. Chem. Phys., pp. 299–369 (1985) | 10.1002/9780470142844.ch7 | Crossref |
| 10 | L. Brewer, J. Tellinghuisen, J. Chem. Phys. 56, 3929 (1972) | 10.1063/1.1677797 | Crossref |
| 11 | J. B. Koffend, R. Bacis, R. W. Field, J. Chem. Phys. 70, 2366 (1979) | 10.1063/1.437744 | Crossref, OpenAlex |
| 12 | M. Lamrini et al., J. Chem. Phys. 100, 8780 (1994) | 10.1063/1.466732 | Crossref, OpenAlex |
| 13 | J. I. Steinfeld et al., J. Chem. Phys. 42, 25 (1965) | 10.1063/1.1695685 | Crossref |
| 14 | R. H. Tipping, J.-P. Bouanich, JQSRT 71, 99 (2001) | 10.1016/s0022-4073(01)00014-0 | Crossref |
| 15 | A. Saiz-Lopez et al., Atmos. Chem. Phys. 4, 1443 (2004) | 10.5194/acp-4-1443-2004 | Crossref; full text read |
| 16 | P. Spietz, J. C. Gómez Martín, J. P. Burrows, Atmos. Chem. Phys. 6, 2177 (2006) | 10.5194/acp-6-2177-2006 | Publisher PDF read; DOI follows the journal's pattern but was not resolved |
| 17 | R. J. Le Roy, Comput. Phys. Commun. 52, 383 (1989) | 10.1016/0010-4655(89)90113-6 | Crossref |
| 18 | R. J. Le Roy, G. T. Kraemer, BCONT 2.2, Univ. Waterloo report CP-650R2 (2004) | none | uwaterloo.ca page |
| 19 | M. S. Child, H. Essén, R. J. Le Roy, J. Chem. Phys. 78, 6732 (1983) | 10.1063/1.444673 | Crossref |
| 20 | G. A. Capelle, H. P. Broida, J. Chem. Phys. 58, 4212 (1973) | 10.1063/1.1678977 | OpenAlex abstract |
| 21 | J. A. Paisner, R. Wallenstein, J. Chem. Phys. 61, 4317 (1974) | 10.1063/1.1681737 | OpenAlex |
| 22 | M. Broyer, J. Vigué, J. C. Lehmann, J. Chem. Phys. 63, 5428 (1975) | 10.1063/1.431275 | OpenAlex |
| 23 | J. Vigué, M. Broyer, J. C. Lehmann, J. Physique 42, 937 / 949 / 961 (1981) | 10.1051/jphys:01981004207093700; 10.1051/jphys:01981004207094900; 10.1051/jphys:01981004207096100 | Crossref; HAL full text of part II |
| 24 | J. P. Pique et al., J. Physique 44, 347 (1983) | 10.1051/jphys:01983004403034700 | OpenAlex abstract |
| 25 | K. C. Shotton, G. D. Chapman, J. Chem. Phys. 56, 1012 (1972) | 10.1063/1.1677205 | Crossref |
| 26 | W.-Y. Cheng et al., Opt. Lett. 27, 571 (2002), plus erratum | 10.1364/ol.27.000571; erratum 10.1364/ol.27.001076 | Crossref, OpenAlex |
| 27 | A. Goncharov et al., Metrologia 44, 275 (2007) | 10.1088/0026-1394/44/5/003 | Crossref; arXiv:0706.4005 |
| 28 | J. Hrabina et al., Meas. Sci. Rev. 14, 213 (2014) | 10.2478/msr-2014-0029 | Crossref; full text read |
| 29 | H.-M. Fang, S. C. Wang, J.-T. Shy, Opt. Commun. 257, 76 (2006) | 10.1016/j.optcom.2005.07.016 | Crossref |
| 30 | J. Ye, L. Robertsson, S. Picard, L.-S. Ma, J. L. Hall, IEEE Trans. Instrum. Meas. 48, 544 (1999) | 10.1109/19.769654 | Crossref, OpenAlex |
| 31 | A. Yu. Nevsky et al., Opt. Commun. 192, 263 (2001) | 10.1016/s0030-4018(01)01190-7 | Crossref |
| 32 | M. L. Eickhoff, J. L. Hall, IEEE Trans. Instrum. Meas. 44, 155 (1995) | 10.1109/19.377797 | Crossref |
| 33 | BIPM mise en pratique, I₂ 532 nm (2007; updated 2012) and 633 nm (2003) | none | PDFs read at bipm.org |
| 34 | D. G. Fletcher, J. C. McDaniel, JQSRT 54, 837 (1995) | 10.1016/0022-4073(95)00105-t | Crossref |
| 35 | E. N. Wolf, PhD thesis (2009), arXiv:0910.5053 | 10.48550/arXiv.0910.5053 | Full text read |
| 36 | G. T. Phillips, G. P. Perram, JQSRT 109, 1875 (2008) | 10.1016/j.jqsrt.2007.12.011 | Crossref |
| 37 | A. Titov, I. Malinovsky, M. Erin, Opt. Commun. 136, 327 (1997) | 10.1016/s0030-4018(96)00685-2 | Crossref |
| 38 | J. Bordé, Ch. J. Bordé, J. Mol. Spectrosc. 78, 353 (1979) | 10.1016/0022-2852(79)90063-8 | Crossref |
| 39 | Ch. J. Bordé, G. Camy, B. Decomps, Phys. Rev. A 20, 254 (1979) | 10.1103/physreva.20.254 | OpenAlex abstract |
| 40 | Ch. J. Bordé et al., J. Physique 42, 1393 (1981) | 10.1051/jphys:0198100420100139300 | Crossref |
| 41 | J. H. Shirley, Opt. Lett. 7, 537 (1982) | 10.1364/ol.7.000537 | Crossref |
| 42 | G. Camy, Ch. J. Bordé, M. Ducloy, Opt. Commun. 41, 325 (1982) | 10.1016/0030-4018(82)90406-0 | Crossref |
| 43 | J. L. Hall, L. Hollberg, T. Baer, H. G. Robinson, Appl. Phys. Lett. 39, 680 (1981) | 10.1063/1.92867 | Crossref |
| 44 | J. L. Hall, C. J. Bordé, Appl. Phys. Lett. 29, 788 (1976) | 10.1063/1.88949 | Crossref |
| 45 | M. Wakasugi et al., J. Opt. Soc. Am. B 6, 1660 (1989) | 10.1364/josab.6.001660 | OpenAlex abstract |
| 46 | E. J. Salumbides et al., Mol. Phys. 104 (16–17) (2006) | 10.1080/00268970600747696 | OpenAlex abstract |
| 47 | F.-L. Hong et al., J. Opt. Soc. Am. B 21, 88 (2004) | 10.1364/josab.21.000088 | Crossref; full text via NIST |
| 48 | PGOPHER help page "Hyperfine structure in the B-X transition in I₂" | none | Page fetched |
