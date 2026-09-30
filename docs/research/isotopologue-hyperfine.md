# Hyperfine structure of ¹²⁹I₂ and ¹²⁷I¹²⁹I: what Salumbides et al. 2006 gives us

**Sources**

* **S06**: E. J. Salumbides, K. S. E. Eikema, W. Ubachs, U. Hollenstein, H. Knöckel and E. Tiemann, "The hyperfine structure of ¹²⁹I₂ and ¹²⁷I¹²⁹I in the B³Π₀u⁺–X¹Σg⁺ band system", *Mol. Phys.* **104**, 2641–2652 (2006), doi:10.1080/00268970600747696. Local copy: `papers/Salumbides2006_MolPhys104_2641.pdf` (VU repository copy; PDF pp. 3–14 = journal pp. 2641–2652).
* **S08**: Salumbides et al., *Eur. Phys. J. D* **47**, 171 (2008). It has no new hyperfine information (§7).
* **BKT02**: Bodermann, Knöckel & Tiemann, *Eur. Phys. J. D* **19**, 31 (2002). This is what `src/i2spec/hfs_params.py` implements now.

**Transcription.** Equations and Tables 1–4 were read from the page images. `pdftotext` drops Greek letters, minus signs and exponents in this PDF, so its output was used only to cross-check digits. Items marked **[flag]** are ambiguities or misprints in the source. Items marked **[inference]** are reasoning made here, not statements in the paper.

---

## 0. Summary

1. **The Hamiltonian is the same four-term effective Hamiltonian as BKT02.** S06 does not write out the operators; it names the terms and defers to Broyer et al. (1978) for the homonuclear case and to Freed (1966) for the heteronuclear case. Its ¹²⁷I₂ eqQ_X is −2452.285 MHz, the same sign and size as BKT02, so our sign and normalization conventions carry over unchanged.
2. **For ¹²⁷I¹²⁹I, each nucleus has its own eqQ and its own spin–rotation constant C:** eqQ⁽¹²⁷⁾, eqQ⁽¹²⁹⁾, C⁽¹²⁷⁾, C⁽¹²⁹⁾ in each electronic state. The nuclear spin–spin terms (δ, d) are left out entirely. There is no exchange symmetry. The paper's basis is |Ωv(JI₁)F₁I₂FM_F⟩.
3. **Parameters for any isotopologue are the ¹²⁷I₂ formulae times a nuclear-moment ratio** (eq. 8): γ_Q = 0.701213 for eqQ of a ¹²⁹I nucleus, γ_μ = 0.6655 for C, and γ_μ² for δ and d in ¹²⁹I₂. The χ_lk get no reduced-mass scaling. The energies E_v,J in the denominators are **¹²⁷I₂ energies at the same (v, J)** for every isotopologue.
4. **S06 also refits the ¹²⁷I₂ formulae**, in Tables 2–4, for eqQ_X, eqQ_B, C_B, δ_B and d_B, extending validity from v′ ≤ 43 to v′ ≤ 53. C_X, δ_X and d_X are unchanged from BKT02. Our checks (§8) found three things:
   * **[flag] The printed δ_B and d_B formula has a sign error in the Ω = 0 perturber term.** As printed, δ_B(v′ = 32) = −78 kHz, against −6 kHz from BKT02. Flipping that term in both δ and d gives −6.8 and −41.3 kHz.
   * **E_v,J must be the J = 0 energy**, as in BKT02. Using E at J, R(87)33-0 misses BIPM by about 4 MHz rms.
   * **For ¹²⁷I₂ at 532 nm, S06 eqQ is as good as BKT02, but S06 C_B is worse.** The BIPM rms goes from 24 to 57 kHz for R(56)32-0 and from 89 to 266 kHz for R(87)33-0. Keep BKT02 for ¹²⁷I₂ at v′ ≤ 43.
5. **Neither S06 nor S08 tabulates any ¹²⁹I₂ or ¹²⁷I¹²⁹I component positions or per-line parameters.** The per-line ΔeqQ and ΔC values are in the S06 electronic supplement [27], which we do not have. The only targets inside S06 are Figs. 3–5, whose vector traces we digitized (§4.1b). The high-precision isotopologue splittings S06 relies on are in Quinn, *Metrologia* **40**, 103 (2003) (the CIPM 2001 list, 633 nm).
6. **Our implementation of the S06 recipe reproduces S06's own calculated spectra.**
   * Fig. 4 (¹²⁹I₂ R34(11–3)) and Fig. 5 (¹²⁷I¹²⁹I R37(11–3)): every resolved peak of the calculated profile agrees within 0.3 MHz (rms 0.1 MHz).
   * Fig. 3 control (¹²⁷I₂ with BKT02): agrees to 0.09 MHz rms, which validates the digitization.
   * The match needs per-nucleus spin–rotation. With one C for both nuclei, the Fig. 5 peaks miss by up to 1.8 MHz (0.9 MHz rms).

---

## 1. The Hamiltonian

### 1.1 What the paper writes (§3.1, p. 2643)

> H_hfs,eff = H_NEQ + H_SR + H_SS + H_TS. (1)
>
> "H_NEQ, H_SR, H_SS and H_TS represent the (effective) nuclear electric quadrupole, the nuclear spin–rotation, the scalar nuclear spin–spin and the tensorial nuclear spin–spin interactions. The matrix elements of each of these terms can be separated into a product of a geometrical factor g_i and a hyperfine parameter, i.e. eqQ, C, δ or d:"
>
> ⟨α′, F | H_hfs,eff | α, F⟩ = eqQ·g_eqQ + C·g_SR + δ·g_SS + d·g_TS. (2)
>
> "The g_i are functions of appropriate angular momentum quantum numbers of the system, and can be calculated with spherical tensor algebra."

No explicit g_i are given anywhere in S06.

**Homonuclear ¹²⁷I₂ and ¹²⁹I₂:**

> "The homonuclear molecules ¹²⁷I₂ and ¹²⁹I₂ are characterized by a wave function |Ωv(I₁I₂)IJFM_F⟩ [17]. The two nuclear spins I₁ and I₂ couple to a total nuclear spin I, which is then coupled with the rotational angular momentum J to give a total angular momentum F. … Here the computer code, described in more detail in [18], was applied."

[17] is M. Broyer, J. Vigué and J. C. Lehmann, *J. Phys. (France)* **39**, 591 (1978). [18] is H. Knöckel, S. Kremser, B. Bodermann and E. Tiemann, *Z. Phys. D* **37**, 43 (1996).

**Heteronuclear ¹²⁷I¹²⁹I** (quoted in full, because it is the whole specification):

> "For the spectra of the heteronuclear species a coupling scheme corresponding to the notation |Ωv(JI₁)F₁I₂FM_F⟩ is used, and all simplifications due to symmetry like in the homonuclear case are absent. Moreover, for each nucleus appropriate hyperfine parameters must be used, which means eqQ⁽¹²⁷⁾, eqQ⁽¹²⁹⁾, C⁽¹²⁷⁾, C⁽¹²⁹⁾, … values for each electronic state. An appropriate computer code was available. The corresponding matrix elements are taken from [19]. The nuclear spin–spin interactions are not included. Due to its small contribution of not more than 400 kHz to the hyperfine splitting and in view of the limited signal-to-noise ratio in our spectra this neglect does not change the results on the molecular parameters." (pp. 2643–2644)

[19] is K. F. Freed, *J. Chem. Phys.* **45**, 1714 (1966).

**Diagonalization and intensities** (p. 2644):

> "…the computer codes calculate the hyperfine splitting frequencies and also relative intensities of the hyperfine lines, using for the dipole matrix element the eigenvectors determined during the diagonalization of the hyperfine interaction matrix for each state, extending the rotational space up to ΔJ = 0, ±2."
>
> "It was verified that the code for heteronuclear species gives the same results, within numerical accuracy, as the code for the homonuclear species for situations where the results should coincide."

**Nuclear spins and component counts** (p. 2646):

> "Due to the nuclear spins I₁ = 5/2 for ¹²⁷I and I₁ = 7/2 for ¹²⁹I, the full hyperfine pattern of a transition with selection rule ΔF = ΔJ of ¹²⁷I₂ exhibits 15 hyperfine components when J″ is even and 21 when J″ is odd; for ¹²⁹I₂ we have 28 resp. 36 components. Transitions of the mixed isotopomer ¹²⁷I¹²⁹I always have 48 hyperfine components, independent of J″ even or odd. The hyperfine components in the spectra are assigned by labels a₁ to a_n from lower to higher frequencies…"

Check against our `allowed_spins`. For ¹²⁹I₂ with the X-state "g" rule (I ≡ J mod 2), even J″ gives I ∈ {0, 2, 4, 6}, so Σ(2I+1) = 1 + 5 + 9 + 13 = 28. Odd J″ gives I ∈ {1, 3, 5, 7}, so Σ(2I+1) = 3 + 7 + 11 + 15 = 36. For ¹²⁷I¹²⁹I, I = 1, …, 6 with no restriction gives 3 + 5 + 7 + 9 + 11 + 13 = 48 (for J ≥ 6). Our homonuclear parity rule works unchanged for ¹²⁹I₂, which is also a pair of identical fermions. ¹²⁷I¹²⁹I needs `symmetry=None`.

### 1.2 Conventions compared with ours

| Item | S06 | i2spec now | Status |
|---|---|---|---|
| Terms | eqQ, C, δ, d (eq. 2) | eqQ·H_EQ + C·(I·J) + d·H_TSS + δ·(I₁·I₂) | Same set |
| eqQ sign | eqQ_X(¹²⁷I₂) = −2452.285 MHz (Table 4); BKT02 has −2452.2916 | Townes–Schawlow form with eqQ_X ≈ −2452 MHz, matched to BIPM 532 nm | **Consistent.** S06 is a refit in the BKT02 framework, and that framework reproduces BIPM with our operator |
| Per-nucleus quadrupole normalization | Not written. Eq. (8) scales eqQ by the bare ratio Q⁽ⁱ⁾/Q⁽ʳᵉᶠ⁾ even though I changes from 5/2 to 7/2 | Σₙ …/[2iₙ(2iₙ−1)(2J−1)(2J+3)] | **Consistent [inference].** The Q-ratio scaling only makes sense if the spin dependence is in g_eqQ, i.e. the standard 1/[2i(2i−1)] per nucleus, which is what we have |
| Spin–rotation, homonuclear | One C | One C acting on I = I₁ + I₂ | Same |
| Spin–rotation, heteronuclear | Separate C⁽¹²⁷⁾ and C⁽¹²⁹⁾ | One C only | **Difference. Code change needed** (§6). The operator must be C⁽¹²⁷⁾ I₁₂₇·J + C⁽¹²⁹⁾ I₁₂₉·J [inference: the natural reading of "C⁽¹²⁷⁾, C⁽¹²⁹⁾ values for each electronic state"; the paper does not write the operator] |
| Magnetic scaling | γ_μ = (μ/I)⁽ⁱ⁾ / (μ/I)⁽ʳᵉᶠ⁾, i.e. a g-factor ratio | — | C ∝ g_I, which is consistent with H_SR = C I·J per nucleus |
| Scalar spin–spin | δ; scaled by γ_μ² in ¹²⁹I₂; **omitted for ¹²⁷I¹²⁹I** | δ I₁·I₂ | Same operator [inference; not written]. For the heteronuclear case set δ = 0 to reproduce S06 |
| Tensor spin–spin | d; scaled by γ_μ² in ¹²⁹I₂; **omitted for ¹²⁷I¹²⁹I** | Our H_TSS (no factor 5), fixed against BIPM | Same normalization [inference: the δ_B, d_B values are a refit in the BKT02 framework]. Set d = 0 for the heteronuclear case |
| Rotational mixing | ΔJ = 0, ±2 | ΔJ = 0, ±2 | Same |
| Homonuclear basis | \|(I₁I₂)I J F⟩ | \|J,(i₁i₂)I; F⟩ | Same |
| Heteronuclear basis | \|(JI₁)F₁ I₂ F⟩, no symmetry restriction | Coupled \|J,(i₁i₂)I; F⟩ | **Equivalent** if all I = \|i₁−i₂\|, …, i₁+i₂ (1–6) are kept for every J, and ΔI = ±1 couplings are kept. Eigenvalues are basis-independent. Which nucleus is "I₁" in the heteronuclear scheme is not stated; it does not affect eigenvalues |
| Energies in denominators | E_v,J "calculated from the Dunham parameters in [23] for ¹²⁷I₂" | E(v, J = 0) from our Hannover 2008 model | **Resolved numerically: use J = 0**, as we do now (§2.3, §8c) |

`hfs_params.py` already departs from the printed BKT02 text by using BKT02 eqs. (12)–(13) as δ_B and d_B, not Δδ and Δd. S06 confirms that reading explicitly (p. 2649): *"remark: erroneously the formulas (12) and (13) in [4] read Δδ and Δd, they should, however read δ_B and d_B"*.

---

## 2. How the parameters are obtained

### 2.1 Isotope scaling (eq. 8, p. 2645)

> χ⁽ⁱ⁾(v, J) = χ⁽ʳᵉᶠ⁾ · γ⁽ⁱ⁾_Q,μ, (8)
>
> "where χ⁽ʳᵉᶠ⁾ is the corresponding parameter for the reference isotopomer ¹²⁷I₂ and γ_Q^i = Q⁽ⁱ⁾/Q⁽ʳᵉᶠ⁾ and γ_μ^i = (μ⁽ⁱ⁾/I₁⁽ⁱ⁾)/(μ⁽ʳᵉᶠ⁾/I₁⁽ʳᵉᶠ⁾) are the ratios of the nuclear electric quadrupole moments and of the nuclear dipole moments referred to the nuclear spin quantum numbers I₁. … In this way the numbers given later can be used directly for ¹²⁷I₂ with γ = 1, while for the other isotopomers the factors γ are different from one. These modifications are applied to the nuclear electric quadrupole and the nuclear spin–rotation interaction, where the nuclear moments enter linearly. The nuclear spin–spin interaction is a product of two matrix elements containing a nuclear magnetic moment. So here the square of the ratio, (γ_μ⁽ⁱ⁾)², must be used."

**Table 1** (p. 2645), "Ratios of nuclear moments used in the fits of the hyperfine parameters":

| Nuclear moment | Value | Reference |
|---|---|---|
| μ⁽¹²⁷⁾ ᵃ | 2.8090(4) | [24] |
| μ⁽¹²⁹⁾ ᵃ | 2.6173(3) | [24] |
| ¹²⁹Q / ¹²⁷Q | 0.701213(15) | [25] |
| γ_μ⁽¹²⁷⁾ | 1 | |
| γ_μ⁽¹²⁹⁾ | 0.6655 | |
| γ_Q⁽¹²⁷⁾ | 1 | |
| γ_Q⁽¹²⁹⁾ | 0.701213 | |

ᵃ In units of μ_N (nuclear magneton). [24] is H. Walchli, R. Livingston and G. Hebert, *Phys. Rev.* **82**, 97 (1951). [25] is R. Livingston and H. Zeldes, *Phys. Rev.* **90**, 609 (1953).

Check: (2.6173/3.5)/(2.8090/2.5) = 0.747800/1.123600 = 0.665540, which rounds to the printed 0.6655. Use the printed 0.6655, since the fits used it; the difference is 6×10⁻⁵ relative, or about 0.002 kHz on C. Then γ_μ² = 0.6655² = 0.44289.

### 2.2 Mass scaling: none (p. 2645)

> "The Dunham theory for the rovibrational motion also reveals dependences of the expansion parameters on the reduced mass. For the hyperfine parameters this is not straightforward. We checked the importance of mass corrections and found that within the present data set the results are not influenced by a Dunham-like dependence of the χ_lk on the reduced masses. Thus, it is omitted in the formulas given earlier. The energies E_v,J are calculated from the Dunham parameters in [23] for the case of ¹²⁷I₂ and are used without mass relations for all isotopomers. Ē_p, E_v,J and E_c are referred to E_v″=0,J″=0 of the ground state."

[23] is S. Gerstenkorn and P. Luc, *J. Phys. (France)* **46**, 867 (1985).

**Recipe.** For isotopologue *i* and level (v, J), evaluate the ¹²⁷I₂ formulae at the same integers v and J. Use E_v,J of **¹²⁷I₂** (v, J), not the isotopologue's own level energy. Then multiply by the γ factor. We do not have the Gerstenkorn–Luc 1985 Dunham set in the code. Using ¹²⁷I₂ energies from our Hannover 2008 model instead should change the denominators by far less than their size (hundreds to thousands of cm⁻¹ for v′ ≤ 53) [inference].

### 2.3 [flag] E_v,J: with or without rotational energy?

* BKT02 defines the energies at J = 0: "E(v″) is the energy of the vibrational level (v″, J″ = 0), referred to the lowest vibrational level (v″ = 0, J″ = 0), E(v″) = G_X(v″) − G_X(v″ = 0)", and likewise "E(v′) is the energy of the vibrational level (v′, J′ = 0) … E(v′) = G_B(v′) − G_X(v″ = 0)". Our `model.hyperfine_components` follows that.
* S06 writes E_v,J throughout, in eqs. (5)–(7), and says the E_v,J "are calculated from the Dunham parameters". It never says J = 0. The explicit J(J+1) terms in the numerators (χ⁽ᵖ⁾₀₁, χ⁽ᵖ⁾₁₁ for C_B; χ₂₁ for δ_B, d_B) would also absorb a J dependence of the denominator. So the text alone does not settle it.
* It matters for C_B. At J′ ≈ 100 the B-state rotational energy is about 280 cm⁻¹, while E − Ē_p is only about −3000 cm⁻¹ at v′ ≈ 11. §8 checks this numerically against the BIPM tables.
* **Resolution (§8c): use E(v, J = 0).** With S06 parameters and E at J, R(56)32-0 misses BIPM by 440–500 kHz rms and R(87)33-0 by about 3.9 MHz rms. With E at J = 0 the misses are 57 and 255 kHz. At v′ = 11 (Figs. 4–5) the choice moves the ¹²⁹I₂ and ¹²⁷I¹²⁹I peaks by less than 0.02 MHz, so the figures cannot tell the two apart.

### 2.4 The ¹²⁷I₂ formulae (reference isotopologue), §3.2 and Tables 2–4

Notation: u = v + 1/2, x = J(J+1). χ_lk multiplies u^l x^k. E ≡ E_v,J in cm⁻¹ above X(v″ = 0, J″ = 0) (see §2.3).

**General forms.** Eq. (3), the Broyer form:

χ(v, J) = χ⁽¹⁾(v, J) + Σ_p [χ⁽ᵖ⁾(v, J)/ΔE^J_v,vp]·⟨v|v⁽ᵖ⁾⟩. (3)

Eq. (4), a Dunham series: χ(v, J) = Σ_l,k χ_lk u^l x^k.

For χ = eqQ and χ = C, eq. (5):

χ(v, J) = Σ_lk χ⁽¹⁾_lk u^l x^k + Σ_p [Σ_lk χ⁽ᵖ⁾_lk u^l x^k] / (E_v,J − Ē_p), (5)

"Ē_p stands for Ē_p(3/2, 3/2) or Ē_p(1/2, 3/2)", i.e. an average perturber energy near the X-state asymptote ²P₃/₂ + ²P₃/₂ or the B-state asymptote ²P₁/₂ + ²P₃/₂.

For δ and d, eqs. (6) and (7):

δ = Σ_lk δ_lk u^l x^k + [Σ_lk χ⁽Ω⁼⁰⁾_lk u^l x^k]/(E_v,J − Ē_Ω=0) + [Σ_lk χ⁽Ω⁼¹⁾_lk u^l x^k]/(E_v,J − Ē_Ω=1) + η exp(−(E_v,J − E_c)²/W), (6)

d = Σ_lk d_lk u^l x^k − [Σ_lk χ⁽Ω⁼⁰⁾_lk u^l x^k]/(E_v,J − Ē_Ω=0) + (1/2)[Σ_lk χ⁽Ω⁼¹⁾_lk u^l x^k]/(E_v,J − Ē_Ω=1) + (η/2) exp(−(E_v,J − E_c)²/W). (7)

The Ω = 0 perturber term enters δ and d with opposite signs; the Ω = 1 term enters d with a factor 1/2. The Gaussian is the local perturbation of B by the 1u(¹Π) state crossing at "E_c ≈ 16 800 cm⁻¹".

**Table 2** (p. 2649), δ_B and d_B, eqs. (6)–(7). "For the ground state δ_X = 3.705 kHz, d_X = 1.524 kHz are used."

| Parameter | Value | Unit (as printed) |
|---|---|---|
| δ₀₀ | 20.88 | kHz |
| d₀₀ | 25.93 | kHz |
| χ⁽Ω⁼⁰⁾₀₀ | 39698 | kHz/cm⁻¹ |
| χ⁽Ω⁼⁰⁾₂₁ | 0.2475 × 10⁻³ | kHz/cm⁻¹ |
| χ⁽Ω⁼¹⁾₀₀ | 106923 | kHz/cm⁻¹ |
| χ⁽Ω⁼¹⁾₂₁ | 0.1712 × 10⁻³ | kHz/cm⁻¹ |
| η | −22.41 | kHz |
| Ē_Ω=0 | 19973 | cm⁻¹ |
| Ē_Ω=1 | 20534 | cm⁻¹ |
| E_c | 16787 | cm⁻¹ |
| W | 260267 | (cm⁻¹)² |

The unit "kHz/cm⁻¹" is as printed. Dimensionally it is kHz·cm⁻¹, because the term is divided by an energy in cm⁻¹. Written out, in kHz with E in cm⁻¹:

```
# (i) AS PRINTED (eqs. 6-7 + Table 2) -- gives unphysical values, see the flag below
δ_B = 20.88 + (39698 + 0.2475e-3·u²·x)/(E − 19973) + (106923 + 0.1712e-3·u²·x)/(E − 20534)
            − 22.41·exp(−(E − 16787)²/260267)
d_B = 25.93 − (39698 + 0.2475e-3·u²·x)/(E − 19973) + 0.5·(106923 + 0.1712e-3·u²·x)/(E − 20534)
            − 11.205·exp(−(E − 16787)²/260267)

# (ii) RECOMMENDED: sign of the Ω = 0 term flipped in both (equivalently χ⁽Ω⁼⁰⁾_lk → −χ⁽Ω⁼⁰⁾_lk)
δ_B = 20.88 − (39698 + 0.2475e-3·u²·x)/(E − 19973) + (106923 + 0.1712e-3·u²·x)/(E − 20534)
            − 22.41·exp(−(E − 16787)²/260267)
d_B = 25.93 + (39698 + 0.2475e-3·u²·x)/(E − 19973) + 0.5·(106923 + 0.1712e-3·u²·x)/(E − 20534)
            − 11.205·exp(−(E − 16787)²/260267)

δ_X = 3.705,  d_X = 1.524
```

**[flag] The printed Ω = 0 sign disagrees with BKT02 and BIPM.** A 400 dpi render confirms the print: eq. (6) has "+" before the Ω = 0 term, eq. (7) has "−", and Table 2 gives χ⁽Ω⁼⁰⁾₀₀ = 39698 with no minus sign.

* **As printed**, at v′ = 32 (E = 18 834 cm⁻¹, J′ = 57): δ_B = −78.0 and d_B = +29.9 kHz. BKT02 gives −6.0 and −45.1 kHz; that set fits BIPM R(56)32-0 to 24 kHz rms.
* **Flipped**, as in (ii): δ_B = −6.75 and d_B = −41.3 kHz.
* **BIPM test (§8c).** The printed form puts R(56)32-0 at 245 kHz rms. The flipped form gives 57 kHz, identical to S06 eqQ/C combined with BKT02 spin–spin (58 kHz). The flip also slightly improves the match to the Fig. 4 calculated curve (§8d).
* **BKT02 has the same structure** but with a negative Ω = 0 amplitude in δ (−32452/(E − 19896)) and a positive one in d.

The paper most likely changed the sign convention in the equations but not in the table, or vice versa. **Use (ii).**

**[flag]** The running text on p. 2649 says "δ_X = 3.705 kHz, d_B = 1.524 kHz are used". Given the Table 2 caption and BKT02, "d_B" there is a misprint for d_X.

Fit notes: the δ/d fit has a standard deviation of 1.6 kHz. δ_B of R30(42–0) [Chen et al.] was excluded ("due to its different sign … it seems to be perturbed"). The only isotopologue input was ¹²⁹I₂ P69(12–6), refitted from the Quinn 2003 frequency tables (see §3). "No parameters from the present measurements were included."

**Table 3** (p. 2649), C_X and C_B, eq. (5), "from the fit of nuclear spin–rotation-interaction parameters":

| State | Parameter | Value | Unit |
|---|---|---|---|
| X¹Σg⁺ | χ⁽¹⁾₀₀ | 1.9245 | kHz |
| | χ⁽¹⁾₁₀ | 0.01356 | kHz |
| | χ⁽ᵖ⁾₀₀ | −15098 | kHz/cm⁻¹ |
| | Ē_p | 12340 | cm⁻¹ |
| B³Π₀u⁺ | χ⁽¹⁾₀₀ | 28.89 | kHz |
| | χ⁽¹⁾₁₀ | 0.9234 × 10⁻¹ | kHz |
| | χ⁽¹⁾₀₁ | 0.2241 × 10⁻² | kHz |
| | χ⁽¹⁾₁₁ | 0.1027 × 10⁻⁴ | kHz |
| | χ⁽ᵖ⁾₀₀ | 36672 | kHz/cm⁻¹ |
| | χ⁽ᵖ⁾₁₀ | −3388 | kHz/cm⁻¹ |
| | χ⁽ᵖ⁾₀₁ | 9.7985 | kHz/cm⁻¹ |
| | χ⁽ᵖ⁾₁₁ | −0.2201 | kHz/cm⁻¹ |
| | Ē_p | 20140 | cm⁻¹ |

```
C_X = 1.9245 + 0.01356·u − 15098/(E − 12340)                                          [kHz]
C_B = 28.89 + 0.9234e-1·u + 0.2241e-2·x + 0.1027e-4·u·x
      + (36672 − 3388·u + 9.7985·x − 0.2201·u·x)/(E − 20140)                          [kHz]
```

C_X is identical to BKT02 eq. (10): "For the X ground state no new data exist, so those parameters are unchanged compared to [4], and were kept fixed in the fit." C_B is a new functional form: "For the upper state we had to introduce two new parameters … describing the perturbation from the electronic states sharing the same asymptote with the B state". The fit used about 640 C_B and ΔC values of ¹²⁷I₂, ¹²⁹I₂ and ¹²⁷I¹²⁹I.

**Table 4** (p. 2650), eqQ_X and eqQ_B, eq. (5) with no perturber term. "For the ratios of nuclear quadrupole moments see table 1."

| State | Parameter | Value | Unit |
|---|---|---|---|
| X¹Σg⁺ | χ⁽¹⁾₀₀ | −2452.285 | MHz |
| | χ⁽¹⁾₁₀ | −0.5474 | MHz |
| | χ⁽¹⁾₂₀ | 0.4487 × 10⁻¹ | MHz |
| | χ⁽¹⁾₀₁ | −0.2089 × 10⁻³ | MHz |
| | χ⁽¹⁾₁₁ | 0.6965 × 10⁻⁵ | MHz |
| B³Π₀u⁺ | χ⁽¹⁾₀₀ | −488.086 | MHz |
| | χ⁽¹⁾₁₀ | −1.83777 | MHz |
| | χ⁽¹⁾₃₀ | 0.99774 × 10⁻⁴ | MHz |
| | χ⁽¹⁾₅₀ | 0.80818 × 10⁻⁸ | MHz |
| | χ⁽¹⁾₀₁ | −0.14020 × 10⁻³ | MHz |
| | χ⁽¹⁾₁₁ | −0.32217 × 10⁻⁵ | MHz |
| | χ⁽¹⁾₂₁ | 0.2716 × 10⁻⁷ | MHz |
| | χ⁽¹⁾₀₂ | −0.3377 × 10⁻⁹ | MHz |

```
eqQ_X = −2452.285 − 0.5474·u + 0.4487e-1·u² − 0.2089e-3·x + 0.6965e-5·u·x                   [MHz]
eqQ_B = −488.086 − 1.83777·u + 0.99774e-4·u³ + 0.80818e-8·u⁵
        − 0.14020e-3·x − 0.32217e-5·u·x + 0.2716e-7·u²·x − 0.3377e-9·x²                     [MHz]
```

Fit notes: more than 650 eqQ_X, eqQ_B and ΔeqQ values for all three isotopologues. The standard deviation is 13 kHz, "mainly determined by the high precision data". The ground-state parameters were also refitted, which "improved the χ² of the fit by a factor of two". P78(1–9) and R113(3–10) of ¹²⁷I₂ from Dubé & Trinczek [32] were excluded. In Fig. 6, the high-precision 532 nm data (points 76–97) "have the highest weight in the fit".

For comparison, BKT02 as in `hfs_params.py` now:

* eqQ_X = −2452.2916 − 0.542u + 0.4534e-1u² − 0.1927e-3x + 0.694e-5ux.
* eqQ_B = −487.879 − 1.8621u + 0.12511e-3u³ − 0.1281e-3x − 0.225e-5ux − 0.308e-9x².
* C_B = −4.016 − 0.1501u − 3.957e-4x − 1.767e-5ux − (110704 + 1.862x)/(E − 19986).
* δ_B = 31.57 − 32452/(E − 19896) + 126257/(E − 20687) − bump.
* d_B = 19.56 + 32452/(E − 19896) + ½(126257/(E − 20687) − bump).

As printed, the S06 Ω = 0 term has the opposite sign to BKT02's: +39698/(E − Ē) in δ, where BKT02 has −32452/(E − Ē). Given the resulting magnitudes, this is a misprint (see the flag above), not a new convention.

### 2.5 Recipe for each isotopologue

Let P₁₂₇(v, J) denote the S06 ¹²⁷I₂ formulae of §2.4, with E from ¹²⁷I₂.

| | ¹²⁷I₂ | ¹²⁹I₂ | ¹²⁷I¹²⁹I |
|---|---|---|---|
| Spins (i₁, i₂) | (5/2, 5/2) | (7/2, 7/2) | (5/2 for ¹²⁷I, 7/2 for ¹²⁹I) |
| Allowed I | I ≡ J (mod 2) in X (g); I ≢ J′ (mod 2) in B (u) | Same rule | All I = 1…6 for every J |
| eqQ (MHz) | eqQ₁₂₇ | 0.701213 · eqQ₁₂₇ (both nuclei) | nucleus ¹²⁷I: eqQ₁₂₇; nucleus ¹²⁹I: 0.701213·eqQ₁₂₇, so eqQ_ratio(129/127) = 0.701213 |
| C (kHz) | C₁₂₇ | 0.6655 · C₁₂₇ | C⁽¹²⁷⁾ = C₁₂₇; C⁽¹²⁹⁾ = 0.6655·C₁₂₇ (acting on the respective nuclei) |
| δ, d (kHz) | δ₁₂₇, d₁₂₇ | 0.6655² · (δ₁₂₇, d₁₂₇) = 0.44289·(…) | **0 (omitted in S06)** |
| E in formulae | ¹²⁷I₂ E_v,J | ¹²⁷I₂ E_v,J at the same (v, J) | ¹²⁷I₂ E_v,J at the same (v, J) |
| Components (ΔF = ΔJ) | 15 (J″ even) / 21 (odd) | 28 / 36 | 48 |

In this recipe:

* δ_B and d_B take the Ω = 0 sign flip (§2.4, form (ii)).
* All E are E(v, J = 0) from ¹²⁷I₂ (§2.3).
* For ¹²⁷I₂ itself, keep BKT02 at v′ ≤ 43, since S06 C_B is worse at 532 nm (§8c). The S06 set is needed only for v′ = 44–53.
* For the isotopologues, S06 C_B is the formula the isotopologue data were fitted with, and it reproduces S06's calculated spectra (§8d).

For ¹²⁷I¹²⁹I, if spin–spin were wanted, the physical scaling would be γ_μ⁽¹²⁷⁾γ_μ⁽¹²⁹⁾ = 0.6655 for both δ and d. That is **not in S06** [inference]; S06 bounds the whole spin–spin effect at ≤ 400 kHz.

In the S06 per-line fits (§4 of the paper), the X-state parameters were held fixed at BKT02 values with eq. (8) scaling. "The nuclear spin–spin parameters also for the upper state [were calculated from the formulas in [4] with isotopic corrections] and kept … fixed for the homonuclear case in all fits, while in the heteronuclear case this interaction was not taken into account at all." This matters only if someone refits the supplement data.

---

## 3. Tables of measured or fitted parameters

* **S06 main text: none per line.** Quote (p. 2648): "Altogether, more than 120 eqQ_B parameters of ¹²⁷I₂, more than 200 for ¹²⁹I₂ and more than 170 for ¹²⁷I¹²⁹I have been determined. The same numbers apply to new C_B values. The tables of the individual hyperfine parameters will not be included here due to their length. The eqQ data are given in part 1 of the electronic supplement [27], while the new data of the spin–rotation parameter C are collected in part 2. The values in the lists are given as ΔeqQ and ΔC according to equations (9), respectively."
  * Eq. (9): ΔeqQ = eqQ_B − eqQ_X, ΔC = C_B − C_X, Δδ = δ_B − δ_X, Δd = d_B − d_X.
  * **[flag]** Reference [27] says "Part 1 contains information on the nuclear spin–rotation coupling constants C_B and C_X; Part 1 contains information on the electric quadrupole coupling constants eQq_B and eQq_X". Both are "Part 1", which contradicts the text. The text says eqQ is in part 1 and C in part 2.
  * Stated uncertainties of the per-line values (p. 2648): repeated lines use the standard deviation of the average; single measurements get estimated uncertainties. "Typical uncertainties for ΔeqQ are 2 MHz and for ΔC are 2 kHz."
  * **The supplement is not in `papers/`.** It is the only source of per-line ¹²⁹I₂ and ¹²⁷I¹²⁹I parameters. It is listed, with the S08 supplement, in `docs/research/data-availability.md`.
* Model-level fitted parameters: Tables 2–4, reproduced in full in §2.4. No uncertainties are given for the individual coefficients. Only the fit standard deviations are given: 13 kHz for eqQ, 1.6 kHz for δ/d, and none stated for C.

---

## 4. Test cases

**Neither S06 nor S08 tabulates any hyperfine component position or splitting for ¹²⁹I₂ or ¹²⁷I¹²⁹I.** The available targets are these.

### 4.1 Figure-only lines in S06 (611 nm, all from the 11–3 band region)

| Figure | Line | Isotopologue | J″ parity → components | Horizontal axis (as drawn) |
|---|---|---|---|---|
| Fig. 3 | P34(11–3) | ¹²⁷I₂ | even → 15 | −800 to +200 MHz "relative frequency" |
| Fig. 4 | R34(11–3) | ¹²⁹I₂ | even → 28 ("not all are resolved") | about −330 to +340 MHz |
| Fig. 5 | R37(11–3) | ¹²⁷I¹²⁹I | 48 ("overlap to a large degree") | about −490 to +330 MHz |

Fig. 2 (upper trace, ¹²⁹I₂ cell, absolute axis near 490 245 000–490 250 000 MHz) identifies these lines from left to right: P144(13–3) ¹²⁹I₂, R37(11–3) ¹²⁷I¹²⁹I, R110(12–3) ¹²⁹I₂, P20(9–2) ¹²⁹I₂, R34(11–3) ¹²⁹I₂, P14(9–2) ¹²⁹I₂, R25(9–2) ¹²⁷I¹²⁹I and R149(13–3) ¹²⁹I₂. The lower trace (¹²⁷I₂ cell) shows P34(11–3) ¹²⁷I₂ and R29(9–2) ¹²⁷I₂. The labels were read from a 300 dpi render.

Figs. 3–5 are vector graphics, so their traces can be digitized precisely; §4.1b gives the resulting peak tables. The figure zero is arbitrary and is not the hyperfine-free line centre: our fitted offsets are −255.0 MHz (Fig. 3), +54.8 MHz (Fig. 4) and −23.9 MHz (Fig. 5).

### 4.1b Digitized peak positions of the S06 calculated profiles

**These are not tabulated values.** They are maxima of S06's *calculated* (red) profiles, read from the PDF vector paths and calibrated on the axis tick marks. Caveats:

* **Blended peaks.** The profiles are Lorentzian (eq. 10 with b = 0; fitted a ≈ 1) with FWHM ≈ 12 MHz, so one peak can contain several components.
* **Per-line fits.** The profiles come from S06's per-line fits (fitted ΔeqQ and ΔC), not from the Table 2–4 formulas.
* **Accuracy.** Digitization is good to about 0.1 MHz (the Fig. 3 control).

Column meanings:

* "Interval" is measured from the first resolved peak.
* "Labels" lists our main components (ΔF = ΔJ, a1… by increasing frequency) within FWHM/2 of the matching model peak.
* "Model − S06" is our S06-recipe profile peak minus the digitized peak, after one common shift fitted to the whole profile. The recipe uses the Ω = 0 flip and E at J = 0.

**Recommended regression test:** convolve our components (strength^1) with a Lorentzian of the listed FWHM, find the peaks, and compare intervals from the first peak. Use a tolerance of about 0.5 MHz.

**Fig. 4: ¹²⁹I₂ R34(11–3)** (v′ = 11, v″ = 3, J″ = 34 even, 28 components). Fitted FWHM 11.9 MHz.

| Interval (MHz) | Labels | Model − S06 (MHz) |
|---|---|---|
| 0.00 | a1 | +0.05 |
| 136.75 | a2, a3 | +0.00 |
| 156.20 | a4, a5 | −0.01 |
| 217.41 | a6 | +0.08 |
| 238.39 | a7 | −0.09 |
| 247.08 | a8 | +0.18 |
| 269.61 | a9, a10 | −0.22 |
| 292.54 | a12, a13 | −0.06 |
| 388.20 | a16, a17 | −0.11 |
| 407.88 | a18 | −0.22 |
| 448.08 | a21 | −0.13 |
| 460.07 | a22 | −0.04 |
| 484.24 | a23 | −0.04 |
| 524.82 | a24, a25 | −0.01 |
| 543.62 | a26, a27 | +0.00 |
| 582.81 | a28 | +0.11 |

**Fig. 5: ¹²⁷I¹²⁹I R37(11–3)** (J″ = 37, 48 components). Fitted FWHM 12.2 MHz.

| Interval (MHz) | Labels | Model − S06 (MHz) |
|---|---|---|
| 0.00 | a1 | +0.21 |
| 21.99 | a2, a3 | −0.13 |
| 44.59 | a4 | −0.14 |
| 136.44 | a5 | +0.21 |
| 170.54 | a6, a7 | +0.29 |
| 196.05 | a8 | −0.08 |
| 236.42 | a9 | +0.05 |
| 265.81 | a10, a11 | +0.01 |
| 291.95 | a12, a13, a14 | −0.15 |
| 304.25 | a15, a16 | −0.03 |
| 323.94 | a17, a18 | +0.07 |
| 344.28 | a20 | −0.03 |
| 415.16 | a21 | +0.14 |
| 459.13 | a23–a26 | −0.10 |
| 504.70 | a28 | +0.02 |
| 518.48 | a29 | −0.08 |
| 551.13 | a30 | −0.08 |
| 562.74 | a31 | +0.05 |
| 586.89 | a33, a34, a35 | −0.04 |
| 618.42 | a37, a38 | −0.16 |
| 630.52 | a39, a40 | −0.14 |
| 682.00 | a41 | +0.11 |
| 695.49 | a42 | +0.01 |
| 706.09 | a43 | +0.04 |
| 726.34 | a44 | −0.13 |
| 745.87 | a46 | +0.07 |

**Fig. 3: ¹²⁷I₂ P34(11–3)** (control, 15 components). Fitted FWHM 11.6 MHz. The BKT02 model is used for the "Model − S06" column.

| Interval (MHz) | 0.00 | 281.59 | 301.94 | 413.84 | 438.66 | 458.63 | 581.90 | 715.25 | 740.99 | 873.84 |
|---|---|---|---|---|---|---|---|---|---|---|
| Labels | a1 | a2, a3 | a4, a5 | a6 | a8 | a9 | a10 | a11, a12 | a13, a14 | a15 |
| Model − S06 (MHz) | −0.14 | +0.00 | −0.10 | −0.05 | −0.13 | +0.10 | −0.01 | +0.11 | −0.13 | +0.02 |

### 4.2 Consistency test stated by S06

The heteronuclear code must reproduce the homonuclear code "within numerical accuracy … for situations where the results should coincide". For us: with i₁ = i₂, eqQ_ratio = 1, C⁽¹⁾ = C⁽²⁾ and `symmetry=None`, the level set must be the union of the "g" and "u" (or I-even and I-odd) homonuclear level sets.

### 4.3 External numeric targets named by S06 (not local)

* **[11] T. J. Quinn, *Metrologia* 40, 103 (2003)** (the CIPM 2001 *mise en pratique*). S06 says: "High precision splittings of few lines around 633 nm are tabulated in [11]." Also: "From the frequency tables in [11] hyperfine parameters for ¹²⁹I₂ were fitted, but only the data for the P69 (12–6) line of ¹²⁹I₂ are sufficiently accurate and complete for fitting reliable nuclear spin–spin parameters." **The ¹²⁹I₂ P(69) 12–6 components at 633 nm are the best regression target.** It is open access (doi:10.1088/0026-1394/40/2/316); the 633 nm BIPM tables hold the only CIPM-level isotopologue anchors.
* [10] M. Tesic and Y. H. Pao, *J. Mol. Spectrosc.* 57, 75 (1975): "the hyperfine structure near 633 nm had been unravelled" (¹²⁷I¹²⁹I context).
* [9] G. W. King et al., *Chem. Phys.* 50, 291 (1980): two B-state vibrational levels of ¹²⁷I¹²⁹I.
* [12] M. Klug et al., *Opt. Commun.* 184, 215 (2000): "some accurate information on the level structure of ¹²⁹I₂ and ¹²⁷I¹²⁹I" from Raman lasing at 532 nm.
* S06 supplement [27]: per-line ΔeqQ and ΔC for more than 200 ¹²⁹I₂ and more than 170 ¹²⁷I¹²⁹I lines (v′ = 8–20, v″ = 0–5).

---

## 5. Validity ranges and caveats

* **Model range** (§5, p. 2648): "we restrict here the ranges of our description to vibrational levels 0 ≤ v″ ≤ 17 for the X¹Σg⁺ state and 0 ≤ v′ ≤ 53 for the B³Π₀u⁺ state." Perturbations appear above v′ = 55, and there are "some irregularities" for v′ ≥ 42 attributed to the 1g(¹Πg) state. Our `E_B_MAX` (v′ ≤ 43, from BKT02) can move up to v′ = 53 with the S06 formulae. Check the pole distances: Ē_Ω=0 = 19973, Ē_p(C_B) = 20140, Ē_Ω=1 = 20534 cm⁻¹.
* **New isotopologue data:** "8 ≤ v′ ≤ 20 in the upper state and 0 ≤ v″ ≤ 5 in the lower state. The range of the rotational quantum number J″ is from 10 to 150 for the ¹²⁹I₂ isotopomer, while the ranges for ¹²⁷I₂ and ¹²⁷I¹²⁹I are a bit smaller" (p. 2647). **[flag]** §4 (p. 2646) says "v′ from 9 to 20 for all three isotopomers". S08 gives J″ 23–128 for its own data.
* **Isotopologue data are sparse outside v′ = 8–20.** At other v′ the isotopologue parameters rest entirely on ¹²⁷I₂ data plus eq. (8).
* **Prediction uncertainty:**
  * ¹²⁷I₂: "less than 60 kHz … for the prediction of component a10 (even J″) resp. a13 (odd J″)".
  * ¹²⁹I₂ and ¹²⁷I¹²⁹I: "we expect a prediction uncertainty of hyperfine splittings of ¹²⁹I₂ and ¹²⁷I¹²⁹I not larger than 1 MHz" (p. 2651).
  * **[flag]** The conclusion says "For the isotopomers ¹²⁹I₂ and ¹²⁷I¹²⁹I the prediction uncertainty is less due to only a few high precision measurements". From context, "less" means *less good*.
* **Spin–spin neglect for ¹²⁷I¹²⁹I:** at most 400 kHz on the splittings in their range. For ¹²⁷I₂ the spin–spin contribution to a2 is 200 kHz at v′ = 11 and 1.1 MHz for P89(53–0). It is similar for a20 and about half for a13 (odd J″); for even J″ the corresponding components are a1, a10 and a15. Since spin–spin in ¹²⁷I¹²⁹I would scale by about 0.67, omitting it is an error of up to a few hundred kHz on some components at high v′ [inference].
* **Sensitivities for ¹²⁷I₂** (§5.2–5.3):
  * The components least sensitive to C are a1 (a2), a10 (a13) and a15 (a20) for even (odd) J″. **[flag]** One sentence on p. 2650 writes "a10 (a12)" and p. 2651 writes "a10 (a12)" once, then "a13 (odd J″)"; a13 is meant.
  * A 1 MHz change in eqQ shifts a1 (a2) by 0.25 MHz, a10 (a13) by 0.05 MHz and a15 (a20) by 0.2 MHz.
  * Only ΔeqQ and ΔC are well determined from ΔF = ΔJ lines with J ≥ 10 (the B and X parameters are strongly correlated). The low-J, ΔF = 0 and ΔF = −ΔJ or crossover data are what separate X from B.
* **Reference line for ¹²⁹I₂:** "Only for ¹²⁹I₂ and even J″, one hyperfine component, a1, is sufficiently separated from the other ones to be used as a reference line free from overlap."
* **Energies:** see §2.3. The ¹²⁷I₂ Gerstenkorn–Luc 1985 Dunham energies are used for all isotopologues.

---

## 6. Code changes needed in `src/i2spec`

1. **Per-nucleus spin–rotation (required for ¹²⁷I¹²⁹I).** Replace C·I·J by C₁ I₁·J + C₂ I₂·J, e.g. with a `C_ratio = C(nucleus 2)/C(nucleus 1)` argument mirroring `eqQ_ratio` (0.6655 when nucleus 1 is ¹²⁷I). This term is diagonal in J but **not in I** when i₁ ≠ i₂ or C₁ ≠ C₂. In the coupled basis:

   ⟨J (i₁i₂)I′; F | Iₙ·J | J (i₁i₂)I; F⟩ = (−1)^(J+I′+F) {F I′ J; 1 J I} · √(J(J+1)(2J+1)) · ⟨I′‖Iₙ‖I⟩,

   where ⟨I′‖Iₙ‖I⟩ is exactly `_one_nucleus(n, I′, I, i1, i2, 1, _vec_reduced(i_n))`. This is the same pattern as the rank-2 terms, with rank 1 and the rotational factor ⟨J‖J‖J⟩. The phase was verified numerically (§8) against a brute-force product-space construction. With C₁ = C₂ it reduces to the present diagonal C[F(F+1) − I(I+1) − J(J+1)]/2.
2. **ΔI = ±1 couplings.** The quadrupole already produces them through `_one_nucleus` when eqQ_ratio ≠ 1 or i₁ ≠ i₂. The module docstring ("ΔI = ±2") should say so. `allowed_spins(..., symmetry=None)` already returns I = 1…6.
3. **Spin–spin for ¹²⁷I¹²⁹I:** set δ = d = 0 to match S06. The existing scalar formula and the 9j tensor element already handle i₁ ≠ i₂ if an extension is wanted later.
4. **New parameter set** (e.g. `hfs_params.salumbides2006`) with these contents:
   * Tables 2–4, with the Ω = 0 sign flip in δ_B and d_B (§2.4, form (ii)).
   * E = E(v, J = 0) of **¹²⁷I₂** above X(0, 0) (§2.3). So `RovibronicModel("129I2")` needs the ¹²⁷I₂ vibrational energies G_X(v) and G_B(v); a small cached table is enough.
   * The eq. (8) scaling.

   Keep BKT02 as the default for ¹²⁷I₂ at v′ ≤ 43, because its C_B is better at 532 nm (§8c). Use S06 for ¹²⁹I₂, for ¹²⁷I¹²⁹I, and for ¹²⁷I₂ at v′ = 44–53.
5. **`model.hyperfine_components`:** remove the `NotImplementedError`. Pass i₁ and i₂, symmetry ("g"/"u" for the homonuclear molecules, `None` for ¹²⁷I¹²⁹I), eqQ_ratio, C_ratio and the scaled parameters. Raise the v′ ceiling to 53 with the S06 set.
6. **`line_components`** needs no change. The dipole acts with I as a spectator (`Ia == Ib`), and the ΔF = ΔJ strongest-first labelling gives 48 labels for ¹²⁷I¹²⁹I.
7. **Tests** to add:
   * the S06 homonuclear/heteronuclear consistency test (§4.2);
   * the component counts (15/21, 28/36, 48);
   * the digitized S06 Fig. 4 and Fig. 5 peak intervals (§4.1b, tolerance about 0.5 MHz);
   * the single-C variant as a negative control (it misses Fig. 5 by up to 1.8 MHz);
   * the Quinn 2003 ¹²⁹I₂ P(69) 12–6 intervals once that table is transcribed.

---

## 7. Salumbides et al. 2008 (EPJD 47, 171): hyperfine content

Nothing new. S08 fits the hyperfine structure of each line only to extract hyperfine-free line centres. "The details of the hyperfine structure models have been discussed in [12] and will not be repeated here", where [12] is S06. For literature lines with too few components, "the hyperfine splitting was calculated taking the models from [12] and subtracted to yield the hyperfine-structure-free frequencies. The uncertainty introduced by this prediction of hyperfine splitting is less than 60 kHz [12]." The absolute calibration used the a1 ("t") component of ¹²⁷I₂ P(72)(11–3). S08 has no tables of isotopologue hyperfine components.

---

## 8. Numerical checks (i2spec, this study)

Energies come from our Hannover 2008 model and are referred to X(0, 0). The line builder used here is the package's `matrix_element` plus the per-nucleus spin–rotation element of §6. With BKT02 parameters it reproduces `RovibronicModel.hyperfine_components` exactly (0.000 kHz on R(56)32-0).

**(a) Per-nucleus spin–rotation.** Setup:

* The §6 coupled-basis element plus the existing quadrupole.
* (i₁, i₂) = (5/2, 7/2), eqQ_ratio = 0.701213, C₂/C₁ = 0.6655.
* Compared with brute-force diagonalization in the product space |m_J m₁ m₂⟩ at J = 3, 4 and 7.

All 336, 432 and 720 eigenvalues agree to 2×10⁻¹². With C₁ = C₂ the element reduces to C[F(F+1) − I(I+1) − J(J+1)]/2 on the diagonal and 0 off the diagonal (to 3×10⁻¹³).

**(b) ¹²⁷I₂ parameter values.** E = E(v, J = 0). BKT02 B-state values above v′ = 43 are frozen at `E_B_MAX`, as the package does.

| Level | eqQ (MHz) BKT02 / S06 | C (kHz) BKT02 / S06 | δ (kHz) BKT02 / S06 printed / S06 flipped | d (kHz) BKT02 / S06 printed / S06 flipped |
|---|---|---|---|---|
| X v=0, J=56 | −2453.155 / −2453.203 | 3.155 / 3.155 | 3.705 | 1.524 |
| B v=11, J=35 | −509.297 / −509.286 | 31.397 / 30.736 | −10.23 / −41.53 / −14.77 | −18.16 / +14.79 / −11.96 |
| B v=32, J=57 | −544.771 / −544.814 | 89.333 / 89.935 | −6.01 / −77.97 / −6.75 | −45.06 / +29.92 / −41.30 |
| B v=33, J=88 | −547.168 / −547.285 | 98.468 / 98.826 | −6.49 / −84.42 / −6.45 | −48.29 / +31.76 / −46.21 |
| B v=43, J=60 | −559.413 / −559.400 | 197.965 / 198.318 | 1.03 / −154.16 / −0.59 | −101.00 / +53.59 / −99.98 |
| B v=53, J=88 | −570.308 / −569.446 | (211.2, frozen) / 490.376 | (2.27) / −349.27 / +100.41 | (−104.24) / +178.12 / −271.57 |

At v′ = 53, E = 19 772 cm⁻¹ is only 201 cm⁻¹ below Ē_Ω=0 = 19 973 cm⁻¹. The δ_B and d_B there are dominated by that pole and are very sensitive to E. S06 does fit data up to v′ = 53, and says the spin–spin contribution to a2 of P89(53–0) is 1.1 MHz.

**(c) BIPM 532 nm intervals.** Data are from `tests/test_bipm_532.py`: R(56)32-0 has 13 intervals from a10 (u = 1.5 kHz), and R(87)33-0 has 21 intervals from a1 (u = 2 kHz). All runs use ΔJ = ±2 mixing.

| Parameter set | E | R(56)32-0 rms / max (kHz) | R(87)33-0 rms / max (kHz) |
|---|---|---|---|
| BKT02 (package now) | J = 0 | 23.7 / 53.7 | 88.9 / 146.6 |
| S06 as printed | J = 0 | 244.8 / 480.7 | 247.9 / 447.1 |
| S06 as printed | J | 498.9 / 1042.8 | 3880.8 / 7016.5 |
| S06, Ω = 0 flipped | J = 0 | 57.0 / 110.9 | 255.0 / 463.7 |
| S06, Ω = 0 flipped | J | 443.8 / 780.5 | 3990.7 / 7036.1 |
| S06 eqQ and C, BKT02 δ and d | J = 0 | 57.6 / 112.5 | 267.1 / 464.1 |
| S06 eqQ only (rest BKT02) | J = 0 | 23.7 / 54.4 | 89.9 / 146.7 |
| S06 C only (rest BKT02) | J = 0 | 57.5 / 111.7 | 266.1 / 464.0 |

Conclusions:

1. E must be E(v, J = 0).
2. The printed Ω = 0 sign must be flipped. The flipped δ_B and d_B are as good as BKT02's.
3. S06 eqQ is equivalent to BKT02 at 532 nm.
4. All the remaining loss comes from S06 C_B. Relative to BKT02 it is +0.60 kHz at v′ = 32, J′ = 57 and +0.36 kHz at v′ = 33, J′ = 88, and those shifts cost 30–180 kHz rms. S06 did not claim kHz-level C_B predictions: it quotes about 1 kHz uncertainty on ΔC and less than 15 kHz only for the C-insensitive components.

**(d) S06 figures: calculated profiles reproduced.** Method:

* Figs. 3–5 are vector paths, not raster images (`pdfimages` finds no images on pp. 2646–2648).
* The calculated (red) and observed (black) polylines were extracted with PyMuPDF.
* The axis tick marks give 4.834, 3.225 and 4.000 MHz/pt, with linear-fit residuals ≤ 0.18 MHz.
* Our components were convolved with the S06 eq. (10) profile (b = 0). One global shift, the FWHM, a, an intensity exponent p (strength^p) and a quadratic background were fitted to each red curve, and the resolved peaks compared.

| Figure / line | Parameters | Fit rms / peak height | Peaks matched | Model − red, rms / max (MHz) |
|---|---|---|---|---|
| Fig. 3, ¹²⁷I₂ P34(11–3) (control) | BKT02 | 0.0045 | 10/10 | 0.09 / 0.14 |
| Fig. 3 | S06, Ω = 0 flipped | 0.0048 | 10/10 | 0.13 / 0.23 |
| Fig. 4, ¹²⁹I₂ R34(11–3) | S06, flipped, E at J = 0 | 0.0052 | 16/16 | 0.11 / 0.22 |
| Fig. 4 | S06 as printed | 0.0054 | 16/16 | 0.13 / 0.27 |
| Fig. 5, ¹²⁷I¹²⁹I R37(11–3) | S06, flipped, per-nucleus C | 0.0024 | 26/26 | 0.12 / 0.29 |
| Fig. 5 | S06 as printed (no heteronuclear δ, d, so same) | 0.0024 | 26/26 | 0.12 / 0.29 |
| Fig. 5 | one C for both nuclei (C⁽¹²⁹⁾ = C⁽¹²⁷⁾) | 0.0154 | 24/26 | 0.89 / 1.84 |

The fits give FWHM 11.6–12.2 MHz, a = 0.9–1.0 (i.e. Lorentzian) and p = 1.0 (intensity proportional to our line strength). At the matched peaks the observed (black) traces scatter about the red curves by 0.36, 0.64 and 1.41 MHz rms (Figs. 3, 4 and 5), so the calculated curve is the more useful target. Using E at J instead of J = 0 changes Figs. 4 and 5 by at most 0.01 MHz.

Conclusions:

1. The digitization is accurate to about 0.1 MHz (Fig. 3 control).
2. Our implementation of the S06 recipe reproduces S06's calculated ¹²⁹I₂ and ¹²⁷I¹²⁹I profiles to 0.1 MHz rms. That confirms together:
   * i = 7/2;
   * γ_Q = 0.701213 and γ_μ = 0.6655;
   * per-nucleus eqQ and C;
   * no heteronuclear spin–spin;
   * the equivalence of the (JI₁)F₁I₂F and (I₁I₂)IJF bases.
3. **Per-nucleus spin–rotation is required** and clearly resolved at this level.
4. The remaining 0.1–0.3 MHz is consistent with the difference between S06's per-line fitted parameters (behind the red curves) and the Table 2–4 formula values we used. It is well inside S06's stated prediction uncertainty of ≤ 1 MHz for the isotopologue splittings.

---

## 9. Implementation status (2026-09-15)

§6 is implemented.

* **`hyperfine.py`.** A `C_ratio` argument (C of nucleus 2 over C of nucleus 1) sits alongside `eqQ_ratio`.
  * The operator is written as C·I·J + (C_ratio − 1)·C·I₂·J, with the §6.1 element for the second term.
  * With C_ratio = 1 the code path is unchanged.
* **`hfs_params.py`.**
  * S06 Tables 2–4 are the `*_s06` functions, with the Ω = 0 sign flipped in δ_B and d_B.
  * Table 1 is in `GAMMA_Q` and `GAMMA_MU`; `nucleus_ratios` derives the per-nucleus ratios from them.
  * `line_states` chooses the parameter set per line: BKT02 for ¹²⁷I₂ at v′ ≤ 43 (`E_B_MAX`), S06 otherwise. Above v′ = 53 the B-state energy is frozen at `E_B_MAX_S06` = 19 785 cm⁻¹.
  * For ¹²⁷I¹²⁹I, nucleus 1 is ¹²⁷I, and δ = d = 0.
* **`model.py`.** `hyperfine_components` handles all three isotopologues. The formula energies are ¹²⁷I₂ E(v, J = 0) above X(0, 0), taken from a ¹²⁷I₂ model with the same parameters, grids and solver, built on first use.
* **Effect on ¹²⁷I₂.**
  * Lines with v′ ≤ 43 are bit-identical to before; checked on R(56)32-0, R(87)33-0, P(60)43-0 and P(34)11-3.
  * Lines with v′ = 44–53 now use S06 instead of frozen BKT02.
* **Tests** (`tests/test_isotopologue_hyperfine.py`):
  * the heteronuclear Hamiltonian against a brute-force product space, for (i₁, i₂) = (5/2, 7/2) in both orders, with ΔJ = 0 and ±2 and all four terms;
  * the §4.2 consistency test;
  * the component counts;
  * the §8(b) S06 B-state values;
  * the Fig. 4 and Fig. 5 peak intervals, with tolerance 0.5 MHz;
  * the single-C negative control.
* **Achieved.** These are peak intervals from the first peak, with a pure Lorentzian at the fitted FWHM.

  | Figure | From the first peak, rms / max | After removing the mean shift (as in §8d), rms / max |
  |---|---|---|
  | Fig. 4 | 0.13 / 0.27 MHz | 0.11 / 0.21 MHz |
  | Fig. 5 | 0.24 / 0.37 MHz | 0.12 / 0.29 MHz |

  With one C, the residuals reach 18.5 MHz, because some tabulated peaks have no model counterpart.
* **Not done.** The Quinn 2003 ¹²⁹I₂ P(69) 12–6 test, because that table is not yet transcribed, and the S06 supplement data.
* **BIPM tables (2026-09-15).** ¹²⁹I₂ P(69) 12–6 is now tested, together with five other ¹²⁹I₂ lines and ¹²⁷I¹²⁹I P(33) 6–3, using the BIPM MEP 2003 tables, which supersede Quinn 2003. The tests are in `tests/test_bipm_other.py`; residuals, and the 5–14 MHz isotope-shift error they revealed, are in `bipm-hyperfine-tables.md`.
