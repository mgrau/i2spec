# BIPM hyperfine tables for the other iodine standards (515–640 nm)

*2026-09-15. Tests: `tests/test_bipm_other.py`; the 532 nm standard stays in `tests/test_bipm_532.py`. Model: Hannover 2008 potentials, with Bodermann 2002 hyperfine (BKT02) for ¹²⁷I₂ at v′ ≤ 43 and Salumbides 2006 (S06) otherwise (`isotopologue-hyperfine.md`).*

## Sources

BIPM, *Recommended values of standard frequencies* (mise en pratique). The URLs are in `data/catalog/precision.yaml`. All were retrieved 2026-09-14.

| λ | Edition | URL | Hyperfine tables used |
|---|---|---|---|
| 515 nm | MEP 2005 | https://www.bipm.org/documents/20126/41549496/M-e-P_I2_515.pdf/63b35d2c-379d-aa88-8480-7a44395637b5 | ¹²⁷I₂ P(13) 43-0, R(15) 43-0, R(98) 58-1 (kHz) |
| 543 nm | MEP 2003 | https://www.bipm.org/documents/20126/41549533/M-e-P_I2_543.pdf/6791a64e-6a7b-bcd9-308a-bf682306ba89 | ¹²⁷I₂ R(12) 26-0, R(106) 28-0 |
| 576 nm | MEP 2003 | https://www.bipm.org/documents/20126/41549542/M-e-P_I2_576.pdf/bd97f28c-f030-c38f-f1be-a380eea5129d | ¹²⁷I₂ P(62) 17-1 |
| 612 nm | MEP 2003 | https://www.bipm.org/documents/20126/41549551/M-e-P_I2_612.pdf/e3257540-b5f1-596f-9626-0531add4912a | ¹²⁷I₂ R(47) 9-2, P(48) 11-3, R(48) 15-5; ¹²⁹I₂ P(110) 10-2, R(113) 14-4 |
| 633 nm | MEP 2003 | https://www.bipm.org/documents/20126/41549560/M-e-P_I2_633.pdf/c4c25f25-ae65-e05d-402a-9bfc84c715c3 | ¹²⁷I₂ R(127) 11-5, P(33) 6-3; ¹²⁹I₂ P(54) 8-4, P(69) 12-6, R(60) 8-4, P(33) 6-3; ¹²⁷I¹²⁹I P(33) 6-3 |
| 640 nm | MEP 2003 | https://www.bipm.org/documents/20126/41549593/M-e-P_I2_640.pdf/89d9d1e5-aa97-3bff-d9d1-49e138f1fc13 | ¹²⁷I₂ P(10) 8-5, R(16) 8-5 |
| 531 nm | CIPM 2015 | https://www.bipm.org/documents/20126/41549505/127I2_531nm_2015.pdf/3f08c121-ed26-3268-89e9-75b9eb9fc663 | none (one recommended value only) |

**Quinn 2003 is not needed.** The MEP 2003 files say their tables "replace those published in BIPM Com. Cons. Long., 2001, 10, … and Metrologia, 2003, 40, 121–128". That includes the ¹²⁹I₂ P(69) 12-6 table that `isotopologue-hyperfine.md` §4.3 wanted from Quinn 2003.

## Transcription notes

* Values were taken from the PDF text layer and checked against the table layout. Minus signs are present in the text layer.
* Every value and u_c is entered exactly as printed. The 515 nm tables are in kHz and were converted to MHz, which is exact.
* **Labels.**
  * The prefix (a, b, c, d, e, m) only names the line. The index counts the ΔF = ΔJ components by increasing frequency, which is how i2spec numbers its labels, so each entry is stored as aN.
  * The mapping is confirmed by the residuals, including the partial tables: c1–c9 with c7 absent (15 components), b19–b36 (36), d23–d28 (28), a1–a10 plus a15 (576 nm), and e1–e17. For e1–e17, a block match against the 36 model components gives 0.02 MHz rms for e1–e4 only when e1 = a1.
* **Blends at 633 nm.**
  * The x-column letters form one sequence across the mixed ¹²⁹I₂/¹²⁷I¹²⁹I-cell spectrum. A letter that appears in both isotopologue tables is one observed feature assigned twice.
  * Those features are r″ and q″ (with explicit "see" notes), k″ and d″, plus n (¹²⁷I¹²⁹I m8 = ¹²⁹I₂ P(54) a24/a25, also with a note).
  * They are excluded from the tests: ¹²⁹I₂ P(54) a24–a25; ¹²⁹I₂ P(69) a5, a6, a8, a12; ¹²⁷I¹²⁹I a8, a28, a29, a34, a35, a40.
* Reference components with no u_c (the zero of each table) are entered as u_c = 0.

## Method

1. Compute the components with `hyperfine_components(..., dJ=2)`.
2. For each component, r = model − table.
3. Remove the common offset: the mean of r over the well-measured components. These are u_c ≤ 0.1 MHz, or ≤ 0.5 MHz if that leaves fewer than three, or else all components.
4. s is the rms of r over those same components.

**The test** requires |r| ≤ 3·√(s² + u_c²) for every unblended component. Each line's s in the test file is the value below, rounded up.

## Residuals per line

"n" is the number of well-measured components over the total, excluding blends. The rms and max are over the well-measured ones.

| Line | Parameters | n | rms (kHz) | max (kHz) | Notes |
|---|---|---|---|---|---|
| ¹²⁷I₂ R(127) 11-5, 633 nm | BKT02 | 18/20 | 28.5 | 99.5 | a2, a3 (u_c 0.5 MHz) within 0.7σ |
| ¹²⁷I₂ P(33) 6-3, 633 nm | BKT02 | 21/21 | 37.1 | 77.0 | |
| ¹²⁷I₂ P(13) 43-0, 515 nm | BKT02 | 17/21 | **2.8** | 5.9 | Blends a10, a14, a18, a19 (u_c 1 MHz) within 0.4σ |
| ¹²⁷I₂ R(15) 43-0, 515 nm | BKT02 | 17/21 | **2.2** | 4.1 | |
| ¹²⁷I₂ R(98) 58-1, 515 nm | S06, frozen at v′ = 53 | 13/15 | 34 000 | 54 000 | **Fails; xfail.** Still 6 MHz rms without the freeze. See finding 4 |
| ¹²⁷I₂ R(12) 26-0, 543 nm | BKT02 | 15/15 | 22.3 | 48.5 | |
| ¹²⁷I₂ R(106) 28-0, 543 nm | BKT02 | 15/15 | 40.0 | 74.5 | |
| ¹²⁷I₂ P(62) 17-1, 576 nm | BKT02 | 11/11 | 96.8 | 207.6 | A fitted ΔeqQ_B, ΔC_B does not help (92 kHz) |
| ¹²⁷I₂ R(47) 9-2, 612 nm | BKT02 | 21/21 | 100.5 | 197.0 | ΔC_B = −0.80 kHz gives 6 kHz |
| ¹²⁷I₂ P(48) 11-3, 612 nm | BKT02 | 13/15 | 69.6 | 135.9 | |
| ¹²⁷I₂ R(48) 15-5, 612 nm | BKT02 | 8/8 | 24.4 | 40.8 | |
| ¹²⁷I₂ P(10) 8-5, 640 nm | BKT02 | 7/15 | 81.4 | 175.1 | Older data, u_c up to 1.3 MHz |
| ¹²⁷I₂ R(16) 8-5, 640 nm | BKT02 | 3/3 (u_c ≤ 0.5 MHz) | 26.7 | 37.0 | Partial: b1–b3 |
| ¹²⁹I₂ P(54) 8-4, 633 nm | S06 | 11/25 | 42.9 | 84.9 | Fitted ΔeqQ_B, ΔC_B: 22 kHz |
| ¹²⁹I₂ P(69) 12-6, 633 nm | S06 | 24/29 | 160.4 | 516.4 | a21 (q′) −0.52 MHz, 5σ of its u_c. See finding 2 |
| ¹²⁹I₂ R(60) 8-4, 633 nm | S06 | 6/6 (u_c 2–5 MHz) | 725 | 1297 | Limited by the data (every point within 0.7σ) |
| ¹²⁹I₂ P(33) 6-3, 633 nm | S06 | 4/17 | **17.3** | 26.9 | e5–e17 (u_c 2–6 MHz) within 1.8σ |
| ¹²⁹I₂ P(110) 10-2, 612 nm | S06 | 11/28 | 86.8 | 139.5 | Fitted ΔeqQ_B, ΔC_B: 12 kHz |
| ¹²⁹I₂ R(113) 14-4, 612 nm | S06 | 15/18 (u_c ≤ 0.5 MHz) | 214.5 | 470.9 | Partial: b19–b36 |
| ¹²⁷I¹²⁹I P(33) 6-3, 633 nm | S06, per-nucleus C, no spin–spin | 31/42 | 107.4 | 422.6 | a47 (m47, r′) +0.42 MHz, 8σ of its u_c: listed as an outlier, checked to < 1 MHz |

## Absolute and cross-line checks (model − table, MHz)

**CIPM recommended components** (tolerance 4 MHz, about 3× the 1.3 MHz rms):

| Component | Model − CIPM |
|---|---|
| R(127) 11-5 a16, 633 nm | +1.18 |
| P(13) 43-0 a3, 515 nm | +1.40 |
| R(106) 28-0 b10, 543 nm | −1.23 |
| P(62) 17-1 a1, 576 nm | −1.90 |
| R(47) 9-2 a7, 612 nm | +0.81 |
| P(10) 8-5 a9, 640 nm | −0.86 |
| (R(56) 32-0 a10, 532 nm, from `test_bipm_532.py`) | (+2.00) |

**Intervals between lines of one isotopologue** (tolerance 5 MHz, about 3× the 1.6 MHz rms):

| Interval | Model − table |
|---|---|
| ¹²⁷I₂ P(33) b21 − R(127) a16 | −0.87 |
| R(15) b1 − P(13) a1 (u_c 5 MHz) | +0.33 |
| R(12) a4 − R(106) b10 | +0.51 |
| P(48) b15 − R(47) a7 | +0.05 |
| R(48) 15-5 c1 − R(47) a7 | −4.21 |
| R(16) b1 − P(10) a9 | −0.06 |
| ¹²⁹I₂ P(69) b27 − P(54) a28 | −1.40 |
| ¹²⁹I₂ R(60) d28 − P(54) a28 | −0.38 |

**Isotope shifts** (`test_isotope_shifts`, xfail, strict):

| Interval | Table (u_c) | Model − table |
|---|---|---|
| ¹²⁹I₂ P(54) a28 − ¹²⁷I₂ R(127) a16 | −42.99 (0.04) | **−10.95** |
| ¹²⁹I₂ P(33) e2 − ¹²⁷I₂ R(127) a16 | 849.4 (0.2) | **−11.9** |
| ¹²⁹I₂ P(110) a1 − ¹²⁷I₂ R(47) a7 | −376.29 (0.05) | **−7.39** |
| ¹²⁹I₂ R(113) b36 − ¹²⁷I₂ R(47) a7 | −187.0 (0.3) | **−13.91** |
| ¹²⁷I¹²⁹I P(33) m10 − ¹²⁹I₂ P(54) a28 | 12.04 (0.03) | **+4.51** |

## Findings

1. **¹²⁷I₂ hyperfine at v′ ≤ 43 is good to 2–100 kHz rms.** That is consistent with the uncertainties BKT02 states. The two 515 nm lines at v′ = 43 reach 2–3 kHz.
   * The 612 nm R(47) 9-2 residual is almost entirely a C_B offset. A fitted ΔC_B of −0.80 kHz, within BKT02's ±2 kHz (2σ), leaves 6 kHz.
   * The 576 nm P(62) 17-1 residual is not an eqQ or C error: fitting both leaves 92 kHz. Its residuals compress each group of components spaced about 12 MHz apart by 0.35–0.55%. [inference] That could be line-shape pulling in the 1979/1984 third-harmonic measurements.
2. **The isotopologue implementation (S06 recipe) is structurally right.**
   * With a fitted ΔeqQ_B and ΔC_B, ¹²⁹I₂ P(54) 8-4 and P(110) 10-2 reach 22 and 12 kHz. ¹²⁹I₂ P(33) 6-3 (e1–e4) is at 17 kHz without any fit.
   * The isotopologue data independently confirm the Ω = 0 sign flip in δ_B and d_B. With the printed sign, P(54) and P(110) go to 79 and 104 kHz, against 43 and 87 kHz with the flip.
   * **¹²⁹I₂ P(69) 12-6 stays at 160 kHz rms**, and fitting ΔeqQ_B and ΔC_B barely changes it. That was checked before k″ and d″ were excluded.
     * Its worst components are a21–a23 (q′, o′, n′: −0.52, +0.27, −0.31 MHz). They sit within a few MHz of ¹²⁷I¹²⁹I r′, s′ and p′ in the mixed-cell spectrum.
     * The ¹²⁷I¹²⁹I outlier m47 (r′, +0.42 MHz) is the neighbour of ¹²⁹I₂ q′. The table puts the pair 4.1 MHz apart, the model 3.2 MHz. [inference] That pattern would fit 3f line shapes interfering between overlapping components of the two isotopologues.
     * All of this is within the ≤ 1 MHz that S06 quotes for isotopologue predictions.
   * **¹²⁷I¹²⁹I P(33) 6-3 is at 107 kHz rms.**
     * Adding the heteronuclear spin–spin terms, with δ and d scaled by γ₁₂₇γ₁₂₉ = 0.6655, lowers it to 95 kHz. S06 omits these terms. a12 goes from +91 to +16 kHz and a5 from +120 to +24 kHz; a44 and a47 do not move.
     * Not adopted, because it is not the S06 recipe. It is a candidate for the refit.
3. **Isotope shifts are off by 5–14 MHz.** ¹²⁷I₂ absolute positions agree within 2 MHz, so this is specific to the isotopologues. It is in the rovibronic model (potentials and Born–Oppenheimer corrections), not in the hyperfine code.
   * ¹²⁹I₂ lines come out 7–14 MHz low relative to ¹²⁷I₂, and ¹²⁷I¹²⁹I about 5 MHz high relative to ¹²⁹I₂ (about −5 MHz relative to ¹²⁷I₂).
   * For scale: switching off the B-state adiabatic correction V_corr moves the ¹²⁹I₂ − ¹²⁷I₂ differences by +93 and +124 MHz. So the errors match a V_corr that is 6–12% too large, or an equivalent isotope-dependent term. [inference] This needs checking against the Salumbides 2008 BOC conventions: `vad_power` = 5 for B, the reference mass, and the X state, which has no BOC terms in `hannover2008.json`.
   * The masses used are NIST: ¹²⁷I 126.9044719 u, ¹²⁹I 128.9049837 u.
4. **v′ = 58 is out of reach.** ¹²⁷I₂ R(98) 58-1 at 515 nm is measured 2.10 GHz below P(13) 43-0 a3. The model puts it 43.38 GHz below, a 1.38 cm⁻¹ error in B(v′ = 58, J′ = 99).
   * This is not numerical: the B-spline and DVR solvers agree to 10⁻⁴ cm⁻¹, and the tabulated and continuity-matched outer extensions give identical results.
   * The v′ = 58 outer turning point (≈ 5.5 Å) lies beyond the B-state R_O = 4.9 Å.
   * The hyperfine intervals miss by up to 54 MHz with the v′ = 53 freeze, and by 6 MHz rms without it.
   * This matters for any spectrum above about 19 800 cm⁻¹.
5. **All 20 lines at 532 nm** (`data/observations/bipm2012a`, 329 observations; added 2026-09-15):
   * **Hyperfine intervals:** 18 lines at 21–135 kHz rms.
   * **R(145) 37-0 is at 1.24 MHz rms** (max 2.1 MHz):
     * every interval from a1 is larger than the model, by 0.1–1%;
     * a2 and a3 are off by +1.1 and +2.1 MHz.
   * **P(142) 37-0** is at 0.33 MHz.
   * The neighbouring high-J lines of v′ = 35 and 36, R(122) 35-0 and P(132) 36-0, are at 46 and 62 kHz.
   * ΔJ = ±2 mixing is converged: adding J ± 4 changes nothing.
   * [inference] This looks like a local perturbation of B(v′ = 37) near J′ ≈ 141–146, which the smooth formulae cannot represent. Down-weight or exclude these two lines in the refit until it is understood.
   * **Intervals between lines:** −6 to +6 MHz, smooth in v′ (all relative to R(56) 32-0):

     | v′ | Offset |
     |---|---|
     | 32 | within 0.3 MHz |
     | 33 | about −4.7 MHz |
     | 34 | about −5.9 MHz |
     | 35 | about −2.2 MHz |
     | 36 | about +3.9 MHz |
     | 37 | about +5.6 MHz |

   * **Absolute f(a10):** 2.0 MHz high.
   * The inter-line offsets and the f(a10) offset are MHz-level errors in the B levels of the published potentials.
6. **Follow-ups:**
   * **The isotope-shift conventions in item 3 have been checked.** Our code follows Salumbides 2008 as printed, and Salumbides 2006 Fig. 2 shows the same offsets (`hannover-model-reproduction.md`).
   * **The v′ = 58 failure in item 4 is part of a systematic error** in the published B levels above v′ ≈ 44, measured across the Salami & Ross atlas (`spectra-validation.md`).
