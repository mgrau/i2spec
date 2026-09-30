# Inputs for intensities and spectra (M3)

*Read from the full texts, 2026-09-15.*

## Line strengths: Tellinghuisen, J. Chem. Phys. 134, 084301 (2011)

**Absorption formula** (eqs. 1–3):

ε_ν = 108.862 · S_L(ν) · ν · s_J′J″/(2J″+1) · |⟨v′J′|μ_e|v″J″⟩|² · f(v″, J″)

| Quantity | Meaning and units |
|---|---|
| ε_ν | Molar absorptivity, L mol⁻¹ cm⁻¹ |
| S_L(ν) | Unit-area line profile, cm |
| ν | Wavenumber, cm⁻¹ |
| μ_e | Transition moment, D |
| s_J′J″ | Hönl–London factor: J″ for P lines, J″ + 1 for R lines |
| f(v″, J″) | Fractional population of the lower level, including nuclear-spin weights (7/6 for odd J″ and 5/6 for even J″, relative to no symmetry, i.e. 21:15) |

**Absorbance and cross section.**
* Absorbance (base 10): A = ε l C, with l in cm and C in mol/L.
* Cross section: σ = 3.8235×10⁻²¹ ε cm². Hence the integrated line cross section is 4.1624×10⁻¹⁹ · ν · s/(2J″+1) · |μ|² · f, in cm.

**Transition moment.**
* Linear fit, 520–640 nm (R ≈ 2.63–2.88 Å), eq. (9): μ_e(R) = 1.1123(32) D + 0.712(24) D/Å · (R − 2.85 Å). Relative standard errors are about 1%.
* Recommended global function, eq. (10): μ_e²(R) = c₁ exp[−c₂(R − c₃)²] R⁻², with c₁ = 21.5(1.9) D² Å², c₂ = 1.18(12) Å⁻², c₃ = 3.65(9) Å.
  * It fits both the absorption data and the radiative rates A_T (Vigué et al. 1981) over R ≈ 2–6 Å.
  * Check values: 0.956 D² at 2.66 Å, and 1.71 D² at 3.3 Å (near the maximum).
* R-centroid for v″ = 0–3: R̄(Å) = 2.6322 + 0.001548(λ − 500) − 9.1×10⁻⁷(λ − 500)², with λ in nm. The R-centroid approximation is only a check; i2spec computes the full matrix elements.

**Line profile and cell conditions.**
* Doppler FWHM (I₂): Δν = 4.494×10⁻⁸ · ν · √T cm⁻¹.
* Vapour pressure, eq. (8): log₁₀(P/Torr) = 12.1891 − 0.001301 T − 0.3523 log₁₀T − 3410.71/T. This is 0.2% above NIST-JANAF 1998, and gives 0.201 Torr at 20 °C.

**Continuum.** Molar absorptivity ε_c (Table II), in L mol⁻¹ cm⁻¹. The cell body was at 34–38 °C and the side arm, which sets the pressure, at 21–27 °C, so these are about 309 K values, not room temperature. The physical continuum model that replaces this table is in `continuum-model.md` and `src/i2spec/continuum.py`.

| λ (nm) | ε_c |
|---|---|
| 520 | 120–133 |
| 530 | 85 |
| 540 | 58–64 |
| 555 | 38–41 |
| 570 | 26–28 |
| 580–600 | 19–23 |
| 610 | 23–27 |
| 625 | 31 |
| 633–635 | 32–35 |

This is the continuum of the A←X and C←X transitions plus the B←X continuum from excited v″. The C←X strength is 25% below the 1982 estimate.

## Transition moment over a wide range: Lamrini et al., J. Chem. Phys. 100, 8780 (1994)

Relative |μ_e|² from laser-induced fluorescence, for R_c = 2.633–6.035 Å:

| R_c range | Fit | Scatter |
|---|---|---|
| 2.633–3.9 Å | log₁₀|μ|² = −134.0416 + 166.1989 R − 75.81764 R² + 15.12918 R³ − 1.12537 R⁴ | 15% |
| 3.9–6.035 Å | log₁₀|μ|² = 2.6642 − 1.1635 R | 55% |

As written, the fit gives |μ|² = 1.00 D² at 2.66 Å. Their final, lifetime-based calibration adds +0.30103 to the constant, giving |μ|²(2.66 Å) = 2.00(11) D². That is twice the absorption-based value. Tellinghuisen 2011 eq. (10) reconciles the absorption data with the lifetimes, so **i2spec uses eq. (10)** and keeps Lamrini's fit only as a check of the shape.

## Validation spectrum: Salami & Ross, J. Mol. Spectrosc. 233, 157 (2005)

**Instrument.**
* Bomem DA3 FTS, 14 250–20 100 cm⁻¹, instrumental resolution 0.02 cm⁻¹.
* Observed linewidths 0.035–0.05 cm⁻¹ (unresolved hyperfine structure).

**Cell.**
* 50 cm long, with a 6 cm sidearm.
* **The sidearm temperature, which sets the I₂ pressure, is not given**, so the column density has to be fitted for each segment.

**Segments** (Table 1; temperature of the main cell body):

| Range (cm⁻¹) | Filters | Cell body |
|---|---|---|
| 19 100–20 500 | CS 4-97 + CS 5-56 | 20 °C |
| 16 500–19 500 | CS 4-97 | 20 °C |
| 15 400–17 400 | CS 3-67 + LS 650 | 50 °C |
| 14 500–16 500 | CS 2-61 + LS 700 | 190 °C |

**Calibration.**
* The segments were spliced and averaged, so the absolute transmission level is not reliable. The baseline is about 80% in their Fig. 1.
* Frequencies were calibrated against Katô and Gerstenkorn–Luc, with a correction of about 4×10⁻⁷ σ.
* Stated accuracy: ±0.003 cm⁻¹, except at 15 797.7–15 798.5 cm⁻¹ (HeNe spike).
