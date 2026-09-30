# Absolute I₂ cross sections for validation (visible B←X)

*2026-09-15. Data are in `data/external/saiz_lopez_2004/`, `spietz_2006/` and `mpi_mainz_i2/`. File formats and terms are in `data/external/README.md`, and checksums in `data/external/SHA256SUMS`. Anything marked **inferred** is arithmetic made here or analysis; everything else is quoted from the cited source.*

## 1. Bottom line

* **Measurement conditions.** Two modern data sets give absolute visible cross sections together with enough measurement conditions to simulate them without free parameters:
  * **Spietz et al. 2006** gives two grating spectra, at 0.25 nm and 0.59 nm FWHM, with the I₂ column, path, bath gas and wavelength scale printed in the file headers. It also gives the tightest absolute anchor: σ(500 nm, air) = (2.186 ± 0.021) × 10⁻¹⁸ cm².
  * **Saiz-Lopez et al. 2004** is an FTS spectrum at 4 cm⁻¹ with 760 Torr of air. Its I₂ column (3.1 × 10¹⁶ cm⁻²) is stated in Spietz et al. 2006, not in the Saiz-Lopez paper. Its absolute scale is ±12 %.
* **Wavelength scales.** The Spietz spectra are in vacuum wavelengths, but the Spietz σ(500 nm) is in air. The Saiz-Lopez paper doesn't say which. Cross-correlating against the Salami & Ross FTS atlas (vacuum wavenumbers) puts Saiz-Lopez on a **vacuum** scale to within 0.04 nm, and confirms that the Spietz 0.25 nm spectrum is vacuum. The Spietz 0.59 nm spectrum, however, runs **0.09–0.17 nm long**, which is about one pixel (§4).
* **Tellinghuisen 2011.** The MPI-Mainz file `Tellinghuisen(2011)_294-300K_520-635nm` is **not** a 0.05 nm-step B–X cross section: it has 38 rows at 18 nominal wavelengths. Its values are 7–60 % of the total cross section, so they can't be the total B–X absorption. The papers are behind a Cloudflare challenge, so the cell conditions for all Tellinghuisen 2011 files are unknown. §4a of `02a-broadband-data.md` describes this file inaccurately.
* **Licences.** Both ACP articles are **CC BY-NC-SA 2.5**, not CC BY. The live MPI-Mainz site states no licence; its 2013 Zenodo snapshot is CC BY 4.0. None of these files should ship with the MIT-licensed package.

## 2. Conditions of the absolute data sets

| Data set | T | Range (nm) | Resolution, instrument function | Sampling | λ scale | I₂ column, path | Bath gas | Stated uncertainty |
|---|---|---|---|---|---|---|---|---|
| Saiz-Lopez 2004 FTS | 295 K | 182–750 | 4 cm⁻¹ Bruker IFS/66; apodization not stated. Spietz: "approx. 0.12 nm FWHM at 550 nm" | 1 nm below 500 and above 630 nm; 0.1/0.2 nm for 500–630 nm | not stated; **vacuum** (inferred, §4) | **3.1 × 10¹⁶ cm⁻²** (per Spietz 2006, §5); path not stated | 760 Torr air | ±0.50 × 10⁻¹⁸ at 533 nm (12 %); ±0.27 × 10⁻¹⁸ at 500 nm (Spietz Table 2) |
| Spietz 2006, 0.25 nm (`I2DOASref_1200.TXT`) | room T (MPI: 298 K) | 543.2–578.0 | Czerny–Turner grating spectrometer: 1200 grooves/mm, 170 µm slit → 0.25 nm FWHM (from Hg–Cd lines; shape not given) | 0.035 nm/pixel | vacuum (stated), ±0.04 nm; confirmed | **(6.86 ± 0.20) × 10¹⁵ cm⁻²**, 26.4 cm | 1000 mbar N₂ | larger of 4 % of σ and 7 × 10⁻²⁰ cm² |
| Spietz 2006, 0.59 nm (`I2DOASref_0300.TXT`) | room T (MPI: 298 K) | 428.5–588.3 | same spectrometer: 300 grooves/mm, 50 µm slit → 0.59 nm FWHM | 0.154 nm/pixel | vacuum (stated), ±0.11 nm; **appears 0.09–0.17 nm long** | **(1.42 ± 0.038) × 10¹⁶ cm⁻²**, 26.4 cm | 1000 mbar N₂ | larger of 2.7 % of σ and 7 × 10⁻²⁰ cm² (header); 3 % (paper) |
| Spietz 2006, σ(500 nm) | room T, no stabilization | 500.0 | 150 grooves/mm, 25 µm slit → 0.58 nm FWHM, 0.32 nm/pixel | — | **air** (stated) | I₂ partial pressure < 0.35 mbar, 26.4 cm, so N < 2.2 × 10¹⁷ cm⁻² (inferred, at 298 K) | none (total < 0.7 mbar) | (2.186 ± 0.021) × 10⁻¹⁸ cm², a maximum-error estimate (~1 %) |
| Tellinghuisen 2011, B–X paper (MPI file) | 294–300 K | 520–635 | 0.1 nm (abstract) | "0.05 nm intervals" (MPI); file has 38 rows | not stated | not available | not stated | per-row s.d. |
| Tellinghuisen 2011, C-state paper (MPI files) | 273, 308, 337 K | 390–900 | 2 nm, double-beam spectrophotometer (MPI) | 5 nm; 2 nm in the files with s.d. | not stated | not available | not stated | s.d. in the 308 K and 337 K `400-500,600-850nm` files |
| Tellinghuisen 1973 (MPI file) | 295–300 K | 420–800 | trapezoidal slit, 2.9 nm base and 2.3 nm top (per Spietz); "~2.5 nm" (MPI); "2.6 nm" (Saiz-Lopez) | 10 nm | not stated | extrapolated to zero I₂ pressure, no bath gas (per Saiz-Lopez) | none | σ(500) = (2.20 ± 0.07) × 10⁻¹⁸ (Spietz Table 2) |
| Bauer 1998 (MPI file) | 295 K | 436, 500 | not checked | two points | not stated | not checked | not checked | 436 nm: (1.41 ± 0.05) × 10⁻¹⁹; 500 nm: (2.25 ± 0.09) × 10⁻¹⁸ |

**Pressures implied by the stated columns (inferred, 298 K).** For the 0.25 nm spectrum, 6.86 × 10¹⁵ / 26.4 cm = 2.6 × 10¹⁴ cm⁻³, or 1.1 Pa of I₂. For the 0.59 nm spectrum, 5.4 × 10¹⁴ cm⁻³, or 2.2 Pa. For Saiz-Lopez, the measured vapour density is (7.4 ± 0.7) × 10¹⁵ cm⁻³ (0.225 ± 0.020 Torr at 295 K). With 3.1 × 10¹⁶ cm⁻², that implies a ≈4.2 cm path if the cell was at saturation; neither the paper nor Spietz states the path.

## 3. Data sets

### 3.1 Saiz-Lopez et al. 2004 (`saiz_lopez_2004/`)

**Measurement (paper, §2.1).**
* **Instrument and sample.** Bruker IFS/66 FTS at 4 cm⁻¹ ("0.1 nm at λ = 500 nm"), CaF₂ beamsplitter. Iodine crystals equilibrated in an optical cell at 295 K and 760 Torr of air.
* **Three overlapping regions:**
  * 182–500 nm: GaP diode, D₂ lamp;
  * 260–555 nm: GaP diode, W lamp;
  * 500–1100 nm: Si diode, W lamp.
* **Absolute scale.** The FTS spectrum was "scaled in the continuum region to the average of a series of spectra recorded using a grating spectrometer (Acton SpectraPro SP-556-I, grating 1200 grooves mm⁻¹, resolution 0.2 nm)" at 295 K and 1 atm of air. The I₂ vapour pressure was measured separately as 0.225 ± 0.020 Torr at 295 K.
* **Joining the regions.** According to Spietz 2006 (§5, citing Saiz-Lopez, private communication), the FTS regions were "joined by linear weighting from 500 to 555 nm". Spietz also states that "the column density in the recording of their FTS spectrum was 3.1·10¹⁶ cm⁻²".
* **Saturation.** The paper argues that 760 Torr of air broadens the lines enough to prevent line-centre saturation, and so gives "true extinction values". Spietz compared spectra at 6.9 × 10¹⁵ and 3.1 × 10¹⁶ cm⁻² and found "effects due to different column density are minor".

**Not stated in any source found.** The FTS apodization or instrument-function shape, the FTS cell path length, whether wavelengths are in air or vacuum, and how the FTS wavenumber grid was resampled to the tabulated wavelengths. Spietz notes a "picket-fence effect" at 0.1 nm steps: band heads are levelled off because the step is close to the resolution.

**File.** `Iodine_xsection.xls` has one sheet with 1474 rows from 182 to 750 nm, in cm² molecule⁻¹. The MPI copy is identical to ≤5 × 10⁻⁶ relative.
* **Steps:** 1 nm for 182–500 and 630–750 nm. 500–513.4 nm is 0.1 nm. Above that, points are increasingly skipped: 7 of the 93 steps between 510 and 520 nm are 0.2 nm, rising to 37 of the 63 steps between 620 and 630 nm.
* **Checks against the paper:**
  * The file maximum is 4.238 × 10⁻¹⁸ at 533.0 nm (paper: 4.24 ± 0.50).
  * σ(500.0) = 2.293 × 10⁻¹⁸, which is 4.9 % above Spietz's 2.186.
* **Noise:** above 630 nm the values scatter point to point; 748–750 nm, for example, jumps between 1.6 × 10⁻²⁰ and 4.3 × 10⁻²⁰.

### 3.2 Spietz, Gómez Martín & Burrows 2006 (`spietz_2006/`)

**Spectrometer.** Spietz et al. did **not** run an FTS themselves. All three products come from one setup:
* a 500 mm Acton Research Czerny–Turner spectrograph with a Roper Scientific 1024 × 1024 CCD (26 µm pixels);
* a 150 W Xe arc lamp;
* a (26.4 ± 0.2) cm glass vessel with fused-silica windows.

Two FTS spectra appear in the paper, both from other groups: the Marcy & Butler (1992) atlas and the Saiz-Lopez spectrum.
* **Marcy & Butler (1992):** 0.04 cm⁻¹, room-temperature vapour pressure, 10 cm path, no bath gas. It was used as the "true" spectrum in their simulations (Gaussian instrument function, with a 0.3 cm⁻¹ Lorentzian standing in for atmospheric broadening) and as the wavelength reference.
* **Saiz-Lopez:** used for comparison only.

**DOAS reference spectra (the supplement).** Both were recorded in flow mode at 1000 mbar (N₂ grade 4.8).
* **Wavelength calibration.** Pt–Cr–Ne and Hg–Cd line lamps, refined against the convolved Marcy & Butler spectrum. The fit used 8 band maxima for the 1200 grooves/mm spectrum and 33 for the 300 grooves/mm spectrum.
* **Resolution.** "FWHM was determined from the apparent shape of isolated emission lines from the mercury-cadmium line source." No line shape is published.
* **Vertical axis.** "optical density, corrected for deposit absorption and drift", scaled to absolute σ through σ(500 nm) = 2.191 × 10⁻¹⁸ cm². For the 0.25 nm spectrum, the concentration was transferred from 500 nm measurements taken before and after, with about 4 % uncertainty. So the files are **apparent** cross sections, σ_app = OD/N, at the stated column.
* **Signal to noise.** About 70:1 for the 0.25 nm spectrum and 50:1 for the 0.59 nm spectrum.
* **Range discrepancy.** The paper says the 0.59 nm spectrum covers "445 nm to 600 nm", but the file spans 428.5–588.3 nm.

**σ(500 nm) (paper §6).**
* **Method.** Optical density at 500 nm was plotted against I₂ concentration measured with an MKS Baratron 627B (0.001–1 mbar) and leak-corrected. The I₂ partial pressure was < 0.35 mbar, with no bath gas.
* **Result.** The two series gave 2.194 ± 0.010 and 2.158 ± 0.018 (×10⁻¹⁸ cm²), for a final σ(500 nm) = (2.186 ± 0.021) × 10⁻¹⁸ cm².
* **Weighted mean.** Combined with Tellinghuisen 1973 (2.20 ± 0.07), Saiz-Lopez 2004 (2.29 ± 0.27) and Bauer et al. (2.25 ± 0.09; described in Spietz as "Bauer et al. (2004)", previously unpublished), the weighted mean is **(2.191 ± 0.02) × 10⁻¹⁸ cm² (1σ)**.
* **Wavelength.** "Wavelength with respect to the 500 nm absorption cross section of I₂ is given in air". That is 500.14 nm vacuum, 19 994 cm⁻¹ (inferred). This lies about 46 cm⁻¹ below the B-state dissociation limit (~20 040 cm⁻¹, see `02a-broadband-data.md`).
* **Linearity.** At 0.58 nm FWHM, the OD at 500 nm was linear in concentration, while the OD at 543.8 and 572.2 nm was not (Fig. 8).

**Other copies.** The MPI-Mainz and IUP Bremen copies carry the same numbers. The MPI metadata for the Spietz entries contain typos: "2.9191×10⁻¹⁸" where the files and paper give 2.191, and "0.25 m FWHM".

### 3.3 Tellinghuisen (MPI-Mainz files)

**2011, B–X paper** (J. Chem. Phys. 134, 084301). The abstract (from PubMed) gives 520–640 nm, recorded "with high quantitative precision at moderate resolution (0.1 nm)", and analysed by least-squares spectral simulation into B–X plus continuum.

The MPI file `294-300K_520-635nm` has 38 rows at 18 nominal wavelengths, most of them in pairs. MPI describes it this way: "measured … at 0.05 nm intervals; the first column gives the nominal short wavelength limits, the second column gives the absorption cross sections for the B-X system, the third column the estimated standard deviations". The values are 4.6 × 10⁻¹⁹ at 520 nm and 1.3 × 10⁻¹⁹ at 635 nm.

Compared with the mean Saiz-Lopez σ over [λ, λ + 5 nm], they are 7–17 % of the total at 520–600 nm and ≈60 % at 635 nm (inferred), so they cannot be the total B–X absorption. They may be the continuum components fitted per analysis window. **The meaning is unverified until the paper is read.**

**2011, C-state paper** (J. Chem. Phys. 135, 054301). The abstract says "Absorption spectra are recorded at low resolution but high quantitative precision for I₂ vapor at 35 °C and 64 °C", and that these are analysed together with literature spectra. MPI gives 2 nm resolution and a double-beam spectrophotometer.
* **Files:** 273 K and 308 K at 390–900 nm (5 nm steps, no s.d.); 308 K and 337 K at 400–500 and 600–850 nm (2 nm steps, with s.d.).
* **The 273 K set:** the abstract doesn't mention it. It may be one of the literature spectra (unverified).

**1973** (J. Chem. Phys. 58, 2821). 420–800 nm at 10 nm steps, 295–300 K. The slit function is given above. Per Saiz-Lopez, "made with no added bath gas but extrapolated to zero vapour pressure".

### 3.4 Every I₂ data set in the MPI-Mainz atlas (listing checked 2026-09-15)

The atlas lists 32 files from 15 author–year sources. The 30 that reach the visible are stored in `mpi_mainz_i2/` under their original names. **No Kiefer data set is listed for I₂.** The Kiefer & Bernstein 1972 work catalogued in `02a` is resonance Raman.

"Digitized" means MPI read the values off a figure and converted molar decadic extinction to σ with the factor 3.8235 × 10⁻²¹.

| Author(s), year | T (K) | Range (nm) | Resolution, method (MPI comments) | Stored |
|---|---|---|---|---|
| Bauer et al. 1998 | 295 | 436, 500 | 2 points with error limits; the 436 nm value calibrated an HOI spectrum | yes |
| IUPAC 2007 | 298 | 400–725 | 5 nm averages of Saiz-Lopez 2004 (derived) | yes |
| JPL-2010 (2011) | 295 | 185–700 | 5 nm averages of Saiz-Lopez 2004 (derived) | yes |
| JPL-2010 (2011) | 295 | 500.6–629.8 | maxima and minima picked from Saiz-Lopez 2004 (derived) | yes |
| Kortüm & Friedheim 1947 | 353 | 452.7–589.7 | digitized, Fig. 3 | yes |
| Kortüm & Friedheim 1947 | 613 | 223.7–601.2 | digitized, Fig. 3 | yes |
| Mathieson & Rees 1956 | 393 | 622.5–848.5 | **calculated** from potential curves | yes |
| Mathieson & Rees 1956 | 393 | 634.6–848.5 | measured points, digitized | yes |
| McMillan 1966 | 343–353 | 400–650 | 10 cm cell, Bausch & Lomb 505; pure I₂; digitized | yes |
| McMillan 1966 | 343–353 | 400–650 | as above, I₂ in air | yes |
| Myer & Samson 1970 | 298 | 175.6–202.0 | VUV | no |
| Rabinowitch & Wood 1936 | 293 | 440.8–605.8 | Hilger spectrograph, selenium photocell; digitized | yes |
| Roxlo & Mandl 1980 | 298 | 170–229 | UV | no |
| Saiz-Lopez et al. 2004 | 295 | 182–750 | 4 cm⁻¹ FTS (§3.1) | yes |
| Spietz et al. 2006 | 298 | 428.46–588.29 | 0.59 nm (§3.2) | yes |
| Spietz et al. 2006 | 298 | 500 | weighted mean, 2.191 ± 0.02 | yes |
| Spietz et al. 2006 | 298 | 543.21–577.95 | 0.25 nm (§3.2) | yes |
| Sulzer & Wieland 1952 | 423, 873, 1323 | 357–741 | tabulated on a 250–500 cm⁻¹ grid and converted | yes |
| Sulzer & Wieland 1952 | 573, 723, 1023, 1173 | 363.6–555.6 | as above | yes |
| Tellinghuisen 1973 | 295–300 | 420–800 | ~2.5 nm | yes |
| Tellinghuisen 2011 | 273, 308 | 390–900 | 2 nm, double-beam spectrophotometer | yes |
| Tellinghuisen 2011 | 294–300 | 520–635 | see §3.3 | yes |
| Tellinghuisen 2011 | 308, 337 | 400–500, 600–850 | 2 nm, with s.d. | yes |
| Vogt & Koenigsberger 1923 | 321, 361 | 441–598, 461.6–589.4 | digitized from Rabinowitch & Wood Fig. 1 | yes |

## 4. Wavelength scale check (inferred)

**Method.** The reference is the Salami & Ross transmission atlas (`salami_ross_2005/`), on a vacuum-wavenumber scale calibrated to ±0.003 cm⁻¹, converted to vacuum wavelengths as λ = 10⁷/ν̃. Each data set was compared against it:
1. Convolve the atlas with a Gaussian of the data set's FWHM: 0.25 nm or 0.59 nm for Spietz, and λ² × 4 cm⁻¹ for Saiz-Lopez.
2. Take −ln of the result.
3. High-pass both the reference and the data set by subtracting a 3 nm Gaussian smooth.
4. Find the shift that maximizes the Pearson correlation, searching ±0.40 nm in 0.0025 nm steps.

Each fit was repeated using the convolved optical depth in place of the convolved transmission. The air–vacuum difference uses the Edlén-type dispersion formula for standard air (Morton 2000).

| Data set | Window (nm) | Shift to **add** to the file's λ to reach vacuum (nm): transmission / optical depth | Shift expected if the file were in air |
|---|---|---|---|
| Spietz 0.25 nm | 546–575 | −0.013 / −0.003 | +0.156 |
| Spietz 0.59 nm | 512–572 | −0.155 / −0.143 | +0.151 |
| Spietz 0.59 nm | 576–586 | −0.168 / −0.143 | +0.161 |
| Saiz-Lopez | 505–540 | −0.008 / −0.003 | +0.145 |
| Saiz-Lopez | 540–572 | +0.008 / +0.020 | +0.154 |
| Saiz-Lopez | 577–600 | +0.018 / +0.040 | +0.163 |
| Saiz-Lopez | 600–625 | +0.025 / +0.038 | +0.169 |

Correlation coefficients were 0.77–0.99. Comparing the 0.59 nm spectrum directly with the 0.25 nm spectrum smoothed to 0.59 nm gives −0.090 nm over 547–574 nm. Windows above 575 nm fall in the atlas's 50 °C segment, so hot-band differences may bias those windows slightly.

**Conclusions.**
1. Saiz-Lopez wavelengths are **vacuum**: the offset is ≤0.04 nm, where an air scale would need +0.15 nm.
2. The Spietz 0.25 nm scale is correct.
3. The Spietz 0.59 nm scale is too long by 0.09–0.17 nm, about one 0.154 nm pixel. The offset has the opposite sign to an air/vacuum mix-up, and at the upper end it exceeds the stated ±0.11 nm.

Also as a consistency check, the mean σ of the 0.59 nm spectrum over 547–574 nm is 5.4 % above that of the 0.25 nm spectrum (inferred). That is within their combined 3 % and 4 % uncertainties.

## 5. What a forward simulation needs

**Forward model.** For every data set the quantity to compute is the apparent cross section, with the logarithm taken *after* instrumental smoothing:

σ_app(λᵢ) = −(1/N) ln ∫ Pᵢ(λ) exp[−N σ(λ; T, p)] dλ

Here Pᵢ is the normalized response of sample i, and σ contains two parts:
* **Lines:** the line-by-line B←X spectrum, with Voigt profiles made of Doppler broadening at T plus pressure broadening by the bath gas;
* **Continua:** A←X, C←X and B←X bound–free (see `intensity-inputs.md`; all of these data sets include the continua).

The line-by-line σ is Doppler-limited, so it has to be computed on a fine grid before the integral is taken.

| Input | Saiz-Lopez 2004 | Spietz 0.25 nm | Spietz 0.59 nm | Spietz σ(500) |
|---|---|---|---|---|
| T | 295 K | ~298 K (not in file) | ~298 K | room T |
| N (cm⁻²) | 3.1 × 10¹⁶ (Spietz, private communication) | 6.86 × 10¹⁵ | 1.42 × 10¹⁶ | ≤2.2 × 10¹⁷; linear regime shown |
| Bath gas | 760 Torr air | 1000 mbar N₂ | 1000 mbar N₂ | none (Doppler only) |
| Pᵢ | FTS instrument function at 4 cm⁻¹ in ν̃; apodization unknown (bracket boxcar-sinc against a strong apodization); then point-sample at the tabulated λ | 0.25 nm FWHM (Gaussian assumed; shape unpublished) ⊗ 0.035 nm pixel boxcar | 0.59 nm FWHM ⊗ 0.154 nm pixel boxcar (pixel = 26 % of the FWHM) | 0.58 nm FWHM at 500 nm |
| λ scale | vacuum (inferred) | vacuum | vacuum after a −0.09 to −0.17 nm correction | air → 500.14 nm vacuum |
| Scale uncertainty | ±12 % (vapour pressure ±9 %) | ±4 % | ±3 % | ±1 % |

**Broadening inputs.** The Spietz simulations used a 0.3 cm⁻¹ Lorentzian for 1 atm; the width convention wasn't stated. Measured broadening coefficients for air near 543 nm (Fletcher & McDaniel 1995) and N₂ at 675 nm (Wolf 2009) are catalogued in `02a-broadband-data.md`.

Note that the Spietz instrument FWHM was measured on the same CCD, so it already includes the pixel width to first order. The extra boxcar matters mainly for the 0.59 nm spectrum.

## 6. Recommended comparisons (in order of discriminating power)

1. **σ(500 nm air = 500.14 nm vacuum) against 2.186 ± 0.021 × 10⁻¹⁸ cm² (Spietz).** This is the tightest absolute test, and saturation plays no part in it. It probes the near-dissociation B levels and the continua, not isolated bands. The Saiz-Lopez value at 500 nm (2.29 ± 0.27) is a looser check.
2. **Spietz 0.25 nm, 546–575 nm, at N = 6.86 × 10¹⁵ cm⁻² and 1000 mbar N₂.** Compare:
   * the window-mean σ_app, an absolute test at ±4 %;
   * the differential band structure (band-head positions and contrast), which doesn't depend on the absolute scale;
   * σ_app at the stated N against σ_app as N → 0, to measure how much saturation a pyodine simulation predicts at this resolution.
3. **Spietz 0.59 nm, 445–588 nm, at N = 1.42 × 10¹⁶ cm⁻².** Determine the wavelength shift against the model and report it (a λ-only nuisance parameter, expected at about −0.1 to −0.15 nm), then compare absolute σ_app at ±3 %. Treat values below ~485 nm with caution: in the σ(500) experiment, window deposit contributed 3 % at 485 nm and 20 % at 460 nm. The reference spectra are stated to be deposit-corrected.
4. **Saiz-Lopez, 500–630 nm, at 4 cm⁻¹, N = 3.1 × 10¹⁶ cm⁻², 760 Torr air and 295 K.** This is the highest-resolution absolute data set. Show it both as published and rescaled by the fixed factor 2.186/2.293 = 0.953, which puts it on the Spietz σ(500). Neither version involves a fitted parameter. Expect band-correlated differences: Spietz found ratios of up to 10 % between their grating spectra and this FTS spectrum, with a jump near 560 nm that they attributed to resolution, step size, drift or deposit.
5. **The red continuum and temperature dependence.** Use Saiz-Lopez at 630–750 nm and Tellinghuisen 2011 at 308 and 337 K (600–850 nm, with s.d.). These can only be compared semi-quantitatively until the Tellinghuisen cell conditions are known.
6. **Coverage checks only.** Sulzer & Wieland (423–1323 K) and the other digitized historical sets provide sanity checks at hot-cell temperatures. They are low-resolution, digitized and have no stated conditions.

## 7. Caveats and gaps

* **No rotationally resolved absolute data.** Every data set here is at ≥0.1 nm resolution, so they test band-integrated intensities, the continua and saturation handling, not individual line strengths.
* **Saiz-Lopez conditions are second-hand.** The I₂ column comes from Spietz, not from the Saiz-Lopez paper. The path, apodization and resampling method are unknown, and the air/vacuum scale is inferred.
* **Air vs vacuum is mixed within Spietz 2006.** The spectra are vacuum and σ(500) is air. The historical, Tellinghuisen and Bauer scales are unknown, though at 2–10 nm resolution this rarely matters.
* **Derived sets are not independent.** The JPL-2010 and IUPAC-2007 files are derived from Saiz-Lopez, so they add no independent information. JPL 19-5 was not checked.
* **Not obtained:**
  * the Tellinghuisen 2011 papers (pubs.aip.org Cloudflare challenge);
  * the Bauer 1998 measurement conditions (ACS not consulted);
  * the MPI-Mainz 2013 Zenodo snapshot (HTTP 403, "unusual traffic"), so it is unverified whether the live files match the CC BY 4.0 snapshot.
* **Search gap.** The web-search budget ran out, so no further search was made for other visible absolute data sets or for a Kiefer data set.
* **Licences.** ACP 2004 and 2006 are CC BY-NC-SA 2.5. That is NonCommercial and ShareAlike, so the files must not be bundled with the MIT-licensed package. Download them on demand, as `data/external/README.md` already requires.

## References

* A. Saiz-Lopez et al., Atmos. Chem. Phys. 4, 1443–1450 (2004), doi:10.5194/acp-4-1443-2004. Supplement: https://acp.copernicus.org/articles/4/1443/2004/acp-4-1443-2004-supplement.zip
* P. Spietz, J. C. Gómez Martín, J. P. Burrows, Atmos. Chem. Phys. 6, 2177–2191 (2006), doi:10.5194/acp-6-2177-2006. Supplement: https://acp.copernicus.org/articles/6/2177/2006/acp-6-2177-2006-supplement.zip. Bremen copies: https://www.iup.uni-bremen.de/gruppen/molspec/databases/referencespectra/i2spectra/index.html
* H. Keller-Rudek et al., Earth Syst. Sci. Data 5, 365–373 (2013), doi:10.5194/essd-5-365-2013. I₂ page: https://uv-vis-spectral-atlas-mainz.org/uvvis/cross_sections/Halogens+mixed%20halogens/I2.spc. Snapshot: doi:10.5281/zenodo.6951 (CC BY 4.0).
* J. Tellinghuisen, J. Chem. Phys. 134, 084301 (2011), doi:10.1063/1.3555623; J. Chem. Phys. 135, 054301 (2011), doi:10.1063/1.3616039; J. Chem. Phys. 58, 2821 (1973), doi:10.1063/1.1679584.
* D. Bauer et al., J. Phys. Chem. A 102, 2857–2864 (1998), doi:10.1021/jp9804300.
* H. Salami, A. J. Ross, J. Mol. Spectrosc. 233, 157–159 (2005), doi:10.1016/j.jms.2005.06.002 (vacuum wavelength reference for §4).
