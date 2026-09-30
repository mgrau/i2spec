# Stage 2a: Broadband, Doppler-limited and supporting data for I₂

*pyodine research notes, compiled September 2026. Companion catalog: [`data/catalog/broadband.yaml`](../../data/catalog/broadband.yaml).*

This document covers the non-sub-Doppler I₂ data a global B–X model needs:

1. Fourier-transform (FTS) absorption atlases and other broadband spectra.
2. Fluorescence and term-value data, and the other electronic states that matter for B–X.
3. Data on the isotopologues ¹²⁹I₂ and ¹²⁷I¹²⁹I.
4. Absolute intensities, the transition moment, lifetimes and predissociation, and vapor pressure.
5. Pressure broadening and shift.

Sub-Doppler, hyperfine-resolved and comb-referenced absolute frequencies of ¹²⁷I₂ B–X lines belong to the companion Stage 2b catalog. They appear here only where they recalibrate the atlases or anchor term values.

## 0. Method, conventions and summary

**Search method.**
- **Sources.** Crossref, OpenAlex, Semantic Scholar, HAL, DataCite, arXiv, NCBI and publisher pages, plus general web search. The full text was also read of Knöckel, Bodermann & Tiemann, *Eur. Phys. J. D* 28, 199 (2004).
- **Forward citations (OpenAlex).**

  | Starting work | Citing works found |
  |---|---|
  | Salami & Ross 2005 | 94 |
  | Knöckel et al. 2004 | 87 |
  | Gerstenkorn & Luc 1979 (atlas correction) | 507 |
  | Gerstenkorn & Luc 1985 (Dunham description of the atlases) | 179 |
  | Martin et al. 1986 (X-state LIF-FTS) | ≈130 |
  | Salumbides et al. 2006 and 2008 | 24 and 31 |
  | Selected line-shape papers | – |

  The Orsay atlases are books and are not indexed as citable works in OpenAlex or Crossref. Their 1979 calibration paper and 1985 Dunham paper therefore served as proxies.
- **Limits reached.** The OpenAlex daily budget and the web-search budget ran out late in the work. A few forward-citation searches in Category 5 were not completed; they are listed in §5.

**Integrity conventions.**
- Every catalog entry was checked against at least one of: a DOI record (Crossref), an OpenAlex, HAL, ADS or DataCite record, a library/ISBN record, a publisher or arXiv page, or the reference list of a peer-reviewed paper we read. All entries have `verified: true`.
- Numbers seen only in a search-engine rendering of an abstract are marked as such in `notes`.
- Our own arithmetic and inferences are labelled "computed here" or "inferred".
- Unknown quantities are `null`.
- `uncertainty_MHz` is filled only when the source states a 1σ (or standard-error) value. Conversions use 1 cm⁻¹ = 29 979.2458 MHz.

**Catalog content.**

| Category | Entries | High | Medium | Low | Published 2004 or later | Machine-readable | Likely in Knöckel 2004 |
|---|---|---|---|---|---|---|---|
| 1. FTS atlases & broadband | 33 | 9 | 9 | 15 | 18 | 2 | 7 |
| 2. fluorescence/term values/other states | 29 | 5 | 10 | 14 | 3 | 0 | 3 |
| 3. isotopologues | 20 | 4 | 1 | 15 | 5 | 0 | 1 |
| 4. intensity/lifetimes/vapor pressure | 56 | 9 | 18 | 29 | 11 | 17 | 0 |
| 5. line shape | 21 | 2 | 11 | 8 | 12 | 0 | 0 |
| **Total** | **159** | 29 | 49 | 81 | 49 | 19 | 11 |

**The most valuable items for the global fit:**

1. **Knöckel et al. 2004 supplementary table** (`knockel2004a`). This is the Hannover fit data set: 1513 assigned rovibronic frequencies, v″ ≤ 17, v′ ≤ 43, J ≲ 238. We could not retrieve it. The paper points to a supplementary table at edpsciences.org, but the current Springer page carries no electronic supplementary material and the Internet Archive was offline. Recovering it (for example from E. Tiemann, or from B. Bodermann, now at PTB) is the single most useful action.
2. **Knöckel 2004 Tables 2–3** (`knockel2004b`). These recalibrate the Gerstenkorn atlases to about 20–30 MHz (1σ).
3. **Salami & Ross 2005** (`salami2005`). A freely downloadable 1.16-million-point ASCII transmission spectrum, 14 250–20 050 cm⁻¹. It is the only machine-readable broadband B–X spectrum.
4. **Nölleke et al. 2018** (`nolleke2018`). About 10 000 lines at 915–985 nm, with 50 MHz absolute accuracy, as CSV on Mendeley Data. This is new territory beyond all the atlases.
5. **Reiners et al. 2024** (`reiners2024`). A comb-calibrated FTS spectrum of an iodine cell over 5150–6300 Å. It reveals a systematic ≈2 m s⁻¹ (≈3.6 MHz) band-correlated error in the Knöckel model. E. Tiemann is a co-author.
6. **Orsay atlases** (`gerstenkorn1978a`, `gerstenkorn1978b`, `gerstenkorn1982`, with the correction `gerstenkorn1979a`). Still the backbone for 11 000–20 000 cm⁻¹, but printed only.
7. **Supporting categories.** The Lyon LIF-FTS X-state work (`martin1986`, `bacis1986`, `cerny1986`), the Novosibirsk sub-MHz X-state anchors (`matyugin2012`, `nesterenko2019`), the isotopologue data of Salumbides et al. (`salumbides2006`, `salumbides2008`), and the intensity references of Tellinghuisen 2011 and Saiz-Lopez 2004 (via the MPI-Mainz atlas).

**A naming note for the project.** "Pyodine" is already the name of an open-source Python package for iodine-cell radial velocities: Heeren et al., *A&A* 674, A164 (2023), doi:10.1051/0004-6361/202244441, code at gitlab.com/Heeren/pyodine. The PyPI name `pyodine` returned no project when queried (September 2026), but the name clash within the iodine-cell community is worth considering.

---

## 1. Fourier-transform absorption atlases and broadband spectra

### 1.1 The Orsay FTS atlases (Laboratoire Aimé Cotton)

The Gerstenkorn–Luc (GL) atlases, recorded with the Laboratoire Aimé Cotton FTS in the late 1970s, remain the only complete, calibrated survey of Doppler-limited B–X absorption from 11 000 to 20 000 cm⁻¹. Gerstenkorn & Luc (1985) state that the atlases hold about 100 000 recorded transitions over 11 000–20 040 cm⁻¹. These are reproduced by 46 Dunham-type constants with a standard error of 0.002 cm⁻¹, covering X v″ ≤ 19 and B v′ ≤ 80.

Recording conditions below are from Knöckel 2004 Table 3; the corrections and 1σ values are from Knöckel 2004 Table 2.

| Catalog id | Range (cm⁻¹) | Publication | Knöckel label | Cell temperature (°C) | Knöckel correction | Estimated 1σ even / odd J″ (10⁻³ cm⁻¹) |
|---|---|---|---|---|---|---|
| `gerstenkorn1982` | 11 000–14 000 | Gerstenkorn, Vergès & Chevillard, LAC Orsay 1982 | atlas I, part 1 (11 000–13 000) | 790 | none | 0.65 / 1.1 |
| | | | atlas I, parts 2–4 (13 000–14 000) | 500 | −2.5×10⁻³ cm⁻¹ | 1.0 / 1.5 |
| `gerstenkorn1978b` | 14 000–15 600 | Gerstenkorn & Luc, LAC Orsay 1978 | atlas II (14 000–15 800) | 250 | ν(1 + 1.7×10⁻⁷) − 8.1×10⁻³ | 0.71 / 0.99 |
| `gerstenkorn1978a` | 14 800–20 000 | Gerstenkorn & Luc, Éditions du CNRS 1978 (546 pp.); issued at LAC as 15 600–17 600 and 17 500–20 000 (1977) | atlas III (14 800–20 000) | 250 | same | 0.71 / 0.99 |
| `gerstenkorn1980b` | 14 800–20 000 | *Complément: identification des transitions du système (B–X)*, Éditions du CNRS 1980 | – | – | – | assignments |

The cell temperatures are printed in a column headed °C. The 790 figure for part 1 is reproduced as printed.

**Status and use.**
- All volumes are printed only. No official electronic line list was found.
- The atlas gives wavenumbers of line maxima of hyperfine-blended Doppler profiles. For odd J″ these differ from the rovibronic (hyperfine-free) frequency by −6 to −8 × 10⁻⁴ cm⁻¹ (Knöckel Table 3).
- Knöckel et al. used atlas lines only to fill gaps: 514–526 nm and 667–776 nm, the latter mostly for 8 ≤ v″ ≤ 11. They also used them to link bands in their near-IR (778–815 nm) local fits.
- The heated-cell 11 000–14 000 cm⁻¹ atlas is the only broadband source for hot bands up to at least v″ ≈ 14. Its band analysis is in Gerstenkorn, Luc & Sinzelle (1980; `gerstenkorn1980a`).
- A–X lines also fall in this region. They were analysed from Orsay FTS data by Appadoo et al. (1996; see §2) and must be modeled or masked.

**Related Orsay papers** (all catalogued, mostly low priority):
- 5350 Å band analysis (1977; `gerstenkorn1977`).
- Assignments of line groups (1979; `gerstenkorn1979b`).
- Calibration and assignments near 5915 Å (1981; `gerstenkorn1981a`).
- A laser-excitation "Complement to the Atlas" at 15 780–15 815 cm⁻¹, with more than 750 transitions, J = 4–167 (1981; `gerstenkorn1981b`).
- Early Dunham constants (Luc 1980; `luc1980`).
- The atlas-wide Dunham description (1985; `gerstenkorn1985a`).

Several of these papers are open access on HAL.

**Not found.** We did not confirm the existence of any Orsay volume below 11 000 cm⁻¹. Nölleke et al. (2018) compare their 915–985 nm lines with 3769 "reference" lines, but they cite the 14 800–20 000 cm⁻¹ atlas, which cannot cover that range. The identity of their reference list is therefore open (see §1.4).

### 1.2 Recalibrations of the Gerstenkorn–Luc wavenumber scale

1. **Gerstenkorn & Luc (1979), Rev. Phys. Appl. 14, 791** (`gerstenkorn1979a`). The published correction was derived by the atlas authors themselves, using an internal-standard method: *subtract 0.0056 cm⁻¹ from all wavenumbers of the 14 800–20 000 cm⁻¹ atlas*, giving an absolute accuracy of ±0.002 cm⁻¹ (≈ ±60 MHz; confidence level not stated). This is the correction usually meant by "the GL correction".
2. **Knöckel, Bodermann & Tiemann (2004), Tables 2–3** (`knockel2004b`). They simulated each atlas line profile: the interpolated hyperfine structure, convolved with the Doppler profile at the recording temperature, the rectangular aperture function and the instrumental sinc. They compared the results with high-precision frequencies and derived the part-by-part rules in the table above, plus the odd-J″ maximum-to-rovibronic shifts. Knöckel also refers to a frequency-dependent calibration effect across atlases II and III "already mentioned by Palmer et al." (Los Alamos informal report LA-8251-MS, 1980). That report is not an I₂ data set and is not catalogued.
   - For the 14 000–20 000 cm⁻¹ volumes the rule is ν_corr = ν_atlas(1 + 1.7×10⁻⁷) − 8.1×10⁻³ cm⁻¹.
   - This gives −5.7, −5.2 and −4.7 × 10⁻³ cm⁻¹ at 14 000, 17 000 and 20 000 cm⁻¹ (arithmetic made here), bracketing the original −5.6 × 10⁻³ cm⁻¹.
   - After correction, unblended lines carry an estimated 1σ of 0.7–1.5 × 10⁻³ cm⁻¹ (20–45 MHz).
3. **Other checks of the atlas scale.** Rakowsky, Zimmermann & Ernst, *Appl. Phys. B* 48, 463 (1989), doi:10.1007/BF00694680, measured accurate wavenumbers of iodine lines in the red. The Amsterdam 1 MHz grids of Velchev et al. 1998 (571–596 nm, doi:10.1006/jmsp.1997.7480) and Xu et al. 2000 (595–655 nm, doi:10.1006/jmsp.2000.8085) supplied part of the precise reference data against which Knöckel recalibrated the atlas. These are sub-Doppler and belong to Stage 2b.
   - The Salami & Ross spectrum was, per its abstract as rendered by a search engine, put on a scale matched to other calibrated atlases.
   - Nölleke et al. report "excellent agreement" with their reference list at 915–985 nm, whose accuracy they quote as 0.45 pm (≈150 MHz).
4. **Test of the model rather than the atlas.** Reiners et al. 2024 (§1.6) tested the Knöckel model with a comb-calibrated FTS. Their residual pattern is the first external check at the MHz level across 515–630 nm.

**Recommendation.** Any atlas line used in pyodine should carry a provenance flag for its correction (none, GL1979, or K2004), plus the K2004 odd-J″ shift.

### 1.3 The Salami & Ross (2005) digital atlas

Salami & Ross, *J. Mol. Spectrosc.* 233, 157 (2005) (`salami2005`), published a Doppler-limited FTS transmission spectrum of I₂ as a single ASCII file. The Elsevier supplementary file `mmc1.txt` is still freely downloadable. We downloaded it in September 2026:

- **Header:** "Transmission spectrum of iodine, B-X system."
- **Columns:** two — wavenumber (cm⁻¹) and transmission (%).
- **Grid:** 1 160 001 points on a uniform 0.005 cm⁻¹ grid from 14 250.100 to 20 050.100 cm⁻¹ (498.8–701.8 nm), with no gaps.
- **Transmission:** 3.1 % to 93.9 %.

The file does not record cell temperature, path length, pressure or calibration. According to the abstract as rendered by a search engine (the full text was not accessible here), the spectrum was taken at 0.02 cm⁻¹ instrumental resolution, is reliable to ±0.003 cm⁻¹, and has its wavenumber scale matched to other calibrated atlases.

Its value for pyodine:
1. It is the only machine-readable broadband B–X absorption spectrum covering most of the GL range.
2. With cell conditions from the paper, it is an end-to-end test of simulated transmission (positions, intensities, Doppler and instrumental widths).
3. Profile fitting against a model could re-derive line positions at the ~10⁻³ cm⁻¹ level without digitizing the printed atlas.

It is heavily cited (94 citing works) as a calibration reference, for example in dual-comb and laser-stabilization papers.

### 1.4 Other laboratory Doppler-limited broadband data (post-2003)

- **Nölleke, Raab, Neuhaus & Falke (TOPTICA), J. Mol. Spectrosc. 346, 19 (2018)** (`nolleke2018`).
  - **Method:** a mode-hop-free diode laser, wavemeter-stabilized, with the wavemeter recalibrated every 5 min against a Cs D1-locked laser, scanned 915–985 nm (10 152–10 929 cm⁻¹).
  - **Cell:** 75 cm, triple-passed (225 cm effective); body at 300 °C, cold finger at 39 °C (127 Pa).
  - **Result:** 9970 lines at about 50 MHz absolute (20 MHz frequency, 33 MHz step and 30 MHz line-centre contributions). About 6000 of them are not in the earlier reference list.
  - **Data:** on Mendeley Data (doi:10.17632/cvr6vv2zps.1, CC BY 4.0): a line list (`iodine_atlas.csv`, 10 162 rows; 3664 of them with a matched reference wavelength, as counted here) and the raw transmission scan (`scan_data.csv`, 6 MB).
  - **Why it matters:** lines at these wavelengths come from high v″ hot bands, so they constrain the X potential well above v″ = 17 (inferred; lines are unassigned).
- **Marcassa group (São Carlos), J. Mol. Spectrosc. 387, 111668 (2022) and 395, 111789 (2023)** (`torres2022`, `rodriguezfernandez2023`). Doppler-limited laser absorption at 14 600–14 710 and 14 400–14 600 cm⁻¹. The abstract of the 2022 paper, as rendered by a search engine, reports 648 lines: 410 agreeing with earlier data and 238 new. Calibration and uncertainties could not be read (paywall). This is the 667–776 nm zone where Knöckel's model is weakest (2σ up to 60 MHz), so these papers deserve a full-text check.
- **Williamson, Kuntzleman & Kafader, J. Chem. Educ. 90, 383 (2013)** (`williamson2013`). An educational B–X data set reported to contain 7381 lines. Not verified; probably derivative.

### 1.5 FTS scans of astronomical iodine cells

Iodine cells are the classic wavelength reference for precise radial velocities. Every such cell needs an ultra-high-resolution FTS "template". These scans are a large, mostly unpublished body of Doppler-limited I₂ spectra at 5000–6300 Å, recorded at known cell temperatures.

| Cell / programme | FTS scans (what we verified) | Availability | Catalog id |
|---|---|---|---|
| Lick Hamilton, Keck HIRES | Kitt Peak McMath 1-m FTS templates (Valenti et al. 1995; Butler et al. 1996); HIRES cell rescanned with the NIST FTS in March 2009 at R ≥ 500 000 (Keck web page) | "available upon request" (Keck); McMath spectra largely public in the NSO archive, but iodine files not identified | `valenti1995`, `butler1996` |
| HET HRS | KPNO FTS (≈ two decades old) and NIST FTS on 15 Nov 2011 at 65, 70 and 75 °C; compared with Tull TS12 échelle spectra using IodineSpec5 | not public | `wang2019` |
| Four ESPRESSO cells (AS-31 to AS-34) | NIST 2-m FTS, March 2023, 14 200–21 300 cm⁻¹ at 0.018, 0.010 and 0.007 cm⁻¹, S/N 600–1700 | not stated; ESPRESSO data in ESO archive | `nave2026` |
| ESPRESSO validation | forward model uses Knöckel + comb-calibrated FTS; model "reliable only for λ ≳ 5150 Å" | ESO archive | `schmidt2025` |
| IAG Göttingen FTS (Bruker IFS 125HR) | 632 spectra of a Hamburg cell at 23–66 °C (2015); iodine + Fabry–Pérot calibration of the solar FTS; comb-referenced iodine spectrum (2023) | not public | `perdelwitz2018`, `debus2023`, `reiners2024` |
| SONG (Tenerife), Lick (reanalysis) | I₂ templates from FTS scans used by the "pyodine" RV code (Heeren et al. 2023); scan details not extracted | via the pyodine RV project (not checked) | – |
| Subaru HDS, Okayama HIDES | cells described by Kambe et al., PASJ 54, 865 (2002), doi:10.1093/pasj/54.6.865; Sato et al., PASJ 54, 873 (2002), doi:10.1093/pasj/54.6.873; Kambe et al., PASJ 60, 45 (2008), doi:10.1093/pasj/60.1.45; FTS details not extracted | – | – |
| Magellan PFS, Lick APF, CTIO CHIRON | instrument papers: Crane et al., SPIE 6269 (2006), 7014 (2008), 7735 (2010); Radovan et al., SPIE 7735 (2010), 9145 (2014); Tokovinin et al., PASP 125, 1336 (2013), doi:10.1086/674012; no FTS template details verified | – | – |

**Assessment.** The astronomical scans are high-S/N, high-resolution and absolutely calibrated (FTS), at controlled temperatures. They are therefore excellent validation data for a line-by-line intensity and line-shape model at 500–630 nm, and a potential source of new line positions (the NIST 0.007 cm⁻¹ scans). Almost none are public. The NIST scans (G. Nave), the IAG spectra (A. Reiners) and the Keck 2009 scan are the obvious requests.

The NSO FTS archive (https://nispdata.nso.edu/ftp/FTS_cdrom/) is organized only by date (FTS01–FTS57). Its README directs users to the NSO Digital Library query tool, so we could not confirm whether any iodine-cell scans are among the public files.

### 1.6 Comb-based and dual-comb broadband spectra

- **Reiners et al., A&A 690, A210 (2024)** (`reiners2024`).
  - **Setup:** a Bruker IFS 125HR (max OPD 136 cm, no apodization) recorded an iodine cell (5200–6200 Å) and a 1 GHz laser frequency comb (8000–8400 Å) simultaneously. 17 807 comb lines fixed the FTS zero point to about 1 cm s⁻¹ (≈ 0.02 MHz).
  - **Comparison with the model:** they compared 19 spectra with IodineSpec-type synthetic spectra from the Knöckel model in 302 chunks over 5150–6300 Å. They found a long-period residual wave of ≈2 m s⁻¹ amplitude (≈3.6 MHz, conversion made here) plus an oscillation that follows the band structure.
  - **Interpretation:** the authors attribute this to the energy scale (potentials and hyperfine), not intensities. They provide a spline correction good to 0.5–1 m s⁻¹ over 5300–6150 Å.
  - **Why it matters:** this is the best modern evidence of residual systematic error in the Hannover model in its best-covered region.
- **Ideguchi et al., Opt. Lett. 37, 4847 (2012)** (`ideguchi2012`, with CLEO 2012 CW1J.3). Visible dual-comb spectrum near 19 240 cm⁻¹ (400 cm⁻¹ span) at Doppler-limited resolution in 12 ms. A demonstration without an absolute frequency scale.
- **Eber et al., Opt. Express 33, 35314 (2025) and Laser Photonics Rev. 20, e02713 (2026)** (`eber2025`, `eber2026`). Comb-resolved visible dual-comb spectra of I₂ agreeing with the literature over 10 THz. The 80 MHz-repetition-rate system shows I₂ hyperfine structure.
- **Nishiyama, Ishikawa & Misono, JOSA B 30, 2107 (2013).** A comb-assisted system over 570–612 nm with about 100 kHz hyperfine accuracy. It is sub-Doppler and belongs to Stage 2b.

No comb-referenced broadband data set with released line positions exists yet. The technology (visible dual-comb, comb-calibrated FTS) is mature enough that a dedicated campaign over 500–700 nm and on heated cells is realistic.

### 1.7 Cavity-enhanced broadband spectra

The cavity-enhanced literature on I₂ is incoherent broadband cavity-enhanced absorption spectroscopy (IBBCEAS) for atmospheric monitoring:
- Vaughan et al., PCCP 10, 4471 (2008) (`vaughan2008`).
- Dixneuf et al., ACP 9, 823 (2009), doi:10.5194/acp-9-823-2009.
- Johansson et al., Appl. Phys. B 114, 421 (2014) (`johansson2013`).

All work at resolutions of order 0.1–1 nm and use literature cross sections (Category 4). We found no high-resolution cavity-enhanced broadband I₂ spectra (CRDS/CEAS line surveys, or comb-cavity spectra).

### 1.8 Iodine-filter absorption models

Iodine vapour filters (filtered Rayleigh scattering, Doppler global velocimetry, HSRL/wind lidar) need accurate Doppler-limited transmission near 532 nm (and 514.5 nm). Relevant work:
- **Forkey, Lempert & Miles, Appl. Opt. 36, 6729 (1997)** (`forkey1997`). A line-by-line absorption model compared with measured cell profiles over the Nd:YAG-SHG tuning range.
- **Earlier and applied filter work** (narrative only): Chan et al., Meas. Sci. Technol. 6, 784 (1995); Shibata et al., Jpn. J. Appl. Phys. 48, 032401 (2009).

These are few in number but provide absolute absorption profiles at known temperature and pressure. That makes them useful end-to-end checks of pyodine intensities and line shapes at 532 nm. Metrology cell-purity studies (Hrabina et al.) are in §5.

### 1.9 What the forward-citation searches added

- **From Salami & Ross (94 citing works):** Nölleke 2018; the Marcassa group papers (2022, 2023); the dual-comb papers (Ideguchi 2012; Eber 2025, 2026); IBBCEAS (Johansson 2013); iodine-cell purity studies (Hrabina 2014–2017); and a paper on frequency-scale correction of visible FTS (Serdyukov, J. Appl. Spectrosc. 83, 307 (2016), doi:10.1007/s10812-016-0287-0; content not accessed). The rest are non-I₂ spectroscopy that used the file for calibration.
- **From Knöckel 2004 (87):** Reiners 2024; Schmidt 2025 (ESPRESSO); Wang 2019; Nölleke 2018; the isotopologue work of Salumbides 2006/2008; the Novosibirsk emission-line frequencies (Category 2); the 129 I LIF-detection papers (Category 3); and many sub-Doppler absolute-frequency papers (Stage 2b).
- **From Gerstenkorn & Luc 1979 and 1985 (686 combined):** after 2003, mostly sub-Doppler metrology plus the Marcassa papers, Perdelwitz 2018, Williamson 2013, and St Petersburg/Edinburgh work on weakly bound and ion-pair states (Category 2 context).

No new FTS absorption atlas of I₂ after 2005 was found.

---

## 2. Fluorescence, term-value, and other-state data

This category covers data that fix the energy levels of the X and B states outside the range probed by room-temperature absorption. It also covers the other electronic states that perturb or predissociate B. For a potential-based global fit, these data do three things:

1. They extend the X potential above v″ = 17. Knöckel et al. (2004) fitted only v″ ≤ 17 and v′ ≤ 43; their data field is shown in their Fig. 2.
2. They pin the long-range parts and asymptotes of both potentials.
3. They supply the physics behind line widths and anomalous structure near the B dissociation limit.

All statements below were checked against a Crossref, HAL, OSTI or Semantic Scholar record, or against the article PDF. Where the only source of a number is a publisher abstract seen through a search-engine snippet, the catalog `notes` field says so.

### 2.1 X-state term values from laser-induced-fluorescence Fourier-transform spectroscopy (Lyon/Orsay, 1980–1986)

Nearly all precise knowledge of the X state above v″ ≈ 20 comes from one programme: the Lyon (Bacis, Churassy) and Laboratoire Aimé Cotton (Vergès) LIF-FTS work in the 1980s.

- **Bacis et al. 1980 (J. Chem. Phys. 72, 34).** B→X fluorescence excited at 514.5 and 501.7 nm and recorded at 1 mK precision for v″ = 10–100. The paper also reports eQq″ from linewidths, quasibound X levels above the rotationless limit, and perturbations of v″ ≥ 92 by two long-range states.
- **Martin, Bacis, Churassy & Vergès 1986 (J. Mol. Spectrosc. 116, 71).** The definitive X-state analysis. Rotational constants are given for every v″ ≤ 107, and v″ = 108–113 are extrapolated. Knöckel 2004 describes it as covering 7400–18350 cm⁻¹ with an uncertainty of "several 10⁻⁷". Knöckel cites the paper but did not fit these data.
- **Cerny, Bacis & Vergès 1986 (J. Mol. Spectrosc. 116, 458).** The same analysis for ¹²⁷I¹²⁹I and ¹²⁹I₂. It is the only extensive X-state term-value set for the minor isotopologues.
- **Bacis, Cerny & Martin 1986 (J. Mol. Spectrosc. 118, 434).** Long-range tests for all three isotopologues. Knöckel 2004 (Table 4) takes the X-state C6, C8 and C10 used in its potential extension from this paper.
- **Martin et al. 1983 (J. Chem. Phys. 79, 3725) and Churassy et al. 1981 (J. Chem. Phys. 75, 4863).** Near-IR LIF-FTS of the gerade states that converge on ²P₃/₂ + ²P₃/₂ (X 0g⁺, a′ 0g⁺, a 1g). Martin et al. obtain D₀(X) = 12 440.083 ± 0.145 cm⁻¹ and describe the X–a′–a interactions near R ≈ 5 Å that perturb the last X levels.

**Caveat on J coverage.** A fluorescence progression from one laser-excited level (v′, J′) samples only J″ = J′ ± 1 at each v″. The J coverage at high v″ is therefore a sparse set fixed by the available laser coincidences, mainly Ar⁺ and Kr⁺ lines. It is not a uniform grid. We could not determine from the available abstracts whether the full line lists were printed or deposited. None is known to exist in machine-readable form.

### 2.2 Sub-MHz anchors for high v″: three-level "emission" spectroscopy (Novosibirsk, 2008–2019)

The forward-citation searches found no new broadband X-state term-value data after 2003. The searches covered the 130 works citing Martin 1986 and the 87 citing Knöckel 2004 (OpenAlex).

The only modern high-v″ data come from the Institute of Laser Physics, Novosibirsk, using three-level laser spectroscopy. A 532 nm Nd:YAG pump and a near-IR diode-laser probe are each locked to hyperfine components that share an upper level, and a femtosecond comb measures the probe frequency.

| Paper | Range | Transitions | Items | Rel. uncertainty | X level reached |
|---|---|---|---|---|---|
| Matyugin et al. 2012, Quantum Electron. 42, 250 | 982–985 nm | R56, P58 (32–48); P85, R86 (all 15 comps.), R87, R88 (33–48) | 20 HFS comps. | 7×10⁻¹⁰ | v″ = 48 |
| Nesterenko et al. 2019, Quantum Electron. 49, 633 | 1053–1068 nm | six lines, bands 32–53 and 32–54 | 18 HFS comps. | 8×10⁻¹⁰ | v″ = 53, 54 |
| Matyugin et al. 2008, Quantum Electron. 38, 755 | not retrieved | HFS of emission lines (method paper) | – | – | – |

The upper levels (v′ = 32, 33) are shared with the well-characterized 532 nm lines, so each emission frequency gives an X term value at v″ = 48, 53 or 54 to about 0.2 MHz. That is roughly 10³ times better than LIF-FTS. These are the most valuable X-state points found in this category, few as they are.

A related source is the Hannover cw I₂ Raman laser (Klug et al. 2000). Pumping at 532 nm, it lases on more than 50 lines between 544 and 1330 nm for all three isotopologues. We did not verify whether absolute line frequencies were reported.

### 2.3 Older X-state sources (historical and cross-check only)

- **UV resonance series (Verma 1960).** D₀ = 12 452.5 ± 1.5 cm⁻¹. LeRoy (1970) reassigned the series, which shows that the near-dissociation assignments are unreliable. LeRoy's own value is D₀ = 12 440.9 ± 1.1 cm⁻¹.
- **Resonance Raman (Holzer, Murphy & Bernstein 1970; Kiefer & Bernstein 1972).** ωe and ωexe from overtone progressions with Q, S and hot-band structure. Much less precise than LIF-FTS.
- **cw optically pumped I₂ laser (Koffend, Bacis & Field 1979).** Lasing from B v′ = 43 to X v″ = 42, 64 and 76–83, per the 1980 follow-up. It is a demonstration of access to high v″ rather than a term-value source.
- **Representations to dissociation.** Ashmore & Tellinghuisen 1986 applied mixed polynomial/NDE fits to I₂(X). They are useful as benchmarks for extrapolating the X potential.

### 2.4 B state near dissociation and long range

Knöckel's fit stops at v′ = 43. Above that, the level structure depends on a patchwork of sources:

| Source | B levels | Nature | Key numbers |
|---|---|---|---|
| Luc (FTS constants, via Tromp 1983) | v′ = 1–62 | Doppler-limited FTS (Category 1) | – |
| Barrow & Yee 1973 | v′ ≈ 63–77 | band rotational analysis | long-range B analysis |
| Danyluk & King 1977 | v′ = 77–82 | two-photon sequential absorption | Bᵥ for v′ > 78 shown erroneous by Tromp 1983 |
| Steinfeld et al. 1969 | up to ~24 cm⁻¹ below limit | absorption | 93 ± 2 bound levels estimated |
| Tromp, Le Roy, Gerstenkorn & Luc 1983 | analysis to v′ ≈ 82 | NDT re-analysis | limit 20 043.16(2) cm⁻¹, v_D = 87.32(4), **C5 = 2.88(3)×10⁵ cm⁻¹ Å⁵**, D₀(X) = 12 440.18(2) cm⁻¹ |
| Gerstenkorn, Luc & Amiot 1985 | IPA potential, v′ = 0–80 | long-range multipole fit | De(B) = 4381.2492(10) cm⁻¹, **C5 = 3.161(33)×10⁵**, C6 = 1.506(43)×10⁶, C8 = 0.248(14)×10⁸, C10 ≤ 0.042(3)×10¹⁰ |
| Gerstenkorn, Luc & Le Roy 1991 | – | constants for ¹²⁷I¹²⁹I and ¹²⁹I₂; asymptotes | Knöckel uses De(X) = 12 547.340 and De(B) = 20 150.317 cm⁻¹ from this reference |
| Pique et al. 1984 (PRL), 1986 (J. Phys.) | near the ²P₃/₂+²P₁/₂ limit; v′ = 76–78 strongly perturbed | sub-Doppler hyperfine (>10 000 lines) | B–a 1g hyperfine coupling breaks u–g symmetry |

Knöckel 2004 fixes the B-state Cn at the Gerstenkorn–Luc–Amiot values. The C5 of Tromp et al. is about 9 % smaller. This tension should be carried into our fit, either as a free parameter or as a prior informed by theoretical C5.

Near the limit, B levels are mixed with a 1g through the hyperfine interaction. This produces extra lines and "superhyperfine" structure that a single-channel B potential cannot reproduce. Modeling B–X absorption below about 505 nm therefore needs a coupled-channel or effective-Hamiltonian treatment, or else the lines have to be excluded.

### 2.5 Other electronic states that matter for B–X modeling

- **Predissociating states.** Natural linewidths come from gyroscopic, hyperfine and magnetic predissociation of B by the states that cross it.
  - Vigué, Broyer & Lehmann (1981; papers I–III, open access on HAL) give the complete experimental and theoretical treatment, with v′ and J dependence of the predissociation and radiative rates. Broyer, Lehmann & Vigué (1975) give B-state gJ factors and lifetimes. Knöckel 2004 compared its nonadiabatic BOC function α(R) with these gJ values but did not fit them.
  - Tellinghuisen (1985) derived potentials for the 1Πu (B″ 1u) and a 1g states, including their repulsive branches where they cross B at R < 3.3 Å.
  - Inard et al. (1999) gave the first rotational analysis of B″ 1u through E→B″ FTS fluorescence: Te = 12 344.804(40) cm⁻¹, De = 202.53(4) cm⁻¹, with long-range C5, C6 and C8. These numbers come from the publisher abstract and should be checked.
  - Lifetime and linewidth data proper are in §4.
- **A 3Π1u.** Appadoo et al. (1996) analysed 9552 FTS absorption lines in 79 A–X bands (v′ = 0–35, v″ = 3–17). A–X lines are interleaved with B–X hot bands in the near IR, so they must be modeled or masked when fitting the 11 000–14 000 cm⁻¹ atlas region. They also give independent X v″ ≤ 17 combination differences. Tellinghuisen (2003) re-represented these data.
- **A′ 2u.** Koffend et al. 1982 (collisionally induced double resonance: De(A′) = 2505.7 ± 1.9 cm⁻¹, T₀ = 9988.7 ± 1.6 cm⁻¹) and Cerny et al. 1997 (D′→A′ LIF-FTS). These are relevant only to energetics and collisional context.
- **Ion-pair and weakly bound (bb) valence states.** The Edinburgh (Ridley, Lawley, Donovan 2007) and St Petersburg (Baturo, Lukashov, Poretsky, Pravilov 2013–2016) groups characterized weakly bound states correlating with I*+I*. They lie far above the B–X region and do not supply X- or B-state term values, so they are not catalogued.

### 2.6 What Knöckel 2004 did and did not use (from the article PDF)

- **Fitted data.** About 1500 rovibronic frequencies with v″ ≤ 17 and v′ ≤ 43. These were the Hannover NIR measurements (778–815 nm), sub-Doppler grids, Kato-atlas-derived frequencies (526–667 nm), and recalibrated Gerstenkorn–Luc atlas lines. A few ¹²⁹I₂ and ¹²⁷I¹²⁹I lines were added for Born–Oppenheimer corrections.
- **Fixed long-range parameters.** X: C6, C8 and C10 from Bacis et al. 1986. B: C5–C10 from Gerstenkorn–Luc–Amiot 1985. Asymptotes De(X) and De(B) from Gerstenkorn–Luc–Le Roy 1991.
- **Cited but not fitted.** The Martin et al. 1986 X-state fluorescence data, and the Broyer et al. 1975 g-factors (used for comparison).

### 2.7 Gap notes (Category 2)

1. **X state, v″ = 18–47.** The only data are 1980s LIF-FTS (about 0.001–0.005 cm⁻¹, or 30–150 MHz), plus a few sub-Doppler NIR hot-band measurements catalogued elsewhere. There are no sub-MHz anchors apart from v″ = 48, 53 and 54 (Novosibirsk). This range governs hot-band absorption in heated cells and all near-IR lines below about 12 000 cm⁻¹.
2. **X state, v″ = 55–113 (to dissociation).** Only the Lyon LIF-FTS (v″ ≤ 107) and long-range analyses exist. The highest levels are perturbed by a′ 0g⁺ and a 1g, and quasibound levels are broadened. No data were found after 1986.
3. **J coverage at high v″.** LIF gives sparse J″ sets tied to laser coincidences. We could not quantify the coverage without the full papers. It is probably the largest structural weakness when extrapolating centrifugal distortion at high v″ and high J″.
4. **B state, v′ = 44–62.** The data are the Doppler-limited atlas (Category 1) plus sub-Doppler lines near 500–515 nm. None of it was in Knöckel's fit.
5. **B state, v′ ≥ 63.** The only sources are 1970s photographic and two-photon work (Barrow & Yee; Danyluk & King, partly erroneous) and the Lyon hyperfine studies. v′ = 76–78 are strongly perturbed by a 1g. The long-range C5 is inconsistent between Tromp 1983 (2.88×10⁵) and Gerstenkorn 1985 (3.161×10⁵).
6. **Isotopologues.** The X-state term values for ¹²⁷I¹²⁹I and ¹²⁹I₂ rest on a single 1986 LIF-FTS study (plus the long-range test).
7. **Machine readability.** Not one Category-2 dataset was found in electronic form. Every one will need digitization from printed tables or a request to the authors. Several Journal de Physique papers are open access on HAL.
8. **Opportunity.** A modern comb-referenced or high-resolution FTS LIF campaign, pumped at several well-known sub-Doppler lines (such as 532, 515 and 633 nm), would put the whole X ladder on a common absolute scale. The Novosibirsk three-level technique shows that such anchors can reach below 1 MHz.

---

## 3. Isotopologues: ¹²⁹I₂ and ¹²⁷I¹²⁹I

**Summary.** The data for the heavier isotopologues are sparse and fall into three groups:

1. Sub-Doppler hyperfine spectra from 1970s–80s laser metrology at the He-Ne lines (633, 612, 640 nm).
2. One laser-induced-fluorescence Fourier-transform (LIF-FTS) study of the ground state (Lyon/Orsay, 1986).
3. One modern systematic study from Amsterdam, Zürich and Hannover (Salumbides et al. 2006, 2008).

No broadband FTS absorption atlas of ¹²⁹I₂ or ¹²⁷I¹²⁹I was found. Searches covered Crossref, OpenAlex, Semantic Scholar, web search, and forward citations of Salumbides 2006/2008, Cerny 1986, Gläser 1981 and Wu 1985.

### What Knöckel et al. (2004) used

Section 2 of Knöckel et al. (2004) lists exactly one ¹²⁷I¹²⁹I line and six ¹²⁹I₂ lines, all taken from the CIPM compilation (Quinn 2003):

| Isotopologue | Lines | Region |
|---|---|---|
| ¹²⁷I¹²⁹I | P(33) 6-3 | 633 nm |
| ¹²⁹I₂ | P(33) 6-3, P(54) 8-4, R(69) 8-4, P(69) 12-6 | 633 nm |
| ¹²⁹I₂ | P(110) 10-2, R(113) 14-4 | 612 nm |

The authors call this basis "still very poor" for determining the adiabatic Born–Oppenheimer correction (BOC). The BIPM mise en pratique (MeP) 633 nm tables list **R(60) 8-4** rather than R(69) 8-4, so the Knöckel list probably has a typo. This is unresolved and noted in the catalog.

The MeP 2003 tables themselves hold 180 tabulated component intervals for these seven lines (entry `quinn2003`), with 1σ uncertainties of 0.03–5 MHz. They are anchored to ¹²⁷I₂ by beat measurements, for example f(a28, ¹²⁹I₂ P(54) 8-4) − f(a16, ¹²⁷I₂ R(127) 11-5) = −42.99(4) MHz.

### Metrology-era hyperfine data (1971–1985)

**633 nm**
- Knox & Pao (1971): 38 ¹²⁹I₂ saturation features.
- Tesic & Pao (1975): assignments for ¹²⁹I₂ and ¹²⁷I¹²⁹I.
- Schweitzer et al. (1973, NBS): ¹²⁹I₂-stabilized lasers, with Kr-referenced wavelengths and pressure broadening.
- Magyar & Brown (1980): ¹²⁹I₂, ¹²⁹I¹²⁷I and ¹²⁷I₂.

**612 nm**
- Ciddor & Brown (1980).
- Dschao, Gläser & Helmcke (1980): 25 components of ¹²⁹I₂ P(110) 10-2.
- Gläser, Dschao & Foth (1981): enriched ¹²⁹I₂, hyperfine structure plus fluorescence analysis.
- Wu, Gaida & Bialas (1985): one ¹²⁷I¹²⁹I line, including its absorption coefficient.

**640 nm**
- Gläser (1985): identification, with calculated positions for all three isotopologues.
- Khodakovskiy et al. (2010): FM spectroscopy of a ¹²⁷I¹²⁹I line, with ¹²⁹I₂ P(90) 5-3 also in the scan.

Most of these were compiled into the BIPM tables. Their value now is mainly as cross-checks of an isotopologue hyperfine model.

### Ground-state term values

Cerny, Bacis & Vergès (1986) is the only extensive X-state dataset found for ¹²⁷I¹²⁹I and ¹²⁹I₂ (LIF-FTS). Bacis, Cerny & Martin (1986) tested X-state long-range behaviour across all three isotopologues.

Abstracts were not available through any API, so v″/J coverage and uncertainties still have to be read from the papers. They are rated **high priority**: Salumbides et al. (2008) could not separate X- and B-state BOC contributions, and independent X-state isotopologue data are exactly what is needed to do so.

### Salumbides et al. (2006, 2008): the key post-2003 data

The method was double saturation spectroscopy, with a ¹²⁹I₂/¹²⁷I¹²⁹I cell and a ¹²⁷I₂ reference recorded simultaneously.

- **Coverage:** 573–583 nm, 610–615 nm and near 633 nm; v″ = 0–5, v′ = 9–20, J′ = 23–128.
- **Molecular Physics (2006):** more than 200 eqQ(B) and C(B) values for ¹²⁹I₂ and more than 170 for ¹²⁷I¹²⁹I (typically 2 MHz and 2 kHz). Also gives isotopically scaled hyperfine formulae valid for v″ ≤ 17 and v′ ≤ 53. The tables are in an electronic supplement.
- **EPJ D (2008):** more than 290 (¹²⁹I₂) plus more than 90 (¹²⁷I¹²⁹I) difference frequencies to ¹²⁷I₂ lines. Uncertainty is 1 MHz (1σ) for the absolute calibration and ≤ 1.5 MHz for the differences.
  - A joint direct-potential fit gives effective adiabatic and non-adiabatic BOC functions for all three isotopologues, and Tₑ(¹²⁷I₂) − Tₑ(¹²⁹I₂) = 94(11) MHz.
  - The data set and potentials are in electronic supplementary material at epj.org (not downloaded, so whether it is machine-readable is unconfirmed).
  - Stated prediction accuracy is about 1.5 MHz within the observed range.

These two papers are the obvious isotopologue backbone for pyodine's global fit. Forward citations (24 and 31) turned up no later dedicated ¹²⁹I₂/¹²⁷I¹²⁹I spectroscopy apart from the Kireev LIF-detection papers.

### Other isotopologue items (low priority)

- **Kireev et al. (2012–2016):** isotope-selective LIF detection of ¹²⁹I with Cu-vapor (510.6 and 578.2 nm), Nd:YAG-SHG and diode lasers. These are useful only as consistency checks of predicted absorption at fixed laser lines.
- **Wieland et al. (1972):** E→B and F→X emission of ¹²⁷I₂ and ¹²⁹I₂.
- **King & McLean (1989):** three-photon absorption of ¹²⁷I₂ and ¹²⁹I₂.
- **Chartier et al. (1992):** found ¹²⁹I₂ contamination in some commercial ¹²⁷I₂ cells. This matters for modeling real cells, since weak ¹²⁹I₂ and ¹²⁷I¹²⁹I lines can appear in "pure" ¹²⁷I₂ spectra.

### Gap notes (isotopologues)

- **Spectral and vibrational coverage.** Precision data for ¹²⁹I₂ and ¹²⁷I¹²⁹I exist only at 573–640 nm, with v′ ≤ 20 and v″ ≤ 5, apart from the single ¹²⁹I₂ P(90) 5-3 feature at 640 nm.
  - Nothing sub-Doppler exists in the green (500–560 nm), near dissociation (v′ > 40), in the red and near-IR hot bands (> 640 nm), or for J″ < 23.
  - There is no broadband FTS absorption or cross-section measurement of either isotopologue, and no B-state lifetime or predissociation data.
- **Hyperfine parameters.** Nuclear spin–spin parameters for ¹²⁹I₂ come from one line, P(69) 12-6.
- **Consequence for the fit.** Predictions for ¹²⁹I₂ outside 573–640 nm rest entirely on mass-scaling of ¹²⁷I₂ potentials plus effective BOC functions constrained over a narrow range. A broadband FTS spectrum of an enriched ¹²⁹I cell, or comb-referenced sub-Doppler lines near 515–532 nm, would be the highest-value new measurements.
- **Commercial ¹²⁹I₂ cells.** Cells are made for 633 nm metrology, so a lab measurement is feasible.

---

## 4. Intensity data: cross sections, transition moment, lifetimes, vapor pressure

This category covers everything that sets the **strength and width** of B–X absorption, as opposed to line positions: absolute cross sections, the B–X electronic transition dipole moment (TDM) function and absolute line strengths, B-state radiative lifetimes and predissociation rates (these set natural linewidths), and I₂ vapor pressure (sets the number density in a cell). Catalog ids are given in brackets. "Inferred" marks an interpretation made here, as opposed to what the source states.

### 4a. Absolute absorption cross sections

**Modern reference data (post-2000).** Three datasets matter most:

- **Saiz-Lopez et al. 2004** [saizlopez2004] used FTS at 4 cm⁻¹ resolution (0.1 nm at 500 nm) to measure the room-temperature cross section from 182 to 750 nm (295 K, 760 Torr). The visible maximum is σ(533.0 nm) = (4.24 ± 0.50) × 10⁻¹⁸ cm². The spectrum is an ACP supplement (zip), and it underlies both the JPL (2006/2010) and IUPAC (2007) recommendations, which are 5-nm averages [sander2011, atkinson2007].
- **Spietz, Gómez Martín & Burrows 2006** [spietz2006a, spietz2006b] measured σ(500 nm) = (2.186 ± 0.021) × 10⁻¹⁸ cm² by a method independent of the iodine vapor pressure. That makes it the most precise absolute anchor found (~1 %). Combined with earlier values, their weighted mean is (2.191 ± 0.02) × 10⁻¹⁸ cm². The same paper gives DOAS reference spectra at 0.25 nm and 0.59 nm FWHM (±4 % and ±3 %). It also shows that retrieved column densities depend strongly on the column density of the reference spectrum at ≥1 nm resolution. This is a useful test of whether a line-by-line model handles optically thick, instrument-convolved spectra correctly.
- **Tellinghuisen 2011** [tellinghuisen2011a, tellinghuisen2011b] gives quantitative absorption at 0.1 nm resolution over 520–640 nm, analyzed by least-squares spectral simulation into discrete B–X and continuum components. It yields the B–X |μₑ|² with **<2 % relative standard error** over most of that region and lowers the C(¹Πᵤ)←X strength by 25 %. The MPI-Mainz `Tellinghuisen(2011)_294-300K_520-635nm` file is **not** a B–X spectrum. It is the paper's Table II continuum absorptivity ε_c, converted to σ, at 18 wavelengths (38 rows, with standard deviations). This correction comes from the stage-3 cross-section data review (2026-09-15); see `cross-section-data.md`. The 0.1 nm spectra themselves were not released. The companion C-state paper [tellinghuisen2011c] provides 2 nm resolution spectra at 273, 308 and 337 K (390–900 nm), the only modern temperature-dependent set.

**Historical data** [tellinghuisen1973a, tellinghuisen1973b, tellinghuisen1982, rabinowitch1936, vogt1923, kortum1947, mathieson1956, mcmillan1966, bauer1998] are low resolution. Their main value is coverage. Tellinghuisen 1973 split the 420–800 nm absorption into ¹Πᵤ←X, B←X and A←X contributions, and his 1982 paper added "between-the-lines" measurements corrected for residual B–X absorption. Mathieson & Rees 1956 is one of the few absolute data sets reaching 848 nm.

**High temperature.** Sulzer & Wieland 1952 [sulzer1952] measured extinction coefficients over 13,500–28,000 cm⁻¹ at seven temperatures from 423 to 1323 K. This is the only data set found above ~613 K; Kortüm & Friedheim 1947 [kortum1947] went to 613 K. Nothing modern exists between 340 and 1300 K.

**Machine-readable access.** The MPI-Mainz UV/VIS Spectral Atlas [kellerrudek2013] lists 21 I₂ datasets (1923–2011; 273–1323 K; 170–900 nm) as ASCII files. Digitized historical curves were converted with the factor 3.8235 × 10⁻²¹. This is the practical ingestion route for all cross-section data. The JPL 19-5 evaluation (Burkholder et al. 2019) could not be reached: the site was unreachable from this environment, so whether its I₂ entry changed is unchecked.

| id | T (K) | range (nm) | resolution | stated uncertainty | machine-readable |
|---|---|---|---|---|---|
| saizlopez2004 | 295 | 182–750 | 4 cm⁻¹ FTS | ±0.50e-18 at 533 nm (~12 %) | yes (ACP supplement, MPI) |
| spietz2006a | 298 | 500 | — | ±1 % | yes (MPI) |
| spietz2006b | 298 | 428–588 | 0.25 / 0.59 nm | ±4 % / ±3 % | yes (MPI) |
| tellinghuisen2011a | 294–300 | 520–635 | Table II only: continuum ε_c at 18 wavelengths (spectra not released) | per-point s.d. | yes (MPI) |
| tellinghuisen2011c | 273–337 | 390–900 | 2 nm | not extracted | yes (MPI) |
| tellinghuisen1973a | 295–300 | 420–800 | low | not extracted | yes (MPI) |
| sulzer1952 | 423–1323 | 357–741 | low | not extracted | yes (MPI, digitized) |
| mathieson1956 | 393 | 635–849 | low | not extracted | yes (MPI, digitized) |

### 4b. Transition dipole moment, Franck–Condon factors, absolute line strengths

- **Experimental TDM functions:**
  - Koffend, Bacis & Field 1979 [koffend1979a] used gain in an optically pumped I₂ laser.
  - Bhale et al. 1985 [bhale1985] used laser-excited fluorescence.
  - Lamrini et al. 1994 [lamrini1994] cover R-centroids **2.633–6.035 Å**, over which |μₑ| falls by about four orders of magnitude. This is the widest experimental R range and is needed for transitions to high v″, i.e. hot bands and the red/near-IR.
  - Tellinghuisen 1997 [tellinghuisen1997] inverted B-state radiative decay rates versus v′ with a sum rule. The strength peaks at 1.85 ± 0.1 D² near R = 3.3 Å, but he rated the quantitative consistency of earlier TDM studies as poor.
  - Tellinghuisen 2011 [tellinghuisen2011b] is now the precision benchmark within its R range.
- **Ab initio:** Zaitsevskii et al. 2000 [zaitsevskii2000] computed B–X, A–X and 1ᵤ–X transition intensities. With the relativistic potential calculations of Teichteil & Pelissier 1994 [teichteil1994] and de Jong et al. 1997 [dejong1997], these are the only guides for extrapolating the TDM beyond the measured R range and for the A–X and 1ᵤ–X continua.
- **Absolute single-line strengths:** Dubé & Trinczek 2004 [dube2004] measured integrated absorption of P(78)1–9, R(86)1–9 and R(113)3–10 near 718 nm. With literature FC factors this gives μₑ = 1.10(3) D at R-centroid 0.293 nm. It is the only modern absolute line-strength measurement on individual B–X lines found in this search. Suwaiyan et al. 1992 [suwaiyan1992] covers a 2 cm⁻¹ window near 588 nm.
- **FC / intensity factors:** Zare 1964 [zare1964] and Tellinghuisen 1978 [tellinghuisen1978] are historical. pyodine will recompute these from its own potentials.

### 4c. B-state lifetimes, predissociation and natural linewidths

The natural width of a B–X hyperfine component is the radiative rate plus predissociation. Predissociation comes from coupling to repulsive states, mainly B″(1ᵤ). The couplings are gyroscopic (∝ J′(J′+1)), hyperfine, and magnetic in an external field. The key sources:

- **Radiative / total lifetimes vs v′:**
  - Capelle & Broida 1973 [capelle1973] is the broadest survey: v′ ≈ 5–70 excited at 640–499.5 nm, lifetimes from <0.4 µs to >7 µs (inferred natural widths ≈ 20–400 kHz), and self-quenching cross sections of 47–90 × 10⁻¹⁶ cm².
  - Rotationally selected lifetimes and self-quenching: Sakurai et al. 1971 [sakurai1971], Paisner & Wallenstein 1974 [paisner1974], and Klein et al. 1994 [klein1994].
- **Predissociation (ENS Paris series):**
  - Broyer, Vigué & Lehmann 1975 [broyer1975b] demonstrated the J′(J′+1) dependence of 1/τ (gyroscopic predissociation).
  - Related ENS papers: [broyer1975a, vigue1975, broyer1976].
  - Vigué, Broyer & Lehmann 1981, parts I–III [vigue1981a–c], give the full theory plus experiments on natural, hyperfine and magnetic predissociation.
  - Hyperfine-sublevel-resolved lifetimes: Pique et al. 1983 [pique1983] and Tench & Ezekiel 1983 [tench1983].
  - Later lifetime analyses: Martínez, Castaño et al. [martinez1988].
  - Tellinghuisen 1972 and 1985 [tellinghuisen1972, tellinghuisen1985] derived potentials for the weakly bound and repulsive states from diffuse spectra and predissociation data.
- **Sub-Doppler linewidths near dissociation:** Cheng et al. 2002 [cheng2002] mapped hyperfine linewidth versus wavelength over 523–498 nm, with the narrowest ≈ 4 kHz near 508 nm. This overlaps the sub-Doppler catalog. The measured widths include transit-time and pressure contributions.
- **Model implication (inferred):** a natural-width model Γ(v′, J′, F) = A_rad(v′) + k_gyro(v′) J′(J′+1) + hyperfine term can be parameterized from these data. The data are almost all from 1971–1994, mostly without machine-readable tables.

### 4d. I₂ vapor pressure and cell modeling

Converting a cold-finger temperature (typically −15 to +40 °C in metrology cells) or a cell temperature into I₂ number density requires the vapor pressure of solid iodine. Candidate sources:

- Gillespie & Fraser 1936 [gillespie1936]: the standard crystalline-iodine reference (gas-current method).
- Baxter et al. 1907 and Baxter & Grose 1915 (50–95 °C) [baxter1907, baxter1915].
- Berkenblit & Reisman 1966 (43–80 °C) [berkenblit1966].
- Stull 1947, via the NIST WebBook Antoine fit: A = 3.36429, B = 1039.159, C = −146.589, valid only 311.9–456 K [stull1947].
- Shirley & Giauque 1959 [shirley1959]: heat capacity (13–327 K) and heat of sublimation, which allow thermodynamically consistent extrapolation below 0 °C.
- Honig & Hook 1960 [honig1960]: existence confirmed, iodine content not checked.

Hrabina et al. 2014 [hrabina2014] studied cells filled to a specified saturation pressure, with no cold finger. Fredin-Picard 1989 [fredinpicard1989] addresses cell impurity effects.

### Gaps (Category 4)

1. **No modern absolute, rotationally resolved B–X intensities.** Absolute data are at ≥0.1 nm resolution (Saiz-Lopez, Tellinghuisen) or are single lines (Dubé & Trinczek, near 718 nm). There is no absolute line-strength dataset at Doppler-limited resolution anywhere in 500–680 nm. FTS scans of iodine cells from the Category 1 astronomical-cell data could fill this if the I₂ column (cell temperature and fill) were known.
2. **Temperature dependence** is covered only at 273–337 K (Tellinghuisen 2011, 2 nm resolution) and 423–1323 K (Sulzer & Wieland 1952, low resolution, digitized). There is nothing modern at the 320–600 K temperatures of heated iodine filters.
3. **Red/near-IR (>650 nm to ~900 nm)** absolute cross sections rest on Mathieson & Rees 1956 and Tellinghuisen's 2 nm spectra. This is the region of the 11,000–14,000 cm⁻¹ atlas and of A–X overlap.
4. **TDM** is measured only for R-centroid 2.63–6.04 Å (Lamrini 1994); the precise 2011 determination covers a narrower range. Literature consistency is poor (Tellinghuisen 1997), and the only ab initio function is from 2000.
5. **Lifetimes and predissociation** are mostly from 1971–1994 and printed only. Coverage in J′ is sparse, hyperfine-dependent rates exist for only a few levels, and there is no modern systematic re-measurement. The v′ > 70 region is covered only by sub-Doppler linewidths (Cheng 2002).
6. **Vapor pressure** below ~0 °C (cold-finger range) is not directly measured in the sources checked; the NIST Antoine fit starts at 311.9 K. A thermodynamic fit combining Gillespie & Fraser with Shirley & Giauque is needed. Temperature ranges and uncertainties of the older papers were not extracted (full texts paywalled).

**Handled elsewhere in this document:** Williamson et al., J. Chem. Educ. 90, 383 (2013), doi 10.1021/ed300455n, is reported to be a 7,381-line B–X data set with no intensity content; it is catalogued in §1.4 (`williamson2013`). Johansson et al. 2013 (IBBCEAS) is in §1.7 (`johansson2013`). The J. Phys. B 47, 055101 (2014) bound–free cross sections are for ion-pair-related excitation, not B–X absorption, so they are not catalogued. Myer & Samson 1970 and Roxlo & Mandl 1980 (VUV/UV) are available through MPI-Mainz but not cataloged.

---

## 5. Line-shape data: pressure broadening and shift

**Summary.** No Doppler-limited *self*-broadening measurement of I₂ B–X lines was found. Foreign-gas broadening and shift at Doppler-limited resolution exist for only two spectral spots:

- 543 nm in air (Fletcher & McDaniel 1995, with temperature scaling).
- 675 nm with 12 buffer gases (Wolf 2009 thesis, with limited temperature dependence).

The sub-Doppler metrology literature adds many *self*-pressure-shift coefficients (−1 to −8 kHz/Pa is typical) and a few self-broadening numbers. These belong to specific lock schemes, and cell impurities often dominate them.

No speed-dependent, Dicke-narrowing or line-mixing study of I₂ B–X was found. Wolf (2009) discusses asymmetric line shapes but fits Voigt profiles.

### Doppler-limited foreign-gas broadening and shift

The table uses Wolf's six-line averages at 292 K (FWHM) and Fletcher's air values converted to 293 K.

| Perturber | Broadening (MHz/Torr) | Shift (MHz/Torr) | Line(s) | Source |
|---|---|---|---|---|
| He | 8.34(31) | −0.201(50) | 675 nm, J″ = 69–133 | Wolf 2009 |
| Ne | 6.17(39) | −0.71(3) | 675 nm | Wolf 2009 |
| Ar | 7.70(58) | −1.32(11) | 675 nm | Wolf 2009 |
| Kr | 6.86(33) | −1.40(6) | 675 nm | Wolf 2009 |
| Xe | 6.47(29) | −1.71(10) | 675 nm | Wolf 2009 |
| H₂, D₂, N₂, CO₂, N₂O, air, H₂O | see thesis Tables 6.6–6.7 | see thesis | 675 nm | Wolf 2009 |
| air | ≈ 10.0 (0.531 ± 0.009 GHz K^0.7/Torr) | ≈ −1.09 (−0.058 ± 0.004 GHz K^0.7/Torr) | R(97) 28-0 and P(60) 27-0, 543 nm | Fletcher & McDaniel 1995 |

**Temperature dependence**
- Fletcher & McDaniel found both shift and broadening scale as P/T^0.7, not the T^0.5 assumed in the Hiller & Hanson (1990) review. They also state that the review's shift value is wrong.
- Wolf measured Ar at 292, 348 and 388 K, and He and CO₂ at 292 and 388 K.

**Sub-Doppler foreign-gas broadening.** Phillips & Perram (2008) measured Ar broadening of hyperfine-resolved P(10) and P(70) of the (17,1) band. The coefficient values were not accessible to us.

### Self-broadening (sub-Doppler, iodine pressure)

| Source | Line | Value | Comment |
|---|---|---|---|
| Schweitzer et al. 1973 | ¹²⁹I₂ at 633 nm | ≈ 13 MHz/Torr (98 kHz/Pa), 2.6 MHz intercept | Saturated-resonance width |
| Eickhoff & Hall 1995 | 532 nm | 148 kHz/Pa (≈ 19.7 MHz/Torr) | Quoted second-hand by Kobayashi et al. 2015 |
| Fang, Wang & Shy 2006 | R(56) 32-0 a10 | — | Dedicated pressure/power broadening study; values not extracted |
| Sakurai et al. 1979; Cérez et al. 1974 | 633 nm | — | Values not extracted |
| Weller et al. 2026 | R(41) 51-0 at 508 nm | — | Linewidth vs vapor pressure |
| Hrabina et al. 2014, 2017 | R(56) 32-0 a10 | Δν = 177 kHz + 242 kHz/Pa × K₀ | Linewidth vs Stern–Volmer contamination K₀ |

Converting sub-Doppler resonance widths to homogeneous impact widths requires care: roughly a factor of 2, plus power and transit broadening and modulation effects.

### Self-pressure-shift coefficients

| λ (nm) | Line / component | Shift (kHz/Pa) | Pressure range | Source |
|---|---|---|---|---|
| 501.7 | R(26) 62-0 | "less than −3.8" | op. 0.33 Pa | Goncharov et al. 2007 |
| 508 | R(41) 51-0 | (widths only) | — | Weller et al. 2026 |
| 520.2 | a1 of P(34), R(36), P(33), R(35) 39-0 | −5.00(16), −4.80(13), −5.03(47), −5.23(54) | 1.44–3.36 Pa | Nishiyama et al. 2024 |
| 531.5 | R(36) 32-0 a1 | −1.3 | ≈ 41 Pa | Kobayashi et al. 2015 |
| 532.2 | R(56) 32-0 a10 | −4.2 (MeP, from Nevsky et al. 2001); −3.4 ± 0.6 (Hrabina et al. 2017) | 0.8–4 Pa | BIPM MeP; Hrabina 2017 |
| 554 | P(69) 25-1 b18 | −78.2 ± 8.8 (outlier?) | 2.6–3.9 Pa | Chen et al. 2024 |
| 578 | R(37) 16-1 a1 | −7.7 | — | Kobayashi et al. 2016 |
| 633 | R(127) 11-5 | −15 kHz/°C of cold finger (≈ −9 kHz/Pa, estimate made here) | ≈ 17 Pa | BIPM MeP 633 |
| 633 (¹²⁹I₂) | several components | ≈ −7 (conversion made here of "2×10⁻⁹ ν/Torr") | — | Schweitzer et al. 1973 |
| 640 | P(10) 8-5 a9 | −7.8 | 16–23 Pa | BIPM MeP 640 |

**Values seen but not verified or not cataloged**
- Web-search summaries quoted coefficients for 548 nm (P(28) 24-0; Hsiao et al., JOSA B 30, 328, 2013), 534 nm (Cheng et al., JOSA B 36, 1816, 2019) and 532 nm (Hong et al., Opt. Commun. 235, 377, 2004; "−1.3 kHz/Pa over 0.4–4.0 Pa"). None could be checked against the papers.
- Zhang et al. (Appl. Opt. 48, 5629, 2009) and Gläser (IEEE TIM 36, 604, 1987) studied pressure broadening and pressure effects at 560 and 612 nm, but no values were extracted.

**Cell purity.** Gläser (1982) found non-linear low-pressure shifts from a condensable impurity. Chartier et al. (1992), Lazar et al. (2009) and Hrabina et al. (2014, 2017) all show that residual foreign gas in sealed cells often dominates measured shifts and widths at the 10 kHz level. This is directly relevant if pyodine is to predict real-cell line shapes.

### Gap notes (line shape)

1. **No self-broadening at Doppler-limited resolution.** There is no measurement of the I₂–I₂ collisional (Lorentzian) width of B–X lines in linear absorption. The metrology numbers (≈ 100–150 kHz/Pa for sub-Doppler resonances, ≈ 10–20 MHz/Torr) are the only proxies.
2. **Foreign-gas coverage is two points.** Only 543 nm (air) and 675 nm (v′ = 4–6, J″ = 69–133, 12 gases) are covered. There is nothing for 500–530 nm, where most astronomical-cell and filter work happens, nothing near dissociation, and nothing for N₂/air at 532 nm apart from what may sit inside iodine-filter models (Category 1).
3. **J and v dependence.** This is essentially unmeasured: Wolf's low-J Ar data are inconclusive, and Phillips & Perram compare only J″ = 10 and 70.
4. **Temperature exponents.** These exist only for air (T^−0.7) and for Ar, He and CO₂ up to 388 K. Heated cells (hot-band atlases, iodine filters at 330–350 K, astronomical cells at 50–70 °C) need exponents at higher temperatures.
5. **Self-shift coefficients.** These have scattered values that depend on the lock method, and there is no systematic study versus v′ and J. The 554 nm outlier (−78 kHz/Pa) needs independent confirmation.
6. **Line-shape physics beyond Voigt.** No speed-dependence, Dicke-narrowing, line-mixing or non-impact data exist. This matters only for comb-resolved or cavity-enhanced Doppler-limited spectra at the 10⁻³ level.
7. **Natural widths.** Hyperfine-dependent predissociation widths (Category 4) are not addressed by any of these studies. Near-dissociation linewidth data (e.g. Goncharov 2007, HWHM 45 kHz at 0.33 Pa) could be combined with pressure extrapolation to separate the two.

### Method caveats

- The OpenAlex daily budget ran out mid-task and Semantic Scholar rate-limited us. As a result, forward-citation searches from Phillips & Perram 2008, Hrabina 2014 (Appl. Opt.), Wolf 2009 and Gläser 1982 were **not** completed.
- Forward citations were completed for Salumbides 2006/2008, Cerny 1986, Gläser 1981, Wu 1985, Fletcher & McDaniel 1995, Hiller & Hanson 1990, Fang 2006 and Sakurai 1979.
- Several 1970s–80s abstracts (Elsevier Opt. Commun. and J. Mol. Spectrosc.) were not retrievable; those entries carry "not extracted" fields.

---

## 6. Gap analysis

This section combines the category-level gap notes above into one view. The question it answers: where would a global potential fit of B–X (with X, B and hyperfine models) be poorly constrained, if it used everything catalogued here plus the Stage 2b sub-Doppler data?

### 6.1 Spectral coverage

| Region (cm⁻¹) | λ (nm) | What exists | Main weakness | Severity |
|---|---|---|---|---|
| < 10 150 | > 985 | No broadband absorption. B–X only via emission: Lyon LIF-FTS down to 7400 cm⁻¹ (1980s); Novosibirsk three-level sub-MHz lines at 982–985 and 1053–1068 nm. A–X absorption (Appadoo 1996). | Absorption is weak here; matters mainly as X-state information at high v″ | medium |
| 10 150–10 930 | 915–985 | Nölleke 2018 only: ≈10 000 unassigned lines at 50 MHz, machine-readable | No assignments, no second source, no intensities | medium–high |
| 10 930–11 000 | 909–915 | Nothing found | Small gap between Nölleke and the atlas | low |
| 11 000–14 000 | 714–909 | Gerstenkorn–Vergès–Chevillard heated-cell atlas (printed; 20–45 MHz after the Knöckel recalibration); sub-Doppler islands (Hannover 778–815 nm, other near-IR metrology in Stage 2b) | Only v′ = 0 bands (0–12 … 0–17) are known at the sub-MHz level; hot bands with v′ > 0 rest on the printed atlas; A–X overlap | high |
| 14 000–15 000 | 667–714 | GL atlas II (printed); Salami & Ross from 14 250 cm⁻¹; Marcassa-group laser absorption (14 400–14 710); Wolf 2009 line shapes at 675 nm | Knöckel model 2σ up to 60 MHz in 667–776 nm; no absolute sub-Doppler grid | high |
| 15 000–19 000 | 526–667 | Best-covered region: GL atlas and Salami & Ross; Kato Doppler-free atlas; Amsterdam 1 MHz grids; many comb-referenced lines (Stage 2b); astronomical FTS (Reiners 2024 comb-calibrated; NIST 0.007 cm⁻¹ scans) | Knöckel 2σ ≈ 3 MHz, but Reiners 2024 shows ≈3.6 MHz systematic, band-correlated residuals | low (positions), medium (MHz-level consistency) |
| 19 000–19 450 | 514–526 | GL atlas and Salami & Ross; dual-comb demonstration (Ideguchi 2012); isolated sub-Doppler lines at 514/515/520 nm (Stage 2b) | Knöckel 2σ ≈ 30 MHz; few systematic data | medium–high |
| 19 450–20 050 | 499–514 | GL atlas and Salami & Ross; sub-Doppler lines near dissociation (Cheng 2002; 501.7 nm; 508 nm, Stage 2b/§5) | Outside the Knöckel fit (B v′ > 43). ESPRESSO team: model "reliable only for λ ≳ 5150 Å". Near-limit B–a 1g hyperfine mixing (Pique 1984/86) needs coupled channels | high |
| > 20 040 | < 499 | Bound–free continuum only (Tellinghuisen 1973; Saiz-Lopez 2004 cross sections) | Intensity-only; no line positions exist | n/a |

### 6.2 Vibrational and rotational coverage

- **X state (term values).**
  - v″ = 0–7: dense (atlas plus sub-Doppler grids).
  - v″ = 8–11: in Knöckel 2004, mostly Gerstenkorn atlas data (≈30 MHz), with no high-precision link to v″ ≤ 7 or to v″ = 12–17.
  - v″ = 12–17: Hannover near-IR sub-Doppler data (v′ = 0 only) plus the heated-cell atlas.
  - v″ = 18–47: only 1980s LIF-FTS (≈30–150 MHz) and unassigned hot-band absorption (Nölleke 2018; atlas I).
  - v″ = 48, 53, 54: a few Novosibirsk sub-MHz anchors, at few J.
  - v″ = 55–107: Martin et al. 1986 only. v″ = 108–113 are extrapolated.
- **B state.**
  - v′ = 0–43: in the Knöckel fit. Dense for v′ ≲ 30; sparse above (Knöckel Fig. 2).
  - v′ = 44–62: atlas and Salami & Ross plus sub-Doppler lines near 500–515 nm. Not in any modern global fit.
  - v′ = 63–82: 1970s–80s long-range analyses. C5 values disagree by 9 % (2.88 vs 3.161 × 10⁵ cm⁻¹ Å⁵).
  - v′ > 82: near-dissociation hyperfine-mixed levels.
- **Rotation.**
  - Knöckel 2004 reaches J ≈ 238, but residuals grow with J (up to 7 MHz at J″ = 238 in the 0–15 band).
  - High-J coverage for v′ > 30 is sparse.
  - LIF data at high v″ give only J′ ± 1 for each pumped level, so J″ coverage there is a sparse set of laser coincidences.
  - Isotopologue data cover J′ = 23–128 only.
  - Inferred: low-J (J ≲ 10) lines in hot bands and near dissociation are rarely measured precisely, because they are weak and blended at room temperature.
- **Isotopologues.** ¹²⁹I₂ and ¹²⁷I¹²⁹I:
  - Precision data exist only at 573–640 nm (v′ ≤ 20, v″ ≤ 5).
  - The X-state isotopologue term values come from one 1986 study.
  - There is no broadband FTS of either species, so the Born–Oppenheimer-correction functions are constrained over a narrow R range.

### 6.3 Intensity, width and line-shape gaps

(See §4 and §5.)

- **Absolute line strengths.** None at Doppler-limited resolution anywhere in 500–680 nm. The Salami & Ross spectrum and the astronomical FTS scans could supply them if cell conditions (I₂ column, temperature) were documented.
- **Temperature dependence.** Nothing modern between 340 and 1300 K, although heated cells and filters operate at 320–600 K.
- **Transition moment.** Measured over R-centroids 2.63–6.04 Å, with poor mutual consistency outside the 520–640 nm window of Tellinghuisen 2011.
- **Natural widths.** Lifetimes and predissociation rates are 1971–1994 data, printed only, sparse in J′ and hyperfine level.
- **Self-broadening.** No measurement in linear absorption (Doppler-limited).
- **Foreign-gas broadening and shift.** Only at 543 nm (air) and 675 nm. No temperature exponents above 388 K. No speed-dependent or line-mixing studies.
- **Vapor pressure.** No direct measurements below about 0 °C (cold-finger range) in the sources checked. A thermodynamic fit will be needed.

### 6.4 Machine readability

- **Directly ingestible now:**
  - Salami & Ross 2005 (spectrum).
  - Nölleke 2018 (line list plus spectrum; CC BY 4.0).
  - The MPI-Mainz UV/VIS atlas files (21 I₂ cross-section sets, including Saiz-Lopez 2004, Spietz 2006 and Tellinghuisen 2011), and the Saiz-Lopez ACP supplement.
- **Electronic supplements that should exist but were not downloaded:** Salumbides et al. 2006 and 2008.
- **Lost or unretrieved:** the Knöckel 2004 supplementary table.
- **Everything else** is printed (atlases, Lyon LIF-FTS term values, lifetimes) or unpublished (astronomical FTS scans, IAG and NIST spectra).

### 6.5 Recommended actions, in priority order

1. **Recover the Knöckel 2004 supplementary table**: ask E. Tiemann (Hannover), B. Bodermann (PTB) or TOPTICA (IodineSpec). If it cannot be recovered, rebuild it from the cited sources, most of which are in this catalog or Stage 2b.
2. **Ingest Salami & Ross 2005 and Nölleke 2018 now**, and write a profile-fitting line extractor that turns the Salami spectrum into positions using a trial model. This gives a machine-readable surrogate for the printed atlas over 14 250–20 050 cm⁻¹.
3. **Digitize** the 11 000–14 000 cm⁻¹ atlas and atlas II (14 000–15 600 cm⁻¹). These are the regions where the model is weakest. Apply the Knöckel recalibration.
4. **Request unpublished FTS spectra**: Reiners et al. 2024 (comb-calibrated, IAG); Nave & Butler 2026 (NIST, 0.007 cm⁻¹); Wang et al. 2019 (NIST temperature series); the Keck HIRES 2009 NIST scan. Use them for position residual maps (as in Reiners 2024) and for intensity-versus-temperature validation.
5. **Full-text checks still needed**:
   - Salami & Ross 2005: cell conditions and calibration.
   - Marcassa group 2022/2023: calibration and uncertainty.
   - Martin et al. 1986: term-value tables.
   - Vigué et al. 1981 (open on HAL): predissociation rates.
   - Tellinghuisen 2011: supplementary data.
6. **Candidate new measurements**, where pyodine could lead:
   - A comb-referenced FTS or dual-comb survey of a heated cell over 10 000–15 000 cm⁻¹.
   - Comb-referenced lines over 19 000–20 050 cm⁻¹ (v′ = 44–80).
   - An FTS spectrum of an enriched ¹²⁹I₂ cell.
   - Doppler-limited self- and N₂-broadening near 500–530 nm versus temperature.

---

## 7. References

References are grouped by the section that discusses them. DOIs were checked against Crossref, OpenAlex or the publisher record.

### Category 1: FTS atlases, broadband and supporting references

- S. Gerstenkorn, P. Luc, *Atlas du spectre d'absorption de la molécule d'iode, 14 800–20 000 cm⁻¹* (Éditions du CNRS, Paris, 1978), 546 pp. ADS 1978adsa.book.....G. Also issued as Laboratoire Aimé Cotton volumes 15 600–17 600 cm⁻¹ (1977) and 17 500–20 000 cm⁻¹ (1977), per Knöckel 2004 ref. [1].
- S. Gerstenkorn, P. Luc, *Atlas du spectre d'absorption de la molécule d'iode, 14 000–15 600 cm⁻¹* (Laboratoire Aimé Cotton, CNRS II, Orsay, 1978; per Knöckel 2004 ref. [1]).
- S. Gerstenkorn, J. Vergès, J. Chevillard, *Atlas du spectre d'absorption de la molécule d'iode, 11 000–14 000 cm⁻¹* (Laboratoire Aimé Cotton, CNRS II, Orsay, 1982).
- S. Gerstenkorn, P. Luc, *Atlas … 14 800–20 000 cm⁻¹. Complément: identification des transitions du système (B–X)* (Éditions du CNRS, Paris, 1980). ISBN 978-2-222-03881-8.
- S. Gerstenkorn, P. Luc, "Absolute iodine (I₂) standards measured by means of Fourier transform spectroscopy," Rev. Phys. Appl. 14, 791–794 (1979). doi:10.1051/rphysap:01979001408079100
- S. Gerstenkorn, P. Luc, "Description of the absorption spectrum of iodine recorded by means of Fourier Transform Spectroscopy: the (B-X) system," J. Phys. (Paris) 46, 867–881 (1985). doi:10.1051/jphys:01985004606086700
- S. Gerstenkorn, P. Luc, J. Sinzelle, "Study of the iodine absorption spectrum by means of Fourier spectroscopy in the region 12 600–14 000 cm⁻¹," J. Phys. (Paris) 41, 1419–1430 (1980). doi:10.1051/jphys:0198000410120141900
- S. Gerstenkorn, P. Luc, R. Perrin, "Rotational analysis of the 5350 Å band of iodine by means of Fourier transform spectroscopy," J. Mol. Spectrosc. 64, 56–69 (1977). doi:10.1016/0022-2852(77)90339-3
- S. Gerstenkorn, P. Luc, "Assignments of several groups of iodine (I₂) lines in the B-X system," J. Mol. Spectrosc. 77, 310–321 (1979). doi:10.1016/0022-2852(79)90111-5
- S. Gerstenkorn, P. Luc, "The iodine absorption spectrum around 5915 Å. Calibration and assignments," Opt. Commun. 36, 322–326 (1981). doi:10.1016/0030-4018(81)90384-9
- S. Gerstenkorn, P. Luc, A. Vetter, "Excitation spectrum of the iodine molecule induced by laser radiation in the 15 780–15 815 cm⁻¹ region: Complement to the Atlas," Rev. Phys. Appl. 16, 529–533 (1981). doi:10.1051/rphysap:01981001609052900
- P. Luc, "Molecular constants and Dunham expansion parameters describing the B-X system of the iodine molecule," J. Mol. Spectrosc. 80, 41–55 (1980). doi:10.1016/0022-2852(80)90269-6
- P. Luc, S. Gerstenkorn, "Fourier transform spectroscopy in the visible and ultraviolet range," Appl. Opt. 17, 1327 (1978). doi:10.1364/AO.17.001327 *(instrument; context only)*
- H. Knöckel, B. Bodermann, E. Tiemann, "High precision description of the rovibronic structure of the I₂ B-X spectrum," Eur. Phys. J. D 28, 199–209 (2004). doi:10.1140/epjd/e2003-00313-4
- B. Bodermann, H. Knöckel, E. Tiemann, "Widely usable interpolation formulae for hyperfine splittings in the ¹²⁷I₂ spectrum," Eur. Phys. J. D 19, 31–44 (2002). doi:10.1140/epjd/e20020052 *(Stage 2b; cited for context)*
- B. A. Palmer, R. A. Keller, Los Alamos Scientific Laboratory Informal Report LA-8251-MS (1980), as cited by Knöckel 2004 *(not an I₂ data set)*
- H. Kato et al., *Doppler-Free High Resolution Spectral Atlas of Iodine Molecule 15 000 to 19 000 cm⁻¹* (Japan Society for the Promotion of Science, 2000) *(Stage 2b)*
- S. Rakowsky, D. Zimmermann, W. E. Ernst, "Accurate determination of wavenumbers for iodine molecular lines in the red spectral region," Appl. Phys. B 48, 463–466 (1989). doi:10.1007/BF00694680
- I. Velchev, R. van Dierendonck, W. Hogervorst, W. Ubachs, "A dense grid of reference iodine lines for optical frequency calibration in the range 571–596 nm," J. Mol. Spectrosc. 187, 21–27 (1998). doi:10.1006/jmsp.1997.7480 *(Stage 2b)*
- S. C. Xu, R. van Dierendonck, W. Hogervorst, W. Ubachs, "A dense grid of reference iodine lines for optical frequency calibration in the range 595–655 nm," J. Mol. Spectrosc. 201, 256–266 (2000). doi:10.1006/jmsp.2000.8085 *(Stage 2b)*
- H. Salami, A. J. Ross, "A molecular iodine atlas in ascii format," J. Mol. Spectrosc. 233, 157–159 (2005). doi:10.1016/j.jms.2005.06.002. Data: https://ars.els-cdn.com/content/image/1-s2.0-S002228520500130X-mmc1.txt
- C. Nölleke, C. Raab, R. Neuhaus, S. Falke, "Absolute frequency atlas from 915 nm to 985 nm based on laser absorption spectroscopy of iodine," J. Mol. Spectrosc. 346, 19–22 (2018). doi:10.1016/j.jms.2017.12.013. Data: Mendeley Data, doi:10.17632/cvr6vv2zps.1
- M. A. Lefrán Torres, D. Rodríguez Fernández, M. R. Cardoso, L. G. Marcassa, "High resolution laser spectroscopy of iodine molecule in the 14 600–14 710 cm⁻¹ range," J. Mol. Spectrosc. 387, 111668 (2022). doi:10.1016/j.jms.2022.111668
- D. Rodríguez Fernández, M. A. Lefrán Torres, M. R. Cardoso, J. D. M. Kondo, L. G. Marcassa, "High resolution laser spectroscopy of iodine molecule in the 14 400–14 600 cm⁻¹ range," J. Mol. Spectrosc. 395, 111789 (2023). doi:10.1016/j.jms.2023.111789
- J. C. Williamson, T. S. Kuntzleman, R. A. Kafader, "A molecular iodine spectral data set for rovibronic analysis," J. Chem. Educ. 90, 383–385 (2013). doi:10.1021/ed300455n
- A. Reiners, M. Debus, S. Schäfer, E. Tiemann, M. Zechmeister, "Accurate calibration spectra for precision radial velocities. Iodine absorption referenced by a laser frequency comb," Astron. Astrophys. 690, A210 (2024). doi:10.1051/0004-6361/202451389 (arXiv:2409.02631)
- M. Debus, S. Schäfer, A. Reiners, "Toward 10 cm s⁻¹ radial velocity accuracy on the Sun using a Fourier transform spectrometer," J. Astron. Telesc. Instrum. Syst. 9(4), 045003 (2023). doi:10.1117/1.JATIS.9.4.045003
- V. Perdelwitz, P. Huke, "A systematic approach to determining the properties of an iodine absorption cell for high-precision radial velocity measurements," Mon. Not. R. Astron. Soc. (2018). doi:10.1093/mnras/sty1523 (arXiv:1806.03040)
- T. Schmidt, A. Reiners, M. Murphy, G. Lo Curto, C. Martins, P. Huke, "Validation of the ESPRESSO wavelength calibration using iodine absorption cell spectra," Mon. Not. R. Astron. Soc. 539, 3301–3318 (2025). doi:10.1093/mnras/staf588
- G. Nave, R. P. Butler, "Design and commissioning of an iodine cell for the ESPRESSO spectrograph," arXiv:2606.13816 (2026)
- S. Wang, J. Wright, P. MacQueen, W. Cochran, D. Doss, C. Gibson, J. Schmitt, "Calibrating iodine cells for precise radial velocities," Publ. Astron. Soc. Pac. 132, 014503 (2020). doi:10.1088/1538-3873/ab5021
- R. Butler, G. Marcy, E. Williams, C. McCarthy, P. Dosanjh, S. Vogt, "Attaining Doppler precision of 3 m s⁻¹," Publ. Astron. Soc. Pac. 108, 500 (1996). doi:10.1086/133755
- J. Valenti, R. Butler, G. Marcy, "Determining spectrometer instrumental profiles using FTS reference spectra," Publ. Astron. Soc. Pac. 107, 966 (1995). doi:10.1086/133645
- G. Marcy, R. Butler, "Precision radial velocities with an iodine absorption cell," Publ. Astron. Soc. Pac. 104, 270 (1992). doi:10.1086/132989
- W. M. Keck Observatory, HIRES iodine cell page: https://www2.keck.hawaii.edu/realpublic/realpublic/inst/hires/iodine_cell.html (accessed Sept 2026)
- NSO McMath FTS archive: https://nispdata.nso.edu/ftp/FTS_cdrom/ (accessed Sept 2026)
- P. Heeren et al., "Pyodine: an open, flexible reduction software for iodine-calibrated precise radial velocities," Astron. Astrophys. 674, A164 (2023). doi:10.1051/0004-6361/202244441
- H. Kambe et al., "Development of iodine cells for the Subaru HDS and the Okayama HIDES. I," Publ. Astron. Soc. Jpn. 54, 865–871 (2002). doi:10.1093/pasj/54.6.865; B. Sato et al., "… II," PASJ 54, 873–882 (2002). doi:10.1093/pasj/54.6.873; H. Kambe et al., "… III," PASJ 60, 45–53 (2008). doi:10.1093/pasj/60.1.45
- J. D. Crane, S. A. Shectman, R. P. Butler, "The Carnegie Planet Finder Spectrograph," Proc. SPIE 6269, 626931 (2006). doi:10.1117/12.672339; status report, Proc. SPIE 7014, 701479 (2008). doi:10.1117/12.789637; integration and commissioning, Proc. SPIE 7735, 773553 (2010). doi:10.1117/12.857792
- M. V. Radovan et al., "A radial velocity spectrometer for the Automated Planet Finder Telescope at Lick Observatory," Proc. SPIE 7735, 77354K (2010). doi:10.1117/12.857726; "The automated planet finder at Lick Observatory," Proc. SPIE 9145, 91452B (2014). doi:10.1117/12.2057310
- A. Tokovinin et al., "CHIRON—A fiber fed spectrometer for precise radial velocities," Publ. Astron. Soc. Pac. 125, 1336–1347 (2013). doi:10.1086/674012
- T. Ideguchi, A. Poisson, G. Guelachvili, T. W. Hänsch, N. Picqué, "Adaptive dual-comb spectroscopy in the green region," Opt. Lett. 37, 4847–4849 (2012). doi:10.1364/OL.37.004847; "Real-time dual-comb spectroscopy of iodine in the visible," CLEO 2012, CW1J.3. doi:10.1364/CLEO_SI.2012.CW1J.3
- A. Eber, C. Gruber, M. Schultze, B. Bernhardt, M. Ossiander, "Streaming self-corrected dual-comb spectrometer," Opt. Express 33, 35314 (2025). doi:10.1364/OE.569404
- A. Eber, M. Pal, L. Fürst, E. Hruska, C. Gruber, M. Ossiander, B. Bernhardt, "NIR/VIS dual-comb spectroscopy comparing high and low repetition rate regimes," Laser Photonics Rev. 20, e02713 (2026). doi:10.1002/lpor.202502713
- A. Nishiyama, D. Ishikawa, M. Misono, "High resolution molecular spectroscopic system assisted by an optical frequency comb," J. Opt. Soc. Am. B 30, 2107 (2013). doi:10.1364/JOSAB.30.002107 *(Stage 2b)*
- S. Vaughan, T. Gherman, A. Ruth, J. Orphal, "Incoherent broad-band cavity-enhanced absorption spectroscopy of the marine boundary layer species I₂, IO and OIO," Phys. Chem. Chem. Phys. 10, 4471 (2008). doi:10.1039/b802618a
- F. Dixneuf, A. A. Ruth, S. Vaughan, R. M. Varma, J. Orphal, "The time dependence of molecular iodine emission from Laminaria digitata," Atmos. Chem. Phys. 9, 823–829 (2009). doi:10.5194/acp-9-823-2009
- O. Johansson, H. Mutelle, A. Parker, S. Batut, P. Demaux, C. Schoemaecker, C. Fittschen, "Quantitative IBBCEAS measurements of I₂ in the presence of aerosols," Appl. Phys. B 114, 421–432 (2014). doi:10.1007/s00340-013-5536-9
- J. N. Forkey, W. R. Lempert, R. B. Miles, "Corrected and calibrated I₂ absorption model at frequency-doubled Nd:YAG laser wavelengths," Appl. Opt. 36, 6729 (1997). doi:10.1364/AO.36.006729
- V. Chan, A. Heyes, D. I. Robinson, J. T. Turner, "Iodine absorption filters for Doppler global velocimetry," Meas. Sci. Technol. 6, 784–794 (1995). doi:10.1088/0957-0233/6/6/015
- Y. Shibata, C. Nagasawa, M. Abo, T. Nagai, "System evaluation of an incoherent wind Doppler lidar using an iodine filter," Jpn. J. Appl. Phys. 48, 032401 (2009). doi:10.1143/JJAP.48.032401
- V. I. Serdyukov, "Frequency scale correction of Fourier spectrometers in the visible," J. Appl. Spectrosc. 83, 307–309 (2016). doi:10.1007/s10812-016-0287-0 *(content not accessed)*

### Category 2: fluorescence, term values, other states


- R. Bacis, S. Churassy, R. W. Field, J. B. Koffend, J. Vergès, "High resolution and sub-Doppler Fourier transform spectroscopy: Iodine molecular fluorescence excited by the 514.5 and 501.7 nm Ar⁺ laser lines," J. Chem. Phys. 72, 34–42 (1980). doi:10.1063/1.438855
- F. Martin, R. Bacis, S. Churassy, J. Vergès, "Laser-induced-fluorescence Fourier transform spectrometry of the XOg⁺ state of I₂: Extensive analysis of the BOu⁺ → XOg⁺ fluorescence spectrum of ¹²⁷I₂," J. Mol. Spectrosc. 116, 71–100 (1986). doi:10.1016/0022-2852(86)90254-7
- D. Cerny, R. Bacis, J. Vergès, "… Extensive analysis of the BOu⁺ → XOg⁺ fluorescence spectrum of ¹²⁷I¹²⁹I and ¹²⁹I₂," J. Mol. Spectrosc. 116, 458–498 (1986). doi:10.1016/0022-2852(86)90140-2
- R. Bacis, D. Cerny, F. Martin, "… Tests of the long-range behavior for three isotopes of iodine," J. Mol. Spectrosc. 118, 434–447 (1986). doi:10.1016/0022-2852(86)90180-3
- F. Martin, S. Churassy, R. Bacis, R. W. Field, J. Vergès, "Long range behavior of the gerade states near the ²P₃/₂ + ²P₃/₂ iodine dissociation limit by LIF-FTS," J. Chem. Phys. 79, 3725–3737 (1983). doi:10.1063/1.446293
- S. Churassy, F. Martin, R. Bacis, J. Vergès, R. W. Field, "Rotation–vibration analysis of the B 0u⁺–a 1g and B 0u⁺–a′ 0g⁺ electronic systems of I₂ by LIF-FTS," J. Chem. Phys. 75, 4863–4868 (1981). doi:10.1063/1.441923
- J. G. Ashmore, J. Tellinghuisen, "Combined polynomial and near-dissociation representations for diatomic spectral data: Cl₂(X) and I₂(X)," J. Mol. Spectrosc. 119, 68–82 (1986). doi:10.1016/0022-2852(86)90202-x
- R. D. Verma, "Ultraviolet Resonance Spectrum of the Iodine Molecule," J. Chem. Phys. 32, 738–749 (1960). doi:10.1063/1.1730793
- R. J. LeRoy, "Spectroscopic Reassignment and Ground-State Dissociation Energy of Molecular Iodine," J. Chem. Phys. 52, 2678–2682 (1970). doi:10.1063/1.1673357
- W. Holzer, W. F. Murphy, H. J. Bernstein, "Resonance Raman Effect and Resonance Fluorescence in Halogen Gases," J. Chem. Phys. 52, 399–407 (1970). doi:10.1063/1.1672699
- W. Kiefer, H. J. Bernstein, "Vibrational-rotational structure in the resonance Raman effect of iodine vapor," J. Mol. Spectrosc. 43, 366–381 (1972). doi:10.1016/0022-2852(72)90048-3
- Yu. A. Matyugin et al., "Study of the hyperfine structure of emission lines of I₂ molecules by the method of three-level laser spectroscopy," Quantum Electron. 38, 755–763 (2008). doi:10.1070/qe2008v038n08abeh013761
- Yu. A. Matyugin et al., "Absolute frequency measurement for the emission transitions of molecular iodine in the 982–985 nm range," Quantum Electron. 42, 250–257 (2012). doi:10.1070/qe2012v042n03abeh014799
- M. I. Nesterenko et al., "Absolute frequency measurements for emission transitions of molecular iodine in the range of 1053–1068 nm," Quantum Electron. 49, 633–640 (2019). doi:10.1070/qel16890
- M. Klug, K. Schulze, U. Hinze, A. Apolonskii, E. Tiemann, B. Wellegehausen, "Frequency stable I₂ Raman laser excited by a cw frequency doubled monolithic Nd:YAG laser," Opt. Commun. 184, 215–223 (2000). doi:10.1016/s0030-4018(00)00921-4
- J. B. Koffend, R. Bacis, R. W. Field, "Continuous wave optically pumped iodine laser," J. Mol. Spectrosc. 77, 202–212 (1979). doi:10.1016/0022-2852(79)90102-4
- J. B. Koffend, F. Wodarczyk, R. Bacis, R. W. Field, "Collisional relaxation of highly excited vibrational levels of the I₂ X state using an I₂ optically pumped laser," J. Chem. Phys. 72, 478–483 (1980). doi:10.1063/1.438874
- R. Barrow, K. Yee (initials as in Crossref), "B ³Π₀₊ᵤ–X ¹Σ⁺g system of ¹²⁷I₂: rotational analysis and long-range potential in the B ³Π₀₊ᵤ state," J. Chem. Soc., Faraday Trans. 2 69, 684–700 (1973). doi:10.1039/f29736900684
- J. I. Steinfeld, J. D. Campbell, N. A. Weiss, "Spectroscopy of I₂ at the dissociation limit," J. Mol. Spectrosc. 29, 204–215 (1969). doi:10.1016/0022-2852(69)90100-3
- M. D. Danyluk, G. W. King, "Energy levels of iodine near the B state dissociation limit," Chem. Phys. 25, 343–351 (1977). doi:10.1016/0301-0104(77)85144-6
- J. W. Tromp, R. J. Le Roy, S. Gerstenkorn, P. Luc, "Reexamination of the I₂ spectrum near the B(³Π₀ᵤ⁺) state dissociation limit," J. Mol. Spectrosc. 100, 82–94 (1983). doi:10.1016/0022-2852(83)90027-9
- J. W. Tromp, R. J. Le Roy, "Near-dissociation expansion representation of large spectroscopic data sets: The B ← X system of I₂," J. Mol. Spectrosc. 109, 352–367 (1985). doi:10.1016/0022-2852(85)90318-2
- S. Gerstenkorn, P. Luc, C. Amiot, "Analysis of the long range potential of iodine in the B ³Π₀ᵤ⁺ state," J. Phys. (Paris) 46, 355–364 (1985). doi:10.1051/jphys:01985004603035500
- S. Gerstenkorn, P. Luc, R. J. Le Roy, "Molecular constants describing the B–X transitions of ¹²⁷,¹²⁹I₂ and ¹²⁹,¹²⁹I₂," Can. J. Phys. 69, 1299–1303 (1991). doi:10.1139/p91-194
- J. P. Pique, F. Hartmann, R. Bacis, S. Churassy, J. B. Koffend, "Hyperfine-induced ungerade-gerade symmetry breaking in a homonuclear diatomic molecule near a dissociation limit: ¹²⁷I₂ at the ²P₃/₂–²P₁/₂ limit," Phys. Rev. Lett. 52, 267–270 (1984). doi:10.1103/PhysRevLett.52.267
- J. P. Pique, F. Hartmann, S. Churassy, R. Bacis, "Hyperfine interactions in homonuclear diatomic molecules and u-g perturbations. I. Theory; II. Experiments on I₂," J. Phys. (Paris) 47, 1909–1916 and 1917–1929 (1986). doi:10.1051/jphys:0198600470110190900; doi:10.1051/jphys:0198600470110191700
- M. Broyer, J.-C. Lehmann, J. Vigué, "g Factors and lifetimes in the B state of molecular iodine," J. Phys. (Paris) 36, 235–241 (1975). doi:10.1051/jphys:01975003603023500
- J. Vigué, M. Broyer, J.-C. Lehmann, "Natural hyperfine and magnetic predissociation of the I₂ B state, I. Theory; II. Experiments on natural and hyperfine predissociation; III. Experiments on magnetic predissociation," J. Phys. (Paris) 42, 937–947; 949–959; 961–978 (1981). doi:10.1051/jphys:01981004207093700; doi:10.1051/jphys:01981004207094900; doi:10.1051/jphys:01981004207096100
- J. Tellinghuisen, "Potentials for weakly bound states in I₂ from diffuse spectra and predissociation data," J. Chem. Phys. 82, 4012–4016 (1985). doi:10.1063/1.448841
- D. R. T. Appadoo et al., "Comprehensive analysis of the A–X spectrum of I₂: An application of near-dissociation theory," J. Chem. Phys. 104, 903–913 (1996). doi:10.1063/1.470814
- J. Tellinghuisen, "On the efficient representation of comprehensive, precise spectroscopic data sets: The A state of I₂," J. Chem. Phys. 118, 3532–3537 (2003). doi:10.1063/1.1539849
- D. Inard, D. Cerny, M. Nota, R. Bacis, S. Churassy, V. Skorokhodov, "E0g⁺→A1u and E0g⁺→B″1u laser-induced fluorescence in molecular iodine recorded by Fourier-transform spectroscopy," Chem. Phys. 243, 305–321 (1999). doi:10.1016/s0301-0104(99)00077-4
- D. Cerny, R. Bacis, S. Churassy, D. Inard, M. Lamrini, M. Nota, "Analysis of the D′2g–A′2u transition in molecular iodine by LIF-FTS," Chem. Phys. 216, 207–226 (1997). doi:10.1016/s0301-0104(97)00008-6
- J. Tellinghuisen, "The D′ → A′ transition in I₂," J. Mol. Spectrosc. 94, 231–252 (1982). doi:10.1016/0022-2852(82)90002-9
- Koffend, Sibai, Bacis, "Collisionally induced optical double resonance in I₂: rotational analysis of the D′(2g)–A′(2u) laser transition," J. Phys. (Paris) 43, 1639–1651 (1982). doi:10.1051/jphys:0198200430110163900
- T. Ridley, K. Lawley, R. Donovan, "Observation of three weakly bound valence states of I₂," J. Chem. Phys. 127 (2007), article number not retrieved. doi:10.1063/1.2795722 *(cited for context only)*
- H. Knöckel, B. Bodermann, E. Tiemann, "High precision description of the rovibronic structure of the I₂ B–X spectrum," Eur. Phys. J. D 28, 199–209 (2004). doi:10.1140/epjd/e2003-00313-4 *(its data field and fixed parameters were read from the article PDF)*

### Categories 3 and 5: isotopologues and line shape


- J. D. Knox, Y.-H. Pao, High-resolution saturation spectra of the iodine isotope ¹²⁹I₂ in the 633-nm wavelength region, Appl. Phys. Lett. 18, 360–361 (1971). doi:10.1063/1.1653696
- J. D. Knox, Y.-H. Pao, Absorption profiles and inverted Lamb dips of I₂ vapor at 633 nm as studied with a He–Ne laser, Appl. Phys. Lett. 16, 129–131 (1970). doi:10.1063/1.1653124
- M. Tesic, Y.-H. Pao, Theoretical assignment of the observed hyperfine structure in the saturated absorption spectra of ¹²⁹I₂ and ¹²⁷I¹²⁹I vapors in the 633 nm wavelength region, J. Mol. Spectrosc. 57, 75–96 (1975). doi:10.1016/0022-2852(75)90043-0
- W. G. Schweitzer Jr. et al., Description, performance, and wavelengths of iodine stabilized lasers, Appl. Opt. 12, 2927–2938 (1973). doi:10.1364/AO.12.002927
- J. A. Magyar, N. Brown, High resolution saturated absorption spectra of iodine molecules ¹²⁹I₂, ¹²⁹I¹²⁷I, and ¹²⁷I₂ at 633 nm, Metrologia 16, 63–68 (1980). doi:10.1088/0026-1394/16/2/001
- P. E. Ciddor, N. Brown, Hyperfine spectra in iodine-129 at 612 nm, Opt. Commun. 34, 53–56 (1980). doi:10.1016/0030-4018(80)90158-3
- K. Dschao, M. Gläser, J. Helmcke, I₂ stabilized He-Ne lasers at 612 nm, IEEE Trans. Instrum. Meas. 29, 354–357 (1980). doi:10.1109/TIM.1980.4314953
- M. Gläser, K. Dschao, H.-J. Foth, Hyperfine structure and fluorescence analysis of enriched ¹²⁹I₂ at the 612 nm wavelength of the He-Ne laser, Opt. Commun. 38, 119–123 (1981). doi:10.1016/0030-4018(81)90212-1
- J. P. Pique, F. Stoeckel, F. Hartmann, Hyperfine structure of ¹²⁷I₂ and ¹²⁹I₂ optical absorption lines, Opt. Commun. 33, 23–25 (1980). doi:10.1016/0030-4018(80)90085-1
- Wu Chengjiu, G. Gaida, J. Bialas, Hyperfine spectrum of a ro-vibrational line of ¹²⁷I¹²⁹I at 612 nm, Metrologia 21, 1–5 (1985). doi:10.1088/0026-1394/21/1/002
- M. Gläser, Identification of hyperfine structure components of the iodine molecule at 640 nm wavelength, Opt. Commun. 54, 335–342 (1985). doi:10.1016/0030-4018(85)90366-9
- V. M. Khodakovskiy et al., Frequency-modulation saturation spectroscopy of molecular iodine hyperfine structure near 640 nm with a diode laser source, Proc. SPIE 7994, 79940L (2010). doi:10.1117/12.882130 (arXiv:0912.3252)
- T. J. Quinn, Practical realization of the definition of the metre, including recommended radiations of other optical frequency standards (2001), Metrologia 40, 103–133 (2003). doi:10.1088/0026-1394/40/2/316
- BIPM, Mise en pratique – recommended radiations, iodine documents (515, 531, 532, 543, 576, 612, 633, 640 nm). https://www.bipm.org/en/publications/mises-en-pratique/standard-frequencies
- D. Cerny, R. Bacis, J. Vergès, LIF Fourier transform spectrometry of the X0g⁺ state of I₂: extensive analysis of the B→X fluorescence spectrum of ¹²⁷I¹²⁹I and ¹²⁹I₂, J. Mol. Spectrosc. 116, 458–498 (1986). doi:10.1016/0022-2852(86)90140-2
- R. Bacis, D. Cerny, F. Martin, LIF Fourier transform spectrometry of the X0g⁺ state of I₂: tests of the long-range behavior for three isotopes of iodine, J. Mol. Spectrosc. 118, 434–447 (1986). doi:10.1016/0022-2852(86)90180-3
- K. Wieland, J. Tellinghuisen, A. Nobs, The band systems E→B and F→X of ¹²⁷I₂ and ¹²⁹I₂…, J. Mol. Spectrosc. 41, 69–83 (1972). doi:10.1016/0022-2852(72)90123-3
- G. W. King, T. D. McLean, The three-photon absorption spectrum of ¹²⁷I₂ and ¹²⁹I₂ in the region 16 500–18 500 cm⁻¹, J. Mol. Spectrosc. 135, 207–222 (1989). doi:10.1016/0022-2852(89)90151-3
- E. J. Salumbides, K. S. E. Eikema, W. Ubachs, U. Hollenstein, H. Knöckel, E. Tiemann, The hyperfine structure of ¹²⁹I₂ and ¹²⁷I¹²⁹I in the B³Π₀u⁺–X¹Σg⁺ band system, Mol. Phys. 104, 2641–2652 (2006). doi:10.1080/00268970600747696
- E. J. Salumbides et al., Improved potentials and Born-Oppenheimer corrections by new measurements of transitions of ¹²⁹I₂ and ¹²⁷I¹²⁹I in the B–X band system, Eur. Phys. J. D 47, 171–179 (2008). doi:10.1140/epjd/e2008-00045-y
- H. Knöckel, B. Bodermann, E. Tiemann, High precision description of the rovibronic structure of the I₂ B–X spectrum, Eur. Phys. J. D 28, 199–209 (2004). (local PDF checked for the isotopologue line list)
- S. V. Kireev et al., Fluorescence of iodine-127 and iodine-129 isotopes excited by radiation of copper vapor lasers (578.2 nm): II, Laser Phys. 23, 105702 (2013). doi:10.1088/1054-660X/23/10/105702; companion papers doi:10.1134/S1054660X1212016X, doi:10.1088/1054-660X/23/7/075701, doi:10.1088/1054-660X/23/10/105701
- S. V. Kireev et al., Fluorescence of iodine-127 and iodine-129 isotopes excited by radiation of copper vapor laser (510.6 nm), Laser Phys. Lett. 11, 075701 (2014). doi:10.1088/1612-2011/11/7/075701
- J.-M. Chartier, S. Picard-Fredin, A. Chartier, International comparison of iodine cells, Metrologia 29, 361–367 (1992). doi:10.1088/0026-1394/29/5/006
- E. N. Wolf, Pressure broadening and pressure shift of diatomic iodine at 675 nm, Ph.D. dissertation, University of Oregon (2009). arXiv:0910.5053
- J. A. Eng, J. L. Hardwick, J. A. Raasch, E. N. Wolf, Diode laser wavelength modulated spectroscopy of I₂ at 675 nm, Spectrochim. Acta A 60, 3413–3419 (2004). doi:10.1016/j.saa.2003.11.044
- D. G. Fletcher, J. C. McDaniel, Collisional shift and broadening of iodine spectral lines in air near 543 nm, J. Quant. Spectrosc. Radiat. Transfer 54, 837–850 (1995). doi:10.1016/0022-4073(95)00105-T
- B. Hiller, R. K. Hanson, Properties of the iodine molecule relevant to laser-induced fluorescence experiments in gas flows, Exp. Fluids 10, 1–11 (1990). doi:10.1007/BF00187865
- B. Hiller, R. K. Hanson, Simultaneous planar measurements of velocity and pressure fields in gas flows using laser-induced fluorescence, Appl. Opt. 27, 33 (1988). doi:10.1364/AO.27.000033
- G. T. Phillips, G. P. Perram, Pressure broadening by argon in the hyperfine resolved P(10) and P(70) (17,1) transitions of I₂ X→B using sub-Doppler laser saturation spectroscopy, J. Quant. Spectrosc. Radiat. Transfer 109, 1875–1885 (2008). doi:10.1016/j.jqsrt.2007.12.011
- P. Cérez, A. Brillet, F. Hartmann, Metrological properties of the R(127) line of iodine studied by laser saturated absorption, IEEE Trans. Instrum. Meas. 23, 526–528 (1974). doi:10.1109/TIM.1974.4314347
- T. Sakurai, S. Iwasaki, T. Oshida, K. Tanaka, Pressure and power broadenings of the saturated absorption lines of iodine at 633 nm, Jpn. J. Appl. Phys. 18, 1199–1200 (1979). doi:10.1143/JJAP.18.1199
- M. Gläser, Frequency shifts at low iodine pressure of ¹²⁷I₂-stabilized He-Ne lasers at 633 nm wavelength, Metrologia 18, 53–58 (1982). doi:10.1088/0026-1394/18/2/002
- M. Gläser, Properties of a He-Ne laser at λ = 612 nm, stabilized by means of an external iodine absorption cell, IEEE Trans. Instrum. Meas. IM-36, 604–608 (1987). doi:10.1109/TIM.1987.6312749
- A. Yu. Nevsky et al., Frequency comparison and absolute frequency measurement of I₂-stabilized lasers at 532 nm, Opt. Commun. 192, 263–272 (2001). doi:10.1016/S0030-4018(01)01190-7
- M. L. Eickhoff, J. L. Hall, Optical frequency standard at 532 nm, IEEE Trans. Instrum. Meas. 44, 155–158 (1995). doi:10.1109/19.377797
- H.-M. Fang, S. C. Wang, J.-T. Shy, Pressure and power broadening of the a10 component of R(56) 32-0 transition of molecular iodine at 532 nm, Opt. Commun. 257, 76–83 (2006). doi:10.1016/j.optcom.2005.07.016
- T. Kobayashi et al., Compact iodine-stabilized laser operating at 531 nm with stability at the 10⁻¹² level and using a coin-sized laser module, Opt. Express 23, 20749 (2015). doi:10.1364/OE.23.020749
- T. Kobayashi et al., Absolute frequency measurements and hyperfine structures of the molecular iodine transitions at 578 nm, J. Opt. Soc. Am. B 33, 725 (2016). doi:10.1364/JOSAB.33.000725
- A. Goncharov, O. Lopez, A. Amy-Klein, F. Du Burck, Absolute frequency measurements for hyperfine structure determination of the R(26) 62-0 transition at 501.7 nm in molecular iodine, Metrologia 44, 275–278 (2007). doi:10.1088/0026-1394/44/5/003
- A. Nishiyama, S. Okubo, T. Kobayashi, A. Kawasaki, H. Inaba, Measurement of transition frequencies and hyperfine constants of molecular iodine at 520.2 nm, J. Opt. Soc. Am. B 41, 2290 (2024). doi:10.1364/JOSAB.531115
- Y. T. Chen et al., Absolute frequency measurement of molecular iodine hyperfine transitions at 554 nm…, Chin. J. Phys. 88, 485–492 (2024). doi:10.1016/j.cjph.2023.12.007
- J. Hrabina et al., Spectral properties of molecular iodine in absorption cells filled to specified saturation pressure, Appl. Opt. 53, 7435 (2014). doi:10.1364/AO.53.007435
- J. Hrabina et al., Comparison of molecular iodine spectral properties at 514.7 and 532 nm wavelengths, Meas. Sci. Rev. 14, 213–218 (2014). doi:10.2478/msr-2014-0029
- J. Hrabina et al., Iodine absorption cells purity testing, Sensors 17, 102 (2017). doi:10.3390/s17010102
- J. Lazar, J. Hrabina, P. Jedlička, O. Číp, Absolute frequency shifts of iodine cells for laser stabilization, Metrologia 46, 450–456 (2009). doi:10.1088/0026-1394/46/5/008
- M. Weller et al., Molecular iodine at 508 nm: hyperfine spectroscopy and operating regimes for optical frequency references, Opt. Express 34, 31567 (2026). doi:10.1364/OE.605121
- Y.-C. Hsiao et al., Absolute frequency measurement of the molecular iodine hyperfine transitions at 548 nm, J. Opt. Soc. Am. B 30, 328 (2013). doi:10.1364/JOSAB.30.000328 (mentioned only)
- F. Cheng et al., Absolute frequency measurement of molecular iodine hyperfine transition at 534 nm, J. Opt. Soc. Am. B 36, 1816 (2019). doi:10.1364/JOSAB.36.001816 (mentioned only)
- J. Zhang, Z. Lu, L. J. Wang, Absolute frequency measurement of the molecular iodine hyperfine components near 560 nm with a solid-state laser source, Appl. Opt. 48, 5629 (2009). doi:10.1364/AO.48.005629 (mentioned only)
- F.-L. Hong et al., Frequency reproducibility of an iodine-stabilized Nd:YAG laser at 532 nm, Opt. Commun. 235, 377–385 (2004). doi:10.1016/j.optcom.2004.02.044 (mentioned only)

### Category 4: intensity, lifetimes, vapor pressure


- Atkinson R. et al., Evaluated kinetic and photochemical data for atmospheric chemistry: Volume III – gas phase reactions of inorganic halogens, Atmos. Chem. Phys. 7, 981–1191 (2007). doi:10.5194/acp-7-981-2007
- Bauer D., Ingham T., Carl S. A., Moortgat G. K., Crowley J. N., Ultraviolet–visible absorption cross sections of gaseous HOI and its photolysis at 355 nm, J. Phys. Chem. A 102, 2857–2864 (1998). doi:10.1021/jp9804300
- Baxter G. P., Hickey C. H., Holmes W. C., The vapor pressure of iodine, J. Am. Chem. Soc. 29, 127–136 (1907). doi:10.1021/ja01956a004
- Baxter G. P., Grose M. R., The vapor pressure of iodine between 50° and 95°, J. Am. Chem. Soc. 37, 1061–1072 (1915). doi:10.1021/ja02170a007
- Berkenblit M., Reisman A., The vapor pressure of iodine in the temperature interval 43°–80°C, J. Electrochem. Soc. 113, 93 (1966). doi:10.1149/1.2423875
- Bhale G. L., Ahmad S. A., Reddy S. G., A study of variation of electronic transition moment of the B–X system of I₂ from its laser-excited fluorescence spectrum, J. Phys. B 18, 645–655 (1985). doi:10.1088/0022-3700/18/4/012
- Brewer L., Tellinghuisen J., Quantum yield for unimolecular dissociation of I₂ in visible absorption, J. Chem. Phys. 56, 3929–3938 (1972). doi:10.1063/1.1677797
- Broyer M., Lehmann J. C., Vigué J., g factors and lifetimes in the B state of molecular iodine, J. Phys. (Paris) 36, 235–241 (1975). doi:10.1051/jphys:01975003603023500
- Broyer M., Vigué J., Lehmann J. C., Direct evidence of the natural predissociation of the I₂ B state through systematic measurements of lifetimes, J. Chem. Phys. 63, 5428–5431 (1975). doi:10.1063/1.431275
- Broyer M., Vigué J., Lehmann J. C., Hyperfine predissociation of molecular iodine, J. Chem. Phys. 64, 4793–4794 (1976). doi:10.1063/1.432067
- Capelle G. A., Broida H. P., Lifetimes and quenching cross sections of I₂(B ³Π₀ᵤ⁺), J. Chem. Phys. 58, 4212–4222 (1973). doi:10.1063/1.1678977
- Cheng W.-Y., Chen L., Yoon T. H., Hall J. L., Ye J., Sub-Doppler molecular-iodine transitions near the dissociation limit (523–498 nm), Opt. Lett. 27, 571–573 (2002). doi:10.1364/OL.27.000571; errata 27, 1076 (2002), doi:10.1364/OL.27.001076
- de Jong W. A., Visscher L., Nieuwpoort W. C., Relativistic and correlated calculations on the ground, excited, and ionized states of iodine, J. Chem. Phys. 107, 9046–9058 (1997). doi:10.1063/1.475194
- Dubé P., Trinczek M., Hyperfine-structure splittings and absorption strengths of molecular-iodine transitions near the trapping frequencies of francium, J. Opt. Soc. Am. B 21, 1113–1126 (2004). doi:10.1364/JOSAB.21.001113
- Fredin-Picard S., A study of contamination in ¹²⁷I₂ cells using laser-induced fluorescence, Metrologia 26, 235–244 (1989). doi:10.1088/0026-1394/26/4/004
- Gillespie L. J., Fraser L. H. D., The normal vapor pressure of crystalline iodine, J. Am. Chem. Soc. 58, 2260–2263 (1936). doi:10.1021/ja01302a050
- Hiller B., Hanson R. K., Properties of the iodine molecule relevant to laser-induced fluorescence experiments in gas flows, Exp. Fluids 10, 1–11 (1990). doi:10.1007/BF00187865
- Honig R. E., Hook H. O., Vapor pressure data for some common gases, RCA Review 21, 360–368 (1960). (no DOI)
- Hrabina J. et al., Spectral properties of molecular iodine in absorption cells filled to specified saturation pressure, Appl. Opt. 53, 7435–7441 (2014). doi:10.1364/AO.53.007435
- Keller-Rudek H., Moortgat G. K., Sander R., Sörensen R., The MPI-Mainz UV/VIS Spectral Atlas of Gaseous Molecules of Atmospheric Interest, Earth Syst. Sci. Data 5, 365–373 (2013). doi:10.5194/essd-5-365-2013
- Klein U., Mastromarino J., Suwaiyan A., Rotational level lifetimes and self-quenching of ¹²⁷I₂ by fluorescence demodulation spectroscopy, Chem. Phys. Lett. 217, 86–90 (1994). doi:10.1016/0009-2614(93)E1366-O
- Koffend J. B., Bacis R., Field R. W., The electronic transition moment of the B 0ᵤ⁺–X ¹Σg⁺ system of I₂ through gain measurements of an I₂ optically pumped laser, J. Chem. Phys. 70, 2366–2372 (1979). doi:10.1063/1.437744
- Kortüm G., Friedheim G., Lichtabsorption und Molekularzustand des Jods in Dampf und Lösung, Z. Naturforsch. A 2, 20–27 (1947). doi:10.1515/zna-1947-0107
- Lamrini M., Bacis R., Cerny D., Churassy S., Crozet P., Ross A. J., The electronic transition dipole moment of the B0ᵤ⁺→X0g⁺ transition in iodine, J. Chem. Phys. 100, 8780–8783 (1994). doi:10.1063/1.466732
- Lukashov S. S., Petrov A. N., Pravilov A. M., Electronic states of iodine molecule and optical transitions between them, in The Iodine Molecule (Springer, 2018), pp. 21–56. doi:10.1007/978-3-319-70072-4_3
- Martínez E., Martínez M. T., Castaño F., Iodine B³Π(0⁺) state predissociation: evaluation of the interaction mechanisms for predissociated levels of the B state, J. Mol. Spectrosc. 128, 554–563 (1988). doi:10.1016/0022-2852(88)90170-1
- Mathieson L., Rees A. L. G., Electronic states and potential energy diagram of the iodine molecule, J. Chem. Phys. 25, 753–761 (1956). doi:10.1063/1.1743043
- McMillan V. (1966), personal communication in Calvert J. G., Pitts J. N., Photochemistry (Wiley, 1966), p. 184. (no DOI)
- Paisner J. A., Wallenstein R., Rotational lifetimes and self-quenching cross sections in the B ³Π₀ᵤ⁺ state of ¹²⁷I₂, J. Chem. Phys. 61, 4317–4320 (1974). doi:10.1063/1.1681737
- Pique J. P., Bacis R., Hartmann F., Sadeghi N., Churassy S., Hyperfine predissociation in the B state of iodine investigated through lifetime measurements of individual hyperfine sublevels, J. Phys. (Paris) 44, 347–351 (1983). doi:10.1051/jphys:01983004403034700
- Rabinowitch E., Wood W. C., The extinction coefficients of iodine and other halogens, Trans. Faraday Soc. 32, 540–546 (1936). doi:10.1039/TF9363200540
- Saiz-Lopez A., Saunders R. W., Joseph D. M., Ashworth S. H., Plane J. M. C., Absolute absorption cross-section and photolysis rate of I₂, Atmos. Chem. Phys. 4, 1443–1450 (2004). doi:10.5194/acp-4-1443-2004
- Sakurai K., Capelle G., Broida H. P., Measurements of lifetimes and quenching cross sections of the B ³Π₀ᵤ⁺ state of iodine using a tunable dye laser, J. Chem. Phys. 54, 1220–1223 (1971). doi:10.1063/1.1674958
- Sander S. P. et al., Chemical Kinetics and Photochemical Data for Use in Atmospheric Studies, Evaluation No. 17, JPL Publication 10-6 (2011). (no DOI)
- Shirley D. A., Giauque W. F., The entropy of iodine. Heat capacity from 13 to 327 K. Heat of sublimation, J. Am. Chem. Soc. 81, 4778–4779 (1959). doi:10.1021/ja01527a005
- Spietz P., Gómez Martín J. C., Burrows J. P., Effects of column density on I₂ spectroscopy and a determination of I₂ absorption cross section at 500 nm, Atmos. Chem. Phys. 6, 2177–2191 (2006). doi:10.5194/acp-6-2177-2006
- Stull D. R., Vapor pressure of pure substances. Organic and inorganic compounds, Ind. Eng. Chem. 39, 517–550 (1947). doi:10.1021/ie50448a022
- Sulzer P., Wieland K., Intensitätsverteilung eines kontinuierlichen Absorptionsspektrums in Abhängigkeit von Temperatur und Wellenzahl, Helv. Phys. Acta 25, 653–676 (1952). doi:10.5169/seals-112328
- Suwaiyan A., Mastromarino J., Klein U., Absorption, fluorescence and saturation of ¹²⁷I₂ in the region 16999.100–17001.300 cm⁻¹, Chem. Phys. Lett. 192, 581–589 (1992). doi:10.1016/0009-2614(92)85520-K
- Teichteil C., Pelissier M., Relativistic calculations of excited states of molecular iodine, Chem. Phys. 180, 1–18 (1994). doi:10.1016/0301-0104(93)E0395-C
- Tellinghuisen J., Spontaneous predissociation in I₂, J. Chem. Phys. 57, 2397–2402 (1972). doi:10.1063/1.1678600
- Tellinghuisen J., Resolution of the visible–infrared absorption spectrum of I₂ into three contributing transitions, J. Chem. Phys. 58, 2821–2834 (1973). doi:10.1063/1.1679584
- Tellinghuisen J., Continuous absorption below the band convergence limit in the I₂ B←X transition, J. Chem. Phys. 59, 849–852 (1973). doi:10.1063/1.1680103
- Tellinghuisen J., Intensity factors for the I₂ B↔X band system, J. Quant. Spectrosc. Radiat. Transfer 19, 149–161 (1978). doi:10.1016/0022-4073(78)90074-2
- Tellinghuisen J., Transition strengths in the visible–infrared absorption spectrum of I₂, J. Chem. Phys. 76, 4736–4744 (1982). doi:10.1063/1.442791
- Tellinghuisen J., Potentials for weakly bound states in I₂ from diffuse spectra and predissociation data, J. Chem. Phys. 82, 4012–4016 (1985). doi:10.1063/1.448841
- Tellinghuisen J., The electronic transition moment function for the B0ᵤ⁺(³Π)↔X ¹Σg⁺ transition in I₂, J. Chem. Phys. 106, 1305–1308 (1997). doi:10.1063/1.473971
- Tellinghuisen J., Intensity analysis of overlapped discrete and continuous absorption by spectral simulation: the electronic transition moment for the B–X system in I₂, J. Chem. Phys. 134, 084301 (2011). doi:10.1063/1.3555623
- Tellinghuisen J., Least-squares analysis of overlapped bound-free absorption spectra and predissociation data in diatomics: the C(¹Πᵤ) state of I₂, J. Chem. Phys. 135, 054301 (2011). doi:10.1063/1.3616039
- Tench R. E., Ezekiel S., Precision measurements of hyperfine predissociation in I₂ vapor using a two-photon resonant scattering technique, Chem. Phys. Lett. 96, 253–258 (1983). doi:10.1016/0009-2614(83)80502-8
- Vigué J., Broyer M., Lehmann J. C., Predissociation effects in the B ³Π₀₊ᵤ state of iodine, J. Chem. Phys. 62, 4941–4947 (1975). doi:10.1063/1.430409
- Vigué J., Broyer M., Lehmann J. C., Natural hyperfine and magnetic predissociation of the I₂ B state I–III, J. Phys. (Paris) 42, 937–947; 949–959; 961–978 (1981). doi:10.1051/jphys:01981004207093700; 10.1051/jphys:01981004207094900; 10.1051/jphys:01981004207096100
- Vogt K., Koenigsberger J., Beobachtungen über Absorption in Joddampf und anderen Dämpfen, Z. Phys. 13, 292–311 (1923). doi:10.1007/BF01328221
- Zaitsevskii A., Pazyuk E. A., Stolyarov A. V., Teichteil C., Vallet V., Theoretical spectroscopy of molecular iodine. 1. Ab initio study on the B0ᵤ⁺–X0g⁺, A1ᵤ–X0g⁺ and B′1ᵤ–X0g⁺ radiative transition intensities, Mol. Phys. 98, 1973–1979 (2000). doi:10.1080/00268970009483400
- Zare R. N., Calculation of intensity distribution in the vibrational structure of electronic transitions: the B ³Π₀ᵤ⁺–X ¹Σ₀g⁺ resonance series of molecular iodine, J. Chem. Phys. 40, 1934–1944 (1964). doi:10.1063/1.1725425
