# Software for calculating I₂ spectra: a survey

*Survey date: 2026-09-14. Scope: the B³Π(0u⁺)–X¹Σg⁺ system of ¹²⁷I₂ including hyperfine structure, and other isotopologues.*

**Tagging convention.** **[S]** means a source states it; the source is given inline or in the references. **[I]** means an inference from the evidence. "Unknown" means nothing reliable was found. No any price, license, version or accuracy figure.

---

## 0. Bottom line

* **No suitable open-source tool exists.** Nothing public computes physics-based B–X line positions, hyperfine structure and intensities over the visible/NIR range, with uncertainties.
* The only tool that does this is **IodineSpec** (Hannover: Knöckel, Bodermann, Tiemann). It is **closed source**.
  * IodineSpec 4 was sold by TOPTICA from about 2002 to about 2008–09.
  * IodineSpec 5 has since been distributed informally "on request" from Prof. Tiemann.
  * The last published model updates are from **2006–2010**.
* The other I₂-specific codes are either narrow or data-only:
  * The Forkey (1997) Fortran model covers only the ~532 nm region.
  * Undergraduate teaching scripts.
  * Iodine-cell radial-velocity pipelines, which use *measured* cell spectra rather than computing them.
  * Atlases (Gerstenkorn–Luc, Kato, Salami–Ross), which are data, not models.
* The general diatomic codes supply building blocks and benchmarks, but none ships an I₂ model:
  * Le Roy's LEVEL/dPotFit
  * Duo, which since 2026 can do hyperfine-resolved spectra of *homonuclear* diatomics
  * PGOPHER
  * PyDiatomic
  * MARVEL
* **Recommendation: GO.** Build our own package. Use IodineSpec (via published parameters and line lists), LEVEL and Duo as benchmarks. Borrow algorithms from Knöckel 2004, Le Roy, Pashov and MARVEL.
* The ExoMol/MARVEL group keeps a 330-entry `I2_MARVEL.bib` bibliography (uploaded 2024-09-03), which suggests a MARVEL analysis of I₂ is in preparation [I].
* **Name.** The repository began as `pyodine`, which is already the name of an iodine-cell radial-velocity package (Heeren et al. 2023) and is one letter from Pyodide; the package is therefore named `i2spec`.

---

## 1. Method and coverage

| Channel | What was done | Result |
|---|---|---|
| Forward citations | OpenAlex `cites:` and Semantic Scholar for Bodermann 2002 (53 / 46 citing works) and Knöckel 2004 (87 / 68) | Full lists reviewed by title. Code-relevant items are discussed below. |
| Full-text index | OpenAlex `fulltext.search:iodinespec` | 26 works, 2002–2025, all users of IodineSpec. No re-implementations. |
| GitHub | Repo search: "iodine spectrum", "iodine hyperfine", "I2 hyperfine", "iodine absorption", "IodineSpec", "iodine atlas", "iodine reference", "I2 line list", "iodine cell", "iodine laser", "molecular iodine"; code search: "IodineSpec", "Bodermann", "Knoeckel iodine", "127I2" | ≤4 hits per repo query. Only the Forkey model (Fortran, plus a MATLAB port) computes an I₂ spectrum. The code hits were text corpora and bibliographies, apart from ExoMol's I₂ bib files (§4.6). |
| GitLab | API project search "iodine", "pyodine", "diatomic" | Only `Heeren/pyodine` is relevant. |
| PyPI | Full simple index (891,702 names) grepped for iodine/I2/diatomic/RKR/Numerov/hyperfine | No I₂ spectrum package. Related packages: `diatomic-py`, `numerov`, `pyexocross`, `molspecutils`, `diatomics`. |
| crates.io | Search "iodine", "diatomic" | Nothing relevant. |
| Julia General registry | `Registry.toml` grepped | No iodine or diatomic-rovibronic package. A generic `Numerov` package exists (not evaluated). |
| Zenodo | "iodine spectrum", "iodine hyperfine", "iodine atlas", "IodineSpec", "I2 line list" | 0 hits each. |
| MATLAB File Exchange | Web search only | Nothing I₂-specific found. |
| Vendor history | Wayback Machine CDX for toptica.com | IodineSpec 4 product pages and brochure found, 2002–2009 (§2.2). |

**Limitations.**
* The Internet Archive went offline part-way through, so captures of the TOPTICA pages after 2009 could not be enumerated.
* OpenAlex rate-limited late searches.
* MATLAB File Exchange was covered only through web search.
* Full texts not consulted: Salami & Ross 2005, and the Kato atlas explanation PDF (by e-mail only).

---

## 2. IodineSpec (Hannover / TOPTICA): the baseline

### 2.1 Underlying model — Knöckel, Bodermann & Tiemann, EPJ D 28, 199 (2004)

Read from the full text.

**Global model [S]**
* The X and B states are represented by potentials in the "X-representation":
  * V(R) = Σ aᵢ Xⁱ, with X = (R − R_m)/(R + b R_m).
  * An exponential inner wall and a dispersion-plus-exchange long-range tail are attached outside the region covered by data.
* The B state uses the Herman–Asgharian effective Hamiltonian:
  * R-dependent nonadiabatic α(R) and adiabatic V_corr(R) Born–Oppenheimer corrections.
  * β(R) is neglected.
  * The BOCs are fitted using a handful of ¹²⁹I₂ and ¹²⁷I¹²⁹I lines.
* Eigenvalues come from the Numerov method (Blatt 1967). χ² was minimised with MINUIT.
* The input was about 1500 absolute rovibronic frequencies (1513 observations):
  * v″ ≤ 17 and v′ up to 43; J up to 238.
  * Hyperfine shifts were removed first, using the Bodermann 2002 interpolation formulae.
  * The minimum assigned uncertainty was 3 MHz.
* All potential parameters are printed in Table 4. The Hannover IQO web page notes misprints in Table 4 and links a corrected version.

**Local NIR model [S]**
* A Dunham-type fit to bands 0–12 … 0–17 (778–815 nm, J″ ≤ 242), with σ_fit = 40 kHz.

**Hyperfine [S]**
* Bodermann, Knöckel & Tiemann, EPJ D 19, 31–44 (2002): interpolation formulae for the ¹²⁷I₂ hyperfine parameters. Knöckel 2004 describes them as capable of a 2σ uncertainty of "few 10 kHz".
* Later measurements show much larger deviations for some lines, up to 0.8 MHz near 514 nm (Yoshiki et al. 2023; §2.3).

**Data fitted [S]**

| Source | Coverage and weight |
|---|---|
| Gerstenkorn–Luc FTS atlas | Recalibrated and line-shape-corrected (Tables 2–3). Main source for 514–526 nm and 667–776 nm, at ~30 MHz. |
| Kato Doppler-free atlas | 526–667 nm, ~400 averaged frequencies, ~3 MHz. |
| Published grids and precision lines | ~1 MHz grids at 560–650 nm, lines near 532 nm and 514 nm. |
| Hannover NIR measurements | 778–815 nm, < 80 kHz. |

The paper was received in May 2003 [S], so the data end around 2002–03 [I].

**Stated accuracy, Knöckel 2004 [S]**

| Range | Stated uncertainty |
|---|---|
| 526–667 nm | 2σ < 3 MHz |
| 514–526 nm | about ±30 MHz (2σ) |
| 667–776 nm | "partly worse", with 60 MHz as the 2σ upper limit |
| > 776 nm | better than 1.5 MHz for J″ < 120; residuals up to ~7 MHz at J″ = 238 (0–15 band) |
| Local NIR model | 1σ < 200 kHz, 778–815 nm |
| Hyperfine components | 1σ "few MHz" |

**Program availability [S]**, from the conclusion and ref. [54]:
> "The potential and the hyperfine structure calculations have been built together in a computer program [54] … to calculate the frequencies of iodine transitions and to construct the iodine spectrum in line shape, width and intensity appropriate to the experimental conditions, if Doppler broadened or Doppler reduced."
>
> "[54] Such program is available under the name "IodineSpec" from TOPTICA Corp., www.toptica.com"

### 2.2 Versions and distribution

**IodineSpec 4 (sold by TOPTICA) [S: archived TOPTICA pages]**
* Product pages were captured from 2002-07 to 2004-08 (`/products/iodinespec/…`).
* A later product page shows "Latest update: Nov 17, 2008".
* The brochure `BR-32014A-Iodine-Calculation-Software.pdf` was captured 2005-11 to 2009-01. The archived PDF was corrupt and could not be read.

What the pages state:

| Item | TOPTICA page |
|---|---|
| Computes | "the absorption spectrum of the iodine molecule, taking into account the hyperfine structure of each rotational transition" |
| User settings | Temperature and hyperfine linewidth |
| Displays | Doppler-broadened and hyperfine spectra |
| Output | Text file with quantum-number assignments; HPGL plots |
| Accuracy | "The overall uncertainty of the calculated line frequencies is less than 25 MHz in the range from 515 nm to 815 nm. A high precision prediction range is given from 776 nm to 815 nm with an uncertainty for selected lines of less than 0,2 MHz." |
| Data | Gerstenkorn–Luc atlas, Hannover NIR measurements, literature |
| Platform | Windows 95/98/NT/2000 |

The pages list references only up to Bodermann et al. 2000, and the text did not change between 2002 and 2008. So IodineSpec 4 probably predates the 2002 hyperfine formulae and the 2004 global potentials [I]. However, the 2004 paper points to TOPTICA as the distributor [S].

**Price:** unknown. No archived page shows one.

**Current TOPTICA status:** a site-restricted web search of toptica.com finds no IodineSpec page today. The product appears to be discontinued [I].

**IodineSpec 5 [S]**
* Users report getting it directly from Hannover:
  * Huang et al., Appl. Opt. 57, 2102 (2018), arXiv:1710.09533: "We have used IodineSpec version 5 … For the actual status of the program, contact knoeckel@iqo.uni-hannover.de"; "We also thank Dr. H. Knöckel for providing the IodineSpec5 software". It quotes a 2σ of ±3 MHz at 647 nm.
  * The Hannover IQO page: "A simulation program called "IodineSpec" containing these models was developed for convenient calculation of calibration lines for the natural isotopomer ¹²⁷I₂ and for the other isotopologues [¹²⁹I₂, ¹²⁷I¹²⁹I] … A present version of the program can be made available on request to Prof. Tiemann."
* License terms are unknown. No source code has been found publicly [S/I].

**Newer versions:** none found. Searches for "IodineSpec6", "IodineSpec 6" and "IodineSpec version 6" returned nothing.

### 2.3 Model updates after 2004 (searched 2004–2026)

Sources: OpenAlex author records, the forward-citation lists, and the Tiemann publication list on the IQO site.

| Year | Paper | What changed [S] |
|---|---|---|
| 2006 | Salumbides, Eikema, Ubachs, Hollenstein, Knöckel, Tiemann, *Mol. Phys.* 104, 2641 | Hyperfine structure of ¹²⁹I₂ and ¹²⁷I¹²⁹I, with models "giving reliable predictions for the hyperfine parameters of all isotopomers". |
| 2008 | Same authors, *EPJ D* 47, 171 | Improved potentials and Born–Oppenheimer corrections from new ¹²⁹I₂ and ¹²⁷I¹²⁹I data. |
| 2010 | Liao, Wu, Lien, Knöckel, Chui, Tiemann, Shy, *JOSA B* 27, 1208 | Comb measurements of 27 hyperfine transitions at 750–780 nm. They differ from the 2004 model "larger than expected", and "an improved model is developed for the range from 755 to 815 nm … error limit … less than 0.2 MHz". |
| 2024 | Reiners, Debus, Schäfer, Tiemann, Zechmeister, *A&A* 690, A210 | No new model; an external test of IodineSpec (below). |

**The 2024 test (Reiners et al.) [S].** IodineSpec line lists were extracted with the program: 4,427,241 lines in 5150–6300 Å at 44 °C. They were compared with a laser-frequency-comb-referenced FTS spectrum:
* The model's stated accuracy is "better than 3 MHz (∼2 m s⁻¹) in the wavelength range of 5260–6670 Å".
* The observations match the model within about 1 m s⁻¹ overall.
* There is, however, a systematic pattern: "one long-period wave overlaid by a rapidly oscillating pattern that coincides with the I₂ absorption band structure". The long-period wave has an amplitude of about 2 m s⁻¹ (≈3 MHz).
* The authors conclude that improvements "should mainly focus on the energy scale (potential functions and hyperfine interaction)".

This is direct evidence of residual model systematics that a new global fit could remove [I].

**Other third-party accuracy statement.** Huang et al., Appl. Opt. 52, 1448 (2013), 671 nm, cites IodineSpec as better than 2 MHz for 19,000–15,000 cm⁻¹ and about 30 MHz for 667–755 nm, "because precision data are lacking there".

**Hyperfine interpolation formulae tested near 514 nm.** Yoshiki et al., EPJ D 77, 140 (2023), read from the accepted manuscript in the YNU repository [S]:
* They measured hyperfine splittings near 514 nm to about 1 kHz and fitted a four-term Hamiltonian.
* They compared the measurements with predictions from the Bodermann 2002 interpolation formulae:

| Line | Obs − Pred (Bodermann 2002 formulae) |
|---|---|
| P(34)44-0 | SD 78.2 kHz |
| R(58)45-0 (J′ = 59, v′ = 45) | SD 386.7 kHz; individual components up to −802.8 kHz (b5), +612.6 kHz (b2) and −578.5 kHz (b9) (Table 8) |

* The authors say the formulae "provide only an approximate depiction of the hyperfine structures with an uncertainty of tens of kilohertz", yet the observed deviations are far larger.
* **Implication [I]:** the hyperfine model in IodineSpec can be off by nearly 1 MHz for some lines at high v′ and J′. A modern refit of the hyperfine parameters against the ~2005–2024 comb-era data is part of beating the baseline.

**Did these updates reach IodineSpec?**
* The 2006 and 2008 papers **did** change the model inside the program [S: IQO page]. The page says the 2004 potential model "was recently extended in cooperation with the Laser Centre of the VU Amsterdam to cover also the isotopologues ¹²⁹I₂ and ¹²⁷I¹²⁹I", citing the 2006 and 2008 papers. It also says "A simulation program called "IodineSpec" containing these models was developed … for the natural isotopomer ¹²⁷I₂ and for the other isotopologues".
* The 2008 fit refits the potentials and Born–Oppenheimer corrections using isotopologue data [S: title]. Its full text was not consulted for this survey, so how much the ¹²⁷I₂ B-state parameters changed is unknown.
* Whether the 2010 NIR model is in IodineSpec 5 is unknown.

**Maintenance status [S/I].** Tiemann is the listed contact [S: IQO page]. No global ¹²⁷I₂ refit has been published since 2004; the only updates are the 2008 isotopologue BOCs and the 2010 NIR local model [I, from the publication searches]. It is unknown whether the 2010 NIR model is included in IodineSpec 5.

**Who uses it [S: citation and full-text lists]:**
* Metrology groups in Taiwan (Shy / L.-B. Wang), Japan (Hong group), China and Europe.
* Yb, Cd and Ra⁺ laser-cooling experiments.
* Astronomical spectrograph calibration: Göttingen, ESPRESSO.
* Lidar: Kühlungsborn RMR lidar, found by the full-text match only.

### 2.4 Fit to our goals

**What IodineSpec defines [I].** It sets the feature bar:
* ¹²⁷I₂ plus two isotopologues
* hyperfine structure
* temperature-dependent relative intensities and Doppler/hyperfine line shapes
* 514–815 nm coverage

It also sets the accuracy bar: < 3 MHz (2σ) at 526–667 nm, and sub-MHz for selected NIR lines.

**Why it is not a solution for us [I]:**
* It is closed and has no public license.
* It cannot be redistributed or extended.
* It is Windows-only (v4 at least).
* It has no published uncertainty per line.
* It depends on one retired group.

**How to use it.** Treat it as a benchmark, not a component:
1. Reproduce Knöckel 2004 Table 4 (with the erratum) in our code.
2. Compare with IodineSpec 5 line lists, if Prof. Tiemann agrees to share them.

---

## 3. Comparison table

Status is the last release, commit or update seen on 2026-09-14.

| Tool | Maintainers | Language | License | Status | What it computes | I₂-specific? | Fit to our goals |
|---|---|---|---|---|---|---|---|
| **IodineSpec 4/5** | Knöckel† (2024), Bodermann, Tiemann (LU Hannover); v4 sold by TOPTICA | Unknown (Windows program) | Proprietary; v5 "on request"; price unknown | Model updates 2006–2010; no new version found | B–X line positions, hyperfine, relative intensities, Doppler/sub-Doppler spectra; ¹²⁷I₂, ¹²⁹I₂, ¹²⁷I¹²⁹I | Yes | **Baseline to beat.** Not reusable (closed). |
| **Forkey I₂ absorption model** (`i2lines2`, `i2spec4`) | Forkey, Lempert, Miles (Princeton 1997); GitHub copy (Lightning469) and MATLAB port (PriybratS) | Fortran; MATLAB port | No license stated | Copy 2017; port 2026 | Line list plus cell transmission vs. wavenumber (T, p, length). README covers 18786–18791 cm⁻¹ near 532 nm and notes a hand-corrected line position | Yes | Low. Narrow range, legacy constants. Useful only as a 532 nm sanity check. |
| **Kato et al. Doppler-free atlas** (book + CD-ROM) | H. Katô et al. (JSPS 2000) | Data plus viewer CD (PC/Mac) | Commercial book; terms unknown | 2000 | Doppler-free spectra with etalon marks, 15,000–19,000 cm⁻¹; absolute wavenumbers of all hyperfine components | Yes (data) | **Key validation data** (the data compilation), not software. |
| **Salami & Ross ASCII atlas** | H. Salami, A. J. Ross (Lyon) | ASCII data | Unknown (JMS article) | 2005 | Doppler-limited FTS absorption spectrum | Yes (data) | Validation data. Range and resolution not verified (paper paywalled). |
| **Gerstenkorn–Luc atlases**; **Simmons & Hougen NBS atlas** | CNRS Orsay (1977–82); NBS (1977) | Printed | n/a | Historical | FTS line lists 11,000–20,000 cm⁻¹; 18,000–19,000 cm⁻¹ | Yes (data) | Validation data; they underpin the Hannover fit. |
| **Nölleke et al. NIR atlas** | Nölleke, Raab, Neuhaus, Falke (JMS 2018) | Data | Journal | 2018 | "almost 10,000 … reference lines with an uncertainty of 50 MHz", 915–985 nm | Yes (data) | Extends the range beyond 815 nm. |
| **pyodine (Heeren et al.)** | P. Heeren, R. Tronsgaard, F. Grundahl | Python 3 | MIT | GitLab active to 2024-02; forks 2025–26 | Precise radial velocities from I₂-cell échelle spectra (Butler et al. method) | Uses I₂ spectra but does not compute them | Not a spectrum calculator. **Name collision.** |
| **JOKARUS "pyodine"** | F. Gutsch (HU Berlin) | Python | Apache-2.0 | 2018 | Control software for an I₂ frequency-reference experiment | No | Irrelevant, except as a name collision. |
| **LEVEL16** (also RKR1-16, dParFit16, dPotFit16, betaFIT16, BCONT) | R. J. Le Roy† (Waterloo); reproduction hosted by U. Waterloo ACE | Fortran (not verified) | Free, but "may not be sold or have any other commercial use made of it without the express written permission of the author" — not OSI | LEVEL last updated 2016-08-11 | Bound/quasibound eigenvalues, expectation values, FCFs and matrix elements; RKR; potential fitting (MLR etc.) | No | **Benchmark and algorithm source.** Do not vendor (license). |
| **PGOPHER** | C. M. Western (Bristol) | Not verified | LGPL (Bristol dataset record) | Last formal release 10.1.180 (Nov 2018); "further updates … will no longer be available" (2022) | Effective-Hamiltonian simulation and fitting of rotational, vibrational and electronic spectra, including nuclear hyperfine | No | Cross-check for hyperfine Hamiltonians and line strengths; frozen. |
| **Duo** | S. Yurchenko, J. Tennyson et al. (UCL) | Fortran 2003 | "Standard CPC licence" (CPC paper); no LICENSE file on GitHub | Active (push 2026-09-08) | Coupled rovibronic Schrödinger equation for diatomics, intensities, refinement fits; hyperfine (2022 heteronuclear; **2026 homonuclear**) | No | **Strongest general cross-validation engine.** License unclear. |
| **ExoCross / PyExoCross** | Yurchenko et al.; Zhang et al. | Fortran 2003 / Python | MIT / GPL-3.0 | Active (2026-09) | Cross sections and spectra from ExoMol-format line lists (Voigt etc.) | No | Could reuse ExoCross (MIT) for cross-section generation, or adopt the ExoMol format. |
| **PyDiatomic** | S. Gibson | Python | GPL-3.0 | v0.3 (2025-02); last default-branch commit 2023-01 | Coupled-channel Schrödinger equation (Johnson renormalized Numerov); photodissociation cross sections | No | Algorithm reference. GPL limits code reuse. |
| **IPA** (Pashov, Jastrzębski, Kowalczyk) | Pashov et al. | Not verified | CPC licence (CPC Library ADLV_v1_0) | 2000 | Pointwise potentials fitted by the inverted perturbation approach | No | Algorithm reference; Knöckel 2004 cites it as an alternative representation. |
| **MARVEL / MARVEL Online** | Furtenbacher, Császár, Tennyson | Web app | Unknown | Online app (unreachable from here on 2026-09-14) | Spectroscopic-network inversion: transitions → empirical energy levels with uncertainties | No, but an **I₂ bibliography exists** (§4.6) | Very relevant to the data compilation and fitting; coordinate with the group. |
| **SPFIT/SPCAT** (Pickett) and wrappers (Pyckett, PySpecTools, C++ port) | H. M. Pickett (JPL); community | C; Python; C++ | Original unknown; wrappers MIT | Wrappers active 2026 | Effective-Hamiltonian fitting with spin and hyperfine terms | No | Marginal: rotational-spectroscopy oriented. |
| **diatomic-py** | Blackmore, Gregory, Hutson, Cornish (Durham) | Python | BSD-3-Clause | v2.1.0 (2026-05) | Rotational plus hyperfine plus Stark/Zeeman structure of ¹Σ molecules | No | Permissive reference for hyperfine matrix elements (X state). |
| **numerov** (PyPI) | R. Bast | Python | MPL-2.0 | 0.5.0 (2017) | Numerov–Cooley vibrational levels | No | Trivial; our own solver will supersede it. |
| **MPI-Mainz UV/VIS Spectral Atlas** (I₂ entry) | Keller-Rudek, Moortgat, Sander, Sörensen | Data | Terms on site | Maintained | Low-resolution I₂ cross sections, 170–800 nm plots | Data | Atmospheric context and a low-resolution intensity check only. |
| **HITRAN / ExoMol databases** | — | — | — | Checked 2026-09-14 | HITRAN has HI and CH₃I but **no I₂**; the ExoMol master file (v20260605, 102 molecules) has **no I₂ dataset** | — | Gap confirmed. |

---

## 4. Other I₂-specific tools and digital resources

### 4.1 Digital atlases

**Kato et al. (2000), *Doppler-Free High Resolution Spectral Atlas of Iodine Molecule 15000–19000 cm⁻¹*** [S: atlas web page]
* Publisher and format: Japan Society for the Promotion of Science, 4000 pp., ISBN 4-89114-000-3.
* Each 0.52 cm⁻¹ panel gives etalon marks, the Doppler-broadened excitation spectrum and the Doppler-free spectrum.
* "The absolute wavenumber of a transition line in this region can be obtained with a standard deviation of 0.000054 cm⁻¹", which is ≈ 1.6 MHz.
* "This spectral data is available on an accompanying CDROM (CDROM of the PC version and Mac version are attached). Using this CDROM, it is possible to obtain the absolute wavenumbers of all the hyperfine components … approximately 0.0001 cm⁻¹", which is ≈ 3 MHz.
* The CD's software format and license are unknown. An explanatory PDF is available by e-mail from the authors.
* **It did ship with software/data on CD-ROM [S].**

**Salami & Ross (2005), "A molecular iodine atlas in ascii format"**, J. Mol. Spectrosc. 233, 157–159.
* A Doppler-limited FTS absorption spectrum distributed as ASCII.
* The article is not open access and was not consulted. Search-engine summaries give roughly 14,250–20,000 cm⁻¹ at 0.02 cm⁻¹; these figures are **unverified**.
* No accompanying software found.
* The same author later published a Te₂ digital atlas (Ross & Cardon, JMS 384, 111589, 2022) [S].

**Other atlases:**
* Gerstenkorn & Luc (CNRS, 1977–1982), with their 1979 absolute standards paper (Rev. Phys. Appl. 14, 791).
* Simmons & Hougen, NBS atlas 18,000–19,000 cm⁻¹ (J. Res. NBS 81A, 25, 1977; open via PMC).
* Nölleke et al. (2018), 915–985 nm NIR atlas.
* All are data; none has an associated modeling code [I].

### 4.2 Forkey–Lempert–Miles model (filtered Rayleigh scattering community)

*Appl. Opt.* 36, 6729 (1997): "a computer model for accurately predicting absorption profiles for molecular iodine cells over the tuning range of frequency-doubled Nd:YAG lasers" [S].

| Item | Detail |
|---|---|
| Fortran source | Posted on GitHub in 2017 by a third party (`Lightning469/Forkey-Iodine-Cell-Transmission-Code`), with no license |
| Files | `i2lines2.f`, `i2spec4.f`, Franck–Condon file `fcfioded`, a pre-computed line file for 18786–18791 cm⁻¹ |
| Known fix | README notes the predicted position of one line had to be shifted by 0.0815 cm⁻¹ |
| MATLAB port | `PriybratS/Iodine_absorption_Forky_Matlab` (2026) |
| Similar code | `eckeratvt/Iodine-Cell-Transmission-Spectra` (2012; MATLAB; no README or license) |
| Hyperfine treatment | Not verified |

Fit [I]: not a general B–X model. It is the de facto open tool for 532 nm iodine filters, so a 532 nm cell-transmission mode in our package would serve that community.

### 4.3 In-house codes from metrology groups

Every precision-measurement group computes hyperfine patterns. Typically they fit measured splittings to the four-term effective hyperfine Hamiltonian (eqQ, C, d, δ) and/or use IodineSpec for predictions. Representative work:

| Group | Papers | What they do [S] |
|---|---|---|
| JILA | Chen & Ye, *Chem. Phys. Lett.* 381, 777 (2003); Chen, Cheng & Ye, *JOSA B* 21, 820 (2004); Chen, de Jong & Ye, *JOSA B* 22, 951 (2005) | Global description of B-state hyperfine interactions and their link to potentials and wavefunctions |
| NMIJ/AIST and Yokohama (F.-L. Hong) | Kobayashi et al., *JOSA B* 33, 725 (2016), 578 nm; Ikeda et al., *JOSA B* 39, 2264 (2022), 514 nm; Yoshiki et al., *EPJ D* 77, 140 (2023), rotation dependence of v′ = 44 hyperfine constants; Nishiyama et al., *JOSA B* 41, 2290 (2024), 520.2 nm | Hyperfine constants from fits to the four-term Hamiltonian |
| Taiwan (NTHU; J.-T. Shy, L.-B. Wang) | Huang et al., *Appl. Opt.* 52, 1448 (2013) and 57, 2102 (2018); Fan et al., *PRA* 89, 032513 (2014) | Use IodineSpec 5 for predictions; comb recalibration of iodine reference lines |
| China | Cheng et al., *JOSA B* 36, 1816 (2019), 534 nm | Absolute frequency measurement |

No public code from any of these groups was found on GitHub, GitLab, Zenodo or PyPI [I]. No Czech or Korean modeling code was found; the Czech work found, e.g. Balling & Křen, *EPJ D* 48, 3 (2007), is measurement-focused.

**Fit [I]:** these papers are sources of data and hyperfine constants (the data compilation), not reusable software.

### 4.4 Iodine cells in astronomy (precision radial velocities)

**The standard method.** The I₂-cell technique (Butler et al., PASP 108, 500, 1996) forward-models the stellar spectrum through a laboratory FTS scan of the specific cell. It does not compute the I₂ spectrum from molecular physics [I].

**Open pipeline.** **pyodine** (Heeren et al., A&A 674, A164, 2023) [S]:
* MIT license; GitLab `Heeren/pyodine`, GitHub mirror `pepeheeren/pyodine`, docs at pyodine.readthedocs.io.
* Adapted to SONG and Lick/Hamilton data.
* Forks for APO (`jaklusmeyer/APO_pyodine`, 2026) and `holtzmanjon/pyodine` (2025).

**Cell characterisation and calibration papers:**
* Wang et al., PASP 132, 014503 (2019): FTS calibration of cells.
* Debus et al., JATIS 9, 045003 (2023).
* Reiners et al., A&A 690, A210 (2024): IodineSpec plus comb-FTS (§2.3).
* Schmidt et al., MNRAS 539, 3301 (2025): ESPRESSO, "full forward modeling approach of the I₂ spectrum". The abstract does not say which model.

**Fit [I].** This community is a natural user group. It needs absolute, physics-based line lists with honest uncertainties across 500–630 nm, as Reiners et al. show. It would not supply code.

### 4.5 Atmospheric chemistry and other databases

* **MPI-Mainz UV/VIS Spectral Atlas** (Keller-Rudek et al., ESSD 5, 365, 2013): has an I₂ entry with cross-section files and plots from 170 to 800 nm [S: live page]. These are data compilations, e.g. Saiz-Lopez et al., ACP 4, 1443 (2004), not tools.
* **HITRAN:** the molecule list has no I₂ [S].
* **ExoMol:** no I₂ dataset [S].
* **BIPM:** the CIPM *mise en pratique* "standard frequencies" page lists recommended values for eight I₂-stabilized lines:
  * 582, 564, 563, 552, 520, 490, 474 and 468 THz.
  * The most recent I₂ entry was updated in 2015 (564 THz).
  * The list is available as PDF and via an API; no software is offered [S].
* **NIST:** the Chemistry WebBook has I₂ diatomic constants; the NBS atlas is noted above. No NIST I₂ line-calculation tool was found [I].

### 4.6 A likely parallel effort: ExoMol/MARVEL I₂ bibliography

The public GitHub repository `ExoMol/bib` contains two I₂ files [S: GitHub API]:

| File | Size | History |
|---|---|---|
| `MARVEL/I2_MARVEL.bib` | 330 entries | Uploaded 2024-09-03 by `jonnytennyson` |
| `exomol/I2.bib` | Small | Updated 2024-05 to 2026-01 |

`I2_MARVEL.bib` lists the classic and comb-era I₂ frequency measurements, e.g. 501.7 nm, 532 nm grids, 535, 543, 548, 578 and 647 nm, alongside MARVEL methodology papers. No MARVEL or ExoMol I₂ paper was found in Crossref or by web search.

**Inference [I]:** a MARVEL analysis of I₂ transitions (and possibly an ExoMol line list) is being compiled at UCL/ELTE. It would overlap with the data compilation and fitting here; their bibliography is a useful checklist of I₂ data.

### 4.7 Teaching and unrelated items (for completeness)

* Undergraduate B–X lab notebooks: `act-cms/vibronic-spectrum-of-iodine` (a template) and `Schmiddy2000/MolecularIodine`.
* `Lazeau/iodine_hfs` concerns atomic I II hyperfine structure, not I₂.
* None is relevant.

---

## 5. General-purpose diatomic codes

### 5.1 Le Roy suite: LEVEL16, RKR1-16, dParFit16, dPotFit16, betaFIT16, BCONT

**Papers:** JQSRT 186 (2017) — LEVEL pp. 167–178, RKR1 pp. 158–166, dPotFit pp. 179–196, dParFit pp. 197–209, betaFIT pp. 210 ff.

**Distribution [S].** The University of Waterloo ACE site hosts "a reproduction of the late Robert J. LeRoy's computer programs and documentation". The original `scienide2.uwaterloo.ca/~rleroy` host did not respond during this survey. The site states:
* LEVEL16 was last updated on 11 August 2016.
* Source code and manuals are also JQSRT supplementary files.
* "All references in the code/documentation for support are currently deprecated."

**Terms [S]:** "distributed free of charge … may not be sold or have any other commercial use made of it without the express written permission of the author." Because of this non-commercial clause it is **not OSI open source** [I].

**Relevance [I]:**
* LEVEL is the community reference for eigenvalues, FCFs and centrifugal-distortion constants of arbitrary potentials. It is ideal for cross-validating our radial solver on the Knöckel 2004 potentials.
* dPotFit's direct-potential-fit strategy and MLR potential forms are candidate algorithms for the global fit.
* We should reimplement rather than vendor the code, because of the license and the lack of maintenance.

### 5.2 PGOPHER

**Paper:** JQSRT 186, 221 (2017).

**License and status [S]:**
* The Bristol data repository record for version 10.1.180 (made available 28 Nov 2018) states: "The program is open source, and covered by the LGPL." Documentation falls under a non-commercial government license.
* The website states (23 Mar 2022): "Regrettably, further updates to PGOPHER will no longer be available."
* The website notes that large hyperfine structure "such as found with iodine" requires `AllowComplex`. That note refers to microwave spectra.

**Relevance [I]:** an effective-Hamiltonian tool, not potential-based. Useful to cross-check hyperfine patterns and line strengths for individual bands. Frozen and single-author.

### 5.3 Duo (and ExoCross)

**Paper:** CPC 202, 262 (2016). It lists "Licensing provisions: Standard CPC licence" and Fortran 2003 [S].

**Repository [S].** The active repository is `Trovemaster/Duo`:
* Last push 2026-09-08.
* No LICENSE file.
* Includes `F1_hyperfine.f90` and `refinement.f90`.

**Hyperfine capability [S]:**
* Qu, Yurchenko & Tennyson, JCTC 18, 1808 (2022): variational hyperfine-resolved rovibronic spectra.
* Yin, Yurchenko & Tennyson, JCP 165, 092501 (2026), published 2026-09-04: extends Duo "to calculate hyperfine-resolved spectra of homonuclear diatomic molecules". It includes five magnetic-dipole terms and one electric-quadrupole term, with explicit nuclear-exchange symmetry. It was validated on H₂, D₂, ¹⁴N₂ and ¹⁴N₂⁺ against experiment and PGOPHER.

**Relevance [I]:**
* Duo is now the only open general code found that could, in principle, produce a hyperfine-resolved, potential-based ¹²⁷I₂ B–X line list with intensities.
* It has not been applied to I₂, as far as could be found.
* Its homonuclear-hyperfine implementation is brand new.
* Its license is ambiguous: the CPC licence is non-OSI, and GitHub has no LICENSE file.
* Best use: an independent cross-validation engine.

**ExoCross** (A&A 614, A131, 2018; `Trovemaster/exocross`, MIT, active) and **PyExoCross** (RASTI 3, 257, 2024; GPL-3.0) turn line lists into temperature-dependent cross sections. ExoCross's MIT license allows reuse [S/I].

### 5.4 PyDiatomic (S. Gibson)

[S] `stggh/PyDiatomic`, GPL-3.0, Python. Tags: v0.1-alpha (2016) and v0.3 (2025-02-03). From its README:
> "solves the time-independent coupled-channel Schroedinger equation using the Johnson renormalized Numerov method … directed to the computation of photodissociation cross sections for diatomic molecules"

**Relevance [I]:** a good algorithmic reference and a Python-native cross-check, e.g. for B-state predissociation and coupled channels. GPL-3.0 means we should not copy code into a permissively licensed package.

### 5.5 Pashov IPA and coupled-channel fitting codes

[S] Pashov, Jastrzębski & Kowalczyk, CPC 128, 622 (2000): "Construction of potential curves for diatomic molecular states by the IPA method". It is in the CPC Program Library (ADLV_v1_0, programs "IPA, SCHROED") under the CPC licence.

Knöckel 2004 refers to Pashov's spline-interpolated pointwise potentials as an alternative representation [S]. Pashov's later coupled-channel fitting codes (used with the Hannover group on alkali dimers) were not found publicly [I].

**Relevance [I]:** algorithms (IPA and regularised IPA, pointwise potentials) are worth borrowing; the code is not.

### 5.6 MARVEL

[S] JMS 245, 115 (2007) and JQSRT 113, 929 (2012). Available as MARVEL Online; the `I2_MARVEL.bib` file also cites "MARVEL 4.1" and "MARVEL ONLINE 2.1". The online app did not respond from this machine on 2026-09-14. License unknown.

**Relevance [I]:** the spectroscopic-network approach checks the consistency of heterogeneous I₂ data and assigns empirical level uncertainties. That is directly useful for the global fit. See also §4.6.

### 5.7 Other Python and Julia packages

| Package | Details [S] | Use [I] |
|---|---|---|
| `diatomic-py` | Durham; CPC 282, 108512 (2023); BSD-3; PyPI 2.1.0 (2026-05). Hyperfine plus field structure of ¹Σ states. | Permissively licensed reference for coupled-basis hyperfine matrix elements. |
| `numerov` | R. Bast; MPL-2.0; 2017 | Trivial. |
| `molspecutils` | MIT, 2022 | HITRAN-oriented. |
| `diatomics` | Wrapper for the FHI Diatomic Molecular Spectroscopy Database | Constants database only. |
| `jj2223/diatomic_py` | GitHub, 2026, no license. Hyperfine and Zeeman levels for homo- and heteronuclear diatomics. | Unevaluated. |
| Julia General registry | Only a generic `Numerov` package | Nothing relevant. |

**Coverage gap [I].** Nothing in the Python, Julia or Rust ecosystems provides a fitted, potential-based B–X model with hyperfine structure.

---

## 6. Conclusion

### 7.1 Does a suitable open-source tool exist?

**No.**

**What the only complete tool lacks.** IodineSpec is complete but closed. It has no public source, no license, no per-line uncertainties and informal distribution. Its global model dates from 2004, with partial updates in 2008 and 2010. An independent 2024 comb-FTS test found ~2 m s⁻¹ (~3 MHz) systematic structure in its energy scale [S: Reiners et al. 2024]. Its hyperfine interpolation formulae miss some 514 nm components by up to ~0.8 MHz [S: Yoshiki et al. 2023].

**Why the open alternatives fall short:**
* The Forkey model is narrow (532 nm) and outdated.
* The atlases are data.
* The astronomy pipelines use measured cell spectra.
* The general codes (LEVEL, Duo, PGOPHER, PyDiatomic) can compute pieces, but none ships an I₂ model, fitted parameters or a data compilation.

### 7.2 Go / no-go

**GO: build our own package.** The niche is real: users in metrology, laser cooling, astronomy and lidar currently rely on a closed program obtained by e-mail. We should also coordinate early with:
1. **E. Tiemann (Hannover)** — permission to use and compare against IodineSpec 5 line lists, the Table 4 erratum, and any unpublished data.
2. **The ExoMol/MARVEL group** — their I₂ bibliography suggests an I₂ MARVEL compilation.

### 7.3 What to reuse, benchmark against, or borrow

**Reuse as code or dependencies (licenses compatible with a permissive package) [I]:**
* The SciPy stack.
* ExoCross (MIT), or just the ExoMol file format, for cross-section and line-list interoperability.
* `diatomic-py` (BSD-3), as a reference for hyperfine matrix elements.
* `numerov` (MPL-2.0) is optional; our own solver will replace it.

**Benchmark against:**

| Benchmark | How to use it |
|---|---|
| **Knöckel 2004 Table 4** (with the IQO erratum) | Milestone 1: reproduce the published potentials and predictions. |
| **IodineSpec 5 line lists** | Obtain from Tiemann if possible. Reiners et al. extracted such lists. |
| **Liao 2010 NIR model** and Salumbides 2006/2008 | Isotopologue and BOC checks. |
| **LEVEL16** | Eigenvalues and FCFs for identical potentials. |
| **Duo** (2026 homonuclear hyperfine) | Independent variational hyperfine-resolved check. |
| **PGOPHER** | Band-by-band hyperfine patterns and intensities. |
| **Validation data** | Kato atlas (CD-ROM), Salami–Ross ASCII atlas, Gerstenkorn–Luc, Nölleke 915–985 nm, the Reiners 2024 comb-referenced FTS spectrum, BIPM recommended I₂ frequencies. |

**Borrow algorithms (reimplement; do not copy code):**
* From Knöckel 2004: the X-representation potential with inner and outer extensions, and R-dependent adiabatic and nonadiabatic BOCs.
* From Bodermann 2002: the hyperfine interpolation formulae.
* The four-term hyperfine Hamiltonian, as used across the metrology literature.
* From Le Roy: MLR potentials and direct-potential-fit strategy (dPotFit); LEVEL's quasibound-level handling.
* From Pashov: IPA and pointwise-spline potentials.
* From PyDiatomic: Johnson renormalized Numerov (coupled channels, predissociation).
* From MARVEL: network-based consistency checking and empirical level uncertainties.
* From Duo: variational treatment of hyperfine structure with exchange symmetry.

### 7.4 Caveats

1. The IodineSpec price and license terms are **unknown**. The "discontinued by TOPTICA" conclusion is an inference from archives and a site search, not a statement from TOPTICA.
2. It is unknown which post-2004 model updates (2008 BOC, 2010 NIR) are inside IodineSpec 5.
3. The Salami–Ross range and resolution and the format of the Kato CD-ROM are unverified.
4. The Duo license is ambiguous: "Standard CPC licence" in the paper, no LICENSE file on GitHub. Le Roy's codes carry a non-commercial clause. PyDiatomic and PyExoCross are GPL-3.0. PGOPHER is LGPL. None of this blocks benchmarking, but it does restrict copying code.
5. The ExoMol/MARVEL I₂ effort is inferred from bibliography files only.

---

## References

All DOIs below were checked to resolve at doi.org on 2026-09-14. Metadata comes from Crossref or OpenAlex.

### Hannover model and IodineSpec
1. B. Bodermann, H. Knöckel, E. Tiemann, "Widely usable interpolation formulae for hyperfine splittings in the ¹²⁷I₂ spectrum", *Eur. Phys. J. D* **19**, 31–44 (2002). https://doi.org/10.1140/epjd/e20020052
2. H. Knöckel, B. Bodermann, E. Tiemann, "High precision description of the rovibronic structure of the I₂ B–X spectrum", *Eur. Phys. J. D* **28**, 199–209 (2004). https://doi.org/10.1140/epjd/e2003-00313-4
3. E. J. Salumbides, K. S. E. Eikema, W. Ubachs, U. Hollenstein, H. Knöckel, E. Tiemann, "The hyperfine structure of ¹²⁹I₂ and ¹²⁷I¹²⁹I in the B³Π0u⁺–X¹Σg⁺ band system", *Mol. Phys.* **104**, 2641–2652 (2006). https://doi.org/10.1080/00268970600747696
4. E. J. Salumbides et al. (same authors), "Improved potentials and Born–Oppenheimer corrections by new measurements of transitions of ¹²⁹I₂ and ¹²⁷I¹²⁹I in the B–X band system", *Eur. Phys. J. D* **47**, 171–179 (2008). https://doi.org/10.1140/epjd/e2008-00045-y
5. C.-C. Liao, K.-Y. Wu, Y.-H. Lien, H. Knöckel, H.-C. Chui, E. Tiemann, J.-T. Shy, "Precise frequency measurements of ¹²⁷I₂ lines in the wavelength region 750–780 nm", *JOSA B* **27**, 1208–1214 (2010). https://doi.org/10.1364/JOSAB.27.001208
6. Institut für Quantenoptik, LU Hannover, "Secondary Frequency Standard" (IodineSpec description, Table 4 erratum, contact). https://www.iqo.uni-hannover.de/de/arbeitsgruppen/molecular-quantum-gases/research-projects/secondary-frequency-standard (accessed 2026-09-14)
7. TOPTICA Photonics, "Iodine Spectrum Calculating Software: IodineSpec 4", archived product pages:
   * http://www.toptica.com/page/Iodine_Spectrum_Calculating_Software_IodineSpec4.php (capture showing "Latest update: Nov 17, 2008")
   * http://www.toptica.com/products/iodinespec/iodinespec_1.htm (captures 2002–2004)
   * Brochure `BR-32014A-Iodine-Calculation-Software.pdf` (captures 2005–2009)

   All via the Internet Archive Wayback Machine, https://web.archive.org/

### Independent tests and users of IodineSpec
8. A. Reiners, M. Debus, S. Schäfer, E. Tiemann, M. Zechmeister, "Accurate calibration spectra for precision radial velocities — Iodine absorption referenced by a laser frequency comb", *A&A* **690**, A210 (2024). https://doi.org/10.1051/0004-6361/202451389 ; arXiv:2409.02631
9. Y.-C. Huang, H.-C. Chen, S.-E. Chen, J.-T. Shy, L.-B. Wang, "Precise frequency measurements of iodine hyperfine transitions at 671 nm", *Appl. Opt.* **52**, 1448 (2013). https://doi.org/10.1364/AO.52.001448
10. Y.-C. Huang, Y.-C. Guan, T.-H. Suen, J.-T. Shy, L.-B. Wang, "Absolute frequency measurement of molecular iodine hyperfine transitions at 647 nm", *Appl. Opt.* **57**, 2102 (2018). https://doi.org/10.1364/AO.57.002102 ; arXiv:1710.09533
11. I. Fan, C.-Y. Chang, L.-B. Wang, S. L. Cornish, J.-T. Shy, Y. Liu, "Refined determination of the muonium-deuterium 1S–2S isotope shift through improved frequency calibration of iodine lines", *Phys. Rev. A* **89**, 032513 (2014). https://doi.org/10.1103/PhysRevA.89.032513
12. M. Debus, S. Schäfer, A. Reiners, "Toward 10 cm s⁻¹ radial velocity accuracy on the Sun using a Fourier transform spectrometer", *J. Astron. Telesc. Instrum. Syst.* **9**, 045003 (2023). https://doi.org/10.1117/1.JATIS.9.4.045003
13. T. M. Schmidt, A. Reiners, M. T. Murphy, G. Lo Curto, C. J. A. P. Martins, P. Huke, "Validation of the ESPRESSO wavelength calibration using iodine absorption cell spectra", *MNRAS* **539**, 3301–3318 (2025). https://doi.org/10.1093/mnras/staf588

### Atlases and I₂ data resources
14. H. Katô et al., *Doppler-Free High Resolution Spectral Atlas of Iodine Molecule 15000 to 19000 cm⁻¹* (Japan Society for the Promotion of Science, 2000), ISBN 4-89114-000-3. Atlas page: http://web1.kcn.jp/kansha-kansha/AtlasofI2.html ; OSU Knowledge Bank abstract: http://hdl.handle.net/1811/19890
15. H. Salami, A. J. Ross, "A molecular iodine atlas in ascii format", *J. Mol. Spectrosc.* **233**, 157–159 (2005). https://doi.org/10.1016/j.jms.2005.06.002
16. S. Gerstenkorn, P. Luc, "Absolute iodine (I₂) standards measured by means of Fourier transform spectroscopy", *Rev. Phys. Appl.* **14**, 791–794 (1979). https://doi.org/10.1051/rphysap:01979001408079100 (The atlases themselves: Laboratoire Aimé Cotton, CNRS, 1977–1982; no DOI.)
17. J. D. Simmons, J. T. Hougen, "Atlas of the I₂ spectrum from 19 000 to 18 000 cm⁻¹", *J. Res. Natl. Bur. Stand.* **81A**, 25–80 (1977). https://doi.org/10.6028/jres.081A.006
18. C. Nölleke, C. Raab, R. Neuhaus, S. Falke, "Absolute frequency atlas from 915 nm to 985 nm based on laser absorption spectroscopy of iodine", *J. Mol. Spectrosc.* **346**, 19–22 (2018). https://doi.org/10.1016/j.jms.2017.12.013 ; arXiv:1801.01342
19. D. Rodríguez Fernández et al., "High resolution laser spectroscopy of iodine molecule in the 14400–14600 cm⁻¹ range", *J. Mol. Spectrosc.* (2023). https://doi.org/10.1016/j.jms.2023.111789
20. A. J. Ross, P. Cardon, "Te₂ absorption spectrum from 19000 to 24000 cm⁻¹", *J. Mol. Spectrosc.* **384**, 111589 (2022). https://doi.org/10.1016/j.jms.2022.111589
21. BIPM, *Mise en pratique* — standard frequencies. https://www.bipm.org/en/publications/mises-en-pratique/standard-frequencies (accessed 2026-09-14)
22. H. Keller-Rudek, G. K. Moortgat, R. Sander, R. Sörensen, "The MPI-Mainz UV/VIS Spectral Atlas of Gaseous Molecules of Atmospheric Interest", *Earth Syst. Sci. Data* **5**, 365–373 (2013). https://doi.org/10.5194/essd-5-365-2013 ; I₂ entry: https://uv-vis-spectral-atlas-mainz.org/uvvis/cross_sections/Halogens+mixed%20halogens/I2.spc
23. A. Saiz-Lopez, R. W. Saunders, D. M. Joseph, S. H. Ashworth, J. M. C. Plane, "Absolute absorption cross-section and photolysis rate of I₂", *Atmos. Chem. Phys.* **4**, 1443–1450 (2004). https://doi.org/10.5194/acp-4-1443-2004
24. HITRAN molecule list: https://hitran.org/docs/molec-meta/ ; ExoMol master file: https://www.exomol.com/db/exomol.all ; NIST Chemistry WebBook, I₂: https://webbook.nist.gov/cgi/cbook.cgi?ID=C7553562&Mask=1000 (all accessed 2026-09-14)

### Metrology hyperfine work (in-house analyses)
25. L. Chen, J. Ye, "Extensive, high-resolution measurement of hyperfine interactions: precise investigations of molecular potentials and wave functions", *Chem. Phys. Lett.* **381**, 777–783 (2003). https://doi.org/10.1016/j.cplett.2003.10.052
26. L. Chen, W.-Y. Cheng, J. Ye, "Hyperfine interactions and perturbation effects in the B0u⁺(³Πu) state of ¹²⁷I₂", *JOSA B* **21**, 820 (2004). https://doi.org/10.1364/JOSAB.21.000820
27. L. Chen, W. A. de Jong, J. Ye, "Characterization of the molecular iodine electronic wave functions and potential energy curves through hyperfine interactions in the B0u⁺(³Πu) state", *JOSA B* **22**, 951 (2005). https://doi.org/10.1364/JOSAB.22.000951
28. T. Kobayashi et al., "Absolute frequency measurements and hyperfine structures of the molecular iodine transitions at 578 nm", *JOSA B* **33**, 725 (2016). https://doi.org/10.1364/JOSAB.33.000725
29. K. Ikeda, T. Kobayashi, M. Yoshiki, D. Akamatsu, F.-L. Hong, "Hyperfine structure and absolute frequency of ¹²⁷I₂ transitions at 514 nm for wavelength standards at 1542 nm", *JOSA B* **39**, 2264 (2022). https://doi.org/10.1364/JOSAB.465499
30. M. Yoshiki, S. Matsunaga, K. Ikeda, D. Akamatsu, F.-L. Hong, "Rotation dependence of v′ = 44 excited-state hyperfine constants …", *Eur. Phys. J. D* **77**, 140 (2023). https://doi.org/10.1140/epjd/s10053-023-00712-7 ; accepted manuscript (open access): https://ynu.repo.nii.ac.jp/records/2001411
31. A. Nishiyama, S. Okubo, T. Kobayashi, A. Kawasaki, H. Inaba, "Measurement of transition frequencies and hyperfine constants of molecular iodine at 520.2 nm", *JOSA B* **41**, 2290 (2024). https://doi.org/10.1364/JOSAB.531115
32. Cheng et al., "Absolute frequency measurement of molecular iodine hyperfine transition at 534 nm", *JOSA B* **36**, 1816 (2019). https://doi.org/10.1364/JOSAB.36.001816

### Other I₂-specific code
33. J. N. Forkey, W. R. Lempert, R. B. Miles, "Corrected and calibrated I₂ absorption model at frequency-doubled Nd:YAG laser wavelengths", *Appl. Opt.* **36**, 6729–6738 (1997). https://doi.org/10.1364/AO.36.006729
    * Code copy: https://github.com/Lightning469/Forkey-Iodine-Cell-Transmission-Code
    * MATLAB port: https://github.com/PriybratS/Iodine_absorption_Forky_Matlab
    * Similar: https://github.com/eckeratvt/Iodine-Cell-Transmission-Spectra
34. R. P. Butler, G. W. Marcy, E. Williams, C. McCarthy, P. Dosanjh, S. S. Vogt, "Attaining Doppler precision of 3 m s⁻¹", *PASP* **108**, 500 (1996). https://doi.org/10.1086/133755
35. P. Heeren, R. Tronsgaard, F. Grundahl, S. Reffert, A. Quirrenbach, E. Pallé, "Pyodine: an open, flexible reduction software for iodine-calibrated precise radial velocities", *A&A* **674**, A164 (2023). https://doi.org/10.1051/0004-6361/202244441
    * Code: https://gitlab.com/Heeren/pyodine , https://github.com/pepeheeren/pyodine
    * Docs: https://pyodine.readthedocs.io
36. S. X. Wang, J. T. Wright, P. MacQueen, W. D. Cochran, D. R. Doss, C. A. Gibson, J. R. Schmitt, "Calibrating iodine cells for precise radial velocities", *PASP* **132**, 014503 (2019). https://doi.org/10.1088/1538-3873/ab5021
37. K. Döringshoff et al., "Iodine frequency reference on a sounding rocket", *Phys. Rev. Applied* **11**, 054068 (2019). https://doi.org/10.1103/PhysRevApplied.11.054068 ; also *IFCS 2018*, https://doi.org/10.1109/fcs.2018.8597507
    * Control software ("pyodine"): https://github.com/trimitri/jokarus

### General-purpose diatomic codes
38. R. J. Le Roy, "LEVEL: A computer program for solving the radial Schrödinger equation for bound and quasibound levels", *JQSRT* **186**, 167–178 (2017). https://doi.org/10.1016/j.jqsrt.2016.05.028
39. R. J. Le Roy, "RKR1: A computer program implementing the first-order RKR method …", *JQSRT* **186**, 158–166 (2017). https://doi.org/10.1016/j.jqsrt.2016.03.030
40. R. J. Le Roy, "dPotFit: A computer program to fit diatomic molecule spectral data to potential energy functions", *JQSRT* **186**, 179–196 (2017). https://doi.org/10.1016/j.jqsrt.2016.06.002
41. R. J. Le Roy, "dParFit: …", *JQSRT* **186**, 197–209 (2017). https://doi.org/10.1016/j.jqsrt.2016.04.004
42. R. J. Le Roy, A. Pashov, "betaFIT: A computer program to fit pointwise potentials to selected analytic functions", *JQSRT* **186**, 210–220 (2017). https://doi.org/10.1016/j.jqsrt.2016.03.036
43. LeRoy Programs (University of Waterloo reproduction). https://uwaterloo.ca/atmospheric-chemistry-experiment/leroy-programs (accessed 2026-09-14)
44. C. M. Western, "PGOPHER: A program for simulating rotational, vibrational and electronic spectra", *JQSRT* **186**, 221–242 (2017). https://doi.org/10.1016/j.jqsrt.2016.04.010
    * PGOPHER v10.1 dataset: https://doi.org/10.5523/bris.3mqfb4glgkr8a2rev7f73t300c
    * Website: https://pgopher.chemistry.bristol.ac.uk/
45. S. N. Yurchenko, L. Lodi, J. Tennyson, A. V. Stolyarov, "Duo: A general program for calculating spectra of diatomic molecules", *Comput. Phys. Commun.* **202**, 262–275 (2016). https://doi.org/10.1016/j.cpc.2015.12.021
    * Code: https://github.com/Trovemaster/Duo
46. Q. Qu, S. N. Yurchenko, J. Tennyson, "A method for the variational calculation of hyperfine-resolved rovibronic spectra of diatomic molecules", *J. Chem. Theory Comput.* **18**, 1808–1820 (2022). https://doi.org/10.1021/acs.jctc.1c01244
47. Yin, S. N. Yurchenko, J. Tennyson, "Hyperfine-resolved variational nuclear motion spectra of homonuclear diatomic molecules", *J. Chem. Phys.* **165**, 092501 (2026). https://doi.org/10.1063/5.0346970
48. S. N. Yurchenko, A. F. Al-Refaie, J. Tennyson, "ExoCross: a general program for generating spectra from molecular line lists", *A&A* **614**, A131 (2018). https://doi.org/10.1051/0004-6361/201732531
    * Code: https://github.com/Trovemaster/exocross
49. Zhang, Tennyson, Yurchenko, "PyExoCross: a Python program for generating spectra and cross-sections from molecular line lists", *RAS Tech. Instrum.* **3**, 257–287 (2024). https://doi.org/10.1093/rasti/rzae016
    * Code: https://github.com/ExoMol/PyExoCross
50. S. Gibson, PyDiatomic. https://github.com/stggh/PyDiatomic
51. A. Pashov, W. Jastrzębski, P. Kowalczyk, "Construction of potential curves for diatomic molecular states by the IPA method", *Comput. Phys. Commun.* **128**, 622–634 (2000). https://doi.org/10.1016/S0010-4655(00)00010-2
    * CPC Library entry: https://elsevier.digitalcommonsdata.com/datasets/7v33vxjxjv/1
52. T. Furtenbacher, A. G. Császár, J. Tennyson, "MARVEL: measured active rotational–vibrational energy levels", *J. Mol. Spectrosc.* **245**, 115–125 (2007). https://doi.org/10.1016/j.jms.2007.07.005
53. T. Furtenbacher, A. G. Császár, "MARVEL: … II. Algorithmic improvements", *JQSRT* **113**, 929–935 (2012). https://doi.org/10.1016/j.jqsrt.2012.01.005
    * MARVEL Online: http://kkrk.chem.elte.hu/marvelonline
54. ExoMol bibliography repository: `MARVEL/I2_MARVEL.bib` and `exomol/I2.bib`. https://github.com/ExoMol/bib (accessed 2026-09-14)
55. H. M. Pickett, "The fitting and prediction of vibration-rotation spectra with spin interactions", *J. Mol. Spectrosc.* **148**, 371–377 (1991). https://doi.org/10.1016/0022-2852(91)90393-O
    * Repositories: https://github.com/laserkelvin/Pickett , https://github.com/Ltotheois/Pyckett , https://github.com/jtbr/spfit-spcat
56. J. A. Blackmore, P. D. Gregory, J. M. Hutson, S. L. Cornish, "Diatomic-py: A Python module for calculating the rotational and hyperfine structure of ¹Σ molecules", *Comput. Phys. Commun.* **282**, 108512 (2023). https://doi.org/10.1016/j.cpc.2022.108512
    * Code: https://github.com/durham-qlm/diatomic-py
57. R. Bast, `numerov`. https://github.com/bast/numerov
