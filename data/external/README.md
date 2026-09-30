# External datasets (third-party; do not redistribute without checking terms)

These are raw data files downloaded for the global fit (stage 5). Catalog entries are in
`data/catalog/broadband.yaml` and `data/catalog/precision.yaml`, and SHA-256 checksums are in
`SHA256SUMS`.

**These files belong to their authors or publishers.** Keep them out of any public repository or
package until the terms have been checked. Add `data/external/*/` to `.gitignore` when the repo is
initialized. The package should download them, or ask for them, not ship them.

## `salami_ross_2005/salami_ross_2005_mmc1.txt`

- **Source:** H. Salami & A. J. Ross, "A molecular iodine atlas in ascii format", *J. Mol. Spectrosc.* **233**, 157–159 (2005), doi:10.1016/j.jms.2005.06.002. Elsevier supplementary file mmc1.
- **URL:** https://ars.els-cdn.com/content/image/1-s2.0-S002228520500130X-mmc1.txt
- **Downloaded:** 2026-09-15.
- **Contents:** Doppler-limited FTS transmission spectrum of the B–X system, in % transmission against wavenumber (cm⁻¹).
  - 14250.100–20050.100 cm⁻¹ with a 0.005 cm⁻¹ step: 1,160,001 points after a 5-line header.
  - The file doesn't record the cell conditions (temperature, path length, pressure) or the instrument resolution. Get these from the paper before fitting line shapes.
- **Terms:** Elsevier supplementary material, license not stated.

## `nolleke_2018/`

- **Source:** C. Nölleke et al., *J. Mol. Spectrosc.* **346**, 19 (2018), doi:10.1016/j.jms.2017.12.013. Data on Mendeley Data, "Laser absorption spectroscopy of iodine between 915 and 985 nm", doi:10.17632/cvr6vv2zps.1, published 2017-12-15.
- **Downloaded:** 2026-09-15. The `.xls` and `.xlsx` duplicates were not downloaded.
- **`iodine_atlas.csv`:** line list, with about 10,000 lines stated at about 50 MHz uncertainty.
- **`scan_data.csv`:** the underlying laser absorption scans.
- **Terms:** **CC BY 4.0** (Mendeley `data_licence`). Redistribution and modification are allowed with attribution and a note of any changes.

## `saiz_lopez_2004/`

- **Source:** A. Saiz-Lopez, R. W. Saunders, D. M. Joseph, S. H. Ashworth & J. M. C. Plane, "Absolute absorption cross-section and photolysis rate of I₂", *Atmos. Chem. Phys.* **4**, 1443–1450 (2004), doi:10.5194/acp-4-1443-2004. ACP supplement.
- **URL:** https://acp.copernicus.org/articles/4/1443/2004/acp-4-1443-2004-supplement.zip
- **Downloaded:** 2026-09-15.
- **`acp-4-1443-2004-supplement.zip`:** the supplement as downloaded. It holds a single file, `Iodine_xsection/Iodine_xsection.xls`.
- **`Iodine_xsection.xls`:** extracted from the zip. Excel 97 format, one sheet named `I2_295 K`, two columns: `wavelength / nm` and `cross-section / cm^2molecule^-1`.
  - Row 2 holds a stray value, 2003.0, in the cross-section column with an empty wavelength. Skip it.
  - 1474 data rows from 182 to 750 nm. The step is 1 nm for 182–500 and 630–750 nm; 500–630 nm is mixed 0.1/0.2 nm, with 0.2 nm steps becoming more common toward the red.
  - Cross section at 295 K and 760 Torr of air, from an FTS at 4 cm⁻¹. The paper doesn't say whether wavelengths are in air or vacuum. For conditions, and a cross-correlation that puts them on a vacuum scale, see `docs/research/cross-section-data.md`.
  - The MPI-Mainz copy (`mpi_mainz_i2/I2_Saiz-Lopez(2004)_295K_182-750nm.txt`) has the same points, agreeing to ≤5 × 10⁻⁶ relative.
- **Terms:** the article is under **CC BY-NC-SA 2.5**, according to the ACP landing page. The supplement's licence isn't stated separately.

## `spietz_2006/`

- **Source:** P. Spietz, J. C. Gómez Martín & J. P. Burrows, "Effects of column density on I₂ spectroscopy and a determination of I₂ absorption cross section at 500 nm", *Atmos. Chem. Phys.* **6**, 2177–2191 (2006), doi:10.5194/acp-6-2177-2006. ACP supplement.
- **URL:** https://acp.copernicus.org/articles/6/2177/2006/acp-6-2177-2006-supplement.zip
- **Downloaded:** 2026-09-15.
- **`acp-6-2177-2006-supplement.zip`:** the supplement as downloaded. It holds the two TXT files below.
- **`I2DOASref_1200.TXT`:** the 0.25 nm FWHM spectrum, 543.217–577.952 nm, 1024 points.
  - Recording: 1200 grooves/mm grating, 170 µm slit, 0.035 nm/pixel.
  - Column and cell: I₂ column (6.86 ± 0.20) × 10¹⁵ cm⁻², 26.4 cm path, 1000 mbar of N₂.
  - Wavelengths are vacuum (±0.04 nm).
  - Uncertainty: the larger of 4 % of σ and 7 × 10⁻²⁰ cm².
- **`I2DOASref_0300.TXT`:** the 0.59 nm FWHM spectrum, 428.463–588.287 nm, 1024 points.
  - Recording: 300 grooves/mm grating, 50 µm slit, 0.154 nm/pixel.
  - Column and cell: I₂ column (1.42 ± 0.038) × 10¹⁶ cm⁻², 26.4 cm path, 1000 mbar of N₂.
  - Wavelengths are vacuum (±0.11 nm) as stated. **But the scale appears to run 0.09–0.17 nm long;** see the research doc.
  - Uncertainty: the larger of 2.7 % of σ and 7 × 10⁻²⁰ cm².
- **Format of both files:** ISO-8859-1 text with CRLF line endings, and a ~30-line header giving the conditions above. Data are two whitespace-separated columns, `lambdaVac` (nm) and `sigma` (cm² molecule⁻¹), scaled to σ(500 nm) = 2.191 × 10⁻¹⁸ cm².
- **Other copies:** these are the same data as the MPI-Mainz copies. They also match the IUP Bremen copies (https://www.iup.uni-bremen.de/gruppen/molspec/I2/RefSpec_0.25nmFWHM_I2_IUPBremen.TXT and `RefSpec_0.59nmFWHM_I2_IUPBremen.TXT`), which have LF line endings and an extra citation note. The Bremen copies were not stored.
- **Terms:** the article is under **CC BY-NC-SA 2.5**, according to the ACP landing page. The supplement's licence isn't stated separately.

## `mpi_mainz_i2/`

- **Source:** H. Keller-Rudek, G. K. Moortgat, R. Sander & R. Sörensen, "The MPI-Mainz UV/VIS Spectral Atlas of Gaseous Molecules of Atmospheric Interest", *Earth Syst. Sci. Data* **5**, 365–373 (2013), doi:10.5194/essd-5-365-2013. I₂ page: https://uv-vis-spectral-atlas-mainz.org/uvvis/cross_sections/Halogens+mixed%20halogens/I2.spc
- **URL pattern:** `https://uv-vis-spectral-atlas-mainz.org/uvvis_data/cross_sections/Halogens+mixed%20halogens/<file>`, with the parentheses and commas URL-encoded.
- **Downloaded:** 2026-09-15. The page lists 32 files. All 30 covering the visible are stored under their original names. The two VUV/UV-only sets, Myer & Samson 1970 and Roxlo & Mandl 1980, were skipped.
- **Contents:** the dataset list, with each set's T, resolution and notes, is in `docs/research/cross-section-data.md`.
- **Format:** CRLF line endings. Mostly two tab-separated columns, wavelength (nm) and σ (cm² molecule⁻¹), with no header. The exceptions are listed below.
  - `Bauer(1998)`, `Spietz(2006)_298K_500nm` and the three `Tellinghuisen(2011)` files with standard deviations have one to three header lines and a third column: the error limit or s.d.
  - The two `Tellinghuisen(2011)_..._400-500,600-850nm` files also have a blank line between the 400–500 and 600–850 nm blocks.
  - `JPL-2010(2011)_..._(max,min,rec)` has four columns, λ and σ at the band maxima and then at the minima.
  - In `JPL-2010(2011)_295K_185-700nm(rec)`, the 225 nm value is written `2.60e-18*)` and is followed by a `*) Corrected value` footnote. Strip the marker before parsing.
  - Some files end with a blank line.
  - Digitized historical curves were converted from decadic molar extinction coefficients with the factor 3.8235 × 10⁻²¹.
  - The atlas doesn't state whether wavelengths are in air or vacuum.
- **Terms:**
  - The live site states no licence. It asks users to cite Keller-Rudek et al. 2013.
  - The ESSD article is CC BY 3.0.
  - A frozen 2013 snapshot of the atlas is on Zenodo under **CC BY 4.0** (doi:10.5281/zenodo.6951). Zenodo refused the download with HTTP 403 ("unusual traffic"), so I have **not** checked that these files match the snapshot.
  - The underlying data remain the work of their original authors. The Saiz-Lopez and Spietz files, for example, come from CC BY-NC-SA 2.5 articles.

## Not obtained (see catalogs for details)

- **Knöckel et al. 2004 supplementary table (the Hannover fit data set):** unavailable. The paper points to `http://www.edpsciences.org`, which no longer resolves, and the SpringerLink page has no supplementary files (checked 2026-09-15). The data set has to be rebuilt from the sources the paper cites.
- **Salumbides et al. 2008 supplementary data:** also missing from SpringerLink (checked 2026-09-15). The potential parameters are printed in its Table 1.
- **Reiners et al. 2024 comb-calibrated FTS spectrum**, and the NIST FTS scans of astronomical iodine cells: not public, available on request.
- **Tellinghuisen 2011 papers** (J. Chem. Phys. 134, 084301 and 135, 054301): read from the published papers. The MPI-Mainz `Tellinghuisen(2011)_294-300K_520-635nm` file is the continuum absorptivity ε_c of Table II of the 134, 084301 paper (converted to σ); see `docs/research/intensity-inputs.md`.
- **MPI-Mainz 2013 Zenodo snapshot** (doi:10.5281/zenodo.6951, CC BY 4.0): the download was refused with HTTP 403 (checked 2026-09-15).

## `heeren_pyodine_atlas/` — two FTS iodine-cell atlases from the `pyodine` RV code

- **Source:** the `iodine_atlas/` directory of the pyodine repository, https://gitlab.com/Heeren/pyodine (branch `master`), fetched 2026-09-20 through the GitLab API. Paper: P. Heeren, R. J. Tronsgaard, F. Grundahl, S. Reffert, A. Quirrenbach, T. Marquart, "Pyodine: an open, flexible reduction software for iodine-calibrated precise radial velocities", *Astron. Astrophys.* **674**, A164 (2023), doi:10.1051/0004-6361/202244441. The repository's README asks that this paper be cited by anyone using the software.
- **`song_iodine_cell_01_65C.h5`:** the SONG (Tenerife) iodine cell at 65 °C, FTS transmission. Datasets `wavelength` (vacuum, Å), `wavelength_air`, `wavenumber` (cm⁻¹), `flux`, `flux_normalized`. 4990.0–6500.0 Å vacuum = 15 384.6–20 040.1 cm⁻¹, 618 003 points, 0.00244 Å ≈ 0.0075 cm⁻¹ per point.
- **`Fischer_Cell_May2022_downsampled3.h5`:** D. Fischer's cell (Lick/Yale), scanned May 2022, downsampled ×3. Same datasets. 4980.0–6250.0 Å = 16 000.0–20 080.3 cm⁻¹, 1 551 029 points, 0.00082 Å ≈ 0.0026 cm⁻¹ per point.
- **Calibration** (`docs/research/data-access-broadband.md` §1.1): both carry a ≈ −2×10⁻⁷ relative scale offset against Salami & Ross 2005, i.e. they were resampled/shifted for RV use; fit the scale out before using positions.
- **Terms:** **MIT** (`LICENSE`, copied here). Cite Heeren et al. 2023.

## `apo_pyodine_atlas/` — a NIST 2-m FTS scan of the APO iodine cell

- **Source:** the `iodine_atlas/` directory of https://github.com/jaklusmeyer/APO_pyodine (branch `main`, Git LFS), fetched 2026-09-20. The repository README describes it as "module for apo song pyodine — small changes from pyodine"; it builds on Heeren et al. 2023 (above), which should be cited for the software. **No licence file is present in that repository**, so the h5 files are used here for research only and are not redistributed.
- **`APO_I2_from_FITS_cleaned_newscan_vnorm.h5`:** a cleaned, normalised copy of an FTS scan that keeps its Xgremlin/NSO header as HDF5 attributes: `ORIGIN = NIST 2m IR-vis-UV FT`, `DATE = Jul 21, 2009`, `SAMPLE = I2, AS-17`, `RESOLUTN = 0.018 cm⁻¹`, `AIRCORR = No`, band 14 200–21 300 cm⁻¹. Datasets `wavelength_air` (Å, 4980.0–6500.0) and `flux_normalized`; 2 492 660 points. The FTS scale is HeNe-referenced: against Salami & Ross the line minima agree to ≲ 5×10⁻⁸ (§1.1 of the same doc), so this is an independent absolute position source once converted from air to vacuum (the header says `AIRCORR = No`, but the dataset name says air; the calibration check above was done with an air→vacuum conversion and is consistent).
- **`APO_I2_cleaned_normalized_shifted.h5`:** the same cell, cleaned, normalised and **shifted** for RV use: `wavelength_vac`, `wavelength_air`, `flux_raw`, `flux_normalized`; 346 647 points. Not a position reference.
- **Attribution for the scan itself:** NIST 2-m Fourier transform spectrometer (Gaithersburg), scan of 21 July 2009, sample AS-17. If positions from it are published, acknowledge NIST and the APO_pyodine repository; the scan has no paper of its own that we have found.

## `simmons_hougen_1977/jresv81An1p25.pdf`

- **Source:** J. D. Simmons, J. T. Hougen, "Atlas of the I₂ Spectrum from 19 000 to 18 000 cm⁻¹", *J. Res. Natl. Bur. Stand. A* **81A**(1), 25–80 (1977), doi:10.6028/jres.081A.005. Fetched 2026-09-20 from the Internet Archive, https://archive.org/download/jresv81An1p25/jresv81An1p25.pdf (56 pages, page images, 55 MB).
- **Contents:** a photographic-plate atlas of the B–X system, 18 000–19 000 cm⁻¹ in 20 cm⁻¹ pages, each with a densitometer trace and a table of line number, observed and calculated wavenumber and rotational/vibrational assignment. FWHM ≈ 0.055 cm⁻¹, absolute accuracy ±0.015 cm⁻¹ (thorium standards). Turning the tables into a line list means digitising the page images; the archive's OCR is lossy on them.
- **Terms:** a US Government work, **public domain**.

## `hal_open_access/` — the Gerstenkorn–Luc papers on HAL

Fetched 2026-09-20 from the HAL open archive (`https://hal.science/<id>/document`). EDP Sciences journals in HAL are free to read; no reuse licence is stated, so these are working copies for transcription only.

| file | paper |
|---|---|
| `jpa-00210032.pdf` | S. Gerstenkorn, P. Luc, "Description of the absorption spectrum of iodine recorded by means of Fourier Transform Spectroscopy: the (B-X) system", *J. Phys. (Paris)* **46**, 867–881 (1985), doi:10.1051/jphys:01985004606086700. **The 46 Dunham constants fitted to ~100 000 atlas transitions, 11 000–20 040 cm⁻¹: the starting point of our own fit.** |
| `jpa-00209975.pdf` | S. Gerstenkorn, P. Luc, C. Amiot, "Analysis of the long range potential of iodine in the B ³Π₀ᵤ⁺ state", *J. Phys. (Paris)* **46**, 355–364 (1985), doi:10.1051/jphys:01985004603035500. B-state IPA potential and long-range analysis. |
| `jpa-00244661.pdf` | S. Gerstenkorn, P. Luc, "Absolute iodine (I₂) standards measured by means of Fourier transform spectroscopy", *Rev. Phys. Appl.* **14**, 791–794 (1979), doi:10.1051/rphysap:01979001408079100. The −0.0056 cm⁻¹ atlas correction and its ±0.002 cm⁻¹ accuracy. |
| `jpa-00208967.pdf` | S. Gerstenkorn, P. Luc, J. Sinzelle, "Study of the iodine absorption spectrum by means of Fourier spectroscopy in the region 12 600–14 000 cm⁻¹", *J. Phys. (Paris)* **41**, 1419–1430 (1980), doi:10.1051/jphys:0198000410120141900. Band analysis of the hot-band region. |
| `jpa-00244946.pdf` | S. Gerstenkorn, P. Luc, A. Vetter, "Excitation spectrum of the iodine molecule induced by laser radiation in the 15 780–15 815 cm⁻¹ region: 'Complement to the Atlas of the iodine molecule spectrum'", *Rev. Phys. Appl.* **16**, 529–534 (1981), doi:10.1051/rphysap:01981001609052900. |

## `orsay_19700_20035/photos/` — Orsay atlas, Partie IV (19 700–20 035 cm⁻¹)

- **Source:** S. Gerstenkorn, P. Luc, *Atlas du spectre d'absorption de la molécule d'iode, Partie IV: 19 700–20 035 cm⁻¹* (Laboratoire Aimé Cotton, CNRS II, Orsay, 1983). Printed volume, library copy (QC462.I1 G47).
- **Photographed:** 2026-09-29, 159 phone photographs (`PXL_20260929_*.jpg`); the order of the file names is the order of the pages.
- **Contents:** introduction (pp. I–VII: 50 cm cell, p ≈ 0.25 torr; the scale is set by R(26) 62-0 at 19 926.0009(10) cm⁻¹ from Pique *et al.* 1981; estimated absolute accuracy 0.002 cm⁻¹; the wavenumbers are 0.0095(10) cm⁻¹ below those of Partie III), 33 plates with line tables (N, σ, ε, I; 3 613 lines with S/N > 1.5), and a classification table (N, σ_mes, I, J, branch, (v′, v″), σ_cal, σ_mes − σ_cal in mK).
- **Terms:** copyrighted print; photographs kept for transcription only, not redistributed.

## Derived: `data/atlas_lines/` (in the repository)

Line positions extracted from `salami_ross_2005/` and `apo_pyodine_atlas/` by `prototypes/atlas_lines.py`
and `prototypes/atlas_dataset.py` are committed as `data/atlas_lines/*.csv` with their own `.toml`
provenance. They are measurements derived here, not copies of the source files; each `.toml` carries the
citation of the spectrum it came from, which must accompany any use of them.
