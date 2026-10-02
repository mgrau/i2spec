# Line centres from the 1993 Kitt Peak cell scans and the Rodríguez Fernández laser spectra

*2026-09-30. `prototypes/kitt_peak_lines.py`, `prototypes/rodriguez_lines.py` (both reuse
`prototypes/atlas_lines.py`). Products: `data/atlas_lines/kitt_peak_1993.{csv,toml}` (5 297 lines) and
`data/external/rodriguez_fernandez_marcassa/rodriguez_fernandez_2023.{csv,toml}` (825 lines; kept out of the
public repository, since the spectra were shared privately by their authors). Model: i2spec2026m.*

## The data

**Kitt Peak** (`data/external/kitt_peak_fts/`). McMath-Pierce 1-m FTS, 17 March 1993, observer Marcy:
twelve scans of radial-velocity iodine cells and one of a bromine cell (`930317R0.010`, skipped). Every
file has the same axis, WSTART + i DELW = 14 951.965 + 0.01409414 i cm⁻¹ (397 312 points, big-endian
float32). The headers give RESOLUTN 0.037 cm⁻¹, 10 co-added scans, a DZE lamp (not divided out), an
8 mm aperture, an evacuated tank, AIRCORR = No (vacuum wavenumbers) and **WAVCORR = 0**: the archive
applied no scale correction. The useful band is 15 875–20 000 cm⁻¹. `930505R0.015` (Xe lamp, 0.15 cm⁻¹
resolution) is too coarse for line centres and was not used.

The published set is the **Keck cell at 50.0 °C, 100 mm** (`930317R0.003`), fitted with full window
coverage. It is the best-known of the cells (the HIRES template of Butler et al. 1996), and at 50 °C
its lines are the least saturated of the Keck series. Three other scans were fitted with windows
every 6 cm⁻¹, only to test the scale: Keck 70 °C (`.005`), ESO 50 °C (`.008`) and Lick backup 50 °C
(`.001`). The cells hold the same iodine and were measured the same day with the same instrument, so
a second full set would repeat the same lines with the same template error and add nothing
independent.

**Rodríguez Fernández / Lefrán Torres** (`data/external/rodriguez_fernandez_marcassa/`).
Doppler-limited diode-laser transmission, 14 400–14 710 cm⁻¹, read on a WS7 wavemeter, in 0.5 cm⁻¹ scans
joined into two files. The first file's axis steps backwards at 207 joins, by up to 0.0045 cm⁻¹. It is
sorted here, and those steps mark scan boundaries. The other joins are invisible, because they happen
to be monotonic. The second file is monotonic throughout.

## How the centres were fitted

As for the NIST/APO scan (`atlas-line-positions.md`): the i2spec model transmission is the template,
and every line deep enough to matter gets its own wavenumber shift. Column density, baseline, zero
offset and a Gaussian width are fitted per window, with hyperfine structure in the template and a
0.01 cm⁻¹ ridge prior on each shift. Repeats are combined by inverse variance. Differences from the
APO run:

- **Kitt Peak:** each scan is divided by its upper envelope (99.5th percentile in 4 cm⁻¹ blocks). The
  windows are 3.2 cm⁻¹ wide, stepped by 2.4 cm⁻¹, because the extractor needs 200 points and the
  sampling is 0.0141 cm⁻¹. The cell is at T = 323.15 K. The axis is divided by (1 − 1.33 × 10⁻⁶) before
  the fits and restored after, because a first window put the lines 700 MHz low, and the 300 MHz ridge
  prior would have pulled the weak lines toward the model. Kept: σ_fit < 30 MHz, window rms < 4 %,
  depth > 0.05. 5 899 lines passed the cuts, of which 5 297 remain after leaving out blends.
- **Rodríguez:** the windows are 1.0 cm⁻¹ wide, stepped by 0.8 cm⁻¹, so each spans only a few scans.
  The axis is shifted by −131 MHz before the fits and restored after. The cell temperature is not in
  the files. A scan of the mean window rms over 295–500 K gives a minimum at **350 K** (0.79 %, against
  0.83 % at 375 K and 1.3 % at 295 K). Kept: σ_fit < 30 MHz, rms < 4 %, depth > 0.03. 966 lines passed
  the cuts, of which 825 remain after leaving out blends.
- **Blends.** A line is left out if another model line within ±0.02 cm⁻¹ (Kitt Peak, half the
  resolution) has at least 0.1 of its strength, or within ±0.01 cm⁻¹ (Rodríguez) has at least 0.3.
  The covariance does not catch these, because a neighbour below the fit threshold is held at its
  model position. The thresholds come from the residuals. At Kitt Peak, lines with such a neighbour
  scattered 28–34 MHz about the scale curve (against 21 MHz for clean lines) and were biased by −4 to
  −14 MHz. In the Rodríguez set they read +6 to +9 MHz (against −1 MHz for clean lines).
- **Uncertainties.** Kitt Peak: σ = √(σ_fit² + 6.3²) MHz, from the scatter against APO in bins of
  σ_fit after removing the scale curve. This gives a median σ of 24 MHz, and z rms 1.00 against APO.
  Rodríguez: σ = √((2.69 σ_fit)² + 14.3²) MHz, fitted to the scatter about the model, which its
  precision sets pin to about 1 MHz here. This gives a median σ of 16 MHz, and z rms 1.03. Salami &
  Ross cannot calibrate this set: in this range they scatter 33 MHz about the model, more than their
  own σ, and would give a 27 MHz floor.

## Kitt Peak: the scale is not a single number

Against i2spec2026m (5 297 lines):

| region (cm⁻¹) | n | median (MHz) | as a scale (ppb) | spread (MHz) |
|---|---|---|---|---|
| 15 500–16 000 | 145 | −580 | −1212 | 22 |
| 16 000–16 500 | 520 | −606 | −1240 | 21 |
| 16 500–17 000 | 510 | −645 | −1285 | 26 |
| 17 000–17 500 | 571 | −693 | −1340 | 25 |
| 17 500–18 000 | 800 | −746 | −1401 | 27 |
| 18 000–18 500 | 817 | −808 | −1475 | 24 |
| 18 500–19 000 | 1047 | −873 | −1553 | 34 |
| 19 000–19 500 | 812 | −931 | −1618 | 32 |
| 19 500–20 000 | 75 | −996 | −1698 | 66 |

The scale, obs = model (1 + s(ν)), is not constant: it runs from −1210 to −1700 ppb across the band.
The median is −784 MHz (robust spread 145 MHz).

| fit | parameters (ppb) | residual spread | 250 cm⁻¹ bin medians rms (max) |
|---|---|---|---|
| constant | −1454 ± 16 | 90 MHz | 80 MHz (143) |
| linear | c0 −1436, c1 −129 | 22 MHz | 9 MHz (23) |
| **quadratic** | **c0 −1422.0 ± 1.6, c1 −133.5 ± 1.2, c2 −14.3 ± 1.2** | **20.5 MHz** | **3.9 MHz (11)** |
| cubic | c3 +0.4 ± 1.4, otherwise the same | 20.5 MHz | 4.0 MHz |

Here s = c0 + c1 x + c2 x² with x = (ν − 17 900 cm⁻¹)/1000 cm⁻¹. The errors come from a bootstrap over
50 cm⁻¹ blocks. Once the quadratic is removed, the medians per 500 cm⁻¹ region are within ±2 MHz
(except −11 MHz at the 75 lines above 19 500 cm⁻¹), and the spread is 16–27 MHz.

**The curve belongs to the instrument, not to the model.** Kitt − APO shows the same curve, measured
without the model: c0 −1388.6, c1 −135.1, c2 −13.8 ppb. The 33 ppb difference in c0 is APO's own
−18 MHz offset from the model. The other three cells give the same shape, c1 = −134 to −138 ppb and
c2 = −12 to −14 ppb, but c0 moves from scan to scan: −1437 (Keck 70 °C), −1404 (ESO) and −1452 ppb
(Lick). That is ±25 ppb, or ±15 MHz. A misset reference laser or the aperture would only make a
constant scale. The aperture alone predicts −(d/2f)²/4 = −560 ppb, much less than what is found.
Where the ν-dependence comes from is not known. The phase correction was off (MODEPHZC = 0) and the
line shape is not modelled, so a line-shape asymmetry that varies with wavenumber is one candidate.
This is inference only.

**Against the other atlases** (after the quadratic):

| pair | n | median | robust spread |
|---|---|---|---|
| Kitt Peak − `apo_nist_2009` | 4 893 | +17.7 MHz | 19.6 MHz (16–24 by 1000 cm⁻¹ region) |
| Kitt Peak − `salami_ross_2005` | 5 096 | −11.8 MHz | 25.7 MHz |
| `apo_nist_2009` − model, same lines | 4 893 | −17.6 MHz | 4.0 MHz |
| `salami_ross_2005` − model, same lines | 5 096 | +12.3 MHz | 14.9 MHz |

After its calibration, Kitt Peak agrees with both atlases to within their known offsets, but it scatters
five times more than APO: 0.037 cm⁻¹ resolution against 0.018 cm⁻¹, and most lines are over 90 % deep.
Its coverage is v′ = 7–51, v″ = 0–4, 77 bands. APO holds 4 893 of its 5 297 lines. Only **158 are in
neither atlas**, and 404 are not in APO.

## Rodríguez Fernández: one offset, and scan offsets in the first file

Against i2spec2026m (825 lines, v′ = 2–5, v″ = 6–8, J = 9–139, 9 bands):

| region (cm⁻¹) | n | median | robust spread |
|---|---|---|---|
| 14 400–14 500 | 276 | +113.9 MHz | 21.3 MHz |
| 14 500–14 600 | 148 | +107.0 MHz | 18.2 MHz |
| 14 600–14 710 | 401 | +114.2 MHz | 14.2 MHz |
| all | 825 | **+113.0 ± 2.0 MHz** | 18.3 MHz |

The ± 2.0 MHz comes from a bootstrap over 5 cm⁻¹ blocks. The weighted mean is +110.2 MHz. A linear
term fits +0.04 MHz per cm⁻¹, i.e. +12 MHz across the range: small, and not needed. The earlier
parabola estimate (+131 MHz, spread 65 MHz, with +148/+113/+138 by region; `data/external/README.md`)
took the transmission minimum of unresolved hyperfine clusters. The template fit removes that bias
and cuts the spread by more than a factor of three. The sign and size agree with the authors' own
comparison: 105–150 MHz above Gerstenkorn & Luc and Salami & Ross. Here, this set − Salami & Ross is
+129.7 MHz (478 lines, spread 38.6 MHz), and Salami & Ross − model is −18.8 MHz on the same lines.

**Per-scan offsets.** In the first file (the 2023 paper, 14 400–14 600 cm⁻¹), the segments between
backward steps of the axis differ in offset beyond their line noise: χ²/dof = 6.5 between the 52
segments that hold ≥ 3 lines, and about 16 MHz of extra scatter. The second file (the 2022 paper)
shows none: χ²/dof 0.72 over 0.5 cm⁻¹ blocks. The offsets were not fitted. Most segments hold one or
two lines, the other joins are invisible, and an offset fitted per scan against the model would make
those lines model-relative. The scatter is in the per-line floor instead, and each line's note names
its file and segment, so a fit can add per-segment offsets.

Of the 825 lines, 347 are not in `salami_ross_2005`. This range is v″ = 6–8, where Salami & Ross is
our only other atlas and scatters 33 MHz about the model.

## Recommendation

- **`rodriguez_fernandez_2023`: worth adding to the fit and the explorer if the authors agree**
  (it stays private until then, and is used only for the comparison above), with one group offset
  (+113 MHz), or one per first-file segment if the fit can carry them. It is better than Salami &
  Ross in this range (18 against 33 MHz about the model), has one clean calibration, and adds 347
  lines in hot bands (v″ = 6–8) that no other data set here measures.
- **`kitt_peak_1993`: not worth adding to the fit.** It needs a three-parameter scale curve, and after
  the curve it is five times noisier than APO. APO already holds 92 % of its lines, and only 158 are
  new. At most it is a low-weight explorer layer, or an independent check of APO on the 4 893 shared
  lines (it confirms APO to 18 ± 20 MHz). If any Kitt Peak scan is used, its constant term must be
  fitted separately: the scans differ by ±25 ppb.

*Caveat: the model comparisons used i2spec2026m as it stood on 2026-09-30. Other work was changing
`src/i2spec/model.py` and `level_corrections_2026m.json` at the same time, so the model residuals may
move by a few MHz. The line centres themselves do not depend on that, beyond the template start.*

## The archive catalogue: every iodine spectrum at Kitt Peak (2026-10-01)

To close out the archive, the header of every full-resolution spectrum in the NSO FTS archive
(`https://nispdata.nso.edu/ftp/FTS_cdrom/`, volumes FTS01–57) was read. NSO describes the archive as not
searchable and the old Digital Library query tool is gone, so the headers were taken from the files:
volumes FTS01–54 by streaming the archive's one tar file (`FTS_01_55.tar`, 31.8 GB, which despite its name
ends at FTS54) once and keeping only the first 5 760 bytes of each member, and FTS55–57 one header per
request, 15 s apart. 24 431 full-resolution spectra in all; none was refused.

136 headers mention iodine. Removing quartz–iodine lamp calibrations, solar and hollow-cathode runs with an
I₂ cell in the beam, and runs outside 11 000–20 100 cm⁻¹ leaves the 1993 scans used above and the 41
laboratory iodine-cell spectra below. All are in `data/external/kitt_peak_fts/` (97 MB with the 1993 scans,
not redistributed). Ranges in cm⁻¹ and resolution (RESOLUTN, cm⁻¹) are from the headers; the ID is the
observer's.

| file | volume | range | resolution | ID |
|---|---|---|---|---|
| `770611R0.003` | FTS01 | 14432–18473 | — | CALIBRATION RUN - KR 86 + I2 CELL NEAR 6300A |
| `770309R0.004` | FTS01 | 15874–22117 | — | IODINE 18000.-19000. |
| `810326R0.012` | FTS05 | 18933–27026 | 0.0283 | I2 CELL NO. 1, 18O2 FOR IO IN BLUE |
| `810914R0.016` | FTS06 | 15269–21273 | — | I2 CELL WITH 5000-6000A INTEGRATED LITE |
| `820620R0.052` | FTS07 | 15422–20426 | 0.0089 | INTEGRATED LITE, I2 CELL |
| `820619R0.009` | FTS07 | 15422–20426 | 0.0089 | INTEGRATED LITE, IODINE CELL |
| `840611R0.002` | FTS12 | 13786–20612 | 0.0138 | IODINE 1.5 INCH CELL, 5050-6650A |
| `840531R0.001` | FTS12 | 13786–20612 | 0.0182 | IODINE 1.5 INCH CELL, 5050-6650A |
| `840531R0.002` | FTS12 | 13786–20612 | 0.0182 | IODINE 1.5 INCH CELL, 5050-6650A |
| `840629R0.001` | FTS12 | 8918–17802 | 0.0200 | IODINE, ARGON-ION LASER, 5145 A, 3.2 WATTS |
| `841018R0.006` | FTS13 | 15511–23085 | 0.0772 | IODINE ABSORPTION, 2 INCH CELL,DZE AT 150 WATTS |
| `841018R0.007` | FTS13 | 15511–23085 | 0.1554 | IODINE ABSORPTION, 2 INCH CELL,DZE AT 150W, LOW RES. |
| `881012R0.002` | FTS20 | 13991–20945 | 0.0098 | I2 CELL AT 23.6 C, 15000-20400 CM-1 |
| `881012R0.003` | FTS21 | 13991–20945 | 0.0098 | I2 CELL, 25CM., 3.6-4.3 C 15000-20400 CM-1 |
| `881015R0.018` | FTS21 | 15191–16309 | 0.0098 | I2 CELL 25 CM,. 23.9 C, 15750-15850 CM-1 |
| `881015R0.017` | FTS21 | 16882–17618 | 0.0098 | I2 25 CM. CELL, 23.9 C, 17180-17380 CM-1 |
| `881012R0.004` | FTS21 | 17482–19718 | 0.0098 | I2 CELL, 3.5-3.8 C, NARROW BAND18400-18900 CM-1 |
| `881014R0.012` | FTS21 | 16500–17127 | 0.0098 | I2 25 CM CELL, 16.1-16.4 C, 16700 - 16900 CM-1 |
| `881015R0.014` | FTS21 | 16500–17127 | 0.0098 | I2 CELL, 25 CM., 23.9 C, 16700-16900 CM-1 |
| `881013R0.006` | FTS21 | 16678–17823 | 0.0098 | I2 CELL 13.9-14.2 C, 17180-17380 CM-1 FOR STD RUN |
| `881015R0.015` | FTS21 | 16145–17754 | 0.0098 | I2 CELL, 25 CM., 23.9 C, 16890-17100 CM-1 |
| `881014R0.008` | FTS21 | 16145–17754 | 0.0098 | I2 CELL 7.8-8.0 C, 16890-17100 CM-1 |
| `881014R0.009` | FTS21 | 16500–17127 | 0.0098 | I2 25 CM CELL, 7.8-8.0 C, 16700- 16900 CM-1 |
| `881014R0.007` | FTS21 | 16800–17727 | 0.0098 | I2 CELL 7.8-8.0 C, 17180-17380 CM-1 |
| `881014R0.013` | FTS21 | 16145–17754 | 0.0098 | I2 25 CM CELL, 16.1-16.4 C, 16890-17100 CM-1 |
| `881014R0.011` | FTS21 | 15682–16827 | 0.0098 | I2 25 CM CELL, 16.1-16.4 C, 16250-16450 CM-1 |
| `881014R0.010` | FTS21 | 17700–19527 | 0.0098 | I2 25 CM CELL, 7.7-7.9 C, 18550-18800 CM-1 |
| `881015R0.016` | FTS21 | 15682–16827 | 0.0098 | I2 25 CM CELL, 23.9 C, 16250-16450 CM-1 |
| `901126R0.002` | FTS27 | 14952–21014 | 0.0370 | IODINE CELL 10CM .01 TORR 50C, 5000 - 6200 A |
| `901126R0.001` | FTS27 | 14952–21014 | 0.0370 | IODINE CELL 10CM .01 TORR 50C, 5000 - 6200 A |
| `951009R0.001` | FTS38 | 13760–20639 | 0.0204 | Lick I2 10cm cell 50C 5000 - 6300A |
| `951009R0.002` | FTS38 | 13760–20639 | 0.0204 | Lick I2 10cm cell 50C 5000 - 6300A |
| `950626R0.001` | FTS38 | 15191–20611 | 0.0204 | Optronics Lamp @ 15Amps long I2 cell, 500-600nm 2x8mm |
| `950626R0.002` | FTS38 | 15191–20611 | 0.0204 | Optronics Lamp @ 15Amps short I2 Cell 500-600 nm 2x8mm |
| `950626R0.003` | FTS38 | 15191–20611 | 0.0204 | Optronics Lamp @ 15Amps no cell 500-600nm 2x8mmm 2x8mm |
| `960405R0.022` | FTS40 | 8355–19833 | 0.0300 | O2 299torr 21.1C 2.4m, I2 10cm, 9000 - 19000 cm-1 |
| `960405R0.017` | FTS40 | 8355–19833 | 0.0351 | O2 572torr 21.4C 2.4m, I2 10cm, 9000 - 19000 cm-1 |
| `960405R0.016` | FTS40 | 8355–19833 | 0.0351 | O2 404.6torr 20.3C 2.4m, I2 10cm, 9000 - 19000 cm-1 |
| `010712R0.027` | FTS52 | 14995–29991 | 0.0499 | Empty 4.01M JPL NO2 cooled cell 17000 - 26000 cm-1 + Iodine cell |
| `010713R0.044` | FTS52 | 14995–29991 | 0.0499 | Empty 4.01M JPL NO2 cooled cell 17000 - 26000 cm-1 I2 spectrum |
| `010713R0.054` | FTS52 | 14995–29991 | 0.0499 | 4.01M,0.05T NO2 + air: 100 T, T=-57 C; 17000-26000 cm-1 + I2 Spec. |

The 1993 scans (`930317R0.001`–`.013`, FTS33) are described in the sections above.

**What they could add.**

- **The October 1988 series (FTS20/21, 16 spectra).** One 25 cm cell at 3.6, 7.9, 14.1, 16.2 and 23.9 °C,
  a broadband scan (15 000–20 400 cm⁻¹) and narrow windows, all at 0.0098 cm⁻¹. With a known path and
  several temperatures this is the one set in the archive that can test the line strengths and their
  temperature dependence (`bx-band-strength.md`), which no other data set has tested.
- **The 1995 Lick cell (FTS38).** 10 cm at 50 °C at 0.020 cm⁻¹, with long and short cells and a scan without
  a cell for the baseline: a second epoch of the radial-velocity cells, at better resolution than 1993.
- **The 1984–1990 cells (FTS12, FTS13, FTS27)** at 0.014–0.077 cm⁻¹: more of the same band, useful mainly as
  checks on the scale.
- The rest are lamp or calibration runs, or iodine scanned alongside other gas cells (O₂, NO₂), with little
  to add.

None has been fitted yet.
