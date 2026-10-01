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
