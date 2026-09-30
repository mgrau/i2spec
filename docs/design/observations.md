# Observation data sets (M4)

Every measurement used by the global fit (stage 5) is stored in one plain-text format, with its provenance. The code is in `src/i2spec/observations.py`.

## Layout

```
data/observations/<id>/meta.toml   source, provenance, unit, conditions
data/observations/<id>/data.csv    one observation per row
```

Where the source has a `data/catalog/*.yaml` entry, `<id>` is the catalog id (for example `bipm2012a`).

Why these formats:
- **TOML** can be read with the standard library (`tomllib`), so there is no new dependency.
- **CSV** diffs cleanly and opens anywhere.

The files are not shipped in the Python package. The fit runs from a checkout.

## `meta.toml`

| Key | Required | Meaning |
|---|---|---|
| `id` | yes | Must equal the directory name |
| `citation` | yes | Full reference |
| `source` | yes | URL or DOI that the numbers were taken from |
| `retrieved` | yes | Date the source was obtained |
| `unit` | yes | `"MHz"` or `"cm-1"`, for both `value` and `uncertainty` |
| `transcription` | yes | How the numbers were obtained (typed, parsed, digitized), what was checked, and every correction or relabelling applied |
| `catalog` | no | Pointer to the catalog entry |
| `partial` | no | What is still missing, if the set is incomplete |
| `uncertainty` | no | How the stated uncertainties were read (1σ or not, statistical or total) |
| `[conditions]` | no | Cell and laser conditions, as printed |
| `exclude` | no | Reason to leave the set out: `load_all()` skips it (and so the fit, the explorer and the validation scripts), `load_all(include_excluded=True)` keeps it |
| `license`, `notes` | no | |

## `data.csv`

Columns, in this order: `line,component,kind,value,uncertainty,ref_line,ref_component,group,note`.

| Column | Meaning |
|---|---|
| `line` | `"<isotopologue> <P or R>(<J''>) <v'>-<v''>"`, e.g. `127I2 R(56) 32-0`. The isotopologue is `127I2`, `129I2` or `127I129I` |
| `component` | Hyperfine component label as published (`a10`, `b12`, `m48`). Empty means the hyperfine-free line centre |
| `kind` | `frequency`: the value is f(line, component). `interval`: the value is f(line, component) − f(ref_line, ref_component) |
| `value`, `uncertainty` | In the data set's unit. The uncertainty is the standard (1σ) uncertainty, read as described in `meta.toml` |
| `ref_line`, `ref_component` | The reference of an interval. An empty `ref_line` means the same line; an empty `ref_component` means the reference line's hyperfine-free centre |
| `group` | Label shared by observations that have a common systematic, such as one spectrum's calibration offset. The fit gives each group its nuisance parameters |
| `note` | Free text |

**Component labels.**
- The number is the component's rank by frequency among the main (ΔF = ΔJ) components, with 1 the lowest.
- The letter only names the line in the source; BIPM uses a, b, d, e and m for different lines.
- i2spec matches components by rank (`component_rank`), using the strongest-first ΔF = ΔJ assignment in `hyperfine.line_components`.
- Crossover resonances and ΔF ≠ ΔJ components cannot be represented yet.

**Frequencies** are for the unperturbed transition as the source publishes them. Any correction made during transcription, such as pressure or power shifts or a calibration update, must be stated in `transcription`.

## Python

```python
from i2spec.observations import Predictor, load_all, residuals

predictor = Predictor()                  # one RovibronicModel per isotopologue, hyperfine with ΔJ = ±2
for ds in load_all():
    r = residuals(predictor, ds)         # observed − model, in the data set's unit
```

`load_dataset` validates on read:
- the column order;
- line keys, component labels and kinds;
- that intervals have a reference;
- that uncertainties are positive.

Errors name the file and the row. `write_observations` writes the CSV back, and the tests check the round trip.

## Data sets

| id | Content | Rows | Status |
|---|---|---|---|
| `bipm2005a` | BIPM 515 nm (MEP 2005): f(a3, P(13) 43-0); hyperfine tables of P(13) 43-0, R(15) 43-0 and R(98) 58-1; the two intervals between lines given in the notes | 57 | complete |
| `bipm2012a` | BIPM 532 nm (CIPM 2007, file 2012): f(a10, R(56) 32-0); a1/a21 intervals of 19 lines to it; hyperfine tables of all 20 lines | 329 | complete; R(121) 35-0 a13 (a misprint) left out |
| `bipm2003b` | BIPM 543 nm (MEP 2003): f(b10, R(106) 28-0); hyperfine tables of R(106) 28-0 and R(12) 26-0 (relative to b10) | 30 | complete |
| `bipm2003c` | BIPM 576 nm (MEP 2003): f(a1, P(62) 17-1); its hyperfine table (a2–a10, a15) | 11 | complete |
| `bipm2003d` | BIPM 612 nm (MEP 2003): f(a7, R(47) 9-2); hyperfine tables of ¹²⁷I₂ R(47) 9-2, P(48) 11-3, R(48) 15-5 and ¹²⁹I₂ P(110) 10-2, R(113) 14-4, all relative to a7 | 90 | complete |
| `bipm2003a` | BIPM 633 nm (MEP 2003): f(a16, R(127) 11-5); hyperfine tables of ¹²⁷I₂ R(127) 11-5, P(33) 6-3; ¹²⁹I₂ P(54) 8-4, P(69) 12-6, R(60) 8-4, P(33) 6-3; ¹²⁷I¹²⁹I P(33) 6-3; one interval between ¹²⁷I₂ lines and two isotope shifts from the notes | 158 | complete; 14 entries for 6 features blended between ¹²⁹I₂ and ¹²⁷I¹²⁹I left out (see `meta.toml`) |
| `bipm2003e` | BIPM 640 nm (MEP 2003): f(a9, P(10) 8-5); hyperfine tables of P(10) 8-5 and R(16) 8-5 (b1–b3) | 18 | complete |
| `bodermann2000a` | Hannover NIR (778–795 nm): 29 absolute hyperfine-component frequencies of 21 lines in bands 0-12, 0-13, 0-14, 0-15, 1-14, 2-15 and 3-16, from the Rb D1, Rb D2, Rb two-photon and four-wave-mixing calibrations; the 3 FWM interval chain steps of §6 | 32 | complete; 4 rows (group `rb_2photon_derived`) are the paper's own re-derivations of other rows, see `meta.toml` |
| `liao2010a` | Liao 2010 (750–780 nm): 27 hyperfine components of nine (0-12)/(0-13) lines (Table 2), three (0-14) cell-calibration lines (Table 1), and the R(141) 0-16 a2 comb value from §2 | 31 | complete; Table 2's "P(100) 0-12" entered as R(100) 0-12 (see `meta.toml`) |
| `velchev1998a` | Velchev 1998 (571-596 nm): 115 measured "t" hyperfine components, ¹²⁷I₂, v' = 13-18, v'' = 1, J'' = 6-137, 1.1-4.5 MHz | 115 | complete for Table 1; the 1584 *predicted* components were supplementary data and are lost |
| `xu2000a` | Xu 2000 (595-655 nm): 473 "t" hyperfine components, ¹²⁷I₂, 15 bands with v' = 5-13, v'' = 2-5, all at 1.0 MHz | 473 | 473 of 481; 4 footnote-b rows excluded (different components), about 4 unreadable in the scan |
| `cornish2000a` | Cornish 2000 (730-733 nm): 2 absolute a15 frequencies, ¹²⁷I₂ R(26) 5-13 and P(258) 7-11, 0.14 and 1.2 MHz | 2 | 2 of 3; R(137) 5-12 is a blended a19-21 and cannot be given a rank |
| `reinhardt2007a` | Reinhardt 2007 (735/772/780 nm): 3 comb-referenced a1 frequencies at 0.30 MHz, v'' = 11 and 14, plus 4 hyperfine intervals at 42 kHz | 7 | complete |
| `huet2013a` | Huet 2013 (716 nm): 14 hyperfine intervals of R(90) 3-10 relative to a1, 0.28-0.38 MHz | 14 | complete; the paper measures no absolute frequency |
| `reinhardt2006a` | Reinhardt 2006 (565/576/585 nm): 56 hyperfine intervals of five lines at 20-30 kHz (v' = 14, 15, 17, 21; v'' = 1) plus one comb-referenced absolute | 57 | complete; Table 4 is referenced to a3, not a1 |
| `dube2004a` | Dubé 2004 (715-719 nm): 62 hyperfine intervals of four lines, v'' = 9 and 10, 0.31-0.48 MHz | 62 | complete; sigma_calib is a per-line scale systematic and is carried as the group label, not as a row uncertainty |
| `morinaga1989a` | Morinaga 1989 (657 nm): 20 intervals of R(69) 3-4 relative to component a, 1-3 MHz | 20 | complete for Table 1; Table 2 and the Table 3 hyperfine constants are not transcribed |
| `yoshiki2023a` | Yoshiki 2023 (514 nm): a1/b1 absolute frequencies of 7 lines at 5.6 kHz and 72 hyperfine splittings of five v′ = 44 lines at 0.4 kHz | 79 | complete; 4 starred (overlapping) rows left out, as the authors did |
| `matsunaga2024a` | Matsunaga 2024 (514 nm): a1/b1 absolute frequencies of 6 lines, v′ = 45–50, at 5.6 kHz and 91 splittings at 0.6 kHz | 97 | complete; 5 starred rows left out |
| `kobayashi2016a` | Kobayashi 2016 (578 nm): a1 absolute frequencies of 4 lines, bands 16-1, 17-1, 18-1, at 7 kHz and 72 splittings at 2 kHz | 76 | complete; 3 unobserved and 5 zero-weight components left out |
| `nishiyama2024a` | Nishiyama 2024 (520 nm): a1 absolute frequencies of 4 lines, v′ = 39, at 5.1–5.5 kHz and 68 splittings at 0.6 kHz | 72 | complete |
| `sansonetti1997a` | Sansonetti 1997 (560–656 nm): one hyperfine component of each of 102 lines, v′ = 5–22, v″ = 0–5, at 0.3–0.8 MHz, from the NIST/LLNL report version | 102 | complete; three misprinted classifications corrected by the model (see `meta.toml`) |
| `bodermann1998b` | Bodermann 1998b (815 nm): 3 absolute frequencies at 14-17 kHz and one 15 kHz splitting, v'' = 16 and 17 | 4 | complete; the paper's derived a10 absolute is omitted as redundant |
| `nesterenko2019` | Nesterenko 2019 (1053-1068 nm): 18 B->X **emission** frequencies at 2-235 kHz, bands 32-53 and 32-54 | 18 | complete; **v'' = 53 and 54**, far beyond anything else held |
| `matyugin2012` | Matyugin 2012 (982-985 nm): B->X **emission** frequencies to v'' = 48 at 3-65 kHz, from B v' = 32 and 33 | 18 | 18 of 20; assignment confirmed by the 2008 companion paper (below) |
| `ikeda2022a` | Ikeda 2022 (514 nm): a1/b2 absolute frequencies of P(57) 45-0, R(73) 46-0 and P(91) 48-0 at 5.4 kHz and 56 splittings at 0.4 kHz | 59 | complete; 2 zero-weight blended and 2 unmeasured components left out; the only absolute at v′ = 48 |
| `cheng2019a` | Cheng 2019 (534 nm): R(53) 31-0 a21 absolute (corrected to zero pressure and power) at 14 kHz and 18 intervals to it | 19 | complete; a5/a6 overlap and were not measured |
| `shie2013a` | Shie 2013 (535 nm): a1, a10, a15 absolute frequencies of P(28) 30-0 at 11 kHz (zero-pressure values) | 3 | complete |
| `hsiao2013a` | Hsiao 2013 (548 nm): a1, a10, a15 absolute frequencies of P(28) 24-0 at 12 kHz (zero-pressure values) | 3 | **excluded** since i2spec2026n (`exclude` in `meta.toml`): 1.1–1.4 MHz from one B v′ = 24 correction that holds hauden2024a, tanabe2022a and yang2011a to 0.1–0.2 MHz, and its own a10−a1, a15−a1 intervals 0.20–0.24 MHz from a hyperfine model that fits hauden2024a's 21 components at the same v′ to 0.03 MHz |
| `zhang2009a` | Zhang 2009 (560 nm): a1, a10, a15 of R(34) 20-0 and a1 of P(144) 23-0, absolute, 4–71 kHz (printed 2σ halved) | 4 | complete |
| `yang2011a` | Yang 2011 (561 nm): P(58) 22-1 a15 and R(130) 24-1 a1 absolute at 4 and 12 kHz; 63 splittings of R(62)/P(58) 22-1, R(103) 23-1 and R(130) 24-1 | 65 | R(103) 23-1 a1 absolute left out: 840 MHz from both the model and the NIST/APO atlas |
| `edwards1999a` | Edwards 1999 (633 nm): P(33) 6-3 b21 absolute at 17 kHz and b1–b20 relative to b21 at 5–14 kHz, 2 MHz modulation | 21 | Table 2 (6 MHz modulation, same intervals remeasured) not entered |
| `huang2013a` | Huang 2013 (671 nm): a1, a10, a15 absolute frequencies of R(78) 4-6 at 24–35 kHz (zero-pressure values); the only precision line between 667 and 716 nm | 3 | complete; the three printed splittings are differences of these and are not entered |
| `hong2009a` | Hong 2009 (578 nm): R(37) 16-1 a1 absolute at 2 kHz (2.4 Pa, 2.1 mW; not extrapolated to zero pressure) | 1 | complete; Kobayashi 2016 finds it 5 kHz above their value moved to the same conditions |
| `ye1999a` | Ye 1999 (532 nm): gaps from R(56) 32-0 a10 to 8 lines (Fig. 3), at 5 kHz (not printed; judgement); hyperfine intervals of P(54) 32-0 and R(57) 32-0 to a1 (Table I) at 1 kHz | 40 | complete; **not independent of `bipm2012a`**, whose Tables 9–10 are these intervals rounded |
| `holzwarth2001a` | Holzwarth 2001 (532 nm): comb-measured a1 absolutes of 16 lines, P(51)–P(58)/R(54)–R(61) 32-0 and P(83) 33-0, plus R(56) 32-0 a10, at 5.2 kHz, at −5 °C (2.42 Pa, uncorrected) | 17 | complete; misprinted P(51) 32-0 value corrected (digit groups swapped); **a source of `bipm2012a`** |
| `jones2002a` | Jones 2002 (515/532 nm): P(13) 43-0 a3 absolute at 2.38 Pa (not pressure-corrected), 1.5 kHz; R(56) 32-0 a10 from the 750 MHz comb check (Fig. 2), repeatability only | 2 | complete; the "378.8 kHz as measured" in BIPM MEP 2005 is a misprint of Fig. 5, the value is 441.8 kHz; **a source of `bipm2005a`** |
| `goncharov2004a` | Goncharov 2004 (515 nm): P(13) 43-0 a3 absolute at 0.12 Pa (as measured), 0.75 kHz | 1 | complete; the paper's extrapolations to 2.38 Pa are derived and not entered; **a source of `bipm2005a`** |
| `simonsen2000a` | Simonsen 2000 (633 nm): hyperfine intervals of ¹²⁷I₂ P(33) 6-3, R(60) 8-4, R(125) 9-4, P(54) 8-4, R(39) 6-3, R(59) 8-4, P(53) 8-4 at 1–10 kHz, and 3 links between lines across ±20 GHz at 5–6 kHz | 122 | complete; 14 blended entries and the MP-97-derived link to R(127) left out; the abstract's "P(39)" is R(39) 6-3 |
| `huang2018a` | Huang 2018 (647 nm): a1, a10, a15 absolute frequencies of P(46) 5-4 at 21 kHz (zero-pressure values) | 3 | complete; the two printed splittings are differences of these and are not entered |
| `manzoor2024a` | Manzoor 2024 (652 nm): 9 comb-referenced absolute frequencies of components of P(63) 4-4 (a1–a4, a7, a12–a14, a18) at 0.19–1.3 MHz | 9 | complete; the fitted centre of gravity and ΔeqQ/ΔC (Table 2) are derived and not entered |
| `fan2014a` | Fan 2014 (730 nm, PRA): R(26) 5-13 a15 absolute at 0.1 MHz | 1 | 1 of 2; the R(137) 5-12 a19–21 reference is a blend |
| `kobayashi2015a` | Kobayashi 2015 (531 nm): R(36) 32-0 a1 absolute at 8 kHz, at 41 Pa and 12.7 mW, not extrapolated to zero pressure | 1 | complete |
| `hauden2024a` | Hauden 2025 (556 nm): all 21 components of P(49) 24-1 at 5–6 kHz, corrected to zero pressure and asymmetry by the authors | 21 | complete; Table 3 hyperfine constants not entered |
| `grieser1994a` | Grieser 1994 (549/585 nm): absolute R(99) 15-1 a13 and R(85) 25-0 a1 at 67–71 kHz via wavelength ratios to the 633 nm He-Ne standard, plus P(13) 43-0 a3 computed from the printed λ514/λ633 ratio | 3 | complete; the 549 nm line is printed as 26-0 and re-assigned to 25-0 (as in Knöckel 2004); CIPM 1992 He-Ne basis, not rescaled (+8 kHz) |
| `badr2006a` | Badr 2006 (661 nm, Ag two-photon paper): P(62) 4-5 a1 absolute from the BIPM comb calibration of the reference laser, 15 kHz (authors' enlarged value; 6 MHz FM lock at a 13 °C cold point, not corrected to zero pressure) | 1 | complete; the only v″ = 5 precision line |
| `hong2001b` | Hong 2001b (532 nm): hyperfine intervals of R(56) 32-0 (a2, a5–a15) to a1 at 0.3 kHz (Table 1, main lines) | 12 | complete for main lines; 23 crossovers not representable; P(54) 32-0 main lines (Table 2) are a reprint of `ye1999a` and are left out; probably **a source of `bipm2012a`** Table 16 (agrees to 0.56 kHz rms) |
| `sakagami2020a` | Sakagami 2020 (531.5 nm): a1/b1 absolute frequencies of 7 lines, bands 32-0 to 35-0 (J″ = 34–112), at 5.7 kHz (2.4 Pa, not extrapolated) and 114 splittings at 1 kHz | 121 | complete; R(75) 33-0 a18 and R(37) 32-0 b2 not measured; a11–a14 of R(38)/P(34) 32-0 kept at 10 kHz; abstract's "P(98)34-0" is R(98)34-0 |
| `yoshii2019a` | Yoshii 2019/2020 (531.5 nm): R(38) 32-0 a10 absolute at 6 kHz (2.5 Pa, 3.0 mW, not extrapolated) | 1 | complete; **not independent of `sakagami2020a`** (same apparatus and budget), agrees within 1 kHz |
| `tanabe2022a` | Tanabe 2022 (556 nm): a1 absolute frequencies of R(53) 24-1, P(49) 24-1, R(95) 25-1 at 7 kHz (4.0 Pa, 2 mW, not extrapolated) and 60 splittings at 2 kHz | 63 | complete; Table 4 hyperfine constants not entered; independent of `hauden2024a`, which it matches at a1 to 3 kHz and on average to −8 kHz |
| `sharma2023a` | Sharma 2023 (739 nm): wavelength-meter frequencies of R(78) 1-11 a1, a10, a15 at 60 MHz | 3 | **excluded**: inconsistent with the known hyperfine gaps and with the paper's own Fig. 3; 3 blended peaks left out |
| `yoshii2021a` | Yoshii 2021 (531 nm, Negative Results): R(36) 32-0 a1 absolute at 9.9 kHz and 20 P(35) 32-0 splittings at 10 kHz, from a ~1 MHz-linewidth laser | 21 | **excluded**: the authors show the splittings are compressed (up to −50 kHz); the a1 absolute is 239 kHz above `kobayashi2015a` |
| `arie1993a` | Arie & Byer 1993 (532 nm): a1 gaps of 7 lines to R(56) 32-0 a1 at 0.2 MHz (incl. P(103) 34-0, held by no other set); hyperfine intervals of P(53)/R(56) 32-0, P(83)/R(86) 33-0, R(106) 34-0, R(134) 36-0 at the per-line fit sd (2.3–9.8 kHz); P(119) 35-0 a21−a1 | 85 | 4 blended/shifted entries left out; errata (JOSA B 11, 866) not seen |
| `arie1994a` | Arie & Byer 1994 (532 nm): P(119) 35-0 a2–a21 relative to a1, 3.7 kHz (fit sd) | 20 | complete; erratum (Opt. Commun. 127, 382) not seen; **probably a source of `bipm2012a`** Table 11 |
| `cheng2001a` | Cheng & Shy 2001 (543 nm): R(12) 26-0 a9–a15 and R(106) 28-0 b1–b15 relative to b10, printed 0.5–1.8 kHz (statistical only) | 21 | complete; a12 misprint (−746.3475 → −476.3475) corrected |
| `zhang2001a` | Zhang 2001 (532 nm): gaps from R(56) 32-0 a10 to 16 lines (a1 of 15, a21 of P(83) 33-0), comb generator, 1 kHz, at −15 °C | 16 | complete; **a source of `bipm2012a`** (column [8]; the BIPM intervals of six lines are these values alone) |
| `hong2001a` | Hong 2001a (532 nm): hyperfine intervals to a1 of R(58) 32-0, P(55) 32-0 and P(104) 34-0, 1–2 kHz | 44 | complete; R(58) a3/a4 and P(55) a5/a6 not measured; **the source of the `bipm2012a` Tables 2–4** |
| `hong2002a` | Hong 2002a (532 nm): hyperfine intervals to a1 of R(87) 33-0, R(145) 37-0 and P(132) 36-0, 1–2 kHz | 54 | complete; **the source of the `bipm2012a` Tables 1, 7, 8**; BIPM R(145) a5 is a misprint (−50 kHz) |
| `hong2000a` | Hong & Ishikawa 2000 (532 nm): hyperfine intervals a2–a15 of R(122) 35-0 and P(84) 33-0 relative to a1, 0.85 and 0.12 kHz (from relock reproducibility; none printed) | 28 | complete; **source of `bipm2012a` Tables 5–6**, whose R(122) a7 (398.2113 MHz) is a misprint of 398.2213 |
| `hong2004a` | Hong 2004a (532 nm): R(85) 33-0 a1 absolute (comb, uncorrected, 0.52 kHz), its interval to R(56) 32-0 a10, and a2–a21 relative to a1 at 0.1 kHz | 22 | complete; Table 2 Δ constants not entered; **source of `bipm2012a` Table 20 and the CI-2007 interval** |
| `hong2004b` | Hong 2004b (532 nm): R(56) 32-0 a10 absolute of laser Y3 at −10 °C, uncorrected, 0.52 kHz | 1 | complete; laser-to-laser offsets and shift slopes not entered |
| `chen2023a` | Chen 2023/2024 (554 nm, arXiv version): one absolute component of each of 9 lines, v′ = 22–26, v″ = 0–1, at 33–47 kHz | 9 | **excluded** (`exclude` in `meta.toml`): all 9 disagree with the model and the NIST/APO atlas by 30–840 MHz, and no relabelling explains it |

Two sources were checked and hold nothing to transcribe: Chen & Ye 2003 (Chem. Phys. Lett. 381, 777) is a
figures-only summary of `chen2004a`, and Chen, de Jong & Ye 2005 (JOSA B 22, 951) is theory.

**Overlapping sources.** `bipm2012a` is built from primary papers that are also here: `ye1999a`,
`holzwarth2001a` and `zhang2001a` (the intervals to R(56) a10), `hong2001a` (its Tables 2–4), `hong2002a`
(Tables 1, 7, 8), `hong2000a` (Tables 5–6), `hong2004a` (Table 20 and the 2007 interval), and probably
`arie1994a` (Table 11) and `hong2001b` (Table 16). `bipm2005a`'s a3 is the mean of `jones2002a` and
`goncharov2004a`, and `yoshii2019a` shares its apparatus with `sakagami2020a`. A fit that uses both sides
counts those measurements twice; the group labels do not yet express this. Transcribing the sources found
two misprints in the MEP 532 nm tables, corrected in `bipm2012a` (R(122) 35-0 a7, R(145) 37-0 a5).

## matyugin2012: how the assignment was settled

Matyugin et al. 2012 (Quantum Electron. 42, 250) measures 20 emission frequencies at v'' = 48 to
3-65 kHz. Its table labels name lines whose branch (R or P into v'' = 48) the paper does not state, and
taken one way the implied model error swings by 15 cm⁻¹ across J = 85-88. The companion paper, Matyugin
et al., Quantum Electron. 38, 755 (2008), states the scheme: pump J'' = 56, v'' = 0 -> J' = 57, v' = 32,
emission J' = 57, v' = 32 -> J'' = 58, v'' = 48. That fixes the convention: the 2012 labels are emission
lines in standard notation, and pairs of rows share an upper level. 18 rows are kept, consistent with
this assignment (`data/observations/matyugin2012/meta.toml`).

## Not yet implemented

- A second component label per line: both sets record only the "t" component, so the hyperfine
  *splittings* of these lines are untouched by them.
- Kinds for intensities, cross sections and band-averaged absorption (the Salami & Ross and Spietz comparisons in `docs/research/spectra-validation.md`).
- Per-line hyperfine parameters (ΔeqQ, ΔC), e.g. the Salumbides 2006 supplement.
- Correlations between observations beyond the `group` label, e.g. a shared reference frequency.
- A MARVEL-style network check of hyperfine-free line centres (combination differences through shared levels).
- A declaration of each group's nuisance model (offset, scale) in `meta.toml`.

## Level constants: `data/x_levels/`

Sources that publish per-level constants rather than lines. `martin1986/` holds Table I of Martin *et al.*,
J. Mol. Spectrosc. 116, 71 (1986): G, B, D, H, L (and M) of 93 X levels, v″ = 8–108, from 14 820 B→X
fluorescence lines, relative to X(0,0) through Luc's G(9) (column definitions in its `meta.toml`). The model
agrees with it to 15–170 MHz at v″ = 26–47 (J ≤ 120), where nothing else measures the X state; above
v″ ≈ 60 the extended X potential of i2spec2026k is off by up to 70 cm⁻¹. Not yet used by any fit.

## Derived hyperfine parameters: `data/hyperfine_parameters/`

Some sources publish the hyperfine constants they fitted to each line rather than the component
frequencies. Those numbers depend on the source's own model choices, such as which X-state values were
held fixed and how many terms were fitted. So they are kept apart from `data/observations/`, with the same
layout (`<id>/meta.toml` and `<id>/data.csv`) and the same required `meta.toml` keys. `meta.toml` must
also record those model choices; for a B-state set, that means the X-state values held fixed.

`data.csv` columns, in this order: `line,eqQ,eqQ_unc,C,C_unc,d,d_unc,delta,delta_unc,fit_sd_kHz,note`,
all in MHz except `fit_sd_kHz`. Uncertainties are 1σ.

```python
from i2spec.observations import load_all_hyperfine_parameters
for s in load_all_hyperfine_parameters():
    for r in s.rows:              # MeasuredHyperfine: r.line, r.J_upper, r.eqQ, r.C_unc, ...
        ...
```

| id | Content | Rows | Status |
|---|---|---|---|
| `chen2004a` | Chen 2004 (500–517 nm): eqQ_B, C_B, d_B, δ_B of 74 lines, v′ = 42–70 with 1–7 J′ per v′ (J′ ≈ 9–112); X state held at BKT02, as in i2spec | 74 | complete for Table 1; the authors left out five lines perturbed by a 1g state |
| `hong2001b` | Hong 2001b (532 nm): eqQ′, C′, d′, δ′ of R(56) and P(54) 32-0 from main + crossover fits; X eqQ″, C″ fitted (in `note`), d″/δ″ held at 1.524/3.705 kHz | 2 | complete for Table 3 |
