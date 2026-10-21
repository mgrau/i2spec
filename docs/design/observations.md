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
| `[overlaps]` | no | In a compilation: the sets whose measurements it prints, and how (see "One record per measurement") |
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

## One record per measurement

Every measurement a paper reports is kept, including those that a later paper or a compilation prints again.
A measurement reported more than once is one record with several citations: the paper that reported it
first, then the compilations (the CIPM's MEP documents) that print it. Each set is still transcribed whole
from its own paper, so the repetitions stay on disk and check each other; the compilation declares them.

A compilation's `[overlaps]` table names each set it draws on, with a relation:

```toml
[overlaps]
hong2000a = "same"         # its rows are this paper's measurements printed again, perhaps rounded or re-referenced
arie1994a = "mean"         # its rows are its own combination of this paper's measurement with others
reinhardt2006a = "superseded"   # its rows are an older measurement of what this paper remeasured
borde1981a = { relation = "mean", lines = ["127I2 P(13) 43-0", "127I2 R(15) 43-0"] }   # only these lines
hong2001b = { relation = "same", groups = ["picard2003a_jilanist2"] }   # only the rows of this group
```

`groups` is for a paper that prints other papers' measurements beside its own on the same quantities, such as
`picard2003a`, whose Table 5 sets three new lasers next to columns reprinted from `hong2001b` and `arie1993a`:
each reprinted column is its own group, and only that group's rows print the source.

A source "determines" a compilation row when its own rows (absolute frequencies or intervals) fix the same
quantity, so a table printed relative to a10 matches a paper's table relative to a1. Then:

- **same**: the row is no record of its own. The source's rows get the compilation as a further citation.
- **mean**: the row stays a record of its own, the compilation's value. The source's rows get the compilation
  as a further citation, since it prints or uses them. A row that two or more sources determine is a mean
  whatever the entries say: one value cannot be two measurements. A source that supersedes the row does not
  count here, since it does not print the row (2026-10-21: `bipm2003c` prints Barwood's P(62) 17-1 table, which
  `razet1993b` holds to the kHz and `reinhardt2006a` supersedes, so the row is `razet1993b`'s measurement
  printed again, not a mean).
- **superseded**: the row stays a record of its own. The source's rows get no citation.

A source that measured the same quantity more than once (two chains, two lasers) has every such row cited.
A later paper of the same group can also be the source of a "superseded" entry in the earlier set:
`sansonetti1990a` names `sansonetti1997a` and `ikeda2020a` names `ikeda2022a`, so the earlier values stay in the
record and the fit takes the later ones.

Every observation carries `citations`, a tuple of set ids with its own set first. `load_all()` gives the
fit: each measurement enters once, from its source, so it leaves out every compilation row a source
determines. `load_all(record=True)` gives the record: it leaves out only the "same" rows, and keeps the
compilations' own values, such as the CIPM recommended frequencies, which the explorer shows under the
compilation. A set with `exclude` (a pending set) has no say until it is admitted. The tests check that every
"same" row equals its source's value to within 1 kHz or half its uncertainty, and that the fit uses only
rows of the record.

The table replaced `defer_to` (a list of sources whose rows the compilation gave way to) on 2026-10-20. The
fit is unchanged by the change.

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
| `arie1993a` | Arie & Byer 1993 (532 nm): a1 gaps of 7 lines to R(56) 32-0 a1 at 0.2 MHz (incl. P(103) 34-0, held by no other set); hyperfine intervals of P(53)/R(56) 32-0, P(83)/R(86) 33-0, R(106) 34-0, R(134) 36-0 at the per-line fit sd (2.3–9.8 kHz); P(119) 35-0 a21−a1 | 85 | 4 blended/shifted entries left out; errata (JOSA B 11, 866) read 2026-10-20: they correct only Table 4's footnote and Table 9, whose standard deviations are the ones used, so no row changed |
| `arie1994a` | Arie & Byer 1994 (532 nm): P(119) 35-0 a2–a21 relative to a1, 3.7 kHz (fit sd) | 20 | complete; erratum (Opt. Commun. 127, 382) not seen; **probably a source of `bipm2012a`** Table 11 |
| `cheng2001a` | Cheng & Shy 2001 (543 nm): R(12) 26-0 a9–a15 and R(106) 28-0 b1–b15 relative to b10, printed 0.5–1.8 kHz (statistical only) | 21 | complete; a12 misprint (−746.3475 → −476.3475) corrected |
| `zhang2001a` | Zhang 2001 (532 nm): gaps from R(56) 32-0 a10 to 16 lines (a1 of 15, a21 of P(83) 33-0), comb generator, 1 kHz, at −15 °C | 16 | complete; **a source of `bipm2012a`** (column [8]; the BIPM intervals of six lines are these values alone) |
| `hong2001a` | Hong 2001a (532 nm): hyperfine intervals to a1 of R(58) 32-0, P(55) 32-0 and P(104) 34-0, 1–2 kHz | 44 | complete; R(58) a3/a4 and P(55) a5/a6 not measured; **the source of the `bipm2012a` Tables 2–4** |
| `hong2002a` | Hong 2002a (532 nm): hyperfine intervals to a1 of R(87) 33-0, R(145) 37-0 and P(132) 36-0, 1–2 kHz | 54 | complete; **the source of the `bipm2012a` Tables 1, 7, 8**; BIPM R(145) a5 is a misprint (−50 kHz) |
| `hong2000a` | Hong & Ishikawa 2000 (532 nm): hyperfine intervals a2–a15 of R(122) 35-0 and P(84) 33-0 relative to a1, 0.85 and 0.12 kHz (from relock reproducibility; none printed) | 28 | complete; **source of `bipm2012a` Tables 5–6**, whose R(122) a7 (398.2113 MHz) is a misprint of 398.2213 |
| `hong2004a` | Hong 2004a (532 nm): R(85) 33-0 a1 absolute (comb, uncorrected, 0.52 kHz), its interval to R(56) 32-0 a10, and a2–a21 relative to a1 at 0.1 kHz | 22 | complete; Table 2 Δ constants not entered; **source of `bipm2012a` Table 20 and the CI-2007 interval** |
| `hong2004b` | Hong 2004b (532 nm): R(56) 32-0 a10 absolute of laser Y3 at −10 °C, uncorrected, 0.52 kHz | 1 | complete; laser-to-laser offsets and shift slopes not entered |
| `chen2023a` | Chen 2023/2024 (554 nm, arXiv version): one absolute component of each of 9 lines, v′ = 22–26, v″ = 0–1, at 33–47 kHz | 9 | **excluded** (`exclude` in `meta.toml`): all 9 disagree with the model and the NIST/APO atlas by 30–840 MHz, and no relabelling explains it |
| `powel2021a` | Powel 2021 (704 nm): comb-referenced zero-pressure absolutes of R(118) 2-8 a10, R(54) 3-9 a1, P(48) 3-9 a15 at 0.36–0.42 MHz | 3 | **pending** (`exclude`, see below); wavemeter column and Table 3 fit not entered |
| `balling2007a_633` | Balling & Křen 2008 (633 nm): R(127) 11-5 a12–a18 (MeP d–j) of two He-Ne/I₂ lasers, SMU-1 and PLO3, at 0.25–0.4 kHz (printed k = 2 halved) | 14 | **pending**; one paper split into three sets, one `[shift_correction]` each |
| `balling2007a_543` | Balling & Křen 2008 (543 nm): 22 components of R(106) 28-0 and R(12) 26-0, laser PLG2, at 1–9 kHz | 22 | **pending**; the printed CCL value of R(12) 26-0 a15 is 2 MHz off (858 for 856) and is corrected |
| `balling2007a_532` | Balling & Křen 2008 (532 nm): R(56) 32-0 a1, a2, a5–a15, laser YAG1, at 0.5 kHz, 0.85 Pa | 13 | **pending** |
| `borde1981a` | Bordé 1981 (515/502 nm): hyperfine intervals of P(13) and R(15) 43-0 (1 kHz), R(98) 58-1 (5 kHz) and R(26) 62-0 (10 kHz), read from the scanned tables | 60 | **pending**; one of the sources averaged into `bipm2005a` Tables 1–2 (declared in its `[overlaps]`); its R(98) 58-1 is independent of BIPM's |
| `demelo2024a` | de Melo 2024 (556 nm): R(158) 25-0 a10 relative to P(49) 24-1 a21, 6 MHz | 1 | **pending**; J″ = 158, the highest J held; +18.7 MHz from the model, which the two FTS atlases also put 14 MHz low |
| `diddams2000a` | Diddams 2000: R(56) 32-0 a10 at 5 kHz; R(127) 11-5 a13 of He-Ne SN-126 at 2.1 kHz | 2 | **pending**; a source of the CIPM 532 nm value |
| `yoon2001a` | Yoon 2001 (633 nm): R(127) 11-5 a13 of JILA-145 at 1.06 kHz (CIPM conditions) and f − i | 2 | **pending**; behind the MeP 633 nm value |
| `weller2026a` | Weller 2026 (508 nm): all 20 intervals a_k − a1 of R(41) 51-0 at 0.5 kHz | 20 | **pending**; the first hyperfine data at B v′ = 51 |
| `sansonetti1990a` | Sansonetti 1990, LLNL report UCRL-CR-107777: one component of each of 52 lines, 560–656 nm, at 0.55 MHz | 52 | **pending**; superseded by `sansonetti1997a` (same 52 lines re-evaluated, −0.96 to +0.21 MHz) |
| `ye2000a` | Ye 2000 (633 nm): f(a16, R(127) 11-5) of JILA-145 by the JILA comb (0.34 kHz) and the NRC chain (0.69 kHz) | 2 | **pending**; a source of the MeP 633 nm f(a16), declared in `bipm2003a` |
| `ye2001a` | Ye 2001 (532 nm): R(56) 32-0 a10 from a year of comb measurements, SD 0.24 kHz | 1 | **pending**; partly the same data as `diddams2000a` (not declared, see `meta.toml`) |
| `ma2006a` | Ma 2006 (543 nm): ten comb measurements of eight He-Ne/I₂ lasers on R(106) 28-0 b10, b15 and R(12) 26-0 a9, a11, a14, a15, as measured, at 5–9 kHz (120 kHz for one) | 10 | **pending**; the weighted means are not entered; DFM1 (outlier) and 08A01/G10 (cell at 21 °C) flagged |
| `merimaa2008a` | Merimaa 2008 (543 nm): R(106) 28-0 b2 of laser MGI2 at 4 MHz modulation, 5 kHz | 1 | **pending** |
| `hackel1980a` | Hackel 1980 (515 nm): R(15) 43-0, 13 components from b13, at 2.6 kHz; P(13) 43-0, 15 components from a7 (the 1975 data), at 5.4 kHz | 28 | **pending**; the P(13) 43-0 rows fall under the standing exclusion of that line from all fits |
| `picard2003a` | Picard 2003 (532 nm): Table 5, R(56) 32-0 a1–a15 (not a3, a4) from a10 for three new lasers, the BIPM laser and four earlier works, at the per-column fit SD (0.4–6.6 kHz) | 83 | **pending**; the JILA–NIST 2 and A&B columns reprint `hong2001b` and `arie1993a` (declared with `groups`) |
| `eickhoff1995a` | Eickhoff & Hall 1995 (532 nm): R(106) 34-0 a2–a15 from a1 at 0.52 kHz, and a1 − P(83) 33-0 a20 | 15 | **pending**; `bipm2012a` Table 13 is this table rounded (declared "same") |
| `jungner1995a` | Jungner 1995 (532 nm): R(56) 32-0 a10 via the Rb D2 line, 70 kHz | 1 | **pending** |
| `nishiyama2013a` | Nishiyama 2013 (576–612 nm): comb-marked absolutes of P(62) v′-1 a1 (v′ = 12–19) and R(47) v′-2 a7 (v′ = 9–15), 0.03–0.12 MHz | 15 | **pending**; ±0.7 MHz from the shipped model at these v′ |
| `ikeda2020a` | Ikeda 2020 (514 nm): R(73) 46-0 a1 at 2.5 Pa, 6 kHz | 1 | **pending**; superseded by `ikeda2022a` (68 kHz apart; the later paper suspects a servo offset) |
| `smith2007a` | Smith 2007 (633 nm): R(127) 11-5 a13 of a Winters He-Ne/I₂ laser, 0.66 kHz, corrected to CIPM conditions | 1 | **pending** |
| `fang2007a` | Fang 2007 (532 nm, in Chinese): R(56) 32-0 a10 of the NIM laser, one 32 min comb run, 25 Hz (measurement only) | 1 | **pending**; no iodine budget; a fit must add the laser's realization uncertainty |
| `cheng2021a` | Cheng 2021 (532 nm): R(56) 32-0 a10, corrected for RAM, impurity and comb, 5.3 kHz | 1 | **pending**; own pressure and power coefficients applied |
| `macfarlane1999a` | Macfarlane 1999 (532 nm, NPL): P(54) 32-0 a10 and R(57) 32-0 a21 absolutes at 60–68 kHz; separations from a1 of both lines | 30 | **pending**; R(57) a8–a10 relabelled a9–a11; R(57) a14 (a misprint) and two blends left out; separations off by tens to hundreds of kHz from `ye1999a` |
| `schweitzer1973a` | Schweitzer 1973 (633 nm, NBS): ¹²⁹I₂ P(54) 8-4 a28 against ⁸⁶Kr (2 MHz) and beats to ¹²⁷I₂ R(127) 11-5 a13 (9.8 kHz) and ¹²⁹I₂ P(33) 6-3 e2 (83 kHz) | 3 | **pending**; the 1973 letters k and B assigned here from the MEP |
| `cerez1974a` | Cérez 1974 (633 nm): R(127) 11-5 a12–a21 from a13 (peaks a–j), 0.5 kHz (statistical) | 9 | **pending**; 7–30 kHz from the MEP |
| `dschao1980a` | Dschao 1980 (612 nm): ¹²⁷I₂ R(47) 9-2 a1–a20 and 20 ¹²⁹I₂ P(110) 10-2 components, from R(47) a7, at 0.05 MHz | 39 | **pending**; three unresolved ¹²⁹I₂ doublets left out; a source of `bipm2003d` (mean) |
| `zhao1985a` | Zhao 1985 (640 nm): P(10) 8-5 and R(16) 8-5 components from P(10) a9, two lasers, 8–220 kHz | 20 | **pending**; a source of `bipm2003e` (mean); a1 (estimated) and the a4/a5 blend left out |
| `bertinetto1987a` | Bertinetto 1987 (605 nm): R(161) 9-0 b5–b21 from P(62) 11-2 a15, 0.1 MHz | 17 | **pending**; J″ = 161; P(62) blends and unassigned peaks left out |
| `razet1991a` | Razet 1991 (612 nm): R(47) 9-2, P(48) 11-3, R(48) 15-5 intervals at 3–15 kHz, the links between them, and R(47) 9-2 a13 absolute (their ref. 7) | 36 | **pending**; column headings print the reversed sign; a source of `bipm2003d` (mean) |
| `picard1996a` | Picard 1996 (637 nm): all 21 components of R(39) 7-4 from a4, 0.03–0.09 MHz, measured by Edwards et al. 1996 (not held) | 20 | **pending**; 8 partially resolved components starred in the notes |
| `pollock1983a` | Pollock 1983 (576 nm, NBS chain): P(62) 17-1 a1 (o), 83 kHz | 1 | **pending**; the "NBS value" behind `bipm2003c`'s a1 (mean) |
| `baird1979a` | Baird 1979 (576 nm): P(62) 17-1 a1 at 60 MHz and a2–a10 from a1 at 1.4 MHz | 10 | **pending**; the first absolute frequency in the visible |
| `sugiyama1988a` | Sugiyama 1988 (612 nm): P(48) 11-3 b1–b10 from b15 at 0.4–1 MHz and b15 − R(47) 9-2 a7 | 11 | **pending**; two unresolved pairs left out |
| `shiner1993a` | Shiner 1993 (597–612 nm): ten single-component absolutes at 0.25 MHz (four in MHz, six in cm⁻¹) | 10 | **pending**; R(9) 11-2 "rs" (a blend) left out |
| `zhang1999a` | Zhang 1999 (633 nm, in Chinese): R(39) 6-3 a1 − R(127) 11-5 a21 and all 21 R(39) 6-3 components at 0.5 MHz; four Fabry–Perot absolutes at 30 MHz | 25 | **pending**; the same measurements as `chen1998a` (declared "same"); two of the 30 MHz values are 57 and 113 MHz from the model |
| `ruben1973a` | Ruben 1973 (515 nm, molecular beam): all 21 components of P(13) and R(15) 43-0 at 1.3–1.4 MHz, and the common zero | 41 | **pending**; Fabry–Perot scale 0.55 % large; P(13) rows under the standing exclusion |
| `landsberg1976a` | Landsberg 1976 (refit paper): Hanes et al. 1973's eleven R(127) 11-5 d–j separations (20 kHz) and Hackel et al. 1975's P(13) 43-0 table | 26 | **pending**; the only copy of the Hanes 1973 data; Table 2 declared "same" as `hackel1980a` |
| `foth1979a` | Foth & Spieweck 1979 (515 nm): 12 R(98) 58-1 splittings from a6 at 30 kHz and a15 − P(13) 43-0 a1 = −1458 MHz | 13 | **pending**; OCR'd scan repaired with the printed obs − calc column |
| `rong1993a` | Rong 1993 (548.6 nm): 16 beat-frequency spacings chaining 17 lines of bands 24-0 to 28-1 over 150 GHz, 0.1 MHz | 16 | **pending**; assignments made here (atlas numbers only in the paper); one row a 100 MHz misprint, corrected from the printed sums |
| `kremser1994a` | Kremser 1994 (793 nm): R(92) 0-15 a2–a10, a15 from a1 at 12–42 kHz | 10 | **pending**; overlap corrections applied (sign settled on the a3/a4 pair) |
| `razet1993a` | Razet 1993 (633 nm): all 21 components of P(33) 6-3 from b21 at 10 kHz; b21 − R(127) 11-5 d; d − f | 22 | **pending**; MEP ref. [30], declared in `bipm2003a` |
| `razet1993b` | Razet 1993 (576 nm): Barwood & Rowley 1984's P(62) 17-1 intervals to the kHz | 10 | **pending**; superseded by `reinhardt2006a`; `bipm2003c`'s table declared "same" |
| `balling1999a` | Balling & Blabla 1999 (612 nm): R(47) 9-2 a1–a21 from a7 for five laser/cell/harmonic set-ups (30 kHz), R(48) 15-5, P(48) 11-3, R(34) 17-6 by 2f spectra (0.1 MHz), h−i and i−j checks | 118 | **pending**; crossovers not representable; R(34) 17-6 is 0.87 MHz from the model |
| `chen1998a` | Chen 1998 (633 nm, SPIE): the measurements later printed as `zhang1999a` | 25 | **pending**; the earlier publication, so it holds the record |
| `chen2000a` | Chen & Wang 2000 (634–635 nm): hyperfine patterns of P(94), R(80), P(74) 8-4 at MHz level | 23 | **record only** (excluded): scale fixed to the calculated a − o, 2–21 MHz off |
| `rovera2002a` | Rovera 2002 (532 nm): BNM-INM R(56) 32-0 a10 by comb, 1.0 kHz | 1 | **pending** |
| `magyar1980a` | Magyar & Brown 1980 (633 nm): 46 peaks of ¹²⁹I₂ P(54) 8-4, P(69) 12-6 and ¹²⁷I¹²⁹I P(33) 6-3 from k, 30 kHz | 46 | **pending**; assignments from the MEP letters; two misprints (m″, e″) corrected; MEP ref. [41] |
| `knox1971a` | Knox & Pao 1971 (633 nm): 30 adjacent-peak separations of ¹²⁹I₂/¹²⁷I¹²⁹I, 0.5–1.5 MHz | 30 | **pending**; letters mapped through the MEP; MEP ref. [39] |
| `tesic1975a` | Tesic & Pao 1975: Knox & Pao's table reprinted | 30 | **pending**; declared "same" as `knox1971a` |
| `jungner1995b` | Jungner 1995 (532 nm, SPIE): R(56) 32-0 hyperfine table from a10 (1.3 kHz), P(83) 33-0 a20 offset, a10 absolute 40 kHz | 13 | **pending**; supersedes `jungner1995a`'s a10; `picard2003a`'s JILA-NIST 1 column declared "same" |
| `lea2003a` | Lea 2003 (633 nm): f of NPL-R18 and NPL-S by comb, 0.6–0.7 kHz | 2 | **pending**; `[shift_correction]` as `bipm2003a` |
| `acef1993a` | Acef 1993 (633 nm): INM12 on f by the LPTF harmonic chain, 5 kHz | 1 | **pending**; 17 kHz below the CIPM 2001 value |
| `bernard2001a` | Bernard 2001 (633 nm): NRC INMS-3 on f via the Sr⁺ standard, 0.72 kHz at MeP conditions | 1 | **pending**; the NRC chain of the CIPM 2001 value (declared in `bipm2003a`) |
| `madej2004a` | Madej 2004 (633 nm): INMS 3 f (0.57 kHz) and INMS 2 e, f in three sessions at NRC and BIPM | 7 | **pending**; INMS 3 mean includes `bernard2001a` (declared "mean") |
| `hamid2006a` | Hamid 2006 (633 nm): UME-L3 on f, BIPM and UME combs, 1.1 kHz | 2 | **pending** |
| `samoudi2012a` | Samoudi 2012 (633 nm): CEM2 on d–i (CEM comb) and f (BEV comb) | 7 | **pending**; printed k = 2 halved |
| `morinaga1990a` | Morinaga 1990 (657 nm): R(69) 3-4 a15 from the Ca intercombination line, 20 kHz | 1 | **pending**; "R(67)" in the result is a misprint; Ca value added here |
| `gill1979a` | Gill & Bennett 1979 (515–568 nm): four wavemeter line centres at 1e-7 | 4 | **record only** (excluded): atlas precision |
| `rakowsky1989a` | Rakowsky 1989 (720–770 nm): the lowest component of 17 lines, v′ = 0–2, v″ = 10–13, at 0.0011 cm⁻¹ | 17 | **pending**; unit cm⁻¹ |
| `wallerand2006a` | Wallerand 2006 (514.7 nm): absolutes of P(49) 44-0 a1, a3–a21 and P(13) 43-0 a3 at 1.5 kHz | 21 | **pending**; the printed a1/a2 rows conflict with Fig. 2 (a1 relabelled, a2 left out); P(13) under the standing exclusion |
| `juncar1981a` | Juncar 1981 (589 nm): vacuum wavelengths of P(114) 14-1 and P(38) 15-2 components, 0.17–0.39 MHz | 2 | **pending**; ranks a5, a6 assigned here (the printed letters do not fit); rescaled from the 1973 He-Ne basis |
| `zhang2024a` | Zhang 2024 (532 nm, PGMF): R(56) 32-0 a10, 5.1 kHz | 1 | **pending**; same group and systematics as `cheng2021a` |
| `sorem1971a` | Sorem 1971 (501.7 nm): R(26) 62-0, 14 components from h at 2 % | 14 | **pending** |
| `ciddor1980a` | Ciddor & Brown 1980 (612 nm): all 28 components of ¹²⁹I₂ P(110) 10-2 (30 kHz) and recorder values for ¹²⁹I₂ R(113) 14-4 and ¹²⁷I₂ P(48) 11-3 | 51 | **pending**; MEP ref. [17] |
| `glaser1981a` | Gläser, Dschao & Foth 1981 (612 nm): eight P(110) 10-2 components Dschao left as pairs, 16 of R(113) 14-4, and Dschao's table again | 44 | **pending**; MEP ref. [18]; Table 2 declared "same" as `dschao1980a` |
| `fredinpicard1990a` | Fredin-Picard & Razet 1990 (543 nm): Chartier et al. 1989's R(12) 26-0 and R(106) 28-0 intervals | 13 | **pending**; the only copy of Chartier 1989; MEP ref. [7] |
| `lin1994a` | Lin 1994 (543 nm): R(12) 26-0 and R(106) 28-0 from b8, with Brand 1993's and Gläser's PTB-report columns | 56 | **pending**; MEP ref. [8] |
| `blabla1993a` | Blabla 1993 (633 nm): CMI R(127) 11-5 a8–a21 from a13, third- and fifth-harmonic, corrected for overlaps | 26 | **pending** |
| `nevsky2001a` | Nevsky 2001 (532 nm): R(56) a10 and P(54) 32-0 a1, a10, a15 absolutes (1.2 kHz) and the three a10 gaps (2 kHz) | 7 | **pending** |
| `glaser1987a` | Gläser 1987 (612 nm): R(47) 9-2 a1–a15 and P(48) 11-3 b11–b15 from a7, external cell | 19 | **pending**; MEP ref. [12]; 3.6 kHz rms from the model |
| `thang2012a` | Thang et al. 2012 (633 nm): the VMI He-Ne/I₂ laser on f by the KRISS comb at three output powers, 1.1 kHz | 3 | **pending**; `[shift_correction]` as `bipm2003a` |
| `ikeda2018a` | Ikeda et al. 2018 (532 nm, arXiv version): R(86) 33-0 a10 of the YNU Nd:YAG laser by comb (1.0 kHz) and of a compact Nd:YVO₄ laser by beat (3.2 kHz) | 2 | **pending**; the printed Nd:YAG value has "322" for "302" (20 MHz), corrected |
| `klein2002a` | Klein & Arie 2002 (778 nm): nine Doppler-limited line centers of the 0-14 and 0-13 bands by wavemeter, 0.01 cm⁻¹ | 9 | **record only** (excluded): five "R" lines are P lines, corrected |

**Pending sets.** The 36 sets transcribed on 2026-10-20 (`powel2021a` to `zhang1999a`) and the 35 of
2026-10-21 (`ruben1973a` to the end, except the three record-only sets) carry an `exclude` line that starts with
"pending": 71 sets, 1278 rows. They load and validate, but `load_all()` leaves them out of every
fit and benchmark. Matt decided on 2026-10-21 that every new set stays out of training and serves as the blind
test of balanced-expanded (`research-20261020-balanced-expanded`). How a later model is trained on the complete
data is decided then. The 36 sets of 2026-10-20 were frozen and scored once in
`src/i2spec/data/research/20261020/validation/pending_blind_v1/`, and the explorer shows them as "blind test, not
fitted" (`load_all(include=...)`). Later sets stay out of the explorer until a registry, frozen from their identity
columns only, has scored them. Each transcription is checked against `Predictor()` (i2spec2026o) only, never
against balanced-expanded or the apps that now prefer it. Deleting the `exclude` line admits a set. The compilations already name the pending sets they print in
`[overlaps]` (`bipm2003a` `ye2000a`, `bernard2001a`, `landsberg1976a`, `cerez1974a`, `razet1993a`, `knox1971a`,
`magyar1980a`; `bipm2003b` `fredinpicard1990a`, `lin1994a`; `bipm2003c` `pollock1983a`, `razet1993b`; `bipm2003d`
`razet1991a`, `dschao1980a`, `glaser1987a`, `ciddor1980a`, `glaser1981a`; `bipm2003e` `zhao1985a`; `bipm2005a`
`borde1981a`; `bipm2012a` `diddams2000a`, `eickhoff1995a`). A pending set has
no say, so admitting one also moves the compilation rows it determines out of the fit (admitting `borde1981a`
moves the 32 rows of the BIPM 515 nm P(13) and R(15) 43-0 tables). Those rows stay in the record.

The three **record-only** sets (`chen2000a`, `gill1979a`, `klein2002a`) have an `exclude` that starts with "record only": they
are kept for the archive and are not candidates for the fits.

**Suspected misprint in `goncharov2007a`.** Its R(26) 62-0 a4 − a1 = 295 233.9 kHz (arXiv version, Table 1) is
89 kHz from `borde1981a` and 81 kHz from the model, where its other 13 intervals agree with both to 10–30 kHz.
Goncharov's own column "frequencies deduced from [11] − ours" ([11] = Bordé 1981) gives −4.4 kHz at a4, which
requires a4 − a1 = 295 333.9 kHz: "2" printed for "3". Not corrected; check the published Metrologia table.

Sources checked that hold nothing to transcribe: Chen & Ye 2003 (Chem. Phys. Lett. 381, 777) is a
figures-only summary of `chen2004a`; Chen, de Jong & Ye 2005 (JOSA B 22, 951) is theory; Khodakovskiy et al.
2009/2010 (arXiv:0912.3252, Proc. SPIE 7994) gives one 30.8 MHz splitting of an unassigned ¹²⁷I¹²⁹I line;
Guo et al. 2004 (Opt. Lett. 29, 1733) and Nomura et al. 2019 (JOSA B 36, 631) report stabilities, not
frequencies.

Checked on 2026-10-20 with nothing entered:
- No iodine numbers: the CPEM digests Hong 2002 and 2008, Cheng 2000, Hu 1998, Jungner 1994 (results promised for
  the talk) and Schnatz 2000; Chepurov 2007; Barbarat 2019 (Zeeman sensitivity only); Ko 2012 (an unassigned line
  65 MHz from the Ca-48 resonance); Madej 2003 (P(33) 6-3 offsets only plotted; its other results are He-Ne lasers
  referenced to iodine); Riis 1994 (prints only the Li⁺ frequencies); Piqué 1984 (theory matrix elements);
  Hänsch 1971 (P(117) 21-1 components only in a figure).
- Unassigned lines: Hodges 2025 (seven 689 nm peaks at 30 MHz, hot-band guesses only) and Xie 2019 (twenty
  650 nm lines relative to the Ba⁺ resonance, MHz-level wavemeter). Picqué 2001 gives only nominal offsets (the
  R(34) 27-0 a1–a15 spacing "440 MHz" and the He 1083 nm offset subtracted in its Fig. 5).
- Below the atlas: Chui 2005 (R(55) 1-13 by wavemeter, 0.02 cm⁻¹).
- Hyperfine constants only: Wallerand 2000 (B v′ = 6, J′ = 32), Wakasugi 1988 (B v′ = 13, J′ = 0–10) and
  Levenson 1972 (eleven lines). They belong in `data/hyperfine_parameters/`, whose loader has no `exclude`, so
  adding them now would change the fits that read it. Bacis 1980 (X-state eQq) is in `data/x_hyperfine/bacis1980/`,
  entered the same day by the research work on second-order hyperfine.

Checked on 2026-10-21 with nothing entered (the rest of the 106-paper download list):
- No iodine numbers, or figures only: Hackel 1976 (a reply to Landsberg: one octupole constant), Sorem & Schawlow
  1972 and Döpel 1978 (techniques, linewidths), Bennett & Cérez 1978 (spectra only), Churassy 1979 (12-0 band head
  in a supersonic jet: figures and fitted constants), the Tanabe 2022 CLEO-PR digest (its numbers are
  `tanabe2022a`'s), Matyugin 2008 (its table reprints R(56) 32-0 intervals from another paper; the emission
  hyperfine is in figures), Robertsson 2001 (laser-to-laser comparisons only), Fang 2005 (widths and broadening
  coefficients), Razet & Picard 1997 and its 1999 erratum (empirical formulas fitted to other papers' constants).
- Theory or compilation of others' measurements: Gläser 1985 (calculated 640 nm patterns; its "meas." column
  combines `zhao1985a`'s two lasers with another paper), Piqué 1980 (constants of three v′ = 15 lines of ¹²⁷I₂ and
  ¹²⁹I₂, spectra in figures).
- Not representable: Cheng 1998 (543 nm crossover resonances and unidentified peaks); Yokozeki & Muenter 1980
  (radio-frequency transitions between hyperfine sublevels of X(v=0, J=13) and B(v=43, J=12) in a molecular beam;
  its fitted X constants are the standard values the whole field holds fixed); Wu 1985 (a ¹²⁷I¹²⁹I hyperfine
  pattern at 612 nm that the paper leaves among six candidate transitions).
- Needs a reference that is not independent: Watanabe 2003 (Na D2 F=2 → F′=3 to P(38) 15-2 f, 1047.2(3) MHz;
  the Na line's absolute value is `juncar1981a`'s, which rests on this same iodine line).
- Below the atlas: Suwaiyan 1992 (Doppler-limited absorption peaks at 16 999–17 001 cm⁻¹, ±0.0025 cm⁻¹, many
  blends) and Cheng 2005 (wavemeter wavelengths to 0.001 nm near 767 nm).
- Hyperfine constants only: Knöckel 1996 (Delta eqQ and Delta C of 13 lines of the 0-15 band, J′ = 0–100) and
  Churassy 1979, Piqué 1980 above; with Wallerand 2000, Wakasugi 1988 and Levenson 1972 they wait for an `exclude`
  in `data/hyperfine_parameters/` and for optional d and delta columns (most give only Delta eqQ and Delta C).

Checked on 2026-10-21 in the literature search (OpenAlex, Crossref, arXiv and the WANTED list) with nothing entered:
- Wu et al. 2022 (JOSA B 39, 1457; arXiv:2203.00883): one hyperfine component of P(70) 1-11 at 739 nm, which one
  not said, from a HighFinesse WS7 wavemeter (about 30 MHz); +38 MHz from the model's line center.
- Döringshoff 2018 (thesis, HU Berlin): Table 6.1's "measured" R(56) 32-0 intervals equal the MEP 532 nm Table 16
  (`bipm2012a`, taken to a1) to the kHz at all eleven components both give, which an independent beat measurement
  would not. As printed it is a reprint, not a measurement.
- No line positions: Gillot et al. 2024 (Appl. Phys. B 130, 100; Zeeman slopes of four lines), Goïcoechéa et al. 2023
  (Opt. Express 31, 15316; microcell linewidths), Philippe 2017 (thesis, PSL; stability on P(34) 44-0 a1), Şahin
  2006 (MSc thesis, ITU; repeats `hamid2006a`'s two values and gives shift coefficients), Kovalchuk et al. 2005
  (Opt. Lett. 30, 3141; Nd:YAG/I₂ against CH₄, no iodine value), Wang et al. 2025 (arXiv:2510.20057) and Wu et al.
  2026 (arXiv:2603.00424) (stabilities on R(56) 32-0 a1), Rerucha et al. (arXiv:1905.00795; a 633 nm diode laser).
- A spectrum, not a line list: Mashburn et al. 2026 (arXiv:2602.07239), a comb-referenced dual-comb I₂ spectrum over
  500-600 nm at 100 MHz resolution, compared with the Salami-Ross atlas; the data would have to be requested.
- Already in: the Sansonetti 1996 LLNL report (UCRL-CR-125216) is the source `sansonetti1997a` was typed from.

The 2026-10-21 search found about 150 papers with iodine measurements that are not held; only those above were open
to a script. The rest (publisher bot walls: IOP, AIP, APS, Elsevier, Optica, IEEE, RSC, Wiley, Springer and
Canadian Science Publishing; some marked open access but walled to scripts) are listed for download in
`papers/inbox/DOWNLOADS_2026-10-21.md`, grouped by value.

Papers these sources print or cite that are not held (candidates to fetch): Hanes, Baird & DeRemigis 1973 (Appl.
Opt. 12, 1600; in `landsberg1976a`), Hackel et al. 1975 (PRL 35, 568; in `landsberg1976a`, `hackel1980a`), Barwood
& Rowley 1984 (Metrologia 20, 19; in `razet1993b`), Chartier, Fredin-Picard & Robertsson 1989 (Opt. Commun. 74, 87;
in `fredinpicard1990a`), Brand 1993 (Opt. Commun. 100, 361; in `lin1994a`), Gläser, PTB-Bericht PTB-Opt-25 (1987;
in `lin1994a`, a source of the MEP 543, 612 and 640 nm tables), Gläser 1986 (Metrologia 23, 45), Cérez & Bennett
1979 (Appl. Opt. 18, 1079), Hanes et al. 1971 (J. Mol. Spectrosc. 39, 506), Morinaga & Tanaka 1978 (Appl. Phys.
Lett. 32, 114), Bergquist & Daniel 1984 (Opt. Commun. 48, 327), Simonsen 1997 (IEEE TIM 46, 141), Brand, Mensing
& Helmcke 1989 (Appl. Phys. B 48, 343), Simonsen & Poulsen 1990 (Appl. Phys. B 50, 7), Bouchareine et al. 1991
(IEEE TIM 40, 200; the absolute behind `razet1991a`'s a13).

**Overlapping sources.** Declared in `[overlaps]`:
- `bipm2012a` (MEP 532 nm) reprints `hong2000a` (Tables 5–6), `hong2001a` (2–4), `hong2001b` (16),
  `hong2002a` (1, 7, 8), `hong2002b` (18–19), `hong2004a` (20 and the 2007 interval), `ye1999a` (9–10) and six
  of `zhang2001a`'s intervals to R(56) a10. That is 231 rows, all within 1 kHz of their papers.
- `bipm2012a`'s own values are its f(a10) (the mean of `holzwarth2001a`, `jones2002a`, `hong2004b` and the
  pending `diddams2000a`), ten intervals to a10 (means of `ye1999a`, `zhang2001a` and `holzwarth2001a`),
  R(56) a1 − a10, and Table 11 (a mean that includes `arie1994a`).
- `bipm2005a`'s a3 is the mean of `jones2002a` and `goncharov2004a`.
- `bipm2003c`'s P(62) 17-1 table (Barwood 1984) is superseded by `reinhardt2006a`.
- Pending sources, declared now: `bipm2012a` Table 13 reprints `eickhoff1995a` (rounded to 1 kHz);
  `bipm2003a`'s f(a16) is a mean including `ye2000a`; `bipm2003c`'s a1 includes `pollock1983a` (the "NBS value");
  the `bipm2003d` tables average `razet1991a` and `dschao1980a`, and the `bipm2003e` tables `zhao1985a`.
- `picard2003a` reprints `hong2001b` and `arie1993a` in two of its Table 5 columns (declared with `groups`).
- Earlier values a later paper replaces: `sansonetti1990a` (by `sansonetti1997a`), `ikeda2020a` (by `ikeda2022a`).
- Declared 2026-10-21, all pending: the MEP 633, 612 and 543 nm documents name their sources table by table, and
  each held source is declared as a "mean" in its compilation, restricted to the lines that table covers
  (`bipm2003a`: `landsberg1976a` and `cerez1974a` for Table 1, `razet1993a` for Table 2, `knox1971a` and
  `magyar1980a` for Tables 3, 4, 7, `bernard2001a` for f(a16); `bipm2003b`: `fredinpicard1990a`, `lin1994a`;
  `bipm2003d`: `glaser1987a`, `ciddor1980a`, `glaser1981a`). `bipm2003c`'s P(62) 17-1 table is Barwood 1984 as
  printed in `razet1993b` ("same"); `razet1993b` is superseded by `reinhardt2006a`.
- Reprints and re-evaluations between papers: `zhang1999a` prints `chen1998a`'s measurements ("same");
  `tesic1975a` prints `knox1971a` ("same"); `landsberg1976a` prints `hackel1980a`'s 1975 P(13) table ("same");
  `glaser1981a` Table 2 prints `dschao1980a` ("same", by group); `picard2003a`'s JILA-NIST 1 column is
  `jungner1995b`'s table ("same", by group); `madej2004a`'s INMS 3 mean includes `bernard2001a`'s run ("mean");
  `jungner1995a`'s a10 is re-evaluated in `jungner1995b` ("superseded").
- Found but not declared, because it would change the training data: the MEP 633 nm Table 2 (P(33) 6-3) cites
  Edwards et al. 1999 (ref. [34]), i.e. `edwards1999a`, which is in the fit. Declaring it as a "mean" source would
  take the 21 MEP P(33) 6-3 rows out of the fit. Matt's decision.

Transcribing the sources found two misprints in the MEP 532 nm tables, corrected in `bipm2012a` (R(122) 35-0
a7, R(145) 37-0 a5).

Not declared:
- `cheng2001a`, `simonsen2000a`, `fang2006a` and `arie1993a` measure the same components as BIPM tables but
  independently (their `meta.toml` says why). Each is its own record and enters the fit.
- `edwards1999a` is one of the sources averaged into `bipm2003a` Table 2 (MEP ref. [34], found 2026-10-21); see
  above.
- `yoon2001a`, `diddams2000a` and `smith2007a` (pending) measured a13 of the 633 nm standard. Yoon's own f − i
  would also claim the a13 − a16 row of the MEP's Table 1, which is not its measurement, so none is declared.
- `ma2006a` (pending) is probably the basis of the MEP 543 nm value (`bipm2003b`), but the MEP's source list was
  not compared; its absolutes would also claim the MEP's R(12) a15 − R(106) b10 interval.
- `ye2001a`'s year-long a10 record includes the comb runs `diddams2000a` reports; neither "same" nor "mean" fits.
- `yoshii2019a` shares its apparatus with `sakagami2020a`, but the two report different measurements.

`bodermann1998c` (the thesis) prints three of `bodermann2000a`'s intervals with the same digits. They are not
entered there, because the thesis's uncertainties for them were not transcribed, so the thesis is not yet
cited on those rows.

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

`data/b_levels/` (new 2026-10-21) holds the same kind of per-level constants for the B state. `king1980/` is King
et al., Chem. Phys. 50, 291 (1980): G_v and B_v of ¹²⁹I₂ B v = 71–79 (and v = 69, unused by the authors) and of
¹²⁷I¹²⁹I v = 70 and 73, from two-photon sequential spectroscopy within 2–25 cm⁻¹ of the dissociation limit
(columns in its `meta.toml`). It carries a "pending" `exclude` like the observation sets. The shipped model gives
these G_v to 0.028 cm⁻¹ rms. No code reads the directory yet.

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
