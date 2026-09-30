# The 778–795 nm absolute anchors: Bodermann *et al.*, EPJ D 11, 213 (2000)

*Written 2026-09-15. The paper was obtained and transcribed in full into
`data/observations/bodermann2000a/`. This file records the sources, the cell and reference
conditions, exactly what was transcribed and checked, every ambiguity found in the paper, and the
residuals of our model against these data.*

**Tagging.** **[S]** = stated by the paper. **[I]** = inference or a number computed here.

---

## 0. Bottom line

1. **32 observations** are now in `data/observations/bodermann2000a/`: 29 absolute hyperfine-component
   frequencies and 3 line-to-line intervals, over **21 distinct lines** in **777.9–795.0 nm**, bands
   **0-12, 0-13, 0-14, 0-15, 1-14, 2-15 and 3-16**. Stated uncertainties are **34–243 kHz** [S].
2. **Our model (Hannover 2008 potentials) is off by a few MHz, as expected** — mean **+1.4 MHz**,
   **rms 2.9 MHz** over the 29 frequencies, spanning **−6.2 to +7.0 MHz**. That is ~100× the
   measurement uncertainties (χ ≈ 104) [I].
3. **The Knöckel 2004 local Dunham model reproduces these same data to 93 kHz rms** (24 of 25 v′ = 0
   rows within 200 kHz) [I]. This is the paper's own fit data, so the agreement is expected — but it
   is also the strongest possible check that **the transcription is correct**, and it localises the
   2.9 MHz to the global potentials, not to our code, our hyperfine treatment, or these numbers.
4. **The ~50 MHz question near 755 nm is not touched by this paper.** Its two 0-12 lines are
   **J″ = 240 and 243**, which sit at **778.0 and 780.3 nm** — the high-J end of that band, not the
   low-J band head at 754.8 nm. There our model is good (−0.27 and +0.39 MHz) [I]. The 755 nm
   extrapolation remains open and still needs Liao 2010.
5. **The worst band is 1-14 at +7.0 MHz**, consistently on both of its lines. Bands with v′ > 0
   (1-14, 2-15, 3-16) are **outside the reach of either Hannover local model**, which are v′ = 0 only,
   so these four rows are the most valuable new constraints in the set [I].
6. Two genuine defects were found in the paper's tables (§4): the **+42 kHz cell correction is
   present in Table 2 and absent in Table 4**, and **Table 2's uncertainties (74 kHz) are smaller than
   the 80 kHz reference they rest on**. Nothing was corrected; both are recorded in `meta.toml`.

---

## 1. Sources

| Field | Value |
|---|---|
| Citation | B. Bodermann, M. Klug, U. Winkelhoff, H. Knöckel, E. Tiemann, *Eur. Phys. J. D* **11**, 213–225 (2000) |
| DOI | [10.1007/s100530070086](https://doi.org/10.1007/s100530070086) |
| Catalog id | `bodermann2000a` (`data/catalog/precision.yaml`, line 872) — the existing entry was used |
| Data set | `data/observations/bodermann2000a/` |
| Retrieved | 2026-09-15 |

This is reference **[17] of Knöckel 2004** and **R22 of Liao 2010**. It is the only absolute-frequency
data we hold between 671 nm and 815 nm, and it is the measurement set behind both Hannover
near-infrared local models.

## 2. Conditions, as printed [S]

**Iodine.** Five cells I1–I5, filled at PTB; I1/I2 are 50 cm (1994), I3–I5 are 80 cm (1996). The cell
body is heated to **600(30) °C** with one end stabilised at **20(1) °C**, giving a vapour pressure of
**27.0(2.4) Pa**. Collinear Doppler-free saturation spectroscopy with a **3f** lock; the lock-in phase
is switched 0°/180° to cancel electronic offsets. FWHM ≈ 10 MHz. Allan deviation 5 × 10⁻¹² at 100 s;
single-cell reproducibility 10 kHz.

**The cell correction — the only extrapolation in the paper.** Background gas red-shifts the lines,
and cell I5 gave the highest frequencies, so I5 is taken as purest and **all quoted frequencies refer
to I5** [S]. Cells I1 and I2 were measured **−42(5) kHz** below I5, so **+42 kHz** was added to refer
them to I5. The cell-to-cell spread over I1–I5 is **26 kHz (1σ)**, which the paper calls the floor on
absolute accuracy.

**There is no pressure or power extrapolation.** The frequencies are *not* reduced to zero pressure or
zero power; they are cell I5's values at 27.0(2.4) Pa with only that fixed +42 kHz offset. Power and
geometry shifts were characterised only for the *rubidium* references. This is recorded in
`meta.toml` under `[conditions]` as `pressure_extrapolation` / `power_extrapolation`.

**The rubidium references** the absolute scale rests on [S], all recorded in `meta.toml` under
`[references]`:

| Reference | Value (kHz) | Source |
|---|---|---|
| Rb D1, component c′ | 377 106 271 488 (80) | Barwood *et al.*, uncertainty enlarged 60 → 80 kHz |
| Rb D2, crossover d/f | 384 227 981 877 (22) | Ye *et al.* |
| ⁸⁵Rb 5s–5d ²D₃/₂, F_g = 2–F_e = 4 | 385 242 216 338 (23) | Nez / Felder, with a −25(10) kHz AC-Stark correction |
| ⁸⁷Rb 5s–5d ²D₅/₂, F_g = 2–F_e = 4 | 385 284 566 348 (24) | Nez / Felder, same correction |

A caveat the paper states plainly: the Rb D1 realisation carries a **~100(10) kHz residual
first-order Doppler shift** from its 10.0(5) mrad pump/probe angle. It cancels against Barwood's
identical geometry in the *transfer*, but not in the true I₂ frequency — a systematic sitting under
the whole `rb_d1` group.

## 3. What was transcribed, and what was checked

Typed from the `pdftotext` layer, then **checked digit by digit against 400 dpi page-image crops**
rendered with `pdftoppm`. All values are printed in kHz and divided by 1000 to give MHz (exact).

| Page | Table | Content | Rows |
|---|---|---|---|
| 216 | Table 2 | 4 lines on Rb D1, component c′ | 4 |
| 219 | Table 4 | 7 lines on Rb D2, crossover d/f | 7 |
| 222 | Table 7 | 5 lines on the Rb two-photon transition | 5 |
| 223 | Table 8 | 9 four-wave-mixing rows, 4 of which repeat a two-photon value | 9 + 4 |
| 223 | §6 | the three-step FWM interval chain | 3 |

**Checked on the page images, not only the text layer:** every digit of all 25 frequency values and
all 16 beat frequencies; every bracketed uncertainty; the sign of every negative entry (Table 2 rows
1–2, Table 4 rows 1–6, Table 7's P(19) and R(139)/⁸⁷Rb entries, and the −79/−61 kHz differences in
Table 8); the **blank cells** of Table 7 (P(19) and R(26) have only an ⁸⁵Rb beat, R(18) and R(16) only
an ⁸⁷Rb beat, and R(139) 1-14 alone has both); Table 8's layout, where `transition` spans several
`comp.` rows, and its two `*` footnote marks; the four Rb reference values on pp. 216, 219 and 221;
and Figure 5 with its caption.

**Arithmetic cross-checks, all of which passed** — these confirm the digits independently of the OCR:

* Table 2 frequency = Rb D1 + beat **+ 42 kHz**, exactly, all four rows.
* Table 4 frequency = Rb D2 + beat, exactly, all seven rows, with **no** +42 kHz (§4.1).
* Table 7 = reference + beat exactly for the four single-reference rows; for R(139) 1-14 the printed
  value is the **mean** of the ⁸⁵Rb and ⁸⁷Rb determinations (…722.1 and …727.9).
* The two R(139) 1-14 beats reproduce the paper's Rb fine-structure interval 42 350 004.2(88) kHz.
* The §6 chain sums to the printed total 61 215 408 kHz, and its uncertainties combine in quadrature
  to the printed 42 kHz.
* Table 8's `difference` column equals FWM − 2-photon for all four rows (+72, +75, −79, −61).

**Kinds.** Every published value is a **hyperfine component** — the paper prints no hyperfine-free
line centres, so no row has an empty `component`. The 29 table values are `kind = frequency`; the 3
§6 values are `kind = interval` with an explicit `ref_line` and `ref_component`. The beat frequencies
are I₂ *minus rubidium*, which the schema cannot express (a reference must be an I₂ line), so each is
preserved verbatim in its row's `note`.

**Groups**, one per independent calibration route, so the fit can give each its own offset:

| group | n | Tied to |
|---|---|---|
| `rb_d1` | 4 | Barwood's Rb D1 realisation (80 kHz) |
| `rb_d2` | 7 | Ye's Rb D2 d/f value (22 kHz) |
| `rb_2photon` | 5 | the Nez/Felder two-photon values (23–24 kHz) |
| `fwm` | 9 | `rb_d2` plus the ref. [13] wavelength calibration |
| `rb_2photon_derived` | 4 | Table 8's second column — **not independent**, see below |
| `fwm_chain` | 3 | line-to-line only; carries no Rb reference |

**Not entered**, and why: the §6 *total* 61 215 408(42) kHz (it is the sum of the three chain
intervals already entered); the Rb reference and fine-structure values (atomic, unrepresentable);
R(96) 0-14 and P(87) 0-14 from Figure 7 (their frequencies come from ref. [13], Bodermann 1998, not
from this paper); and Figure 7's 6-digit THz axis labels, which are locators, not measurements.

**Redundancy [I].** The four `rb_2photon_derived` rows are re-derivations, not new measurements:

* P(164) 0-13 a10 = 385 363 916 113(55) is *exactly* the R(16) 0-14 a10 row plus the three chain
  intervals (385 302 700 705 + 61 215 408), its uncertainty being 35 ⊕ 42 kHz.
* P(164) 0-13 a1 is that value minus an unpublished a10 − a1 interval of 588 330 kHz.
* P(19) 0-14 a14 and a16 are the Table 7 a13 row plus unpublished intervals of 60 573 and
  116 075 kHz.

Dropping the whole `rb_2photon_derived` group leaves an uncorrelated set; keeping it adds only those
three unpublished intra-line hyperfine intervals. The `fwm` rows for the same components *are* a
genuinely different route and are kept.

## 4. Ambiguities and apparent misprints

Nothing below was "corrected" in the data. All are recorded in `meta.toml`.

### 4.1 The +42 kHz cell correction is in Table 2 but not in Table 4 [I]

The text says the Table 2 measurements used cell I1 and the Table 4 measurements used cell I2, and in
**both** places that "the absolute frequencies are corrected by +42 kHz". The arithmetic disagrees:

* Table 2: frequency − (Rb D1 + beat) = **+42.0 kHz**, exactly, for all four rows.
* Table 4: frequency − (Rb D2 + beat) = **0.0 kHz**, exactly, for all seven rows.

Either Table 4's beat column is already corrected, or the correction was dropped from Table 4. At
42 kHz this **exceeds the 34–55 kHz uncertainties** of those rows, so it is not a rounding artefact.
Entered as printed. If it is an omission, the whole `rb_d2` group is 42 kHz low.

### 4.2 Table 2's uncertainties are smaller than their own reference [I]

Table 2 quotes 74–75 kHz, but every one of those values rests on the Rb D1 reference whose
uncertainty the same page gives as **80 kHz** (60 ⊕ 15 ⊕ 51 kHz, stated explicitly). An absolute
frequency cannot be more certain than its reference. The 74 kHz could not be reconstructed from any
stated combination. Entered as printed; the `rb_d1` group is the least trustworthy in the set.

### 4.3 Figure 5's arrow for P(19) 0-14 points at a16, while Table 7 tabulates a13 — both are right

The "calibrated:" arrow in Fig. 5 sits over **a16**, but Table 7 lists **a13** for that line. This is
not a misprint: a16 is one of the two P(19) components in **Table 8's** two-photon column, so the
figure marks a Table 8 component rather than the Table 7 one. The model settles the labelling
decisively [I]: Tables 7 and 8 imply a14 − a13 = 60.573 MHz, our model predicts **60.557 MHz**, and
**a13 is the only one of the 21 components within 25 MHz** of that spacing. No relabelling was needed.

### 4.4 Minor

Table 7's caption names the reference only as "the Rb 5s ²S₁/₂–²D two photon transition" without
saying which fine-structure component each column uses; the column headers do say, and are used here.

### 4.5 A positive check on the labels [I]

All nine intra-line hyperfine spacings that the tables imply agree with our model to **15–85 kHz**:

| Line | Spacing | Observed (MHz) | Model (MHz) | o − m (kHz) |
|---|---|---|---|---|
| P(19) 0-14 | a14 − a13 | 60.573 | 60.557 | +16 |
| P(19) 0-14 | a16 − a13 | 116.075 | 116.094 | −19 |
| P(19) 0-14 | a16 − a14 | 55.502 | 55.536 | −34 |
| P(164) 0-13 | a10 − a1 | 588.330 | 588.366 | −36 |
| P(164) 0-13 | a15 − a10 | 294.236 | 294.254 | −18 |
| R(56) 0-14 | a15 − a10 | 294.447 | 294.362 | +85 |
| P(31) 0-14 | a7 − a4 | 119.011 | 119.026 | −15 |

This confirms every component label, and confirms that the paper's numbering convention (Fig. 5
caption: "numbered 1…15 for even J″ and 1…21 for odd J″", increasing with ν) is **i2spec's own**
rank-by-frequency convention. Labels were therefore carried over unchanged. Every label is also
consistent with the parity of its J″ (a1–a21 / b1–b13 on odd-J″ lines, a1–a15 on even-J″ lines).

## 5. Residuals of our model against these data

Computed with `i2spec.observations.Predictor()` (Hannover 2008 potentials, B-spline solver, ΔJ = 2)
and `residuals()`. **No model parameter was tuned to these data.**

### 5.1 Full table

`o − m` in kHz. "local" is the Knöckel 2004 local Dunham model (§5.4), which is v′ = 0 only.

| Line | comp | group | λ (nm) | observed (MHz) | u (kHz) | o − m (kHz) | local (kHz) |
|---|---|---|---|---|---|---|---|
| P(105) 0-15 | a14 | rb_d1 | 795.021 | 377 087 649.355 | 74 | **+2155** | +10 |
| R(205) 0-14 | a2 | rb_d1 | 794.985 | 377 104 695.663 | 75 | **−6219** | +37 |
| R(113) 0-15 | a12 | rb_d1 | 794.957 | 377 117 726.169 | 74 | **+2672** | +124 |
| P(104) 0-15 | a10 | rb_d1 | 794.929 | 377 131 003.300 | 74 | **+2231** | +86 |
| P(70) 0-14 | a10 | rb_d2 | 780.293 | 384 205 192.7886 | 36 | **+2816** | −41 |
| R(188) 0-13 | a10 | rb_d2 | 780.286 | 384 208 473.5014 | 36 | **−2011** | +12 |
| R(117) 2-15 | b13 | rb_d2 | 780.285 | 384 208 927.420 | 52 | −231 | — |
| P(43) 3-16 | a2 | rb_d2 | 780.265 | 384 218 608.088 | 45 | **−3213** | — |
| P(243) 0-12 | b1 | rb_d2 | 780.265 | 384 218 812.773 | 55 | −268 | +361 |
| R(78) 0-14 | a10 | rb_d2 | 780.257 | 384 222 471.5179 | 34 | **+3387** | +1 |
| P(148) 1-14 | a1 | rb_d2 | 780.231 | 384 235 614.6238 | 36 | **+7002** | — |
| R(56) 0-14 | a10 | fwm | 779.120 | 384 783 612.796 | 117 | +2903 | −45 |
| R(56) 0-14 | a15 | fwm | 779.119 | 384 783 907.243 | 118 | +2988 | +40 |
| P(31) 0-14 | a4 | fwm | 778.532 | 385 074 117.864 | 138 | +2144 | +99 |
| P(31) 0-14 | a7 | fwm | 778.532 | 385 074 236.875 | 144 | +2129 | +84 |
| P(19) 0-14 | a13 | rb_2photon | 778.240 | 385 218 524.993 | 38 | **+1893** | +21 |
| P(19) 0-14 | a14 | rb_2photon_derived | 778.240 | 385 218 585.566 | 40 | +1909 | +36 |
| P(19) 0-14 | a14 | fwm | 778.240 | 385 218 585.638 | 140 | +1981 | +108 |
| P(19) 0-14 | a16 | rb_2photon_derived | 778.240 | 385 218 641.068 | 40 | +1875 | +2 |
| P(19) 0-14 | a16 | fwm | 778.240 | 385 218 641.143 | 140 | +1950 | +77 |
| R(26) 0-14 | a10 | rb_2photon | 778.210 | 385 233 438.181 | 34 | **+2155** | −21 |
| R(139) 1-14 | a13 | rb_2photon | 778.143 | 385 266 729.725 | 35 | **+7048** | — |
| R(18) 0-14 | a10 | rb_2photon | 778.091 | 385 292 072.466 | 35 | **+2009** | −7 |
| R(16) 0-14 | a10 | rb_2photon | 778.070 | 385 302 700.705 | 35 | **+1998** | +16 |
| P(164) 0-13 | a1 | fwm | 777.948 | 385 363 327.704 | 242 | +63 | −61 |
| P(164) 0-13 | a1 | rb_2photon_derived | 777.948 | 385 363 327.783 | 56 | +142 | +18 |
| P(164) 0-13 | a10 | fwm | 777.946 | 385 363 916.052 | 242 | +45 | −79 |
| P(164) 0-13 | a10 | rb_2photon_derived | 777.946 | 385 363 916.113 | 55 | +106 | −18 |
| P(164) 0-13 | a15 | fwm | 777.946 | 385 364 210.288 | 243 | +27 | −97 |

Intervals (`kind = interval`):

| Interval | observed (MHz) | u (kHz) | o − m (kHz) |
|---|---|---|---|
| P(164) 0-13 a10 − R(240) 0-12 a10 | 16 324.846 | 27 | −279 |
| R(240) 0-12 a10 − R(138) 1-14 a10 | 24 845.370 | 30 | **−6536** |
| R(138) 1-14 a10 − R(16) 0-14 a10 | 20 045.192 | 12 | **+4924** |

### 5.2 Summary

| Subset | n | mean (kHz) | rms (kHz) | range (kHz) |
|---|---|---|---|---|
| all rows | 32 | +1244 | **3126** | −6536 … +7048 |
| `kind = frequency` | 29 | +1437 | **2911** | −6219 … +7048 |
| `kind = interval` | 3 | −630 | 4727 | −6536 … +4924 |
| independent rows (no `rb_2photon_derived`) | 28 | +1277 | 3303 | −6536 … +7048 |

By group: `rb_d1` rms 3723, `rb_d2` rms 3442, `rb_2photon` rms 3631, `fwm` rms 1950,
`rb_2photon_derived` rms 1341, `fwm_chain` rms 4727 kHz.

**χ = √⟨(o − m)²/u²⟩ = 104.** The residuals are ~100× the measurement uncertainties, so these data
constrain the model overwhelmingly more tightly than it currently satisfies them.

### 5.3 The pattern with v′, v″ and J

**By band** (frequency rows):

| Band | n | J″ range | mean (kHz) | sd about the mean (kHz) |
|---|---|---|---|---|
| 0-12 | 1 | 243 | −268 | — |
| 0-13 | 6 | 164–188 | −271 | 779 |
| 0-14 | 15 | 16–205 | +1728 | 2175 |
| 0-15 | 3 | 104–113 | +2353 | 228 |
| **1-14** | 2 | 139–148 | **+7025** | **23** |
| 2-15 | 1 | 117 | −231 | — |
| 3-16 | 1 | 43 | −3213 | — |

**By upper state:** v′ = 0 → mean +1243, rms 2338 kHz (n = 25); **v′ = 1 → +7025 kHz** (n = 2);
v′ = 2 → −231 kHz (n = 1); v′ = 3 → −3213 kHz (n = 1).

**By lower state:** v″ = 12 → −268; v″ = 13 → −271; v″ = 14 → +2351; v″ = 15 → +1707;
v″ = 16 → −3213 kHz.

Three things stand out [I]:

1. **The 1-14 band is offset by +7.0 MHz, and remarkably consistently** — its two lines, R(139) and
   P(148), agree with each other to 23 kHz while both sit 7 MHz from the model. That is the signature
   of a **band-origin error**, not a rotational one. It is the largest discrepancy in the set, and it
   is confirmed independently by the FWM chain (§5.5).
2. **Within 0-14 the residual grows smoothly with J″ and then turns over hard.** From J″ = 16 to 78 it
   climbs monotonically, +1998 → +2009 → ~+1900 → +2155 → +2137 → +2946 → +2816 → +3387 kHz, a slope
   of **+2.3 MHz per 100 J″**. But at **J″ = 205 it is −6219 kHz** — 8 MHz below where that trend
   points. The same turnover appears in 0-13 (+100 kHz at J″ = 164, −2011 kHz at J″ = 188). The model
   errs in the **rotational** description at high J″, and with a sign change.
3. **The v′ = 0 bands closest to the middle of the range are the best.** 0-13 at J″ = 164 and 0-12 at
   J″ = 240–243 are within 0.4 MHz. The error is not a uniform offset over the region.

### 5.4 The same data against the Knöckel 2004 local Dunham model

The local model (`tests/test_hannover2008.py`, `S_B`/`S_X`, Knöckel 2004 Table 1) is v′ = 0 only, so it
covers 25 of the 29 frequency rows. Its centres were combined with **our** hyperfine offsets, which
§4.5 shows are good to < 100 kHz.

| Model | n | mean (kHz) | rms (kHz) | max abs (kHz) |
|---|---|---|---|---|
| Hannover 2008 potentials (ours) | 25 | +1243 | **2338** | 6219 |
| Knöckel 2004 local Dunham | 25 | +31 | **93** | 361 |

22 of 25 rows are within 100 kHz of the local model and 24 of 25 within 200 kHz — matching the 2004
paper's own claim of "1σ accuracies well below 200 kHz". **The single outlier is P(243) 0-12 at
+361 kHz, whose J″ = 243 is just past the model's stated J″ ≤ 242 validity limit** [I] — a pleasingly
self-consistent failure.

This is circular as a test of the *local model* (these are its fit data), but it is decisive as a test
of the *transcription*: a mis-typed digit, a wrong component label or a wrong branch/J″ assignment
would show up here as a multi-MHz outlier, and none does. It also localises our 2.3 MHz cleanly to
the **global potentials**, and not to the solver, the hyperfine code, or these numbers.

### 5.5 The FWM chain is internally consistent with the absolute residuals

Chaining our residual for R(16) 0-14 a10 through the three §6 intervals gives implied absolute
residuals for the two lines that have no absolute entry [I]:

| Line | implied o − m (MHz) | independent check |
|---|---|---|
| R(16) 0-14 a10 | +1.998 | measured directly |
| R(138) 1-14 a10 | **+6.921** | the 1-14 band directly gives +7.002 and +7.048 |
| R(240) 0-12 a10 | +0.385 | — |
| P(164) 0-13 a10 | **+0.106** | measured directly: +0.045 (FWM), +0.106 (derived) |

The chain closes to ~100 kHz at both ends. This independently confirms the 1-14 band's +7 MHz
offset, confirms the three interval values, and confirms that the large interval residuals (−6.5 and
+4.9 MHz) are just the band offsets differencing — not a defect in the interval data.

### 5.6 What about the ~50 MHz near the 755 nm band head?

**These data do not test it.** The expectation recorded in `nir-model-2010.md` §6.1 is that our model
and the 2004 local model diverge to −48.9 MHz at the **0-12 band head**, which is at **754.8 nm and
low J″**. This paper's only 0-12 lines are **P(243) and R(240)** — J″ = 243 and 240, at **780.3 and
778.0 nm**, the far high-J end of that band. Our model there is **−0.27 and +0.39 MHz**, and
potentials − local is **+0.63 MHz** at P(243), not tens of MHz.

So the 755 nm question is untouched and still open; it needs Liao 2010's (0–12) and (0–13) low-J
measurements. What these data *do* newly show is a **different** large error, +7 MHz in the 1-14
band, in a region everyone assumed was well covered.

---

## 6. What this changes

1. **The red gap now has a real anchor.** Before this, our only absolute check between 671 and 815 nm
   was Huang 2013 at 671 nm. We now have 29 absolute frequencies at 34–243 kHz across 778–795 nm.
2. **Our model's accuracy here is quantified: a few MHz, and up to 7 MHz.** The catalog's
   "uncertainty_MHz: 0.08" describes the *measurements*, not our agreement with them.
3. **The case for promoting the local NIR model is now data-backed**, not inferred. `nir-model-2010.md`
   §7 Step 0 proposed promoting the Knöckel 2004 local Dunham model out of test-only code. §5.4 above
   shows it is **25× better than the potentials on real data** in its domain (93 kHz vs 2338 kHz rms).
4. **The v′ > 0 rows are uniquely valuable.** 1-14, 2-15 and 3-16 lie outside both Hannover local
   models (v′ = 0 only), so the +7.0 and −3.2 MHz they expose can only be fixed by the global
   potentials or by a Phase B refit. They are the four most informative rows in the set.
5. **A high-J rotational defect is now visible**, from the sign change between J″ ≈ 80 and J″ ≈ 200
   within 0-14 and 0-13 (§5.3). This is testable against other high-J data we already hold.

## 7. Suggested changes outside this file's scope

Reported, not made — these files belong to other work.

* **`data/catalog/precision.yaml`, entry `bodermann2000a`:** set `machine_readable: true`,
  `n_items: 32`, and add a pointer to `data/observations/bodermann2000a/`. Three fields are wrong as
  written: `range_nm: [778.0, 815.0]` should be **[777.9, 795.1]** (778–815 nm is the *2004 model's*
  validity range, not this paper's data); `transitions: "NIR hyperfine components (bands 0-12 ...
  0-17)"` should be **bands 0-12, 0-13, 0-14, 0-15, 1-14, 2-15, 3-16** — the paper has no v″ = 17, and
  it does have v′ = 1, 2 and 3, which the entry does not mention; and `uncertainty: "mostly much
  smaller than 80 kHz"` is optimistic, the true spread being **34–243 kHz**.
* **`docs/research/nir-model-2010.md`:** §4 still says Bodermann 2000 is "worth requesting", §5 lists
  it under "Not held", and §8 under "Worth adding once obtainable". It is now held and transcribed;
  those three places should point here. §6.1's −48.9 MHz expectation should be annotated with §5.6
  above: the 2000 data do **not** probe the 755 nm head.
* **`src/i2spec/lookup.py` and a new `src/i2spec/local_nir.py`:** §5.4 is the missing empirical
  justification for Step 0 of `nir-model-2010.md` §7. `lookup.uncertainty(line)` currently reports a
  global-potential estimate in this window; on this evidence it understates the error for v′ > 0 and
  overstates it for v′ = 0 where the local model applies.
* **A regression test.** `bodermann2000a` is a good acceptance test for any NIR model change: assert
  the local Dunham model stays within 200 kHz rms of the v′ = 0 rows, and record the current
  potential-model rms of 2.9 MHz so a refit's improvement is visible.
