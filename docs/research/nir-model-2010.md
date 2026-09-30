# The 2010 near-infrared model: Liao *et al.*, JOSA B 27, 1208

*Implemented 2026-09-15 from the published paper (journal pp. 1208–1214). Everything below was typed from the `pdftotext` layer and then **checked
digit by digit against 300 dpi page renders** (`pdftoppm`) of Table 1, Table 2 (both halves),
Table 3 and the surrounding text. Knöckel 2004 Table 1 was re-checked the same way.*

**Tagging.** **[S]** = stated by the paper, with page. **[I]** = computed or inferred here.

---

## 0. Bottom line

1. **The model is what §2 of the earlier draft predicted [S, p. 1212 §4]:** a **refit of the local
   Dunham model of Knöckel 2004** (its eqs. 1–2), *not* new potentials, *not* a correction layer on
   top of them. Same functional form, **same 19 parameters**, a few indices changed, fitted to the
   27 new frequencies plus the 2004 NIR data. Normalised standard deviation **σ = 0.82**.
2. **It is published in a form that cannot be used.** Table 3 prints 3–9 significant digits where
   Knöckel 2004 printed 7–13. The parameters are strongly correlated, so the rounding destroys the
   model: as printed it misses **its own paper's measurements by 37.7 MHz rms (74.9 MHz worst)**,
   against the 0.2 MHz it claims (§3).
3. **The good news is much bigger than the bad.** The 31 comb-referenced frequencies are now in
   `data/observations/liao2010a/`, and against them **our Phase A potentials are within 4.3 MHz
   (2.4 MHz rms)** across 755–810 nm (§6) — not the 30–60 MHz the literature assumes for this
   region. The paper's own potential-model column says the same thing (≤4.4 MHz) [S, Table 2].
4. **The 755 nm band-head question is settled, and the old answer was wrong.** The −48.9 MHz that
   the earlier draft reported at the (0-12) head was the **2004 local Dunham model extrapolating**,
   not an error in our potentials. At R(48) 0-12 (755.6 nm) our potentials are **−0.13 MHz** from
   the measurement while the 2004 local model is **+45.1 MHz** off (§6.1).
5. **What we now ship:** `src/i2spec/local_nir.py` with both parameter sets as data files; the 2004
   set is the usable one and is **10× better than the potentials inside its 775–815 nm window**
   (0.23 MHz rms vs 2.31 MHz on 34 measured anchors). `lookup.uncertainty` is rewired accordingly
   (§7).
6. **Two misprints in the paper** [I]: Table 2's "P(100) 0–12" is **R(100) 0–12** (§4.3), and its
   caption cites "[29]" for both model columns where the text means **[16]** (§5).

---

## 1. The paper

| Field | Value [S] |
|---|---|
| Authors | Chun-Chieh Liao, Kuo-Yu Wu, Yu-Hung Lien (NTHU); Horst Knöckel (Hannover); Hsiang-Chen Chui (NTHU, then NCKU); Eberhard Tiemann (Hannover); Jow-Tsong Shy (NTHU) |
| Dates | received 3 Dec 2009; revised 19 Mar 2010; accepted 22 Mar 2010; published 12 May 2010 (Doc. ID 120901) |
| Abstract claim | 27 hyperfine transitions of the (0–12) and (0–13) bands, 750–780 nm, relative uncertainty a few × 10⁻¹⁰; deviations from the 2004 model "larger than expected"; "**an improved model is developed for the range from 755 to 815 nm** … error limit of the absolute frequency less than **0.2 MHz**" |
| Motivation | p. 1208: in 667–776 nm the only references are the Gerstenkorn atlases, "**not better than 30–90 MHz**". The new work is "an improvement of the prediction uncertainty by **at least a factor of 100** in the wavelength range between 755 and 778 nm" |

---

## 2. The model

### 2.1 Form [S, p. 1212 §4]

> "The Dunham parameters are defined in a power series of quantum numbers ν + 1/2 and J(J + 1)
> according to formulas (1) and (2) of [16]. We use the letter S instead of the conventional
> nomenclature Y to indicate clearly that the derived parameters do not have their normal meaning,
> but are applicable strictly only in a limited interval of quantum numbers."

So the model is Knöckel 2004 eq. (1), with the upper state collapsed by eq. (2) because only
v′ = 0 occurs:

**ν(v′ = 0, v″, J″, branch) = Σ_k S₀ₖ [J′(J′+1)]ᵏ − Σ_{l,k} S_lk (v″+½)^l [J″(J″+1)]ᵏ**,
with J′ = J″ + 1 for R and J″ − 1 for P, all coefficients in cm⁻¹.

### 2.2 Table 3, transcribed exactly [S, p. 1213]

*"Improved Dunham Parameters for High Precision Description of the Wavelength Range 755–815 nm for
the Bands (0–12) to (0–17) of ¹²⁷I₂." Footnote: "For Details See Text. All Values are Given in
cm⁻¹."* Uncertainties are as printed, in the last digits.

| Upper state | Value (cm⁻¹) | | Lower state | Value (cm⁻¹) |
|---|---|---|---|---|
| S₀₀ | 1.583 177 77(95) × 10⁴ | | S₁₀ | 2.145 488(26) × 10² |
| S₀₁ | 2.892 599 03(26) × 10⁻² | | S₂₀ | −6.1206(26) × 10⁻¹ |
| S₀₂ | −6.233 13(29) × 10⁻⁹ | | S₃₀ | −9.57(12) × 10⁻⁴ |
| S₀₃ | −2.3460(82) × 10⁻¹⁵ | | S₄₀ | −1.653(20) × 10⁻⁵ |
| S₀₄ | −3.551(71) × 10⁻²¹ | | S₀₁ | 3.737 090 3(98) × 10⁻² |
| | | | S₁₁ | −1.145 46(21) × 10⁻⁴ |
| | | | S₂₁ | −2.326(15) × 10⁻⁷ |
| | | | S₃₁ | −6.115(34) × 10⁻⁹ |
| | | | S₀₂ | −4.522 49(52) × 10⁻⁹ |
| | | | S₁₂ | −2.7825(41) × 10⁻¹¹ |
| | | | S₃₂ | −1.9829(65) × 10⁻¹⁴ |
| | | | S₀₃ | −2.658(82) × 10⁻¹⁶ |
| | | | S₁₃ | −4.536(10) × 10⁻¹⁷ |
| | | | S₀₄ | −6.99(72) × 10⁻²² |

**19 parameters**, 5 upper and 14 lower — exactly the count of Knöckel 2004 Table 1 [S, p. 1212:
"a good fit with a normalized standard deviation of σ = 0.82 with the same number of parameters as
in [16], but few parameters have different indices than before"]. The index changes are **S₃₀ and
S₃₂ in place of S₂₂ and S₁₄** [I, by comparing the two tables].

Stored verbatim in `src/i2spec/data/liao2010_nir.json`; the 2004 set is in
`src/i2spec/data/knockel2004_nir.json`.

### 2.3 Validity [S]

| Limit | Value | Source |
|---|---|---|
| Electronic/vibrational | v′ = 0; v″ = 12…17 | Table 3 title |
| Wavelength | **755–815 nm** | Table 3 title, abstract, §5 |
| Claimed accuracy | **< 0.2 MHz** absolute; below 100 kHz for most (0-12)/(0-13) lines | abstract; p. 1213 |
| Rotational | Fig. 4 spans J″ = 20…240; **"it is recommended to use lines with rotational quantum numbers lower than 140 for the calibration"** (above that the (0-12)/(0-13) uncertainty rises to a few × 100 kHz) | p. 1213 |
| Hard edge | **"For the bands of lower vibrational quantum numbers in the ground state, such as (0–11) and below … An extrapolation with the improved model presented here is not recommended."** | p. 1213 |

We carry J″ ≤ 242 in the domain metadata, from the 2004 fit's own limit; the 2010 paper states no
J″ ceiling of its own.

### 2.4 Relation to the 2004 and 2008 models — confirmed

* **It replaces the 2004 *local* model, and only that.** The 2004 local model is the one "also used
  in the program 'IodineSpec'" [S, p. 1212].
* **It does not touch the potentials.** The potential model is discussed separately [S, p. 1212] and
  is left as it stands; the paper reports its performance rather than changing it.
* The paper quantifies where the 2004 local model may be used: it "used high precision data of
  (0–17) to (0–14) bands and additionally **three high precision values at high rotational quantum
  numbers J″ > 160 for the bands (0–13) and (0–12)**… The model was assumed to be valid in the range
  from **776 to 815 nm**" [S, p. 1212]. For (0-12) the fitted lines are **R(240) and P(243)**, for
  (0-13) **P(164) and R(188)** [S, p. 1212] — i.e. the new measurements at J″ = 48–126 in (0-12) are
  far outside anything the old fit saw, which is exactly why it fails there.
* **The paper's verdict on its predecessor** [S, p. 1212]: "the Dunham parameter model should not be
  used beyond the limits given in [16], and one should select the calibration lines from the band
  (0–14) to (0–16) whenever possible."

### 2.5 Hyperfine treatment [S, p. 1211 §3]

> "For this comparison, the hyperfine structure was subtracted from the measured frequencies
> according to the models in [15,33]."

[15] = Salumbides *et al.* 2006, [33] = Bodermann *et al.* 2002 — the same formulae
`src/i2spec/hfs_params.py` implements. **The paper publishes no new hyperfine parameters**; the
local model is hyperfine-free and the components come from those formulae. Our own offsets agree
with the paper's to **15–137 kHz** per line (§6.3), which is what makes the comparison below
meaningful.

---

## 3. The defect: the printed table cannot reproduce the model

### 3.1 The symptom [I]

Evaluated as printed, Table 3 misses the paper's own 31 measured components by **37.7 MHz rms,
74.9 MHz at worst** (§6.2) — 200× its claimed error limit, and **worse than the 2004 model it
replaces** (21.6 MHz rms on the same rows). The fit itself is fine: σ = 0.82 and Fig. 4 shows
sub-100 kHz prediction uncertainties. The loss is in the printing.

### 3.2 Three independent demonstrations that rounding is the cause [I]

1. **Digit counts.** Table 3 gives 3–9 significant digits (S₃₀ has three: −9.57 × 10⁻⁴). Knöckel
   2004 Table 1 gives 7–13, and says so deliberately: "The number of digits given in the table
   corresponds to those necessary for a round off smaller than few kHz" [S, 2004 p. 201]. Liao 2010
   carries no such statement.
2. **Monte Carlo over the printed rounding intervals.** Sampling each parameter uniformly within
   ±½ of its last printed digit moves the predicted frequencies of these lines by
   **23–33 MHz (1σ)** — the size of the observed miss. The same exercise on the 2004 table moves
   nothing (< 0.005 MHz).
3. **Round the 2004 table to the 2010 table's digit count** and it breaks the same way: its
   predictions shift by **−14.6 to −20.0 MHz** on these lines. So the mechanism is the printing, not
   the 2010 fit.

**Confirmation:** a bounded least-squares fit that is allowed to move every parameter only *within
its own printed rounding interval* brings the residuals from 36.5 MHz rms down to **1.9 MHz rms**.
The printed digits are consistent with a model that fits — they just cannot be reconstructed from
what is on the page. Recovering it exactly would need the unrounded parameters or the covariance
matrix, neither of which is published; the paper has no supplement.

### 3.3 What we do about it

`load_local_nir("liao2010")` **raises** unless the caller passes `trust_printed_digits=True`. The
numbers ship for the record and for the domain/uncertainty metadata, and `tests/test_local_nir.py`
pins the defect so that a future, better parameter source is noticed immediately.

---

## 4. The measurements

### 4.1 Method and conditions [S, §2, pp. 1208–1210]

| Item | Value |
|---|---|
| Laser | Ti:sapphire (Coherent MBR-110), locked by Doppler-free FM saturation spectroscopy |
| Cell | 60 cm, Opthos Inc., degassed three times before filling; body ≈ **850 K**; cold finger **293.0 ± 0.1 K** → **27.0 ± 0.2 Pa** (Gillespie & Fraser formula) |
| Optics | pump ≈ 40 mW, probe phase-modulated (EOM New Focus 4002), pump chopped at 45 kHz by an AOM, beam ⌀ 2.8 mm |
| Line shape | a10 linewidth ≈ 5 MHz, S/N ≈ 200 in 125 Hz; residual lock jitter < 25 kHz |
| Comb | 1 GHz Kerr-lens mode-locked Ti:sapphire (Gigaoptics Gigajet20), PPLN f–2f, GPS-disciplined Rb (SRS PRS10); comb accuracy better than 1 × 10⁻¹² at 1000 s |
| Counting | Agilent 53132A, 0.2 s gate, 5000 points per measurement, binned into 100 and fitted with a Gaussian; σ is the uncertainty |
| Reproducibility | R(78) 0-14 a10 measured 8× over 3 days: per-run SD 10–13 kHz, SD of the eight means **8 kHz** |

**Uncertainty budget** [S, pp. 1210–1211]: the parenthesised numbers are "the summation of the
standard deviation of the frequency measurement and the possible shift caused by the offset of
electronics (5 kHz)". Neglected: comb (< 1 kHz), residual Doppler background (< 1 kHz), cold-finger
temperature (≈ 1 kHz, via the ≈ 5 kHz/Pa NIR pressure shift). But the cell's own impurity shift is
not characterised, so: *"For the worst case (121 kHz for the pressure shift in Table 1 and 50 kHz
for the measurement uncertainty in Table 2), the accuracy of the iodine transition frequency is
171 kHz … **an overall accuracy of 200 kHz is assigned to the transition frequencies in Table 2**."*
That 200 kHz is what `data/observations/liao2010a/data.csv` carries as the uncertainty; the printed
per-row values (8–50 kHz) are kept in the note column.

### 4.2 The three data groups, all now in `data/observations/liao2010a/` (31 rows)

* **Table 2 (p. 1211), 27 rows.** a1, a10 and a15 of nine lines with even J″, bands (0-12) and
  (0-13). Published *without* the pressure correction (footnote b).
* **Table 1 (p. 1211), 3 rows.** a10 of R(26), P(70), R(78) of (0-14) — the cell-calibration lines,
  measured against Bodermann *et al.* 2000. They come out **−121, −118, −102 kHz** below the
  five-cell PTB values, "a nearly constant redshift (average of −114 kHz)".
* **§2 text (p. 1210), 1 row.** a2 of R(141) 0-16 = **370 007 070.341(10) MHz**, the 2007 comb
  remeasurement of a line the group had determined in 2002 as 370 007 070.355(40) MHz through a
  CH₄-referenced difference-frequency chain. The 14 kHz agreement is the paper's evidence that the
  cell has not aged.

### 4.3 The mislabelled line [I]

**Table 2's "P(100) 0–12" is R(100) 0–12.** Our potential model puts P(100) 0-12 at
394 874 166 MHz — **347 GHz** from the printed value — while R(100) 0-12 lands **2.6 MHz** away, in
family with every other row (−3.5 to +1.2 MHz). A scan of both branches of both bands over
J″ = 40–180 finds no other candidate within 9 GHz. The hyperfine spacings confirm it: measured
a10 − a1 = 588.665 and a15 − a1 = 883.117 MHz, against our R(100) values 588.686 and 883.160 MHz
(+21 and +43 kHz). The row also keeps the table's monotonic frequency ordering. Entered as
R(100) 0-12, flagged in the row's note and documented in `meta.toml`.

---

## 5. What the paper reports about the two Hannover models [S, §3]

The 3σ prediction uncertainties it quotes for the **potential** model: **4.5 MHz** in 526–667 nm and
776–815 nm, "while it is expected to be less accurate with 3σ of **45 MHz or even larger** in the
range between 667 and 776 nm", because there the fit was filled with atlas lines of 30–90 MHz.

Its verdict, which our own numbers reproduce: *"Obviously, all predictions are fairly close to the
measurements; the uncertainty range of 4.5 MHz is never exceeded. However, with the model approach
and the data used here, the prediction is actually no longer an extrapolation as it was with the
local Dunham approach. Here it is an interpolation between two ranges where precise data exist."*

Measured − model, from Table 2 [S]: the **potential** model runs −4441 … +709 kHz; the **2004 local
Dunham** model runs −22 kHz at P(176) 0-13 up to **+37 924 kHz at R(48) 0-12**. The conclusion
(p. 1213) puts it plainly: the refit "reduces the deviations between measurements and predictions
from ≈38 MHz to less than 90 kHz".

*Caption misprint [I]:* Table 2 is titled "…the Predictions from the Dunham Parameter **[29]** and
Potential **[29]** Model Approaches". [29] is Bodermann *et al.* 2000; §3 makes clear both columns
come from **[16]**, Knöckel *et al.* 2004.

---

## 6. Our numbers

Computed with the committed code — `RovibronicModel("127I2")` (Hannover 2008 potentials, B-spline
solver, ΔJ = ±2 hyperfine) and `i2spec.local_nir` — against
`data/observations/liao2010a` and `data/observations/bodermann2000a`. **Sign convention here is
model − measurement**, the opposite of the paper's Table 2. The measurement is the published value
**plus the 114 kHz** pressure correction that the paper's own fit applied.

### 6.1 The 755 nm band head: the earlier draft's −48.9 MHz explained

The previous version of this file predicted a ~49 MHz problem at the (0-12) band head, from the
disagreement between our potentials and the 2004 local Dunham model. Liao's four (0-12) lines land
exactly there and settle it:

| Line | λ (nm) | potentials − obs | K2004 local − obs |
|---|---|---|---|
| R(48) 0-12 | 755.61 | **−0.13 MHz** | **+45.08 MHz** |
| P(80) 0-12 | 757.69 | −1.25 | +36.72 |
| R(100) 0-12 | 758.54 | −2.64 | +29.95 |
| P(126) 0-12 | 761.68 | −2.68 | +20.86 |

**The 49 MHz was the local model extrapolating, not the potentials.** Our potentials are the
accurate model at the band head; the 2004 local Dunham model — used blind, outside its stated
775–815 nm window — is the one that is wrong, by up to 45 MHz. This also resolves the open item in
`hannover-model-reproduction.md` ("0-12 at low J … differ by tens of MHz"): it is the local model's
extrapolation, and the potentials are fine there.

As `nir-anchor-data-2000.md` §5.6 notes, the Bodermann 2000 anchors could not test this: their only
(0-12) lines are J″ = 240 and 243, at 778–780 nm, the far high-J end of the band.

### 6.2 Three-way residuals, Liao's 31 components

Model − measurement in kHz, a1 unless noted, ordered by wavelength. Within a line the three
components agree to **4–137 kHz**, so one row per line tells the whole story.

| Line | λ (nm) | observed (MHz) | potentials | K2004 local | Liao 2010 as printed |
|---|---|---|---|---|---|
| R(48) 0-12 | 755.614 | 396 752 906.602 | **−126** | +45 076 | −23 087 |
| P(80) 0-12 | 757.689 | 395 666 265.779 | −1 245 | +36 720 | −25 787 |
| R(100) 0-12 | 758.543 | 395 220 756.623 | −2 641 | +29 953 | −28 126 |
| P(126) 0-12 | 761.685 | 393 590 751.059 | −2 678 | +20 855 | −31 890 |
| P(82) 0-13 | 769.321 | 389 683 863.049 | −2 646 | +6 203 | −33 058 |
| P(92) 0-13 | 770.067 | 389 306 731.793 | −2 816 | +5 374 | −34 366 |
| R(106) 0-13 | 770.522 | 389 076 805.543 | −3 516 | +3 966 | −36 484 |
| P(152) 0-13 | 776.319 | 386 171 334.557 | −1 170 | +587 | −45 305 |
| R(26) 0-14 (a10) | 778.210 | 385 233 438.060 | −2 148 | **+28** | −35 967 |
| P(176) 0-13 | 779.705 | 384 494 408.214 | +1 163 | −383 | −51 306 |
| R(78) 0-14 (a10) | 780.258 | 384 222 471.416 | −3 399 | **−14** | −40 551 |
| P(70) 0-14 (a10) | 780.293 | 384 205 192.670 | −2 811 | **+46** | −39 524 |
| R(141) 0-16 (a2) | 810.233 | 370 007 070.341 | −4 280 | **+63** | −74 905 |
| **all 31 components** | | | **mean −1.91, rms 2.39, max 4.28 MHz** | mean +14.38, rms 21.55, max 45.09 MHz | mean −36.08, rms 37.65, max 74.91 MHz |

### 6.3 Both NIR data sets together (Liao 2010 + Bodermann 2000)

| Subset | n | potentials | Knöckel 2004 local |
|---|---|---|---|
| v′ = 0, **775–815 nm** (the 2004 model's own window) | 34–35 | rms **2.31**, max 6.22 MHz | rms **0.23**, max 0.72 MHz |
| v′ = 0, 755–775 nm (outside it) | 21 | rms **2.46**, max 3.52 MHz | rms 26.18, max 45.09 MHz |
| v′ > 0 (bands 1-14, 2-15, 3-16) | 4 | rms **5.22**, max 7.05 MHz | no model exists |

Three conclusions [I]:

1. **Inside its window the local model is 10× better than the potentials** and should be used there.
2. **Outside it the local model is worthless and the potentials are unchanged in quality** — our
   accuracy at 755–775 nm is the same few MHz as at 775–815 nm.
3. **v′ > 0 is the weak spot**, at 5–7 MHz, and no published local model covers it.

Our hyperfine offsets are validated along the way: the measured a10 − a1 and a15 − a1 spacings of
the nine Liao lines agree with ours to **4–137 kHz**, and Bodermann's nine implied spacings to
15–85 kHz (`nir-anchor-data-2000.md` §4.5).

---

## 7. What was implemented

| File | Change |
|---|---|
| `src/i2spec/local_nir.py` | **new.** `LocalNIRModel` + `Domain` (v′, v″ range, J″ limit, wavelength window, claimed uncertainty, recommended J″) and `load_local_nir(name, trust_printed_digits=False)`. `transition`, `transition_MHz`, `covers`, `covers_wavelength`, `covers_wavenumber` |
| `src/i2spec/data/knockel2004_nir.json` | **new.** Knöckel 2004 Table 1, with uncertainties and the measured validation numbers of §6.3 |
| `src/i2spec/data/liao2010_nir.json` | **new.** Liao 2010 Table 3, with uncertainties, the domain, and the §3 defect recorded as `printed_digits_reproduce_the_fit: false` |
| `src/i2spec/lookup.py` | `uncertainty()` rebuilt below 15 000 cm⁻¹ (§7.1); new `NU_NIR_ANCHORS` (815 nm), `local_nir_covers()`; flag "beyond 755 nm" replaced by "beyond 815 nm" and "local NIR model" |
| `src/i2spec/__init__.py` | exports `LocalNIRModel`, `load_local_nir` |
| `data/observations/liao2010a/` | **new**, 31 rows (§4.2), validated by `load_dataset` |
| `tests/test_local_nir.py` | **new**, 7 tests (§8) |
| `tests/test_hannover2008.py` | the inline copy of Knöckel 2004 Table 1 now loads from the JSON, so there is one copy |
| `tests/test_lookup.py` | the NIR expectations follow the new bands |

### 7.1 `lookup.uncertainty`, before and after

| Range | Before | After | Basis for the new value |
|---|---|---|---|
| ≥ 15 000 cm⁻¹ (≤ 667 nm) | 3 MHz | 3 MHz | unchanged (BIPM) |
| 13 250–15 000 (667–755 nm) | 5 MHz | 5 MHz | unchanged (the 671 nm anchor) |
| 12 270–13 250 (**755–815 nm**), v′ = 0 | 50 MHz | **3 MHz** | 56 comb-referenced components (Liao 2010 + Bodermann 2000): rms 2.4 MHz, worst 6.2 MHz. The basis string also points to `i2spec.local_nir`, which does 0.23 MHz above 775 nm |
| 12 270–13 250, **v′ > 0** | 50 MHz | **8 MHz** | the 1-14 band sits 7.0 MHz off in Bodermann 2000, and no local model covers v′ > 0 |
| < 12 270 (> 815 nm) | 50 MHz | 50 MHz | unchanged, but the basis is now honest: nothing has ever been measured absolutely beyond 815 nm |

The old 50 MHz rested on the model-vs-model disagreement that §6.1 has now shown to be the local
model's fault. Note that `uncertainty()` describes **the position i2spec reports**, which always
comes from the potentials; where a local model does better, the basis string says so and the line
carries the flag `local NIR model`.

**Not done, deliberately:** switching `Catalog`/`lookup` to *use* the local model for positions. That
means a second position source inside the master line list and touching the CLI and TUI, which other
work owns. The hook is in place (`local_nir_covers`), and the 0.23 MHz is one call away for any
caller that wants it.

---

## 8. Tests and remaining work

`tests/test_local_nir.py` (all tolerances set just outside what the models achieve today, so drift
in either direction fails):

| Test | Asserts |
|---|---|
| `test_liao2010_is_refused_unless_the_caller_insists` | the 2010 set raises without `trust_printed_digits` |
| `test_domain` | v′ = 0 only, v″ 12–17, J″ ≤ 242, the two wavelength windows, branch validation |
| `test_knockel2004_beats_the_potentials_inside_its_own_window` | local < 0.30 MHz rms, potentials > 2.0 MHz, on the 34 anchors above 775 nm |
| `test_knockel2004_fails_where_liao_measured` | local > 20 MHz rms below 775 nm while the potentials stay < 3 MHz |
| `test_the_755_nm_band_head` | R(48) 0-12: potentials within 0.5 MHz, local model beyond 40 MHz |
| `test_potentials_against_the_liao_measurements` | rms < 2.6 MHz, max < 4.5 MHz over all 31 components |
| `test_liao2010_as_printed_does_not_reproduce_its_own_measurements` | the §3 defect, pinned |
| `test_model_metadata` | 19 parameters in each set, the two wavelength windows, the claimed 0.2 MHz, and that the refit moved S₂₂ to S₃₂ |

**For the Phase B refit.** These 31 components plus Bodermann's 32 are the only precise anchors
between 671 and 815 nm. Weight them at their true uncertainties (0.2 MHz for Liao — the cell shift,
not the 8–50 kHz statistics), give the `liao2010` group its own nuisance offset (the 114 kHz
pressure shift is exactly that), and the refit should reach the local model's 0.1 MHz across
755–815 nm while keeping v′ > 0 honest. If it does, the local layer becomes redundant — the
desirable outcome.

**Still open in this region:** 680–715 nm and 740–750 nm have no absolute anchor at all; v′ > 0 in
the NIR rests on four rows; and beyond 815 nm nothing has been measured. The unrounded 2010
parameters are not recoverable from the literature, and there is no supplement — but §6.3 shows we
do not need them: the 2004 local model plus our potentials already cover 755–815 nm at 0.23 and
2.4 MHz respectively, and the Phase B refit is the real fix.

---

## 9. Data set row

For the "Data sets" table of `docs/design/observations.md` (owned elsewhere):

| `liao2010a` | Liao 2010 (750–780 nm): 27 hyperfine components of nine (0-12)/(0-13) lines (Table 2), three (0-14) cell-calibration lines (Table 1), and the R(141) 0-16 a2 comb value from §2 | 31 | complete; Table 2's "P(100) 0-12" entered as R(100) 0-12 (see `meta.toml`) |

---

## 10. Addendum, 2026-09-19: the corrections are in the default

§7's "not done, deliberately" is done, in a different form. Rather than a second position source,
`i2spec2026b` carries measured **level corrections** on top of the potentials — polynomials in J(J+1)
per X level v″ = 11–17 and per B level v′ = 1–3, 5, 7, fitted to all 71 comb-referenced NIR lines
(`prototypes/nir_corrections.py`, `docs/design/parameter-sets.md`). In-sample 0.25 MHz over the 71
lines (Liao 0.34, Bodermann 2000 0.16), 0.38 MHz median leave-one-line-out; the local NIR module stays
as the published reference. The reason for the potentials' 2–6 MHz is now visible: a smooth, strongly
curved J-dependence within each band, and a B v′ = 1 rotational constant off by ~0.6 kHz.
