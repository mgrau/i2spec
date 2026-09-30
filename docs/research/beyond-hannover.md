# Beyond Hannover: where the model can still improve

*2026-09-18. Numbers from `prototypes/mlr_x.py`, `prototypes/hyperfine_stage2.py`,
`prototypes/hyperfine_chen2004.py` and a per-set error budget (below). "Published" is `hannover2008`
with BKT02/S06 hyperfine, i.e. what IodineSpec implements; "MLR X" swaps in `mlr_x_2026b`.*

## What "better than Hannover" would mean

IodineSpec is the reference because nothing else predicts the B–X spectrum to MHz. It has known limits,
and each is a place we can be measurably better:

| where | published model today | what the data could support | status |
|---|---|---|---|
| ¹²⁷I₂ line positions, v″ ≤ 17, v′ ≤ 44 | 2–5 MHz rms on comb-referenced sets (Hannover claims 1.5 MHz 1σ) | 0.02–0.3 MHz | **0.02–0.3 MHz in-sample, 0.1–0.8 held out** by the level corrections of `i2spec2026c/d` |
| X levels v″ = 48–54 | 7.8–21 cm⁻¹ out | 0.1 MHz | **7–19 MHz in the default** (`i2spec2026d`, the MLR X; the rest is X hyperfine) |
| X levels v″ = 18–47 | unmeasured; MLR predicts −21 cm⁻¹ at v″ = 42 | — | needs the Gerstenkorn atlases |
| B levels v′ = 45–62 | up to 2.26 cm⁻¹ out; R(98) 58-1 at 5.4 GHz | 0.005 cm⁻¹ (atlas) | **in the default** (`i2spec2026d`): 0.03 cm⁻¹ on the atlas, yoshiki 0.6 MHz, matsunaga 1.5 MHz in-sample |
| hyperfine, ¹²⁷I₂ v′ ≤ 53 | 25 kHz best, 330 kHz rms at 532 nm | 2–5 kHz | open: lever 3 |
| hyperfine, v′ = 55–70 | 63 MHz median, 470 MHz worst | ~5 kHz | per-v′ tables validated, not built |
| isotopologues | ¹²⁹I₂ 7–14 MHz low, ¹²⁷I¹²⁹I 4.5 MHz high | 1 MHz | open: lever 4 |
| per-line uncertainties | none published | — | rule-based only; lever 5 |

## Measured against the real thing (2026-09-19)

IodineSpec5 was run on every observation line (`docs/research/iodinespec5.md`).
It reproduces `hannover2008` to 0.48 MHz rms for ¹²⁷I₂, uses the Salumbides 2006 hyperfine formulae
throughout, and carries unpublished isotopologue corrections and a NIR local Dunham model. Against it
the default set is better on six sets (`bipm2003a` 2.2×, `reinhardt2006a` 2.1×, `kobayashi2016a`
1.8×, `bipm2003c`, `nishiyama2024a`, `bipm2003b`), equal on nine, worse on six: the four NIR sets
(its in-sample local model, 0.09–2.1 MHz against our 2.3–3.4), the 612 nm isotopologue set (1.5 vs
2.8) and the v′ = 43 hyperfine rows (0.16 vs 0.26). Levers 1 and 4 below are what those six need. **Answered the same day:** `i2spec2026c` carries measured level corrections for every comb-referenced line with v′ ≤ 43, v″ ≤ 17 (`docs/design/parameter-sets.md`); against IodineSpec5 it is better on 16 sets, equal on 3, worse on 2 (one ¹²⁹I₂ line; Liao's two bands, 0.34 vs 0.14 in-sample).

## The error budget

For each data set: the rms of the stated uncertainties (the floor), and the model's rms residual, split
into rows that constrain the potentials (frequencies, inter-line intervals) and rows that constrain only
the hyperfine Hamiltonian (intra-line intervals). MHz throughout.

| set | rows | potential rows: floor / published / MLR X | hyperfine rows: floor / published | what limits it |
|---|---|---|---|---|
| `bipm2012a` | 20 + 309 | 0.005 / 3.82 / 4.04 | 0.004 / 0.33 | B levels v′ 32–37 (±6 MHz, smooth in v′); C_B |
| `bipm2003b` | 16 + 14 | 0.014 / 0.57 / 1.45 | 0.011 / 0.04 | potentials |
| `bipm2003c` | 1 + 10 | 0.20 / 1.90 / 4.00 | 0.06 / 0.10 | potentials; a3/a4 data conflict |
| `bipm2003e` | 4 + 14 | 0.14 / 0.44 / 0.77 | 0.41 / 0.16 | at the floor |
| `bipm2005a` | 3 + 54 | 2.9 / 23 674 / 23 675 (MLR X+B: 638) | 0.8 / 21.3 | v′ = 58: B potential and B hyperfine |
| `bodermann1998b` | 3 + 1 | 0.015 / 3.92 / 1.89 | 0.015 / 0.04 | potentials, v″ 16–17 |
| `bodermann2000a` | 32 | 0.10 / 3.13 / 2.47 | — | potentials, v″ 12–16 |
| `liao2010a` | 31 | 0.20 / 2.30 / 1.35 | — | potentials, v″ 12–16 |
| `reinhardt2007a` | 3 + 4 | 0.30 / 4.89 / 5.56 | 0.04 / 0.01 | potentials |
| `reinhardt2006a` | 1 + 56 | 0.07 / 1.48 / 3.71 | 0.026 / 0.12 | C_B (per line reaches 23 kHz) |
| `cornish2000a` | 2 | 0.85 / 5.61 / 4.06 | — | potentials |
| `nesterenko2019` | 18 | 0.10 / 544 836 / 20.6 | — | was X; now hyperfine at v″ 53–54 |
| `matyugin2012` | 18 | 0.035 / 233 738 / 15.9 | — | was X; now hyperfine at v″ 48 |
| `bipm2003a` | 80 + 78 (117 iso) | 1.22 / 3.82 / 3.22 | 0.72 / 0.72 | isotopologue BO corrections |
| `bipm2003d` | 70 + 20 (46 iso) | 1.06 / 8.69 / 8.52 | 0.19 / 0.19 | isotopologue BO corrections |
| `velchev1998a` | 115 | 2.16 / 1.95 / 2.70 | — | at the floor |
| `xu2000a` | 473 | 1.00 / 0.94 / 1.24 | — | at the floor |
| `dube2004a`, `huet2013a`, `morinaga1989a` | 96 | — | 0.3–1.5 / 0.3–1.0 | at the floor |

Three things stand out.

1. **Every comb-referenced position set sits 10–250× above its floor.** `bodermann1998b` is measured to
   15 kHz and reproduced to 3.9 MHz; `bodermann2000a` 0.1 → 3.1; `liao2010a` 0.2 → 2.3; the BIPM
   532 nm inter-line intervals 5 kHz → 3.8. Hannover's own claim, 1.5 MHz 1σ, sits between. That is
   the potentials, and it is where a joint fit has the most to gain.
2. **The hyperfine rows sit 4–80× above theirs**, and nothing else touches them: `bipm2012a`'s 309
   intra-line intervals are 4 kHz data reproduced to 330 kHz. This is a separate problem with its own
   data (656 rows plus Chen 2004) and its own answer (Stage 2 found per-line parameters reach the noise).
3. **The MLR potentials removed the two catastrophic failures without spending anything.** X went from
   544 GHz to 21 MHz at v″ = 53–54 with De and C₆–C₁₀ untouched; the B MLR, fitted to atlas bands only,
   takes the *held-out* R(98) 58-1 from 5431 to 148 MHz. Both were a functional form extrapolating past
   its data, not physics or measurement. Neither is at its floor yet: 16–21 MHz at v″ = 48–54 is the
   hyperfine model's own error there (26–45 MHz spread among one line's components), and the B MLR is
   204 MHz rms below v′ = 45 until the joint fit holds both ends.

## The levers, ranked by what the data say they are worth

### 1. Finish the joint MLR fit — positions from 2–5 MHz to ~1 MHz, and both high-v ends

*Gain:* every position set above; 5.4 GHz → ≲150 MHz at v′ = 58; the 2.2 cm⁻¹ atlas failure gone.
*Evidence:* the separate fits already reach these numbers at each end (`docs/design/mlr-x.md`).
*What it takes:* the joint fit exists and its Jacobian is verified; what it lacks is weighting. The data
run from 5 kHz to 150 MHz in stated precision while the model's own error is tens of MHz, so `soft_l1`
at 3 MHz writes off everything it cannot nearly fit. A model-error floor in quadrature (as
`hyperfine_fit.MODEL_FLOOR`) and per-group calibration offsets with the ±3 MHz bound from Stage 1 —
Stage 1 found them all under 0.2 MHz, so calibration is not what limits the fit — plus multi-start
against the local minima seen in every MLR run so far. Days, not weeks.
*Ceiling:* the hyperfine offsets frozen in the fit are 25–330 kHz wrong (lever 3), so sub-100 kHz
positions need lever 3 first.

**Status 2026-09-18:** the joint least squares stalls in a stiff landscape, but alternating single-state
level fits work: B fitted to the 787 measured levels with X held reaches 1.7 MHz at v′ ≤ 30 and
17 MHz at v′ = 45–50, and the X+B pair is within 2× of the published model on the precision sets and
10²–10⁴ better beyond them (`docs/design/mlr-x.md`). Done since: intervals folded in, one alternation (`2026c` pair): 2–4× of the published model on the
precision sets, 10²–10⁴ better beyond. Not the default until it is better on every set — more B terms and
another alternation are the next cycle.

### 2. Transcribe the kHz-level data we already hold — the cheapest gain of all

*Gain:* about 380 comb-referenced hyperfine components at 5–7 kHz, none of which the Hannover fit had:
`yoshiki2023a` (103, v′ = 44–45), `matsunaga2024a` (97, v′ = 45–50), `kobayashi2016a` (81, 578 nm),
`nishiyama2024a` (520 nm, v′ = 39), plus `sansonetti1997a` (102 components, 560–656 nm at 1 MHz, and
the labelling authority for the whole component convention). v′ = 44–50 is exactly where the B
potential starts to go wrong and where the atlas is our only constraint; these replace 0.005 cm⁻¹ band
shifts with 5 kHz frequencies. Yoshiki also fits eqQ′ to 4 kHz, which is hyperfine data at several J′
within v′ = 44–45, the thing Stage 2 said it needed.
*What it takes:* the PDFs are in hand. The same script-plus-page-image transcription as Chen 2004 and
the Velchev/Xu sets, validated line by line against the model.

**Done 2026-09-18:** all five sets are in `data/observations/`, 426 rows (`docs/design/observations.md`).
Per-line four-parameter hyperfine fits reproduce every new line at the papers' own 1–3 kHz, and the two
514 nm papers agree on the v′ = 45 B-state offset (+300 and +306 MHz). Three misprinted classifications
in Sansonetti's table were caught and corrected by the model.

### 3. A hyperfine model that predicts — 330 kHz to ~30 kHz, and 63 MHz to kHz above v′ = 53

*Gain:* the largest in ratio terms, and the only lever that improves the BIPM lines themselves.
*Evidence:* per-line B parameters reach the measurement noise on 25 of 29 precise lines; the error is
mostly C_B, smooth in J′ within a v′ and stepping between v′ (Chen 2004). Above v′ = 53 interpolating
measured parameters across v′ beats the frozen formula 20×.
*Two routes, in order:*
- **Tables.** Per-v′ parameter tables from Chen 2004 and the Stage 2 per-line fits, a line in J′(J′+1)
  under BKT02 and a quadratic under S06 and above, with the formulae as fallback at unmeasured v′.
  **Done 2026-09-18** (`i2spec.hfs_table`, default on): held-out lines at a measured v′ go from 270 to
  151 kHz overall and from 25–190 to 1–10 kHz where several J′ exist; R(98) 58-1 from 42 MHz to ~4.
- **Physics.** The Phase C plan (`03a-theory-physics.md` §3.2): eqQ(R) as an R-dependent function with
  ⟨vJ|eqQ(R)|vJ⟩ expectation values, and C, d, δ as explicit second-order sums over a few model
  perturber states with fitted strengths. The v′-to-v′ steps Chen's data show are what energy
  denominators to specific perturbing levels produce, so this is the form that could predict an
  unmeasured (v′, J′). It needs lever 2's data to fit, and it is the only route to a global model.

### 4. Isotopologues — 4–14 MHz to ~1 MHz

*Gain:* ¹²⁹I₂ and ¹²⁷I¹²⁹I lines, `bipm2003d` from 8.7 MHz to its 1 MHz floor.
*Evidence:* the code follows Salumbides 2008 as printed and reproduces the paper's own test shift, but
every BIPM isotope shift and Salumbides 2006's own spectra come out 5–14 MHz off, in the same direction.
*What it takes:* refit the Born–Oppenheimer functions (B has α with 6 terms and V_ad with 4; X has
none in the published set) to the 163 isotopologue rows, or replace them with Le Roy/Watson per-atom
forms. With two isotopes only, the adiabatic and field-shift terms are degenerate, so this needs priors
and will not separate physics from convention. A modest, bounded job once lever 1 is done.

**Done 2026-09-18, in the smallest form that the data support:** one constant in V_ad(B), +0.02154 cm⁻¹,
takes the 124 isotopologue rows from 7.06 to 2.36 MHz rms (floor 1.1–1.2). Seven independent numbers
cannot fix a shape, and fits that try are degenerate. `docs/design/parameter-sets.md`.

### 5. Uncertainties — something IodineSpec does not have at all

*Gain:* a per-line uncertainty with a stated basis, instead of the rule in `lookup.uncertainty` and
instead of nothing.
*What it takes:* M6 as planned — the joint fit's covariance through Hellmann–Feynman sensitivities,
a discrepancy term for the model error the budget above measures, and blocked cross-validation for
the coverage claim. Depends on lever 1.

### 6. Measure the gap — the test of the MLR prediction

*Gain:* v″ = 18–47 has never been measured. The MLR says the published potential is 21 cm⁻¹ high at
v″ = 42, and two fits that agree there to 0.002 cm⁻¹ diverge by 118 cm⁻¹ at v″ = 58. Data in the gap
turn a prediction into a validated potential, or kill it.
*What it takes:* the Gerstenkorn atlases at 11 000–15 600 cm⁻¹, or any fluorescence data reaching v″ > 20.
*(Since answered for v″ = 18–25 by the Orsay atlas: `orsay-atlas-11000-14000.md`.)*

### 7. Intensities — already at the level of the reference disagreement

The relative intensities hold the implied I₂ column to ±7 % across 2300 cm⁻¹ of the atlas, and the
absolute scale sits −0.2 % from Tellinghuisen and +3 % from Spietz, which disagree with each other by
that 3 %. A refit of μ(R) to the atlas could narrow the 7 % (much of it is the 190 °C segment's
temperature), but there is no external reference to beat below 3 %. Lower priority.

### Already better, for the record

- Level energies converged to under 1 kHz (B-spline solver with knots at the joins); the sinc-DVR
  was 0.24 MHz off at 532 nm, and the line list carried that until 17 September.
- Line positions and hot-band coverage to v″ = 54 and v′ = 62 in the lookup, with honest flags.
- Physical constants pinned to CODATA 2022 (a 40 kHz effect); the mass convention is worth ±1 MHz in
  the NIR when comparing with IodineSpec, and is documented.

## The order that makes sense

2 → 1 → 3 (tables, then physics) → 5 → 4, with 6 running in parallel as the library delivers. Lever 2
is a week of transcription that feeds everything after it; lever 1 is the fit we already have with its
weighting fixed; lever 3's table route is a day and pays immediately at 500–510 nm; levers 5 and 4
need the joint fit to exist. Lever 7 waits for a better reference.
