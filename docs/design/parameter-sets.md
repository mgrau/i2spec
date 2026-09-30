# Parameter sets

`src/i2spec/data/<name>.json`, loaded by `potentials.load_potentials(name)`; `constants.DEFAULT_PARAMETERS`
names the one the model uses unless told otherwise.

| name | what it is | status |
|---|---|---|
| `hannover2008` | Salumbides et al. 2008, Table 1, as printed: X and B X-representation potentials, the B-state Born–Oppenheimer corrections | the published model; kept unchanged for comparison |
| `i2spec2026a` | `hannover2008` plus the changes i2spec has validated, listed in the file's `changes` | **the default** since 2026-09-18 |

Two things the model applies on top of whichever set is loaded, and which are not parameter-set
entries: the measured B-state hyperfine corrections (`i2spec.hfs_table`, on by default for ¹²⁷I₂,
`docs/design/hyperfine-fit.md`) and, when asked for, the trial MLR potentials in `data/potentials/`
(`docs/design/mlr-x.md`), which are not the default because their joint fit is unfinished.

## i2spec2026a

### B-state V_ad constant: +0.02154 cm⁻¹ (2026-09-18)

Every ¹²⁹I₂ line in the BIPM 612 and 633 nm tables sat 7–14 MHz above the published model, and
¹²⁷I¹²⁹I 4.6 MHz below ¹²⁹I₂, each with sub-MHz scatter inside the line: a shift of the line, not a
hyperfine error. The published set carries all Born–Oppenheimer terms in the B state — X has none —
as α(R), six terms, in the centrifugal part and V_ad(R), four terms, entering the levels as
(1 − μ₁₂₇/μ) (2Rm/(R+Rm))⁵ Σ vad_i x^i. A constant in V_ad therefore shifts every level of an
isotopologue by (1 − μ₁₂₇/μ) vad₀: 0.0155 vad₀ for ¹²⁹I₂ and 0.0085 vad₀ for ¹²⁷I¹²⁹I. One number
moves the two species in the ratio the data show, and it is the same thing as an electronic isotope
shift of the B–X origin.

`prototypes/bo_refit.py` fits it to the 124 isotopologue rows that constrain the potentials (the
intra-line hyperfine intervals are left out), with the hyperfine offsets frozen:

| fitted | rms on the 124 rows |
|---|---|
| published | 7.06 MHz |
| vad₀ + 0.02154 | **2.36 MHz** |
| vad₀, vad₁ | 2.35–2.48, depending on the start |
| vad₀, vad₁, α₀ | 1.39, with α₀ ≈ 0 and run-dependent vad |

Per line after the constant: P(33) 6-3 +2.3, P(54) 8-4 +1.5, P(110) 10-2 −2.0, R(113) 14-4 +5.1,
¹²⁷I¹²⁹I P(33) 6-3 (relative to ¹²⁹I₂) +0.1, the two intervals within ¹²⁹I₂ +1.6 and +0.2 MHz. The
data floor is 1.1–1.2 MHz. Seven independent numbers cannot determine a shape for V_ad, and the fits
that try find degenerate, run-dependent answers, so only the constant is adopted. What the shape
would change is the prediction for isotopologue lines at v′ far from the measured 6–14: about
±10 MHz between the candidate shapes. R(113) 14-4 remaining at +5.1 MHz says the true correction does
have some v′ or J′ dependence.

`tests/test_bipm_other.py::test_isotope_shifts` holds this, and its neighbour keeps the record that
`hannover2008` as printed misses the same shifts by more than 7 MHz.

## i2spec2026b

`i2spec2026a` plus measured level corrections in the near infrared (`level_corrections_2026a.json`,
named by the set's `level_corrections` key and applied by `RovibronicModel.levels` to ¹²⁷I₂ only).

### NIR level corrections: X v″ = 11–17, B v′ = 1–3, 5, 7 (2026-09-19)

The comb-referenced lines at 755–815 nm (`liao2010a`, `bodermann2000a`, `bodermann1998b`,
`reinhardt2007a`, `cornish2000a`: 71 lines measured to 0.02–0.3 MHz) sit 3.0 MHz rms from the
potentials, and IodineSpec5 answers exactly this with a local Dunham model fitted to the same data
(`docs/research/iodinespec5.md`). The residual is smooth in J within each band — 0-13 runs from +3.3 MHz
at J″ = 106 through 0 at 164 to −2.0 at 188 — and the sets agree with each other to 0.1 MHz where they
overlap, so the error is the potentials', not the data's.

The correction is a polynomial in y = J(J+1)/10⁴ per level, δE(state, v, J) = Σ cₖ yᵏ, degree 2–3 for
v″ = 12–14 (the well-covered bands), 1 for v″ = 15–16, 0 for v″ = 11 and 17 and for v′ = 2, 3, 5, 7,
and 1 for v′ = 1 (the 1-14 band spans 12 MHz over J″ = 42–148: the B v′ = 1 rotational constant is
off by ~0.6 kHz, and nothing but these lines constrains it). δB(v′ = 0) = 0 by convention, since almost
every line comes from v′ = 0. It is linear in the coefficients, so `prototypes/nir_corrections.py`
fits it as a robust (Huber) ridge regression on the residuals, floor 0.3 MHz in quadrature, prior 5 MHz.

| set | before | in-sample | leave-one-line-out, median |
|---|---|---|---|
| `liao2010a` (31) | 2.30 | 0.34 | 0.59 |
| `bodermann2000a` (32) | 3.13 | 0.16 | 0.16 |
| `bodermann1998b` (3) | 3.92 | 0.02 | — |
| `reinhardt2007a` (3) | 4.89 | 0.33 | 1.2 |
| `cornish2000a` (2) | 5.61 | 0.02 | — |
| all 71 | 3.04 | **0.25** | **0.38** |

Leave-one-*set*-out does not work, and says why: Liao measured the low and middle J of the 0-12 and
0-13 bands and Bodermann the high J, so each set is the other's extrapolation. The corrections are an
interpolation inside the measured J range, and `lookup.uncertainty` gives 0.5 MHz only within 15 in J
of it (`level_corrections.J_MARGIN`); outside, the old 3 and 8 MHz stand. `tests/test_level_corrections.py`.

## i2spec2026c

`i2spec2026b`'s level corrections extended from the near infrared to every comb-referenced ¹²⁷I₂ line
with v′ ≤ 43, v″ ≤ 17 (`level_corrections_2026b.json`, `prototypes/level_corrections_fit.py`).

### Level corrections across the visible (2026-09-19)

The residual structure of the default against the precise sets is level-shaped, not scatter:
`nishiyama2024a` is a pure +5.57 MHz offset of B v′ = 39 (scatter 0.11), `bipm2012a` a set of per-v′
offsets at v′ = 32–37 (−2 to +4 MHz, the "drift" the uncertainty rule carried as 5 MHz), `bipm2005a`
a −1.4 MHz offset of B v′ = 43, while the Doppler-limited sets (`xu2000a`, `velchev1998a`,
`sansonetti1997a`) are scatter at their own floor. The NIR fit's form applies unchanged: a polynomial in
J(J+1)/10⁴ per level, its degree set by the J coverage of the precise lines (σ ≤ 0.3 MHz) that reach it,
a robust ridge regression on the residuals against `i2spec2026a`, floor 0.3 MHz, prior 5 MHz. Only
precise lines create a level's parameters; the Doppler-limited lines weigh in. δX = 0 below v″ = 11
and δB(v′ = 0) = 0 remain the conventions. 839 position rows, 605 touching one of 32 levels, 46
coefficients. Where the potentials are off by GHz with a J-shape one polynomial cannot follow
(B v′ ≥ 44, X v″ ≥ 18) the fit is not attempted: that region is the MLR pair's (`mlr-x.md`).

| set | before | in-sample | leave-one-line-out rms / median |
|---|---|---|---|
| `bipm2012a` (532 nm, v′ 32–37) | 4.52 | 0.25 | 0.77 / 0.30 |
| `nishiyama2024a` (v′ 39) | 5.57 | 0.11 | 0.15 / 0.14 |
| `bipm2003b` (v′ 26–28) | 0.57 | 0.02 | 0.75 / 0.73 |
| `bipm2003d` (¹²⁷I₂ rows, v′ 9–15) | 2.34 | 0.24 | 1.51 / 0.15 |
| `kobayashi2016a` (v′ 16–18) | 2.14 | 0.62 | 1.12 / 0.81 |
| `liao2010a`, `bodermann2000a` | 2.30, 3.13 | 0.34, 0.16 | 0.68 / 0.59, 4.99 / 0.15 |
| `xu2000a`, `velchev1998a`, `sansonetti1997a` | 0.94, 2.22, 1.67 | 0.85, 2.05, 1.19 | at their floors |
| all 605 | 1.88 | 1.08 | 1.70 / 0.58 |

Each level's held-out rms over the precise lines that reach it is stored (`held_out_MHz`) and is what
`lookup.uncertainty` reports for a corrected line — 0.06 MHz for B v′ = 34, 0.15 for v′ = 39, 0.68 for
v′ = 32 — with 1 MHz assumed for a level fixed by a single line (v′ = 3, 5, 10, 22, 26, 43; X v″ = 11,
17), where nothing validates it. Against IodineSpec5 this set is better on 16 of the 21 sets inside its
range, equal on 3, worse on 2 (`docs/research/iodinespec5.md`).

## i2spec2026d

`i2spec2026c` plus the extended-range potentials: the MLR pair `mlr_x_2026c` / `mlr_b_2026c`
(`docs/design/mlr-x.md`) supplies every X level from v″ = 18 and every B level from v′ = 44 up
(`"extended"` key; `RovibronicModel.levels` splices the two solvers), and the level corrections are
refitted on that hybrid over the whole range (`level_corrections_2026c.json`).

### Why a hybrid, and not the pair everywhere (2026-09-19)

Inside the published range the two bases are equivalent once corrected — leave-one-line-out 1.68 MHz
(pair) against 1.70 (published) over the same 605 rows, the pair better on `bipm2003d`, `sansonetti`,
`velchev` and worse on `bipm2003b` (3.55 vs 0.75) and `cornish` — so neither passes the ship rule
against the other there, and the published curves stay where they are validated. Beyond it there is
no contest: the published X is +7.8 cm⁻¹ at v″ = 48 and −21 cm⁻¹ at v″ = 54, the published B 2.2 cm⁻¹
against the atlas and 25 GHz on `matsunaga2024a`, while the pair tracks Martin 1986's X levels to
0.1–0.3 GHz through the unmeasured gap (`docs/research/iodinespec5.md`) and is 12 MHz on Matsunaga
before its corrections. The seam is continuous: X v″ = 18 moves 0.001 cm⁻¹, B v′ = 44 0.001 cm⁻¹.

The corrections then reach the extended region too. With the prior widened where the potentials are
known to be off (100 MHz for B v′ = 44–50 and X v″ ≥ 18, 10⁵ MHz above v′ = 50) and any line measured
to 5 MHz allowed to fix a level above v′ = 50:

| set | before (hybrid) | in-sample | held-out rms / median |
|---|---|---|---|
| `yoshiki2023a` (v′ 44–45) | 3.26 | 0.58 | 0.74 / 0.71 |
| `matsunaga2024a` (v′ 45–50) | 12.5 | 1.51 | 12.2 / 5.2 (one or two lines per level) |
| `bipm2005a` R(98) 58-1 | 3266 | 0.00 | — (one line; the 1g region) |
| `matyugin2012` (v″ 48) | 6.5 | 6.5 | X hyperfine, not the levels |
| `nesterenko2019` (v″ 53–54) | 18.7 | 18.7 | X hyperfine, not the levels |

Inside the range every number of `i2spec2026c` is reproduced (`bipm2012a` 0.25 / 0.77, `nishiyama2024a`
0.11 / 0.15, `liao2010a` 0.34 / 0.68). B v′ = 44 and 45 are validated at 0.7 MHz held out; v′ = 46–50
and 58 are fixed by one line each. `lookup.uncertainty` follows: corrected lines carry their held-out
figure, and the extended region without them 15 MHz (v′ 44–50), 2 GHz (v′ > 50, the 1g crossing),
0.3 GHz (v″ 18–28, checked against Martin), 1 GHz (v″ 29–47) and 20 MHz (v″ 48–54). The line list
matches its DVR intensity levels to the position model by index below each asymptote, since the two
now differ by cm⁻¹ at high v″, and drops the few quasi-bound lines at the B asymptote it cannot match.

One rule that the fit needed and that is worth remembering: its base must never carry corrections.
`BarePredictor` strips them whatever the set names, after two runs that silently fitted corrections on
top of corrections and wrote zeros.

## i2spec2026e

`i2spec2026d` with one addition, `level_corrections_2026d.json`: the same level corrections plus a
`band` section — a correction to the *lines* of v′ = 0 → v″ = 12–17, applied within 15 in J″ of the
data at each v″ and nowhere else (`i2spec.level_corrections.LevelCorrections.band_shift`, used by
`RovibronicModel.transition` and `intensity.master_line_list`).

### NIR band correction (2026-09-23)

What the NIR lines needed after the level corrections was a J′-dependence of B v′ = 0 (R-high / P-low
at the same J″). As a level correction it moved 0-9/0-10 by up to +47 MHz and the v′ ≥ 1 bands to
v″ = 11–17 by up to −130 MHz, unmeasured; as a band correction it moves nothing outside the bands.
`liao2010a` 0.344 → 0.145 MHz raw (0.074 with Liao's 114 kHz pressure shift), 0.22 held out; the 60
precise NIR rows 0.274 → 0.248 held out, worst 0.72 MHz. `prototypes/nir_band_fit.py --degree=3,2,3
--floor=0.05`, chosen by a rule fixed in advance; the full account, including why the more flexible
fits that beat IodineSpec5 in-sample were rejected, is in `docs/research/iodinespec5.md`.

## i2spec2026f

`i2spec2026e` refitted with Bodermann's 1998 thesis data (`data/observations/bodermann1998c`, 54 rows:
the unpublished beat frequencies behind the Hannover NIR local model, and the Metrologia 1998
wavelength-comparison lines). `level_corrections_2026e.json` = level corrections refitted with them,
plus a band section of degree (5, 3, 4). `liao2010a` 0.119 MHz raw against IodineSpec5's 0.138 (0.027
with Liao's pressure shift, 0.053 held out); all 108 precise NIR rows 0.132 MHz held out. Details and the
selection grid: `docs/research/iodinespec5.md`, "With Bodermann's thesis".

## i2spec2026g

`i2spec2026f` plus X-state level corrections for v″ = 18–25 from the Orsay atlas part I
(`level_corrections_2026f.json`, its `orsay` section; `data/atlas_lines/orsay1982_part1*.csv`). 3 700 of
the atlas's 4 826 lines assigned (4 % expected by chance), fitted with the atlas's own pure-scale
calibration (+24 ppb) and a linear J(J+1) correction per level: −13 MHz (v″ = 18) to −45 MHz (v″ = 24) at
J″ = 100, level uncertainty 15–52 MHz where `lookup` quoted 300. No measured data set moves. Details:
`docs/research/orsay-atlas-11000-14000.md`, "Assignment and the X levels".

## i2spec2026h

`i2spec2026g` with the extended-range X refitted to the Orsay-atlas levels: `mlr_x_2026d`
(`prototypes/mlr_x_atlas.py`), and the level corrections re-expressed on it (`level_corrections_2026g`), so
every atlas-covered level is unchanged to < 1 MHz. Leave-one-level-out over v″ = 18–25 beats the previous
X on every level; v″ = 26–47 move by up to −100 MHz. Details: `docs/research/orsay-atlas-11000-14000.md`,
"The potential refitted".

## i2spec2026i

`i2spec2026h` with `level_corrections_2026h.json`: beyond each level's measured J range ± 15, the quadratic
and cubic terms of a level correction are held at their value at the edge (`clamp_J`), while the constant
and J(J+1) terms keep extrapolating. Extrapolated curvature had moved lines by tens of MHz (1-13 +57,
1-14 −35 MHz at J″ = 300); holding the whole polynomial was tried first and fails on `velchev1998a`'s
v′ = 16 lines beyond the fitted J. No measured data set moves.

## i2spec2026j

`i2spec2026i` with the B section of the level corrections refitted with seven more comb-referenced data
sets (`level_corrections_2026i.json`; the X, NIR band and Orsay sections are unchanged): `ikeda2022a`,
`cheng2019a`, `shie2013a`, `hsiao2013a`, `zhang2009a`, `yang2011a` and `edwards1999a`, 174 rows. The fit
is the one behind `level_corrections_2026c`–`h` (`prototypes/level_corrections_fit.py --base=i2spec2026d
--extended`), with its degree rule unchanged. Run without the new sets, it reproduces the B section of
2026h to 10⁻⁶ MHz, so only the new data change it. v′ > 50 stays as in 2026h (only v′ = 58 corrected).

rms of observed − model (MHz), absolute frequencies:

| set | 2026i | 2026j |
|---|---|---|
| `shie2013a` (v′ = 30) | 4.69 | 0.02 |
| `hsiao2013a` (v′ = 24) | 4.31 | 0.11 |
| `zhang2009a` (v′ = 20, 23) | 3.35 | 0.07 |
| `cheng2019a` (v′ = 31) | 0.91 | 0.003 |
| `ikeda2022a` (v′ = 45, 46, 48) | 15.5 | 1.65 |
| `yang2011a` (v′ = 22, 24) | 2.13 | 4.93 |
| `sansonetti1997a` (all rows; creates no parameters) | 1.28 | 1.12 |
| `matsunaga2024a` | 1.51 | 1.70 |

Every other set is unchanged to 0.01 MHz. Before the refit the new visible absolutes sat 1–5 MHz from the
model, against a quoted 3 MHz; P(91) 48-0 sat 26 MHz off, against 15.

**B v′ = 24 is inconsistent at one offset.** Hsiao's P(28) 24-0 (J′ = 27) wants −4.3 MHz and Yang's
R(130) 24-1 (J′ = 131) +2.7 MHz. Two J values allow only a constant under the degree rule, and the robust
weighting follows Hsiao's three components, leaving Yang's line 7 MHz off; the level's held-out figure,
and so the quoted uncertainty of its lines, is 7 MHz. v′ = 20 shows a J-dependence of the same sign:
Sansonetti's lines are −3.6 MHz at J′ = 33–43 (confirming Zhang's −2.9 at J′ = 35) and −1.4 at J′ = 85.
The pattern is not one missing rotational term: at v′ = 15 Sansonetti's lines at J′ = 115–123 sit
−3.7 to −4.9 MHz, the opposite sign at high J, and most v′ = 14–18 lines scatter about zero. The
residual J-dependence changes from level to level, which is what a potential slightly off in shape (rather
than in a rotational constant) would do. More lines at high J, especially v′ = 15–25, would settle it.

## i2spec2026k

`i2spec2026j` with the B state taken to the dissociation limit, from the Orsay atlas Partie IV
(Gerstenkorn & Luc 1983; `docs/research/orsay-atlas-19700-20035.md`):

- **`mlr_b_2026d`**, the B MLR refitted with 2 009 unblended atlas lines (v′ = 51–79) added to
  `mlr_b_2026c`'s targets, now with the comb lines at every v′ (`prototypes/mlr_b_dissociation.py`):
  20 β, and a constant `q_far` = −6.6 × 10⁻⁴ added to the published α(R) beyond 5 Å. Without it every
  band showed the same J trend (+500 MHz at J′ = 0 to −500 at J′ = 80): α(R) was fitted in the well and
  its extrapolation sets the rotational energy of the last levels. Levels to v′ = 86 are bound at J = 0;
  a 40 Å box with a graded mesh (`BSplineSolver(mesh=...)`, `ADAPTIVE_GRIDS["B"]` to v′ = 90) holds them,
  equal to a uniform 30 Å box to 10⁻⁴ MHz.
- **`level_corrections_2026j`**: B v′ = 44–50 refitted on the new potential (the fit of 2026j, B v′ = 44–50
  only; every other section as in `level_corrections_2026i`), and B v′ = 51–79 corrected level by level
  from the atlas and the comb lines there (`prototypes/orsay4_fit.py`, section `orsay4`). The atlas offset
  is fitted with them, +119.6 MHz, against +82, +150 and +138 MHz on the three atlas lines also measured
  against a comb. Held-out 20–70 MHz per level at v′ = 51–72, 75–250 at v′ = 73–79.

Against the atlas (unblended lines, offset applied), median obs − model per band: 3–390 GHz in 2026j
above v′ = 58, and no bound level above v′ = 75; now within ±15 MHz for every band v′ = 51–78, with
13–80 MHz scatter to v′ = 72 and 70–150 MHz above.

rms of observed − model (MHz), absolute frequencies:

| set | 2026j | 2026k |
|---|---|---|
| `goncharov2007a` (R(26) 62-0) | 469 | 0.34 |
| `sakamoto2024a` (v′ = 52, 53) | 260 | 0.18 |
| `ikeda2022a` (v′ = 45, 46, 48) | 1.65 | 0.61 |
| `matsunaga2024a` (v′ = 45–50) | 1.70 | 0.99 |
| `yoshiki2023a` (v′ = 44) | 0.57 | 0.20 |

Every other set is unchanged to 0.02 MHz. `lookup.uncertainty` quotes each atlas-covered level's
held-out figure (flag `v′ 51-79: atlas-measured`) and 1 GHz for v′ > 50 outside that coverage, in place
of the flat 2 GHz of 2026j, which the atlas showed to be wrong by up to 200× above v′ ≈ 58.
