# Our own potentials: a direct fit from the atlas constants to our data

*Started 2026-09-20. `src/i2spec/rkr.py`, `prototypes/mlr_from_rkr.py`, `prototypes/global_fit.py`,
`prototypes/global_fit_evaluate.py`. This is Phase B/C proper: not the published potentials with
corrections on top, but potentials fitted here, to data we assembled, from a starting point that owes
nothing to Hannover.*

## What makes this possible now

Hannover fitted ~1950 lines, most of them atlas lines at 0.002–0.005 cm⁻¹. Until this week we had
1013 precision positions — a thousand times more accurate, but touching a small part of the (v′, v″)
plane, which is why every attempt to refit the potentials wandered in the gaps. The atlas extraction
(`docs/research/atlas-line-positions.md`) changed that: **67 729 line centres at 10–25 MHz**, spanning
v′ = 2–58 and v″ = 0–8, an order of magnitude better than the atlas lines Hannover had and thirty times
more numerous.

## The chain, and what each link is checked against

1. **Gerstenkorn & Luc 1985 Dunham constants** (`i2spec.dunham`), transcribed from the HAL copy and
   checked against their own Table V and against the published potentials.
2. **RKR inversion** (`i2spec.rkr`) → classical turning points. Checked against the published X
   potential: **0.05 mÅ at every v″ the constants cover**. The B state agrees to 0.6 mÅ to v′ = 30 and
   parts by 4.5 mÅ at v′ = 60, where that curve is independently known to fail.
3. **MLR fitted to the RKR curve** (`prototypes/mlr_from_rkr.py`), long-range coefficients held at
   their published measurements: Bacis, Cerny & Martin 1986 for X (De, C₆, C₈, C₁₀), Gerstenkorn, Luc
   & Amiot 1985 for B (C₅, C₆, C₈, C₁₀).
4. **Each state fitted to its own Gerstenkorn–Luc term values** (`global_fit.py --stage=dunham`).
   A joint transition fit is badly conditioned before either potential is right; a state fit is not,
   because the X levels depend on the X potential alone. T₀₀ enters as one extra row — without it the
   fit is blind to Te while De drifts 44 cm⁻¹ through the asymptote tie, and the pair ends 1.3 GHz out.
   Result: **16 MHz (X, v″ ≤ 19), 172 MHz (B, v′ ≤ 60)** against the atlas model, with
   **Te = 15769.0637, De = 4381.248 cm⁻¹** — within 0.003 cm⁻¹ of the B potential our earlier fit to
   measured levels reached by a completely different route.
5. **The global fit to our data** (`--stage=data`): MLR X (Re, β) and B (Te, Re, β) together, plus one
   calibration offset per atlas, against the atlas centres and the precision positions.

Before seeing any of our data, the pair from step 4 reproduces all 68 516 measured positions at
**186 MHz rms**. That is the honest measure of the 1985 atlas model against modern measurement.

## The fit

Derivatives are Hellmann–Feynman, ⟨ψ|∂V/∂θ|ψ⟩ on the solver's own quadrature grid, verified against
finite differences to 2×10⁻⁵ on every parameter. One B-spline solver per state per iteration (the
basis and the kinetic and potential matrices do not depend on J; only the rotational term does), so an
iteration costs one eigen-solve per J per state and nothing extra for the 36 parameters. A validity
wall returns a large flat residual where the parameters are not a pair of potentials, so the optimiser
can back out instead of crashing on a NaN.

**Loss.** Linear first. Started under `soft_l1` the fit dies at once: with residuals at 27σ every point
is in the down-weighted regime, the gradient flattens and the trust region collapses (two iterations,
step norm 2×10⁻⁴). Robust loss is a polish, not a start.

**Which atlas lines.** The two atlases measured 23 144 lines in common, and that split the error cleanly:
48 MHz is common to both (the model) and 23 MHz is each atlas's own. But the 16 984 lines only APO
detected sit 244 MHz from any model — weak detections whose fitted σ does not reflect their real error.
σ_fit is the discriminator (≤ 14 MHz → 40–67 MHz residual; 22–45 MHz and unconfirmed → 265 MHz), so the
fit takes σ ≤ 20 MHz and depth ≥ 0.15: 39 757 of 67 729 lines.

## Honest limits

- **The atlas positions are template-assisted.** They were measured as shifts of a model template, so
  for a blended line the result is pulled toward that template, which descends from Hannover. The
  0.01 cm⁻¹ ridge prior bounds it and a strong isolated line is unaffected, but the atlas data are not
  fully independent of the Hannover lineage. The Gerstenkorn–Luc starting point and the precision sets
  are.
- **The hyperfine structure is not fitted here.** The precision rows have their hyperfine offsets
  removed with the current model (BKT02/S06 plus our measured table), and the atlas template used the
  same. A fit of the potentials cannot absorb an error in the hyperfine Hamiltonian, and does not try.
- **Born–Oppenheimer corrections are absent**, so this fit is ¹²⁷I₂ only.

## The result, 2026-09-21 (converged)

*The first run stopped at 30 optimiser evaluations a stage and was not converged; with the fit made
5.2x faster an iteration it now gets 200, and the precision sets improve by a factor of four.*

| set | published | our default | **our fit** |
|---|---|---|---|
| `matyugin2012` (v″ = 48) | 233 738 | 6.5 | **1 748** |
| `nesterenko2019` (v″ = 53, 54) | 544 836 | 18.7 | **1 274** |
| `matsunaga2024a` (v′ = 45–50) | 24 905 | 1.5 | **66** |
| `yoshiki2023a` (v′ = 44–45) | 119 | 0.6 | **44** |
| `bipm2005a` (v′ = 58) | 23 662 | 0.2 | **6 612** |
| `xu2000a` (473 lines) | **0.94** | 0.88 | 18 |
| `velchev1998a` | **1.96** | 1.83 | 18 |
| `sansonetti1997a` | **1.63** | 1.28 | 24 |
| `liao2010a` | **2.30** | 0.34 | 16 |
| `bipm2012a` (532 nm) | **3.91** | 0.25 | 123 |
| atlas, Salami–Ross / APO | 2 909 / 6 385 | 37 / 59 | 63 / 174 |

Better by 130–430× wherever the published potentials leave their data, 10–20× worse inside it. The fit
itself reaches 76.7 MHz over its 40 578 lines, with the atlas calibrations at +4.26 MHz and +33.7 ppb
(Salami–Ross) and −2.37 MHz and −10.5 ppb (APO) — the same values every run finds, and consistent with
the atlas-against-atlas measurement made before any fitting.

**A floor on the v″ > 17 rows is part of the answer, not a fudge.** Taken at their own kHz
uncertainties those 36 rows dominate the cost (normalised residual 66) and the fit spends 40 000 lines
it can describe to move them from 1.53 to 1.09 GHz — `apo_nist_2009` goes 69 → 135 MHz. A single smooth
potential cannot reach v″ = 48–54 at kHz while describing everything else, and the honest way to say so
is a model-error floor, as `MODEL_FLOOR` does in the hyperfine fits. Every homotopy stage is now saved,
because the last is not always the best.

## The first run, 2026-09-21 (under-converged, kept for the record)

`data/potentials/i2spec_x_2026a.json` and `i2spec_b_2026a.json`: X with Re = 2.666320 Å, De = 12538.32
(Bacis's value refined by −0.07 %), 14 β; B with Re = 3.026347 Å, Te = 15769.0276, De = 4372.27, 20 β;
X long-range refined within its priors (C₈ −5.5 %, C₁₀ −6.7 %, both 0.2σ). Atlas calibrations fitted
alongside: Salami–Ross +4.3 MHz and +26.7 ppb, APO −1.9 MHz and −10.1 ppb — a 37 ppb difference, which
is the ~20 MHz drift between the two atlases that we measured independently before the fit saw it.

Against every observation (rms, MHz):

| set | published | our default | **our fit** |
|---|---|---|---|
| `matyugin2012` (v″ = 48) | 233 738 | 6.5 | **707** |
| `nesterenko2019` (v″ = 53, 54) | 544 836 | 18.7 | **1 563** |
| `matsunaga2024a` (v′ = 45–50) | 24 905 | 1.5 | **44** |
| `yoshiki2023a` (v′ = 44–45) | 119 | 0.6 | **32** |
| `bipm2005a` (v′ = 58) | 23 662 | 0.2 | **9 536** |
| `xu2000a` (473 lines) | **0.94** | 0.88 | 76 |
| `velchev1998a` | **1.96** | 1.83 | 38 |
| `sansonetti1997a` | **1.63** | 1.28 | 50 |
| `bipm2012a` (532 nm) | **3.91** | 0.25 | 127 |
| atlas, Salami–Ross / APO | 2 909 / 6 385 | 37 / 59 | 105 / 225 |

**What this is.** A genuinely global pair: 300–500× better than the published potentials wherever those
leave their data — v′ > 43 and v″ > 17, two thirds of the B–X system — and 20–100× worse inside it. The
published potentials bought 1–5 MHz in their core with a form that diverges outside; this one buys
validity everywhere at 30–200 MHz. They are different objects, and the published one also had Kato's
Doppler-free atlas, which we cannot obtain.

**The inner number is a bound, not a limit of the form.** Fitting the 787 precision rows *alone*, with
no atlas data at all, reaches 19 MHz and was still improving when it hit its evaluation limit — 100
evaluations is not many for 41 parameters, at 19 s each. And our own earlier pair, `mlr_x_2026c` /
`mlr_b_2026c`, fitted as alternating level fits, reaches 2–4 MHz on those same sets with the same
functional form. So an MLR can do several times better than this fit did; what is missing here is
optimiser iterations, not physics. That is the first thing to spend compute on.

## Bridging the gap, and what it confirms

The X potential's exponent is unconstrained between the atlas constants (v″ ≤ 19, R ≤ 3.07 Å) and the
only measurements above them (36 emission lines at v″ = 48, 53, 54). An MLR fitted to the first alone
reproduces v″ = 17 to **0.001 cm⁻¹** and is **408 cm⁻¹ low at v″ = 48**; with R_ref = R_e it diverges
outright. Handed the 36 levels at their own weight it throws the well away — they outweigh the 260
atlas levels 10⁸ to one. `prototypes/x_bridge.py` walks the weight down instead, from 100 cm⁻¹ to
0.05, and the potential crosses the gap: the implied levels come from 13.4 THz to 1.55 GHz while the
atlas region holds at 11–14 MHz.

What the crossing lands on matters more than the fit:

| v″ | bridged X (atlas constants + 36 lines) | `mlr_x_2026c` (our data fit) | difference |
|---|---|---|---|
| 19 | 3834.291 | 3834.286 | +0.005 cm⁻¹ |
| 30 | 5825.865 | 5825.582 | +0.283 |
| **42** | **7776.312** | **7775.652** | **+0.660** |
| 48 | 8651.335 | 8651.007 | +0.328 |
| 54 | 9450.722 | 9450.513 | +0.209 |

Two routes with no fitted input in common agree to 0.66 cm⁻¹ in the middle of a 28-level stretch nobody
has ever measured, where the published potential is about 21 cm⁻¹ away. With IodineSpec5's Martin-1986
column, which confirmed the same prediction to v″ = 28, the gap now has three independent supports.

## What was learned about fitting this system

- **The long-range tail is not the handle for an MLR.** `fitting.md` found it was the missing orthogonal
  direction for the ξ-form, and there it was. An MLR carries the correct tail by construction, so
  freeing De and C₆–C₁₀ under Bacis's own uncertainties moves them 0.1–0.2σ and fixes nothing. What is
  unconstrained is the exponent across the gap.
- **Thirty-six rows cannot win a joint fit against forty thousand.** They have to enter as what they
  are — direct constraints on one state's levels, through the other state's potential.
- **Robust loss cannot start a fit.** At 27σ every point is in the down-weighted regime, the gradient
  flattens, and the trust region collapses in two iterations. Linear first, robust to polish.
- **A homotopy on the weights is what crosses a gap**, in the same way the hyperfine fits needed a
  floor schedule. Jumping to full weight moved the high-v levels by a factor of 3; walking moved them
  by a factor of 8 600.

## Why the fit is 10-20x worse than Hannover inside his range, and what would fix it

The converged from-zero fit reaches 18-24 MHz on the precision sets where the published potentials
reach 1-5. Four explanations were tested; three are dead, and the fourth is the answer.

| hypothesis | test | verdict |
|---|---|---|
| the 40 000 atlas lines outweigh the 787 precision lines — a weighting problem | fit the precision rows *alone*, 500 evaluations | **no.** 15.7 MHz with no atlas data present at all, so weights cannot be the cause |
| R_ref = 1.25 R_e centres the exponent outside the well | whole chain refitted with R_ref = R_e, the long-range tail added as a constraint so the extrapolation stays physical | **no.** The state fits improve (14.5 → 10.5 MHz) but the joint fit is far worse: 332 MHz against 77 |
| a joint transition fit only sees E_B − E_X, leaving a degenerate direction | three rounds of alternating X-block / B-block fits | **no.** The total improves (80 → 70 MHz), the precision rows do not move at all |
| the from-zero start is in the wrong basin | start the same precision-only fit from `mlr_2026c` | **yes.** 3.58 MHz in 31 evaluations |

So the machinery, the objective and the weights are all sound: given a start near the right solution
they find 3.58 MHz immediately. The Gerstenkorn-Luc RKR chain lands in a different basin — one that
describes the atlas beautifully and the precision lines at 16 MHz — and no amount of gradient descent
leaves it. Note what `mlr_x_2026c` had that this chain does not: the *published Hannover levels* at
v'' <= 17 as shape guidance. Its 2.6 MHz is partly inherited, which is exactly why a from-zero chain is
worth having as an independent check even at 18 MHz.

**To be no worse than Hannover anywhere, in order of cost:**

1. **The hybrid, which is guaranteed and already exists.** Fitted potentials plus measured per-level
   corrections (`i2spec2026d`) are better everywhere by construction, because the corrections absorb
   whatever the potential cannot. IodineSpec5 does the same thing in the near infrared with its local
   Dunham model. This is the answer for the shipped model.
2. **A better start, for the pure potential.** Multi-start from perturbed RKR curves, or a start taken
   from the precision data alone before the atlas is added. The basin exists and is easy to sit in; the
   problem is finding it from 40 parameters away.
3. **Not more weight, and not more iterations.** Both were measured and neither moves it.
