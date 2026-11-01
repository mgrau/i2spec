# Frozen release: October 31, 2026 — the released model (version 1.0.0)

Release 20261030 with the corrections of the final audit and the hyperfine beyond the core in its final form. The
core fit, the potentials beyond it and the isotope shift are those of release 20261030 (and 20261029), unchanged.
This is the model the apps and the documentation present (role `i2spec`):

```python
i2spec.ResearchModel("i2spec", isotopologue)        # or i2spec.preferred.model(isotopologue)
i2spec.ResearchModel("research-20261031-i2spec", isotopologue)
```

Release 20261030 stays as it was: the third and fourth blind tests and the locked test were scored with it.

## What changed from release 20261030

- **Levels above the dissociation limit (a correction).** For J ≤ 200 the levels B v′ 44–59 are solved in a box that
  ends at 8.6 Å. Release 20261030 kept every eigenvalue of that box, so above the limit (from J′ 132 at v′ 59 and
  J′ 197 at v′ 46) it counted standing waves between the centrifugal barrier and the box wall as levels. Here a
  level that a shorter box puts at or above the limit is solved in the long box with its barrier rule (90 % of the
  density inside the barrier; `structure.above_limit` in `expanded/model.json`; the fit engine applies the same
  rule). 485 ¹²⁷I₂ levels change (76 move, by up to 14 cm⁻¹; 409 are gone), 437 of ¹²⁹I₂ and 463 of ¹²⁷I¹²⁹I,
  all B v′ 46–59 at J′ 132–200 above the limit. No other level moves by more than 1 kHz, the v ladder stays
  continuous and ordered at every J, no training row touches an affected level (the continuation's training cost
  is unchanged to 10⁻¹² of its value), and no reference prediction's center moves.
- **The B hyperfine beyond the core.** The second-order form of release 20261030 missed its own training data at
  v′ 63–75 by MHz. The form here frees the 1g mixing angle (θ_g), adds curvatures of the C, d and δ backgrounds past
  the knot (kC, kd, kδ, the shape of eqQ's γ) and raises the anchored background orders of d and δ to 5 and 4 (they
  are fitted to the core fit's own predictions, not to data beyond the core): four fitted terms more, 110 fitted
  molecular coefficients in all. The 1g(¹Πg) C5 scale is chosen by training cost, as before, from 0.6–0.9: 0.85.
  Fitted: θ_u 11.10°, θ_g 109.79°, γ 37.0, kC 314, kd −1755, kδ −185 kHz; out-of-core training cost 40 961 (release
  20261030: 75 614), core anchor cost 5 185 (13 679).
  - *Selected after the blind scores were read.* The form was found among eleven candidates in a study made after
    the third blind test was scored; it is a post hoc choice, and no blind test of the region remains. In the nested
    cross-validation (ten v′ blocks; each of the five laboratories that train beyond the core left out in turn; the
    calibrate sets that never train) it improves eight regions by more than 2 SE and worsens none (with the fixed
    scale 0.65 it would worsen one). The potentials of that study are this release's, so its verdict stands.
  - The angles and curvatures are effective parameters: they move with the curve scale and absorb errors of the
    inherited perturber curves. Held-out misses remain at the MHz level beyond v′ 55.
  - **The B~1g(¹Πg) resonance (B v′ 77–78) keeps the 1g coupling of release 20261030** (`resonance`, its
    `form`): that release's hyperfine fit (1g mixing angle 30°, the separated-atom value; 1g(¹Πg) C5 scale 0.65;
    its backgrounds and θ_u), on this release's levels. The resonance block is release 20261030's (constants to
    10⁻¹³ MHz, overlaps, resonant level and coupling identical), so these levels take a different hyperfine form
    from the other levels beyond the core. The new form's θ_g is an effective parameter fitted away from the
    resonance; in the coupling of the resonant level it would move the spans of these lines' patterns by up to
    2.6 GHz, and no measurement at v′ 77–78 tests either. Above the resonance's last J′ (35 at v′ 77, 31 at v′ 78)
    the ordinary tables apply, as before.
- **Levels whose hyperfine is not predicted.** Where a second-order energy denominator nearly vanishes, a level of
  an inherited perturber curve lies within a few GHz and the sum is not valid. With the new form two B levels show
  it, (v′, J′) = (70, 33) and (71, 20) (eqQ −9.5 GHz and d −0.45 GHz at the latter); the nine of release 20261030
  (worst (74, 42)) are smooth now, since the 1g(¹Πg) scale moved the perturber levels. Unlike at v′ 77–78, no
  measured level places the perturber, so an exact diagonalization would rest on a level position the curves do
  not know: these levels are left out of the tables (`B.not_predicted`, with the reason; `B.spikes` keeps their
  rows as computed). Their lines have a center and no hyperfine components, and a neighbor whose J ± 2 basis would
  include such a level leaves it out.
- **The X spline window (audit, not changed).** A penalty on the curvature ripple of the X spline correction
  (5.6–6.6 Å) was tested: the weakest penalty of a grid fixed in advance that makes V″ monotone (1 cm⁻¹/Å²) is
  worse held out, under the paired rule, at X v″ 26–47 (+0.21 ± 0.03 in z²) and 61–89 (+0.53 ± 0.06). The window
  is kept as it was.
- **Uncertainty.** Rebuilt on this release: the uncertainty build on the cross-validation with the new nested
  hyperfine folds, then release 20261030's step-4 refit (its forms and protocol) with the calibrate residuals
  recomputed with this model. The covariance and the coverage tables are unchanged.

| term | this release | 20261030 |
|---|---|---|
| core σ_HF floor a, scale error κ | 63 kHz, 2.1 × 10⁻⁴ | the same |
| core σ_model, near precision rows / other | 0.162 / 0.288 MHz | the same |
| B σ_model at v′ 44 / 60 / 70 / 80 (J = 50) | 0.79 / 4.6 / 39 / 750 MHz | 0.79 / 4.6 / 39 / 749 MHz |
| B σ_HF floor at v′ 44 / 60 / 70 / 75 | 0.016 / 1.1 / 16 / 63 MHz | 0.063 / 1.8 / 14 / 41 MHz |
| B κ at v′ 44 / 55 / 70 | 4.6 × 10⁻⁵ / 1.0 × 10⁻³ / 5.1 × 10⁻³ | 5.3 × 10⁻⁵ / 1.5 × 10⁻³ / 8.9 × 10⁻³ |
| calibrate rows within 1σ / 2σ | 83.1 % / 96.5 % | 80.2 % / 96.2 % |

- **New terms at the edges** (`expanded/uncertainty.json`; earlier releases read none of them):
  - *X near dissociation* (`model_error["X beyond"].near_dissociation`): for X v″ ≥ 85, 101 MHz × e^(0.078 (v″ − 85))
    in quadrature with the level's share, fitted by maximum likelihood to the 67 cross-validation records there.
    Coverage within 1σ / 2σ 58 / 88 % without it, 66 / 93 % with it, each v″ left out in turn.
  - *Beyond the tables* (`edges`), rules, not fits: B v′ 81–86 take the v′ = 80 share at the same J times 3 per
    quantum (calibrated on three levels, v′ 81–83, of Danyluk & King 1977 with their reduction emulated:
    z = 1.0, 0.75, 0.5); X v″ 109–113 continue the near-dissociation term, capped at √(binding² + σ(D₀)²). The apps'
    level tables hold these levels (¹²⁷I₂ equivalents); their lines are flagged outside the fitted domain, and
    "level identity uncertain" where a level's σ exceeds half its vibrational spacing. The model's own API
    (`ResearchModel.predict`, `uncertainty`) keeps the fitted domain.
  - *Tunneling width* (`edges["B quasi-bound width"]`): a quasi-bound B level adds its width Γ = (E_v − E_v−1) T / 2π
    in quadrature, T the transmission of a parabolic barrier with the top and curvature of the effective potential.
    Below 1.2 kHz more than 5 cm⁻¹ under the barrier top; tens to hundreds of MHz within 2–3 cm⁻¹ of it, where a 60 Å
    box moves 52 of 59 sampled levels by less than Γ (or 1 kHz) and the others by up to 4.6 Γ (the worst, B (61, 125),
    by 246 MHz against Γ = 54 MHz).
- **Recipe text.** The covariance is over 37 tail terms and 46 core potential coefficients (the earlier text said 38
  and 47, the counts of release 20261020).

## Verification

Every model's references are replayed by `ResearchModel(id).verify()`: the core's exactly; the continued model's
fit-engine references (4 847 rows) to 2.4 Hz; the isotopologues' (762) exactly. Against release 20261030 the
references move only beyond the core in B: no center, no core or X line, and in the isotopologues' list only the
components of B v′ 45–70 lines (up to 6.9 MHz) and, in the continued model's, of B v′ ≥ 44 (median 0.05 MHz at
v′ 44, 50 MHz at v′ 74); the v′ 77–78 lines are release 20261030's (to 10 µHz). Fifteen reference rows of release 20261030 (R(30) 70-0/1/2)
are left out: an upper level's J ± 2 partner is not predicted.

## Locked test

Scored once with release 20261030 (`../20261030/validation/locked_test`), not scored again. Of its 164 scored rows,
150 have the same prediction (to 1 Hz) and σ (to 10⁻⁶) here; the other 14 are the hyperfine intervals of R(26) 62-0
(borde1981a), beyond the core (`manifest.json` "locked_test").
