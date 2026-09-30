# The quasi-bound B levels of the Orsay atlas Partie IV

The Orsay atlas Partie IV (`docs/research/orsay-atlas-19700-20035.md`) holds 70 unblended lines whose upper
level lies above the B asymptote (20 150.316 cm⁻¹ above the X minimum) or within 0.3 cm⁻¹ of it. They
cover v′ = 63–77 and J′ = 35–96, up to 11.2 cm⁻¹ above the asymptote. Up to i2spec2026m they were left out
of the fit, on the assumption that the bound-state solver could not represent them. That was only partly
true.

## They are narrow resonances

`prototypes/quasibound_b.py` takes each level on the bare extended B potential (`mlr_b_2026d`) and computes:

- **Box size.** The eigenvalue in the model's 40 Å graded box and in an 80 Å box. For every level
  localised inside the barrier the two agree to better than 10⁻⁴ MHz.
- **Barrier.** The barrier top lies 1.0–19 cm⁻¹ above the level.
- **Width.** A WKB tunnelling width, Γ = (ΔG/2π) exp(−2θ), is below 10⁻⁷ MHz for all of them. The
  largest, 3 × 10⁻⁸ MHz, is (v′, J′) = (74, 52), R(51) 74-0, 1.2 cm⁻¹ below its barrier top.

At this precision these are bound states. The atlas sees them as sharp lines, and so does the solver,
given a large enough box.

## What went wrong in the model

Above the asymptote the 40 Å box also holds continuum states, standing waves between the barrier and
the box wall. Two problems followed:

1. **The v count shifted.** `model.energy(v, J)` counted every box state, so at low J, where the barrier
   is 1–2 cm⁻¹ high, continuum states slipped in below a resonance. Four levels came out 1.1–43 GHz off:
   (v′, J′) = (76, 41), (75, 45), (74, 49) and (74, 52).
2. **The ladder was cut short.** The grid kept only 95 eigenvalues. At J′ ≈ 80–96 the continuum states
   use up that budget before the highest resonances, so 17 levels got an unrelated state.

The fix is `i2spec.model.resonance_ladder`. On the dissociation grid it keeps every state below the
asymptote, and above it only the states with at least 90 % of their probability inside the centrifugal
barrier. v then counts physical levels. The grid keeps 260 eigenvalues. The Partie IV fit
(`prototypes/orsay4_fit.py`) now numbers its levels the same way and fits all 2 079 unblended lines.

| lines, obs − i2spec2026n (MHz) | n | median | robust spread | max |
|---|---|---|---|---|
| near-limit, before (box count; not fitted) | 70 | — | — | 43 000 |
| near-limit, numbered as resonances, not fitted | 70 | −39 | 118 | 441 |
| near-limit, numbered as resonances, fitted | 70 | −6 | 92 | 337 |
| the other Partie IV lines | 2 009 | −3 | 38 | 485 |

The near-limit spread matches the atlas's own at v′ ≥ 70, where the per-level held-out spread is
60–120 MHz.

## Not covered

- **Intensities.** The explorer's line list near dissociation comes from a 12 Å box
  (`intensity.NEAR_DISSOCIATION`), whose physical prefix stops at the first box state. The quasi-bound
  lines are therefore not in the exported list yet.
- **Broad resonances.** Levels within about 0.1 cm⁻¹ of their barrier top would be broad. None of the 70
  lines is that close.
