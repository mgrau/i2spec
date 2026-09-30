# Near the B dissociation limit

*2026-09-26. Preparation for the Orsay 19 700–20 035 cm⁻¹ volume.*

## Where the model stopped, and why

The B state dissociates 20 043.2 cm⁻¹ above X(0, 0). Two limits kept the model short of it:

- **Levels.** The default B-spline grid solves 60 levels (v′ ≤ 59, 19 899 cm⁻¹ at J = 0). Any request
  above that raised an IndexError, so R(26) 62-0 — the best-measured line in the region — could not be
  evaluated.
- **Line list.** The position and intensity grids end at 7.0 Å. Above v′ ≈ 62 the outer turning point
  passes that, the levels become box states, and `master_line_list` drops them: no lines above
  ~19 955 cm⁻¹.

## Grid convergence

At the validated knot spacing h = 0.01 Å, levels of the extended-range B (`mlr_b_2026c`) against a 30 Å box:

| box | converged to < 0.1 MHz (J = 0 / 26 / 100) |
|---|---|
| 8 Å | v′ ≤ 69 / 69 / 61 (all bound) |
| **12 Å** | **v′ ≤ 73 / 72 / 61** — every level bound by more than 0.3 cm⁻¹ |
| 22 Å | v′ ≤ 75 / 72 / 61 |

So the 8 Å box was nearly enough; the default grid was limited by its level count. h = 0.02 Å is *not*
converged (errors to 1.3 GHz) and must not be used to save time here.

## What changed

- `ADAPTIVE_GRIDS["B"]` gains a 12 Å, 80-level entry (v′ ≤ 73). `RovibronicModel.energy` takes any level
  beyond its default grid from the cheapest validated grid that reaches it, built on first use; the
  default solve is untouched. `_reference_term` (the hyperfine formulae's energy) does the same on the
  raw published curve.
- `intensity_model(grid=NEAR_DISSOCIATION)` puts the DVR box at 12 Å, and `position_model` follows it.
  `prototypes/near_dissociation_lines.py` builds the 19 600–20 050 cm⁻¹ list with it.
- Two comb-referenced data sets are entered: `goncharov2007a` (R(26) 62-0, all 15 components, 250 Hz)
  and `sakamoto2024a` (P(40) 52-0, P(52) 53-0, 8 kHz absolute, all 15 components each).

## Where the model stood (i2spec2026i)

*Superseded by i2spec2026k, refitted to the Orsay atlas Partie IV: R(26) 62-0 now 0.34 MHz, P(40) 52-0 and
P(52) 53-0 0.18 MHz rms (`docs/research/orsay-atlas-19700-20035.md`).*

| line | obs − model |
|---|---|
| R(26) 62-0 | −469 MHz |
| P(40) 52-0 | +271 MHz |
| P(52) 53-0 | +248 MHz |
| hyperfine intervals, these three lines | 108–260 kHz rms |

All inside the 2 GHz `lookup` quotes for v′ > 50. They are anchors for the refit with the atlas volume,
not yet fitted. The 1g perturbation at v′ ≈ 57–60 (R(98) 58-1 is 3.3 GHz off before its own one-line
correction) needs a coupled-channel treatment that a single potential cannot give.

## The line list to the limit

`prototypes/near_dissociation_lines.py`, 19 600–20 050 cm⁻¹, lines above 10⁻²⁶ cm at 295 K:

| grid | lines | highest v′ | highest line | lines above 19 955 cm⁻¹ |
|---|---|---|---|---|
| default (7 Å) | 3 009 | 62 | 19 943.80 cm⁻¹ | 0 |
| 12 Å | **6 435** | **73** | **20 041.11 cm⁻¹** | **738** |

The 12 Å build takes ~17 minutes once and is cached. On that grid the DVR intensity model (the published
B curve) holds a few more bound levels near the limit than the extended MLR that supplies positions
(77 against 74 at J = 10); those have no position and are treated as continuum, like box states. The cap
applies only on the 12 Å grid, so the default list is unchanged bit for bit.
