# The X state to its dissociation limit: Martin et al. 1986

F. Martin, R. Bacis, S. Churassy and J. Vergès, *J. Mol. Spectrosc.* **116**, 71–100 (1986), Table I: G, B,
D, H, L (and M) of 93 X levels of ¹²⁷I₂, v″ = 8–108, from 14 820 B→X fluorescence lines (LIF Fourier
transform), transcribed in `data/x_levels/martin1986/`. Cerny, Bacis & Vergès, *ibid.* 458 (1986), gives
the same for ¹²⁷I¹²⁹I and ¹²⁹I₂ (`data/x_levels/cerny1986_*`).

## What the previous X showed

`mlr_x_2026d`, fitted only to v″ ≤ 25 and the emission levels v″ = 48, 53, 54, followed Martin to 15–170 MHz
at v″ = 26–47 (J ≤ 120) — the first measurement of that range, and a confirmation of the interpolation — but
departed above v″ ≈ 60: +0.24 cm⁻¹ at v″ = 70, +9 at 80, +46 at 89, about +70 cm⁻¹ near v″ = 95.

## The refit (`prototypes/mlr_x_martin.py`, `mlr_x_2026e`)

Targets: the v″ ≤ 17 levels of the previous X at 3 MHz, the corrected atlas and emission levels, and
Martin's levels at J inside their coverage at the paper's precision of recomputed lines (0.005 cm⁻¹ to
J = 80, 0.01 to 100, 0.02 to 120), with one origin offset (Martin's G(9) is Luc's value, ~1 mK from the
comb-referenced levels) held to 0(2) mK. Martin's levels enter in stages, v″ ≤ 60, 75, 89, 108; 18 β.

| stage | anchor | atlas | emission | Martin 26–47 | 49–60 | 61–75 | 76–89 | offset |
|---|---|---|---|---|---|---|---|---|
| start (2026d) | 0.00 | 16.4 | 3.2 | 78 | 174 | 2.9 × 10⁴ | 6.8 × 10⁵ | 0 |
| v″ ≤ 60 | 0.39 | 8.6 | 4.3 | 27 | 36 | — | — | −60 |
| v″ ≤ 75 | 0.45 | 8.7 | 8.0 | 34 | 63 | 67 | — | −63 |
| v″ ≤ 89 | 0.42 | 13.3 | 28.0 | 60 | 55 | 73 | 174 | −102 |

(rms, MHz.) The v″ ≤ 89 stage is `mlr_x_2026e`. The v″ ≤ 108 stage diverged: the last levels, within a few
cm⁻¹ of the asymptote and at J ≤ 20, pulled the whole curve GHz away. The emission and atlas levels, which
lose some accuracy in the refit, are restored by their level corrections.

**Beyond v″ = 89 `mlr_x_2026e` is not physical.** It dips up to 170 cm⁻¹ below its own dispersion limit
𝔇ₑ − u(R) at 7–9 Å, which gives it 137 bound J = 0 levels in a 120 Å box, 15 more than the 2008 curve,
and its levels v″ = 91–108 sit 70 cm⁻¹ rms from Martin's. Three ways of adding those levels were tried
(`--from=89`, 2026-09-30), each continuing from the v″ ≤ 89 stage:

| run | anchor | atlas | emission | 26–47 | 49–60 | 61–75 | 76–89 | 91–108 | offset |
|---|---|---|---|---|---|---|---|---|---|
| v″ ≤ 89 (`mlr_x_2026e`) | 0.42 | 13.3 | 28 | 60 | 55 | 73 | 174 | 2.1×10⁶ | −102 |
| 91–108 at σ = 1 cm⁻¹ (`--soft=1.0`) | 2.85 | 35.7 | 84 | 110 | 95 | 121 | 353 | 34 531 | −215 |
| then at Martin's σ, 18 β | 3.65 | 59.7 | 190 | 125 | 155 | 272 | 763 | 3 449 | −245 |
| then 22 β | 9.42 | 83.8 | 103 | 251 | 333 | 278 | 611 | 1 540 | −529 |
| soft, then σ ≥ 0.05 cm⁻¹ above 89 (`--high=0.05`) | 5.36 | 57.0 | 101 | 168 | 190 | 159 | 426 | 5 584 | −342 |

(rms, MHz.) Every curve that follows the last levels gives up 2–10× on the measured ones below, and drives the
origin offset 3–9σ from its 0(2) mK prior: a single-channel curve cannot hold both. Martin *et al.* 1983 place
the X–a′–a interactions that perturb the last X levels near 5 Å. `mlr_x_2026e` is kept, and the lookup
flags v″ > 89 as not physical (2.1 THz). No B–X line is affected: the highest B level (20 043 cm⁻¹) minus
X v″ = 90 (≈ 12 390 cm⁻¹) is below 7 700 cm⁻¹, and the line list starts at 9 361 cm⁻¹. The partition
function is not affected either: those levels lie above 12 300 cm⁻¹.

## Isotopologues

Once the part of the model error the three isotopologues share is removed, Cerny's ¹²⁹I₂ and ¹²⁷I¹²⁹I levels
follow pure mass scaling of the X curve to 2–5 mK, inside the stated origin error: no X-state
Born–Oppenheimer term is detectable at this precision.

## Next

A joint fit of the three isotopologues near dissociation, which sample the long range at different
energies; the levels above v″ = 89 with a better-conditioned start (fit C6, C8 and De to v″ ≥ 90 alone first).
