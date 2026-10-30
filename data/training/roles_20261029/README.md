# Roles of the pending data sets (decision of 2026-10-28, refined 2026-10-29)

`roles.json` gives every pending data set one role, by set and lab, never by row: **train** (added to the fits),
**calibrate** (sets the uncertainty model), **locked test** (scored once, at the end) or **record** (no fit or
uncertainty uses it now). Its rules, the seed and the draw order of the locked test are in the file; each set's data
files are hashed, and `manifest.json` holds the file's own SHA256. The roles were fixed after the third blind test
(`src/i2spec/data/research/20261029/validation/pending_blind_v3`) was scored, by region, lab and kind only: no set
was placed by its residuals.

| role | sets | rows |
|---|---|---|
| train | 22 line and hyperfine sets (every set beyond the core) and 5 atlas line lists | 1,890 and 37,852 lines |
| calibrate | 95 sets in the core | 2,374 |
| locked test | 17 sets of 15 independent laboratories | 287 (16.4 % of the non-compilation core and isotopologue rows) |
| record | level and constants tables, widths, lifetimes, other states, theory, duplicates | — |

## Outcome of the training blocks (2026-10-30)

`roles.json` is kept as it was fixed. `outcome.json` records what the training blocks decided for the 27 train-role
sets: each block was added to the fits on its own and judged by the paired rule (added only if some region's
held-out error improves by more than 2 SE and none worsens by more than 2 SE, one value per held-out level); a block
not added goes to calibrate. Every variant tried is listed with the regions it flagged.

| block | sets | outcome |
|---|---|---|
| near dissociation (B v' 61-80) | gerstenkorn1983a, pique1981a, duburck2005a | calibrate; B v' 73-80 is a limit of the model |
| B hyperfine beyond the core | 8 sets (6 of lines, 2 of constants) | calibrate |
| other isotopologues | 8 sets | train (core fit) |
| atlas line lists | gerstenkorn1977, gerstenkorn1978, gerstenkorn1978b | train, with model assignments |

wallerand2006a and koffend1978a, in no block, go to calibrate; ikeda2020a (superseded by ikeda2022a), gerstenkorn1979
and gerstenkorn1980b (1978 atlas lines printed again) are record. Two groups of otherwise-used sets are record:
gerstenkorn1983a's Fourier rows (the Orsay Partie IV atlas printed again) and pique1981a's lambdameter rows (an
instrument offset of about -120 MHz). A measurement printed in two places enters the fits once, from the more current
record; separate measurements both enter.

| role | sets |
|---|---|
| train | 11 |
| calibrate | 108 |
| locked test | 17 |
| record | 72 |
