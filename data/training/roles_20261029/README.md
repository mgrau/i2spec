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
