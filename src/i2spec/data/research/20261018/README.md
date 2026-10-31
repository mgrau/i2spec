# Frozen release: October 18, 2026 — the core fit

The exact-J core fit of ¹²⁷I₂ on which the model is built (role `balanced`; the directory names of the frozen
files are their artifact identifiers). Its domain is X v'' ≤ 17, B v' ≤ 43, J'' ≤ 200; it has 61 fitted and 18
inherited molecular coefficients and 4 calibration coefficients, and no line or level offsets. One line,
P(13) 43-0, is held out of every fit as a diagnostic.

`manifest.json` records the release's inventory, domain, benchmark scores and the SHA256 of every frozen file.
The files under `balanced/` are byte-for-byte copies of the verified snapshot, checked against those hashes when
`i2spec.research` reads them; source paths in `provenance.json` are historical records, not runtime dependencies.
`reference_predictions.json` holds 60 absolute and 592 intraline predictions, which the package reproduces within
1 Hz. These are reproduction checks, not independent predictive validation.

Release 20261021 (`../20261021/`) reads these coefficients as its core; the released model, release 20261030, has its own core.
