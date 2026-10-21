# Data

Everything measured that the model is fitted to, tested against or compared with. Each data set sits in
its own folder or file pair, next to a `meta.toml` or `.toml` that gives its source, how it was read or
derived, its uncertainty and its terms. The catalog of every iodine measurement we know of, held here or
not, is in `catalog/`.

## What is here

| folder | contents | how it got here |
|---|---|---|
| `observations/` | precision line positions and hyperfine intervals: absolute frequencies of hyperfine components, BIPM recommended values, intervals between components and lines (67 sources) | transcribed from the published tables, the BIPM *mises en pratique* and one thesis |
| `hyperfine_parameters/` | effective hyperfine parameters (eqQ, C, d, δ) that a source fitted, one row per line | transcribed from the published tables |
| `x_levels/` | ground-state term values from laser-induced fluorescence | transcribed from the published tables |
| `atlas_lines/` | Doppler-limited line positions from Fourier-transform and laser atlases, most of them assigned; for `debus_2023` also line intensities | fitted here to published spectra (Salami & Ross 2005, the NIST scan of the APO cell, Debus *et al.* 2023, the 1993 Kitt Peak Keck cell); transcribed from printed atlases (Orsay 1982, all four parts; Orsay 1983); the authors' own list (Nölleke 2018, unassigned) |
| `spectra/` | measured spectra small enough to keep and free to share | as published (Nölleke 2018 scans) |
| `potentials/` | potential-energy curves fitted in this work | `prototypes/` |
| `catalog/` | the bibliography of iodine data (about 300 entries), with what each measures and where its numbers are | compiled from the literature |
| `external/` | third-party raw files, **not in the repository** (`.gitignore`) | see `external/README.md` |

Line positions and intensities from spectra whose terms do not allow redistribution are kept in a private
companion repository (`i2spec_data`, see `external/README.md`) until they do.

## Terms

The repository's code and the data derived in this work (fits, assignments, transcriptions) are under its
MIT licence. Measured values transcribed from publications are facts reported by their authors: each
`meta.toml` or `.toml` names the publication, which any use of the numbers should cite. Data sets
published under their own licence keep it, as recorded in their `.toml`:

| data | licence | credit |
|---|---|---|
| `atlas_lines/nolleke_2018.*`, `spectra/nolleke_2018_scans.csv.gz` | CC BY 4.0 | Nölleke *et al.*, J. Mol. Spectrosc. 346, 19 (2018), and their Mendeley Data set, doi:10.17632/cvr6vv2zps.1 |
| `atlas_lines/debus_2023.*` (fitted here to their spectra) | the spectra are CC0 1.0 | Debus, Schäfer & Reiners, JATIS 9, 045003 (2023), and doi:10.25625/PONFBL |
| `atlas_lines/pyodine_*` (fitted here to their spectra) | the spectra are MIT | Heeren *et al.*, A&A 674, A164 (2023), and gitlab.com/Heeren/pyodine |
