# i2spec

An open, physics-based model of the molecular iodine (I₂) B–X spectrum, with an uncertainty on every
line. It computes level energies from potential-energy curves, line positions, hyperfine components
from an effective nuclear-spin Hamiltonian, line strengths, absorption spectra and cell transmission,
for ¹²⁷I₂, ¹²⁹I₂ and ¹²⁷I¹²⁹I from the dissociation limit near 499 nm into the near infrared.

The model has one Morse/long-range potential per state, with a fitted rotational function, and radial
B-state hyperfine functions. These are fitted in one problem to every precision measurement of ¹²⁷I₂ and two
Fourier-transform atlases with X v″ ≤ 17 and B v′ ≤ 43, then continued to both dissociation limits
(X v″ ≤ 108, B v′ ≤ 80) by terms that act only beyond that core. It has no correction to any individual level
or line. Every prediction comes with a 1σ from a frozen uncertainty model: the fit's covariance, plus the
error of predicting levels and hyperfine intervals held out of the fit. It predicts held-out comb data to
0.16 MHz (absolute components) and 0.06 MHz (hyperfine intervals), and on 35 data sets scored blind against
frozen predictions 90 % of the measurements lie within 2σ.

¹²⁹I₂ and ¹²⁷I¹²⁹I use the same potentials at their own reduced mass, with a Born–Oppenheimer correction of two
coefficients fitted to the CIPM isotope shifts (each line predicted from the others to 0.08 MHz) and the
hyperfine functions scaled by the moments of each nucleus. Line strengths use the transition moment of
Tellinghuisen (2011) with the model's own wavefunctions.

## Quick start

Everything runs through [uv](https://docs.astral.sh/uv/); the first command creates the environment.

```sh
uv run i2spec tui                                # the terminal browser (also plain `uv run i2spec`)
uv run i2spec web                                # the web app at http://localhost:8777
uv run i2spec lines 532.2 532.3                  # every line in a range, with uncertainties and flags
uv run i2spec line "R(56) 32-0"                  # one line with its hyperfine components
```

The first search for an isotopologue builds its line list (a few minutes, in parallel) and caches it
in `~/.cache/i2spec`. The first `i2spec web` also exports the model for the app (`web/data/`).

## Line lookup

```sh
uv run i2spec lines 18788 18789 --unit cm-1 --sort strength --limit 10
uv run i2spec lines 632.9 633.0 -i 129I2 --json  # any isotopologue; JSON for scripting
```

In the terminal browser the top panel is the absorption spectrum of the whole search range: `+`/`-`
zoom, `[`/`]` pan, `c` centers on the selected line, `l` toggles a log scale, and `p` switches the
lower plot between the selected line's Doppler profile and its sub-Doppler trace (Lamb dips and
crossovers). `s` sorts by strength, `a` shows the weak ΔF ≠ ΔJ components, `/` returns to the search box.

Positions carry a 1σ uncertainty and the reason for it: 0.2–1 MHz in the core, 1–30 MHz for most levels
beyond it, and 0.1–0.2 GHz for the last X levels; lines whose levels lie beyond the data are flagged
*extrapolated*. `i2spec line` also gives each hyperfine component's 1σ. Every measurement of a line is listed
beside the model, with its own uncertainty.

## Web app

The site has two parts: the **explorer** in `web/` (a range search in vacuum or air
wavelength, wavenumber, THz or MHz; the spectrum or the transmission of a cell you specify; a sortable
table with each line's strength, uncertainty and — for measured lines — its sources, linked to their
DOIs; the hyperfine components and every measurement of the selected line), and the **documentation**
in `documentation/`, built with MkDocs Material: getting started, how the calculation works, and
every measurement source. Nothing is solved in
the browser: `i2spec web export` writes the model out as JSON and the app reads it.

```sh
uv run i2spec web                    # build the site (explorer + docs) and serve it locally
uv run mkdocs serve                  # the documentation alone, reloading as you edit
uv run i2spec web export             # re-export after a model change
uv run i2spec web build site         # the static site in site/, as GitHub Pages publishes it
node web/test_physics.mjs            # check the browser physics against the Python model
```

The research progress dashboard runs in the terminal:

```sh
uv run i2spec progress
```

It refreshes every 30 seconds and shows remaining refits, accepted and failed
attempts, worker log activity, finish estimates, and follow-up studies. Press
`r` to refresh, `space` to pause display updates, `1`–`4` to switch views, and
`q` to quit. Fits keep running when the dashboard is paused or closed. Use
`--interval 10` for faster updates, `--out PATH` for another artifact directory,
or `--json` for a single machine-readable snapshot. Counts and ETA cover the
current validation queues; failed attempts count as finished. After the original
queues finish, the overview shows the numerical recoveries, starting-value
checks, potential/hyperfine combinations, joint-data pilots and exact-J numerical
refits from the local finalization registry. Their counts describe the current stages; completed
pilots still need scientific review. The dashboard includes supplemental workers.

Pushing to `main` runs `.github/workflows/pages.yml`, which builds the site and publishes it to GitHub
Pages (Settings → Pages → Source: GitHub Actions). See [`docs/design/web-app.md`](docs/design/web-app.md)
for the data format and the deliberate limits (no sub-Doppler, no continuum, one isotopologue at a time).

## Python

```python
from i2spec.preferred import model

m = model()                                        # the model, for 127I2 (model("129I2") for another)
m.frequency("R(56) 32-0", "a10")                   # ≈ 563 260 223.36 MHz (CIPM: 563 260 223.513 MHz)
u = m.uncertainty("R(56) 32-0", "a10")             # 1σ 0.42 MHz: statistical, model and hyperfine parts
m.uncertainty("R(56) 32-0", "a15", "a1").sigma_MHz # 0.058 MHz for an interval within the line
```

The published 2008 model (the potentials of Salumbides et al., for comparison):

```python
from i2spec import MHZ_PER_CM, RovibronicModel

m = RovibronicModel("127I2", "hannover2008")
nu = m.transition(32, 0, 56, "R")   # R(56) 32-0, hyperfine-free, cm^-1
print(nu, nu * MHZ_PER_CM)          # ≈ 18788.3356 cm^-1 ≈ 563 260 131 MHz (532 nm)

nu0, comps = m.hyperfine_components(32, 0, 56, "R")   # MHz; offsets relative to nu0
a10 = next(c for c in comps if c.label == "a10")
print(nu0 + a10.offset)             # ≈ 563 260 225.51 MHz (CIPM value: 563 260 223.513 MHz)
```

The fit's frozen releases, with every coefficient, table and reference prediction under SHA256, can be listed
and replayed:

```sh
uv run i2spec research list
uv run i2spec research verify --model MODEL      # MODEL: an id from the list
```

See [the release inventory and API](docs/research-models.md) and the
[prospective validation guide](docs/research-model-validation.md), which provides frozen laboratory targets and
measurement-scoring commands.

## Repository layout

| path | what |
|---|---|
| `src/i2spec/` | the package: potentials, radial solvers, level corrections, hyperfine, intensities, spectra, lookup, CLI, TUI, web export |
| `src/i2spec/data/` | the model parameters: potentials, level corrections, hyperfine tables |
| [`src/i2spec/data/research/20261018/`](src/i2spec/data/research/20261018/README.md) | the core fit's frozen release: coefficients, source hashes, validation and reference predictions |
| [`src/i2spec/data/research/20261020/`](src/i2spec/data/research/20261020/README.md) | the fitted model's release: its continuation, hyperfine tables, uncertainty model, audits, validation and the blind tests |
| [`models/research-20261008/`](models/research-20261008/README.md) | earlier frozen research candidate and standalone evaluator |
| `data/observations/` | the measured line positions and hyperfine intervals, one directory per source |
| `data/atlas_lines/` | assigned FTS-atlas line positions (Salami & Ross, NIST/APO, Orsay) |
| `data/catalog/` | the literature catalog: every data source found, with DOIs |
| `web/` | the explorer, a static web app; `web/data/` is generated |
| `documentation/`, `mkdocs.yml` | the documentation site (MkDocs Material) |
| `tests/` | `uv run pytest` |
| `prototypes/` | the research scripts behind each fit and figure; see `prototypes/README.md` |
| `docs/design/`, `docs/research/` | design records and research notes |
| `notebooks/` | an executed walkthrough |

## Walkthrough

[`notebooks/walkthrough.ipynb`](notebooks/walkthrough.ipynb) steps through one calculation: potential curves, vibrational
and rotational levels, line positions, hyperfine structure against the BIPM tables, intensities, a
spectrum against the Salami & Ross atlas, the continuum, and the absolute cross section. It is saved
with its outputs. To run it yourself:

```sh
uv run --group notebook jupyter lab notebooks/walkthrough.ipynb
```

## Development

```sh
uv run pytest                         # the test suite (a few minutes; builds line lists on first run)
uv run --group research python prototypes/<script>.py   # a research script that needs matplotlib or h5py
```

## Sources

- The model: Morse/long-range form of R. J. Le Roy & R. D. E. Henderson, *Mol. Phys.* **105**,
  663 (2007); second-order B hyperfine of L. Chen, W. A. de Jong & J. Ye, *J. Opt. Soc. Am. B* **22**, 951
  (2005). See the documentation's reference list.
- The published 2008 model (for comparison): E. J. Salumbides et al., *Eur. Phys. J. D* **47**, 171 (2008), Table 1. This
  supersedes the misprinted Table 4 of H. Knöckel, B. Bodermann & E. Tiemann, *Eur. Phys. J. D*
  **28**, 199 (2004). See `docs/research/hannover-model-reproduction.md`.
- Hyperfine parameters: B. Bodermann, H. Knöckel & E. Tiemann, *Eur. Phys. J. D* **19**, 31 (2002),
  eqs. (10)–(15). For the conventions we established, see `docs/research/hannover-model-reproduction.md`.
- Research notes, data catalogs and the plan are in `docs/` and `data/catalog/`.

## License

MIT, copyright "i2spec contributors"; see [`LICENSE`](LICENSE).
