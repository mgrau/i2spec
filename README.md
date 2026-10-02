# i2spec

An open, physics-based model of the molecular iodine (I₂) B–X spectrum, with an uncertainty on every
line. It computes level energies from potential-energy curves, line positions, hyperfine components
from an effective nuclear-spin Hamiltonian, line strengths, absorption spectra and cell transmission,
for ¹²⁷I₂, ¹²⁹I₂ and ¹²⁷I¹²⁹I from the dissociation limit near 499 nm into the near infrared.

The model keeps the potentials of Salumbides et al. (2008) where they were fitted (X v″ ≤ 17,
B v′ ≤ 43), uses Morse/long-range curves fitted here beyond them, and adds measured level corrections
from every comb-referenced line, with a Gaussian process for the B levels that have no data of their
own; it reproduces the precision data sets to 0.02–0.3 MHz where they exist. The web app's
"How it works" page summarises the physics.

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
zoom, `[`/`]` pan, `c` centres on the selected line, `l` toggles a log scale, and `p` switches the
lower plot between the selected line's Doppler profile and its sub-Doppler trace (Lamb dips and
crossovers). `s` sorts by strength, `a` shows the weak ΔF ≠ ΔJ components, `/` returns to the search box.

Positions carry a 1σ uncertainty estimated from held-out data, not from a fit, with the reason:
0.05–1 MHz on measured levels, tens of MHz for the atlas-measured X v″ = 18–25, up to GHz where
nothing has been measured.

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

Pushing to `main` runs `.github/workflows/pages.yml`, which builds the site and publishes it to GitHub
Pages (Settings → Pages → Source: GitHub Actions). See [`docs/design/web-app.md`](docs/design/web-app.md)
for the data format and the deliberate limits (no sub-Doppler, no continuum, one isotopologue at a time).

## Python

```python
from i2spec import MHZ_PER_CM, RovibronicModel

m = RovibronicModel("127I2")
nu = m.transition(32, 0, 56, "R")   # R(56) 32-0, hyperfine-free, cm^-1
print(nu, nu * MHZ_PER_CM)          # ≈ 18788.3355 cm^-1 ≈ 563 260 129 MHz (532 nm)

nu0, comps = m.hyperfine_components(32, 0, 56, "R")   # MHz; offsets relative to nu0
a10 = next(c for c in comps if c.label == "a10")
print(nu0 + a10.offset)             # ≈ 563 260 223.57 MHz (CIPM value: 563 260 223.513 MHz)
```

## Repository layout

| path | what |
|---|---|
| `src/i2spec/` | the package: potentials, radial solvers, level corrections, hyperfine, intensities, spectra, lookup, CLI, TUI, web export |
| `src/i2spec/data/` | the model parameters: potentials, level corrections, hyperfine tables |
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

- Potentials: E. J. Salumbides et al., *Eur. Phys. J. D* **47**, 171 (2008), Table 1. This
  supersedes the misprinted Table 4 of H. Knöckel, B. Bodermann & E. Tiemann, *Eur. Phys. J. D*
  **28**, 199 (2004). See `docs/research/hannover-model-reproduction.md`.
- Hyperfine parameters: B. Bodermann, H. Knöckel & E. Tiemann, *Eur. Phys. J. D* **19**, 31 (2002),
  eqs. (10)–(15). For the conventions we established, see `docs/research/hannover-model-reproduction.md`.
- Research notes, data catalogs and the plan are in `docs/` and `data/catalog/`.

## License

MIT, copyright "i2spec contributors"; see [`LICENSE`](LICENSE).
