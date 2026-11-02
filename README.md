# i2spec

An open, physics-based model of the B–X absorption spectrum of molecular iodine (I₂). It provides the position,
strength and hyperfine structure of every rovibronic line of ¹²⁷I₂, ¹²⁹I₂ and ¹²⁷I¹²⁹I from the B-state
dissociation limit (499 nm) into the near infrared. Every position has a 1σ uncertainty. From the line list the
package computes absorption spectra and the transmission of an iodine cell.

The X and B states are each represented by a Morse/long-range potential with a rotational function, and the B state
also by radial functions for its hyperfine parameters. These are fitted together to every precision measurement of
¹²⁷I₂ with X v″ ≤ 17 and B v′ ≤ 43 and to two Fourier-transform atlases. They make up the core fit. Further terms
that act only beyond the core continue the model to both dissociation limits (X v″ ≤ 108, B v′ ≤ 80). No term is
specific to an individual level or line.

¹²⁹I₂ and ¹²⁷I¹²⁹I use the same potentials at their own reduced mass. A Born–Oppenheimer correction of two
coefficients is fitted to the CIPM isotope shifts and predicts each line from the others to 0.08 MHz. The hyperfine
functions are scaled by the moments of each nucleus. Line strengths use the transition moment of Tellinghuisen
(2011) with the model's own wavefunctions.

The uncertainty of each prediction comes from a frozen uncertainty model. It combines the covariance of the fit with
the error of predicting levels and hyperfine intervals held out of the fit. Data sets outside the fit calibrate it.
The core fit predicts held-out comb data to 0.18 MHz for absolute components and to 0.05 MHz for hyperfine
intervals.

A locked test of 17 data sets, used in no fit or calibration and scored once, found 88 % of 164 measurements within
2σ. The absolute frequencies (98 %) and the measurements of the other isotopologues met the target. The hyperfine
intervals within a line (78 %) did not; most of the shortfall comes from one data set
([locked test](https://github.com/mgrau/i2spec/blob/main/src/i2spec/data/research/20261030/README.md)). The released
version then changed the B hyperfine structure beyond the core and the uncertainty at the edges of the domain. These
changes leave 150 of the 164 locked-test rows as they were
([release notes](https://github.com/mgrau/i2spec/blob/main/src/i2spec/data/research/20261031/README.md)).

## Quick start

Everything runs through [uv](https://docs.astral.sh/uv/) from a clone of this repository; the first command creates
the environment.

```sh
uv run i2spec tui                                # the terminal browser (also plain `uv run i2spec`)
uv run i2spec web                                # the web app at http://localhost:8777
uv run i2spec lines 532.2 532.3                  # every line in a range, with uncertainties and flags
uv run i2spec line "R(56) 32-0"                  # one line with its hyperfine components
```

The first search for an isotopologue builds its line list in a few minutes. The line list is stored in
`~/.cache/i2spec` for later searches. The first `i2spec web` also exports the model for the app to `web/data/`.

## Line lookup

```sh
uv run i2spec lines 18788 18789 --unit cm-1 --sort strength --limit 10
uv run i2spec lines 632.9 633.0 -i 129I2 --json  # any isotopologue; JSON for scripting
```

In the terminal browser the top panel shows the absorption spectrum of the whole search range. The keys:

| key | action |
|---|---|
| `+` / `-` | zoom |
| `[` / `]` | pan |
| `c` | center on the selected line |
| `l` | toggle a log scale |
| `p` | switch the lower plot between the Doppler profile and the sub-Doppler trace (Lamb dips and crossovers) |
| `s` | sort by strength |
| `a` | show the weak ΔF ≠ ΔJ components |
| `/` | return to the search box |

Each position has a 1σ uncertainty, listed with its reason. It is 0.2–1 MHz in the core and 1–40 MHz for most
levels beyond it (to B v′ 70 and X v″ 80). For the last B and X levels it reaches 0.1–0.4 GHz. Lines whose levels
lie beyond the data are flagged *extrapolated*. `i2spec line` lists the 1σ of each hyperfine component and every
measurement of the line with its own uncertainty.

## Web app

The site has two parts:

- the **explorer** in `web/`. It searches a range in vacuum or air wavelength, wavenumber, THz or MHz, and plots
  the cross section or the transmission of a cell. A sortable table lists each line's strength and uncertainty,
  and the sources of measured lines with links to their DOIs. The detail panel shows the hyperfine components and
  every measurement of the selected line. `i2spec web export` writes the model as JSON for the app.
- the **documentation** in `documentation/`, built with MkDocs Material: installation and usage, the calculation,
  and every measurement source.

```sh
uv run i2spec web                    # build the site (explorer + docs) and serve it locally
uv run mkdocs serve                  # the documentation alone, reloading as you edit
uv run i2spec web export             # re-export after a model change
uv run i2spec web build site         # the static site in site/, as GitHub Pages publishes it
node web/test_physics.mjs            # check the browser physics against the Python model
```

Pushing to `main` runs `.github/workflows/pages.yml`, which builds the site and publishes it to GitHub
Pages (Settings → Pages → Source: GitHub Actions).

## Python

```python
from i2spec.preferred import model

m = model()                                        # the model, for 127I2 (model("129I2") for another)
m.frequency("R(56) 32-0", "a10")                   # ≈ 563 260 223.41 MHz (CIPM: 563 260 223.513 MHz)
u = m.uncertainty("R(56) 32-0", "a10")             # 1σ 0.39 MHz: statistical, model and hyperfine parts
m.uncertainty("R(56) 32-0", "a15", "a1").sigma_MHz # 0.20 MHz for an interval within the line
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

The model's frozen releases hold every coefficient, table and reference prediction under SHA256:

```sh
uv run i2spec model                  # its domain, coefficient count and held-out accuracy
uv run i2spec model verify           # replay its frozen reference predictions
```

## Repository layout

| path | what |
|---|---|
| `src/i2spec/` | the package: the model, radial solvers, hyperfine, intensities, spectra, lookup, CLI, TUI, web export |
| `src/i2spec/data/` | the published 2008 potentials, for comparison |
| [`src/i2spec/data/research/`](https://github.com/mgrau/i2spec/blob/main/src/i2spec/data/research/20261031/README.md) | the model's frozen releases: coefficients, hyperfine tables, uncertainty model, audits, reference predictions and the blind tests |
| `data/observations/` | the measured line positions and hyperfine intervals, one directory per source |
| `data/atlas_lines/` | assigned FTS-atlas line positions (Salami & Ross, NIST/APO, Orsay) |
| `data/catalog/` | the literature catalog: every data source found, with DOIs |
| `web/` | the explorer, a static web app; `web/data/` is generated |
| `documentation/`, `mkdocs.yml` | the documentation site (MkDocs Material) |
| `tests/` | `uv run pytest` |
| `notebooks/` | an executed walkthrough |

## Walkthrough

[`notebooks/walkthrough.ipynb`](https://github.com/mgrau/i2spec/blob/main/notebooks/walkthrough.ipynb) steps
through one calculation:

- potential curves, vibrational and rotational levels, and line positions;
- the hyperfine structure, compared with the BIPM tables;
- intensities, and a spectrum compared with the Salami & Ross atlas;
- the continuum and the absolute cross section.

It is saved with its outputs. To run it yourself:

```sh
uv run --group notebook jupyter lab notebooks/walkthrough.ipynb
```

## Development

```sh
uv run pytest                         # the test suite (about two minutes; builds line lists on first run)
```

## Sources

- The model: Morse/long-range form of R. J. Le Roy & R. D. E. Henderson, *Mol. Phys.* **105**,
  663 (2007); second-order B hyperfine of L. Chen, W. A. de Jong & J. Ye, *J. Opt. Soc. Am. B* **22**, 951
  (2005). See the documentation's reference list.
- The published 2008 model (for comparison): E. J. Salumbides et al., *Eur. Phys. J. D* **47**, 171 (2008), Table 1. This
  supersedes the misprinted Table 4 of H. Knöckel, B. Bodermann & E. Tiemann, *Eur. Phys. J. D*
  **28**, 199 (2004).
- Published hyperfine formulae (the 2008 model, and lines beyond the model's domain): B. Bodermann, H. Knöckel &
  E. Tiemann, *Eur. Phys. J. D* **19**, 31 (2002), eqs. (10)–(15); E. J. Salumbides et al., *Mol. Phys.* **104**,
  2641 (2006).
- The data catalog, with every source and its DOI, is in `data/catalog/`.

## License

MIT, copyright "i2spec contributors"; see [`LICENSE`](https://github.com/mgrau/i2spec/blob/main/LICENSE).
