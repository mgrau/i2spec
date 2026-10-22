# Installation and usage

i2spec is a Python package managed with [uv](https://docs.astral.sh/uv/). The first `uv run` command
creates the environment and installs the dependencies.

```sh
git clone https://github.com/mgrau/i2spec.git
cd i2spec
# terminal browser (also: uv run i2spec)
uv run i2spec tui
# line explorer and documentation, at http://localhost:8777
uv run i2spec web
```

The first query for an isotopologue solves the model's levels and computes its line list, and caches both
in `~/.cache/i2spec` (230 MB for all three isotopologues). The first `i2spec web` also exports the data of
the line explorer to `web/data/` (590 MB): the line lists of all three isotopologues and the hyperfine
pattern of every line. After that, every query takes under a second. With an empty cache:

| first run | Apple M5 Max (18 cores, 128 GB) |
|---|--:|
| `i2spec lines` (¹²⁷I₂) | 2 min 45 s |
| `i2spec web`, with the export | 16 min 30 s |

About half of the export is solving the levels of the three isotopologues. It needs about 5 GB in its main
process and 1-2 GB in each worker process: with three workers (a 4-core machine, or `I2SPEC_WORKERS=3`), it
peaks at 10 GB.

## Command line

```sh
# all lines in a range, with uncertainties
uv run i2spec lines 532.2 532.3
# in cm-1, the ten strongest
uv run i2spec lines 18788 18789 --unit cm-1 \
    --sort strength --limit 10
# another isotopologue, as JSON
uv run i2spec lines 632.9 633.0 -i 129I2 --json
# one line with its hyperfine components
uv run i2spec line "R(56) 32-0"
```

Ranges can be given in `nm` (vacuum wavelength), `nm-air`, `cm-1`, `THz` or `MHz`. Line strengths
are integrated cross sections in cm per molecule at the chosen temperature (option `-T`, default
293.15 K).

Every position and every hyperfine component comes with its 1σ, and `i2spec line` lists both. The model can
also be used directly from Python, with the uncertainty of any prediction:

```python
from i2spec.preferred import model
# the fitted model for ¹²⁷I₂
m = model()
# component a10 of R(56) 32-0 and its 1σ, in MHz
m.frequency("R(56) 32-0", "a10")
m.uncertainty("R(56) 32-0", "a10").sigma_MHz
# 1σ of an interval within the line, a15 − a1
m.uncertainty("R(56) 32-0", "a15", "a1").sigma_MHz
```

## Terminal browser

`uv run i2spec tui` shows the absorption spectrum of the selected range, the list of lines, and the
hyperfine components of the selected line.

| key | action |
|---|---|
| `+` / `-` | zoom the spectrum in or out |
| `[` / `]` | move the spectrum window |
| `c` | center the spectrum on the selected line |
| `l` | logarithmic intensity scale |
| `p` | line profile: Doppler-broadened or sub-Doppler (Lamb dips and crossover resonances) |
| `s` | sort the list by strength |
| `a` | include the weak ΔF ≠ ΔJ components |
| `/` | return to the search field |

## Line explorer and documentation

```sh
# build the site and serve it on this computer
uv run i2spec web
# serve it on the local network as well
uv run i2spec web --lan
# re-export the model data after a change to the model
uv run i2spec web export
# write the static site to site/
uv run i2spec web build site
# preview this documentation while editing it
uv run mkdocs serve
```

## Python interface

```python
from i2spec import MHZ_PER_CM, RovibronicModel

m = RovibronicModel("127I2")
# R(56) 32-0, hyperfine-free position in cm⁻¹
nu = m.transition(32, 0, 56, "R")
# 18788.3355 cm⁻¹, 563 260 129 MHz
print(nu, nu * MHZ_PER_CM)

# nu0 and the component offsets in MHz
nu0, comps = m.hyperfine_components(32, 0, 56, "R")
a10 = next(c for c in comps if c.label == "a10")
# 563 260 223.52 MHz (CIPM value: 563 260 223.513 MHz)
print(nu0 + a10.offset)
```

`i2spec.lookup.Catalog` returns line lists at any temperature, `i2spec.spectrum` computes cross
sections and cell transmission, and `i2spec.continuum` computes the bound–free absorption.

## Tests

```sh
# Python test suite
uv run pytest
# line explorer against the Python model
node web/test_physics.mjs
# redraw the documentation figures
uv run --group research python \
    documentation/figures/make_figures.py
```
