# Installation and usage

i2spec is a Python package managed with [uv](https://docs.astral.sh/uv/). The first `uv run` command
creates the environment and installs the dependencies.

```sh
git clone https://github.com/mgrau/i2spec.git && cd i2spec
uv run i2spec tui          # terminal browser (also: uv run i2spec)
uv run i2spec web          # line explorer and documentation at http://localhost:8777
```

The first query for an isotopologue computes its line list, which takes a few minutes and is then
cached in `~/.cache/i2spec`. The first `i2spec web` also exports the model data used by the line
explorer to `web/data/`.

## Command line

```sh
uv run i2spec lines 532.2 532.3                          # all lines in a range, with uncertainties
uv run i2spec lines 18788 18789 --unit cm-1 --sort strength --limit 10
uv run i2spec lines 632.9 633.0 -i 129I2 --json          # another isotopologue; JSON output
uv run i2spec line "R(56) 32-0"                          # one line with its hyperfine components
```

Ranges can be given in `nm` (vacuum wavelength), `nm-air`, `cm-1`, `THz` or `MHz`. Line strengths
are integrated cross sections in cm per molecule at the chosen temperature (option `-T`, default
293.15 K).

## Terminal browser

`uv run i2spec tui` shows the absorption spectrum of the selected range, the list of lines, and the
hyperfine components of the selected line.

| key | action |
|---|---|
| `+` / `-` | zoom the spectrum in or out |
| `[` / `]` | move the spectrum window |
| `c` | centre the spectrum on the selected line |
| `l` | logarithmic intensity scale |
| `p` | line profile: Doppler-broadened or sub-Doppler (Lamb dips and crossover resonances) |
| `s` | sort the list by strength |
| `a` | include the weak ΔF ≠ ΔJ components |
| `/` | return to the search field |

## Line explorer and documentation

```sh
uv run i2spec web                    # build the site and serve it on this computer
uv run i2spec web --lan              # serve it on the local network as well
uv run i2spec web export             # re-export the model data after a change to the model
uv run i2spec web build site         # write the static site to site/
uv run mkdocs serve                  # preview this documentation while editing it
```

## Python interface

```python
from i2spec import MHZ_PER_CM, RovibronicModel

m = RovibronicModel("127I2")
nu = m.transition(32, 0, 56, "R")       # R(56) 32-0, hyperfine-free position in cm⁻¹
print(nu, nu * MHZ_PER_CM)              # 18788.3355 cm⁻¹, 563 260 129 MHz

nu0, comps = m.hyperfine_components(32, 0, 56, "R")    # nu0 in MHz; offsets in MHz
a10 = next(c for c in comps if c.label == "a10")
print(nu0 + a10.offset)                 # 563 260 223.57 MHz (CIPM value 563 260 223.513 MHz)
```

`i2spec.lookup.Catalog` returns line lists at any temperature, `i2spec.spectrum` computes cross
sections and cell transmission, and `i2spec.continuum` computes the bound–free absorption.

## Tests

```sh
uv run pytest                                            # Python test suite
node web/test_physics.mjs                                # line explorer against the Python model
uv run --group research python documentation/figures/make_figures.py   # redraw the figures
```
