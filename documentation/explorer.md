# Line explorer

The [line explorer](../index.html) is a static web page that displays the exported model; it solves no equations.
`i2spec web export` writes the temperature-independent line list, the energies needed for the partition function,
the hyperfine components of the lines and every measurement with its source. The page computes line strengths,
Doppler profiles, sub-Doppler spectra and cell transmission for the chosen temperature and cell, and is tested
against the Python model (partition function to 10⁻⁶, line strengths to 0.2 %, sub-Doppler resonances to 10⁻⁹).

It opens in sub-Doppler mode on ¹²⁷I₂ P(53) 32–0, ±2 cm⁻¹ around it. A shared link's range, mode and selection take
precedence.

## Spectrum and line classes

Each line belongs to one of three classes, with the same color in the spectrum and the table:

- <span class="kind precision"></span>**precision**: measured in a sub-Doppler or frequency-comb data set;
- <span class="kind atlas"></span>**atlas**: measured only in a Fourier-transform atlas (Doppler-limited line center);
- <span class="kind model"></span>**model only**: not measured; the position is the model's prediction.

In the cross-section view the filled **model** area is the cross section of every line, measured or not; atlas and
precision lines are drawn over it in their colors as their share of it, and each measured line is marked by a
diamond below the trace. The **show** toggles act on this: **model** includes every line in the spectrum and the
table; **atlas** and **precision** color their measured lines (switched off, those lines stay in the model area
and the table, uncolored). With **model** off, only the measured lines of the kinds switched on remain. The legend
counts the lines in view for each toggle.

The y axis shows the cross section, the transmission of a cell, or the **sub-Doppler** spectrum. The x axis gives the
chosen unit below and wavenumber (or vacuum wavelength) above.

## Hyperfine structure and measurements

Below a span of 0.1 nm each line is drawn as the sum of its main hyperfine components (ΔF = ΔJ), each with the
line's Doppler width and marked by a stick. These patterns omit the ΔJ = ±2 couplings, which move the components by
less than 1 MHz; the detail panel shows the full calculation. At these spans a click selects the component under
the pointer (within 12 pixels), and the table lists every line with a component in view.

Every measured absolute frequency in view is then drawn against the model:

- **a strip under the spectrum** gives measured − model for each measured component (or atlas line center), on its
  own scale in MHz: a dot with the measurement's 1σ bar, in front of a gray box for the model's 1σ. A hollow dot is
  a source that no fit uses. Atlas residuals off the scale of the precision data are drawn as arrows at the edge;
- **on the spectrum, at true scale**, the model's ±1σ is shaded around each measured component and each measurement
  is a dot at its measured frequency, one row per source, once the difference or the σ spans more than a pixel.

Hovering over a dot gives the source, measured − model, and both uncertainties; the values are those of the detail
panel. For example, zooming in on P(13) 43–0 a3 shows a 5.2 MHz miss, well outside the model's 1.0 MHz; that line
is kept out of every fit as a diagnostic.

## Sub-Doppler spectra

The **sub-Doppler** view draws the saturated-absorption spectrum a pump–probe experiment records, with the model
of `i2spec.saturation`, for views of at most 5 cm⁻¹:

- every hyperfine component gives a **Lamb dip** proportional to S² (the weak-saturation limit), S the component
  strength at the cell temperature; the largest feature in view is scaled to 1;
- each resonance is a Lorentzian of the **homogeneous width Γ** (default 2 MHz);
- **detection**: the saturation *signal*, or the *1f* or *3f* harmonic under frequency modulation of a given
  peak-to-peak width (default 1 MHz). The BIPM frequencies are defined at the zero crossings of the 3f signal;
- **crossovers** appear halfway between two components that share a level. They involve a weak ΔF ≠ ΔJ component,
  are large at low J and fade to 0.3 % of a dip by J ≈ 55. Weak components above 10⁻⁴ of their line are exported.

A resonance narrower than about three pixels is drawn at that width: zoom in until three pixels are less than Γ to
see the true line shape.

## Detail panel

The detail panel of the selected line gives its position, strength, uncertainty and data status, its hyperfine
components (with each component's 1σ where the full calculation is exported: the strongest lines and every measured
line), and every measurement of the line with its source, stated uncertainty and difference from the model. A
difference several times larger than both uncertainties marks a disagreement worth a look. Sources marked “not
fitted” are in no fit: held out for the blind tests, set aside for a reason given on hovering, or atlases shown for
comparison. Data sets still waiting for a blind test are not shown.

Below the measurements, the panel lists what is held about the line and its two levels beyond line positions:
hyperfine constants, X-state hyperfine intervals, level constants, lifetimes and line widths, each with its source
and, where the model has one, the model's value. A button highlights the line's band v′–v″ in two reds (P branch
dark, R branch light) in the spectrum and the table; the **band** toggle hides and restores the highlight, and with
**model**, **atlas** and **precision** off it shows the band alone.

## Export

The **export CSV** buttons download either every line in view that the toggles keep (label, isotopologue,
wavelengths, wavenumber, frequency, strength, E″, 1σ, domain, class and sources, and below 0.1 nm the hyperfine
offsets) or the hyperfine components of the selected line. Each file starts with comment lines giving the model
version, the export date, the temperature, the view and a link that reopens it.

## Controls

| action | effect |
|---|---|
| drag the spectrum | move the window |
| shift-drag across the spectrum | zoom to the selected range |
| mouse wheel over the spectrum | zoom about the cursor position |
| sideways scroll (trackpad) | move the window |
| ← / → or `[` / `]` | move the window by a quarter of its width |
| click on the spectrum | select the nearest line, or, zoomed in, the hyperfine component under the pointer |
| ↑ / ↓ | move through the table; past either end the window moves on |
| `/` | search for a line by its label, e.g. `R(56) 32-0` |

On a phone the spectrum takes the upper half of the screen and the line table the lower half; the settings are under
**Options**, one finger drags, two pinch to zoom, and a tap selects. The address bar always holds the range,
isotopologue, temperature, selected line and component (`comp=a10`), y axis (`y=sigma`, `y=trans`, `y=sub`),
sub-Doppler settings (`gamma`, `det`, `mod`) and highlighted band (`band=32-0`), so a view can be shared as a link.

The line list is loaded in parts as a view needs them; until they arrive, a wide view is drawn from a coarse
spectrum. Line profiles are averaged over each screen pixel, so a line's drawn height varies smoothly while zooming
and its integrated cross section is preserved.

!!! note "Limitations"
    The explorer shows one isotopologue at a time and does not include the bound–free continuum, which the Python
    package has. Its crossovers stop at weak components of 10⁻⁴ of the line. The few lines the model does not hold
    say *outside the fitted domain*.
