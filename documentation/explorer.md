# Line explorer

The [line explorer](../index.html) is a static web page that displays the exported model. `i2spec web export`
writes the line lists, the energies needed for the partition function, the hyperfine components of the lines and
every measurement with its source. From these the page computes line strengths, Doppler profiles, sub-Doppler
spectra and cell transmission for the chosen temperature and cell. These calculations are tested against the Python
model: the partition function agrees to 10⁻⁶, line strengths to 0.2 % and sub-Doppler resonances to 10⁻⁹.

The explorer opens in sub-Doppler mode, ±2 cm⁻¹ around ¹²⁷I₂ P(53) 32–0, unless a shared link sets another view.

## Spectrum and line classes

Each line belongs to one of three classes, with the same color in the spectrum and the table:

- <span class="kind precision"></span>**precision**: measured in a sub-Doppler or frequency-comb data set;
- <span class="kind atlas"></span>**atlas**: measured only in a Fourier-transform atlas (Doppler-limited line center);
- <span class="kind model"></span>**model only**: not measured; the position is the model's prediction.

In the cross-section view the filled **model** area is the cross section of every line. Atlas and precision lines
are drawn over it in their colors, as their share of it, and a diamond below the trace marks each measured line.
The **show** toggles control this:

- **model** includes every line in the spectrum and the table;
- **atlas** and **precision** color their measured lines. When one is off, its lines stay in the model area and the
  table without color;
- with **model** off, the spectrum and table hold only the measured lines of the kinds switched on.

The legend counts the lines in view for each toggle.

The y axis shows the cross section, the transmission of a cell, or the **sub-Doppler** spectrum. The x axis shows the
chosen unit below and wavenumber (or vacuum wavelength) above.

## Hyperfine structure and measurements

Below a span of 0.1 nm each line is drawn as the sum of its main hyperfine components (ΔF = ΔJ), each with the
line's Doppler width and marked by a stick. These patterns leave out the ΔJ = ±2 couplings, which move the
components by less than 1 MHz; the detail panel shows the full calculation. At these spans a click selects the
component under the pointer (within 12 pixels). The table lists every line with a component in view.

Every measured absolute frequency in view is then drawn against the model:

- **a strip under the spectrum** shows measured − model for each measured component (or atlas line center), on its
  own scale in MHz: a dot with the measurement's 1σ bar, in front of a gray box for the model's 1σ. A hollow dot marks
  a source outside the fit. Atlas residuals off the scale of the precision data are drawn as arrows at the edge;
- **on the spectrum, at true scale**, the model's ±1σ is shaded around each measured component and each measurement
  is a dot at its measured frequency, one row per source, once the difference or the σ spans more than a pixel.

Hovering over a dot shows the source, measured − model and both uncertainties, as in the detail panel. For
example, P(13) 43–0 a3, which is kept out of every fit as a diagnostic, lies 5.2 MHz from the model, well outside
the model's 1σ of 1.0 MHz.

## Sub-Doppler spectra

The **sub-Doppler** view draws the saturated-absorption spectrum a pump–probe experiment records, with the model
of `i2spec.saturation`, for views of at most 5 cm⁻¹:

- every hyperfine component gives a **Lamb dip** proportional to S² (the weak-saturation limit), S the component
  strength at the cell temperature; the largest feature in view is scaled to 1;
- each resonance is a Lorentzian of the **homogeneous width Γ** (default 2 MHz);
- **detection**: the saturation *signal*, or the *1f* or *3f* harmonic under frequency modulation of a given
  peak-to-peak width (default 1 MHz). The BIPM frequencies are defined at the zero crossings of the 3f signal;
- **crossovers** appear halfway between two components that share a level. They involve a weak ΔF ≠ ΔJ component.
  They are large at low J and fade to 0.3 % of a dip by J ≈ 55. The export includes weak components down to 10⁻⁴ of
  their line.

A resonance narrower than about three pixels is drawn three pixels wide. To see the true line shape, zoom in until
three pixels span less than Γ.

## Detail panel

The detail panel of the selected line shows:

- its position, strength, uncertainty and data status;
- its hyperfine components, with each component's 1σ for the lines whose full calculation is exported (the
  strongest lines and every measured line);
- every measurement of the line, with its source, stated uncertainty and difference from the model. A difference
  several times larger than both uncertainties marks a disagreement worth a look.

Sources marked “not fitted” are held out for the blind tests, set aside for a reason shown on hovering, or atlases
shown for comparison. Data sets that are still waiting for a blind test are hidden.

Below the measurements, the panel lists the other data on the line and its two levels: hyperfine constants, X-state
hyperfine intervals, level constants, lifetimes and line widths, each with its source and the model's value where
the model has one. A button highlights the line's band v′–v″ in the spectrum and the table, the P branch in dark red
and the R branch in light red. The **band** toggle hides and restores the highlight. With **model**, **atlas** and
**precision** off, the spectrum shows the band alone.

## Export

The **export CSV** buttons download one of two tables:

- every line in view that the toggles keep: label, isotopologue, wavelengths, wavenumber, frequency, strength, E″,
  1σ, domain, class and sources, plus the hyperfine offsets below a span of 0.1 nm;
- the hyperfine components of the selected line.

Each file starts with comment lines that record the model version, the export date, the temperature, the view and a
link that reopens it.

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

On a phone the spectrum takes the upper half of the screen and the line table the lower half. The settings are under
**Options**. One finger drags, two fingers pinch to zoom, and a tap selects.

The address bar holds the whole view, so a view can be shared as a link: the range, isotopologue, temperature,
selected line and component (`comp=a10`), y axis (`y=sigma`, `y=trans`, `y=sub`), sub-Doppler settings (`gamma`,
`det`, `mod`) and highlighted band (`band=32-0`).

The line list loads in parts as each view needs them. Until they arrive, a wide view is drawn from a coarse
spectrum. Line profiles are averaged over each screen pixel. A line's drawn height therefore varies smoothly while
zooming, and its integrated cross section is preserved.

!!! note "Limitations"
    The explorer shows one isotopologue at a time. The bound–free continuum is available in the Python package. The
    crossovers include weak components down to 10⁻⁴ of the line. The few lines whose levels lie outside the model
    are marked *outside the fitted domain*.
