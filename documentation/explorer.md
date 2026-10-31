# Line explorer

The [line explorer](../index.html) is a static web page that displays the exported model. It does
not solve any equations. `i2spec web export` writes the temperature-independent line list, the
energies needed for the partition function, the hyperfine components of selected lines, and all
measurements with their sources. The page computes line strengths, Doppler profiles and cell
transmission for the chosen temperature and cell. These calculations are tested against the Python
model (partition function to 10⁻⁶, line strengths to 0.2 %).

The explorer opens in **sub-Doppler** mode, centered on **¹²⁷I₂ P(53) 32–0** with a
range of **±2 cm⁻¹**. That line is selected in the detail panel. A shared link's
explicit range, display mode and line selection take precedence.

## Classification of lines

Each line belongs to one of three classes, shown in the same color in the spectrum and in the
table:

- <span class="kind precision"></span>**precision**: measured in a sub-Doppler or frequency-comb data set;
- <span class="kind atlas"></span>**atlas**: measured only in a Fourier-transform atlas (Doppler-limited line center);
- <span class="kind model"></span>**model only**: not measured; the position is the model prediction.

A line measured in both kinds of data set is shown as a precision line.

In the cross-section display the filled **model** area is the cross section the model predicts from
every line, measured or not. Over it, the atlas and precision lines are each drawn in their own color
as their share of that cross section, and every measured line is marked by a diamond below the
trace. The three **show** toggles act on this:

- **model** includes every line, so the model area and the table cover all of them;
- **atlas** and **precision** pick out their measured lines in color in the spectrum. With one of them
  off, its lines are still part of the model area and the table, but are not colored or marked in the
  spectrum.

With **model** off, the spectrum and the table keep only the measured lines of the kinds switched on.
The legend gives the number of lines in view for each toggle that is on. While a band is highlighted
(see below), a fourth toggle, **band**, joins them.

The spectrum has two horizontal axes: the chosen unit (normally vacuum wavelength) below, and
wavenumber above (or vacuum wavelength, when the chosen unit is wavenumber). Both carry round tick
values.

## Hyperfine structure in the spectrum

When the displayed range is narrower than 0.1 nm, each line is drawn as the sum of its main hyperfine
components (ΔF = ΔJ), each with the Doppler width of the line, and every component is marked by a
vertical stick whose height is the peak of its own profile. The patterns used for this are computed
for every line without the ΔJ = ±2 couplings of the hyperfine Hamiltonian; this changes the component
offsets by less than 1 MHz (0.9 MHz for R(56) 32–0), which is small compared with the Doppler width.
The detail panel shows the full calculation, including these couplings, where it is available.

Where the hyperfine structure is drawn, a click on the spectrum selects the main component under the pointer
(within 12 pixels), and hovering names it; a click elsewhere selects the line of the nearest component, never a
line whose center lies far from the click. At
these spans the table lists every line with a component in view, even when its center lies outside. The
spectrum marks it, the detail panel gives its wavelength, frequency, offset from the line center, 1σ and
relative strength, and its row is highlighted in the hyperfine table and among the measurements. A row of the
hyperfine table selects its component too. A click on a measurement's dot selects that measurement's component.

## Measurement against the model

Below the same 0.1 nm, every measured absolute frequency in view is drawn against the model, in two ways:

- **A strip under the spectrum** gives measured − model for each measured component (or line center, for an
  atlas line), on its own vertical scale in MHz (GHz when only atlas lines are in view). The strip takes its
  height from the line table, so the spectrum keeps its own. Each measurement is a dot at the model's position of its component,
  with its own 1σ as a vertical bar. The gray box behind it is the model's 1σ for that component. Several
  measurements of one component stand side by side. A hollow dot is a source that no fit uses: a set held
  out for the blind tests, a set set aside, or an atlas shown only for comparison. The vertical scale is set by the precision measurements when there are any. An atlas
  residual (tens of MHz) then lies off the scale and is drawn as an arrow at the edge, but it shows at
  true scale on the spectrum itself.
- **On the spectrum, at true scale**, the model's ±1σ is shaded gray around each measured component, and
  each measurement is a dot at its measured frequency with a ±1σ bar. There is one row per source near
  the top of the plot. Once a dot stands apart from the model's position, a dashed line drops from it
  to the baseline. These marks appear only once the difference, the measurement's 1σ or the model's spans
  more than a pixel: an atlas line's at about 0.1 nm, a comb measurement's typically below 10 MHz across the
  plot. Zooming in on one component shows, for example, the 5.2 MHz miss of P(13) 43–0 a3, well outside
  the model's 1.0 MHz.

The strip needs some height and is left out on a phone held so that the plot is short. Hovering over a
dot gives the line, component and source, measured − model with the measurement's uncertainty, and the
model's 1σ; a click selects the line. The values are those of the detail panel's measurement table.
Both views take the model's 1σ for each component from the full hyperfine patterns (loaded once, about
5 MB). Until they arrive, the line's own 1σ is used.

## Sub-Doppler spectra

The **sub-Doppler** button of the y-axis control draws the saturated-absorption spectrum that a
pump–probe experiment records, using the same model as the Python package (`i2spec.saturation`). It needs the hyperfine patterns, so it is drawn only when
the view spans at most 5 cm⁻¹; a wider view shows the cross section.

- Every hyperfine component gives a **Lamb dip** of amplitude proportional to S², the weak-saturation
  (bilinear) limit, where S is the component strength at the cell temperature. Dips of different lines
  in the view therefore compare as the squares of their strengths. The trace is relative: the largest
  feature in the view is scaled to 1.
- Each resonance is a Lorentzian of the **homogeneous width Γ** (default 2 MHz, as in the terminal
  browser: transit-time and pressure broadening in a typical cell).
- **Detection**: *signal* is the saturation signal itself; *1f* and *3f* are the in-phase first and
  third harmonics under sinusoidal frequency modulation of the given peak-to-peak width (default 1 MHz).
  The BIPM recommended frequencies are defined at the zero crossings of the 3f signal.
- A resonance narrower than about three screen pixels is drawn at that width (the modulation widened in
  proportion, so the shape is kept). Zoom in until three pixels are less than Γ to see the true line shape.
- **Crossovers** appear halfway between two components that share a level. The main (ΔF = ΔJ)
  components never share one, so every crossover involves a weak ΔF ≠ ΔJ component. They are large at
  low J (the strongest is 40% of the strongest dip at J = 2, 5% at J = 10) and fade to 0.3% by J ≈ 55.
  The export carries every weak component above 10⁻⁴ of its line with the numbers of its levels, so
  the explorer draws the same crossovers as the Python package (`saturation.resonances`) except those
  below about 0.3% of a dip, which it leaves out.

When a line is selected in sub-Doppler mode, the hyperfine plot in the detail panel shows its
sub-Doppler spectrum in place of the Doppler-broadened profile. The browser implementation is tested
against the Python model to 10⁻⁹ in the resonance amplitudes (including crossovers, from a full
calculation with level labels) and to 10⁻⁶ in the line shape.

## Highlighting a band

The detail panel of a selected line has a button that highlights its vibrational band v′–v″, or clears the
highlight. The band's lines are then drawn over the rest in two reds: the P branch (ΔJ = −1) dark, the R
branch (ΔJ = +1) light. The B–X system is 0u⁺ ← 0g⁺, so it has no Q branch. This works in the cross section
and the sub-Doppler views, at any span, and the legend counts each branch's lines in view. The band's rows in
the table are tinted to match.

The **band** show toggle turns the highlight off and on without clearing it. With **model** off, it keeps the
band's lines, measured or not: switch off model, atlas and precision to see the band alone, in the spectrum and
the table, or keep atlas or precision on to see the band against those measured lines. Highlighting a band
switches the toggle on.

## Export

The **export CSV** buttons download:

- **lines in view**: every line in the current view that the show toggles keep, not only the 500 listed
  in the table, in the table's order. The columns are the label, isotopologue, vacuum and air
  wavelength, wavenumber, frequency, line strength at the current temperature, E″, the 1σ position
  uncertainty, the domain (fitted, or outside for lines the fit does not yet cover), the class (precision,
  atlas, model only) and the measuring sources; below
  0.1 nm, where they are loaded, the ΔJ = 0 hyperfine offsets and strengths of each line.
- **hyperfine of selected**: the hyperfine components of the selected line as the detail panel shows
  them, with offsets, absolute frequencies, wavelengths and relative strengths.

For ¹²⁷I₂ lines with the full hyperfine calculation (the strongest lines and every measured line), the detail
panel's hyperfine table also gives each component's 1σ: the center's uncertainty combined with the hyperfine
model's own (see [Uncertainty estimates](uncertainty.md)). `i2spec line` gives it for any line.

Each file starts with comment lines (`#`) giving the model version, git revision and
export date of the data, the download time, isotopologue, temperature, view, and a link that
reopens the same view.

## Rendering

Each line profile is averaged over the width of a screen pixel instead of being sampled at the pixel
center, and lines narrower than a pixel are broadened in quadrature to one pixel. The drawn height of
a line therefore varies smoothly as the view is moved or zoomed, and the integrated cross section is
preserved.

The line list is stored in 64 parts per isotopologue, loaded as a view first needs them. A wide view that needs
parts not yet loaded (the whole range, the first time) is drawn at once from a coarse spectrum of every part
(the strength of its lines of each class summed in 64 bins, at seven cell temperatures between which the
explorer interpolates). Each part's own lines replace its coarse bins as it arrives. While the view is zoomed
or moved, the last picture is stretched to the new range and then cross-faded into the redrawn one. The
picture is drawn a quarter of the view's width beyond each side, so a zoom out by one wheel step or a short pan
shows no blank margin before the redraw. A wheel
zoom or an arrow-key move glides to its new range in about a quarter of a second, keeping the point under the
cursor in place.

## Line table

The table lists the lines in order of increasing wavelength. When the table is scrolled beyond its
last line, or the last line is selected and ↓ is pressed, the spectrum window moves to longer
wavelengths and the table continues with the following lines. Scrolling beyond the first line moves
to shorter wavelengths. For measured lines the table gives the sources, each linked to its DOI. The
detail panel lists all measurements of the selected line, each with its own stated uncertainty and its
difference from the model, beside the model's position and 1σ. A difference several times larger than both
uncertainties marks a disagreement worth a look. Sources marked “not fitted” are in no fit: held out for
the blind tests, set aside for a reason the mark gives on hovering, or atlases shown for comparison. The explorer
shows every data set we hold except those still waiting for a blind test, which stay unseen until it is scored.

## Data by line and by level

Below the measurements, the detail panel lists what we hold about the selected line and its two levels beyond line
positions: the effective hyperfine constants a source fitted to this line (of its B level, or B minus X as the
source gives them), X-state hyperfine constants and sublevel intervals, level constants (G or T, B, D, H by v, as
printed, the term value from the source's own origin), and, where held, lifetimes and line widths. Each is given
with its 1σ and source and, where the model has one, its value: its hyperfine constants of the same levels, its
term value and B(v) = [E(v, 2) − E(v, 0)]/6. Data sets still pending a blind test are not shown.

## The model behind each line

The positions, hyperfine patterns, uncertainties and line strengths of all three isotopologues are the model's
(see the [model overview](calculation.md)). Each line's flags give the data status of its levels; the few lines
the model does not hold say *outside the fitted domain*.

## Controls

| action | effect |
|---|---|
| drag the spectrum | move the window (the spectrum follows the pointer) |
| shift-drag across the spectrum | zoom to the selected range |
| mouse wheel over the spectrum | zoom about the cursor position |
| sideways scroll (trackpad) | move the window |
| ← / → or `[` / `]` | move the window by a quarter of its width |
| click on the spectrum | select the nearest line, or, zoomed in far enough, the hyperfine component under the pointer |
| ↑ / ↓ | move through the table |
| `/` | search for a line by its label, e.g. `R(56) 32-0` |

On a phone or narrow window the spectrum takes the upper half of the screen and the line table the
lower half; the settings are in a panel opened with **Options**. On a touch screen, drag with one finger to
move the window, pinch with two fingers to zoom about the point between them, and tap to select the
nearest line. The page does not scroll: the spectrum takes the upper half of the
screen, the line table below it scrolls on its own, and the selected line appears in a bar at the
bottom that opens to the full detail when tapped. The vertical axis is labeled inside the plot, with
round tick values and their power of ten given once above it.

The address bar always contains the current range, isotopologue, temperature and selected line (and component,
`comp=a10`), so that a view can be shared as a link. The range applies as it is typed, on Enter or on leaving
the field. The address also records the y axis (`y=sigma`, `y=trans` or `y=sub`), the sub-Doppler settings
(`gamma` in MHz when not the default, `det=1f` or `det=3f` with `mod` in MHz) and a highlighted band (`band=32-0`).
To go to a wavelength or frequency, type the range or pick a preset; an older link that names a laser (`laser`,
`lunit`) but no range still opens on it.

!!! note "Limitations"
    The line explorer shows one isotopologue at a time and does not include the bound–free continuum,
    which the Python package has. Its crossovers stop at weak components of 10⁻⁴ of the line.
