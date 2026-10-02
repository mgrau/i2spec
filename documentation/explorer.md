# Line explorer

The [line explorer](../index.html) is a static web page that displays the exported model. It does
not solve any equations. `i2spec web export` writes the temperature-independent line list, the
energies needed for the partition function, the hyperfine components of selected lines, and all
measurements with their sources. The page computes line strengths, Doppler profiles and cell
transmission for the chosen temperature and cell. These calculations are tested against the Python
model (partition function to 10⁻⁶, line strengths to 0.2 %).

## Classification of lines

Each line belongs to one of three classes, shown in the same colour in the spectrum and in the
table:

- <span class="kind precision"></span>**precision**: measured in a sub-Doppler or frequency-comb data set;
- <span class="kind atlas"></span>**atlas**: measured only in a Fourier-transform atlas (Doppler-limited line centre);
- <span class="kind model"></span>**model only**: not measured; the position is the model prediction.

A line measured in both kinds of data set is shown as a precision line. In the cross-section display
the contributions of the three classes are shown as stacked areas. The three **show** toggles (model,
atlas, precision) switch each class on or off in both the spectrum and the table.

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

## Sub-Doppler spectra

The **sub-Doppler** button of the y-axis control draws the saturated-absorption spectrum that a
pump–probe experiment records, using the same model as the Python package (`i2spec.saturation`,
derived in `docs/research/sub-doppler.md`). It needs the hyperfine patterns, so it is drawn only when
the view is narrower than 0.1 nm; a wider view shows the cross section with a note to zoom in.

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
  proportion, so the shape is kept); the note in the plot gives the width drawn. Zoom in until it
  equals Γ to see the true line shape.
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

## Laser wavelength and harmonics

The **laser fundamental** field takes the wavelength or frequency of a laser (vacuum or air
wavelength in nm, cm⁻¹, THz or MHz) and a harmonic order n = 1–4, for example 1064.49 nm and ×2 for a
frequency-doubled Nd:YAG laser, 1542 nm and ×3 or ×4 for a telecom laser, or 1319, 1111 or 1156 nm
and ×2. Enter (or **Go**) shows ±0.1 nm of vacuum wavelength around the harmonic (±0.04 nm in
sub-Doppler mode). Harmonic generation multiplies the optical frequency by n, so an air wavelength is
converted to vacuum at the fundamental before it is multiplied; the harmonic of 1064 nm in air is not
532 nm in air.

The harmonic is marked in the spectrum by a dotted line labelled *laser ×n*, and the status bar gives
both values, for example `laser 2 × 281.630263 THz = 563.260526 THz (1064.4895 nm → 532.2447 nm vac)`.
For the selected line the detail panel adds its fundamental-equivalent wavelength and frequency (÷ n)
and its distance from the laser harmonic, both at the harmonic and at the fundamental: the tuning the
fundamental needs to reach the line. Clearing the field removes the laser.

## Export

The **export CSV** buttons download:

- **lines in view**: every line in the current view that the show toggles keep, not only the 500 listed
  in the table, in the table's order. The columns are the label, isotopologue, vacuum and air
  wavelength, wavenumber, frequency, line strength at the current temperature, E″, the 1σ position
  uncertainty, the class (precision, atlas, model only) and the measuring sources. With a laser set,
  the fundamental-equivalent wavelength and frequency and the distance from the laser harmonic are
  added; below 0.1 nm, where they are loaded, the ΔJ = 0 hyperfine offsets and strengths of each line.
- **hyperfine of selected**: the hyperfine components of the selected line as the detail panel shows
  them, with offsets, absolute frequencies, wavelengths and relative strengths.

Each file starts with comment lines (`#`) giving the model version, git revision and
export date of the data, the download time, isotopologue, temperature, view, laser, and a link that
reopens the same view.

## Rendering

Each line profile is averaged over the width of a screen pixel instead of being sampled at the pixel
centre, and lines narrower than a pixel are broadened in quadrature to one pixel. The drawn height of
a line therefore varies smoothly as the view is moved or zoomed, and the integrated cross section is
preserved.

## Line table

The table lists the lines in order of increasing wavelength. When the table is scrolled beyond its
last line, or the last line is selected and ↓ is pressed, the spectrum window moves to longer
wavelengths and the table continues with the following lines. Scrolling beyond the first line moves
to shorter wavelengths. For measured lines the table gives the sources, each linked to its DOI. The
detail panel lists all measurements of the selected line with their uncertainties and their
differences from the model.

## Comparison with the 2008 model

For every line the detail panel also gives the position of the published model of Salumbides *et al.*
(2008), whose potentials are those of the program IodineSpec, as its difference from this model.
It is computed here from the published potentials and Born–Oppenheimer corrections, and the measured
components are compared with the published hyperfine formulae alone, so the measurements table has a
second residual column, *obs − 2008*. Inside the range the 2008 potentials were fitted to (X v″ ≤ 17,
B v′ ≤ 43) the two models agree to a few MHz; outside it the 2008 curves are extrapolations, and the
detail panel says so: there the difference reaches tens of GHz. The line-list CSV export carries the same
difference as a column.

## Controls

| action | effect |
|---|---|
| drag the spectrum | move the window (the spectrum follows the pointer) |
| shift-drag across the spectrum | zoom to the selected range |
| mouse wheel over the spectrum | zoom about the cursor position |
| sideways scroll (trackpad) | move the window |
| ← / → or `[` / `]` | move the window by a quarter of its width |
| click on the spectrum | select the nearest line |
| ↑ / ↓ | move through the table |
| `/` | search for a line by its label, e.g. `R(56) 32-0` |
| `L` | enter a laser fundamental |

On a phone or narrow window the spectrum takes the upper half of the screen and the line table the
lower half; the settings are in a panel opened with **Options**. On a touch screen, drag with one finger to
move the window, pinch with two fingers to zoom about the point between them, and tap to select the
nearest line. The page does not scroll: the spectrum takes the upper half of the
screen, the line table below it scrolls on its own, and the selected line appears in a bar at the
bottom that opens to the full detail when tapped. The vertical axis is labelled inside the plot, with
round tick values and their power of ten given once above it.

The address bar always contains the current range, isotopologue, temperature and selected line, so
that a view can be shared as a link. It also records the y axis (`y=trans` or `y=sub`), the
sub-Doppler settings (`gamma` in MHz when not the default, `det=1f` or `det=3f` with `mod` in MHz) and
the laser (`laser`, its unit `lunit` and the harmonic `n`). A link with a laser but no range opens on
the laser's harmonic.

!!! note "Limitations"
    The line explorer shows one isotopologue at a time and does not include the bound–free continuum,
    which the Python package has. Its crossovers stop at weak components of 10⁻⁴ of the line.
