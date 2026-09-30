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

On a phone or narrow window the spectrum takes the upper half of the screen and the line table the
lower half; the settings are in a panel opened with **Options**. On a touch screen, drag with one finger to
move the window, pinch with two fingers to zoom about the point between them, and tap to select the
nearest line. The page does not scroll: the spectrum takes the upper half of the
screen, the line table below it scrolls on its own, and the selected line appears in a bar at the
bottom that opens to the full detail when tapped. The vertical axis is labelled inside the plot, with
round tick values and their power of ten given once above it.

The address bar always contains the current range, isotopologue, temperature and selected line, so
that a view can be shared as a link.

!!! note "Limitations"
    The line explorer shows one isotopologue at a time and does not include sub-Doppler spectra or
    the bound–free continuum. Both are available in the Python package and the terminal browser.
