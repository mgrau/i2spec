# i2spec

i2spec is an open model of the B³Π(0u⁺) ← X¹Σg⁺ absorption spectrum of molecular iodine. It gives
the position, strength and hyperfine structure of every rovibronic line of ¹²⁷I₂, ¹²⁹I₂ and
¹²⁷I¹²⁹I between the B-state dissociation limit (499 nm) and the near infrared, and an uncertainty for
each position.

All quantities are computed from a Hamiltonian:

- level energies from the radial Schrödinger equation with potential-energy curves for X and B;
- hyperfine components from an effective nuclear-spin Hamiltonian;
- line strengths from an electronic transition-moment function and the vibrational wavefunctions.

[Open the line explorer :material-arrow-right:](../index.html){ .md-button .md-button--primary }

Measured line positions are used in two ways only: as the data to which the potentials and hyperfine
functions are fitted, and to test and calibrate the uncertainty of each predicted position.

<div class="grid cards" markdown>

-   :material-console:{ .lg } **[Installation and usage](getting-started.md)**

    Command line, terminal browser, local web server and Python interface.

-   :material-chart-bell-curve:{ .lg } **[Line explorer](explorer.md)**

    Search a wavelength or frequency range; plot the absorption cross section or the transmission
    of a cell; list the lines with their hyperfine components, uncertainties and measurements.

-   <span class="twemoji lg"><svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 24 24" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"><path fill="none" d="M4 3C4.6 14 6.4 20 9.2 20C12.4 20 13.6 8 21 7.5"/><path fill="none" d="M8.5 15.5H9.9M7.7 11.5H12.2"/></svg></span> **[Model](calculation.md)**

    Hamiltonians, potentials, numerical methods, intensities and uncertainty estimates.

-   :material-book-open-variant:{ .lg } **[References](references.md)**

    All measured data sets used by the model, with links to the original publications.

</div>

## The model

The model has seven parts:

- One Morse/long-range potential for each of X and B, with a fitted rotational function.
- Radial functions for the B-state hyperfine parameters.
- One least-squares fit of the potentials and hyperfine functions, the core fit, made to every precision
  measurement of ¹²⁷I₂ and two Fourier-transform atlases with X v″ ≤ 17 and B v′ ≤ 43. It has no correction to
  any individual level or line.
- A continuation of both curves to their dissociation limits (X v″ ≤ 108, B v′ ≤ 80): a few dozen terms
  that act only beyond the core. They are fitted to the Orsay atlases, the X level constants of Martin
  *et al.*, emission lines and the comb lines above B v′ = 43.
- Second-order hyperfine sums for the B levels near their dissociation limit.
- The other isotopologues on the same potentials at their own reduced mass, with a Born–Oppenheimer correction
  of two coefficients fitted to the measured isotope shifts, and the hyperfine functions scaled by the moments
  of each nucleus.
- An uncertainty model with three parts: the covariance of the fit, the error of predicting a level from
  the data around it (fitted to held-out measurements), and the same for hyperfine intervals.

The core fit predicts 60 absolute component frequencies of newer data sets, held out of it, to 0.16 MHz rms,
and 592 hyperfine intervals to 0.06 MHz. Beyond the core, a level predicted from its neighbors alone reproduces
comb lines to 1.4 MHz (B v′ = 44–50), atlas lines to 22–80 MHz and the X level constants to 44–122 MHz. The
isotope shifts of ¹²⁹I₂ and ¹²⁷I¹²⁹I are predicted, each line from the others, to 0.08 MHz. On 35 data sets
scored once against predictions and uncertainties frozen beforehand, 90 % of 345 measurements lie within 2σ.
See [Uncertainty estimates](uncertainty.md).

The model gives every line of the three isotopologues whose levels it holds (J″ ≤ 445). Line strengths use the
transition moment of Tellinghuisen (2011) with the model's own vibrational wavefunctions, which reproduce
measured cell spectra better than those of earlier potentials (see [Line intensities and spectra](spectra.md)).
The explorer lists each measurement of a line beside the model, with both uncertainties.
