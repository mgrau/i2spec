# i2spec

i2spec is an open model of the B³Π(0u⁺) ← X¹Σg⁺ absorption spectrum of molecular iodine. It gives
the position, strength and hyperfine structure of every rovibronic line of ¹²⁷I₂, ¹²⁹I₂ and
¹²⁷I¹²⁹I between the B-state dissociation limit (499 nm) and the near infrared, and an uncertainty for
each position.

All quantities are computed from a Hamiltonian:

- level energies from the radial Schrödinger equation with potential-energy curves for X and B;
- hyperfine components from an effective nuclear-spin Hamiltonian;
- line strengths from an electronic transition-moment function and the vibrational wavefunctions.

[Open the spectrum :material-arrow-right:](../index.html){ .md-button .md-button--primary }

Measured line positions are used in two ways only: as the data to which the potentials and hyperfine
functions are fitted, and to test and calibrate the uncertainty of each predicted position.

<div class="grid cards" markdown>

-   :material-chart-bell-curve:{ .lg } **[Line explorer](../index.html)**

    Search a wavelength or frequency range; plot the absorption cross section or the transmission
    of a cell; list the lines with their hyperfine components, uncertainties and measurements.

-   :material-console:{ .lg } **[Installation and usage](getting-started.md)**

    Command line, terminal browser, local web server and Python interface.

-   :material-atom:{ .lg } **[Model](calculation.md)**

    Hamiltonians, potentials, numerical methods, intensities and uncertainty estimates.

-   :material-book-open-variant:{ .lg } **[References](references.md)**

    All measured data sets used by the model, with links to the original publications.

</div>

## The model

For ¹²⁷I₂ the model is **balanced-expanded**. It has six parts:

- One Morse/long-range potential for each of X and B, with a fitted rotational function.
- Radial functions for the B-state hyperfine parameters.
- One least-squares fit of the potentials and hyperfine functions, the core fit, made to every precision
  measurement and two Fourier-transform atlases with X v″ ≤ 17 and B v′ ≤ 43. It has no correction to any
  individual level or line.
- A continuation of both curves to their dissociation limits (X v″ ≤ 108, B v′ ≤ 80): a few dozen terms
  that act only beyond the core. They are fitted to the Orsay atlases, the X level constants of Martin
  *et al.*, emission lines and the comb lines above B v′ = 43.
- Second-order hyperfine sums for the B levels near their dissociation limit.
- An uncertainty model with three parts: the covariance of the fit, the error of predicting a level from
  the data around it (fitted to held-out measurements), and the same for hyperfine intervals.

The core fit predicts 60 absolute component frequencies of newer data sets, held out of it, to 0.16 MHz rms,
and 592 hyperfine intervals to 0.06 MHz. Beyond the core, a level predicted from its neighbors alone reproduces comb lines to
1.4 MHz (B v′ = 44–50), atlas lines to 22–80 MHz and the X level constants to 44–122 MHz. On 36 data sets kept out of every fit, the model ties or beats the
earlier i2spec2026o on lines nobody had measured before. i2spec2026o is better where it had been corrected to
measurements of the same lines.

¹²⁹I₂ and ¹²⁷I¹²⁹I, and the few ¹²⁷I₂ lines beyond the domain, use the earlier model, i2spec2026o. It is
built from four pieces:

- the potentials of Salumbides *et al.* (2008), a collaboration of the Amsterdam (VU) and Hannover groups;
- Morse/long-range curves beyond their range;
- empirical corrections to individual levels;
- the published hyperfine formulae.

All line strengths come from i2spec2026o's transition moment and wavefunctions. Every line says which model
gave it, and the explorer lists each measurement beside the model with both uncertainties. See
[Uncertainty estimates](uncertainty.md).
