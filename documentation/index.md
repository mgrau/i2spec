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

Measured line positions are used in two ways only: as the data to which the potentials and level
corrections are fitted, and to estimate the uncertainty of each predicted position.

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

The model combines:

- the X and B potentials of Salumbides *et al.* (2008), a collaboration of the Amsterdam (VU) and
  Hannover groups, referred to here as the 2008 potentials, within the range of their fit
  (X v″ ≤ 17, B v′ ≤ 43);
- Morse/long-range potentials fitted in this work for the levels outside that range;
- empirical corrections to individual levels, fitted to all frequency-comb-referenced measurements, and
  near the B dissociation limit (v′ = 51–79) to the Orsay atlas of Gerstenkorn & Luc (1983);
- a Gaussian process that corrects the B levels without data of their own (v′ = 3–35) from the corrected
  levels around them;
- the published hyperfine parameter formulae, with a table of measured corrections to the B-state
  parameters.

With these, the precision data sets are reproduced to 0.02–0.3 MHz. The uncertainty of positions
elsewhere ranges from below 1 MHz to several GHz, depending on the available data; see
[Uncertainty estimates](uncertainty.md).
