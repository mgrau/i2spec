# i2spec

i2spec is an open model of the B³Π(0u⁺) ← X¹Σg⁺ absorption spectrum of molecular iodine. It provides
the position, strength and hyperfine structure of every rovibronic line of ¹²⁷I₂, ¹²⁹I₂ and
¹²⁷I¹²⁹I from the B-state dissociation limit (499 nm) into the near infrared.

All quantities are computed from a Hamiltonian:

- the energy levels of the X and B potential-energy curves, calculated from the radial Schrödinger equation;
- the hyperfine components, from an effective nuclear-spin Hamiltonian;
- the line strengths, from an electronic transition-moment function and the vibrational wavefunctions.

[Open the line explorer :material-arrow-right:](../index.html){ .md-button .md-button--primary }

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

- The X and B states are each represented by a Morse/long-range potential with a rotational function, and the B
  state also by radial functions for its hyperfine parameters. These are fitted together to the precision
  measurements of ¹²⁷I₂ with X v″ ≤ 17 and B v′ ≤ 43, two Fourier-transform atlases and
  {{ fit.n_iso_centers }} line centers of the other isotopologues, and make up the **core fit**. No term is specific
  to an individual level or line.
- The **continuation** extends both curves to their dissociation limits (X v″ ≤ 108, B v′ ≤ 80) with terms that act
  only beyond the core. Near the B limit, second-order sums over the neighboring electronic states enter the
  hyperfine parameters.
- ¹²⁹I₂ and ¹²⁷I¹²⁹I use the same potentials at their own reduced mass, with a fitted isotope shift. Their
  hyperfine functions are scaled by the moments of each nucleus.
- The **uncertainty model** combines the covariance of the fit with the errors of predicting data the fit did not
  see.

[The fit](fit.md) lists which data determine which parameters. Line strengths use the transition moment of
Tellinghuisen (2011) with the model's own vibrational wavefunctions (see [Line intensities and spectra](spectra.md)).

The core fit predicts the {{ fit.benchmark_absolute }} absolute component frequencies of the newest data sets, which
were held out of it, to {{ fit.benchmark_absolute_rms }} MHz rms, and {{ fit.benchmark_intraline }} hyperfine
intervals to {{ fit.benchmark_intraline_rms }} MHz. A locked test of {{ fit.locked_sets }} data sets, kept out of the
fit and of the uncertainty model and scored once against criteria recorded beforehand, found {{ unc.locked_2 }} of
{{ unc.locked_rows }} measurements within 2σ. The absolute frequencies ({{ unc.locked_abs_2 }}) and the isotopologue
measurements met the target; the hyperfine intervals within a line ({{ unc.locked_hf_2 }}) did not (see
[Uncertainty estimates](uncertainty.md#the-locked-test)). The [worked example](calculation.md#worked-example)
follows the 532 nm line R(56) 32–0 through every step of the model.
