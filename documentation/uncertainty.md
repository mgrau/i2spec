# Uncertainty estimates

Each line position is given with a standard (1σ) uncertainty. The uncertainty is not derived from the
covariance of a fit. It is estimated from the agreement of the model with measurements that were not
used to determine the relevant parameters, and it depends on how the upper and lower levels of the
line are determined (Table 1).

| levels of the line | uncertainty | basis |
|---|---|---|
| both levels corrected by frequency-comb data, within 15 in J of the data | 0.05–1 MHz | held-out rms of each level |
| NIR bands v′ = 0 → v″ = 12–17, within the measured J range | 0.13 MHz | held-out rms of the band correction |
| X v″ = 18–25 within the J range of the Orsay atlas | 15–50 MHz | held-out rms of each level |
| range of the 2008 fit, without level corrections | 3–8 MHz | comparison with BIPM and comb data by region |
| B v′ = 44–50 | 15 MHz | comparison with Matsunaga (2024) and Yoshiki (2023) |
| X v″ = 26–75 | 150 MHz | the X potential fitted to the level constants of Martin *et al.* (1986), 55–73 MHz rms |
| X v″ = 76–89 | 400 MHz | the same, 174 MHz rms |
| X v″ > 89 | 2.1 THz | not physical: the X potential is fitted to v″ = 89 only (no B–X line in the list reaches these levels) |
| B v′ = 51–79 within the J range of the Orsay atlas Partie IV | 20–250 MHz | held-out rms of each level |
| B v′ > 50 otherwise | 1 GHz | the refitted B potential alone, 150–350 MHz from the atlas levels, extrapolated |

*Table 1. Uncertainty of the hyperfine-free line position.*

<figure class="fig" markdown="span">
<div class="svg" data-svg="figures/uncertainty_map.svg" role="img" aria-label="Map of position uncertainty over upper and lower vibrational quantum numbers"></div>
<figcaption markdown="span">**Figure 1.** Uncertainty of the R(50) line of each band v′–v″. Below v″ = 17 and v′ = 43 the 2008 potentials and the level corrections give a few MHz or better; the darkest rows are levels with direct frequency-comb measurements. At v″ = 18–25 the Orsay atlas determines the levels to 15–50 MHz where it observed them, and v″ = 48–54 are determined by emission measurements. Elsewhere the positions are extrapolations of the potentials.</figcaption>
</figure>

The uncertainty of the hyperfine offsets is given on the [Hyperfine structure](hyperfine.md) page.
Where the local near-infrared model of Knöckel, Bodermann & Tiemann (2004) applies, the command line
and the line explorer indicate this.
