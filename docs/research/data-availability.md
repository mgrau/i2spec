# Data: what the model uses, what exists elsewhere, and what is missing

This note describes the measurements behind i2spec, the terms under which they are included, the data
that are known to exist but are not publicly available, and the gaps in the published record.

## What the repository contains

| directory | content | origin |
|---|---|---|
| `data/observations/` | precision line positions and hyperfine intervals, one directory per source | transcribed from the published tables, the BIPM *mises en pratique* and one doctoral thesis; each `meta.toml` gives the source, DOI and transcription notes |
| `data/atlas_lines/` | assigned line positions from three Fourier-transform atlases | Salami & Ross (2005) and the NIST/APO cell spectrum: positions fitted in this work to the published spectra; Orsay atlas (1982): line tables transcribed from the printed volume |
| `data/catalog/` | a bibliography of I₂ data (about 300 entries), with DOIs and the quantity each measures | compiled from the literature |
| `data/external/` | third-party raw files (spectra, supplements) used in fits and tests | **not in the repository**; `data/external/README.md` lists each file, where it comes from and its licence |

The files in `data/observations/` and `data/atlas_lines/` contain measured numbers taken from the
publications cited with each set. Raw spectra and supplementary files are not redistributed; the scripts
that use them expect them in `data/external/`.

The [references page](../../documentation/references.md) of the documentation lists every data set used,
with its DOI.

## Data that exist but are not publicly available

Ranked by what each would change in the model.

1. **The data set of Knöckel, Bodermann & Tiemann (2004).** About 1 500 assigned sub-Doppler frequencies
   (v″ = 0–17, v′ = 0–43, 514–833 nm) to which the published potentials were fitted. The electronic
   supplement named in the paper is no longer available from the publisher. With it, the published fit
   could be reproduced from its input rather than only from its parameters.
2. **The supplement of Salumbides *et al.* (2008).** More than 380 differences between ¹²⁷I₂, ¹²⁹I₂ and
   ¹²⁷I¹²⁹I lines: the only data that constrain the breakdown of the Born–Oppenheimer approximation,
   where the isotope shifts of the model are uncertain by several MHz. Also no longer available from the
   publisher.
3. **The supplement of Salumbides *et al.* (2006).** Per-line hyperfine parameters of the minor
   isotopologues, residuals to the formulae in the paper's tables. Not found on the publisher's site,
   even with institutional access.
4. **Kato *et al.*, *Doppler-Free High Resolution Spectral Atlas of Iodine* (2000).** 526–667 nm at about
   3 MHz, published in print with discs; no digital copy is known.
5. **The comb-referenced Fourier-transform spectrum of Reiners *et al.* (2024).** It reports a
   band-correlated deviation of about 2 m s⁻¹ from the Hannover model over 515–630 nm; the spectrum
   would show whether the deviation lies in the model or the measurement.
6. **The component absorption bands of Tellinghuisen (2011).** The separate A←X, C←X and B←X bound–free
   contributions, published as supplementary material on a host that no longer exists. The model
   reproduces their published sum to 0.5 %.
7. **Fourier-transform spectra of astronomical iodine cells at NIST** (2009–2023) and those of Perdelwitz
   & Huke (2018) and Debus *et al.* (2023): measured, described in publications, not deposited.
8. **The line list of Rodríguez Fernández *et al.* (2023)** (1 204 lines at 14 400–14 600 cm⁻¹): described
   in the paper as supplementary material, but not provided with it.
9. **Scans of the Gerstenkorn & Luc atlases** once served by the Laboratoire Aimé Cotton. They are no
   longer online and were not archived. The printed volumes remain: the 11 000–14 000 cm⁻¹ volume and
   Partie IV (19 700–20 035 cm⁻¹) were photographed and transcribed here
   (`docs/research/orsay-atlas-11000-14000.md`, `orsay-atlas-19700-20035.md`); the assignment tables of
   the 14 800–20 000 cm⁻¹ atlas are still only in print.
10. **IodineSpec line lists.** No output of the program has been published as a data set
    (`docs/research/iodinespec5.md`).

## Gaps in the published record

- **Absolute frequencies at 680–715 nm and 740–750 nm.** Doppler-limited positions now cover 499–702 nm
  (Salami & Ross; NIST/APO) and 709–909 nm (Orsay), but sub-Doppler absolute frequencies are sparse there.
- **Beyond 815 nm.** Apart from the emission lines from v″ = 48–54 at 982–1 068 nm, no absolute
  measurement exists; the Orsay atlas provides Doppler-limited positions to 909 nm.
- **X v″ = 26–47.** No measurement; the levels there are an interpolation of the Morse/long-range
  potential between the atlas levels (v″ ≤ 25) and the emission levels (v″ = 48–54).
- **B v′ > 50.** Three comb-referenced lines (v′ = 52, 53, 62) and the Orsay atlas Partie IV (v′ = 51–79,
  30–100 MHz per line); nothing precise at v′ = 63–79, and the levels above the asymptote are unmeasured.
- **Isotope shifts.** No comb-referenced measurements of ¹²⁹I₂ or ¹²⁷I¹²⁹I beyond the BIPM tables.
- **The absolute cross-section scale.** Spietz *et al.* (2006) and Tellinghuisen (2011) differ by 3 % at
  500 nm, and Saiz-Lopez *et al.* (2004) by 4.7 %; no independent modern measurement exists.

Measurements in any of these regions, or copies of the data sets above, would improve the model
directly. The format for adding a data set is described in `docs/design/observations.md`.
