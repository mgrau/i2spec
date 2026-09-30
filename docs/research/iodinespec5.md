# Comparison with IodineSpec5

IodineSpec5 (H. Knöckel and E. Tiemann, Institut für Quantenoptik, Leibniz Universität Hannover) is the
program most laboratories use to predict iodine lines. It implements the Hannover potential model
(Knöckel *et al.* 2004; Salumbides *et al.* 2008), the hyperfine interpolation formulae of Salumbides
*et al.* (2006), a local Dunham model for the near-infrared bands (Knöckel *et al.* 2004), and Dunham
constants from the literature outside the range of the potentials. Its source and parameter files are
not published; the program is available from its authors on request.

This note compares i2spec with IodineSpec5 V5.1 (range 11 200–19 430 cm⁻¹). The program's output was
not redistributed: the repository contains only the statistics derived from it and the script that
computes them (`prototypes/iodinespec5_compare.py`), which reads the output files of any copy of the
program.

## Method

IodineSpec5 was run on windows of ±0.15 cm⁻¹ around every line in the precision data sets, at 293 K
with an intensity threshold of 10⁻¹⁴ and with hyperfine structure, which gives 115 322 distinct lines.
Its output was compared with i2spec line by line and, for every data set, the rms deviation of the
measurements from each model was computed on the same rows.

## Reproduction of the published potentials

Where IodineSpec5 uses its potentials (¹²⁷I₂, v′ ≤ 43, v″ ≤ 17), our implementation of the published
2008 potentials reproduces it to 0.48 MHz rms (maximum 2.1 MHz) over 82 735 lines. The difference is
smooth in v″ and J″, as expected from the rounding of the printed coefficients, and is not a coding
error. For ¹²⁹I₂ and ¹²⁷I¹²⁹I the difference is tens of MHz: the isotopologue corrections in IodineSpec5
are not the published functions.

IodineSpec5 uses the Salumbides *et al.* (2006) hyperfine formulae at every v′. i2spec uses those of
Bodermann *et al.* (2002) below v′ = 44; the two differ by a term in J(J+1) of 0.3–0.5 MHz in eQq.
Neither choice is better on every data set; over all 743 hyperfine rows the rms deviation is 0.42 MHz
with the 2006 formulae and 0.45 MHz with those of Bodermann *et al.*

## Data sets

rms of observed − model, MHz, on the rows inside the range of IodineSpec5. i2spec is the current default
parameter set, `i2spec2026i`.

| data set | rows | IodineSpec5 | i2spec |
|---|---|---|---|
| `bipm2003a` (633 nm, three isotopologues) | 158 | 1.98 | **0.90** |
| `bipm2003b` (543 nm) | 30 | 0.46 | **0.02** |
| `bipm2003c` (576 nm) | 11 | 0.75 | **0.10** |
| `bipm2003d` (612 nm, ¹²⁷I₂ and ¹²⁹I₂) | 90 | **1.55** | 2.52 |
| `bipm2003e` (640 nm) | 18 | 0.23 | **0.17** |
| `bipm2005a` (515 nm) | 42 | 0.16 | **0.15** |
| `bipm2012a` (532 nm) | 329 | 0.97 | **0.06** |
| `bodermann1998b` (815 nm) | 4 | 0.09 | **0.02** |
| `bodermann1998c` (780–790 nm) | 7 | 4.07 | **0.24** |
| `bodermann2000a` (778–795 nm) | 32 | 2.14 | **0.09** |
| `cornish2000a` (730 nm) | 2 | 4.59 | **3.47** |
| `dube2004a` (718 nm) | 62 | 0.99 | 1.03 |
| `huet2013a` (716 nm) | 14 | 0.37 | **0.30** |
| `kobayashi2016a` (578 nm) | 76 | 0.88 | **0.14** |
| `liao2010a` (756–810 nm) | 31 | 0.14 | **0.12** |
| `morinaga1989a` (657 nm) | 20 | 1.03 | **0.95** |
| `nishiyama2024a` (520 nm) | 72 | 1.58 | **0.03** |
| `reinhardt2006a` (565–585 nm) | 57 | 0.42 | **0.02** |
| `reinhardt2007a` (735–780 nm) | 7 | 3.24 | **0.17** |
| `sansonetti1997a` (560–656 nm) | 102 | 1.71 | **1.28** |
| `velchev1998a` (572–594 nm) | 115 | 1.96 | **1.83** |
| `xu2000a` (596–655 nm) | 473 | 0.88 | 0.89 |

i2spec is closer to the measurements on 19 of the 22 sets, equal within the scatter on two (`xu2000a`,
`dube2004a`) and worse on one (`bipm2003d`). Both columns are in-sample: the level corrections of i2spec
and the local near-infrared model and hyperfine formulae of IodineSpec5 were fitted to these data. The
held-out figures for i2spec are given with each corrected level (`docs/design/parameter-sets.md`,
0.05–0.8 MHz for levels fixed by several lines).

The sets at 514 nm and beyond 890 nm (`matsunaga2024a`, `yoshiki2023a`, `matyugin2012`,
`nesterenko2019`, `goncharov2007a`, `sakamoto2024a`) lie outside the range of IodineSpec5.

### Near infrared

The near-infrared sets compare a local model with a global one. IodineSpec5 replaces its potentials by a
local Dunham model for the bands v′ = 0 → v″ = 12–17, fitted to the measurements of Bodermann *et al.*
and Liao *et al.* i2spec applies level corrections, and for these bands a correction to the lines
themselves, confined to the measured J″ range (`docs/design/parameter-sets.md`, `i2spec2026f`). On
`liao2010a` the two give 0.12 and 0.14 MHz; allowing for the 114 kHz pressure shift that Liao *et al.*
did not apply to their table, 0.03 and about 0.04 MHz. The held-out rms of i2spec over all 108 precise
near-infrared rows is 0.13 MHz.

Outside the measured J″ range the two predictions differ by up to about 30 MHz (4.9 MHz mean, 10.6 MHz
rms over 141 lines). No measurement decides between them there.

### Isotopologues

`bipm2003d` consists mainly of hyperfine intervals of two ¹²⁹I₂ lines, which is why it weighs heavily in
the table. Over the seven independent isotope-sensitive quantities in `bipm2003a` and `bipm2003d`, i2spec
deviates by 2.1 MHz rms and IodineSpec5 by 4.8 MHz. Neither model's deviations follow a smooth function
of the quantum numbers, and four ¹²⁹I₂ lines cannot constrain a refit of the Born–Oppenheimer
corrections.

## Beyond the range of the potentials

Above v″ = 17 IodineSpec5 uses the X-state Dunham constants of Martin *et al.* (1986) with the B-state
constants of Gerstenkorn and Luc. Its predictions there test the extrapolation of the X potential, for
which no other data existed before the Orsay atlas (`docs/research/orsay-atlas-11000-14000.md`). For 113
lines with J″ ≤ 150 in the near infrared:

| v″ | lines | Dunham − 2008 potential | Dunham − MLR X |
|---|---|---|---|
| 18 | 12 | +30 ± 32 MHz | −7 ± 27 MHz |
| 20 | 13 | +0.22 GHz | −0.16 GHz |
| 22 | 6 | +1.37 GHz | −0.27 GHz |
| 24 | 11 | +5.8 GHz | −0.20 ± 0.11 GHz |
| 26 | 12 | +19.0 GHz | −0.12 ± 0.07 GHz |
| 28 | 7 | +41.1 GHz (1.37 cm⁻¹) | −0.21 ± 0.03 GHz |

The published potential departs from the Dunham levels beyond v″ ≈ 20; the MLR X potential follows them
to 0.1–0.3 GHz, about the accuracy of the Dunham constants themselves.

Above v′ = 43 IodineSpec5 gives only Dunham predictions from Gerstenkorn and Luc, with a stated
uncertainty of ±30 MHz (2σ); elsewhere its manual quotes ±3 MHz for 526–667 nm and the measured
near-infrared bands, and ±90 MHz otherwise. It gives no uncertainty for individual lines.
