# Sub-Doppler (saturated-absorption) spectra

*2026-09-15. Code: `src/i2spec/saturation.py`, tests in `tests/test_saturation.py`.*

Every iodine-stabilized laser standard is measured this way, and the BIPM recommended frequencies are defined for it, so the model has to produce Doppler-free spectra, not only Doppler-limited ones. This module is Phase A: it predicts *where* the resonances are and *roughly how strong*, from the hyperfine components we already compute.

## 1. What happens in the cell

A counter-propagating pump and probe interrogate one velocity class at a time. Two kinds of narrow feature sit inside the Doppler profile of a line:

* **Lamb dips** at every hyperfine component, from the v_z = 0 class, which both beams address.
* **Crossovers** halfway between two components that **share a level**, from the two velocity classes v_z = ±(ν_j − ν_i)/2, where the pump is resonant with one component and the probe with the other. Two components that share no level never produce one.

The velocity class a crossover needs is off the centre of the Maxwell distribution, so crossovers fade as the pair separation approaches the Doppler width. That is what makes the list finite.

## 2. The amplitude model

In the weak-saturation rate-equation limit the signal is bilinear in the pump and probe transition strengths (Letokhov & Chebotayev 1977; Demtröder, *Laser Spectroscopy*, §2.3, §7.2): the pump burns a hole proportional to S_pump, and the probe reads it with a cross section proportional to S_probe. With S the relative strengths our `hyperfine.line_components` already returns:

```
dip(i)        = S_i**p
crossover(ij) = 2 · w_ij · (S_i · S_j)**(p/2) · exp(−ln2 · Δ_ij² / Δν_D²)
```

* **p = `exponent`, default 2.** The bilinear limit. Strong saturation flattens the dependence towards p = 1; `exponent=1` reproduces the convention of the Salumbides 2006 calculated traces (§5).
* **The exponential** is the fraction of molecules in the velocity class the crossover uses, taken from the Gaussian Doppler profile of FWHM Δν_D at detuning Δ_ij/2.
* **The factor 2** counts the two velocity classes, ±v_z, that land on the same laser frequency.
* **Amplitudes are relative and sum to 1.** The module gives the shape of a Doppler-free spectrum, not its depth.
* For i = j the crossover formula reduces to the dip formula, so the geometric mean is the natural generalization: the whole set is the outer product S ⊗ S restricted to level-sharing pairs.

### Every assumption, stated

| # | Assumption | Consequence if wrong |
|---|---|---|
| A1 | Two-level rate equations, weak saturation, no coherent (recoil, Dicke, crossover-dip interference) effects | Amplitudes only; positions are unaffected |
| A2 | Levels are identified by (I, F), with I the **dominant** nuclear spin from `HyperfineLevel` | At high J the quadrupole mixes I; two levels of one F with similar I content could be merged or split wrongly, creating or missing a crossover |
| A3 | V-type (shared lower) and Λ-type (shared upper) get the same weight, `lambda_weight=1` | Λ-type should in fact be weaker: it works through population accumulated in the B state, which lives about 1 µs, against a ground-state hole that survives the transit time, about 5 µs for a millimetre beam at 300 K. A weight near 0.2 is the physical expectation; it is **not** applied by default, because it is an estimate, not a measurement |
| A4 | No optical pumping between hyperfine levels: a molecule that absorbs is lost from its ground level, and the B state fluoresces to many v″ | Holds for I₂, unlike alkalis, where optical pumping inverts some crossovers |
| A5 | One homogeneous width for every resonance, supplied by the caller | Crossover widths are the mean of the two component widths, which is the same number here |
| A6 | One rovibronic line at a time | Two different lines that share a level (for example R(J) and P(J+2), which share the upper level) can also produce crossovers. Their separation is normally far beyond the Doppler width, so the velocity factor would kill them, but the module never generates them |
| A7 | No absolute depth, saturation broadening or power dependence | The caller sets the width; depths are relative |

### A structural consequence, worth knowing

**Main components never produce crossovers with each other.** The main (ΔF = ΔJ) components take one upper and one lower level each, so no two of them share a level — `test_main_components_never_produce_crossovers` asserts it. Every crossover therefore involves at least one weak ΔF ≠ ΔJ component, whose strength at these J is 10⁻²–10⁻⁴ of a main component. With p = 2 the crossovers land 10⁻²–10⁻³ below the dips. This is why component-only simulations describe I₂ Doppler-free spectra as well as they do.

## 3. Line shape and 3f detection

`lineshape(x, fwhm, harmonic=n, modulation=m)` returns

S_n(x) = (2/π) ∫₀^π L(x + (m/2) cos θ) cos(nθ) dθ,

the in-phase n-th Fourier component of a Lorentzian under sinusoidal frequency modulation of peak-to-peak width m. For m ≪ FWHM it tends to the n-th derivative. The BIPM mise en pratique specifies **third-harmonic detection** with about 1 MHz peak-to-peak modulation, and defines each recommended frequency at the **zero crossing** of that 3f signal, which is what `harmonic=3` reproduces. `signal()` interpolates the modulated shape from an internal kernel of step FWHM/256, accurate to about 10⁻⁵ of the peak even across that zero crossing.

`transit_time_width()` gives the width floor in a cell: 0.12 MHz for a 1 mm beam at 300 K.

## 4. What the Hannover lineage does

* The software survey lists "Doppler/sub-Doppler spectra" among IodineSpec's features (`01-software-survey.md` §2.4). IodineSpec is closed, and its treatment of crossovers is not documented publicly, so there is nothing to compare against directly.
* Salumbides et al. (2006) plot calculated Doppler-free profiles in their Figs. 3–5. An earlier digitization of those figures (`isotopologue-hyperfine.md` §4.1b, §8d) fitted them as a **plain sum over components** with a Lorentzian-like profile (their eq. 10 with b = 0) and intensity ∝ S¹. The fit reproduced every resolved peak to 0.1 MHz rms, so their calculated curves contain **no crossovers** and use p = 1.

## 5. Validation against Salumbides 2006, Figs. 3–5

The digitized observed ("black") traces were tested against this model, using the shift and profile width from that earlier fit. Positions are matched within 3 MHz, the scatter of the observed traces.

| Figure | Line | Components (main) | Crossovers | Above 1% of the strongest dip |
|---|---|---|---|---|
| 3 | ¹²⁷I₂ P(34) 11-3 | 78 (15) | 428 | 0 |
| 4 | ¹²⁹I₂ R(34) 11-3 | 204 (28) | 1848 | 4 |
| 5 | ¹²⁷I¹²⁹I R(37) 11-3 | 570 (48) | 8748 | 4 |

**Predicted size of the effect**, as the largest change to the simulated profile at each figure's resolution:

| Figure | p = 2 (weak saturation) | p = 1 (the Salumbides plotting convention) |
|---|---|---|
| 3 | 1.0% of peak | 14.4% |
| 4 | 3.0% | 37.0% |
| 5 | 1.9% | 30.9% |

**What the traces show.** Nearly every observed peak sits on a component: 18 of 18 in Fig. 4 and 34 of 37 in Fig. 5. The leftovers are wiggles of prominence 0.021–0.042 of the trace span, i.e. digitization noise.

**The match rate against crossovers is not evidence, and the script says so.** With hundreds to thousands of crossovers spread over about 1 GHz, the ±3 MHz windows around them cover 75%, 97% and 100% of the axis in Figs. 3, 4 and 5. A random position would "match" at that rate. Restricted to crossovers strong enough to be visible (≥ 1% of the strongest dip), the windows cover 0–3% of the axis, and **no leftover peak falls in one**.

**Conclusion.**
* These figures cannot confirm the crossover model: at p = 2 it predicts 1–3% features, below the traces' noise.
* They do disfavour p = 1 combined with equal V and Λ weights, which would distort the profiles by 14–37% — the paper's own component-only calculation would then not fit its own observations as well as it does.
* The positions are what this module is currently good for. The amplitudes need a spectrum with known saturation parameter and cell conditions to test.

## 6. How to validate it properly

1. A published Doppler-free I₂ spectrum with **resolved** crossovers, recorded at a stated saturation parameter, beam diameter and pressure. The BIPM mise en pratique traces are 3f and heavily processed, so they test positions and the lock point rather than amplitudes.
2. Compare at two saturation powers: p should move from 2 towards 1 as the dips saturate.
3. A line at low J, where the ΔF ≠ ΔJ components are relatively strong and crossovers should be largest.
4. A Λ-type-dominated pair, to measure `lambda_weight` instead of assuming it.

## 7. API

```python
from i2spec import RovibronicModel
from i2spec.saturation import line_resonances, signal
import numpy as np

model = RovibronicModel("127I2")
nu0, res = line_resonances(model, 32, 0, 56, "R", 300.0)   # R(56) 32-0, the 532 nm standard
nu = np.linspace(-800, 500, 20001)                          # MHz from the hyperfine-free centre
lamb = signal(nu, res, 1.0)                                 # Lamb dips and crossovers, 1 MHz wide
lock = signal(nu, res, 1.0, harmonic=3, modulation=1.0)     # what a 3f-locked laser sees
```

Each `Resonance` carries its `offset`, relative `amplitude`, `kind` ("dip" or "crossover"), `sharing` ("lower", "upper" or None), the indices of the components involved and their labels, so a caller can pick out, say, the a10 dip or ask which components make a given crossover.

## 8. Not implemented

* Absolute signal depth, saturation broadening and the power dependence of the dips.
* Crossovers between different rovibronic lines that share a level (A6).
* Coherent effects: recoil splitting, which is about 2 kHz at 532 nm and matters at the BIPM's level, and crossover-dip interference.
* Optical pumping between hyperfine levels (A4), which is negligible for I₂ but not for the alkalis this machinery resembles.
