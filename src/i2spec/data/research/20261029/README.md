# Frozen release: October 29, 2026 — the final form

The model rebuilt in its final form on the training data of release 20261021 (`data/training/20261021`, unchanged),
and frozen for the one blind test of the pending data sets. Its three models carry the roles of the earlier ones
(the core fit, the continued model, and the model for every isotopologue and J), but a role resolves to the first
release that has it, so these are reached by id:

```python
i2spec.ResearchModel("research-20261029-i2spec", isotopologue)
```

Release 20261021 stays the default.

## What changed from release 20261021

- **Core fit.** Comb rows are floored at 0.3 MHz (was 0.5). The B limit is tied exactly to D0(X) plus the atomic
  splitting, so De(B) is no free parameter. The six CIPM isotopologue line centers (bipm2003a/d, components reduced
  with the model's own hyperfine) are in the fit, with the B isotope shift (1 − ρ)(1 − y)^5 (a0 + a1 y): they fix the
  split between the potentials' shape and the rotational function that the ¹²⁷I₂ lines leave open. B eqQ has order
  6 (was 5).
- **Continued model.** The B rotational function beyond the core has order 2 (was 3).
- **Named parameters.** The files store Be(X), Be(B), T0 and D0(X) in place of the coefficients they replace; the
  runtime converts them on load (`i2spec.named`, exact to about 1 Hz). `expanded/uncertainty.json` gives their
  values, σ and correlations (`named`).
- **B~1g(1Πg) resonance at B v′ 77–78.** The overlaps of the 1g level with the B levels carry one phase convention
  (positive innermost lobe). The earlier tables took the eigenvectors' signs as they came.
- **Uncertainty.** Rebuilt on this model's own cross-validation (`data/training/20261029`), with the benchmark rows'
  own σ after the model's own hyperfine reduction.

## Inventory

| model | fitted molecular | inherited | calibration |
|---|---|---|---|
| `core` | 63: 46 potential terms, the shift's a0 and a1, 15 B hyperfine terms | 18 | 4 |
| `expanded` | 106: the core's 63, 37 tail terms, 6 hyperfine terms beyond the core | 73 | 5 |
| `i2spec` | 106 | 75: the continued model's 73 and the two nuclear ratios | 5 |

No coefficient is specific to a line or level. 127I2 P(13) 43-0 is left out of every fit.

## Held out, against release 20261021

| check | this release | 20261021 |
|---|---|---|
| comb benchmark, 81 line centers, rms | 0.251 MHz | 0.275 MHz |
| core benchmark, 60 absolute components, rms | 0.184 MHz | 0.156 MHz |
| the same per line (34 lines), rms | 0.162 MHz | 0.174 MHz |
| core benchmark, 578 intraline intervals, rms | 0.049 MHz | 0.059 MHz |
| isotope shift, leave one CIPM line out, rms | 0.086 MHz | 0.078 MHz |

The core benchmark's component rms is higher because of one line, P(49) 24-1, which holds 21 of the 60 components:
B v′ 24 moved by −0.10 MHz in every benchmark row that reaches it (two laboratories), to −0.22 MHz, which is 0.7σ.
Per line, the benchmark improves. The uncertainty model's core terms: σ_HF a = 32 kHz (37), σ_model 0.139 MHz for
lines near precision training rows (0.120) and 0.213 MHz for the rest (0.666).

## Verification

`ResearchModel(id).verify()` replays each model's `reference_predictions.json`. The continued model's 842 centers,
2,586 components and 1,434 intervals come from the fit engine (`i2spec.fit.package`) and are reproduced within
2.4 Hz. The core's and the isotopologue model's references were written by the runtime from these coefficients and
are reproduced exactly.
