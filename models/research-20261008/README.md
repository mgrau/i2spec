# Research model frozen October 8, 2026

This snapshot preserves the completed candidate with the lowest absolute-component
RMS in the current accepted research results: **0.155929392 MHz** on 60 distinct
components from 34 rovibrational transitions. Its architecture is X14/B20 MLR with
anchored radial X hyperfine functions. It has **61 fitted molecular coefficients,
11 inherited molecular constants, and 4 separately fitted atlas-calibration
coefficients**, with no individual line shifts or vibrational-level offsets.

The absolute benchmark has been reused for architecture comparisons and is
exploratory validation. This candidate has not completed the seven band/source
holdouts. Its RMS on 592 held-out intraline intervals is **0.108003070 MHz**,
compared with **0.048595155 MHz** for the existing 62-fitted/18-inherited reference
(absolute RMS 0.164413544 MHz). The latter has completed its seven holdouts.
These models trade absolute error, intraline error, coverage and coefficient
count. “Best” here refers specifically to the lowest absolute-component RMS
among completed, accepted, training-selected candidates at the snapshot time.
The production model and the explorer continue to use `i2spec2026o`.

The domain is **127I2, X v<=17 and B v<=43**. The exact `127I2 P(13) 43-0`
exclusion is retained in potential fitting, hyperfine fitting and headline
scores. Predictions for that line remain diagnostic only.

## Evaluate or verify the snapshot

From the repository root, install the locked environment and run:

```sh
uv sync --locked
uv run python models/research-20261008/predict.py 'P(53) 32-0'
uv run python models/research-20261008/predict.py 'R(56) 32-0' --component a10
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1 \
  uv run python models/research-20261008/predict.py --verify
```

The evaluator uses the repository's B-spline solver and angular hyperfine
Hamiltonian directly. It needs none of the ignored research scripts or output
folders. The verification compares all 60 absolute-component predictions and
592 intraline predictions with the frozen source outputs, requiring agreement
within 100 Hz. It verifies reproduction, not independent predictive accuracy.

## Files and conventions

| File | Contents |
|---|---|
| `model.json` | All fitted and inherited molecular constants, radial grids, units and original optimizer coefficients |
| `predict.py` | Standalone snapshot evaluator using the installed `i2spec` package |
| `reference_predictions.json` | Frozen absolute and intraline predictions, with observation identities |
| `selection.json` | Exact terminal fit tags, training objectives, benchmark errors, coefficient inventory and selection caveats |
| `audits.json` | Original dense shape, radial mesh, actual-J, quadrature and hyperfine checks |
| `provenance.json` | Source artifact and implementation SHA256 hashes and the base Git revision |

Potentials use the standard MLR form from `i2spec.potentials.MLRPotential`,
with a Chebyshev exponent. For X, the exponent coordinate is `(y_q - 0.2)/0.6`;
for B, it is `y_q`. The rotational factor is
`1 + 2 Re/(R+Re) sum_i alpha_i T_i(((R-Re)/(R+Re) - 0.3)/0.4)`.
No adiabatic correction is applied to the reference isotopologue.

Hyperfine polynomials use powers of
`z = ((R^2-Re^2)/(R^2+Re^2))/scale`, with scale 0.15 for X and 0.35 for B,
averaged over the corresponding radial wavefunction. X eqQ and C are anchored
at X(v=0,J=0): each is its inherited anchor plus the fitted polynomial's
expectation minus its ground-state expectation. X d and delta are inherited
constants. B couplings are direct radial polynomial expectations. The angular
Hamiltonian includes J±2 mixing and averages couplings between rotational
levels, as in the fitted implementation.

The evaluator uses full core-domain radial boxes (X 2.10–4.20 Å and B
2.35–6.10 Å). The original source selected smaller grid tiers where possible;
the stored reference comparison checks the resulting numerical agreement.
