# Frozen prospective validation evidence

`prospective_v1/` contains immutable protocol, model fingerprints and 192
point forecasts for 48 component targets on 16 lines. Its manifest SHA256 is
`e0cdfbaba8fc4758a1700edfeceff9b92b0af76a56b737b081a9a56dfa9a32f4`.
`plan_v1.json` is the preparation input; `measurements_template.json` is a
blank input template to COPY before filling it. `lab_targets.csv` lists the
balanced predictions and vacuum wavelengths for laboratory planning.

No new measurements were available. The prospective pending report scores
zero outcomes. `completion_state.json` distinguishes computational completion
from experimental evidence still required.

`literature_v1/` preserves the earlier frozen published-case registry byte for
byte, including exposure chronology. Its companion input and pending report
retain two candidates: one value unavailable and one correction/uncertainty
convention unverified. Neither is currently scored as comparable.

`uncertainty_summary.json` and the byte-preserved full-precision
`envelope_calibration_inputs.json` document reused calibration evidence and
start sensitivity. The 0.6 and 13 MHz half-widths have no nominal probability.
They are based on the balanced architecture and are descriptive comparisons
for the other three roles. No full model covariance is available.

Use `python -m i2spec.validation verify REGISTRY` and
`python -m i2spec.validation score REGISTRY MEASUREMENTS --out NEW_REPORT`.
Scoring never changes model coefficients, forecasts, widths or primary model.
The full workflow and limitations are in `docs/research-model-validation.md`.
