# Per-line uncertainties (M6)

Two halves. The **discrepancy** half — what the model actually misses, region by region, measured
against every observation we hold — is done and is what `lookup.uncertainty` quotes. The **parameter**
half — the covariance of a fitted model propagated to each line — waits for the joint MLR fit
(`docs/design/mlr-x.md`), whose Jacobian is in place but whose fit is not converged.

## The error map

`prototypes/error_map.py` computes obs − model with the default model (`i2spec2026a` potentials, the
hyperfine table) for the 1013 observations that constrain a line position, and pools them by
isotopologue and by v′ and v″ band. MHz; "floor" is the rms of the stated measurement uncertainties.

| isotopologue | v′ | v″ | rows | lines | rms | median | 90 % | floor | rule before |
|---|---|---|---|---|---|---|---|---|---|
| ¹²⁷I₂ | 0–17 | 0–5 | 672 | 618 | **1.17** | 0.66 | 1.87 | 1.09 | 3.0 |
| ¹²⁷I₂ | 0–17 | 6–12 | 16 | 8 | 2.40 | 1.80 | 2.55 | 0.35 | 3.0 |
| ¹²⁷I₂ | 0–17 | 13–17 | 55 | 29 | 3.21 | 2.46 | 5.50 | 0.15 | 3.0 |
| ¹²⁷I₂ | 18–30 | 0–5 | 70 | 56 | 2.30 | 1.31 | 4.13 | 1.79 | 3.0 |
| ¹²⁷I₂ | 31–44 | 0–5 | 32 | 32 | 17.1 | 4.84 | 38.1 | 0.88 | 3.0 / 150 |
| ¹²⁷I₂ | 45–53 | 0–5 | 7 | 7 | 23 000 | 3 960 | 38 900 | 0.01 | 60 000 |
| ¹²⁷I₂ | 54–70 | 0–5 | 1 | 1 | 41 000 | — | — | 1.0 | 60 000 |
| ¹²⁷I₂ | 31–44 | 48–60 | 36 | 10 | 419 000 | 365 000 | 631 000 | 0.07 | 650 000 |
| ¹²⁹I₂ | 0–17 | 0–5 | 54 | 5 | 3.36 | 2.02 | 4.98 | 1.51 | 15.3 |
| ¹²⁹I₂ | 0–17 | 6–12 | 29 | 1 | 1.65 | 1.69 | 1.88 | 1.06 | 15.3 |
| ¹²⁷I¹²⁹I | 0–17 | 0–5 | 41 | 1 | 0.25 | 0.10 | 0.38 | 1.01 | 15.3 |

Three things it says:

- **Where the data are, the model is at the data.** 672 rows at v′ ≤ 17, v″ ≤ 5 sit at 1.17 MHz rms
  against a 1.09 MHz floor. The rule's 3 MHz there is conservative, and stays: the cell mixes
  0.4 MHz sets with 2 MHz ones, and the model's error is not uniform inside it (v″ = 13–17 is 3.2).
- **The rule was wrong in two places.** Isotopologues carried a 15 MHz allowance from the era before
  the V_ad constant (`docs/design/parameter-sets.md`); ¹²⁹I₂ is now 3.4 MHz rms and ¹²⁷I¹²⁹I 0.3, so
  the allowance is 3 MHz. And v′ = 31–43 in the BIPM range claimed 3 MHz where the B levels drift by
  ±6 MHz across v′ = 32–37 (the BIPM 532 nm intervals) and sit 5.5 MHz off at v′ = 39 (Nishiyama 2024):
  it is now 5 MHz.
- **Above v′ = 44 the default model is not usable for positions**, 4–40 GHz, exactly as the rule says.
  The MLR B potential brings that to ≲150 MHz (`docs/design/mlr-x.md`), but it is not the default until
  the joint fit holds the low-v′ precision at the same time.

**Coverage.** Over the 1013 rows, 95 % of residuals fall within one rule σ and 100 % within two
(a calibrated Gaussian would give 68 % and 95 %). The rule is a bound, not a standard deviation, and
reads as one.

## What is still missing

- **Parameter uncertainty.** A line the data do not touch (v″ = 18–47, J″ above the measured range,
  high v′) has an error the map cannot see. That is the joint fit's covariance, through
  ∂ν/∂θ = ⟨ψ|∂V/∂θ|ψ⟩ (`prototypes/mlr_joint_fit.py` computes exactly these), plus the priors on the
  long-range constants. Not before the fit converges.
- **Coverage by held-out set.** The honest test of any uncertainty is to hold out one data set, fit,
  and check that its residuals divided by the predicted σ have unit variance. Same dependency.
- **Hyperfine components.** The map is for line centres. With the table, components at a measured v′
  are good to 1–10 kHz; elsewhere the formulae's 25 kHz–1 MHz applies (`docs/design/hyperfine-fit.md`).
