# Frozen release: October 20, 2026 — the continued model

The core fit of the [October 18 release](../20261018/README.md) continued to X v'' ≤ 108 and B v' ≤ 80, with its
uncertainty model (artifact identifier `research-20261020-balanced-expanded`). The released model
(`../20261021/`) extends it to the other isotopologues and to J'' = 445; select that one, as
`i2spec.preferred.model(isotopologue)`.

`balanced_expanded/model.json` holds the new terms and settings only. The base coefficients are
read at load time from `../20261018/balanced/model.json`, and their SHA256 must match the one
recorded here. `hyperfine.json` holds the out-of-core hyperfine constants per level and the
B v'77/78 resonance data. `reference_predictions.json` holds 842 line centers, 2,586 component
frequencies and 1,434 intraline intervals from the prototype code, which `verify()` reproduces
within 2 Hz.

`audits.json` covers:

- reproduction of the final fit's training rows (within 1.4 kHz);
- the port's agreement with the prototype potentials;
- the Jacobian, physical and grid audits, and the restarts;
- the hyperfine tables' coverage;
- how much the source of the rotational energies moves the components.

`validation.json` holds the four-fold cross-validation of the final settings, the nested
hyperfine cross-validation, the held-out comb benchmark, the comparison with an earlier version of the model
and with Hannover 2008 on the same rows, and the isotopologue check. `selection.json` records the
decisions. `provenance.json` lists the SHA256 of every script, result and data file used.

Inventory:

- 105 fitted molecular coefficients: 61 base, 38 new potential terms and 6 new hyperfine
  parameters;
- 73 inherited: 18 base and 55 theory/atomic inputs to the hyperfine;
- 5 calibration terms;
- no line or level terms.

The exact P(13) 43-0 exclusion is retained.

`uncertainty.json` (added October 21, 2026) is the frozen uncertainty model behind
`ResearchModel.uncertainty(...)`. It gives the 1σ of every line center, hyperfine component,
intraline interval and interval between two lines in the domain. It has three parts:

- **Statistical.** The fit's linearized covariance is propagated through each level's Hellmann–Feynman
  derivatives. It is stored as a 97 × 86 factor over the 38 free new terms (B's dDe is X's), the 47 frozen core potential
  coefficients, 7 inherited dispersion factors and 5 calibrations.
- **Model error.** Smooth functions of state, v and J, fitted to held-out residuals of the core and
  extension cross-validation.
- **Hyperfine.** A floor plus a pattern-scale error, fitted to held-out intervals and held-out
  measured constants.

It also holds each level's data status (measured, interpolated or extrapolated), the recipe and its checks,
and the SHA256 of its inputs. No pending-set outcome entered it. `uncertainty_references.json` holds 69 values
that the runtime reproduces exactly (`tests/test_research_expanded.py`).

`validation/pending_blind_v1/` is the frozen forecast registry for the pending data sets. It was
created after this release and before any of their values were read. Its `coverage_posthoc.json` scores
the uncertainty model on those rows. The scores are post hoc: the accuracy results had been read before
the uncertainty model was built, although its recipe used only the training fit and its cross-validation.
