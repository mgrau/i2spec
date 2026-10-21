# Frozen research release: October 21, 2026

One model, `research-20261021-i2spec` (role `i2spec`): the model the apps and the documentation present. It is
`research-20261020-balanced-expanded` for every isotopologue and every J the line lists reach. Select it as
`i2spec.ResearchModel("i2spec", isotopologue)` or through `i2spec.preferred.model(isotopologue)`.

The release holds only what it adds. `model.json`, `hyperfine.json` and `uncertainty.json` are read from
[`../20261020/balanced_expanded/`](../20261020/README.md), and their SHA256 must match the ones recorded in
`manifest.json` (`extends`). For ¹²⁷I₂ up to J'' = 200 every position, component and uncertainty is the frozen
model's, bit for bit.

`i2spec/isotopologues.json`:

- **Levels of ¹²⁹I₂ and ¹²⁷I¹²⁹I** come from the same potentials at their reduced masses. The rotational function is
  scaled by ρ = μ(¹²⁷I₂)/μ. The B state gets an isotope shift relative to X,
  ΔV_B = (1 − ρ)(1 − y)^p (a0 + a1 y) with y = (R − Re)/(R + Re). a0 and a1 are fitted to the CIPM isotope shifts
  (bipm2003a/d: six ¹²⁹I₂ lines and one ¹²⁷I¹²⁹I line, each tied to a measured ¹²⁷I₂ frequency). The power p = 5 has
  the lowest leave-one-line-out rms among 0, 3, 5 and 8 (0.078 MHz). The pending sets of the same lines stay blind.
- **Hyperfine of the other isotopologues** uses the ¹²⁷I₂ functions per nucleus, scaled by the Q and g ratios of
  Salumbides et al. (2006). In the core, B uses the radial functions over the isotopologue's own wavefunction and X
  the formulae at the ¹²⁷I₂-equivalent level. Beyond the core, the ¹²⁷I₂ constants are interpolated at that
  level. Nothing is fitted: the CIPM intervals give 0.19 MHz rms, against 0.17 MHz for the 2006 formulae fitted
  to them. A floor of 0.10 MHz per component is added to the uncertainty.
- **Uncertainty of the other isotopologues.** The fit's covariance and the shift's two coefficients are
  conditioned on the seven CIPM lines (`joint_factor`), because the ¹²⁷I₂ levels do not separate a potential's
  shape from its rotational function, which the isotopologues see scaled by ρ. Two further terms are added: the
  spread of the p = 3 and 8 forms, and the leave-one-line-out floor. A B level more than 2 in ¹²⁷I₂-equivalent v
  from the lines that fix the shift (v′ = 6–14) is *extrapolated*.
- **Above J = 200** (`j_split`), levels are solved on the long grid with the barrier rule, up to J'' = 445.
  Hyperfine constants beyond the tables' J = 203 are held at that J. A level's model σ grows as
  J(J + 1)/J_r(J_r + 1) beyond the J the training rows reach at its v. This is a rule, not a fit, and it also
  applies to ¹²⁷I₂ levels beyond their data below J = 200.

`i2spec/reference_predictions.json` holds centers, components and one intraline interval per line for 10 ¹²⁷I₂,
10 ¹²⁹I₂ and 5 ¹²⁷I¹²⁹I lines. `ResearchModel("i2spec", iso).verify()` reproduces them exactly.

Built by `prototypes/explore/physical_search/package_i2spec.py`
(stages `fit write condition references verify`), with `isotopologues.py` (in the gitignored `prototypes/explore/`).
