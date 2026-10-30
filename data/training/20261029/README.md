# Training inputs of release 20261029

Release 20261029 is fitted to the training rows of release 20261021 (`../20261021/`, unchanged). This folder holds
the two inputs of its fits that are not there.

| file | what | read by |
|---|---|---|
| `isotopologue_centers.json` | the six CIPM isotopologue line centers (bipm2003a/d: five ¹²⁹I₂ lines and one ¹²⁷I¹²⁹I line), components reduced with the model's hyperfine, each with its σ; core training rows (floor 0.1 MHz) | `i2spec.fit.isotopologues.core_rows` |
| `cross_validation.json.gz` | outputs of the cross-validation refits the uncertainty model is built from, in the layout of `../20261021/cross_validation.json.gz`: the continuation's four v folds, the core fit's seven physical folds and its future fold, the nested hyperfine folds, and the benchmark predictions with the rows' own σ | `i2spec.fit.uncertainty.build(cv=...)` |

The fits' settings: comb floor 0.3 MHz, the B limit tied to D0(X) plus the atomic splitting, the B isotope shift of
power 5 in the core fit, B eqQ order 6 in the core hyperfine, B rotation order 2 in the continuation. The exact line
127I2 P(13) 43-0 is left out of every fit and calibration.
