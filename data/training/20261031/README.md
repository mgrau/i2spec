# Training inputs of release 20261031

Release 20261031 is rebuilt from the repository by `i2spec.fit.final.rebuild()` (stage by stage, each on the previous
one's output): its fits read the training rows of release 20261021 (`../20261021/`), the six CIPM isotopologue line
centers and the cross-validation of release 20261029 (`../20261029/`), the roles of the pending sets
(`../roles_20261029/`) and the scored third blind registry (`src/i2spec/data/research/20261029/validation/pending_blind_v3`)
for the calibrate residuals. This folder holds the three inputs that are not there.

| file | what | read by |
|---|---|---|
| `hyperfine_folds.json.gz` | the nested cross-validation of the release's hyperfine form beyond the core (`hyperfine_beyond.FORM_20261031`): at each 1g(1Pig) C5 scale 0.6-0.9, each of the ten B v' blocks held out in turn, the fit's background, angles and curvatures, training costs and held-out eqQ (70 refits; `fit.final.hyperfine_fold` refits one) | `i2spec.fit.final.cross_validation` |
| `s2_terms.json.gz` | second-order terms the release's hyperfine stage read from its cache that no computation on a committed B potential reproduces: 390 levels of the 1g(1Pig) state at the fitted C5 scale 0.85 (every level the fit reads among them) and B v' 43 J' 12 of six states | `i2spec.fit.final.s2_fixed` |
| `edges.json` | the rules for the levels beyond the tables (B v' 81-86, X v'' 109-113), the J rule and X C slope status, and the tunneling width of quasi-bound B levels, each with the evidence it was set on; rules, not fits | `i2spec.fit.final.edges` |

None holds a pending value: the folds are fitted to the training hyperfine rows (`../20261021/hyperfine.csv.gz`), the
second-order terms are computed from the potentials alone, and the edge rules' evidence quotes scored registries only.

The release read its second-order terms from a cache keyed without the B potential, so each came from the first
computation that filled it. The other terms are computed again by the rule `fit.final.s2_frozen`: every state on
release 20261021's hyperfine potential (`../20261021/hyperfine_potential.json`) except the 1g(1Pig) at the fitted
scale, on the release's continuation. The cache key now holds the potential.
