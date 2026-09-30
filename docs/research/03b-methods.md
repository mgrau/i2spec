# Stage 3b — Global fitting and uncertainty quantification methods (plus implementation-stack check)

Status: complete draft, 2026-09-15.

Verification limits:
- Paywalled full texts (ScienceDirect) and Le Roy's program site could not be read.
- Claims about Watson (2003), Le Roy (1998) and MARVEL come from metadata and from passages quoted in
  open-access citing papers. Unverified details are flagged in place.
- No local wasm32 compile test was possible; the wasm findings come from the crates' repositories and CI.
Scope: methods for a global fit of I₂ B–X (and X) potentials plus hyperfine Hamiltonian to heterogeneous
data (FTS line positions, sub-Doppler hyperfine components, relative hyperfine splittings, fluorescence
term values), and for uncertainty quantification (UQ) of predicted line positions across the whole
interpolation range. Radial-solver benchmarking is out of scope (handled elsewhere).

Conventions used in this document:
- **[Source]** marks statements attributed to a cited reference (verified to exist; content from abstracts,
  publisher/repository pages or documentation unless noted).
- **[Analysis]** marks our own reasoning, proposals or estimates.
- "unknown" means we could not verify the point.

---

## 1. Global fitting strategies in diatomic spectroscopy

### 1.1 DPF vs Dunham/term-value fits vs MARVEL-style networks

**What the sources establish**
- [Source] Knöckel, Bodermann & Tiemann (2004) fit X-representation potentials for X and B directly to I₂ B–X
  data, i.e. a direct potential fit (DPF). The details in the brief (32 B + 14 X + 8 BO parameters, Numerov,
  MINUIT, 3 MHz floor, ~1500 lines) are taken from the project brief and were not re-read here. The hyperfine
  structure was handled with the separately fitted interpolation formulae of Bodermann, Knöckel & Tiemann
  (2002).
- [Source] Le Roy's program suite (JQSRT vol. 186, 2017) separates the three classical strategies into
  separate tools:
  - *dParFit*, "fitting … to parameterized level energy expressions" (Dunham/band-constant type);
  - *DPotFit*, "to potential energy functions" (DPF);
  - *RKR1* (semiclassical inversion) and *LEVEL* (radial Schrödinger solver).

  Published DPFs have used the same approach for N₂, CO, Rb₂, Li₂, the hydrogen halides and others (Le Roy et al.
  2006, 2009; Coxon & Hajigeorgiou 2004; Seto et al. 2000, whose Rb₂ paper is subtitled "Nothing else will do!").
- [Source] MARVEL (Furtenbacher, Császár & Tennyson 2007; algorithmic update 2012) inverts a network of
  assigned transitions into empirical energy levels with uncertainties. Tóbiás et al. (2020) state that
  network theory "based on the generalized Ritz principle" can be used to design and validate kHz-level
  precision-spectroscopy experiments and to derive accurate energy differences. ExoMol's 2024 release
  describes extensive "MARVELization" of line lists for high-resolution use (Tennyson et al. 2024).

**Assessment for I₂ B–X** [Analysis]

| Strategy | Parameters | Predicts unobserved lines? | Role for us |
|---|---|---|---|
| Dunham / band constants | Many Yₖₗ, or per-band constants | Only by polynomial extrapolation, which is poor near the B-state dissociation limit and at J ≈ 300 | Diagnostic only |
| Term-value fit (every level free) or MARVEL network | One per level: ~10⁴ X+B levels, a sparse *linear* LSQ | No (observed levels only) | Model-free validation layer: finds misassignments, estimates per-dataset calibration offsets independently of the potential model, gives empirical level uncertainties |
| DPF (V_X(R), V_B(R), BOB/adiabatic corrections, hyperfine radial functions) | ~50–150 | Yes, physically constrained, including isotopologues through mass scaling | **Primary model** |

The recommended architecture is hybrid:
1. Run a MARVEL-like linear network/term-value fit first, to detect outliers and misassignments and to
   estimate calibration nuisance parameters without potential-model bias.
2. Then run the DPF on the cleaned data. The network step is cheap: sparse linear least squares with ~10⁵
   rows.
3. The network also gives an independent check of DPF residual structure. Systematic DPF-minus-network
   differences in a (v, J) region show that the model, not the data, is failing there.

### 1.2 Combining data with uncertainties spanning ~6 orders of magnitude

**Sources**
- [Source] Knöckel et al. 2004 used a fixed 3 MHz floor on all uncertainties (project brief).
- [Source] When data are mutually inconsistent, metrology practice uses the Birge ratio, √(χ²/ν), to diagnose
  under-stated uncertainties (Birge 1932). CODATA adjustments apply expansion factors to input uncertainties
  (Tiesinga et al. 2021).

**Analysis** [Analysis]
- A plain 1/σ² weighting is optimal only if (i) the stated σᵢ are correct and (ii) the model is adequate at
  the σᵢ level. Point (ii) is the real problem.
  - kHz hyperfine data require the model to be right at ~3×10⁻⁸ cm⁻¹. A ~50-parameter potential plus an
    effective hyperfine Hamiltonian may not be that good everywhere.
  - If the kHz points are taken at face value, they dominate χ², and small local model deficiencies are
    "fixed" by distorting potential parameters that the FTS data cover globally.
  - Knöckel's 3 MHz floor is a crude but effective guard against this.
- We recommend replacing the fixed floor with estimated per-dataset (or per-data-class) error-model
  hyperparameters:

  σ²ᵢ,eff = (s_k σᵢ)² + τ_k² + σ²_disc(vᵢ, Jᵢ)

  - s_k is a scale factor, a per-dataset Birge ratio.
  - τ_k is a floor or "jitter" term.
  - σ_disc is a model-discrepancy term (§2.4).

  They are estimated by maximum marginal likelihood or sampled in the Bayesian fit, then frozen for
  production. This makes the floor data-driven and reportable. It also stops a single optimistic dataset from
  dominating.
- Numerics: σ ratios of ~10⁶ mean weight ratios of ~10¹².
  - Never form normal equations (condition number squared). Solve the weighted Jacobian problem with
    QR/SVD, with column scaling.
  - Compute residuals as differences of eigenvalues with absolute accuracy ≪ 1 kHz. This is handled by the
    solver work stream.
- Most relative hyperfine splittings (kHz) constrain hyperfine parameters almost exclusively. The absolute
  sub-Doppler frequencies (kHz–MHz) are what tie the rovibronic term values to the kHz scale. This block
  structure (§1.6) limits how much damage the high weights can do to the potential. It is still worth
  enforcing through the discrepancy term.

### 1.3 Weighting and robust (outlier-resistant) fitting

**Sources** [Source]
- Watson (2003), "Robust weighting in least-squares fits" (J. Mol. Spectrosc. 219, 326), is the standard
  reference in molecular spectroscopy. MARVEL uses it: "The uncertainties of the transition frequencies weight
  this inversion process using a robust reweighting procedure advocated by Watson (2003) allowing Marvel to
  yield the uncertainty of each extracted energy level" (McKemmish et al. 2017, ApJS 228, 15, quoting the
  MARVEL methodology of Furtenbacher & Császár 2012).
- The scheme is iterative. Each point's weight depends on its residual from the previous iteration, which
  lowers the weight of outliers without hard rejection.
- *Unverified detail:* we recall Watson's weight as wᵢ = 1/(σᵢ² + δᵢ²/3), where δᵢ is the previous-iteration
  residual. The full text is paywalled and was not re-read in this session. **Check before implementing.**

**Assessment** [Analysis]
- With 10⁴–10⁵ legacy lines, some misassignments, unresolved blends and transcription errors are certain.
  Candidate approaches:
  - Watson-type IRLS: a smooth reweighting that approaches a Cauchy-like influence for large residuals.
  - A Student-t likelihood per data class, with ν estimated. This is its Bayesian equivalent.
  - Hard sigma-clipping. Brittle, and its outcome depends on the order of operations.
- **Caveat:** robust weights computed from an *inadequate* model down-weight the very data that reveal the
  inadequacy, e.g. a high-J or near-dissociation region, or a local perturbation.
  - So down-weighting must be audited. Every line whose effective weight drops by more than ×10 is logged
    with (dataset, v′, J′, v″, J″).
  - If down-weighted lines cluster in (v, J) space, that is treated as a model problem (§1.8, §2.4), not
    noise.
  - The model-free network check (§1.1) separates the two cases. An outlier in the network fit is a data
    problem. A line consistent in the network but deviant in the DPF points to the model.
- Uncertainty propagation must use the *final* weights. Where misspecification is suspected, the
  sandwich (Huber–White) covariance (JᵀWJ)⁻¹ JᵀW Ŝ W J (JᵀWJ)⁻¹, with Ŝ the diagonal of squared residuals, is a
  cheap robustness check against the model-based (JᵀWJ)⁻¹.

### 1.4 Sequential rounding and refitting (Le Roy)

**Sources** [Source]
- Le Roy (1998), "Uncertainty, Sensitivity, Convergence, and Rounding in Performing and Reporting
  Least-Squares Fits", defines a *parameter sensitivity*. Later papers restate it:
  - it "essentially measures how fast the observables change when one parameter is varied with all others
    held fixed", and "to reproduce the results from a set of parameters, it is often necessary to specify
    many more digits than implied by the uncertainty" (Berninger et al. 2013);
  - "the parameters need to be specified to within their sensitivities, not their confidence limits, in
    order to reproduce the results" (Takekoshi et al. 2012).
- The *sequential rounding and refitting* procedure is:
  1. round the parameter with the largest relative uncertainty at the level of its sensitivity;
  2. fix it, and refit the remaining parameters;
  3. repeat.

  It yields a compact table that still reproduces the data. The algorithm description follows Le Roy's
  programs and papers as we know them; its exact steps were not re-read in the 1998 full text.

**Assessment** [Analysis]
- For a software package the rounding problem mostly disappears. We ship full-precision (float64)
  parameters together with the covariance or its factor (§2.7), in a machine-readable file.
- Le Roy's procedure remains useful for (i) a printed table in the paper and (ii) as a test. The rounded set
  must reproduce all predicted lines to < 0.1σ_pred everywhere on the prediction grid, not only at data
  points.
- An orthogonalized basis (§1.8) greatly reduces the number of digits needed. Parameter sensitivities
  become comparable to parameter uncertainties when the correlation is removed.

### 1.5 Per-dataset calibration offsets and scale factors as nuisance parameters

**Source** [Source] Gerstenkorn & Luc (1979), abstract: "the wavenumbers previously published in the iodine-atlas
extending from 14 800 to 20 000 cm⁻¹ must be corrected by subtracting 0.0056 cm⁻¹ from all wavenumbers. The
accuracy of the absolute wavenumbers obtained in this way is estimated ± 0.002 cm⁻¹." Salami & Ross (2005)
distribute the atlas in ASCII form. Which correction their file includes is not verified here.

**Proposed treatment** [Analysis]
- **Observation model per dataset k:** ν_obs = (1 + κ_k) ν_model + δ_k, optionally piecewise by atlas
  segment or plate.
  - A pure FTS scale error grows with ν. Over 14 800–20 000 cm⁻¹, a single additive correction and a scale
    factor are not equivalent, so the data should decide.
  - κ_k and δ_k enter linearly. They can be profiled out exactly at each iteration (variable projection,
    Golub & Pereyra 1973) or appended to the Jacobian. The cost is negligible.
- **Anchoring:** comb-referenced absolute frequencies are the anchor (κ = δ = 0, fixed). Without an anchor,
  the δ_k are degenerate with the B–X term offset.
- **Priors / penalties:** use each dataset's stated calibration uncertainty, e.g. N(0, 0.002 cm⁻¹) on δ for
  the atlas after the 1979 correction. A per-dataset prior keeps the nuisance parameters from absorbing
  model error.
- **Other nuisance types:**
  - Fluorescence series: each series' upper-level energy is a free parameter. This is common DPF practice; we
    did not verify DPotFit's specific implementation.
  - Sub-Doppler data sets: pressure/power/modulation shifts become per-set offsets, with priors taken from the
    reported shift coefficients.
  - Isotopologue data: isotope-specific shifts are handled physically through mass scaling and BOB terms,
    not by free offsets.
- **Blend-aware forward model (important):** a Doppler-limited FTS "line position" of ¹²⁷I₂ is the centre of
  a blend of hyperfine components. The manifold width is comparable to the Doppler width (hundreds of MHz),
  while the FTS precision is 30–300 MHz.
  - The forward model should therefore compute the observed FTS position from the hyperfine components,
    weighted by intensity and convolved with the instrument and Doppler profile. It should not equate the
    FTS position with the hyperfine-free transition.
  - Otherwise, a J- and v-dependent bias of order tens of MHz is absorbed into the potential. This is an
    order-of-magnitude estimate that should be checked in simulation.
  - The same applies to blended lines of different (v′, J′) within one FTS feature. Flag them or model them
    as a sum.

### 1.6 Hyperfine and rovibronic parameters: simultaneous vs sequential

**Sources** [Source]
- Bodermann, Knöckel & Tiemann (2002) give "widely usable interpolation formulae for hyperfine splittings in
  the ¹²⁷I₂ spectrum". Knöckel et al. (2004) then fit the rovibronic structure. This is a sequential
  approach.
- Effective-Hamiltonian fitting of hyperfine structure is classical (e.g. Pickett 1991, SPFIT/SPCAT).
- Recent variational treatments compute hyperfine-resolved rovibronic spectra of diatomics directly from
  R-dependent coupling curves:
  - Qu, Yurchenko & Tennyson (2022), JCTC;
  - Yin et al. (2026), JCP, for homonuclear diatomics.

**Assessment** [Analysis]
- Observable structure:
  - an absolute component frequency = hyperfine-free rovibronic transition + E_hfs(upper) − E_hfs(lower);
  - a relative splitting depends only on the hyperfine terms;
  - an FTS position is the blend centroid (§1.5).
- The joint Jacobian is therefore block-structured and nearly block-diagonal. Rovibronic parameters couple to
  hyperfine ones only through vibrational averaging of the hyperfine radial functions and through the
  blend-centroid model.
- The problem with a *sequential* fit is not bias. The problem is that the first stage's covariance is
  usually discarded, so the second stage's uncertainties are too small and the cross-covariances are lost.
- Recommendation: one **joint likelihood**. Hyperfine coefficients are parameterized as smooth functions of R
  (eqQ(R), C_I(R), δ(R), …) and vibrationally averaged with the same wavefunctions, or as a (v, J) surface
  where the R-picture is inadequate, as in the Bodermann-style formulae.
  - Solve it with block elimination (Schur complement): the result equals the simultaneous fit, but the
    cost is close to that of two separate fits.
  - A sequential fit remains useful as the initializer.

### 1.7 Model-order selection

**Sources** [Source]
- Information criteria: AIC (Akaike 1974), BIC (Schwarz 1978).
- Cross-validation for model selection is surveyed by Arlot & Celisse (2010).
- Bayesian evidence (the "Occam factor") under a Laplace/Gaussian approximation is developed by MacKay (1992).

**Assessment** [Analysis]
- With N ≈ 10⁵ and under-stated or non-Gaussian errors, AIC/BIC penalties (2k, k ln N) are tiny compared
  with the Δχ² from fitting noise or model deficiencies. They therefore over-select parameters unless the
  error model (§1.2) is calibrated first, i.e. the dimensionless standard error ≈ 1 per data class. Use them
  only after that calibration.
- The criterion that matters for us is **predictive accuracy where there are no data**. Random K-fold CV
  leaves near-duplicates (neighbouring J) in the training set and is over-optimistic. Use **blocked CV**:
  - leave-one-v′-band-out;
  - leave-J-range-out, e.g. drop J > 250 and predict it;
  - leave-one-dataset-out;
  - leave-one-isotopologue-out.

  Blocked CV errors estimate interpolation and extrapolation error directly and feed the discrepancy model
  (§2.4).
- An **order scan** gives an estimate of model-form uncertainty. Refit with increasing numbers of B- and
  X-potential terms (and different b or R_m in the X-representation), and plot the spread of predictions over
  a (v, J) grid.
- Evidence (Laplace) is a reasonable automatic criterion for order, and for choosing Tikhonov/prior
  strengths, once the likelihood is calibrated.

### 1.8 Strongly correlated potential parameters

**Sources** [Source]
- Transtrum, Machta & Sethna (2010, 2011) explain why multi-parameter nonlinear fits are hard by viewing
  fitting "as a generalized interpolation procedure" on the model manifold. Gutenkunst et al. (2007) show that
  "sloppy" parameter spectra are generic.
- Pashov, Jastrzębski & Kowalczyk (2000), the IPA method, construct point-wise diatomic potentials. The
  claim that IPA uses SVD with singular-value truncation is our recollection; the full text was not re-read.
- Hansen (1998) is the standard reference for truncated SVD and Tikhonov regularization of rank-deficient
  problems.
- Le Roy & Henderson (2007) introduced the MLR form with built-in long-range behaviour.
- Tellinghuisen (1990, 2001) treats error propagation from correlated fit parameters, including to RKR
  potentials.

**Assessment** [Analysis]
- A 32-term power series in X is Vandermonde-like, and its parameters are strongly correlated. This is
  harmless for *predictions* inside the data-supported region, provided:
  1. the full covariance is kept, and parameters are stored with enough digits (see Le Roy's "sensitivity",
     §1.4);
  2. predictions are propagated with the full covariance (§2.1).

  Reporting individual parameter uncertainties is not meaningful; reporting prediction uncertainties is.
- Conditioning fixes that do not change the function space:
  - Replace Xⁱ by orthogonal polynomials, e.g. Chebyshev or Legendre on X ∈ [−1/b, 1), since X → 1 as
    R → ∞ and X = −1/b at R = 0.
  - Scale columns.
  - Use SVD-based Levenberg–Marquardt / trust-region steps (with geodesic acceleration, Transtrum 2011, if
    convergence is slow).

  These reduce the parameter condition number by orders of magnitude and make rounding/refitting less
  fragile. Predictions are unchanged.
- Regularization, which does change the solution: truncated SVD or Tikhonov with a smoothness penalty on
  V″(R) outside the classically sampled region. In the Bayesian view this is a Gaussian prior (§2.3). Its
  strength is chosen by evidence or blocked CV (§1.7), never by eye.
- Near the B-state dissociation limit, a functional form with the correct long-range tail (MLR-like, with
  theory-fixed C₅/C₆/C₈ where available) matters more for extrapolation than the number of polynomial
  terms. The high-v′ extrapolation region is where the model form, not the parameter covariance, dominates
  the uncertainty (§2.6).

---

## 2. Uncertainty quantification across the interpolation range

### 2.1 Linearized covariance propagation (Jacobian + Hellmann–Feynman)

**Physics** [Source + Analysis] The Hellmann–Feynman theorem (Feynman 1939) gives, for a normalized
eigenstate, ∂E_n/∂p = ⟨ψ_n|∂H/∂p|ψ_n⟩.
- For a potential V(R; p), ∂E_{vJ}/∂p_k = ⟨vJ|∂V/∂p_k|vJ⟩. This is one quadrature per level and parameter,
  using wavefunctions the solver already has.
- The same holds for BOB/adiabatic-correction and centrifugal-correction parameters, through their radial
  functions.
- For the hyperfine matrix H_hfs (finite-dimensional), ∂E/∂c = uᵀ(∂H_hfs/∂c)u. This is first-order,
  non-degenerate perturbation theory. Exactly degenerate hyperfine levels need the degenerate form.
- For a transition, ∂ν/∂p = ∂E′/∂p − ∂E″/∂p.

This gives exact analytic Jacobians at negligible extra cost. No finite differences are needed, and it
avoids differentiating through the eigensolver.

**Propagation** [Analysis]
- Σ_p = (JᵀWJ + Λ_prior)⁻¹, computed by QR/SVD, with W from the calibrated error model (§1.2).
- For any set of predicted lines with sensitivity matrix G (rows gᵢ = ∂νᵢ/∂p):
  - Cov(ν_pred) = G Σ_p Gᵀ;
  - Var(νᵢ) = gᵢᵀ Σ_p gᵢ.
- Correlations between predicted lines are large. For example, all lines of one band share v′-dependent
  errors. They matter for any downstream use that combines lines, such as a frequency reference chain or a
  fit of an unknown spectrum to the model.
- [Source] Error propagation with correlated parameters is treated by Tellinghuisen (1990) for RKR
  potentials and by Tellinghuisen (2001) in general.
- **Validity:**
  - The linearization holds if the forward model is close to linear over the posterior width of the
    predictions. For sloppy models this usually holds for predictions even when parameter marginals look
    strange (Transtrum et al. 2010, 2011).
  - It must still be *tested*: against Laplace-vs-MCMC (§2.3) and against parametric bootstrap coverage
    (§2.2).
  - Linear propagation reflects only parametric uncertainty within the chosen model form. It says nothing
    about model-form error, hence §2.4.

### 2.2 Bootstrap and jackknife

**Source** [Source] Efron (1979) introduced the bootstrap, as "another look at the jackknife".

**Assessment** [Analysis]

| Variant | What it captures | Cost | Use |
|---|---|---|---|
| Case (line) resampling | Random scatter only. With 10⁵ lines it ≈ linear covariance | B full refits | Low value |
| **Cluster/block bootstrap by dataset or band** | Dataset-level systematics (calibration, line-shape conventions) that per-line σ's miss | B refits, warm-started | Worthwhile offline, B ≈ 100–200 |
| Leave-one-dataset-out jackknife | Influence of each dataset; same as blocked CV (§1.7) | K refits | **Recommended** diagnostic |
| Parametric bootstrap (simulate from fitted model + error model, refit) | Estimator bias; coverage of the linear/Laplace intervals | B refits | **Recommended** validation of §2.1 |

Each refit is a warm-started nonlinear LSQ and should converge in a few Gauss–Newton steps. Feasibility
depends on the forward-model cost, which is owned by the solver work stream. The bootstrap is a
*validation* tool here, not the production UQ.

### 2.3 Bayesian inference: MCMC/HMC and the Laplace approximation

**Tools** [Source]
- emcee is an affine-invariant ensemble MCMC sampler (Foreman-Mackey et al. 2013).
- NUTS is an adaptive-path-length Hamiltonian Monte Carlo (Hoffman & Gelman, arXiv:1111.4246).
- NumPyro is a JAX-based probabilistic programming library with NUTS (Phan et al., arXiv:1912.11554).
- BlackJAX provides composable samplers in JAX (Cabezas et al., arXiv:2402.10797).
- MacKay (1992) develops the Laplace/Gaussian approximation to the posterior and the evidence, used for model
  comparison and for setting regularization constants.

**Assessment** [Analysis]
- **Dimension.** There are ~100–150 potential/hyperfine parameters, plus 2–3 nuisances per dataset, plus
  error-model hyperparameters. Every likelihood evaluation is a full forward model: all levels of X and B for
  every J, plus hyperfine diagonalizations.
- **Samplers:**
  - *emcee:* ensemble samplers need many walkers (≳2× the dimension). Affine invariance handles linear
    correlations but is inefficient in ~150 dimensions with curved, sloppy posteriors. Not recommended as the
    main tool.
  - *NUTS (NumPyro/BlackJAX):* needs gradients, which Hellmann–Feynman supplies analytically (§2.1). It can
    be wrapped as a custom JAX primitive with a custom JVP/VJP, so the eigensolver itself need not be
    differentiable.
  - With the Laplace covariance as a dense mass matrix (whitening), NUTS is the right tool for *validating*
    Gaussianity of predictions and for marginalizing hyperparameters.
- **Laplace approximation:**
  - Take a Gaussian at the MAP with the Gauss–Newton Hessian JᵀWJ + Λ_prior. It costs the same as §2.1 and
    also gives the log-evidence used for model order and prior strength (§1.7).
  - If §2.3 sampling shows predicted-line marginals are Gaussian to within ~10% in width, Laplace is the
    production UQ.
- **Hyperparameters** (s_k, τ_k, discrepancy-GP parameters) are best fitted by maximizing the Laplace
  evidence (type-II maximum likelihood), then held fixed. Alternatively they can be sampled jointly in the
  NUTS validation run.

### 2.4 Model discrepancy (Kennedy–O'Hagan, Gaussian-process residual models)

**Sources** [Source]
- Kennedy & O'Hagan (2001) is the reference framework for calibrating a computer model with an explicit
  model-inadequacy (discrepancy) function, modelled as a Gaussian process (GP).
- Brynjarsdóttir & O'Hagan (2014), abstract: "In order to make appropriate use of observations of the physical
  system it is important to recognize model discrepancy, the difference between reality and the simulator
  output."
- Nuclear physics has the closest analogue to our extrapolation problem:
  - Neufcourt et al. (2018) apply "Bayesian Gaussian processes and neural networks" to extrapolate
    two-neutron separation energies, and find that GPs "deliver a more stable performance".
  - Melendez et al. (2019) propose a GP-based truncation-error model "for predictions that are correlated
    across the independent variables" (see also Furnstahl et al. 2015).
- GP extrapolation of quantum observables: Vargas-Hernández et al. (2018).
- Standard GP reference: Rasmussen & Williams (2005).
- We found no verified example of a Kennedy–O'Hagan-type discrepancy model in a diatomic DPF. This appears to
  be new for the field (our search was limited).

**Proposed discrepancy model for I₂** [Analysis]
- **Put the discrepancy on level energies, not on lines:** δ_X(v″, J″) and δ_B(v′, J′), plus separate terms
  per isotopologue if needed.
  - A line inherits δ_B − δ_X. Lines sharing a level then get correctly correlated, with the
    combination-difference structure, and the number of GP inputs is the number of levels (~10⁴), not lines.
  - Inputs in physically scaled coordinates, e.g. (v + ½) or G(v)/D_e, and J(J+1). Anisotropic Matérn
    kernel.
  - Hyperparameters: amplitude and length scales, estimated from blocked-CV residuals or by
    evidence maximization.
- **Implementation as a linear basis.** Approximate the GP by a smooth 2-D tensor B-spline (or
  Hilbert-space/Fourier) basis in (v, J), with a Gaussian prior on its coefficients.
  - The discrepancy then enters the Jacobian as extra *linear* parameters with a prior. The same Laplace/QR
    machinery (§2.1, §2.3) yields posterior predictive covariances that include discrepancy, with no separate
    GP code.
  - The rank equals the number of basis functions, which helps storage (§2.7).
- **Identifiability.** A flexible δ can absorb what the potential should explain. Mitigations:
  - (i) keep the prior amplitude small, i.e. at the level of the blocked-CV residual scatter;
  - (ii) use long length scales;
  - (iii) optionally project the discrepancy basis onto the orthogonal complement of the potential Jacobian's
    column space, so that δ only models what the potential *cannot*.
- **Effect.**
  - In data-dense regions the posterior δ is pinned by the residuals, so the added uncertainty is small.
  - In gaps it reverts towards the prior amplitude, giving an honest, larger uncertainty.
  - Beyond the data (high v′ near dissociation, J > J_max) a stationary GP simply reverts to the prior
    amplitude. This probably *under*-states real extrapolation error, so let the amplitude grow with distance
    from data support (§2.6).
- This formalizes Knöckel's per-region "rule of thumb" into a computed, per-line quantity.

### 2.5 Prior art on uncertainties of predicted line positions (DPF, ExoMol, HITRAN, MARVEL)

**What the sources show** [Source]

| Source | Uncertainty representation | Correlations between lines |
|---|---|---|
| **HITRAN** (documentation page; the table there is attributed to the HITRAN2024 paper; see also HITRAN2020) | Per-line *uncertainty code*, in decade bins in cm⁻¹. Code 0 = "≥1 or unreported", 1 = 0.1–1, 2 = 0.01–0.1, 3 = 0.001–0.01, 4–9 progressively finer down to < 10⁻⁸ cm⁻¹ | None mentioned on the documentation page |
| **ExoMol** | *Per state*, not per line. The dataset definition (`.def.json`) has `uncertainties_available`; when true, the States file carries a field described as "Energy uncertainty in cm-1" (true for ⁴⁸Ti¹⁶O "Toto", false for ¹²C¹⁶O "Li2015", checked in the live def.json files). The 2020 release emphasizes "accurate characterisation of transition frequencies for use in high resolution studies"; the 2024 release describes "extensive MARVELization" | None |
| **MARVEL** | Per-level uncertainties from the weighted (Watson-reweighted) network inversion (McKemmish et al. 2017) | Not reported, as far as we verified |
| **Knöckel et al. 2004** (I₂) | No parameter uncertainties; a rule of thumb per spectral region (project brief) | None |
| **DPF (Le Roy)** | Parameter uncertainties and sensitivities (Le Roy 1998; citing works) | Whether DPotFit reports predicted-level uncertainties: *unknown* (manual not reachable during this session) |

**Assessment** [Analysis]
- No standard line database carries inter-line correlations.
- Per-state uncertainties combined in quadrature assume independent upper and lower levels. For I₂, where
  both levels come from the same fitted potentials, that is wrong in both directions: systematic errors
  cancel within a band and add across bands.
- Our package can do better at negligible cost: exact parametric covariance (§2.1) plus a discrepancy
  covariance (§2.4), available on demand.
- For interchange with HITRAN/ExoMol, export per-line σ (and a HITRAN-style code) and per-state σ as
  lossy summaries.

### 2.6 Quantifying and reporting extrapolation risk

[Analysis]
1. **Support metric per predicted level.** Candidates:
   - the distance, in scaled (v, J) coordinates, to the nearest observed level of the same state and
     isotopologue;
   - the statistical leverage hᵢ = gᵢᵀ Σ_p gᵢ / σ²_local, where σ_local is the uncertainty of nearby data.

   From these, label each prediction *interpolated*, *near-extrapolation* or *far-extrapolation*. The flag
   ships with every line.
2. **Total predicted uncertainty.** σ²_tot = σ²_param (Laplace, §2.1/2.3) + σ²_disc (GP/basis, §2.4) +
   σ²_form, where σ²_form is the spread across alternative model forms: the order scan (§1.7), the X-representation
   vs an MLR-type long-range form, and different BOB parameterizations. The spread gives an honest term where
   the forms disagree, i.e. almost only in extrapolation.
3. **Calibration of σ_tot.** Use blocked CV with held-out v′ bands, J ranges, datasets and isotopologues.
   Report coverage: the fraction of held-out residuals within ±1σ, ±2σ, and the PIT histogram. Coverage
   is a published, tested property of the package, not a claim.
4. **Known risk regions for I₂ B–X**, to be covered by the above:
   - v′ near the B-state dissociation limit (long-range form, hyperfine-induced effects);
   - J beyond the observed maximum;
   - high v″ hot bands;
   - isotopologues with little or no data (mass scaling + BOB only).
5. **Reporting.** Publish maps of σ_tot over (v′, J′) and (v″, J″), with the data support overlaid. These
   replace Knöckel's per-region rule of thumb.

### 2.7 Compact storage of per-line uncertainties and correlations for a web app

**Source** [Source] Randomized low-rank factorizations (Halko, Martinsson & Tropp 2011) give near-optimal rank-r
approximations of large matrices at low cost.

**Options** [Analysis] P ≈ 150 parameters, N_L ≈ 10⁴–10⁵ levels, N_ℓ ≈ 10⁵–10⁶ lines incl. hyperfine
components.

| Option | Stored | Size (order of magnitude) | Exact? | Needs solver in browser? |
|---|---|---|---|---|
| **A. p̂ + factor of Σ_p, Jacobians on the fly** | p̂ (float64) + Cholesky/eigen factor of Σ_p (P×P) + discrepancy coefficients' posterior factor | ~0.2–1 MB | Exact, including all correlations | Yes. We need it anyway; Hellmann–Feynman gradients cost P quadratures per level |
| **B. Per-level low-rank factors** | Energy + factor Lₙ ∈ ℝ^r per level (r ≤ P + K_disc), with a line's factor = L_upper − L_lower | N_L·r·4 B, e.g. 10⁵ × 32 × 4 B ≈ 13 MB | Exact if r = full rank; truncated otherwise | No |
| C. Per-line low-rank factors | Line factor ∈ ℝ^r | ~10× more than B | Same as B | No |
| D. Posterior samples | S parameter vectors (S ≈ 100–200) | ~0.2 MB | Captures nonlinearity | Yes, S forward runs, fine for small windows |

- The parametric covariance G Σ_p Gᵀ has rank ≤ P. With a basis-function discrepancy (§2.4), the total rank
  is ≤ P + K_disc. So B is *exact* at modest rank. Truncation to r ≈ 20–40 (randomized SVD) can be checked
  against A.
- Store absolute positions in float64. 20 000 cm⁻¹ at 10⁻⁹ cm⁻¹ needs 13 significant digits. Factors can be
  float32 or even float16.
- **Recommendation:** A as the primary representation, since it is small, exact and uses the WASM solver we
  ship anyway. B as a precomputed cache for static downloads and the TUI. D as an optional "nonlinear
  check" mode.

---

## 3. Implementation stack (brief)

All repository checks were made against default branches on 2026-09-14/15 via the GitHub API. Versions are
from crates.io.

### 3.1 JAX vs NumPy/SciPy for Jacobians

**What JAX provides** [Source]
- JAX's derivative rule for `eigh` (`_eigh_jvp_rule` in `jax/_src/lax/linalg.py`) is documented in the source
  as "classic nondegenerate perturbation theory". Comments there note that degenerate eigenvalues need
  different treatment.
- 64-bit arithmetic requires the `jax_enable_x64` flag, which "by default … is set to False"
  (`docs/101/default_dtypes.md`).

**Assessment** [Analysis]
- For eigenvalues, autodiff through `eigh` and Hellmann–Feynman give the same derivative, first-order
  perturbation theory. The analytic HF route:
  - is exact, and costs one quadrature per level and parameter;
  - needs no differentiable eigensolver, and does not store the dense-eigh backward pass for ~600 J-blocks;
  - ports verbatim to Rust/WASM.

  It should be the production Jacobian.
- Autodiff is most useful for **∂V(R; p)/∂p on the grid** (MLR/damping/BOB functions, where hand-coded
  derivatives are error-prone). That step involves no eigensolver: `jax.jacfwd` of V(R; p), then an HF
  quadrature. It also serves as a unit test of hand-written derivatives.
- JAX is also the natural host for NumPyro/BlackJAX validation (§2.3), with the forward model wrapped as a
  `custom_jvp` that uses HF gradients.
- Proposal: a NumPy/SciPy core with JAX as an optional extra (`i2spec[bayes]`). This keeps the base install
  light and float64-by-default.

### 3.2 Rust eigensolvers that compile to wasm32-unknown-unknown (faer, nalgebra)

| | **nalgebra 0.35.0** (2026-05-24) | **faer 0.24.4** (2026-06-24) |
|---|---|---|
| Pure Rust | Yes. LAPACK bindings live in a separate `nalgebra-lapack` crate | Yes: README "implements low level linear algebra routines … in pure rust" |
| Symmetric eigensolver | `SymmetricEigen::new/try_new`, plus `symmetric_eigenvalues()`. Householder tridiagonalization (`SymmetricTridiagonal`, a public type) followed by implicit QR with Wilkinson shift (`src/linalg/symmetric_eigen.rs`) | `Mat::self_adjoint_eigen(side)` / `self_adjoint_eigenvalues(side)`, plus the lower-level `faer::linalg::evd::self_adjoint_evd`. Uses tridiagonalization then QR or divide-and-conquer (`tridiag_evd.rs`: `qr_algorithm`, `divide_and_conquer` with a QR fallback threshold) |
| Tridiagonal input | No dedicated tridiagonal-input API verified (*unknown*) | **Yes:** public `faer::linalg::evd::tridiagonal_self_adjoint_evd` |
| Subset of eigenpairs (lowest k) | Not found (*unknown*) | Not found (*unknown*) |
| wasm32-unknown-unknown | **Verified in CI:** the `build-wasm` job runs `rustup target add wasm32-unknown-unknown` and `cargo build --target wasm32-unknown-unknown`. no_std CI builds for `x86_64-unknown-none` and `thumbv7em-none-eabihf` | **Not in CI** (no wasm target in `run-tests.yml` or `code-quality.yml`). CI does run a no_std test crate (`faer-no-std-test`). The x86 GEMM backend (`private-gemm-x86`) is target-gated to `x86_64`. Issue #222, "Sparse cholesky crash on wasm32" (closed 2025-04-25), shows faer is used on wasm32 via wasm-pack; whether a fix landed was not checked |
| Default features | `std`, `macros` | `std`, `rayon`, `sparse-linalg`, `rand`, `npy`. For WASM use `default-features = false, features = ["std"]` (our recommendation, not yet compiled) |

- A local compile test was **not possible**: the Homebrew Rust toolchain here has no `wasm32-unknown-unknown`
  standard library, and rustup is absent.
- **Action item:** add a CI job that builds both crates' symmetric eigensolvers for wasm32-unknown-unknown
  and runs a numerical check in a WASM runtime, e.g. `wasm-bindgen-test` or Node.
- [Analysis] For a dense sinc-DVR Hamiltonian, either crate works. faer is expected to be faster on native;
  this is claimed by faer's benchmarks, which we did not re-run. nalgebra has the stronger wasm guarantee.
- A small hand-written symmetric-tridiagonal solver (implicit QL or bisection plus inverse iteration,
  ~300 lines) is a zero-dependency fallback. It fits only if the solver work stream chooses a banded or
  tridiagonal discretization.
- Recommendation: put the eigen-backend behind a trait in the core crate. Start with nalgebra for WASM
  certainty, benchmark faer natively, and switch per target if it pays.

### 3.3 PyO3/maturin and a shared Rust core

**Facts** [Source]
- Current versions: pyo3 0.29.2, maturin 1.15.0, wasm-bindgen 0.2.128 (crates.io, Aug–Sep 2026).
- maturin's `generate-ci` offers an `emscripten` platform option (`guide/src/distribution.md`), so
  Pyodide-style wheels are within reach. Pyodide specifics were not verified.

**Proposed layout** [Analysis] A Cargo workspace with:
- `i2spec-core`: pure Rust. Holds potentials, the radial solver, hyperfine Hamiltonian, HF Jacobians,
  covariance propagation and line-list generation. It has no PyO3 or wasm-bindgen dependencies, uses float64,
  and threading is an optional `parallel` feature.
- `i2spec-py`: PyO3 + maturin, with numpy arrays at the boundary.
- `i2spec-wasm`: wasm-bindgen for `wasm32-unknown-unknown`.

**When to start** [Analysis]
- The web app needs only the **forward model and uncertainty propagation**, never the fitter.
- Write the fitter in Python (SciPy least-squares/custom LM, optional JAX/NumPyro), since it will change
  constantly during research.
- Write the forward model in Rust **from the start of the forward-model freeze (v0.1)**, not as a later port.
  The existing Python prototype (sinc-DVR) becomes the reference oracle, with golden tests at the
  ≤ 1 kHz (≈3×10⁻⁸ cm⁻¹) level.
- The Python fitter then calls the Rust forward model through PyO3 for speed. This removes the
  two-implementation drift that a later port would create.

---

## 4. Recommendation

[Analysis throughout; it builds on the sourced material above.]

### 4.1 Fitting framework

1. **Data model with provenance.**
   - Each observation carries: type (FTS position | sub-Doppler absolute component | relative hyperfine
     splitting | fluorescence term/series), assignment (isotopologue, v′, J′, v″, J″, hyperfine labels),
     value, σ, dataset ID and conditions.
   - Each dataset declares its observation model (κ, δ, segments, pressure/power shifts) and the priors on
     those nuisance parameters.
2. **Stage 0: model-free network check.** A MARVEL-like sparse linear LSQ of level energies from all
   assigned transitions, with robust weights. It yields misassignment flags, initial calibration nuisances,
   and empirical level energies for later comparison with the DPF (§1.1, §1.3).
3. **Forward model.**
   - X and B potentials in the X-representation, re-expressed in an **orthogonal polynomial basis in X**, and
     joined to analytic long-range tails with theory-constrained Cₙ.
   - An MLR-type form is kept as an alternative, used for model-form spread.
   - BOB/adiabatic corrections for isotopologues.
   - Hyperfine parameters as R-dependent functions averaged with the same wavefunctions, with (v, J)
     correction terms only where needed.
   - A **blend-aware FTS observable** (§1.5).
4. **Estimator: MAP with a calibrated likelihood.**
   - Robust likelihood (Student-t, or Watson IRLS) with error model σ²_eff = (s_kσ)² + τ_k² + σ²_disc.
   - Smoothness/Tikhonov priors on the potentials, and priors on the nuisance parameters.
   - Solved by trust-region Levenberg–Marquardt on the SVD/QR of the weighted Jacobian, using the analytic
     **Hellmann–Feynman Jacobian**.
   - Nuisance parameters profiled by variable projection. Hyperfine–rovibronic blocks coupled through a
     Schur complement, making it one joint fit, not a sequential one.
   - Hyperparameters (s_k, τ_k, discrepancy amplitude and length scales) set by Laplace evidence, alternating
     with the MAP solve.
   - All down-weighting is audited.
5. **Model selection.**
   - **Blocked cross-validation**: leave out v′ bands, J ranges, datasets, isotopologues.
   - Laplace evidence and order scans.
   - Choose the smallest model whose blocked-CV errors agree with its predicted σ_tot.
6. **Publication.** Full-precision parameters, covariance factor and discrepancy posterior, in a
   machine-readable file. Le Roy sequential rounding only for printed tables, with a
   "reproduce every predicted line to < 0.1σ" test.

### 4.2 UQ approach

- **Production UQ:**
  - linear/Laplace parametric covariance, via HF Jacobians;
  - a basis-function GP **discrepancy on level energies** (X and B separately, so lines inherit correct
    correlations);
  - a model-form term from alternative forms and orders.

  Together these give σ_tot per line, the full line–line covariance on demand, and an interpolated /
  near- / far-extrapolation flag per line.
- **Validation before each release:**
  1. NUTS (NumPyro or BlackJAX, dense mass matrix from Laplace, HF gradients through `custom_jvp`), on the
     full problem if the forward model is fast enough, otherwise on a reduced one. It confirms that predicted
     marginals are Gaussian.
  2. Leave-one-dataset-out jackknife and parametric bootstrap.
  3. A published **coverage table** from blocked CV.
- **Storage/delivery:**
  - Option A (p̂ + Σ_p factor + discrepancy factor; Jacobians computed on the fly in WASM) is primary.
  - Option B (per-level low-rank factors, line factor = upper − lower) serves static line lists and the TUI.
  - Export per-line σ and per-state σ in HITRAN/ExoMol style for interchange.

### 4.3 Rationale

- DPF is the only strategy that predicts every (v, J) and isotopologue from a few hundred numbers with
  physical constraints. It is also what the baseline (Knöckel 2004) used, so results are directly
  comparable.
- A fitted error model replaces the fixed 3 MHz floor. This keeps kHz information where the model is
  adequate and down-weights it automatically where it is not, and the floor becomes a reportable, testable
  quantity.
- Laplace/linear UQ costs one Jacobian, reuses the fitting machinery and provides the evidence. It is
  *validated* by sampling rather than assumed.
- A discrepancy term is the only way to make interpolation *and* extrapolation uncertainties honest: purely
  parametric covariance ignores model-form error. Placing it on levels gives the right correlation
  structure and low rank.
- A shared Rust forward core started at v0.1: the web app needs only forward + UQ, and a single
  implementation avoids drift between Python and WASM.

### 4.4 Open questions / follow-ups

1. The exact form of Watson's (2003) robust weight. We recall it as 1/(σ² + δ²/3); verify from the paper.
2. Forward-model cost for one full evaluation (all J, both states, hyperfine), owned by the solver work
   stream. It decides whether full-problem NUTS and bootstrap are routine or offline-only.
3. Hyperfine parameterization near the B-state dissociation limit: is vibrational averaging of R-dependent
   hyperfine functions adequate, or are Bodermann-type (v, J) formulae plus a discrepancy needed?
4. How each legacy FTS dataset defines "line position" (peak vs centroid, handling of blends). Simulate the
   bias before fitting.
5. faer on wasm32-unknown-unknown is not in its CI. nalgebra is, but offers no verified tridiagonal-input or
   subset-eigenpair API. Add our own wasm32 CI test.
6. The anchor set: which comb-referenced lines fix the absolute scale, and how pressure/power shifts in
   sub-Doppler sets are modelled.
7. Local perturbations or couplings in B that a single-channel DPF cannot represent. The discrepancy model and
   the robust-weight audit will *flag* them, not fix them. A coupled-channel extension may be needed.
8. How the discrepancy prior amplitude should grow with distance from the data. Calibrate it on
   leave-J-range-out and leave-v′-band-out tests.

---

## References

All DOIs were checked against Crossref (or DataCite for arXiv) on 2026-09-14/15. Volume/page data are from
Crossref, and "–" marks pages that Crossref did not list. Content claims in the text come from abstracts,
publisher metadata, quoted passages in open-access citing papers, or repository files, as noted in place.

- Akaike, H. (1974). A new look at the statistical model identification. *IEEE Trans. Autom. Control* 19, 716–723. https://doi.org/10.1109/TAC.1974.1100705
- Arlot, S., Celisse, A. (2010). A survey of cross-validation procedures for model selection. *Statistics Surveys* 4. https://doi.org/10.1214/09-SS054
- Berninger, M., et al. (2013). Feshbach resonances, weakly bound molecular states, and coupled-channel potentials for cesium at high magnetic fields. *Phys. Rev. A* 87, 032517. https://doi.org/10.1103/PhysRevA.87.032517
- Birge, R. T. (1932). The calculation of errors by the method of least squares. *Phys. Rev.* 40, 207–227. https://doi.org/10.1103/PhysRev.40.207
- Bodermann, B., Knöckel, H., Tiemann, E. (2002). Widely usable interpolation formulae for hyperfine splittings in the ¹²⁷I₂ spectrum. *Eur. Phys. J. D* 19, 31–44. https://doi.org/10.1140/epjd/e20020052
- Brynjarsdóttir, J., O'Hagan, A. (2014). Learning about physical parameters: the importance of model discrepancy. *Inverse Problems* 30, 114007. https://doi.org/10.1088/0266-5611/30/11/114007
- Cabezas, A., Corenflos, A., Lao, J., Louf, R., et al. (2024). BlackJAX: Composable Bayesian inference in JAX. arXiv:2402.10797. https://doi.org/10.48550/arXiv.2402.10797
- Coxon, J. A., Hajigeorgiou, P. G. (2004). Direct potential fit analysis of the X¹Σ⁺ ground state of CO. *J. Chem. Phys.* 121, 2992–3008. https://doi.org/10.1063/1.1768167
- Efron, B. (1979). Bootstrap methods: another look at the jackknife. *Ann. Statist.* 7. https://doi.org/10.1214/aos/1176344552
- Feynman, R. P. (1939). Forces in molecules. *Phys. Rev.* 56, 340–343. https://doi.org/10.1103/PhysRev.56.340
- Foreman-Mackey, D., et al. (2013). emcee: The MCMC Hammer. *PASP* 125, 306–312. https://doi.org/10.1086/670067
- Furnstahl, R. J., Phillips, D. R., Wesolowski, S. (2015). A recipe for EFT uncertainty quantification in nuclear physics. *J. Phys. G* 42, 034028. https://doi.org/10.1088/0954-3899/42/3/034028
- Furtenbacher, T., Császár, A. G., Tennyson, J. (2007). MARVEL: measured active rotational–vibrational energy levels. *J. Mol. Spectrosc.* 245, 115–125. https://doi.org/10.1016/j.jms.2007.07.005
- Furtenbacher, T., Császár, A. G. (2012). MARVEL: measured active rotational–vibrational energy levels. II. Algorithmic improvements. *JQSRT* 113, 929–935. https://doi.org/10.1016/j.jqsrt.2012.01.005
- Gerstenkorn, S., Luc, P. (1979). Absolute iodine (I₂) standards measured by means of Fourier transform spectroscopy. *Rev. Phys. Appl.* 14, 791–794. https://doi.org/10.1051/rphysap:01979001408079100
- Golub, G. H., Pereyra, V. (1973). The differentiation of pseudo-inverses and nonlinear least squares problems whose variables separate. *SIAM J. Numer. Anal.* 10, 413–432. https://doi.org/10.1137/0710036
- Gordon, I. E., et al. (2022). The HITRAN2020 molecular spectroscopic database. *JQSRT* 277, 107949. https://doi.org/10.1016/j.jqsrt.2021.107949
- Gutenkunst, R. N., et al. (2007). Universally sloppy parameter sensitivities in systems biology models. *PLoS Comput. Biol.* 3, e189. https://doi.org/10.1371/journal.pcbi.0030189
- Halko, N., Martinsson, P. G., Tropp, J. A. (2011). Finding structure with randomness: probabilistic algorithms for constructing approximate matrix decompositions. *SIAM Rev.* 53, 217–288. https://doi.org/10.1137/090771806
- Hansen, P. C. (1998). *Rank-Deficient and Discrete Ill-Posed Problems.* SIAM. https://doi.org/10.1137/1.9780898719697
- HITRAN documentation, "Uncertainty codes", https://hitran.org/docs/uncertainties/ (accessed 2026-09-15).
- Hoffman, M. D., Gelman, A. (2011). The No-U-Turn Sampler: adaptively setting path lengths in Hamiltonian Monte Carlo. arXiv:1111.4246. https://doi.org/10.48550/arXiv.1111.4246
- Kennedy, M. C., O'Hagan, A. (2001). Bayesian calibration of computer models. *J. R. Stat. Soc. B* 63, 425–464. https://doi.org/10.1111/1467-9868.00294
- Knöckel, H., Bodermann, B., Tiemann, E. (2004). High precision description of the rovibronic structure of the I₂ B–X spectrum. *Eur. Phys. J. D* 28, 199–209. https://doi.org/10.1140/epjd/e2003-00313-4
- Le Roy, R. J. (1998). Uncertainty, sensitivity, convergence, and rounding in performing and reporting least-squares fits. *J. Mol. Spectrosc.* 191, 223–231. https://doi.org/10.1006/jmsp.1998.7646
- Le Roy, R. J. (2017). RKR1. *JQSRT* 186, 158–166. https://doi.org/10.1016/j.jqsrt.2016.03.030
- Le Roy, R. J. (2017). LEVEL. *JQSRT* 186, 167–178. https://doi.org/10.1016/j.jqsrt.2016.05.028
- Le Roy, R. J. (2017). dPotFit: a computer program to fit diatomic molecule spectral data to potential energy functions. *JQSRT* 186, 179–196. https://doi.org/10.1016/j.jqsrt.2016.06.002
- Le Roy, R. J. (2017). dParFit: a computer program for fitting diatomic molecule spectral data to parameterized level energy expressions. *JQSRT* 186, 197–209. https://doi.org/10.1016/j.jqsrt.2016.04.004
- Le Roy, R. J., Henderson, R. D. E. (2007). A new potential function form incorporating extended long-range behaviour: application to ground-state Ca₂. *Mol. Phys.* 105, 663–677. https://doi.org/10.1080/00268970701241656
- Le Roy, R. J., Huang, Y., Jary, C. (2006). An accurate analytic potential function for ground-state N₂ from a direct-potential-fit analysis of spectroscopic data. *J. Chem. Phys.* 125, 164310. https://doi.org/10.1063/1.2354502
- Le Roy, R. J., Dattani, N. S., Coxon, J. A., Ross, A. J., et al. (2009). Accurate analytic potentials for Li₂(X) and Li₂(A) from 2 to 90 Å, and the radiative lifetime of Li(2p). *J. Chem. Phys.* 131, 204309. https://doi.org/10.1063/1.3264688
- MacKay, D. J. C. (1992). Bayesian interpolation. *Neural Comput.* 4, 415–447. https://doi.org/10.1162/neco.1992.4.3.415
- McKemmish, L. K., et al. (2017). MARVEL analysis of the measured high-resolution rovibronic spectra of ⁴⁸Ti¹⁶O. *ApJS* 228, 15. https://doi.org/10.3847/1538-4365/228/2/15
- Melendez, J. A., et al. (2019). Quantifying correlated truncation errors in effective field theory. *Phys. Rev. C* 100, 044001. https://doi.org/10.1103/PhysRevC.100.044001
- Neufcourt, L., et al. (2018). Bayesian approach to model-based extrapolation of nuclear observables. *Phys. Rev. C* 98, 034318. https://doi.org/10.1103/PhysRevC.98.034318
- Pashov, A., Jastrzębski, W., Kowalczyk, P. (2000). Construction of potential curves for diatomic molecular states by the IPA method. *Comput. Phys. Commun.* 128, 622–634. https://doi.org/10.1016/S0010-4655(00)00010-2
- Phan, D., Pradhan, N., Jankowiak, M. (2019). Composable effects for flexible and accelerated probabilistic programming in NumPyro. arXiv:1912.11554. https://doi.org/10.48550/arXiv.1912.11554
- Pickett, H. M. (1991). The fitting and prediction of vibration-rotation spectra with spin interactions. *J. Mol. Spectrosc.* 148, 371–377. https://doi.org/10.1016/0022-2852(91)90393-O
- Qu, Q., Yurchenko, S. N., Tennyson, J. (2022). A method for the variational calculation of hyperfine-resolved rovibronic spectra of diatomic molecules. *J. Chem. Theory Comput.* 18, 1808–1820. https://doi.org/10.1021/acs.jctc.1c01244
- Rasmussen, C. E., Williams, C. K. I. (2005). *Gaussian Processes for Machine Learning.* MIT Press. https://doi.org/10.7551/mitpress/3206.001.0001
- Salami, H., Ross, A. J. (2005). A molecular iodine atlas in ascii format. *J. Mol. Spectrosc.* 233, 157–159. https://doi.org/10.1016/j.jms.2005.06.002
- Schwarz, G. (1978). Estimating the dimension of a model. *Ann. Statist.* 6. https://doi.org/10.1214/aos/1176344136
- Seto, J. Y., Le Roy, R. J., Vergès, J., Amiot, C. (2000). Direct potential fit analysis of the X¹Σg⁺ state of Rb₂: Nothing else will do! *J. Chem. Phys.* 113, 3067–3076. https://doi.org/10.1063/1.1286979
- Takekoshi, T., et al. (2012). Towards the production of ultracold ground-state RbCs molecules: Feshbach resonances, weakly bound states, and coupled-channel models. *Phys. Rev. A* 85, 032506. https://doi.org/10.1103/PhysRevA.85.032506
- Tellinghuisen, J. (1990). Statistical uncertainties in RKR potentials: an exercise in error propagation. *J. Mol. Spectrosc.* 141, 258–264. https://doi.org/10.1016/0022-2852(90)90162-J
- Tellinghuisen, J. (2001). Statistical error propagation. *J. Phys. Chem. A* 105, 3917–3921. https://doi.org/10.1021/jp003484u
- Tennyson, J., et al. (2020). The 2020 release of the ExoMol database. *JQSRT* 255, 107228. https://doi.org/10.1016/j.jqsrt.2020.107228
- Tennyson, J., et al. (2024). The 2024 release of the ExoMol database. *JQSRT* 326, 109083. https://doi.org/10.1016/j.jqsrt.2024.109083
- Tiesinga, E., et al. (2021). CODATA recommended values of the fundamental physical constants: 2018. *Rev. Mod. Phys.* 93, 025010. https://doi.org/10.1103/RevModPhys.93.025010
- Tóbiás, R., Furtenbacher, T., Simkó, I., Császár, A. G., et al. (2020). Spectroscopic-network-assisted precision spectroscopy and its application to water. *Nat. Commun.* 11, 1708. https://doi.org/10.1038/s41467-020-15430-6
- Transtrum, M. K., Machta, B. B., Sethna, J. P. (2010). Why are nonlinear fits to data so challenging? *Phys. Rev. Lett.* 104, 060201. https://doi.org/10.1103/PhysRevLett.104.060201
- Transtrum, M. K., Machta, B. B., Sethna, J. P. (2011). Geometry of nonlinear least squares with applications to sloppy models and optimization. *Phys. Rev. E* 83, 036701. https://doi.org/10.1103/PhysRevE.83.036701
- Vargas-Hernández, R. A., Sous, J., Berciu, M., Krems, R. V. (2018). Extrapolating quantum observables with machine learning. *Phys. Rev. Lett.* 121, 255702. https://doi.org/10.1103/PhysRevLett.121.255702
- Watson, J. K. G. (2003). Robust weighting in least-squares fits. *J. Mol. Spectrosc.* 219, 326–328. https://doi.org/10.1016/S0022-2852(03)00100-0
- Yin, Y., Yurchenko, S. N., Tennyson, J. (2026). Hyperfine-resolved variational nuclear motion spectra of homonuclear diatomic molecules. *J. Chem. Phys.* 165, 092501. https://doi.org/10.1063/5.0346970

**Software, repositories and data files** (checked 2026-09-14/15):
- nalgebra: https://github.com/dimforge/nalgebra (`src/linalg/symmetric_eigen.rs`, `.github/workflows/nalgebra-ci-build.yml`); crates.io v0.35.0.
- faer: https://github.com/sarah-quinones/faer-rs (`faer/Cargo.toml`, `faer/src/linalg/evd/`, `faer/src/linalg/solvers.rs`, `.github/workflows/run-tests.yml`, issue #222); crates.io v0.24.4.
- JAX: https://github.com/jax-ml/jax (`jax/_src/lax/linalg.py`, `docs/101/default_dtypes.md`).
- maturin: https://github.com/PyO3/maturin (`guide/src/distribution.md`); crates.io v1.15.0. PyO3: crates.io v0.29.2. wasm-bindgen: crates.io v0.2.128.
- ExoMol dataset definitions: https://www.exomol.com/db/TiO/48Ti-16O/Toto/48Ti-16O__Toto.def.json and https://www.exomol.com/db/CO/12C-16O/Li2015/12C-16O__Li2015.def.json.
