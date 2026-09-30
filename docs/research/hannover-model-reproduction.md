# Reproducing the Hannover I₂ B–X potential model (milestone 1)

*2026-09-15. Code: `prototypes/hannover_baseline.py`.*

## Summary

We rebuilt the Hannover rovibronic model (the potentials inside IodineSpec) from published parameters. It solves correctly with our own radial solver. For the X state, **use the 2008 parameter set** (Salumbides et al., EPJD 47, 171 (2008), Table 1; open copy on the VU Amsterdam repository). The 2004 Table 4 (Knöckel, Bodermann & Tiemann, EPJD 28, 199) is known to contain misprints. The Hannover IQO web page says so, but its "look here" erratum link is dead and the Wayback Machine was offline. Our tests pin down where the misprints are:

* **X-state column of the 2004 Table 4.**
  * The level spacings come out wrong by up to 1.3 cm⁻¹ at v″ = 17: NIR lines are off by 17–44 GHz.
  * The printed extension constants (A_I, B_I, A_O, B_O) don't join the printed power series smoothly.
  * No single-coefficient correction brings the RMS below 1.7 GHz, so more than one entry is wrong (or the column is scrambled).
* **B-state A_O and B_O in 2004** don't join smoothly either. The 2008 values do.
* **Sign of the long-range exponential.** The 2008 paper defines it explicitly: V = De − Σ Cₙ/Rⁿ − A_O·exp(−B_O(R − R_O)).

## Implementation

* **Potentials.** X-representation potentials: V = Σ aᵢXⁱ with X = (R − Rm)/(R + b·Rm), plus the exponential inner wall and the dispersion + exponential outer branch.
* **Rotational term.** B-state nonadiabatic α(R) = (2Rm/(R+Rm)) Σ αᵢXⁱ in the centrifugal term. For the reference isotopologue ¹²⁷I₂, V_corr = 0.
* **Solver.** sinc-DVR (Colbert–Miller), dense, with the lowest 50 (B) or 20 (X) eigenvalues taken per J.
  * Grid: B on 2.35–7.0 Å with h = 0.005 Å; X on 2.10–4.0 Å with h = 0.004 Å.
  * The full test script runs in about 3 s on a laptop.
* **Constants.**
  * ħ²/(2u·Å²) = 16.857629192 cm⁻¹ (CODATA 2018 via SciPy).
  * m(¹²⁷I) = 126.9044719 u (AME2020, atomic mass).
  * c = 29 979.2458 MHz per cm⁻¹.

## Checks (2008 parameters)

**Joins.** The series matches the tabulated extension constants at all four joins.

| Join | Series value (cm⁻¹) | Extension value (cm⁻¹) | Difference |
|---|---|---|---|
| B, R_I = 2.647 Å | 19605.26886 | 19605.26903 | 1.7×10⁻⁴ |
| B, R_O = 4.9 Å | 19810.452297 | 19810.452315 | 2×10⁻⁵ |
| X, R_I = 2.40 Å | 4580.93282 | 4580.93294 | 1.2×10⁻⁴ |
| X, R_O = 3.30 Å | 6950.242940 | 6950.242966 | 3×10⁻⁵ |

**Absolute energies** against Gerstenkorn, Luc & Le Roy (1991):

| Quantity | This model | GLL 1991 |
|---|---|---|
| E_X(v=0, J=0) | 107.0995 cm⁻¹ | D_e − D₀ = 107.101 |
| T₀ | 15724.5872 cm⁻¹ | 15724.587 |

**NIR bands.** Rovibronic R(J″) lines of the potential model minus the independent local Dunham model (2004 Table 1, stated 1σ below 200 kHz), in MHz:

| Band | J″=20 | 50 | 80 | 100 | 120 | 150 | 180 | 210 |
|---|---|---|---|---|---|---|---|---|
| 0-12 | −48.4 | −44.9 | −38.4 | −32.7 | −26.1 | −16.1 | −7.3 | −0.9 |
| 0-13 | −10.9 | −10.6 | −9.4 | −8.1 | −6.2 | −2.7 | 1.0 | 4.1 |
| 0-14 | −2.1 | −2.8 | −3.5 | −3.5 | −2.9 | −0.9 | 2.6 | 6.7 |
| 0-15 | −0.4 | −1.3 | −2.3 | −2.6 | −2.5 | −1.0 | 2.5 | 7.0 |
| 0-16 | 2.6 | 1.6 | −0.1 | −1.5 | −3.0 | −4.9 | −5.8 | −6.6 |
| 0-17 | −1.5 | −3.1 | −6.8 | −10.9 | −16.5 | −27.7 | −41.3 | −57.8 |

* The central bands agree at the few-MHz level, which is what the Hannover papers state for the global model.
* The edges (0-12 at low J, 0-17 at high J) differ by tens of MHz. Two explanations are possible:
  * the local Dunham model extrapolating beyond its data; Table 1 has no S₃₀ term, and its data coverage per band is uneven;
  * the 2008 refit changing the X potential. Liao et al. (JOSA B 27, 1208, 2010) found deviations "larger than expected" at 750–780 nm and made a new local NIR model.
* To resolve this, compare against measured NIR frequencies from the stage 2 catalog.

**Absolute hyperfine components.** These are coarse checks until the hyperfine model exists. The differences below are model rovibronic frequency minus the measured component; they should equal the component's offset from the line centre, which is up to about ±500 MHz.

| Line | Model − component (MHz) | Expectation |
|---|---|---|
| R(56)32-0 a10, 532 nm (BIPM) | −92.4 | a10 near the middle of the 15-component pattern |
| R(127)11-5 a16 (f), 633 nm (BIPM) | −234.8 | Consistent with a16 sitting in the upper part of the 21-component pattern |
| R(37)16-1 a1, 578 nm (Hong et al. 2009) | +509.0 | Consistent with a1 being the lowest component |

**Numerical convergence.** The joins are only C¹: value and slope match, but the second derivative jumps. That destroys the spectral convergence of the sinc-DVR for levels whose wavefunctions reach a join.
* v = 0 levels move by less than 1 kHz when h → 0.8h.
* X v″ = 15 moves by about 20 kHz, and B v′ = 45 by about 60 kHz.
* **B v′ = 32, J′ = 57**, whose inner turning point lies near R_I = 2.647 Å. Relative to h = 0.005 Å, the energy is
  +48, −222, −305, −236, −247 and −239 kHz at h = 0.004, 0.0033, 0.0025, 0.002, 0.0016 and 0.00125 Å.
  The convergence oscillates, and the value settles about −240 ± 10 kHz from the default only for h ≤ 0.002 Å.

**Fixed, 2026-09-15.** `src/i2spec/bspline.py` expands the wavefunctions in B-splines with knots at R_I and R_O (compare Derevianko et al., arXiv:0806.1368):
* At those knots the basis is C³, the regularity of the exact eigenfunctions where V″ jumps.
* Matrix elements are integrated piecewise, so no quadrature interval straddles a join.
* **B(32, 57):** with order 10 and 0.01 Å breakpoint spacing it converges to 1×10⁻⁹ cm⁻¹ (≈ 30 Hz). Orders 8, 10 and 12, at spacings 0.005–0.012 Å, all agree. The converged value is 4.7 kHz from the fine-step DVR estimate, which itself scatters by ±10 kHz.
* **X(15, 100):** converged to 1×10⁻¹² cm⁻¹.
* **Analytic check:** the Morse levels come out to 10⁻¹¹ cm⁻¹.

The B-spline solver is now the default for level energies. The sinc-DVR stays for intensities, which need eigenvectors on a grid shared by X and B.

**Line-list positions, fixed 2026-09-17.** Until then the master line list (`intensity.master_line_list`) took its positions from those DVR eigenvalues, so the lookup, the TUI and the web app showed DVR positions while `RovibronicModel.transition` gave B-spline ones. The error was larger than the 0.24 MHz at R(56) 32-0 suggested:
* **B state:** up to ~3 MHz, the DVR's convergence error with a 0.005 Å step (both solvers used the same 7 Å box).
* **X state, v″ ≳ 23: up to 86 GHz.** The shared grid started at 2.33 Å, a hard wall just inside the inner turning points of high v″ (2.39 Å at v″ = 25, 2.35 Å at v″ = 39). Levels came out 51 MHz high at v″ = 25, 1.6 GHz at 30 and 12.6 GHz at 34. Line strengths at v″ = 30–39 moved by up to 6.8% once the wall was removed.

The fix has two parts:
* **Positions**, lower-level energies and the partition function now come from B-spline levels (`POSITION_GRIDS`), matched to the DVR levels by energy, not by index. The list agrees with `RovibronicModel.transition` to 1 Hz at R(56) 32-0, and to under 1 kHz on a random sample.
* **X is solved on the B grid's points extended inward to 2.10 Å.** B cannot follow it there: its nonadiabatic α(R) is a fitted polynomial that diverges below ~2.35 Å (+19 at 2.33 Å, −19 000 at 2.25), and on a 2.20 Å grid the B levels collapse by 10²² cm⁻¹ at J = 150. B wavefunctions vanish below 2.33 Å, so the overlaps are unchanged.

Against the previous list, positions moved by a median 2 kHz for v″ ≤ 16, 176 kHz for v″ = 17–22, 22 MHz for v″ = 23–29 and 3.6 GHz for v″ = 30–39. Strengths for v″ ≤ 22 moved by under 0.01%.

**A caveat on `RovibronicModel.energy(state, v, J)` that this exposed.** At high J (seen at J″ = 226–330), B levels above the dissociation asymptote are held only by the centrifugal barrier. There, box states of the solver's box interleave with the real quasi-bound levels, and `v` counts them. The level labelled v′ = 32 at J′ = 277 sits at 20 495 cm⁻¹ in the default 8 Å box and 20 337 cm⁻¹ in a 12 Å box, while the real quasi-bound level (20 533.72 cm⁻¹) is index 32, 50 and 88 in 7, 10 and 14 Å boxes. The line list is not affected: it keeps only levels below the first box state and matches energies. `RovibronicModel.transition` is: 5 of 400 randomly sampled lines (v′ ≤ 59) disagreed with the list by 0.6–1.7 THz for this reason, all far above any measured data. Not fixed.

**Physical constants.** SciPy 1.13 uses CODATA 2018 and SciPy 1.18 uses CODATA 2022. The change in the atomic mass constant (1.4×10⁻⁹) shifts B–X lines by about 40 kHz. The package now pins h, c (exact) and u (CODATA 2022) in `constants.py`. To reproduce IodineSpec at the kHz level, we would also need the constants and masses Hannover used.

**Mass convention.** Shifting m(¹²⁷I) by ±4×10⁻⁶ u moves NIR lines by about 1 MHz and leaves T₀ unchanged. The Hannover mass value must be matched to reproduce IodineSpec at the sub-MHz level.

## Milestone 2: hyperfine structure of ¹²⁷I₂

*Code: `src/i2spec/hyperfine.py`, `src/i2spec/hfs_params.py`. Tests: `tests/test_hyperfine.py`, `tests/test_bipm_532.py`.*

**Method.**
* A four-term effective Hamiltonian (eqQ, C, d, δ), set up in the coupled basis |J (i₁i₂) I; F⟩ and diagonalized for each F.
* The basis includes the rotational neighbours J ± 2, with rotational energies from the potential model, and ΔI = ±2 mixing.
* Matrix elements come from general tensor algebra (6j and 9j symbols via SymPy). They are checked against an independent brute-force construction in the uncoupled basis, built from explicit spin matrices and spherical-harmonic quadrature.
* Parameters come from Bodermann, Knöckel & Tiemann, EPJD 19, 31 (2002), eqs. (10)–(15). The PDF's text layer and page images agree on every coefficient.

**Conventions, fixed by the BIPM tables.** The literature doesn't state these explicitly, and each one decides the result:

| Choice | Alternative, R(56)32-0 intervals |
|---|---|
| EQ sign: within ΔJ = 0, H_EQ = −eqQ Σₙ[3(Iₙ·J)² + (3/2)(Iₙ·J) − Iₙ²J²]/[2i(2i−1)(2J−1)(2J+3)] (the Townes–Schawlow sign) | Opposite sign: ~100 MHz errors (mirrored pattern) |
| TSS: d·[3(I₁·J)(I₂·J) + 3(I₂·J)(I₁·J) − 2(I₁·I₂)J(J+1)]/[(2J−1)(2J+3)], i.e. d times the dipolar form I₁·I₂ − 3(I₁·n̂)(I₂·n̂) | With an extra factor of 5: 0.74 MHz rms. Without the d term: 0.19 MHz rms |
| ΔJ = ±2 couplings included | Without them: 0.55 MHz rms. Adding ±4 changes nothing |
| Eqs. (12)–(13) read as the absolute δ_B and d_B; the paper says to add δ_X and d_X | Checked against the paper's Table 4: mean residual +0.1 kHz for our reading vs +3.8 kHz (δ); +0.15 vs +1.7 kHz (d) |

**Results** (Hannover 2008 potentials plus Bodermann 2002 hyperfine formulae):

| Test | Result | Stated accuracy of the formulae |
|---|---|---|
| R(56)32-0: 13 intervals from a10 (BIPM, u = 1.5 kHz) | rms 23.7 kHz, max 53.7 kHz. F−J=0 components: a1 +4.7 kHz, a15 −11.7 kHz | 20–30 kHz for F−J=0 components; up to 1 MHz for the others at high J |
| R(87)33-0: 21 intervals from a1 (u = 2 kHz) | rms 89 kHz, max 147 kHz. The F−J=0 components a2, a13 and a20 agree with each other to 2 kHz but sit +83 kHz from a1 | Same |
| R(56)32-0 a10, absolute | Model − CIPM = **+2.00 MHz** with the B-spline solver (+2.24 MHz with the sinc-DVR, which is 0.23 MHz off at v′ = 32) | Rovibronic ≈ 1.5 MHz (1σ; Salumbides 2008) |

**Limitations.**
* **Speed:** the Wigner symbols are now floating-point Racah formulas (`wigner.py`), not SymPy. A ¹²⁷I₂ line takes a small fraction of a second with ΔJ = ±2 mixing, and an isotopologue line 1.3–1.5 s.
* **Isotopologues:** implemented with Salumbides et al. 2006 (`isotopologue-hyperfine.md`). Intervals agree with the BIPM tables to ≤ 1 MHz (`bipm-hyperfine-tables.md`). The isotope shifts are still open; see the next section.
* **Range:** v′ ≤ 43, v″ ≤ 17 and J″ ≤ 200 for the Bodermann 2002 formulae. Salumbides 2006 extends v′ to 53.

## Isotope shifts of ¹²⁹I₂ and ¹²⁷I¹²⁹I (resolved 2026-09-18 by a constant in V_ad)

*A constant added to the B-state adiabatic correction, +0.02154 cm⁻¹, takes the isotopologue lines
from 7.06 to 2.36 MHz rms against a data floor of 1.1–1.2 MHz; the default parameter set `i2spec2026a`
carries it (`docs/design/parameter-sets.md`). The section below is the analysis as it stood.*

*Investigated 2026-09-15. The scripts were throwaway and are not in the repository.*

**Symptom.** Model − BIPM for the isotope-shift intervals in `tests/test_bipm_other.py` (strict xfail):

| Interval | Model − BIPM (MHz) |
|---|---|
| ¹²⁹I₂ P(54) 8-4 − ¹²⁷I₂ R(127) 11-5 | −10.95 |
| ¹²⁹I₂ P(33) 6-3 − ¹²⁷I₂ R(127) 11-5 | −11.86 |
| ¹²⁹I₂ P(110) 10-2 − ¹²⁷I₂ R(47) 9-2 | −7.39 |
| ¹²⁹I₂ R(113) 14-4 − ¹²⁷I₂ R(47) 9-2 | −13.91 |
| ¹²⁷I¹²⁹I P(33) 6-3 − ¹²⁹I₂ P(54) 8-4 | +4.51 |

Knöckel 2004 lists exactly these BIPM isotopologue lines in its fit data. Salumbides 2008 kept that data base and used a 3 MHz minimum uncertainty. So the published model ought to reproduce these intervals to a few MHz.

**Salumbides 2006's own spectra show the same offsets.** Fig. 2 of Salumbides 2006 is a vector graphic with an absolute MHz axis, calibrated on ¹²⁷I₂. Method:
* Fit our hyperfine patterns (Lorentzian, fitted width 12–13 MHz) to its traces.
* The ¹²⁷I₂ lines sit −1.3 MHz (P34 11-3) and −3.3 MHz (R29 9-2) from the model.
* Correct the upper trace (the ¹²⁹I₂ cell) by +37.5 MHz. Salumbides 2006 and 2008 both state that the pump-offset AOM shifts that setup's resonances 37.5 MHz to lower frequencies. That the plotted trace is uncorrected is inference made here.

Model − data relative to ¹²⁷I₂ P34(11-3) is then:
* ¹²⁹I₂: −13.7 MHz (P144 13-3), −8.7 (R110 12-3), −5.3 (R34 11-3), −11.1 (R149 13-3);
* ¹²⁷I¹²⁹I: −2.9 MHz (R37 11-3), i.e. +6.8 MHz relative to the ¹²⁹I₂ mean.

The sign, size and isotopologue ratio match the BIPM comparison.

**Ruled out:**
* **Numerics.** B-spline and DVR agree to better than 1 kHz on these lines.
* **Masses.**
  * Changing m(¹²⁹I) by 1×10⁻⁵ u moves the intervals by less than 0.4 MHz.
  * Nuclear instead of atomic masses moves them by hundreds of MHz.
  * The ±2 MHz agreement of ¹²⁷I₂ absolute frequencies already pins the atomic-mass convention.
* **BOC conventions.** The code follows Salumbides 2008 eqs. (5)–(9) as printed:
  * α = (μ_ref/μ)(2R_m/(R+R_m)) Σαᵢ Xⁱ;
  * V_corr = (1 − μ_ref/μ)(2R_m/(R+R_m))⁵ Σvᵢ Xⁱ;
  * X as in eq. (2);
  * applied to the B state only, with no extensions.

  Every alternative tried is much worse:

  | Variant | Change in the intervals (MHz) |
  |---|---|
  | Knöckel 2004 form, power 1 in V_ad | −117 to −202 |
  | Power 4 or 6 in V_ad | −58 to +29 |
  | X with b = 0, −0.5 or 1 | hundreds |
  | Prefactor power 0 or 2 in α | ±70 |
  | α not scaled, or scaled by (μ_ref/μ)² | ±60 |
* **The paper's own number.** Our V_corr gives T_e(¹²⁷I₂) − T_e(¹²⁹I₂) = 96.1 MHz, against the 94(11) MHz quoted in Salumbides 2008.

**What fits: a uniformly weaker adiabatic correction.**
* **All of V_ad × 0.919:** all ten comparisons (5 BIPM, 5 from Fig. 2) fit with rms 2.5 MHz and max 4.7 MHz. T_e(¹²⁷I₂) − T_e(¹²⁹I₂) becomes 88.3 MHz, still inside 94(11).
* **v₀ alone** (−0.2097 → −0.1858 cm⁻¹): rms 2.7 MHz, max 5.0 MHz.
* **No single α coefficient** does better than rms 5.2 MHz.

**Conclusion.** No paper documents the cause. Two candidates:
* a misprint in the Table 1 V_ad coefficients of Salumbides 2008 (the 2004 Table 4 had misprints);
* an unstated convention in the Hannover code for the isotope mass factor.

The code keeps the 2008 parameters as printed, and `test_isotope_shifts` stays a strict xfail. The proper fix is the Phase B refit, which will include these data. ¹²⁷I₂ is unaffected, because V_corr = 0 and α is unscaled for the reference isotopologue.

## Next steps

**Done since this list was first written:**
* the join-aware B-spline solver;
* fast Wigner symbols;
* M3 intensities and spectra (`spectra-validation.md`, `continuum-model.md`);
* isotopologue hyperfine structure;
* all BIPM tables (`bipm-hyperfine-tables.md`, and the data sets in `data/observations/`);
* the mass values, which are ruled out as the cause of the isotope-shift error (a change of 10⁻⁵ u moves it by < 0.4 MHz).

**Still to do:**
1. The Phase B refit. It has to fix three things with the published parameters:
   * the B levels above v′ ≈ 44, which are wrong by up to 2.2 cm⁻¹ (`spectra-validation.md`);
   * the isotope shifts (above);
   * the v′-dependent errors of −6 to +6 MHz at 532 nm (BIPM inter-line intervals).
2. Measured NIR and visible data from the stage 2 catalog, entered as observation data sets (`docs/design/observations.md`).

There will be no author contact, so the 2004 Table 4 erratum stays unavailable and the 2008 parameter set is used. The reading of eqs. (12)–(13) rests on the Table 4 evidence above.

## The outer-branch sign: the 2004 and 2008 papers disagree

*Found 2026-09-15, by plotting the potential: a wrong sign shows up there directly.*

The two papers print the outer extension differently:

* **Knöckel 2004, p. 6:** V(R) = D_e − Σ C_n/Rⁿ **+** A_O exp(−B_O(R − R_O))
* **Salumbides 2008, eq. (4):** V(R) = D_e − Σ C_n/Rⁿ **−** A_O exp(−B_O(R − R_O))

Only the minus sign is consistent with the tabulated A_O values, for either state:

| State | R_O (Å) | Series at R_O | Printed 2004 form, +A_O | 2008 form, −A_O |
|---|---|---|---|---|
| X | 3.3 | 6950.2430 | 9057.4710 (**+2107.2**) | 6950.2430 (exact) |
| B | 4.9 | 19810.4523 | 19794.2417 (**−16.2**) | 19810.4523 (exact) |

The slopes agree the same way: at the X join the series gives 12042.211 cm⁻¹/Å, the +A_O form 9389.147, and the −A_O form 12042.211.

* **i2spec implements the minus sign** (`potentials.py`), so this is a defect of the 2004 text, not of our model.
* The damage is asymmetric. A_O is large and positive for X, so the 2004 form breaks the ground state by 2107 cm⁻¹ — obvious in any plot. For B it is small and negative, so the error is only 16 cm⁻¹ and easy to miss.
* `with_continuous_extensions()` rederives A_I, B_I, A_O and B_O from continuity alone. They agree with the printed values to about 1×10⁻⁵, which confirms the **table** is sound and only the 2004 **equation** is wrong.

This is the fourth defect found in this family of papers, after the Salumbides 2006 Ω = 0 sign flip, the Bodermann 2002 Δδ/Δd mislabelling, and the 2004 Table 4 misprints.
