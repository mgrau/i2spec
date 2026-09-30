# The global fit (M5, Phase B)

`src/i2spec/fitting.py`. Refits the potentials to everything in `data/observations/`, replacing the
published Hannover parameters with our own.

```python
from i2spec.fitting import GlobalFit
g = GlobalFit()            # loads every data set, reduces each observation to level energies
result = g.run()           # scipy least_squares, analytic Jacobian, robust loss
```

## What varies, and what does not

**Varies:** the power-series coefficients `a_i` of the X and B X-representation potentials — 14 and
32, so 46 parameters.

**Derived, not fitted:** the inner wall (A_I, B_I) and the outer branch (A_O, B_O). After every
change to `a_i` they are recomputed from C¹ continuity at R_I and R_O, via
`XRepPotential.with_continuous_extensions()`. The potential is therefore continuous *by
construction*; there are no continuity constraints for the optimiser to violate.

**Held fixed — and this is a result, not a convenience:**

- **The hyperfine parameters.** `bodermann1998b` measures an absolute frequency *and* a hyperfine
  splitting for the same lines at v″ = 16–17. The splitting is reproduced to **36 kHz**; the
  positions are out by **2.5 and 5.8 MHz**. Both components of R(42) 0-17 miss by the same +2.53
  MHz. The error is in the term values, not the hyperfine Hamiltonian, so freeing the hyperfine
  constants here would only let them soak up potential error.
- **The Born–Oppenheimer corrections.** Seven observations in the whole collection are not ¹²⁷I₂.

Each observation's hyperfine offset is computed once from the starting model and carried as an
additive constant.

## How an observation becomes a residual

Every datum is reduced to a signed combination of level energies plus a constant:

| kind | model value |
|---|---|
| `frequency` | `[E_B(v′,J′) − E_X(v″,J″)]·c + offset(component)` |
| `interval`, different lines | the same, minus the reference line's |
| `interval`, same line | **the rovibronic centres cancel exactly** |

That last row matters. A splitting measured within one line constrains the hyperfine Hamiltonian and
says nothing about the potentials while the hyperfine parameters are fixed. Of 1510 observations,
**854 constrain the potentials and 656 are hyperfine-only**; the latter are held in
`GlobalFit.hyperfine_only` for reporting and excluded from the fit. Stage 2 fits them: see
`docs/design/hyperfine-fit.md`.

## The Jacobian is analytic

By Hellmann–Feynman, `∂E(v,J)/∂p = ⟨ψ_vJ| ∂V/∂p |ψ_vJ⟩`, so **one eigensolve per (isotopologue,
state, J) gives the derivative with respect to every parameter**. `BSplineSolver.wavefunctions(J)`
returns ψ on the quadrature grid, normalised so that `Σ W ψ² = 1`.

`∂V/∂p` itself is taken by central differences *of the potential* on that grid — two function
evaluations per parameter, microseconds each — which automatically includes the induced change in
the re-derived wall and tail.

Measured: the Jacobian costs **1.3× one forward model** instead of 46×, a **37× speed-up** in situ
(412× on a single state in isolation), and agrees with finite differences to ~10⁻⁵ relative.

This is why the fitter stays in Python. Profiling shows **97–100% of the runtime is inside LAPACK**;
matrix assembly is 1%. A port to Rust or Julia would call the same routines. The remaining wins are
algorithmic: the B-spline Hamiltonian has bandwidth 9 at n = 583 (2.9% dense), so a banded solver
(`dsbgvx`) should give about 10×, and the 293 solves per evaluation are embarrassingly parallel.

## The prior is not optional

The data constrain only **29 of the 46 directions**. v′ = 19–31 and 38–57 and v″ = 7–10 have no data
at all, and the B series alone has 32 coefficients. An unregularised fit would let the unconstrained
directions drift to fit noise, and the potential would go wrong exactly where nothing is watching.

So each coefficient carries a Gaussian prior centred on the published value with
`sigma = prior_rel·|a_i| + prior_floor`. This is a genuine Bayesian prior, not a numerical dodge: the
published values carry the information from the Knöckel 2004 fit set, **which no longer exists**
(see `docs/research/data-availability.md`). The prior is how that lost data still constrains the result.

The optimiser works in scaled units `u = (a − a0)/sigma_prior`, which makes the prior block the
identity. Unscaled, the Jacobian's condition number is 2.5×10¹⁷ and its numerical rank 23/46; scaled,
it is 8.8×10¹⁰ at full rank 46/46, with the smallest singular value exactly 1 — the prior holding up
the directions the data cannot see.

## Robust loss

The default is `soft_l1`. One datum — R(98) 58-1 d6, at v′ = 58 — sits **41 GHz** from the model,
because the published B potential is known to fail above v′ ≈ 44. Its stated uncertainty is 1 MHz.
Under a quadratic loss that single point would contribute ~10⁹ to χ² and drive the entire fit, and
with no data between v′ = 44 and 58 it would distort the B potential across a region nothing else
constrains.

## Result of the first real run: the series coefficients alone cannot reach high v''

With `nesterenko2019` in the set (v'' = 53 and 54, 2-235 kHz, starting 14-21 cm-1 out) the fit
converges and **does essentially nothing** — every data set changes by 0%. That is not a solver
failure; there is no improving direction. Three measurements establish it:

1. **The step that would fix v'' = 53/54 alone is 8x10^9 prior sigmas**, and it would move the bulk
   data by 2x10^7 MHz rms and violate the shape prior by 2x10^10 sigma.
2. **High-v'' levels respond to the coefficients just as strongly as low-v'' ones** — the ratio of
   response magnitudes is 1.095, so it is not a sensitivity problem.
3. **The responses are 94% collinear.** The cosine between the mean response direction of the
   v'' = 53/54 data and that of everything else is +0.940. Raising v'' = 53 by 14 cm-1 necessarily
   drags v'' <= 17 with it, and the orthogonal 6% can only be reached through an enormous,
   near-cancelling excursion.

**The missing handle is the long-range tail.** Stage 1 varies only the series coefficients a_i; the
X-state De and C_6, C_8, C_10 and the join radius R_O are held at their published values. Those act
where the high-v'' wavefunctions have amplitude and the low-v'' ones have none, which is exactly the
orthogonal direction the fit cannot otherwise reach. Knoeckel 2004 Table 4 takes the X long-range
coefficients from Bacis et al. 1986 (`10.1016/0022-2852(86)90180-3`), which is now in hand.

So the next step is not more optimiser tuning; it is to free the X long-range parameters.

## Freeing the long-range tail: it breaks the degeneracy, but cannot supply the correction

De and C6, C8, C10 are now fitted (`long_range=("X",)`), with priors taken from the uncertainties of
the measurement they come from — Bacis, Cerny & Martin 1986 give De = 12 547.335(128) cm-1,
C6 = 1.48(12)e6, C8 = 3.86(1.20)e7, C10 = 1.0(0.5)e8, so C8 is known to 31% and C10 to 50%.

That does what it was meant to. **Collinearity between the high-v'' and low-v'' response directions
falls from +0.940 to +0.559**, with each long-range knob carrying about 0.17 of the high-v''
direction and exactly 0.000 of the low-v'' one. Fitting the four of them to the 38 high-v'' data
alone drops that residual from 384 GHz to 20.5 GHz and leaves every other data set numerically
unchanged, which is the orthogonality working as intended.

**But the magnitudes do not reach.** C6 moved by its own 1 sigma shifts v'' = 53 by only 0.6 cm-1,
and the discrepancy is 14 cm-1. The unconstrained prefit gets there only by going to
C6 = -2.2e7 (negative, so unphysical) at 192 sigma, C8 at 428 sigma and C10 at 550 sigma. Within
anything Bacis 1986 would recognise, the tail supplies perhaps 1-2 cm-1.

**R_O is the one handle that can.** It is a join radius chosen by the modeller, not a measured
quantity, and moving it acts only where the high-v'' levels live:

| R_O (A) | dE(v''=0) | dE(v''=17) | dE(v''=48) | dE(v''=53) | dE(v''=54) |
|---|---|---|---|---|---|
| 3.300 | 0 | 0 | 0 | 0 | 0 |
| 3.310 | -0.0000 | -0.0000 | +8.30 | +13.87 | +15.06 |

A shift of 0.010 A — tiny, and free of any measurement to contradict it — supplies the +14.7 cm-1
that v'' = 53 needs, with no effect whatever below v'' = 24.

## The obstacle is an oscillating extrapolation error

Mapping the published X potential against every absolute-frequency datum — model term value minus
the value the measurement implies, using the well-determined B levels — gives:

| v'' | 0-17 | 48 | 53 | 54 |
|---|---|---|---|---|
| error (cm-1) | -0.016 to +0.003 | **+7.76** | **-14.64** | **-21.12** |

Zero through the fitted range, then positive at v'' = 48, **through zero**, and strongly negative by
v'' = 53. The two high-v'' data sets are not in conflict: both are consistent with one
non-monotonic error, which is what a power series does when extrapolated past the data it was fitted
to. v'' = 18-47 is empty, so the turning point is unobserved.

This is why no single smooth handle repairs it. R_O, De and the dispersion coefficients all produce
**monotonic** corrections across v''; an oscillation needs either many more series terms — which the
prior and the shape constraint rightly resist — or a functional form whose long-range behaviour is
correct by construction. That is the Phase C case, and it is now evidence rather than preference.

## Superseded framing: a contradiction in the data

`nesterenko2019` needs E_X(53) and E_X(54) raised by 14.7 and 21.1 cm-1. `matyugin2012` needs
E_X(48) **lowered** by 7.8. Every smooth handle we have — R_O, the dispersion terms, the series —
moves all three the same way. No single-parameter change can raise 53 and lower 48; their
wavefunctions overlap too much.

**The assignment has since been confirmed from the authors' own companion paper.** Matyugin et al.,
Quantum Electron. 38, 755 (2008) states that the emission studied is J' = 57, v' = 32 -> J'' = 58,
v'' = 48, pumped at J'' = 56, v'' = 0 -> J' = 57, v' = 32. That fixes the convention (the labels are
emission lines in standard notation, and pairs share an upper level) and matches the residuals of the
18 rows kept: +9.80 and +9.93 cm-1 for the J' = 57 pair, +7.51 and +7.31 for the J' = 87 pair.

So the first explanation is now much weaker, and the contradiction looks real:

**The X-representation form appears unable to span v'' = 0 to 54 with one smooth outer branch.** That
is the Phase C argument, and it now rests on confirmed assignments rather than a reconstruction. An
MLR potential with a theoretically anchored long-range tail exists precisely for this case.

**This has now been tested, and it works**: see `docs/design/mlr-x.md`. One MLR fits v'' <= 17 and the
measured v'' = 48, 53 and 54 simultaneously, with De and the dispersion coefficients at their published
values, so the contradiction was in the functional form rather than in the data.

Two caveats remain. Two of matyugin2012's twenty rows (R87 and P85) are still unexplained under any
branch choice or transposition, so something in that table is not fully understood. And the argument
is linear, around the published parameters; a nonlinear path that satisfies both data sets has not
been excluded, only shown to be unreachable by any single smooth handle.

## Two ways the fit went wrong first, both worth recording

**A soft prior on the calibration offsets is no constraint against a large physics error.** Given a
0.3 MHz Gaussian prior and residuals of 5x10^5 sigma, the optimiser found it cheaper to absorb the
14 cm-1 discrepancy into free per-group offsets than to bend the potential: the offsets ran to
**432 GHz** and every data set degraded by two to four orders of magnitude. They now carry a hard
+/- 3 MHz bound (`calibration_bound`), which is all a frequency-chain error could credibly be. With
the bound in place the fitted offsets are all under 0.3 MHz, as they should be.

**A very large first f_scale hands the fit to a handful of points.** Under `soft_l1` the gradient
from a residual is bounded by f_scale, so starting at 1e7 let eighteen nesterenko points with 2-235
kHz uncertainties outvote 872 data. The schedule now starts at 1e4.

## Known gaps

- Calibration nuisance parameters per `group` are not implemented yet. The groups are carried on
  every datum and Dubé's four lines already show why they are needed.
- Stage 2 (hyperfine, potentials frozen) is written and run, with a negative result: `docs/design/hyperfine-fit.md`.
  Stage 3 (BO corrections) is not written.
- The covariance returned is the linearised one at the optimum; M6 replaces it.
