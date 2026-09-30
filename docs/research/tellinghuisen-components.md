# The model's A←X, B←X and C←X against Tellinghuisen's component bands

*2026-09-30. Script: `prototypes/tellinghuisen_components.py`. Output: `prototypes/out/tellinghuisen_components.json` and `docs/figures/tellinghuisen_components.png`. Reference: Table IIS of the [T11C] supplement (J. Tellinghuisen, J. Chem. Phys. 135, 054301 (2011)), in `data/external/tellinghuisen_2011_supplement/` (see `data/external/README.md`). Model: `i2spec2026m`, `continuum_model` plus the line list of `paper/figures/fig9_xsec.py`.*

## What was compared

* **Units.** Table IIS gives decadic molar absorptivity ε (l mol⁻¹ cm⁻¹; A = ε l C, [T11B] eq. 2). σ = ε · 1000 ln 10 / N_A = 3.8235 × 10⁻²¹ ε cm². Below, ratios and differences are given in ε. To get σ, multiply by 3.82 × 10⁻²¹. "z" is (model − Tellinghuisen) divided by his LS standard error at 35 °C.
* **Wavelengths.** Read as air and converted to vacuum. The supplement doesn't say which scale it uses. Every earlier check of the model against [T11C] Table II favoured air (`continuum-model.md` §10), and so does every component here: 400–495 nm means of 0.996 (air) against 0.990 (vacuum) for C, 0.990 against 0.980 for B, and 1.004 against 1.008 for A at 730–850 nm. The means in the JSON are given both ways.
* **Temperatures.** 273.15 and 308.15 K.
* **Components.**
  * A and C: the model's A←X and C←X continua. They use the same Table I potentials and moments as [T11C].
  * B at λ ≤ 495 nm: the model's B←X bound–free. It uses Hannover B with the [T11C] inner wall and μ_B from [T11B] eq. (10).
  * B at λ ≥ 500 nm: [T11C]'s B is a pseudocontinuum, the discrete bands smeared into one smooth band. So the model's B there is the discrete line list plus the B continuum above the line-list cut, averaged over a 4 nm FWHM Gaussian. The 2 nm and 8 nm averages are in the JSON. The 2 nm average leaves band structure of ±10–20 % at 580–640 nm; the 8 nm average is monotonic and within 1 % of the 4 nm one.

## Results (35 °C; 0 °C is the same to within 0.3 % unless stated)

| Component | Range | model / T11C | Mean difference (ε) | z |
|---|---|---|---|---|
| A←X | 530–650 nm | 0.992–0.999 | −0.04 | ≤ 0.3 |
| A←X | 655–725 nm | 0.999–1.002 | +0.01 | ≤ 0.6 |
| A←X | 730–850 nm | 1.002–1.009 (mean 1.004) | +0.03 | 0.2–0.6 |
| C←X | 400–495 nm | 0.991–1.000 (mean 0.996) | −0.11 | ≤ 0.4 |
| C←X | 500–650 nm | 1.000–1.007 (0 °C: up to 1.015) | +0.08 | ≤ 0.15 |
| B←X bound–free | 400–495 nm | 0.984–0.998 (mean 0.990; 0.986 at 430–445 nm) | −0.67 | ≤ 0.5 |
| B←X banded (lines + continuum) | 500–520 nm | 0.989–1.001 | up to −5.9 at 510 nm | −2.5 to +0.3 |
| B←X banded | 525–560 nm | 1.004–1.018 | +3 to +8 | +1.5 to +7.9 |
| B←X banded | 565–625 nm | 1.013–1.047 (0 °C: up to 1.054) | +1 to +8 | +3 to +11 |
| B←X banded | 630–650 nm | 1.010–1.041 | +0.2 to +0.5 | +1 to +4 |
| B←X hot bands | 655–705 nm | 1.05–1.13 | +0.03 to +0.3 | +3 to +7 |
| Total | 400–495 nm | 0.993–0.999 (mean 0.995) | −0.61 | −1.8 to 0 |
| Total | 500–650 nm | 0.992–1.030 (mean 1.012) | +3.1 | −5.6 to +9.6 |
| Total | 655–850 nm | 1.002–1.009 | ≤ +0.1 | ≤ 2.7 |

**Measured ε (Table IS).** These are Tellinghuisen's spectra at 35.4 and 64.0 °C. The model total is lines plus continuum averaged over 2 nm, in the weak-absorption limit.

| Range | Model / measured, 35.4 °C | Model / measured, 64.0 °C |
|---|---|---|
| 420–500 nm | 0.987–0.997 | 0.991–1.002 |
| 600–625 nm | 1.040 (+3.6 ε) | 1.026 (+2.7 ε) |
| 625–775 nm | 1.000–1.009 | 0.993–1.013 |

Tellinghuisen's own total matches the same data to within 1 % at 460–630 nm.

## What it means

1. **A←X and C←X are Tellinghuisen's bands, reproduced.** They agree to ≤ 0.9 % everywhere they exceed 1 ε, and to within 0.6 of his standard errors. This is a check of the implementation, not of the physics: the model uses his Table I potentials and moments unchanged. The A curve is 0.4–0.9 % high at 730–850 nm. That is where the Morse stand-in for the A wall beyond 2.806 Å takes over (`continuum-model.md` §8, item 2). The size is within z ≤ 0.6.
2. **B←X bound–free is 1.0 % low at 400–495 nm, and 1.4 % low at 430–445 nm.** That is within z ≤ 0.5.
   * About 0.7–1.1 % of it comes from μ_B: [T11B] eq. (10), used by the model, has μ² 0.7–1.1 % below the Table I quadratic at R = 2.50–2.55 Å, and the two agree above 2.65 Å.
   * The rest comes from the B wall join (Hannover beyond 2.7154 Å, where [T11C] used the G&L RKR curve).
3. **The total's agreement at 400–495 nm doesn't hide compensating errors.** The model's B and C deviations have the same sign (−1.0 % and −0.4 %). The −0.5 % in the total is their sum, with B carrying about 90 % of it (at 480 nm, −1.4 of −1.6 ε).
   * Tellinghuisen's split is itself compensating. His B and C standard errors (about 3 ε each at 470–490 nm) are more than three times the total's (0.7–0.9). That implies a B–C correlation of about −0.95. This matches his statement that the B/C split depends on the B attachment point while the total does not.
   * The model inherits his C and his split, so this comparison can't test the split. It only shows that the model's B continuum matches his to 1 %.
4. **The discrete B←X in 525–650 nm is the one real disagreement.** The band-averaged line list rises steadily against his pseudocontinuum: +0.5 % at 525 nm, +2 % at 560 nm, +3–5 % at 585–625 nm, and +5–13 % in the hot-band tail at 655–705 nm (about 0.3 ε, under 1 % of the total there).
   * It isn't the transition moment: eq. (10) and the Table I μ_B differ by ≤ 1 % in μ² over 2.65–2.90 Å.
   * It isn't the averaging: the 8 nm average shows the same trend.
   * So it lies between the model's discrete bands (Hannover potentials) and his reflection-type pseudocontinuum on the G&L RKR curve. His table note warns that 500–650 nm is limited by that treatment.
   * The measured data point the same way where they exist. At 600–625 nm the model is 4.0 % (35 °C) and 2.6 % (64 °C) above Table IS, while his total is within 1 %. At 1–2 nm resolution, line-centre saturation in those spectra would also push the data low, so this is suggestive, not conclusive.
   * The 1.2 % mean excess of the total at 500–650 nm is entirely this B excess.

## Recommendation

* **Keep the A←X and C←X transition moments and potentials as they are.** Nothing in Table IIS asks for a change. Their ≤ 1 % differences are inside his errors, and the C split can't be tested against a table the model was built from.
* **Keep the B←X continuum.** The −1 % at 430–450 nm is z ≤ 0.5 and within [T11C]'s stated 0.5 % reliability plus its B/C model error. Switching the continuum to the Table I quadratic μ_B would remove most of it, but lines and continuum would then no longer share one μ_B(R). Don't do that.
* **Don't refit a continuum moment against Table IIS at 500–650 nm.** There it is a pseudocontinuum, not a measurement.
* **The item to follow up is the discrete B←X band strength at 560–630 nm (and its hot bands).** The model is 2–5 % above both Tellinghuisen's pseudocontinuum and his 600–625 nm spectra. The test should use resolved, unsaturated absolute data:
  * the Spietz 0.25 nm and 0.59 nm spectra and Saiz-Lopez (fig. 9, `cross-section-data.md`);
  * the [T11B] 520–640 nm spectra, if they can be obtained.
  * If an excess of 3–5 % at 580–630 nm is confirmed there, refit the slope of μ_B(R) at R ≈ 2.8–2.95 Å (the [T11B] eq. (10) parameters c₂, c₃) against those spectra. Do not touch the continuum.
