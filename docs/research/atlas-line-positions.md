# Line positions from the FTS atlases, and why the Nölleke list stays unassigned

*2026-09-20. `prototypes/atlas_lines.py`, `prototypes/atlas_dataset.py`, `prototypes/nolleke_assign.py`.
The product is `data/atlas_lines/`: 67 729 measured ¹²⁷I₂ B–X line centres.*

## Why

Hannover fitted their potentials to ~1950 lines, of which perhaps two thirds were atlas lines — Kato's
Doppler-free atlas and recalibrated Gerstenkorn–Luc, at 0.002–0.005 cm⁻¹ each. Those are what tie a
global potential down band by band; our 1013 precision rows are a thousand times more accurate but
touch only a fraction of the (v′, v″) plane. To fit our own potentials we need the same kind of broad
coverage, and the two Doppler-limited FTS atlases we hold supply it.

## Method

For each 2.4 cm⁻¹ window, stepped by 2 cm⁻¹: build the model transmission (master line list at the
cell temperature, hyperfine structure for every line that matters, the bound–free continuum), fit the
window as `compare_salami_ross.fit` does — column density, linear baseline, zero offset, instrument
width — then give **every line deep enough to matter its own wavenumber shift** and fit all of them
jointly by linearised least squares with a 0.01 cm⁻¹ ridge prior. The Jacobian is analytic (one
convolution of `N σ_k′(ν) exp(−Nσ)` per line), so a window takes about 3 s. A blend comes out of the
covariance with a large σ instead of a confident wrong number. Repeats from overlapping windows and
overlapping temperature segments are combined by inverse variance.

The measurement this makes is the **hyperfine-free line centre on the atlas's own wavenumber scale**.
The scale itself — each atlas's calibration — stays out of the per-line uncertainty and belongs in the
fit as a group offset.

## Result

| atlas | lines | range | coverage | median σ |
|---|---|---|---|---|
| `salami_ross_2005` | 27 601 | 14 320–19 509 cm⁻¹ | v″ = 0–8, v′ = 2–53 | 19 MHz |
| `apo_nist_2009` | 40 128 | 15 405–19 466 cm⁻¹ | v″ = 0–5, v′ = 5–58 | 14 MHz |

The APO scan is the better of the two, as its provenance suggests (a NIST 2-m FTS scan at 0.018 cm⁻¹
against Salami & Ross's 0.02 cm⁻¹ resolution and coarser sampling), and its cell temperature, 333 K,
was fitted from the window residuals.

**Uncertainties.** The window fit's covariance is optimistic: it knows the noise but not what the
template gets wrong on a blend or a hyperfine pattern. The two atlases are independent measurements of
the same lines, so their difference calibrates the rest: binned by fit σ it gives
σ_eff = √((1.55 σ_fit)² + 7.5²) MHz. Checked against every precision measurement of the same lines,
this is right — normalised residuals 0.67–1.32:

| set | Salami–Ross bias / rms / z | APO bias / rms / z |
|---|---|---|
| `sansonetti1997a` (98, 92 lines) | +8.9 / 18.7 / 1.17 | −18.1 / 10.1 / 1.04 |
| `velchev1998a` (108, 113) | +7.1 / 18.2 / 1.04 | −17.6 / 13.7 / 0.89 |
| `xu2000a` (444, 402) | +5.4 / 24.6 / 1.32 | −18.2 / 11.9 / 0.67 |
| `kobayashi2016a`, `nishiyama2024a`, `yoshiki2023a` | −2 to +10 / 14–19 | −18 to −22 / 6–11 |

So each atlas carries a constant calibration offset — Salami–Ross about +7 MHz, APO about −18 MHz,
their difference +30 MHz — and, once that is removed, gives line centres good to **10–25 MHz**. That is
an order of magnitude inside the ±0.003 cm⁻¹ (±90 MHz) the atlases claim for themselves, because we are
not reading a scale off a chart: we are fitting a known template and letting the calibration float.

*One bug worth recording: the first run wrote its output with `%.6g`, which rounds a wavenumber to
0.1 cm⁻¹. Everything looked sane except a suspiciously uniform ±1025 MHz spread against the model
whatever the cut — the signature of a quantisation, not of physics.*

## The Nölleke 2018 list stays unassigned

Nölleke et al. (2018) published 10 162 absorption lines at 915–985 nm (10 152–10 929 cm⁻¹) with 50 MHz
absolute accuracy and **no quantum numbers**. Those lines come from v″ ≈ 24–45 — the gap nothing else
reaches, and exactly where the published X potential and our MLR disagree. Assigning them would be a
real contribution and would anchor the X state through its unmeasured region.

It does not work with the model we have:

- **Band pattern search.** For each of the 41 candidate bands, slide its predicted lines over ±0.6 cm⁻¹
  and count intensity-weighted matches within 0.01 cm⁻¹: best z = 4.5, against ~20 000 trials — nothing.
- **Two-parameter search** (origin shift ±6 cm⁻¹ and a rotational correction ±4×10⁻⁴ cm⁻¹ in J(J+1)/10⁴)
  on the eight strongest bands: best z = 5.0, still consistent with chance.
- **Control.** The same matcher, on our own APO line centres at 17 000–17 800 cm⁻¹ treated as anonymous
  positions — the same line density, 12.1 per cm⁻¹ against Nölleke's 13.1 — finds the correct bands at
  **z = 21 with a shift of 0.000–0.004 cm⁻¹**, 45 bands above z = 10. The method works; the prediction
  is what fails.
- **Full spectrum cross-correlation** over 30 cm⁻¹ windows: the model predicts 781–866 features where
  they report 411–417, and no shift stands out (z ≈ 3).

Three things would have to be settled before another attempt: their detection threshold and cell
conditions (the model predicts twice as many features as they list, which alone buries the signal);
whether their reported positions are hyperfine-free centres or the peaks of blended hyperfine clusters,
which at 950 nm span 0.02–0.03 cm⁻¹ against a 0.069 cm⁻¹ mean line spacing; and the X levels at
v″ = 24–45, which is the thing the assignment was supposed to measure. The honest reading is that this
is its own project, and that it needs the Gerstenkorn 11 000–14 000 cm⁻¹ atlas (requested) far more
than it needs a better search.

## Second attempt (2026-09-26): the raw scans, at the right temperature

The first attempt modelled the cell at 300 K. Nölleke's cell body was at 300 °C (573 K; cold finger 39 °C,
127 Pa, 225 cm path), and absorption at 10 150–10 930 cm⁻¹ needs X levels at least ~5 000 cm⁻¹ up, so the
temperature sets the band shapes the search relied on. `prototypes/nolleke_transmission.py` models the
transmission for their cell and cross-correlates it with the **raw scans** (`scan_data.csv`), which avoids
their peak-picking threshold and blend splitting, per 10 cm⁻¹ window over ±0.3 cm⁻¹ of shift.

- **The scan contains the 940 nm water band.** Its strongest features (absorbance 0.10–0.18 at
  10 450–10 750 cm⁻¹) are atmospheric H₂O in the open beam; the list flags 192 of them. They are masked
  (±0.05 cm⁻¹, and anything above absorbance 0.03), 2.6 % of the scan.
- **Neither model registers:** i2spec2026i 0 of 78 windows at z > 5, hannover2008 2 of 78 at shifts of
  ~7 GHz (chance) — including the windows near 10 900 cm⁻¹ whose strongest B–X lines come from
  v″ = 24–26, which the Orsay atlas has just measured to ~50 MHz.
- **The intensities do not fit B–X.** At their stated conditions B–X lines from v″ ≈ 25–30 (Boltzmann
  factor ~e⁻¹³) peak near 0.1 % absorption; the listed iodine lines are mostly 1–3 %. A factor of 10–30
  would need a column density or temperature far from what they report.

The likely reading is that these lines are **not B–X hot bands** but another system from low, well-
populated v″ — the A ³Π₁ᵤ–X and A′ ³Π₂ᵤ–X systems lie at 10 000–11 000 cm⁻¹ — which our model does not
contain. A ground-state combination-difference census (`prototypes/nolleke_combdiff.py`), which needs no
upper-state model, is too weak to decide (10 000 lines × ~200 J give ~10⁵ chance pairs against a few
hundred real ones per band). This remains an inference, not a demonstration.

Either way, the conclusion for this project changes: lines from low v″ would not constrain X at
v″ = 26–30, which is what the list was wanted for. It is set aside; an A-state model would be its own
project.
