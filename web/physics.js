// The physics the browser has to do. Everything else — the line positions, the hyperfine
// components, the temperature-independent strengths — is precomputed by `i2spec web export` (i2spec/webapp.py), so this file
// holds only what depends on the viewer's choices: the unit they typed in, the cell temperature,
// the cold finger. Kept separate from the UI so test_physics.mjs can check it against Python.

export const C2 = 1.438776877;          // hc/k, cm K
export const MHZ_PER_CM = 29979.2458;
const KB = 1.380649e-23, U = 1.66053906892e-27, C_LIGHT = 299792458, TORR = 133.322368;

export const MASS = {"127I2": 253.8089438, "129I2": 257.8099674, "127I129I": 255.8094556};  // u

// --- units ------------------------------------------------------------------------------------

/** Refractive index of standard air (Edlén 1966), λ in nm. */
export function airIndex(lamVac) {
  const s2 = (1e3 / lamVac) ** 2;
  return 1 + 8.34254e-5 + 2.406147e-2 / (130 - s2) + 1.5998e-4 / (38.9 - s2);
}
// The dispersion formula wants the vacuum wavelength, so converting the other way iterates once.
export const airToVac = lamAir => lamAir * airIndex(lamAir * airIndex(lamAir));
export const vacToAir = lamVac => lamVac / airIndex(lamVac);

export const UNITS = ["nm", "nm-air", "cm-1", "THz", "MHz"];

export function toWavenumber(v, unit) {
  switch (unit) {
    case "nm": return 1e7 / v;
    case "nm-air": return 1e7 / airToVac(v);
    case "cm-1": return v;
    case "THz": return v * 1e6 / MHZ_PER_CM;
    case "MHz": return v / MHZ_PER_CM;
  }
  throw new Error("unknown unit " + unit);
}

export function fromWavenumber(nu, unit) {
  switch (unit) {
    case "nm": return 1e7 / nu;
    case "nm-air": return vacToAir(1e7 / nu);
    case "cm-1": return nu;
    case "THz": return nu * MHZ_PER_CM / 1e6;
    case "MHz": return nu * MHZ_PER_CM;
  }
  throw new Error("unknown unit " + unit);
}

/** Digits worth showing in each unit, so a 1 MHz step is still visible. */
export const DIGITS = {"nm": 6, "nm-air": 6, "cm-1": 5, "THz": 7, "MHz": 3};

// --- level populations --------------------------------------------------------------------------

/**
 * Rovibrational partition function of the X state.
 *
 * `table` is what the export shipped: `g[J]` is the nuclear-spin weight of that J (it alternates
 * with the parity of J, and differs between the homonuclear and the mixed isotopologue) and
 * `levels[J]` lists the term values of every bound v at that J.
 */
export function partitionFunction(table, T) {
  let q = 0;
  for (let J = 0; J < table.g.length; J++) {
    const levels = table.levels[J];
    let inner = 0;
    for (let i = 0; i < levels.length; i++) inner += Math.exp(-C2 * levels[i] / T);
    q += table.g[J] * (2 * J + 1) * inner;
  }
  return q;
}

/**
 * Line strength in cm/molecule at temperature T.
 *
 * The export stores `s0 = S(300 K) · Q(300 K) · exp(+c2 E″/300)`, which is the part of the strength
 * that does not depend on temperature; undoing that here costs one exponential per line.
 */
export const strengthAt = (line, qT, T) => line.s0 * Math.exp(-C2 * line.el / T) / qT;

/** Doppler FWHM in cm⁻¹. */
export const dopplerFWHM = (nu, T, iso) =>
  nu * Math.sqrt(8 * KB * T * Math.LN2 / (MASS[iso] * U * C_LIGHT ** 2));

/** I₂ vapour pressure over the solid, in Pa (Tellinghuisen 2011 eq. 8). */
export const vapourPressure = T =>
  Math.pow(10, 12.1891 - 0.001301 * T - 0.3523 * Math.log10(T) - 3410.71 / T) * TORR;

/** I₂ molecules per cm³ in a cell whose pressure is set by its cold finger. */
export const numberDensity = (Tcold, Tcell) => vapourPressure(Tcold) / (KB * Tcell) * 1e-6;

// --- lasers and their harmonics -------------------------------------------------------------------

/**
 * The wavenumber of the n-th harmonic of a laser at `value` in `unit`: harmonic generation multiplies
 * the optical frequency by n. An air wavelength is converted to vacuum at the fundamental, where the
 * laser is, so 1064 nm in air doubles to a vacuum wavenumber, not to 532 nm in air.
 */
export function harmonicWavenumber(value, unit, n) {
  if (!(Number.isInteger(n) && n >= 1)) throw new Error("harmonic order must be a positive integer, not " + n);
  return n * toWavenumber(value, unit);
}

/** A value of the n-th harmonic, `nu` in cm⁻¹, expressed at the fundamental in `unit`. */
export const fundamentalOf = (nu, n, unit = "cm-1") => fromWavenumber(nu / n, unit);

/** The view [lo, hi] in cm⁻¹ that spans ± halfNm of vacuum wavelength around the n-th harmonic. */
export function harmonicView(value, unit, n, halfNm = 0.1) {
  const lam = 1e7 / harmonicWavenumber(value, unit, n);
  return [1e7 / (lam + halfNm), 1e7 / (lam - halfNm)];
}

// --- sub-Doppler (saturated-absorption) spectra ------------------------------------------------------
// The same model as i2spec/saturation.py, which derives it (docs/research/sub-doppler.md): a Lamb dip at
// every hyperfine component and a crossover halfway between two components that share a level, with
// weak-saturation amplitudes, each drawn as a Lorentzian of the homogeneous width or as its n-th
// harmonic under sinusoidal frequency modulation.

/** Homogeneous width (MHz) drawn by default: the terminal browser's (tui.SUB_DOPPLER_FWHM). */
export const SUB_DOPPLER_FWHM = 2.0;
/** Crossovers are dropped beyond this multiple of the Doppler FWHM (saturation.REACH). */
export const REACH = 3.0;

/**
 * Lamb dips and crossovers of one line, strongest first, amplitudes summing to 1 (saturation.resonances).
 *
 * `components`: [{offset (MHz), strength, lower, upper}], where `lower` and `upper` identify the levels
 * (saturation.py uses (I, F)); two components share a level when those compare equal. Components without
 * them produce Lamb dips only. `dopplerWidth` is the Doppler FWHM in MHz. Returns {list, total}, `total`
 * being the sum the amplitudes were divided by, so lines of different strength can be combined.
 */
export function resonances(components, dopplerWidth, {exponent = 2.0, lambdaWeight = 1.0, reach = REACH, threshold = 1e-5} = {}) {
  const found = components.map((c, i) => ({offset: c.offset, amplitude: c.strength ** exponent, kind: "dip",
                                           sharing: null, components: [i, i]}));
  const known = x => x !== undefined && x !== null;
  for (let i = 0; i < components.length; i++) {
    const a = components[i];
    for (let j = i + 1; j < components.length; j++) {
      const b = components[j];
      const lower = known(a.lower) && a.lower === b.lower, upper = known(a.upper) && a.upper === b.upper;
      const sharing = lower && upper ? null : lower ? "lower" : upper ? "upper" : null;
      const gap = b.offset - a.offset;
      if (sharing === null || Math.abs(gap) > reach * dopplerWidth) continue;
      const weight = sharing === "lower" ? 1.0 : lambdaWeight;
      const amplitude = 2 * weight * (a.strength * b.strength) ** (exponent / 2) *
        Math.exp(-Math.LN2 * (gap / dopplerWidth) ** 2);
      found.push({offset: (a.offset + b.offset) / 2, amplitude, kind: "crossover", sharing, components: [i, j]});
    }
  }
  const strongest = found.reduce((m, r) => Math.max(m, r.amplitude), 0);
  const kept = found.filter(r => r.amplitude > threshold * strongest);
  const total = kept.reduce((s, r) => s + r.amplitude, 0);
  const list = kept.map(r => ({...r, amplitude: r.amplitude / total})).sort((p, q) => q.amplitude - p.amplitude);
  return {list, total};
}

const lorentzian = (x, fwhm) => 1 / (1 + (2 * x / fwhm) ** 2);

/** The midpoint-rule weights and cosines of lineshape(), which depend only on the harmonic and the points. */
const QUADRATURE = new Map();
function quadrature(harmonic, points) {
  const key = harmonic + "|" + points;
  if (!QUADRATURE.has(key)) {
    const cos = new Float64Array(points), weight = new Float64Array(points);
    for (let k = 0; k < points; k++) {
      const theta = (k + 0.5) * Math.PI / points;
      cos[k] = Math.cos(theta); weight[k] = Math.cos(harmonic * theta) * (2 / points);
    }
    QUADRATURE.set(key, {cos, weight});
  }
  return QUADRATURE.get(key);
}

/**
 * Saturation line shape at detuning x (MHz), peak 1 unmodulated (saturation.lineshape): a Lorentzian, or
 * with modulation m (MHz peak to peak) its in-phase n-th harmonic,
 * S_n(x) = (2/π) ∫₀^π L(x + (m/2) cos θ) cos(nθ) dθ, by the midpoint rule.
 */
export function lineshape(x, fwhm, {harmonic = 0, modulation = 0.0, points = 256} = {}) {
  if (harmonic === 0 && modulation === 0) return lorentzian(x, fwhm);
  const {cos, weight} = quadrature(harmonic, points);
  let s = 0;
  for (let k = 0; k < points; k++) s += weight[k] * lorentzian(x + 0.5 * modulation * cos[k], fwhm);
  return s;
}

/** The modulated shape tabulated as saturation.signal() does: step fwhm/256 out to `wings` half-widths. */
function kernel(fwhm, harmonic, modulation, points, wings, step) {
  const half = wings * (fwhm + modulation);
  const start = -half, n = Math.ceil((half + step - start) / step);    // numpy.arange's length
  const values = new Float64Array(n);
  for (let k = 0; k < n; k++) values[k] = lineshape(start + k * step, fwhm, {harmonic, modulation, points});
  return {start, step, n, values};
}
/** numpy.interp on the kernel's grid, zero outside it. */
function interpKernel(K, x) {
  const f = (x - K.start) / K.step;
  if (f < 0 || f > K.n - 1) return 0;
  const k = Math.min(Math.floor(f), K.n - 2), t = f - k;
  return K.values[k] + t * (K.values[k + 1] - K.values[k]);
}

/**
 * The Doppler-free signal at each of `nu` (MHz, same origin as the resonance offsets), as
 * saturation.signal(): Σ amplitude × lineshape(ν − offset), positive at each resonance unmodulated.
 * Modulated shapes are interpolated from a tabulated kernel exactly as the Python does.
 */
export function signal(nu, list, fwhm, {harmonic = 0, modulation = 0.0, points = 256, wings = 100.0} = {}) {
  const out = new Float64Array(nu.length);
  if (!list.length) return out;
  if (harmonic === 0 && modulation === 0) {
    for (let k = 0; k < nu.length; k++) {
      let s = 0;
      for (const r of list) s += r.amplitude * lorentzian(nu[k] - r.offset, fwhm);
      out[k] = s;
    }
    return out;
  }
  const K = kernel(fwhm, harmonic, modulation, points, wings, fwhm / 256);
  for (let k = 0; k < nu.length; k++) {
    let s = 0;
    for (const r of list) s += r.amplitude * interpKernel(K, nu[k] - r.offset);
    out[k] = s;
  }
  return out;
}

let lastKernel = {key: null, K: null};      // a redraw at the same settings reuses the table
/**
 * The same signal added into `bins`, sampled at x0 + (k + ½)·dx (MHz), for a plot. Each resonance only
 * touches the bins within `wings` half-widths of it, and a modulated shape is tabulated `perWidth` times
 * per FWHM, so a view with a thousand resonances still draws in a frame. `scale` multiplies every amplitude.
 */
export function addSignal(bins, x0, dx, list, fwhm, {harmonic = 0, modulation = 0.0, points = 64, wings = 40, perWidth = 128, scale = 1} = {}) {
  const reach = wings * (fwhm + modulation), n = bins.length;
  const plain = harmonic === 0 && modulation === 0;
  const key = [fwhm, harmonic, modulation, points, wings, perWidth].join("|");
  if (!plain && lastKernel.key !== key) lastKernel = {key, K: kernel(fwhm, harmonic, modulation, points, wings, fwhm / perWidth)};
  const K = plain ? null : lastKernel.K;
  for (const r of list) {
    const k0 = Math.max(0, Math.floor((r.offset - reach - x0) / dx)), k1 = Math.min(n - 1, Math.ceil((r.offset + reach - x0) / dx));
    const a = r.amplitude * scale;
    for (let k = k0; k <= k1; k++) {
      const x = x0 + (k + 0.5) * dx - r.offset;
      bins[k] += a * (plain ? lorentzian(x, fwhm) : interpKernel(K, x));
    }
  }
  return bins;
}
