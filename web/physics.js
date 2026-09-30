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
