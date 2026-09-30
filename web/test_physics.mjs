// Checks physics.js against the Python model. The reference values are regenerated with
//   uv run python web/make_reference.py > web/reference.json
// and the exported data in web/data/ is used directly, so this also catches an export that no
// longer matches what the app expects.
//
//   node web/test_physics.mjs

import {readFileSync} from "node:fs";
import {dirname, join} from "node:path";
import {fileURLToPath} from "node:url";

import * as P from "./physics.js";

const HERE = dirname(fileURLToPath(import.meta.url));
const read = p => JSON.parse(readFileSync(join(HERE, p), "utf8"));
const ref = read("reference.json");
const manifest = read("data/manifest.json");

let failures = 0;
function close(name, got, want, rtol) {
  const err = Math.abs(got - want) / Math.abs(want || 1);
  const ok = err <= rtol;
  if (!ok) failures++;
  console.log(`${ok ? "  ok  " : "  FAIL"} ${name.padEnd(46)} ${got.toPrecision(10)} vs ${want.toPrecision(10)}  (${err.toExponential(1)})`);
}
function check(name, ok, detail = "") {
  if (!ok) failures++;
  console.log(`${ok ? "  ok  " : "  FAIL"} ${name}${detail ? "  " + detail : ""}`);
}

console.log("air refractive index and the air/vacuum conversions");
for (const [lam, n, air, vac] of ref.air) {
  close(`n(${lam} nm)`, P.airIndex(lam), n, 1e-12);
  close(`vac→air(${lam})`, P.vacToAir(lam), air, 1e-12);
  close(`air→vac(${lam})`, P.airToVac(lam), vac, 1e-9);
}

console.log("\nunit round trips");
for (const unit of P.UNITS) {
  for (const nu of [12600.0, 15800.0, 18788.4372]) {
    close(`${unit} round trip at ${nu} cm⁻¹`, P.toWavenumber(P.fromWavenumber(nu, unit), unit), nu, 1e-12);
  }
}

console.log("\nvapour pressure and number density");
for (const [T, p, n] of ref.vapor) {
  close(`p(${T} K)`, P.vapourPressure(T), p, 1e-12);
  close(`N(${T} K cold, 293.15 K cell)`, P.numberDensity(T, 293.15), n, 1e-12);
}

console.log("\nDoppler width");
for (const [nu, T, w] of ref.doppler) close(`FWHM(${nu} cm⁻¹, ${T} K)`, P.dopplerFWHM(nu, T, "127I2"), w, 1e-12);

console.log("\npartition function, from the exported table");
const table = read("data/partition_127I2.json");
for (const [T, q] of ref.Q) close(`Q(${T} K)`, P.partitionFunction(table, T), q, 1e-6);

console.log("\nline strengths, from the exported shards");
// find the exported rows for the reference lines and rebuild S(T) the way the browser does
const shards = manifest.isotopologues["127I2"].shards;
for (const want of ref.strengths) {
  const shard = shards.find(s => want.nu >= s.nu0 && want.nu < s.nu1);
  const d = read("data/" + shard.file);
  let k = -1, best = Infinity;
  for (let i = 0; i < d.nu.length; i++) {
    const e = Math.abs(d.nu[i] - want.nu);
    if (e < best) { best = e; k = i; }
  }
  check(`line at ${want.nu.toFixed(5)} cm⁻¹ found in ${shard.file}`, best < 1e-4, `Δν = ${best.toExponential(1)} cm⁻¹`);
  const line = {s0: d.s0[k], el: d.el[k]};
  for (const [T, S] of Object.entries(want.S)) {
    // 1e-3 because the export rounds s0 to four significant digits and E″ to 0.01 cm⁻¹
    close(`  S(${T} K)`, P.strengthAt(line, P.partitionFunction(table, Number(T)), Number(T)), S, 2e-3);
  }
}

console.log("\nper-line uncertainty, as the Python model computes it");
for (const want of ref.uncertainty) {
  const shard = shards.find(s => want.nu >= s.nu0 && want.nu < s.nu1);
  const d = read("data/" + shard.file);
  let k = -1, best = Infinity;
  for (let i = 0; i < d.nu.length; i++) {
    const e = Math.abs(d.nu[i] - want.nu);
    if (e < best) { best = e; k = i; }
  }
  close(`${want.label}: u`, d.u[k], want.u, 5e-3);           // exported to three significant figures
  const reason = manifest.reasons[d.w[k]];
  check(`${want.label}: reason and flags`, reason && reason.why.length > 40 &&
    JSON.stringify(reason.flags.map(f => f[0])) === JSON.stringify(want.flags), JSON.stringify(reason && reason.flags));
  check(`${want.label}: marked measured`, d.m.includes(k));
}

console.log("\nlaser harmonics");
// 1064 nm doubled is 532 nm exactly in vacuum; 281.63 THz doubled is 563.26 THz
close("2 × 1064 nm (vacuum)", P.fromWavenumber(P.harmonicWavenumber(1064, "nm", 2), "nm"), 532, 1e-15);
close("2 × 281.63 THz", P.fromWavenumber(P.harmonicWavenumber(281.63, "THz", 2), "THz"), 563.26, 1e-15);
close("3 × 1542 nm", P.harmonicWavenumber(1542, "nm", 3), 3e7 / 1542, 1e-15);
close("4 × 9394 cm⁻¹", P.harmonicWavenumber(9394, "cm-1", 4), 37576, 1e-15);
close("1 × 1319 nm", P.harmonicWavenumber(1319, "nm", 1), 1e7 / 1319, 1e-15);
close("2 × 3e8 MHz", P.harmonicWavenumber(3e8, "MHz", 2), 6e8 / P.MHZ_PER_CM, 1e-15);
// an air wavelength doubles through its frequency: the harmonic of 1064 nm air is not 532 nm air
const airH = P.fromWavenumber(P.harmonicWavenumber(1064, "nm-air", 2), "nm-air");
close("2 × 1064 nm air = vac(1064 air)/2 in vacuum", P.fromWavenumber(P.harmonicWavenumber(1064, "nm-air", 2), "nm"),
  P.airToVac(1064) / 2, 1e-12);
check("2 × 1064 nm air differs from 532 nm air by the dispersion", Math.abs(airH - 532) > 1e-3, `${airH.toFixed(6)} nm`);
for (const n of [1, 2, 3, 4]) {
  const nu = P.harmonicWavenumber(1111.6, "nm", n);
  close(`fundamental of ${n} × 1111.6 nm`, P.fundamentalOf(nu, n, "nm"), 1111.6, 1e-14);
}
const [vlo, vhi] = P.harmonicView(1064.49, "nm", 2, 0.1);
close("view around 2 × 1064.49 nm: short edge", 1e7 / vhi, 532.145, 1e-13);
close("view around 2 × 1064.49 nm: long edge", 1e7 / vlo, 532.345, 1e-13);
let threw = false;
try { P.harmonicWavenumber(1064, "nm", 0); } catch { threw = true; }
check("harmonic order 0 is refused", threw);

console.log("\nsub-Doppler spectra, against i2spec.saturation");
const sd = ref.sub_doppler;
const comps = sd.components;
for (const {options, list: want} of sd.resonances) {
  const opts = {exponent: options.exponent, lambdaWeight: options.lambda_weight, threshold: options.threshold};
  for (const k of Object.keys(opts)) if (opts[k] === undefined) delete opts[k];
  const {list: got} = P.resonances(comps, sd.doppler_width, opts);
  const name = `${sd.line} ${JSON.stringify(options)}`;
  check(`${name}: ${want.length} resonances (${want.filter(r => r.kind === "crossover").length} crossovers)`,
    got.length === want.length, `got ${got.length}`);
  // match by the components each one comes from; the order of equal amplitudes is not significant
  const key = r => r.components.join(",");
  const byKey = new Map(got.map(r => [key(r), r]));
  let worst = 0, same = true;
  for (const w of want) {
    const g = byKey.get(key(w));
    if (!g || g.kind !== w.kind || g.sharing !== w.sharing || Math.abs(g.offset - w.offset) > 1e-9) { same = false; continue; }
    worst = Math.max(worst, Math.abs(g.amplitude - w.amplitude) / w.amplitude);
  }
  check(`${name}: same offsets, kinds and shared levels`, same);
  check(`${name}: amplitudes to 1e-9 relative`, worst < 1e-9, `worst ${worst.toExponential(1)}`);
}
const dipsOnly = P.resonances(comps.map(({offset, strength}) => ({offset, strength})), sd.doppler_width).list;
check("without level labels: Lamb dips only", dipsOnly.every(r => r.kind === "dip") && dipsOnly.length > 0);
for (const {fwhm, harmonic, modulation, x, y} of sd.lineshape) {
  const peak = Math.max(...y.map(Math.abs));
  let worst = 0;
  x.forEach((xi, k) => { worst = Math.max(worst, Math.abs(P.lineshape(xi, fwhm, {harmonic, modulation}) - y[k]) / peak); });
  check(`lineshape Γ = ${fwhm}, n = ${harmonic}, m = ${modulation}: to 1e-12 of its peak`, worst < 1e-12, `worst ${worst.toExponential(1)}`);
}
const found = P.resonances(comps, sd.doppler_width).list;
for (const {fwhm, harmonic, modulation, nu, y} of sd.signal) {
  const got = P.signal(nu, found, fwhm, {harmonic, modulation});
  const peak = Math.max(...y.map(Math.abs));
  let worst = 0, worstRel = 0;
  y.forEach((v, k) => {
    worst = Math.max(worst, Math.abs(got[k] - v) / peak);
    if (Math.abs(v) > 1e-3 * peak) worstRel = Math.max(worstRel, Math.abs(got[k] - v) / Math.abs(v));
  });
  check(`signal Γ = ${fwhm}, n = ${harmonic}: to 1e-6 relative (above 1e-3 of the peak)`, worstRel < 1e-6,
    `worst ${worstRel.toExponential(1)}; ${worst.toExponential(1)} of the peak anywhere`);
  // the plot's version: the same curve on a pixel grid, from a coarser kernel
  const x0 = -600, dx = 1200 / 2000, bins = P.addSignal(new Float64Array(2000), x0, dx, found, fwhm, {harmonic, modulation});
  const at = Array.from(bins, (_, k) => x0 + (k + 0.5) * dx), exact = P.signal(at, found, fwhm, {harmonic, modulation});
  const dev = Math.max(...Array.from(bins, (b, k) => Math.abs(b - exact[k]))) / Math.max(...Array.from(exact, Math.abs));
  check(`plot signal Γ = ${fwhm}, n = ${harmonic}: within 1e-3 of the peak`, dev < 1e-3, `${dev.toExponential(1)}`);
}

console.log("\nexported data is self-consistent");
for (const [iso, info] of Object.entries(manifest.isotopologues)) {
  const total = info.shards.reduce((a, s) => a + s.n, 0);
  const d0 = read("data/" + info.shards[0].file);
  check(`${iso}: rows carry u, w and their measurements`, d0.u.length === d0.nu.length && d0.w.length === d0.nu.length &&
    d0.ms.length === d0.m.length && d0.w.every(w => w >= 0 && w < manifest.reasons.length));
  check(`${iso}: shard counts sum to n_lines`, total === info.n_lines, `${total} = ${info.n_lines}`);
  // stored in increasing wavelength: shards, and the rows within each, run to lower wavenumber
  check(`${iso}: shards are contiguous and in wavelength order`,
    info.shards.every((s, i) => s.nu0 < s.nu1 && (i === 0 || s.nu1 <= info.shards[i - 1].nu0 + 1e-6)));
  check(`${iso}: every shard has a ΔJ = 0 hyperfine file with one pattern per row`,
    info.shards.slice(0, 6).every(s => {
      const d = read("data/" + s.file), h = read("data/" + s.hfs);
      return h.o.length === d.nu.length && h.s.length === d.nu.length &&
        h.o.every((o, i) => o.length === h.s[i].length && o.every((v, k) => k === 0 || v >= o[k - 1]));
    }));
  check(`${iso}: rows within a shard are in wavelength order`,
    info.shards.slice(0, 8).every(s => { const d = read("data/" + s.file); return d.nu.every((v, i) => i === 0 || v <= d.nu[i - 1]); }));
  const hfs = read("data/" + info.hfs);
  const bad = Object.entries(hfs).filter(([, c]) => c.o.length !== c.s.length || c.o.length !== c.l.length);
  check(`${iso}: ${Object.keys(hfs).length} hyperfine patterns have matching arrays`, bad.length === 0);
  const named = Object.values(hfs).filter(c => c.l.some(Boolean)).length;
  check(`${iso}: every pattern has labelled main components`, named === Object.keys(hfs).length);
}

console.log(failures ? `\n${failures} FAILED` : "\nall checks passed");
process.exit(failures ? 1 : 0);
