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
