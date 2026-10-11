// The line explorer. All the model physics lives in physics.js; this file is the interface:
// a search over the exported line list, a spectrum, a results table and a line detail pane.

import * as P from "./physics.js";

const $ = s => document.querySelector(s);
const $$ = s => [...document.querySelectorAll(s)];
// toLocaleString builds a new formatter on every call, which is most of the cost of a table of numbers
const FORMATTERS = new Map();
const fmt = (v, d = 0) => {
  let f = FORMATTERS.get(d);
  if (!f) FORMATTERS.set(d, f = new Intl.NumberFormat("en-US", {minimumFractionDigits: d, maximumFractionDigits: d}));
  return f.format(v);
};
const fmtU = v => fmt(v, v < 1 ? 2 : v < 10 ? 1 : 0);        // an uncertainty, to two figures or so
const sentence = t => t ? t.charAt(0).toUpperCase() + t.slice(1) + (/[.)]$/.test(t) ? "" : ".") : "";
const DEFAULT_VIEW = {line: "P(53) 32-0", near: 18788.44, halfSpan: 2};   // cm⁻¹; resolve the exact centre from the export

const state = {
  manifest: null,
  iso: "127I2",
  T: 293.15,
  lo: 0, hi: 0,            // the view, in cm-1
  lines: [],               // lines inside the view
  selected: null,
  mode: "sub",
  sort: {key: "nu", dir: 1},  // by wavelength, shortest first
  shards: new Map(),       // file -> parsed shard
  partitions: new Map(),   // iso -> table
  qCache: new Map(),       // iso|T -> Q
  hfs: new Map(),          // iso -> patterns (ΔJ = ±2, for the lines the export selected)
  hfsShards: new Map(),    // hyperfine file -> ΔJ = 0 patterns of every row of one shard
  index: new Map(),        // iso -> label -> line, built only when someone searches by label
  refs: null,              // the measurement sources: index -> {short, citation, link, ...}
  kinds: new Set(["model", "atlas", "precision"]),   // which kinds of line the plot and table show
  show: "all",             // the same as a key: "all", or the kinds shown, e.g. "atlas,precision"
  laser: null,             // {value, unit, n, nu}: a fundamental laser and the harmonic of it iodine sees
  sub: {fwhm: P.SUB_DOPPLER_FWHM, harmonic: 0, modulation: 1.0},   // the sub-Doppler display
  selectedComps: null,     // the hyperfine components the detail pane shows, for the CSV export
  busy: 0,
};

const PRESETS = [
  ["532 nm", "Nd:YAG second harmonic", 18787.8, 18789.0],
  ["543 nm", "green HeNe", 18409.0, 18410.2],
  ["612 nm", "HeNe", 16338.0, 16339.2],
  ["633 nm", "HeNe", 15797.6, 15798.8],
  ["671 nm", "Li cooling", 14900.0, 14901.2],
  ["780 nm", "Rb D2", 12816.0, 12817.2],
];

// --- data ---------------------------------------------------------------------------------------

const url = p => new URL("data/" + p, import.meta.url).href;

const inflight = new Map();             // path -> promise, so a file asked for twice is fetched once
async function getJSON(path, store, key = path) {
  if (store && store.has(key)) return store.get(key);
  if (!inflight.has(path)) {
    inflight.set(path, (async () => {
      const r = await fetch(url(path));
      if (!r.ok) throw new Error(`${path}: ${r.status}`);
      return r.json();
    })().finally(() => inflight.delete(path)));
  }
  const data = await inflight.get(path);
  if (store) store.set(key, data);
  return data;
}

function partitionAt(iso, T) {
  const key = `${iso}|${T}`;
  if (!state.qCache.has(key)) state.qCache.set(key, P.partitionFunction(state.partitions.get(iso), T));
  return state.qCache.get(key);
}

/** Row -> its measurements, per shard: [{s: source, f: [[component, MHz, u, obs - model]], i: intervals}]. */
function measOf(d) {
  if (!d._meas) d._meas = new Map(d.m.map((i, k) => [i, d.ms ? d.ms[k] : []]));
  return d._meas;
}

/** Decode one shard row into the object the rest of the app uses. */
const rowAt = (d, i, iso) => {
  const meas = measOf(d).get(i) || null;
  return {
    iso, nu: d.nu[i], s0: d.s0[i], el: d.el[i], vu: d.vu[i], vl: d.vl[i],
    J: Math.abs(d.j[i]), branch: d.j[i] > 0 ? "R" : "P", meas, measured: meas !== null, u: d.u[i], w: d.w[i],
    d08: d.d08 ? d.d08[i] : null,          // the published 2008 model's position minus this one, MHz
    hfsFile: d._hfsFile, row: i,           // where this line's ΔJ = 0 hyperfine pattern is
  };
};

/** A line's ΔJ = 0 hyperfine pattern, {o: MHz offsets, s: per mille}, if its shard's file is loaded. */
function hfsOf(line) {
  const h = line.hfsFile && state.hfsShards.get(line.hfsFile);
  if (!h || !h.o[line.row] || !h.o[line.row].length) return null;
  return {o: h.o[line.row], s: h.s[line.row], x: h.x ? h.x[line.row] : null};
}

// Doppler profiles resolve components below 0.1 nm. Sub-Doppler views allow 5 cm⁻¹,
// including the initial ±2 cm⁻¹ window around P(53) 32-0 (~0.113 nm wide).
const HFS_SPAN_NM = 0.1;
const SUB_DOPPLER_SPAN_CM = 5;
const viewSpanNm = () => 1e7 / state.lo - 1e7 / state.hi;
const hyperfineShown = () => state.mode === "sub" ? state.hi - state.lo <= SUB_DOPPLER_SPAN_CM :
  viewSpanNm() < HFS_SPAN_NM * (1 - 1e-6);   // a Doppler view of exactly 0.1 nm is not "below"

const source = k => (state.refs && state.refs.sources[k]) || {short: `source ${k}`, link: null};
const escapeHTML = t => String(t).replace(/[&<>"]/g, c => ({"&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;"}[c]));
/** A source's short name, linked to its DOI (or its document when it has none). */
function sourceLink(k) {
  const s = source(k);
  const name = escapeHTML(s.short);
  return s.link ? `<a href="${escapeHTML(s.link)}" target="_blank" rel="noopener" title="${escapeHTML(s.citation || "")}">${name}</a>`
    : `<span title="${escapeHTML(s.citation || "")}">${name}</span>`;
}
/** "precision", "atlas" or "model": how a line is known. A line in both kinds of source counts as precision. */
function kindOf(l) {
  if (!l.meas) return "model";
  if (l._kind) return l._kind;
  const kinds = l.meas.map(e => source(e.s).kind);
  return (l._kind = kinds.includes("precision") ? "precision" : kinds.includes("atlas") ? "atlas" : "model");
}
/** Which kinds of line are shown: three independent toggles. state.show is a key for caches. */
const showKey = () => state.kinds.size === 3 ? "all" : [...state.kinds].sort().join(",");
const keepLine = () => (state.kinds.size === 3 ? () => true : l => state.kinds.has(kindOf(l)));
let shownCache = {lines: null, show: null, out: null};
/** Lines whose centre lies within this margin (cm⁻¹) outside the view still draw there: their hyperfine
 * components (up to ~1 GHz from the centre) and Doppler wings reach in. About 1.5 GHz. */
const PLOT_PAD_CM = 0.05;
let plotCache = {lines: null, show: null, out: null};
/** The lines the plot draws: the view's, plus those just outside it (state.plotLines), filtered as shown. */
const plotLines = () => {
  const all = state.plotLines || state.lines;
  if (state.show === "all") return all;
  if (plotCache.lines !== all || plotCache.show !== state.show)
    plotCache = {lines: all, show: state.show, out: all.filter(keepLine())};
  return plotCache.out;
};
const shownLines = () => {
  if (state.show === "all") return state.lines;
  if (shownCache.lines !== state.lines || shownCache.show !== state.show)
    shownCache = {lines: state.lines, show: state.show, out: state.lines.filter(keepLine())};
  return shownCache.out;
};
const KIND_LABEL = {precision: "precision", atlas: "atlas", model: "model only"};
const KINDS = ["precision", "atlas", "model"];
const KIND_INDEX = {precision: 0, atlas: 1, model: 2};

/** The smallest uncertainty among a line's absolute measurements, MHz (Infinity for intervals only). */
const bestU = l => l.meas ? Math.min(...l.meas.flatMap(e => e.f.map(f => f[2])), Infinity) : Infinity;

const label = l => `${l.branch}(${l.J}) ${l.vu}-${l.vl}`;

async function loadShards(iso, files, note) {
  const missing = files.filter(f => !state.shards.has(f));
  for (let k = 0; k < missing.length; k++) {
    if (note) setStatus(`${note} ${k + 1}/${missing.length}…`);
    await getJSON(missing[k], state.shards);
  }
}

/** Fetch shards without waiting for them; redraw the view when they arrive if they are in it. */
function loadInBackground(iso, shards) {
  const missing = shards.filter(s => !state.shards.has(s.file));
  if (!missing.length) return;
  Promise.all(missing.map(s => getJSON(s.file, state.shards).catch(() => null))).then(() => {
    if (state.iso === iso && missing.some(s => s.nu1 >= state.lo && s.nu0 <= state.hi)) refreshView();
  });
}

async function linesIn(iso, nuLo, nuHi, {loadedOnly = false} = {}) {
  const info = state.manifest.isotopologues[iso];
  let need = info.shards.filter(s => s.nu1 >= nuLo && s.nu0 <= nuHi);
  if (loadedOnly) {
    // during a gesture: draw what is here now, and let the rest arrive without holding up the frame
    loadInBackground(iso, need);
    need = need.filter(s => state.shards.has(s.file));
  } else {
    await loadShards(iso, need.map(s => s.file), need.length > 3 ? "loading line list" : null);
  }
  // The export stores the lines in the order the table lists them, increasing wavelength, within
  // each shard and across shards (manifest order). The lines in view are a run of each shard, found by
  // bisection, and their concatenation is already in order: nothing is sorted here. The decoded row
  // objects are kept per shard, so zooming does not rebuild them.
  const out = [];
  out.slices = [];                               // which rows of which shard: for per-shard summaries
  for (const s of need) {
    const d = state.shards.get(s.file);
    d._hfsFile = s.hfs;
    if (!d._rows) d._rows = Array.from(d.nu, (_, i) => rowAt(d, i, iso));
    const a = firstBelow(d.nu, nuHi), b = firstBelow(d.nu, nuLo, true);
    for (let i = a; i < b; i++) out.push(d._rows[i]);
    out.slices.push({d, s, a, b});
  }
  return out;
}

/** In a decreasing array, the first index with arr[i] <= x (or < x when strictly is set). */
function firstBelow(arr, x, strictly = false) {
  let lo = 0, hi = arr.length;
  while (lo < hi) {
    const mid = (lo + hi) >> 1;
    if (strictly ? arr[mid] >= x : arr[mid] > x) lo = mid + 1; else hi = mid;
  }
  return lo;
}

/** The ΔJ = 0 hyperfine files of the shards a view overlaps (older exports have none: then nothing). */
async function loadHfs(iso, nuLo, nuHi) {
  const need = state.manifest.isotopologues[iso].shards.filter(s => s.hfs && s.nu1 >= nuLo && s.nu0 <= nuHi);
  for (const s of need) {
    try { await getJSON(s.hfs, state.hfsShards); } catch { /* drawn without hyperfine structure */ }
  }
}

/** Build label -> line for one isotopologue. Needs every shard, so it asks first. */
async function buildIndex(iso) {
  if (state.index.has(iso)) return state.index.get(iso);
  const info = state.manifest.isotopologues[iso];
  await loadShards(iso, info.shards.map(s => s.file), "indexing every line");
  const index = new Map();
  for (const s of info.shards) {
    const d = state.shards.get(s.file);
    for (let i = 0; i < d.nu.length; i++) {
      const line = rowAt(d, i, iso);
      index.set(label(line), line);
    }
  }
  state.index.set(iso, index);
  return index;
}

// --- the spectrum --------------------------------------------------------------------------------

/**
 * Bin the lines onto the pixel grid. A line narrower than one pixel becomes a spike of the right
 * area rather than the right height, so the trace keeps its integral as you zoom; a resolved line
 * is spread as the Gaussian it is.
 */
/** erf(x), Abramowitz & Stegun 7.1.26: absolute error below 1.5e-7, far under a pixel. */
function erf(x) {
  const t = 1 / (1 + 0.3275911 * Math.abs(x));
  const y = 1 - t * (0.254829592 + t * (-0.284496736 + t * (1.421413741 + t * (-1.453152027 + t * 1.061405429)))) * Math.exp(-x * x);
  return x >= 0 ? y : -y;
}

/**
 * Add a Gaussian of area S, centre c and width sig to the bins, as its mean over each bin.
 * Averaging over the bin, rather than sampling at its centre, conserves the area and keeps a line's
 * drawn height from jumping as it moves across the pixel grid.
 */
function addGaussian(bins, lo, dnu, c, sig, S) {
  const nb = bins.length, reach = 5 * sig;
  const k0 = Math.max(0, Math.floor((c - reach - lo) / dnu));
  const k1 = Math.min(nb - 1, Math.floor((c + reach - lo) / dnu));
  if (k1 < k0) return;
  const f = Math.SQRT1_2 / sig;
  let prev = erf((lo + k0 * dnu - c) * f);
  for (let k = k0; k <= k1; k++) {
    const next = erf((lo + (k + 1) * dnu - c) * f);
    bins[k] += 0.5 * S * (next - prev) / dnu;
    prev = next;
  }
}

/**
 * The cross section on the pixel grid, one array per kind of line, and the hyperfine components drawn as
 * sticks when the view is narrow enough. A line narrower than a pixel is widened to one pixel in
 * quadrature before it is binned: the pixel-scale smoothing that stops lines a few pixels wide from
 * flickering as the view moves. Wider lines are unaffected (a 10-pixel line gains 0.5% in width).
 */
/** A line's strength at T, kept on the line until T changes (the exponential is the costly part). */
function strengthOf(line, q, T) {
  if (line._T !== T) { line._S = P.strengthAt(line, q, T); line._T = T; }
  return line._S;
}

//: Bins per shard in the precomputed histograms, and how much coarser than them a view's pixels must be
//: before the view is drawn from them: 2048 bins make about 0.07 cm⁻¹ each.
const HIST_BINS = 2048, HIST_MIN_RATIO = 4;

/**
 * A wide view is drawn from per-shard histograms instead of from the lines: each shard's absorption,
 * binned once per temperature and filter on a fine grid, is resampled onto the pixels. At full range
 * that replaces 270 000 line profiles per frame with about 130 000 bin copies; lines there are far
 * narrower than a pixel, so nothing visible is lost. Shards not loaded, or not covering the view, are
 * simply absent, as their lines would be.
 */
function histogramOf(d, s, T, iso) {
  const key = `${T}|${state.show}`;
  if (d._hist && d._hist.key === key) return d._hist;
  const lo = s.nu0, hi = s.nu1, w = (hi - lo) / HIST_BINS;
  const parts = KINDS.map(() => new Float64Array(HIST_BINS));    // area per bin
  const marks = new Uint8Array(HIST_BINS);                        // 1 atlas, 2 precision measured here
  const q = partitionAt(iso, T), keep = keepLine();
  for (const line of d._rows) {
    if (!keep(line)) continue;
    const x = (line.nu - lo) / w - 0.5, k = Math.floor(x), f = x - k;
    const kind = kindOf(line), S = strengthOf(line, q, T), bins = parts[KIND_INDEX[kind]];
    if (k >= 0) bins[k] += S * (1 - f);                          // split between the two nearest bins
    if (k + 1 < HIST_BINS) bins[k + 1] += S * f;
    const km = Math.min(HIST_BINS - 1, Math.max(0, Math.round(x)));
    if (kind === "precision") marks[km] = 2; else if (kind === "atlas" && !marks[km]) marks[km] = 1;
  }
  return (d._hist = {key, lo, w, parts, marks});
}

function binnedFromHistograms(lo, hi, nbins, T, iso) {
  const parts = KINDS.map(() => new Float64Array(nbins));
  const dnu = (hi - lo) / nbins;
  const info = state.manifest.isotopologues[iso];
  for (const s of info.shards) {
    if (s.nu1 < lo || s.nu0 > hi) continue;
    const d = state.shards.get(s.file);
    if (!d || !d._rows) continue;
    const h = histogramOf(d, s, T, iso);
    for (let p = 0; p < 3; p++) {
      const src = h.parts[p], dst = parts[p];
      for (let j = 0; j < HIST_BINS; j++) {
        const A = src[j];
        if (!A) continue;
        const x = (h.lo + (j + 0.5) * h.w - lo) / dnu - 0.5, k = Math.floor(x), f = x - k;
        if (k >= 0 && k < nbins) dst[k] += A * (1 - f) / dnu;
        if (k + 1 >= 0 && k + 1 < nbins) dst[k + 1] += A * f / dnu;
      }
    }
  }
  return {parts, sticks: []};
}

function binned(lines, lo, hi, nbins, T, iso, withHfs) {
  const dnu = (hi - lo) / nbins;
  const s0 = state.manifest.isotopologues[iso].shards[0];
  if (!withHfs && dnu >= HIST_MIN_RATIO * (s0.nu1 - s0.nu0) / HIST_BINS) return binnedFromHistograms(lo, hi, nbins, T, iso);
  const parts = KINDS.map(() => new Float64Array(nbins));
  const sticks = [];
  const q = partitionAt(iso, T);
  const root2pi = Math.sqrt(2 * Math.PI);
  for (const line of lines) {
    const kind = kindOf(line);
    const bins = parts[KIND_INDEX[kind]];
    const S = strengthOf(line, q, T);
    const sigma = P.dopplerFWHM(line.nu, T, iso) / 2.35482;
    const sig = Math.hypot(sigma, dnu);
    const h = withHfs ? hfsOf(line) : null;
    if (!h) { addGaussian(bins, lo, dnu, line.nu, sig, S); continue; }
    for (let c = 0; c < h.o.length; c++) {
      const nuC = line.nu + h.o[c] / P.MHZ_PER_CM, Sc = S * h.s[c] / 1000;
      addGaussian(bins, lo, dnu, nuC, sig, Sc);
      sticks.push({nu: nuC, peak: Sc / (sigma * root2pi), kind});
    }
  }
  return {parts, sticks};
}

/** The key of a line in the full hyperfine patterns (hfs_<iso>.json). */
const hfsKey = l => `${l.vu}-${l.vl}${l.branch}${l.J}`;

/**
 * Components with the level identities physics.resonances() compares: main component k joins upper and
 * lower level k (no two share a level), and each weak one [offset, strength, upper, lower] of the export
 * carries the numbers of its levels on the same scheme (webapp.weak_links). `scale` puts the main
 * strengths on the weak ones' scale (1/1000 for the per-mille shard patterns).
 */
function withLevels(o, s, x, scale = 1) {
  const comps = o.map((offset, k) => ({offset, strength: s[k] * scale, upper: k, lower: k}));
  for (const [offset, strength, upper, lower] of x || []) comps.push({offset, strength, upper, lower});
  return comps;
}

/**
 * A line's hyperfine components as {offset (MHz), strength, upper, lower}: the full calculation when the
 * detail pane has loaded it, else the ΔJ = 0 pattern of its shard, else the line alone. Each carries the
 * weak ΔF ≠ ΔJ components that share a level with a main one (10⁻² to 10⁻⁴ of the line), so
 * physics.resonances() finds the crossovers as well as the Lamb dips (docs/research/sub-doppler.md §2).
 */
function componentsOf(line) {
  const full = (state.hfs.get(line.iso) || {})[hfsKey(line)];
  if (full) return withLevels(full.o, full.s, full.x);
  const h = hfsOf(line);
  if (h) return withLevels(h.o, h.s, h.x, 1 / 1000);
  return [{offset: 0, strength: 1}];
}

/**
 * The saturated-absorption signal on the pixel grid (physics.resonances and physics.addSignal, the model of
 * i2spec/saturation.py). Each line's components are weighted by its strength at T, so the dips of
 * different lines compare as S² in the weak-saturation limit. A resonance narrower than about three
 * pixels is drawn at three pixels (the modulation widened in proportion, keeping the shape) so that it
 * shows at all; the note on the plot says so.
 */
function subDopplerBins(lines, lo, hi, nbins, T, iso) {
  const {fwhm, harmonic, modulation} = state.sub;
  const q = partitionAt(iso, T);
  const dx = (hi - lo) * P.MHZ_PER_CM / nbins, x0 = lo * P.MHZ_PER_CM;
  const drawn = Math.hypot(fwhm, 3 * dx), widen = drawn / fwhm;
  const bins = new Float64Array(nbins);
  let nDips = 0, nCross = 0;
  for (const line of lines) {
    const S = strengthOf(line, q, T), f0 = line.nu * P.MHZ_PER_CM;
    const comps = componentsOf(line).map(c => ({...c, strength: S * c.strength}));
    const {list, total} = P.resonances(comps, P.dopplerFWHM(line.nu, T, iso) * P.MHZ_PER_CM);
    for (const r of list) if (r.kind === "crossover") nCross++; else nDips++;
    P.addSignal(bins, x0, dx, list.map(r => ({offset: f0 + r.offset, amplitude: r.amplitude})), drawn,
      {harmonic, modulation: harmonic ? modulation * widen : 0, scale: total});
  }
  return {bins, harmonic, fwhm, fwhmDrawn: drawn, nDips, nCross};
}

/**
 * Per shard and filter: running totals of each kind of line, and the rows of the measured lines, so
 * the counts of any run of rows are two subtractions and its measured lines a slice.
 */
function shardTotals(d) {
  if (d._tot && d._tot.show === state.show) return d._tot;
  const n = d._rows.length, keep = keepLine();
  const cum = KINDS.map(() => new Int32Array(n + 1)), meas = [];
  for (let i = 0; i < n; i++) {
    const l = d._rows[i], k = keep(l) ? KIND_INDEX[kindOf(l)] : -1;
    for (let p = 0; p < 3; p++) cum[p][i + 1] = cum[p][i] + (p === k ? 1 : 0);
    if (k >= 0 && l.measured) meas.push(i);
  }
  return (d._tot = {show: state.show, cum, meas});
}

/** Counts by kind of the lines shown, and their measured lines, for the current view. */
let summaryCache = {lines: null, show: null};
function viewSummary() {
  const lines = state.lines;
  if (summaryCache.lines === lines && summaryCache.show === state.show) return summaryCache;
  const counts = {precision: 0, atlas: 0, model: 0};
  let nMarked = 0;
  const slices = lines.slices || [];
  for (const {d, a, b} of slices) {
    const t = shardTotals(d);
    KINDS.forEach((k, p) => { counts[k] += t.cum[p][b] - t.cum[p][a]; });
    nMarked += lowerIndex(t.meas, b) - lowerIndex(t.meas, a);
  }
  // the measured lines themselves, only when few enough to draw one by one
  let marked = null;
  if (nMarked <= 5000) {
    marked = [];
    for (const {d, a, b} of slices) {
      const t = shardTotals(d);
      for (let j = lowerIndex(t.meas, a); j < t.meas.length && t.meas[j] < b; j++) marked.push(d._rows[t.meas[j]]);
    }
  }
  return (summaryCache = {lines, show: state.show, counts, nMarked, marked});
}
function lowerIndex(arr, x) {                     // first index with arr[i] >= x, arr ascending
  let lo = 0, hi = arr.length;
  while (lo < hi) { const mid = (lo + hi) >> 1; if (arr[mid] < x) lo = mid + 1; else hi = mid; }
  return lo;
}

/** Round tick values (1, 2 or 5 times a power of ten) between a and b, about n of them. */
function niceTicks(a, b, n) {
  const raw = (b - a) / n, p10 = 10 ** Math.floor(Math.log10(raw)), f = raw / p10;
  const step = (f < 1.5 ? 1 : f < 3.5 ? 2 : f < 7.5 ? 5 : 10) * p10;
  const decimals = Math.max(0, -Math.floor(Math.log10(step) + 1e-9));
  const out = [];
  for (let v = Math.ceil(a / step - 1e-9) * step; v <= b + step * 1e-9; v += step) out.push(+v.toFixed(decimals));
  return {values: out, decimals};
}

const UNIT_TITLE = {"nm": "vacuum wavelength (nm)", "nm-air": "wavelength in air (nm)", "cm-1": "wavenumber (cm⁻¹)",
                    "THz": "frequency (THz)", "MHz": "frequency (MHz)"};
/** The quantity on the top axis: wavenumber, or vacuum wavelength when the main unit is wavenumber. */
const secondUnit = unit => unit === "cm-1" ? "nm" : "cm-1";

/**
 * Both horizontal axes: the chosen unit below the plot, and wavenumber (or wavelength) above it, each
 * with round tick values placed where they fall. Labels that would overlap a neighbour are left out.
 */
function drawXLabels(g, lo, hi, m, pw, H, narrow, flip) {
  const unit = $("#unit").value;
  const xOf = nu => m.l + pw * (flip ? (hi - nu) : (nu - lo)) / (hi - lo);
  g.fillStyle = ink("--ink-3"); g.strokeStyle = ink("--rule");
  g.font = '11px "IBM Plex Mono", ui-monospace, monospace'; g.lineWidth = 1;
  for (const [u, yText, yTick] of [[unit, H - 16, m.t + m.ph + 16], [secondUnit(unit), m.t - 7, m.t - 4]]) {
    const a = P.fromWavenumber(lo, u), b = P.fromWavenumber(hi, u);
    const {values, decimals} = niceTicks(Math.min(a, b), Math.max(a, b), narrow ? 3 : 6);
    let lastRight = -Infinity;
    g.textAlign = "center"; g.textBaseline = "alphabetic";
    // left to right as drawn (a reversed axis lists its values the other way), so the overlap test works
    const placed = values.map(v => [v, xOf(P.toWavenumber(v, u))]).sort((p, q) => p[1] - q[1]);
    for (const [v, x] of placed) {
      if (x < m.l - 0.5 || x > m.l + pw + 0.5) continue;
      g.beginPath(); g.moveTo(Math.round(x) + 0.5, yTick); g.lineTo(Math.round(x) + 0.5, yTick + 4); g.stroke();
      const text = v.toFixed(decimals), w = g.measureText(text).width;
      const cx = Math.min(Math.max(x, m.l + w / 2), m.l + pw - w / 2);        // keep the ends inside
      if (cx - w / 2 < lastRight + 10) continue;
      g.fillText(text, cx, yText);
      lastRight = cx + w / 2;
    }
  }
  g.textAlign = "center"; g.textBaseline = "bottom";
  g.fillText(UNIT_TITLE[unit], m.l + pw / 2, H - 2);
  g.textAlign = "right"; g.textBaseline = "middle";
  g.fillText(UNIT_TITLE[secondUnit(unit)], m.l + pw, narrow ? 26 : 10.5);   // the top axis is named at the right
  g.textBaseline = "alphabetic";
}

/**
 * The last finished picture, and a preview made from it: the plot area stretched so that each
 * wavenumber lands where it would in the new range (blank where the old picture had nothing), with
 * the axis labels drawn fresh. Returns false when there is nothing to stretch.
 */
let lastFrame = null;
function keepFrame(canvas, meta) { lastFrame = meta; }
function drawPreview(lo, hi) {
  const f = lastFrame, canvas = $("#plot");
  if (!f || f.iso !== state.iso || f.W !== $("#plotwrap").clientWidth) return false;
  const {m, pw, ph, H, dpr, flip, narrow} = f;
  // where a wavenumber was drawn, and where it belongs now: x' = m.l + (x - m.l) sx + tx
  const sx = (f.hi - f.lo) / (hi - lo);
  const tx = pw * (flip ? (hi - f.hi) : (f.lo - lo)) / (hi - lo);
  $("#trace").style.transform = `translateX(${tx}px) scaleX(${sx})`;
  const g = canvas.getContext("2d");
  g.setTransform(dpr, 0, 0, dpr, 0, 0);
  g.clearRect(0, m.t + ph + 14, f.W, H - (m.t + ph + 14));        // both axes' labels, redrawn for the new range
  g.clearRect(0, m.t - 20, f.W, 20);
  drawXLabels(g, lo, hi, m, pw, H, narrow, flip);
  return true;
}

// Reading a CSS variable costs a style lookup, so the colours are read once per draw, not per mark.
let inkCache = null;
const ink = name => {
  if (!inkCache) inkCache = getComputedStyle(document.documentElement);
  return inkCache.getPropertyValue(name).trim();
};

function drawPlot() {
  inkCache = null;
  const kindInk = Object.fromEntries(KINDS.map(k => [k, ink(`--cls-${k}`)]));
  const canvas = $("#plot"), wrap = $("#plotwrap");
  // at most 2 canvas pixels per CSS pixel: a phone's 3x screen would otherwise make every frame
  // fill 2.25 times as many pixels for a difference no one can see in a line plot
  const dpr = Math.min(window.devicePixelRatio || 1, 2);
  const W = wrap.clientWidth, H = wrap.clientHeight;
  if (W < 40 || H < 40) return;
  const narrow = W < 700;
  // on a phone the vertical axis is labelled inside the plot, so the spectrum gets the full width
  const m = {l: narrow ? 8 : 78, r: narrow ? 8 : 18, t: narrow ? 56 : 42, b: 50};
  const pw = W - m.l - m.r, ph = H - m.t - m.b;
  // no room for a plot (a phone with the options panel open): drawing anyway would turn it upside
  // down, so keep the last picture until the plot is resized and redrawn (see wire)
  if (pw < 20 || ph < 20) return;
  // resizing a canvas reallocates its pixels, so only when the size has actually changed
  if (canvas.width !== Math.round(W * dpr) || canvas.height !== Math.round(H * dpr)) {
    canvas.width = Math.round(W * dpr); canvas.height = Math.round(H * dpr);
    canvas.style.width = W + "px"; canvas.style.height = H + "px";
  }
  const gb = canvas.getContext("2d");
  gb.setTransform(dpr, 0, 0, dpr, 0, 0);
  gb.clearRect(0, 0, W, H);
  let g = gb;                  // the axes, grid and labels; the trace itself goes on its own layer (gt)
  m.ph = ph;
  const {lo, hi, T, iso} = state;
  // The trace layer: a canvas the size of the plot, inside a box clipped to the plot area. A zoom or
  // pan in progress moves it with a CSS transform, which the compositor does without redrawing.
  const tc = $("#trace"), clip = $("#traceclip");
  if (tc.width !== canvas.width || tc.height !== canvas.height) {
    tc.width = canvas.width; tc.height = canvas.height;
    tc.style.width = W + "px"; tc.style.height = H + "px";
  }
  Object.assign(clip.style, {left: m.l + "px", width: pw + "px", height: H + "px"});
  Object.assign(tc.style, {left: -m.l + "px", transform: "none", transformOrigin: `${m.l}px 0`});
  const gt = tc.getContext("2d");
  gt.setTransform(dpr, 0, 0, dpr, 0, 0);
  gt.clearRect(0, 0, W, H);
  const lines = plotLines();
  const nb = Math.max(2, Math.round(pw * dpr));        // one bin per canvas pixel
  // one binned cross section per kind of line, stacked precision / atlas / model from the baseline
  const withHfs = hyperfineShown();
  const subOn = state.mode === "sub" && withHfs;
  const {parts, sticks} = binned(lines, lo, hi, nb, T, iso, withHfs);
  const sub = subOn ? subDopplerBins(lines, lo, hi, nb, T, iso) : null;
  const sigma = new Float64Array(nb);
  for (const p of parts) for (let k = 0; k < nb; k++) sigma[k] += p[k];

  let y, ylabel, ticks, top = 1;
  if (sub) {
    // relative: the resonances of each line scale with its strength squared (the bilinear limit)
    let peak = 0;
    for (const v of sub.bins) peak = Math.max(peak, Math.abs(v));
    peak = peak || 1;
    if (sub.harmonic) {
      y = Array.from(sub.bins, v => 0.5 + 0.5 * v / (peak * 1.08));
      ticks = [[0.5 - 0.5 / 1.08, "−1"], [0.5, "0"], [0.5 + 0.5 / 1.08, "1"]];
      ylabel = `${sub.harmonic}f signal (rel.)`;
    } else {
      y = Array.from(sub.bins, v => v / (peak * 1.08));
      ticks = [0, 0.25, 0.5, 0.75, 1].map(v => [v / 1.08, v.toFixed(2)]);
      ylabel = "sat. signal (rel.)";
    }
  } else if (state.mode === "trans") {
    const N = P.numberDensity(Number($("#cold").value) + 273.15, T) * Number($("#path").value);
    y = Array.from(sigma, s => Math.exp(-s * N));
    ylabel = "transmission";
    ticks = [0, 0.25, 0.5, 0.75, 1].map(v => [v, v.toFixed(2)]);
  } else {
    top = 0;
    for (const v of sigma) top = Math.max(top, v);
    top = top > 0 ? top * 1.08 : 1e-20;
    y = Array.from(sigma, s => s / top);
    ylabel = "cross section (cm²)";
    // Round ticks: a step of 1, 2 or 5 times a power of ten, about four of them. In units of that
    // power of ten every tick is an integer, so a phone shows the integers and names the unit once.
    const raw = top / 4, p10 = 10 ** Math.floor(Math.log10(raw)), f = raw / p10;
    const mult = f < 1.5 ? 1 : f < 3.5 ? 2 : f < 7.5 ? 5 : 10;
    const step = mult * p10, e = Math.round(Math.log10(p10)) + (mult === 10 ? 1 : 0);
    const values = [];
    for (let k = 0; k * step <= top * (1 + 1e-9); k++) values.push(k * step);
    if (narrow) {
      ticks = values.map(v => [v / top, String(Math.round(v / 10 ** e))]);
      ylabel = `σ (10${String(e).replace("-", "⁻").replace(/\d/g, d => "⁰¹²³⁴⁵⁶⁷⁸⁹"[d])} cm²)`;
    } else {
      ticks = values.map(v => [v / top, v === 0 ? "0" : `${Math.round(v / 10 ** e)}e${e}`]);
    }
  }

  // a wavelength axis has to run the other way round, or it reads right to left
  const flip = state.unitFlips;
  const xPix = nu => m.l + pw * (flip ? (hi - nu) : (nu - lo)) / (hi - lo);
  const yPix = v => m.t + ph * (1 - v);
  g.strokeStyle = ink("--rule-2"); g.lineWidth = 1;
  g.beginPath();
  for (const [v] of ticks) { g.moveTo(m.l, yPix(v)); g.lineTo(m.l + pw, yPix(v)); }
  g.stroke();

  g.fillStyle = ink("--ink-3"); g.font = '11px "IBM Plex Mono", ui-monospace, monospace';
  if (narrow) {
    // inside the plot, just above each grid line, and the axis name in the free row above the plot
    g.textAlign = "left"; g.textBaseline = "bottom";
    g.strokeStyle = ink("--surface"); g.lineWidth = 3; g.lineJoin = "round";   // a halo, where a line passes behind
    for (const [v, text] of ticks) if (v > 0) { g.strokeText(text, m.l + 2, yPix(v) - 2); g.fillText(text, m.l + 2, yPix(v) - 2); }
    g.fillText(ylabel, m.l, 30);
  } else {
    g.textAlign = "right"; g.textBaseline = "middle";
    for (const [v, text] of ticks) g.fillText(text, m.l - 10, yPix(v));
  }
  g.textAlign = "center"; g.textBaseline = "alphabetic";
  const unit = $("#unit").value, digits = P.DIGITS[unit];
  drawXLabels(g, lo, hi, m, pw, H, narrow, flip);
  if (!narrow) { g.save(); g.translate(15, m.t + ph / 2); g.rotate(-Math.PI / 2); g.fillText(ylabel, 0, 0); g.restore(); }

  g = gt;                      // from here to the marker strip: the trace layer
  const xAt = k => m.l + (k + 0.5) * pw / nb;
  const at = (arr, k) => arr[flip ? nb - 1 - k : k];
  if (sub) {
    if (!sub.harmonic) {
      g.beginPath(); g.moveTo(xAt(0), yPix(0));
      for (let k = 0; k < nb; k++) g.lineTo(xAt(k), yPix(at(y, k)));
      g.lineTo(xAt(nb - 1), yPix(0)); g.closePath();
      g.fillStyle = ink("--trace-fill"); g.fill();
    }
    g.beginPath();
    for (let k = 0; k < nb; k++) k === 0 ? g.moveTo(xAt(k), yPix(at(y, k))) : g.lineTo(xAt(k), yPix(at(y, k)));
    g.strokeStyle = ink("--trace"); g.lineWidth = 1.4; g.lineJoin = "round"; g.stroke();
  } else if (state.mode === "trans") {
    g.beginPath();
    for (let k = 0; k < nb; k++) k === 0 ? g.moveTo(xAt(k), yPix(at(y, k))) : g.lineTo(xAt(k), yPix(at(y, k)));
    g.strokeStyle = ink("--ink-2"); g.lineWidth = 1.6; g.lineJoin = "round"; g.stroke();
  } else {
    // stacked areas, each kind's share of the cross section in its own colour, then the total on top
    let base = new Float64Array(nb);
    for (let p = 0; p < KINDS.length; p++) {
      const upper = base.map((b, k) => b + parts[p][k]);
      if (upper.some((u, k) => u > base[k])) {
        g.beginPath();
        for (let k = 0; k < nb; k++) g.lineTo(xAt(k), yPix(at(upper, k) / top));
        for (let k = nb - 1; k >= 0; k--) g.lineTo(xAt(k), yPix(at(base, k) / top));
        g.closePath();
        g.fillStyle = ink(`--cls-${KINDS[p]}`); g.globalAlpha = KINDS[p] === "model" ? 0.45 : 0.8; g.fill();
        g.globalAlpha = 1;
      }
      base = upper;
    }
    g.beginPath();
    for (let k = 0; k < nb; k++) k === 0 ? g.moveTo(xAt(k), yPix(at(y, k))) : g.lineTo(xAt(k), yPix(at(y, k)));
    g.strokeStyle = ink("--ink-2"); g.lineWidth = 1; g.lineJoin = "round"; g.stroke();
    // hyperfine components: a stick at each, as tall as that component's own Doppler profile
    if (sticks.length && sticks.length < 6000) {
      g.lineWidth = 1;
      for (const st of sticks) {
        const x = xPix(st.nu);
        if (x < m.l || x > m.l + pw) continue;
        g.strokeStyle = kindInk[st.kind];
        g.globalAlpha = st.kind === "model" ? 0.7 : 0.95;
        g.beginPath(); g.moveTo(x, yPix(0)); g.lineTo(x, yPix(Math.min(st.peak / top, 1))); g.stroke();
      }
      g.globalAlpha = 1;
    }
  }

  // measured lines: a diamond each in a strip under the trace, one per pixel column at most
  const strip = m.t + ph + 8, seen = new Map();
  const {marked, counts, nMarked} = viewSummary();
  const nMeasured = nMarked;
  // one diamond per pixel column, precision winning a shared column; one path per colour
  const mark = (nu, kind) => {
    const x = Math.round(xPix(nu));
    if (x < m.l || x > m.l + pw || seen.get(x) === "precision") return;
    seen.set(x, kind);
  };
  if (marked) for (const line of marked) mark(line.nu, kindOf(line));
  else {
    // too many to place one by one: where the per-shard histograms recorded them
    for (const {d, s} of (state.lines.slices || [])) {
      const h = histogramOf(d, s, T, iso);
      for (let j = 0; j < HIST_BINS; j++) if (h.marks[j]) mark(h.lo + (j + 0.5) * h.w, h.marks[j] === 2 ? "precision" : "atlas");
    }
  }
  for (const k of ["atlas", "precision"]) {
    g.beginPath();
    for (const [x, kind] of seen) {
      if (kind !== k) continue;
      g.moveTo(x, strip - 4); g.lineTo(x + 3.5, strip); g.lineTo(x, strip + 4); g.lineTo(x - 3.5, strip); g.closePath();
    }
    g.fillStyle = kindInk[k]; g.fill();
  }
  g = gb;                      // labels and legend: the axes layer

  // legend, top left: what the colours mean and how many of each are in view
  g.font = '11px "IBM Plex Sans", system-ui, sans-serif'; g.textAlign = "left"; g.textBaseline = "middle";
  let lx = m.l + 10;
  for (const k of KINDS) {
    if (!counts[k]) continue;
    const text = narrow ? `${k === "model" ? "model" : KIND_LABEL[k]} ${fmt(counts[k])}` : `${KIND_LABEL[k]} ${fmt(counts[k])}`;
    if (lx + 14 + g.measureText(text).width > W - 4) break;      // no room: leave the rest out
    g.fillStyle = ink(`--cls-${k}`); g.globalAlpha = k === "model" ? 0.6 : 0.9;
    g.fillRect(lx, 5, 10, 10); g.globalAlpha = 1;
    g.fillStyle = ink("--ink-2");
    g.fillText(text, lx + 14, 10.5);
    lx += 14 + g.measureText(text).width + 14;
  }
  if (nMeasured && lx + 14 + g.measureText("measured").width <= W - 4) {
    // the diamonds under the trace
    g.fillStyle = ink("--ink-2");
    g.beginPath(); g.moveTo(lx + 5, 5.5); g.lineTo(lx + 9, 10); g.lineTo(lx + 5, 14.5); g.lineTo(lx + 1, 10); g.closePath(); g.fill();
    g.fillText("measured", lx + 14, 10.5);
  }
  g.textBaseline = "alphabetic"; g.textAlign = "left";
  g.strokeStyle = ink("--rule"); g.lineWidth = 1;
  g.beginPath(); g.moveTo(m.l, m.t + ph + 0.5); g.lineTo(m.l + pw, m.t + ph + 0.5); g.stroke();

  g = gt;                      // the selected line moves with the trace
  if (state.selected) {
    const x = xPix(state.selected.nu);
    if (x >= m.l && x <= m.l + pw) {
      g.strokeStyle = ink("--pick"); g.lineWidth = 2; g.setLineDash([5, 4]);
      g.beginPath(); g.moveTo(x, m.t); g.lineTo(x, m.t + ph); g.stroke();
      g.setLineDash([]);
      g.textAlign = "center"; g.textBaseline = "top";
      g.font = '600 11px "IBM Plex Mono", ui-monospace, monospace';
      const lx = Math.min(Math.max(x, m.l + 40), m.l + pw - 40);
      g.strokeStyle = ink("--surface"); g.lineWidth = 4; g.lineJoin = "round";
      g.strokeText(label(state.selected), lx, m.t + 4);                // a halo over the trace
      g.fillStyle = ink("--pick"); g.fillText(label(state.selected), lx, m.t + 4);
      g.textBaseline = "alphabetic";
    }
  }
  if (state.laser) {
    // the harmonic of the laser: dashed, labelled at the foot of the plot so it clears the selected line's label
    const x = xPix(state.laser.nu);
    if (x >= m.l && x <= m.l + pw) {
      g.strokeStyle = ink("--accent"); g.lineWidth = 1.5; g.setLineDash([2, 3]);
      g.beginPath(); g.moveTo(x, m.t); g.lineTo(x, m.t + ph); g.stroke();
      g.setLineDash([]);
      const text = `laser ×${state.laser.n}`;
      g.font = '500 11px "IBM Plex Mono", ui-monospace, monospace'; g.textBaseline = "bottom";
      g.textAlign = x > m.l + pw - 80 ? "right" : "left";
      const tx = x + (g.textAlign === "right" ? -4 : 4);
      g.strokeStyle = ink("--surface"); g.lineWidth = 4; g.lineJoin = "round";
      g.strokeText(text, tx, m.t + ph - 4);
      g.fillStyle = ink("--accent"); g.fillText(text, tx, m.t + ph - 4);
      g.textBaseline = "alphabetic"; g.textAlign = "left";
    }
  }
  // lines whose strength was not computed (emission lines, X v'' >= 48) draw nothing; when they are all the
  // view holds, a note says why the plot is empty (the measured-line diamonds below still mark them)
  const silent = lines.filter(l => !(l.s0 > 0));
  if (silent.length) {
    if (silent.length === lines.length) {
      const text = "no thermal absorption here: emission lines (X v″ ≥ 48), marked below the plot";
      g.font = '11px "IBM Plex Sans", system-ui, sans-serif'; g.textAlign = "left"; g.textBaseline = "top";
      g.strokeStyle = ink("--surface"); g.lineWidth = 4; g.lineJoin = "round";
      g.strokeText(text, m.l + 8, m.t + 24); g.fillStyle = ink("--ink-2"); g.fillText(text, m.l + 8, m.t + 24);
      g.textBaseline = "alphabetic";
    }
  }
  if (state.mode === "sub") {
    // what the sub-Doppler trace is, or why it is not drawn: inside the plot, with a halo; top left, or on a
    // phone, where the tick labels are inside the plot at the left, top right
    const note = !sub ? `sub-Doppler: zoom to a span of ${SUB_DOPPLER_SPAN_CM} cm⁻¹ or less` :
      `${sub.nDips} Lamb dips${sub.nCross ? `, ${sub.nCross} crossovers` : ""} · Γ ${sub.fwhm} MHz` + (sub.fwhmDrawn > sub.fwhm * 1.05 ?
        ` (drawn ${sub.fwhmDrawn < 10 ? sub.fwhmDrawn.toFixed(1) : sub.fwhmDrawn.toFixed(0)} MHz: zoom in)` : "");
    g.font = '11px "IBM Plex Sans", system-ui, sans-serif'; g.textAlign = narrow ? "right" : "left"; g.textBaseline = "top";
    const nx = narrow ? m.l + pw - 4 : m.l + 8;
    g.strokeStyle = ink("--surface"); g.lineWidth = 4; g.lineJoin = "round";
    g.strokeText(note, nx, m.t + 6);
    g.fillStyle = ink(sub ? "--ink-2" : "--warn"); g.fillText(note, nx, m.t + 6);
    g.textAlign = "left";
    g.textBaseline = "alphabetic";
  }
  state.geom = {m, pw, ph, y, flip};
  keepFrame(canvas, {lo, hi, iso, m, pw, ph, W, H, dpr, flip, narrow});
}

// --- the results table ---------------------------------------------------------------------------

const COLUMNS = [
  {key: "label", name: "line", cell: l => label(l) + (l.measured ? " ◆" : ""), sort: null},
  {key: "nu", name: "λ vac (nm)", cell: l => (1e7 / l.nu).toFixed(6)},
  {key: "cm", name: "ν (cm⁻¹)", cell: l => l.nu.toFixed(6)},
  {key: "f", name: "f (MHz)", cell: l => fmt(l.nu * P.MHZ_PER_CM, 1)},
  // s0 = 0: a measured line outside the computed list (emission to X v'' >= 48), whose strength was not computed
  {key: "S", name: "S (cm)", cell: l => l.s0 > 0 ? P.strengthAt(l, partitionAt(l.iso, state.T), state.T).toExponential(2) : "—"},
  {key: "u", name: "u (MHz)", cell: l => fmtU(unc(l).value)},
  {key: "meas", name: "measured by", html: l => measuredCell(l)},
];

function measuredCell(l) {
  if (!l.meas) return "";
  const links = l.meas.map(e => sourceLink(e.s));
  return links.length > 3 ? links.slice(0, 3).join(", ") + ` +${links.length - 3}` : links.join(", ");
}
/** 1σ position uncertainty (MHz), why, and the flags that go with it: computed by the Python model at export. */
const unc = l => {
  const r = state.manifest.reasons[l.w] || {why: "", flags: []};
  return {value: l.u, why: r.why, flags: r.flags};
};

const byLam = (a, b) => b.nu - a.nu;       // ascending wavelength
const SORTERS = {
  nu: byLam,
  cm: (a, b) => a.nu - b.nu,
  f: (a, b) => a.nu - b.nu,
  S: (a, b) => P.strengthAt(a, partitionAt(a.iso, state.T), state.T) - P.strengthAt(b, partitionAt(b.iso, state.T), state.T),
  u: (a, b) => unc(a).value - unc(b).value,
  label: (a, b) => a.vu - b.vu || a.vl - b.vl || a.J - b.J,
  meas: (a, b) => KINDS.indexOf(kindOf(a)) - KINDS.indexOf(kindOf(b)) || (bestU(a) - bestU(b)),
};

function renderTable() {
  const head = $("#head");
  head.replaceChildren(...COLUMNS.map(c => {
    const th = document.createElement("th");
    th.textContent = c.name + (state.sort.key === c.key ? (state.sort.dir > 0 ? " ▲" : " ▼") : "");
    th.className = [state.sort.key === c.key ? "sorted" : "", c.html ? "src" : ""].join(" ").trim();
    th.tabIndex = 0;
    const flip = () => {
      state.sort = {key: c.key, dir: state.sort.key === c.key ? -state.sort.dir : (c.key === "S" ? -1 : 1)};
      renderTable();
    };
    th.onclick = flip;
    th.onkeydown = e => { if (e.key === "Enter" || e.key === " ") { e.preventDefault(); flip(); } };
    return th;
  }));

  const cmp = tableOrder();
  const pool = shownLines();
  let shown;
  if (state.sort.key === "nu") {
    // the pool is already in wavelength order, as exported: take its head (or, reversed, its tail)
    shown = state.sort.dir > 0 ? pool.slice(0, 500) : pool.slice(-500).reverse();
  } else if (pool.length > 4000) {
    shown = topK(pool, 500, cmp);
  } else {
    shown = [...pool].sort(cmp).slice(0, 500);
  }
  state.shown = shown;                       // the rows on screen, for keyboard navigation
  const body = $("#rows");
  body.replaceChildren(...shown.map((line, i) => {
    const tr = document.createElement("tr");
    tr.className = `k-${kindOf(line)}`;
    tr.line = line;
    const isSelected = state.selected && state.selected.nu === line.nu;
    if (isSelected) {
      tr.setAttribute("aria-selected", "true");
      if (!state.holdScroll) queueMicrotask(() => {
        // scroll the table, never the page: on a phone the page scrolling away from the plot is disorienting
        const w = tr.closest(".tablewrap"), r = tr.getBoundingClientRect(), b = w.getBoundingClientRect();
        const head = w.querySelector("thead").getBoundingClientRect().height;
        if (r.top < b.top + head || r.bottom > b.bottom) w.scrollTop += (r.top + r.bottom) / 2 - (b.top + head + b.bottom) / 2;
        if (state.keyNav) { tr.focus({preventScroll: true}); state.keyNav = false; }
      });
    }
    tr.append(...COLUMNS.map(c => {
      const td = document.createElement("td");
      if (c.html) { td.innerHTML = c.html(line); td.className = "src"; } else td.textContent = c.cell(line);
      return td;
    }));
    tr.tabIndex = isSelected || (!state.selected && i === 0) ? 0 : -1;   // one tab stop into the table
    tr.onclick = () => select(line);
    tr.onkeydown = e => {
      // the arrows walk the list and take the selection with them, rather than scrolling
      const step = {ArrowDown: 1, ArrowUp: -1, PageDown: 10, PageUp: -10}[e.key];
      if (step !== undefined) {
        // past either end of the list in wavelength order: move the window along the spectrum
        if (byWavelength() && ((i === shown.length - 1 && step > 0) || (i === 0 && step < 0))) {
          e.preventDefault(); tablePan(step > 0 ? state.sort.dir : -state.sort.dir, {select: true});
          return;
        }
        const j = Math.min(Math.max(i + step, 0), shown.length - 1);
        if (j !== i) { e.preventDefault(); state.keyNav = true; select(shown[j]); }
        else e.preventDefault();
        return;
      }
      if (e.key === "Home" || e.key === "End") {
        e.preventDefault(); state.keyNav = true;
        select(shown[e.key === "Home" ? 0 : shown.length - 1]);
        return;
      }
      if (e.key === "Enter" || e.key === " ") { e.preventDefault(); select(line); }
    };
    return tr;
  }));
  // in wavelength order the list continues past its ends: a row at each end says so and moves the window
  if (byWavelength()) {
    const edge = (towardLonger, where) => {
      const tr = document.createElement("tr");
      tr.className = "more";
      const td = document.createElement("td");
      td.colSpan = COLUMNS.length;
      const shorter = towardLonger < 0;
      const atEnd = shorter ? state.hi >= state.manifest.isotopologues[state.iso].nu_max - 1e-9
                            : state.lo <= state.manifest.isotopologues[state.iso].nu_min + 1e-9;
      td.textContent = atEnd ? `end of the exported range` :
        `${where === "top" ? "▲" : "▼"} ${shorter ? "shorter" : "longer"} wavelengths: scroll on, or click`;
      tr.append(td);
      if (!atEnd) { tr.onclick = () => tablePan(towardLonger); tr.tabIndex = -1; }
      return tr;
    };
    body.prepend(edge(-state.sort.dir, "top"));
    body.append(edge(state.sort.dir, "bottom"));
  }
  const counts = {precision: 0, atlas: 0, model: 0};
  for (const l of state.lines) counts[kindOf(l)]++;
  setStatus(`${fmt(state.lines.length)} line${state.lines.length === 1 ? "" : "s"} in view: ` +
    `${fmt(counts.precision)} precision, ${fmt(counts.atlas)} atlas, ${fmt(counts.model)} model only` +
    (pool.length > 500 ? `; ${state.sort.key === "S" ? "strongest" : "first"} 500 listed` : "") +
    ` · ${state.manifest.isotopologues[state.iso].label} · ${state.T} K` + (state.laser ? " · " + laserText() : ""));
}

/** "laser 2 × 281.630500 THz = 563.261000 THz (1064.4900 nm → 532.2450 nm)": the fundamental and its harmonic. */
function laserText() {
  const {n, nu} = state.laser;
  const thz = v => P.fromWavenumber(v, "THz").toFixed(6), nm = v => (1e7 / v).toFixed(4);
  return `laser ${n} × ${thz(nu / n)} THz = ${thz(nu)} THz (${nm(nu / n)} nm → ${nm(nu)} nm vac)`;
}

const byWavelength = () => state.sort.key === "nu";
/** The table's sort order, as a comparator (the CSV export uses it too). */
const tableOrder = () => (a, b) => state.sort.dir * ((SORTERS[state.sort.key] || byLam)(a, b) || byLam(a, b));

/**
 * Move the window along the spectrum from the table: towardLonger = +1 continues to longer wavelengths
 * (lower wavenumber), -1 to shorter. The row at the edge of the visible part of the table stays where
 * it is on screen, so the list reads as one continuous list.
 */
async function tablePan(towardLonger, {select: pick = false} = {}) {
  if (state.panning) return;
  const wrap = $(".tablewrap"), box = wrap.getBoundingClientRect();
  const head = wrap.querySelector("thead").getBoundingClientRect().height;
  const rows = [...$("#rows").children].filter(tr => tr.line);
  const down = towardLonger === state.sort.dir;             // continuing below the list, not above it
  const inView = rows.filter(tr => { const r = tr.getBoundingClientRect(); return r.bottom > box.top + head && r.top < box.bottom; });
  const anchorRow = down ? inView[inView.length - 1] : inView[0];
  const span = state.hi - state.lo, margin = span * 0.01;
  let lo, hi;
  if (towardLonger > 0) {             // to lower wavenumber
    hi = anchorRow ? Math.min(anchorRow.line.nu + margin, state.hi - 0.25 * span) : state.hi - 0.5 * span;
    lo = hi - span;
  } else {
    lo = anchorRow ? Math.max(anchorRow.line.nu - margin, state.lo + 0.25 * span) : state.lo + 0.5 * span;
    hi = lo + span;
  }
  const before = anchorRow ? anchorRow.getBoundingClientRect().top - box.top : 0;
  const anchorNu = anchorRow && anchorRow.line.nu;
  state.panning = true; state.holdScroll = true;
  try {
    await show(lo, hi);
    const tr = [...$("#rows").children].find(r => r.line && r.line.nu === anchorNu);
    if (tr) wrap.scrollTop += tr.getBoundingClientRect().top - wrap.getBoundingClientRect().top - before;
    else wrap.scrollTop = down ? 0 : wrap.scrollHeight;
    if (pick) {
      const next = [...$("#rows").children].filter(r => r.line);
      const i = next.findIndex(r => r.line.nu === anchorNu);
      const target = next[i + (down ? 1 : -1)] || next[down ? 0 : next.length - 1];
      if (target) { state.keyNav = true; await select(target.line); }
    }
  } finally {
    state.panning = false; state.holdScroll = false;
  }
}

/** The k smallest of items under cmp, sorted: a bounded max-heap, O(n log k). */
function topK(items, k, cmp) {
  const heap = [];                      // max-heap: the worst kept item at the root
  const up = i => { while (i > 0) { const p = (i - 1) >> 1; if (cmp(heap[i], heap[p]) <= 0) break; [heap[i], heap[p]] = [heap[p], heap[i]]; i = p; } };
  const down = i => {
    for (;;) {
      const l = 2 * i + 1, r = l + 1;
      let m = i;
      if (l < heap.length && cmp(heap[l], heap[m]) > 0) m = l;
      if (r < heap.length && cmp(heap[r], heap[m]) > 0) m = r;
      if (m === i) return;
      [heap[i], heap[m]] = [heap[m], heap[i]]; i = m;
    }
  };
  for (const x of items) {
    if (heap.length < k) { heap.push(x); up(heap.length - 1); }
    else if (cmp(x, heap[0]) < 0) { heap[0] = x; down(0); }
  }
  return heap.sort(cmp);
}

function setStatus(text) { $("#status").textContent = text; }

/** Keep the address bar describing the view, so a range or a line can be sent to someone. */
function pushURL() {
  const unit = $("#unit").value;
  const q = new URLSearchParams({
    unit, iso: state.iso, T: String(state.T),
    from: Math.min(P.fromWavenumber(state.lo, unit), P.fromWavenumber(state.hi, unit)).toFixed(P.DIGITS[unit]),
    to: Math.max(P.fromWavenumber(state.lo, unit), P.fromWavenumber(state.hi, unit)).toFixed(P.DIGITS[unit]),
  });
  if (state.selected) q.set("line", label(state.selected));
  if (state.show !== "all") q.set("show", state.show || "none");
  q.set("y", state.mode);
  if (state.mode === "sub") {
    if (state.sub.fwhm !== P.SUB_DOPPLER_FWHM) q.set("gamma", String(state.sub.fwhm));
    if (state.sub.harmonic) { q.set("det", state.sub.harmonic + "f"); q.set("mod", String(state.sub.modulation)); }
  }
  if (state.laser) {
    q.set("laser", String(state.laser.value)); q.set("lunit", state.laser.unit); q.set("n", String(state.laser.n));
  }
  history.replaceState(null, "", "?" + q);
}

// --- the detail pane -----------------------------------------------------------------------------

function drawHyperfine(canvas, comps, line) {
  const dpr = window.devicePixelRatio || 1;
  const W = canvas.clientWidth, H = 120;
  canvas.width = W * dpr; canvas.height = H * dpr; canvas.style.height = H + "px";
  const g = canvas.getContext("2d"); g.setTransform(dpr, 0, 0, dpr, 0, 0);
  g.clearRect(0, 0, W, H);

  const span = Math.max(...comps.o.map(Math.abs)) * 1.12 || 1;   // MHz
  const pad = 10, base = H - 22, top = 8;
  const xOf = o => pad + (W - 2 * pad) * (o + span) / (2 * span);
  const smax = Math.max(...comps.s) || 1;

  // the Doppler-broadened blend the components add up to, on the same axis, or in sub-Doppler mode the
  // saturation signal of the same components
  const fwhm = P.dopplerFWHM(line.nu, state.T, line.iso) * P.MHZ_PER_CM;
  const sigma = fwhm / 2.35482;
  let prof = new Float64Array(W + 1);
  let pmax = 0, subNote = null;
  const sub = state.mode === "sub";
  if (sub) {
    const {fwhm: gamma, harmonic, modulation} = state.sub;
    const dx = 2 * span / (W - 2 * pad), drawn = Math.hypot(gamma, 3 * dx);
    const {list} = P.resonances(withLevels(comps.o, comps.s, comps.x), fwhm);
    P.addSignal(prof, -span - (pad + 0.5) * dx, dx, list, drawn,
      {harmonic, modulation: harmonic ? modulation * drawn / gamma : 0});
    for (const v of prof) pmax = Math.max(pmax, Math.abs(v));
    subNote = `${harmonic ? harmonic + "f" : "sub-Doppler"} · Γ ${drawn > 1.05 * gamma ? "drawn " + drawn.toFixed(0) : gamma} MHz`;
  }
  for (let px = 0; px <= W && !sub; px++) {
    const o = -span + 2 * span * (px - pad) / (W - 2 * pad);
    let v = 0;
    for (let i = 0; i < comps.o.length; i++) {
      const x = (o - comps.o[i]) / sigma;
      if (Math.abs(x) < 5) v += comps.s[i] * Math.exp(-0.5 * x * x);
    }
    prof[px] = v; pmax = Math.max(pmax, v);
  }
  if (pmax > 0) {
    // a harmonic signal swings both ways: centred on the middle of the box
    const mid = sub && state.sub.harmonic ? (base + top) / 2 : base, amp = mid === base ? base - top : (base - top) / 2;
    g.beginPath();
    for (let px = 0; px <= W; px++) g.lineTo(px, mid - amp * prof[px] / pmax);
    g.strokeStyle = ink(sub ? "--trace" : "--ink-3"); g.lineWidth = sub ? 1.2 : 1.5; g.stroke();
  }
  for (let i = 0; i < comps.o.length; i++) {
    const x = xOf(comps.o[i]), h = (base - top) * comps.s[i] / smax;
    g.strokeStyle = comps.l[i] ? ink("--trace") : ink("--rule");
    g.lineWidth = comps.l[i] ? 2 : 1;
    g.beginPath(); g.moveTo(x, base); g.lineTo(x, base - h); g.stroke();
  }
  g.strokeStyle = ink("--rule"); g.lineWidth = 1;
  g.beginPath(); g.moveTo(pad, base); g.lineTo(W - pad, base); g.stroke();

  g.fillStyle = ink("--ink-3"); g.font = '10px "IBM Plex Mono", ui-monospace, monospace';
  g.textAlign = "left"; g.fillText(`−${(span / 1000).toFixed(2)} GHz`, pad, H - 6);
  g.textAlign = "right"; g.fillText(`+${(span / 1000).toFixed(2)} GHz`, W - pad, H - 6);
  g.textAlign = "center"; g.fillText(subNote || `Doppler FWHM ${fwhm.toFixed(0)} MHz`, W / 2, H - 6);
}

/** The published 2008 model (Salumbides et al.) was fitted to these levels; beyond them it extrapolates. */
const inside2008 = l => l.vu <= 43 && l.vl <= 17;
/** A signed frequency difference in MHz, or GHz once it is large. */
const signedMHz = v => {
  const a = Math.abs(v), s = v >= 0 ? "+" : "−";
  return a >= 1e4 ? `${s}${fmt(a / 1000, a >= 1e6 ? 0 : 1)} GHz` : `${s}${fmt(a, a < 10 ? 2 : 1)} MHz`;
};
function published2008Row(line) {
  if (line.d08 === null || line.d08 === undefined) return "";
  const where = inside2008(line) ? "inside its fitted range" : "outside its fitted range (v′ ≤ 43, v″ ≤ 17): extrapolated";
  return `<dt>2008 model</dt><dd>${signedMHz(line.d08)} from this model <span class="n">· ${where}</span></dd>`;
}

/** The selected line seen from the laser: its fundamental equivalent, and how far the harmonic is from it. */
function laserRows(line) {
  const {n, nu} = state.laser, d = (line.nu - nu) * P.MHZ_PER_CM;
  const signed = v => (v >= 0 ? "+" : "−") + fmt(Math.abs(v), Math.abs(v) < 100 ? 2 : 1);
  return `<dt>÷ ${n}</dt><dd>${P.fundamentalOf(line.nu, n, "nm").toFixed(6)} nm · ${P.fundamentalOf(line.nu, n, "THz").toFixed(7)} THz</dd>
    <dt>− laser</dt><dd>${signed(d)} MHz (${signed(d / n)} MHz at the fundamental)</dd>`;
}

async function select(line) {
  state.selected = line;
  state.selectedComps = null;
  $("#csvhfs").disabled = true;
  renderTable(); drawPlot(); pushURL();
  const {value, why, flags} = unc(line);
  const S = P.strengthAt(line, partitionAt(line.iso, state.T), state.T);
  const lam = 1e7 / line.nu;
  const isoLabel = state.manifest.isotopologues[line.iso].label;
  $("#detail").innerHTML = `
    <h2>${isoLabel} ${label(line)}</h2>
    <div class="chips">${flags.map(([t, c]) => `<span class="chip ${c}">${t}</span>`).join("")}
      <span class="chip kind k-${kindOf(line)}">${line.measured ? "measured: " + KIND_LABEL[kindOf(line)] : "model only"}</span></div>
    <dl>
      <dt>λ vac</dt><dd>${lam.toFixed(6)} nm</dd>
      <dt>λ air</dt><dd>${P.vacToAir(lam).toFixed(6)} nm</dd>
      <dt>ν</dt><dd>${line.nu.toFixed(6)} cm⁻¹</dd>
      <dt>f</dt><dd>${fmt(line.nu * P.MHZ_PER_CM, 1)} MHz</dd>
      <dt>u (1σ)</dt><dd>± ${fmtU(value)} MHz</dd>
      <dt>S</dt><dd>${line.s0 > 0 ? `${S.toExponential(3)} cm at ${state.T} K` : "not computed: no thermal absorption (measured in emission)"}</dd>
      <dt>E″</dt><dd>${line.el.toFixed(2)} cm⁻¹</dd>
      <dt>Doppler</dt><dd>${(P.dopplerFWHM(line.nu, state.T, line.iso) * P.MHZ_PER_CM).toFixed(1)} MHz FWHM</dd>
      ${published2008Row(line)}
      ${state.laser ? laserRows(line) : ""}
    </dl>
    <p class="note">${sentence(why)}</p>
    ${measurementsHTML(line)}
    <div id="hfsbox"><p class="note">loading hyperfine structure…</p></div>`;

  const patterns = await getJSON(state.manifest.isotopologues[line.iso].hfs, state.hfs, line.iso);
  if (state.selected !== line) return;                       // the user moved on while we fetched
  let comps = patterns[`${line.vu}-${line.vl}${line.branch}${line.J}`], approx = false;
  if (!comps && line.hfsFile) {
    // not among the lines with the full calculation: use its ΔJ = 0 pattern from the shard file
    try { await getJSON(line.hfsFile, state.hfsShards); } catch {}
    if (state.selected !== line) return;
    const h = hfsOf(line);
    if (h) { comps = {o: h.o, s: h.s.map(v => v / 1000), l: h.o.map((_, i) => `a${i + 1}`), x: h.x}; approx = true; }
  }
  const box = $("#hfsbox");
  if (!comps) {
    box.innerHTML = `<h3>Hyperfine structure</h3><p class="note">Precomputed for the strongest
      lines of each isotopologue, for the strongest of every part of the spectrum, and for every
      line we hold a measurement of. This one
      is not among them; <code>i2spec line "${label(line)}"${line.iso === "127I2" ? "" : " -i " + line.iso}</code>
      computes any line on demand.</p>`;
    return;
  }
  state.selectedComps = {comps, approx};
  $("#csvhfs").disabled = false;
  box.innerHTML = `<h3>Hyperfine structure</h3><canvas id="hfs"></canvas><div id="hfstable"></div>`;
  drawHyperfine($("#hfs"), comps, line);
  const main = comps.l.map((l, i) => [l, comps.o[i], comps.s[i]]).filter(r => r[0]);
  const f0 = line.nu * P.MHZ_PER_CM;
  $("#hfstable").innerHTML =
    `<table class="sub"><thead><tr><th>comp</th><th>offset (MHz)</th><th>f (MHz)</th><th>rel. strength</th></tr></thead><tbody>` +
    main.map(([l, o, s]) => `<tr><td>${l}</td><td>${o.toFixed(3)}</td><td>${fmt(f0 + o, 1)}</td><td>${s.toFixed(4)}</td></tr>`).join("") +
    `</tbody></table><p class="note">${main.length} main components (ΔF = ΔJ)${
      comps.o.length > main.length ? `, plus ${comps.o.length - main.length} weaker ones drawn faint above` : ""}.
     ${approx ? "Computed without the J ± 2 couplings (ΔJ = 0), which changes the offsets by less than 1 MHz; " +
       `<code>i2spec line "${label(line)}"</code> gives the full calculation.`
       : "The ΔF ≠ ΔJ components are omitted below 2% of the peak, where nothing in a spectrum shows them."}
     Offsets are from the hyperfine-free centre, which is not itself an observable line.</p>`;
}

function measurementsHTML(line) {
  if (!line.meas) return "";
  const f0 = line.nu * P.MHZ_PER_CM;
  const rows = [];
  const signed = r => r === null ? "–" : (r >= 0 ? "+" : "−") + Math.abs(r).toFixed(Math.abs(r) < 1 ? 3 : 1);
  for (const e of line.meas) {
    // one heading row per source, carrying its link, then its values
    rows.push(`<tr class="srcrow"><td colspan="5">${sourceLink(e.s)}${e.i ?
      ` <span class="n">· ${e.i} hyperfine interval${e.i === 1 ? "" : "s"}</span>` : ""}</td></tr>`);
    for (const [comp, value, u, res, res08] of e.f) {
      rows.push(`<tr><td>${comp || "centre"}</td><td>${fmt(value, u < 0.1 ? 4 : u < 10 ? 2 : 1)}</td>
        <td>${u < 0.1 ? u.toFixed(4) : u.toFixed(u < 10 ? 2 : 1)}</td><td>${signed(res)}</td>
        <td>${res08 === undefined || res08 === null ? "–" : Math.abs(res08) >= 1e4 ? signedMHz(res08).replace(" GHz", "&nbsp;G") : signed(res08)}</td></tr>`);
    }
  }
  const shown = rows.slice(0, 60);
  return `<h3>Measurements</h3>
    <div class="measwrap"><table class="sub meas"><thead><tr><th>comp</th><th>f (MHz)</th><th>u (MHz)</th><th>− model</th><th>− 2008</th></tr></thead>
    <tbody>${shown.join("")}</tbody></table></div>
    <p class="note">${rows.length > shown.length ? `First ${shown.length} of ${rows.length} rows. ` : ""}Each source links to
    its DOI; the <a href="docs/references.html">references</a> page lists them all. Frequencies are as published
    (atlas lines are hyperfine-free centres on the atlas's own scale; the Orsay scales are corrected by +23.85 ppb, part I, and +200.8 ppb, Partie IV);
    “− model” is the measured value minus this model's centre or component (MHz), ${fmt(f0, 1)} MHz for the centre;
    “− 2008” minus the published model of Salumbides <i>et al.</i> (2008), computed here from its potentials and
    the published hyperfine formulae (MHz; G: GHz).</p>`;
}

// --- CSV export ------------------------------------------------------------------------------------

const csvCell = v => {
  const t = v === null || v === undefined ? "" : String(v);
  return /[",\n\r]/.test(t) ? `"${t.replace(/"/g, '""')}"` : t;
};
const csvRow = cells => cells.map(csvCell).join(",");

/** Comment lines that say where a file came from: model, parameters, export, download date, T and view. */
function csvHeader(what) {
  const m = state.manifest, unit = $("#unit").value, d = P.DIGITS[unit];
  const a = P.fromWavenumber(state.lo, unit), b = P.fromWavenumber(state.hi, unit);
  const out = [
    `# ${what}, from the i2spec line explorer`,
    `# model: ${m.model.package} ${m.model.version}` +
      `${m.model.git ? ", git " + m.model.git : ""}; data exported ${m.generated}`,
    `# downloaded ${new Date().toISOString()}; isotopologue ${state.iso}; T = ${state.T} K; ` +
      `view ${Math.min(a, b).toFixed(d)} to ${Math.max(a, b).toFixed(d)} ${unit}; shown: ${state.show || "none"}`,
  ];
  if (state.laser) out.push(`# ${laserText()}`);
  out.push(`# link: ${location.href}`);
  return out;
}

function download(name, lines) {
  const blob = new Blob([lines.join("\n") + "\n"], {type: "text/csv;charset=utf-8"});
  const a = document.createElement("a");
  a.href = URL.createObjectURL(blob);
  a.download = name;
  document.body.append(a);
  a.click();
  a.remove();
  setTimeout(() => URL.revokeObjectURL(a.href), 1000);
}
const fileStem = () => {
  const unit = $("#unit").value, d = Math.max(0, P.DIGITS[unit] - 2);
  const a = P.fromWavenumber(state.lo, unit), b = P.fromWavenumber(state.hi, unit);
  return `i2spec_${state.iso}_${Math.min(a, b).toFixed(d)}-${Math.max(a, b).toFixed(d)}${unit.replace("-", "")}`;
};

/**
 * Every line in view that the show toggles keep (not only the 500 the table lists), in the table's order.
 * The ΔJ = 0 hyperfine offsets are added when their files are loaded, which they are below 0.1 nm.
 */
function exportLines() {
  const pool = [...shownLines()].sort(tableOrder());
  const q = partitionAt(state.iso, state.T), laser = state.laser;
  const withHfs = pool.some(l => hfsOf(l));
  const head = ["line", "isotopologue", "lambda_vac_nm", "lambda_air_nm", "wavenumber_cm-1", "frequency_MHz",
    `S_cm_at_${state.T}K`, "E_lower_cm-1", "u_MHz_1sigma", "class", "measured_by", "published_2008_minus_model_MHz"];
  if (laser) head.push(`fundamental_lambda_vac_nm_div${laser.n}`, `fundamental_MHz_div${laser.n}`, "minus_laser_harmonic_MHz");
  if (withHfs) head.push("hfs_offsets_MHz_dJ0", "hfs_strengths_permille_dJ0");
  const rows = pool.map(l => {
    const lam = 1e7 / l.nu;
    const cells = [label(l), l.iso, lam.toFixed(7), P.vacToAir(lam).toFixed(7), l.nu.toFixed(6),
      (l.nu * P.MHZ_PER_CM).toFixed(2), P.strengthAt(l, q, state.T).toExponential(4), l.el.toFixed(2), l.u,
      KIND_LABEL[kindOf(l)], l.meas ? [...new Set(l.meas.map(e => source(e.s).short))].join("; ") : "",
      l.d08 === null || l.d08 === undefined ? "" : l.d08];
    if (laser) cells.push(P.fundamentalOf(l.nu, laser.n, "nm").toFixed(7), (l.nu * P.MHZ_PER_CM / laser.n).toFixed(2),
      ((l.nu - laser.nu) * P.MHZ_PER_CM).toFixed(2));
    if (withHfs) {
      const h = hfsOf(l);
      cells.push(h ? h.o.join(" ") : "", h ? h.s.join(" ") : "");
    }
    return csvRow(cells);
  });
  download(fileStem() + ".csv", [
    ...csvHeader(`${fmt(pool.length)} B-X lines`),
    "# S is the line strength at the cell temperature; u is the model's 1-sigma position uncertainty;" +
      " class is how the line is known (precision, atlas, model only); published_2008_minus_model is the" +
      " position of the published model of Salumbides et al. (2008), computed from its potentials, minus this model's",
    ...(withHfs ? ["# hfs_*: main hyperfine components without the J +- 2 couplings (offsets from the centre, MHz;" +
      " strengths per mille), space separated"] : []),
    csvRow(head), ...rows]);
  setStatus(`saved ${fmt(pool.length)} line${pool.length === 1 ? "" : "s"} to ${fileStem()}.csv`);
}

/** The hyperfine components of the selected line, as the detail pane shows them. */
function exportHyperfine() {
  const line = state.selected, sc = state.selectedComps;
  if (!line || !sc) return;
  const {comps, approx} = sc, f0 = line.nu * P.MHZ_PER_CM;
  const head = ["component", "offset_MHz", "frequency_MHz", "lambda_vac_nm", "wavenumber_cm-1", "relative_strength"];
  if (state.laser) head.push(`fundamental_MHz_div${state.laser.n}`);
  const rows = comps.o.map((o, i) => {
    const f = f0 + o, nu = f / P.MHZ_PER_CM;
    const cells = [comps.l[i] || "", o.toFixed(3), f.toFixed(3), (1e7 / nu).toFixed(8), nu.toFixed(7), comps.s[i]];
    if (state.laser) cells.push((f / state.laser.n).toFixed(3));
    return csvRow(cells);
  });
  const name = `i2spec_${line.iso}_${label(line).replace(/[()\s]+/g, "_").replace(/_+$/, "")}_hyperfine.csv`;
  download(name, [
    ...csvHeader(`hyperfine components of ${state.manifest.isotopologues[line.iso].label} ${label(line)}`),
    `# centre (hyperfine-free, not itself observable): ${f0.toFixed(3)} MHz = ${line.nu.toFixed(7)} cm-1; ` +
      (approx ? "computed without the J +- 2 couplings (offsets good to about 1 MHz)"
              : "full calculation; components other than the main ones (ΔF = ΔJ) below 2% of the peak omitted"),
    csvRow(head), ...rows]);
  setStatus(`saved ${comps.o.length} components to ${name}`);
}

// --- driving it ----------------------------------------------------------------------------------

function syncInputs() {
  const unit = $("#unit").value, d = P.DIGITS[unit];
  const a = P.fromWavenumber(state.lo, unit), b = P.fromWavenumber(state.hi, unit);
  $("#lo").value = Math.min(a, b).toFixed(d);
  $("#hi").value = Math.max(a, b).toFixed(d);
  syncSummary();
}

/** One line under the header on small screens, standing in for the hidden options. */
function syncSummary() {
  if (!state.manifest) return;
  const unit = $("#unit").value, d = P.DIGITS[unit];
  const a = P.fromWavenumber(state.lo, unit), b = P.fromWavenumber(state.hi, unit);
  const u = unit === "cm-1" ? "cm⁻¹" : unit === "nm-air" ? "nm (air)" : unit;
  $("#summary").textContent = `${Math.min(a, b).toFixed(Math.max(0, d - 2))}–${Math.max(a, b).toFixed(Math.max(0, d - 2))} ${u}` +
    ` · ${state.manifest.isotopologues[state.iso].label} · ${state.T} K` +
    (state.show !== "all" ? ` · ${state.show ? state.show.replace(",", " + ") : "nothing"} shown` : "");
}

//: Never zoom tighter than this (cm-1), about 0.3 MHz: below it the view is empty of structure.
const MIN_SPAN = 1e-5;

let tableTimer = null;

/**
 * Clamp a view to the range the export covers, widening rather than sliding: truncating one edge made
 * zooming out past the limit look like a pan.
 */
function clampView(nuLo, nuHi) {
  const info = state.manifest.isotopologues[$("#iso").value];
  const LO = info.nu_min, HI = info.nu_max;
  if (nuHi - nuLo < MIN_SPAN) {
    const mid = 0.5 * (nuLo + nuHi);
    nuLo = mid - MIN_SPAN / 2; nuHi = mid + MIN_SPAN / 2;
  }
  if (nuHi - nuLo >= HI - LO) return [LO, HI];
  if (nuLo < LO) { nuHi += LO - nuLo; nuLo = LO; }
  if (nuHi > HI) { nuLo -= nuHi - HI; nuHi = HI; }
  return [Math.max(nuLo, LO), Math.min(nuHi, HI)];
}

async function show(nuLo, nuHi, {keepInputs = false, deferTable = false} = {}) {
  if (!(nuLo < nuHi)) return;
  state.iso = $("#iso").value;
  const info = state.manifest.isotopologues[state.iso];
  [nuLo, nuHi] = clampView(nuLo, nuHi);
  state.lo = nuLo; state.hi = nuHi;
  state.T = Number($("#temp").value);
  state.unitFlips = $("#unit").value.startsWith("nm");
  if (!keepInputs) syncInputs();
  const token = ++state.busy;
  setStatus("loading…");
  try {
    await getJSON(state.manifest.isotopologues[state.iso].partition, state.partitions, state.iso);
    const lines = await linesIn(state.iso, state.lo, state.hi, {loadedOnly: deferTable});
    const padded = await linesIn(state.iso, state.lo - PLOT_PAD_CM, state.hi + PLOT_PAD_CM, {loadedOnly: deferTable});
    // a wide view: fetch the rest of the list in the background, so zooming out finds it loaded
    if (state.hi - state.lo > 0.2 * (info.nu_max - info.nu_min)) loadInBackground(state.iso, info.shards);
    if (token !== state.busy) return;                    // a newer view won
    if (hyperfineShown()) await loadHfs(state.iso, state.lo - PLOT_PAD_CM, state.hi + PLOT_PAD_CM);
    if (token !== state.busy) return;
    state.lines = lines;
    state.plotLines = padded;
    if (state.selected && (state.selected.iso !== state.iso ||
        state.selected.nu < state.lo || state.selected.nu > state.hi)) {
      state.selected = null;
      $("#detail").innerHTML = EMPTY_DETAIL;
    }
    // while the wheel is turning, redraw only the spectrum; rebuild the table once it stops
    clearTimeout(tableTimer);
    syncSummary();
    if (deferTable) { drawPlot(); tableTimer = setTimeout(() => { renderTable(); pushURL(); }, 160); }
    else { renderTable(); drawPlot(); pushURL(); }
  } catch (err) {
    setStatus("could not load the line data: " + err.message);
  }
}

const EMPTY_DETAIL = `<h2>No line selected</h2>
  <p class="note">Click a row in the table, or click near a line in the spectrum, to see its
  hyperfine components, its strength and how far the model can be trusted there.</p>`;

function readRange() {
  const unit = $("#unit").value;
  let lo = P.toWavenumber(Number($("#lo").value), unit);
  let hi = P.toWavenumber(Number($("#hi").value), unit);
  if (!Number.isFinite(lo) || !Number.isFinite(hi)) return null;
  return lo < hi ? [lo, hi] : [hi, lo];
}

async function refresh() {
  const range = readRange();
  if (range) await show(range[0], range[1], {keepInputs: true});
}

async function findLine(text, {keepView = false} = {}) {
  const m = text.trim().match(/^([PR])\s*\(?\s*(\d+)\s*\)?\s+(\d+)\s*-\s*(\d+)$/i);
  if (!m) { setStatus(`"${text}" is not a line label; write it like R(56) 32-0`); return; }
  const key = `${m[1].toUpperCase()}(${m[2]}) ${m[3]}-${m[4]}`;
  // a shared link names both a view and a line in it: select it without loading the whole list
  const here = keepView && state.iso === $("#iso").value && state.lines.find(l => label(l) === key);
  if (here) { await select(here); return; }
  const index = await buildIndex($("#iso").value);
  const line = index.get(key);
  if (!line) { setStatus(`${key} is not in the exported list (it is weaker than the cutoff and has not been measured)`); return; }
  const w = P.dopplerFWHM(line.nu, state.T, line.iso) * 12;
  await show(line.nu - w, line.nu + w);
  const same = state.lines.find(l => l.nu === line.nu) || line;
  await select(same);
}

/**
 * Ask for a new view without waiting for the last one. A trackpad sends wheel events far faster than
 * the view can be redrawn; each request only updates the target, and one redraw per frame catches up
 * with the latest target, so a burst of events costs a few redraws, not one each.
 */
let target = null, previewPending = false, renderTimer = null, rendering = false, lastReal = 0;
function requestView(lo, hi) {
  target = clampView(lo, hi);
  clearTimeout(tableTimer);                // the table waits until the gesture is over
  // every frame: stretch the last finished picture to the new range, which costs almost nothing ...
  if (!previewPending) {
    previewPending = true;
    requestAnimationFrame(() => {
      previewPending = false;
      if (target && !drawPreview(target[0], target[1])) scheduleReal(0);
    });
  }
  // ... and redraw properly when the gesture pauses, or every 500 ms while it goes on
  scheduleReal(performance.now() - lastReal > 500 ? 0 : 110);
}
function scheduleReal(delay) {
  clearTimeout(renderTimer);
  renderTimer = setTimeout(realRender, delay);
}
async function realRender() {
  if (rendering) { scheduleReal(40); return; }
  if (!target) return;
  const [a, b] = target;
  target = null;
  rendering = true;
  lastReal = performance.now();
  try { await show(a, b, {deferTable: true}); } finally { rendering = false; }
}
/** Redraw the current view properly, e.g. once data for it has arrived. */
function refreshView() { if (!target) target = [state.lo, state.hi]; scheduleReal(0); }
/** The view the next redraw will show: the pending target if there is one. */
const viewNow = () => target || [state.lo, state.hi];

function closeOptions() {
  if (!$("#controls").classList.contains("open")) return;
  $("#controls").classList.remove("open");
  $("#optbtn").setAttribute("aria-expanded", "false");
  $("#optbtn").textContent = "Options ▾";
}

function nuAtClientX(clientX) {
  const {m, pw, flip} = state.geom;
  const rect = $("#plot").getBoundingClientRect();
  let f = Math.min(Math.max((clientX - rect.left - m.l) / pw, 0), 1);
  if (flip) f = 1 - f;
  return state.lo + (state.hi - state.lo) * f;
}

function wire() {
  $("#go").onclick = () => { refresh(); closeOptions(); };
  // on a phone the detail pane is a bottom sheet: its heading opens and closes it
  $("#detail").addEventListener("click", e => {
    if (e.target.closest("h2") && matchMedia("(max-width: 880px)").matches) $("#detail").classList.toggle("expanded");
  });
  $("#optbtn").onclick = () => {
    const open = !$("#controls").classList.contains("open");
    $("#controls").classList.toggle("open", open);
    $("#optbtn").setAttribute("aria-expanded", String(open));
    $("#optbtn").textContent = open ? "Close ▴" : "Options ▾";
  };
  for (const id of ["#lo", "#hi"]) $(id).addEventListener("keydown", e => { if (e.key === "Enter") refresh(); });
  $("#unit").onchange = () => { syncInputs(); state.unitFlips = $("#unit").value.startsWith("nm"); drawPlot(); };
  $("#iso").onchange = () => show(state.lo, state.hi);
  $("#temp").onchange = () => show(state.lo, state.hi);
  for (const id of ["#cold", "#path"]) $(id).oninput = drawPlot;
  $("#find").addEventListener("keydown", e => { if (e.key === "Enter") findLine(e.target.value); });
  $$("#show button").forEach(b => b.onclick = () => {
    const k = b.dataset.kind;
    state.kinds.has(k) ? state.kinds.delete(k) : state.kinds.add(k);
    b.setAttribute("aria-pressed", String(state.kinds.has(k)));
    state.show = showKey();
    renderTable(); drawPlot(); pushURL();
  });

  // scrolling on past either end of the table moves the window along the spectrum
  const wrap = $(".tablewrap");
  let over = 0;
  wrap.addEventListener("wheel", e => {
    if (!byWavelength() || state.panning || !state.shown || !state.shown.length) return;
    const atBottom = wrap.scrollTop + wrap.clientHeight >= wrap.scrollHeight - 2;
    const atTop = wrap.scrollTop <= 0;
    if ((e.deltaY > 0 && atBottom) || (e.deltaY < 0 && atTop)) {
      over += Math.abs(e.deltaY);
      if (over > 160) { over = 0; tablePan(e.deltaY > 0 ? state.sort.dir : -state.sort.dir); }
    } else over = 0;
  }, {passive: true});

  const buttons = PRESETS.map(([name, why, lo, hi]) => {
    const b = document.createElement("button");
    b.textContent = name; b.title = why;
    b.onclick = () => { show(lo, hi); closeOptions(); };
    return b;
  });
  const all = document.createElement("button");
  all.textContent = "all";
  all.title = "the whole exported range; the first time this loads about 10 MB of line data";
  all.onclick = () => {
    const i = state.manifest.isotopologues[$("#iso").value];
    show(i.nu_min, i.nu_max);
  };
  buttons.push(all);
  $("#presets").replaceChildren(...buttons);

  $$("#mode button").forEach(b => b.onclick = () => setMode(b.dataset.mode));
  const readSub = () => {
    const g = Number($("#gamma").value), m = Number($("#mod").value);
    if (g > 0) state.sub.fwhm = g;
    if (m > 0) state.sub.modulation = m;
    state.sub.harmonic = Number($("#det").value);
    $$(".modonly").forEach(el => el.hidden = state.mode !== "sub" || !state.sub.harmonic);
    redrawAll(); pushURL();
  };
  for (const id of ["#gamma", "#mod"]) $(id).oninput = readSub;
  $("#det").onchange = readSub;

  $("#laser").addEventListener("keydown", e => { if (e.key === "Enter") { applyLaser(); closeOptions(); } });
  $("#lasergo").onclick = () => { applyLaser(); closeOptions(); };
  for (const id of ["#lunit", "#harm"]) $(id).onchange = () => { if ($("#laser").value.trim()) applyLaser(); };
  $("#csv").onclick = exportLines;
  $("#csvhfs").onclick = exportHyperfine;

  $("#theme").onclick = () => {
    const dark = matchMedia("(prefers-color-scheme: dark)").matches;
    const now = document.documentElement.dataset.theme || (dark ? "dark" : "light");
    document.documentElement.dataset.theme = now === "dark" ? "light" : "dark";
    try { localStorage.setItem("theme", document.documentElement.dataset.theme); } catch {}
    redrawAll();
  };

  const plot = $("#plot");
  plot.addEventListener("mousemove", e => {
    if (!state.geom) return;
    const nu = nuAtClientX(e.clientX), unit = $("#unit").value;
    const {m, pw, ph, y} = state.geom;
    const k = Math.round((e.clientX - plot.getBoundingClientRect().left - m.l) / pw * y.length);
    const ro = $("#readout");
    ro.textContent = `${P.fromWavenumber(nu, unit).toFixed(P.DIGITS[unit])} ${unit}   ${nu.toFixed(5)} cm⁻¹`;
    ro.style.top = (m.t + 6) + "px";      // inside the plot frame, clear of the top axis and its labels
    ro.style.right = (m.r + 6) + "px";
    const cross = $("#crosshair");
    const x = e.clientX - plot.getBoundingClientRect().left;
    cross.style.left = x + "px";
    cross.style.top = m.t + "px";
    cross.style.height = ph + "px";
    cross.hidden = x < m.l || x > m.l + pw;
  });
  plot.addEventListener("mouseleave", () => { $("#crosshair").hidden = true; $("#readout").textContent = ""; });
  plot.addEventListener("click", e => {
    const pool = shownLines();
    if (!pool.length || dragged) return;
    const nu = nuAtClientX(e.clientX);
    let best = pool[0];
    for (const l of pool) if (Math.abs(l.nu - nu) < Math.abs(best.nu - nu)) best = l;
    select(best);
  });
  plot.addEventListener("wheel", e => {
    e.preventDefault();
    const [lo, hi] = viewNow();
    if (Math.abs(e.deltaX) > Math.abs(e.deltaY)) {
      // a sideways scroll (trackpad, tilt wheel, shift-wheel) pans; the spectrum moves with the fingers
      const dnu = e.deltaX / (state.geom ? state.geom.pw : 800) * (hi - lo) * (state.geom && state.geom.flip ? -1 : 1);
      requestView(lo + dnu, hi + dnu);
      return;
    }
    const frac = state.geom ? Math.min(Math.max((e.clientX - $("#plot").getBoundingClientRect().left - state.geom.m.l) / state.geom.pw, 0), 1) : 0.5;
    const nu = lo + (hi - lo) * (state.geom && state.geom.flip ? 1 - frac : frac);
    // scale the step with the wheel delta, so a trackpad's many small events zoom smoothly
    const f = Math.pow(1.3, Math.max(-3, Math.min(3, e.deltaY / 100)));
    requestView(nu - (nu - lo) * f, nu + (hi - nu) * f);
  }, {passive: false});

  // Dragging selects a range on the wavenumber axis and zooms to it, the way a spectrum viewer
  // should; panning is on [ and ]. A drag shorter than a few pixels is treated as a click.
  let drag = null, dragged = false;
  const sel = $("#selection");

  function drawSelection(x0, x1) {
    const {m, pw, ph} = state.geom;
    const a = Math.min(Math.max(Math.min(x0, x1), m.l), m.l + pw);
    const b = Math.min(Math.max(Math.max(x0, x1), m.l), m.l + pw);
    sel.style.left = a + "px";
    sel.style.width = (b - a) + "px";
    sel.style.top = m.t + "px";
    sel.style.height = ph + "px";
    const unit = $("#unit").value;
    const [lo, hi] = [nuAtClientX(a + plot.getBoundingClientRect().left),
                      nuAtClientX(b + plot.getBoundingClientRect().left)].sort((p, q) => p - q);
    const span = Math.abs(P.fromWavenumber(hi, unit) - P.fromWavenumber(lo, unit));
    sel.dataset.span = `${span.toPrecision(3)} ${unit === "cm-1" ? "cm⁻¹" : unit}`;
    sel.hidden = false;
  }

  // Touch: one finger pans, two fingers pinch-zoom about the point between them, a tap picks a line.
  // Each gesture works from the view it started with, so the spectrum stays under the fingers
  // however many redraws the view falls behind by.
  const touches = new Map();                     // pointerId -> clientX
  let gesture = null;
  let plotLeft = 0;                              // the plot's position, read once per gesture: reading
  const plotFrac = clientX => {                  // it on every move forces a layout each time
    const {m, pw} = state.geom;
    return (clientX - plotLeft - m.l) / pw;
  };
  const nuAtFrac = (f, lo, hi) => lo + (hi - lo) * (state.geom.flip ? 1 - f : f);
  function startGesture() {
    plotLeft = plot.getBoundingClientRect().left;
    const [lo, hi] = viewNow();
    const xs = [...touches.values()];
    if (xs.length >= 2) {
      const [a, b] = xs;
      const f = plotFrac((a + b) / 2);
      gesture = {kind: "pinch", lo, hi, dist: Math.max(Math.abs(a - b), 20), nu: nuAtFrac(f, lo, hi)};
    } else if (xs.length === 1) {
      gesture = {kind: "pan", lo, hi, x: xs[0]};
    } else gesture = null;
  }

  plot.addEventListener("pointerdown", e => {
    if (!state.geom) return;
    if (e.pointerType === "touch") {
      plot.setPointerCapture(e.pointerId);
      touches.set(e.pointerId, e.clientX);
      if (touches.size === 1) dragged = false;
      startGesture();
      return;
    }
    // mouse: a drag pans (the spectrum follows the pointer); a shift-drag selects a range to zoom to
    const x = e.clientX - plot.getBoundingClientRect().left;
    const [lo, hi] = viewNow();
    drag = {x0: x, x1: x, select: e.shiftKey, lo, hi};
    dragged = false;
    plot.setPointerCapture(e.pointerId);
    if (!e.shiftKey) plot.classList.add("grabbing");
  });
  plot.addEventListener("pointermove", e => {
    if (e.pointerType === "touch") {
      if (!touches.has(e.pointerId) || !gesture || !state.geom) return;
      touches.set(e.pointerId, e.clientX);
      const {lo, hi} = gesture, w = hi - lo;
      if (gesture.kind === "pan") {
        const dx = e.clientX - gesture.x;
        if (Math.abs(dx) > 6) dragged = true;           // a moved finger is not a tap
        if (!dragged) return;
        // the spectrum follows the finger: moving right shows what was to the left
        const dnu = dx / state.geom.pw * w * (state.geom.flip ? 1 : -1);
        requestView(lo + dnu, hi + dnu);
      } else {
        const [a, b] = [...touches.values()];
        dragged = true;
        const scale = gesture.dist / Math.max(Math.abs(a - b), 20);   // fingers apart: zoom in
        const nw = w * scale, f = plotFrac((a + b) / 2);
        // the wavenumber that was under the midpoint stays under the (possibly moved) midpoint
        const g = state.geom.flip ? 1 - f : f;
        requestView(gesture.nu - g * nw, gesture.nu + (1 - g) * nw);
      }
      return;
    }
    if (!drag || !state.geom) return;
    drag.x1 = e.clientX - plot.getBoundingClientRect().left;
    if (Math.abs(drag.x1 - drag.x0) < 4 && !dragged) return;
    dragged = true;
    if (drag.select) { drawSelection(drag.x0, drag.x1); return; }
    const w = drag.hi - drag.lo;
    const dnu = (drag.x1 - drag.x0) / state.geom.pw * w * (state.geom.flip ? 1 : -1);
    requestView(drag.lo + dnu, drag.hi + dnu);
  });
  const endTouch = e => {
    touches.delete(e.pointerId);
    startGesture();                                // a lifted second finger turns a pinch back into a pan
  };
  plot.addEventListener("pointerup", e => {
    try { plot.releasePointerCapture(e.pointerId); } catch {}
    if (e.pointerType === "touch") { endTouch(e); return; }
    sel.hidden = true;
    plot.classList.remove("grabbing");
    if (!drag || !state.geom) { drag = null; return; }
    const {x0, x1, select: selecting} = drag;
    drag = null;
    if (Math.abs(x1 - x0) < 4) return;           // a click: the click handler picks the line
    if (!selecting) return;                      // a pan: its redraw is already scheduled
    const rect = plot.getBoundingClientRect();
    const a = nuAtClientX(x0 + rect.left), b = nuAtClientX(x1 + rect.left);
    const [lo, hi] = [Math.min(a, b), Math.max(a, b)];
    if (hi - lo > 1e-9) show(lo, hi);
  });
  plot.addEventListener("pointercancel", e => {
    if (e.pointerType === "touch") { endTouch(e); return; }
    drag = null; sel.hidden = true; plot.classList.remove("grabbing");
  });

  addEventListener("keydown", e => {
    if (e.target.matches("input, select")) return;
    if (e.key === "/") { e.preventDefault(); $("#find").focus(); }
    if (e.key === "l" || e.key === "L") { e.preventDefault(); $("#laser").focus(); $("#laser").select(); }
    if ((e.key === "ArrowDown" || e.key === "ArrowUp") && state.shown && state.shown.length) {
      // not already inside the table: jump into it rather than scrolling the page
      e.preventDefault();
      state.keyNav = true;
      select(state.selected && state.shown.includes(state.selected) ? state.selected : state.shown[0]);
    }
    if (e.key === "[" || e.key === "]" || e.key === "ArrowLeft" || e.key === "ArrowRight") {
      e.preventDefault();
      const [lo, hi] = viewNow();
      // ← and [ move to the left of the plot as drawn; with a wavelength axis that is higher wavenumber
      const right = e.key === "]" || e.key === "ArrowRight";
      const w = (hi - lo) * 0.25 * (right ? 1 : -1) * (state.unitFlips ? -1 : 1);
      requestView(lo + w, hi + w);
    }
  });

  let pending;
  const later = () => { clearTimeout(pending); pending = setTimeout(redrawAll, 80); };
  addEventListener("resize", later);
  // the plot also changes size without the window doing so, e.g. as the options panel opens and closes
  new ResizeObserver(later).observe($("#plotwrap"));
  matchMedia("(prefers-color-scheme: dark)").addEventListener("change", redrawAll);
}

/** The y axis: "sigma" (cross section), "trans" (cell transmission) or "sub" (sub-Doppler). */
function setMode(mode) {
  state.mode = ["sigma", "trans", "sub"].includes(mode) ? mode : "sub";
  $$("#mode button").forEach(o => o.setAttribute("aria-pressed", String(o.dataset.mode === state.mode)));
  $$(".cellonly").forEach(el => el.hidden = state.mode !== "trans");
  $$(".subonly").forEach(el => el.hidden = state.mode !== "sub");
  $$(".modonly").forEach(el => el.hidden = state.mode !== "sub" || !state.sub.harmonic);
  // sub-Doppler needs the hyperfine patterns, which a wider view has not loaded
  if (!(state.manifest && state.hi > state.lo)) return;      // before the first view: show() draws it
  if (state.mode === "sub" && hyperfineShown()) loadHfs(state.iso, state.lo - PLOT_PAD_CM, state.hi + PLOT_PAD_CM).then(redrawAll);
  else redrawAll();
  pushURL();
}

/** ±0.1 nm around a laser's harmonic; ±0.04 nm for a closer sub-Doppler view. */
const laserHalfSpan = () => state.mode === "sub" ? 0.4 * HFS_SPAN_NM : HFS_SPAN_NM;

/**
 * The laser fields: a fundamental in its own unit and a harmonic order. Jumps to ±0.1 nm around the
 * harmonic unless told not to (a shared link keeps its own view), and marks it on the plot. An empty
 * field forgets the laser.
 */
async function applyLaser({jump = true} = {}) {
  const text = $("#laser").value.trim();
  if (!text) {
    state.laser = null;
    redrawAll(); renderTable(); pushURL();
    if (state.selected) select(state.selected);
    return;
  }
  const value = Number(text), unit = $("#lunit").value, n = Number($("#harm").value);
  if (!Number.isFinite(value) || value <= 0) { setStatus(`"${text}" is not a number; give the fundamental in ${unit}`); return; }
  const nu = P.harmonicWavenumber(value, unit, n);
  state.laser = {value, unit, n, nu};
  const info = state.manifest.isotopologues[$("#iso").value];
  if (nu < info.nu_min || nu > info.nu_max) {
    renderTable(); drawPlot(); pushURL();
    setStatus(`${laserText()}: outside the exported range, ${(1e7 / info.nu_max).toFixed(1)}–${(1e7 / info.nu_min).toFixed(1)} nm`);
    return;
  }
  if (jump) {
    const [lo, hi] = P.harmonicView(value, unit, n, laserHalfSpan());
    await show(lo, hi);
  } else { renderTable(); drawPlot(); pushURL(); }
  if (state.selected) select(state.selected);           // the detail pane shows the line from the laser
}

function redrawAll() {
  drawPlot();
  const c = $("#hfs");
  if (c && state.selected) {
    let comps = (state.hfs.get(state.selected.iso) || {})[
      `${state.selected.vu}-${state.selected.vl}${state.selected.branch}${state.selected.J}`];
    const h = !comps && hfsOf(state.selected);
    if (h) comps = {o: h.o, s: h.s.map(v => v / 1000), l: h.o.map((_, i) => `a${i + 1}`), x: h.x};
    if (comps) drawHyperfine(c, comps, state.selected);
  }
}

async function start() {
  try { const t = localStorage.getItem("theme"); if (t) document.documentElement.dataset.theme = t; } catch {}
  state.manifest = await getJSON("manifest.json", null);
  const m = state.manifest;
  try { if (m.references) state.refs = await getJSON(m.references, null); } catch { state.refs = null; }
  $("#iso").replaceChildren(...Object.entries(m.isotopologues).map(([id, info]) => {
    const o = document.createElement("option");
    o.value = id; o.textContent = info.label;
    return o;
  }));
  $("#iso").value = "127I2";
  $("#prov").innerHTML =
    `${m.model.package} ${m.model.version}${m.model.git ? " · " + m.model.git : ""}` +
    `<br>exported ${m.generated} · lines above ${m.cutoff_cm_at_300K.toExponential(0)} cm at 300 K`;
  $("#prov").title = m.model.note;
  $("#detail").innerHTML = EMPTY_DETAIL;
  wire();

  const q = new URLSearchParams(location.search);
  if (q.get("theme") === "dark" || q.get("theme") === "light") document.documentElement.dataset.theme = q.get("theme");
  if (P.UNITS.includes(q.get("unit"))) $("#unit").value = q.get("unit");
  if (m.isotopologues[q.get("iso")]) $("#iso").value = q.get("iso");
  if (q.get("show")) {
    // "atlas,precision", or the older single filters: measured, precision, atlas
    const v = q.get("show"), legacy = {measured: "atlas,precision"};
    state.kinds = new Set((legacy[v] || v).split(",").filter(k => KINDS.includes(k)));
    state.show = showKey();
    $$("#show button").forEach(o => o.setAttribute("aria-pressed", String(state.kinds.has(o.dataset.kind))));
  }
  if (Number(q.get("T")) >= 200 && Number(q.get("T")) <= 600) $("#temp").value = q.get("T");
  // sub-Doppler settings and the y axis
  if (Number(q.get("gamma")) > 0) { state.sub.fwhm = Number(q.get("gamma")); $("#gamma").value = q.get("gamma"); }
  if (/^[13]f$/.test(q.get("det") || "")) { state.sub.harmonic = Number(q.get("det")[0]); $("#det").value = String(state.sub.harmonic); }
  if (Number(q.get("mod")) > 0) { state.sub.modulation = Number(q.get("mod")); $("#mod").value = q.get("mod"); }
  setMode(q.get("y") || state.mode);
  // a laser: its fields, and its marker; the view is the link's own when it has one
  const laser = Number(q.get("laser"));
  if (laser > 0) {
    $("#laser").value = q.get("laser");
    if (P.UNITS.includes(q.get("lunit"))) $("#lunit").value = q.get("lunit");
    if (["1", "2", "3", "4"].includes(q.get("n"))) $("#harm").value = q.get("n");
    const n = Number($("#harm").value);
    state.laser = {value: laser, unit: $("#lunit").value, n, nu: P.harmonicWavenumber(laser, $("#lunit").value, n)};
  }
  const from = Number(q.get("from")), to = Number(q.get("to"));
  let initialLine = null;
  if (Number.isFinite(from) && Number.isFinite(to) && from !== to) {
    const a = P.toWavenumber(from, $("#unit").value), b = P.toWavenumber(to, $("#unit").value);
    await show(Math.min(a, b), Math.max(a, b));
  } else if (state.laser) {
    await show(...P.harmonicView(state.laser.value, state.laser.unit, state.laser.n, laserHalfSpan()));
  } else {
    // Usually one shard contains the target. Resolve its actual model position
    // without indexing the whole list; other isotopes can use the label index.
    const iso = $("#iso").value, {line, near, halfSpan} = DEFAULT_VIEW;
    const nearby = await linesIn(iso, near - halfSpan, near + halfSpan);
    initialLine = nearby.find(l => label(l) === line) || (await buildIndex(iso)).get(line);
    const centre = initialLine ? initialLine.nu : near;
    await show(centre - halfSpan, centre + halfSpan);
  }
  if (q.get("line")) await findLine(q.get("line"), {keepView: true});
  else if (initialLine) await select(initialLine);
  else if (state.lines.length) {
    // open on the strongest line in view, so the hyperfine pane shows what this tool is for
    await select(state.lines.reduce((a, b) =>
      P.strengthAt(a, partitionAt(a.iso, state.T), state.T) >= P.strengthAt(b, partitionAt(b.iso, state.T), state.T) ? a : b));
  }
}

start();
