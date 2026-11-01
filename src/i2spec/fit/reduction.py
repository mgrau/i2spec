"""Hyperfine-free line centers from the measurements: the reduction that builds the training table (rows).

Precision rows (comb, Doppler-free, emission): data.fit_datasets() (excluded and pending sets dropped, a mean giving
way to its primaries and a reprint entering once, the role decision's trained sets in and its row records out, absolute frequencies taken to zero pressure and power with each set's shift corrections,
whose allowance widens sigma); intervals within one line (pure hyperfine) are left to the hyperfine fits. Each
component value becomes the hyperfine-free line center with the model's own hyperfine offsets: center = value -
offset(line, component) [+ offset(ref line, ref component) for intervals]. All components of one line (or pair of
lines) in one set merge into one row: the inverse-variance mean, sigma = the median of the components' sigmas (they
share the set's calibration), in quadrature with the spread of the component-derived centers (n >= 2: the
hyperfine model's error) and a hyperfine-removal floor (5 kHz; 20 kHz at X v'' >= 18 or B v' >= 44).

Atlases: Salami & Ross 2005 and the APO cell, lines with sigma <= 20 MHz and depth >= 0.15, each on a scale
calibrated against the precision rows of the same lines (atlas - precision = a + b (nu - 17000 cm-1)/1000, a
3-sigma-clipped weighted fit; refitted on the core training rows alone in `stages`). Orsay part I (+23.85 ppb) and
Partie IV (+200.8 ppb, unblended lines), sigma = hypot(printed precision, 20 / 45 MHz). Martin et al. 1986 X levels
at J = 0, 20, ... (to 120 for v'' <= 85, 60 to 96, 20 above; sigma 0.005-0.02 cm-1), and Martin et al. 1983's
v'' 93, 97, 102, 103, 106 on the 1986 origin (`martin1983_rows`).

The model's own offsets (ModelOffsets): in the core the core fit's hyperfine (core_hyperfine); beyond it the B model
of hyperfine_beyond and the X model of x_hyperfine; a line beyond the core whose own intervals were measured has its
four B constants fitted to them (per_line_fit), so its offsets rest on its own data. Lines without intervals carry,
in sigma, the spread between the X eqQ models and, for B v' > 43, the per-line-fit scatter of B v' within 3.
"""
from __future__ import annotations

import csv
from collections import defaultdict

import numpy as np

from .. import hfs_params as hp
from ..constants import MHZ_PER_CM
from ..hyperfine import HyperfineParameters, level_structure, line_components
from ..observations import Line, component_rank
from . import core_hyperfine, hyperfine_beyond as hb, hyperfine_rows, x_hyperfine as xh
from .continuation import DATA_V
from . import rows as training
from .data import fit_datasets
from .rows import CORE_V, EXCLUSION, Rows, line_levels, merge

ATLAS = training.ROOT / "data" / "atlas_lines"
EMISSION_SETS = ("matyugin2012", "nesterenko2019")
DOPPLER_SIGMA = .3                 # MHz: a set whose median stated sigma reaches this is "doppler"
HF_FLOOR, HF_FLOOR_FAR = .005, .02
ATLAS_MAX_SIGMA, ATLAS_MIN_DEPTH = 20., .15
ORSAY_SCALE_PPB = {"orsay1": 23.85, "orsay4": 200.8}
ORSAY_FLOOR = {"orsay1": 20., "orsay4": 45.}
MARTIN_DISSOC_V = 90
MARTIN1983_ORIGIN = .0056          # cm-1, Martin 1986 - 1983 at their 11 common levels
MARTIN1983_LEVELS = (93, 97, 102, 103, 106)
IDENTITY = ("source", "kind", "isotopologue", "branch", "J_lower", "v_upper", "v_lower",
            "isotopologue2", "branch2", "J_lower2", "v_upper2", "v_lower2", "state", "v", "J")


# ------------------------------------------------------------------------------------------- the table

def precision_rows(offset, center, datasets=None, raw=None):
    """comb / doppler / emission rows, merged per (set, line[, ref line]), hyperfine-free.

    offset(line, component) -> MHz of the component above the line's center; center(line) -> MHz (only the `nu`
    label of interval rows uses it). Rows whose lines `offset` cannot reduce (it returns None) are left out."""
    datasets = datasets if datasets is not None else fit_datasets()
    raw = raw if raw is not None else fit_datasets(raw=True)
    groups, meta, skipped = defaultdict(list), {}, []
    for ds, ds_raw in zip(datasets, raw):
        assert ds.id == ds_raw.id and len(ds) == len(ds_raw), ds.id
        scale = MHZ_PER_CM if ds.unit == "cm-1" else 1.
        for o, o_raw in zip(ds.observations, ds_raw.observations):
            ref = o.ref_line or o.line
            if o.kind == "interval" and ref == o.line:
                continue
            off = offset(o.line, o.component) if o.component else 0.
            if o.kind == "interval" and o.ref_component and off is not None:
                r = offset(ref, o.ref_component)
                off = None if r is None else off-r
            if off is None:
                skipped.append((ds.id, str(o.line)))
                continue
            key = (ds.id, "abs" if o.kind == "frequency" else "interval", o.line, ref if o.kind == "interval" else None)
            groups[key].append((o.value*scale-off, o.uncertainty*scale, o_raw.uncertainty*scale))
        meta[ds.id] = float(np.median([o.uncertainty*scale for o in ds.observations]))
    out = []
    for (sid, kind, line, ref), comps in groups.items():
        c = np.array(comps)
        w = 1/c[:, 1]**2
        value = float(np.sum(w*c[:, 0])/np.sum(w))
        spread = float(np.std(c[:, 0], ddof=1)) if len(c) >= 2 else 0.
        far = max(line.v_lower, ref.v_lower if ref else 0) >= 18 or max(line.v_upper, ref.v_upper if ref else 0) >= 44
        sigma = float(np.sqrt(np.median(c[:, 1])**2+spread**2+(HF_FLOOR_FAR if far else HF_FLOOR)**2))
        cat = "emission" if sid in EMISSION_SETS else ("doppler" if meta[sid] >= DOPPLER_SIGMA else "comb")
        lv = line_levels(line.isotopologue, line.branch, line.J_lower, line.v_upper, line.v_lower)
        if ref is not None:
            lv += line_levels(ref.isotopologue, ref.branch, ref.J_lower, ref.v_upper, ref.v_lower, -1.)
        nu = (value if kind == "abs" else center(line))/MHZ_PER_CM
        out.append(dict(source=sid, kind=kind, category=cat, isotopologue=line.isotopologue, branch=line.branch,
                        J_lower=line.J_lower, v_upper=line.v_upper, v_lower=line.v_lower,
                        isotopologue2=ref.isotopologue if ref else "", branch2=ref.branch if ref else "",
                        J_lower2=ref.J_lower if ref else -1, v_upper2=ref.v_upper if ref else -1,
                        v_lower2=ref.v_lower if ref else -1, state="", v=-1, J=-1, value=value, sigma=sigma,
                        sigma_stated=float(c[:, 2].min()), n_components=len(c), nu=nu,
                        new=sid in training.NEW_SETS, levels=merge(lv)))
    return out, skipped


def atlas_raw(name):
    rows = []
    with open(ATLAS/f"{name}.csv", newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            sigma = float(r["uncertainty"])*MHZ_PER_CM
            depth = float(r["note"].split()[1].rstrip(";"))
            if sigma > ATLAS_MAX_SIGMA or depth < ATLAS_MIN_DEPTH:
                continue
            ln = Line.parse(r["line"])
            rows.append((ln.branch, ln.J_lower, ln.v_upper, ln.v_lower, float(r["value"])*MHZ_PER_CM, sigma))
    return rows


def calibrate(atlas, precision):
    """atlas - precision on shared lines (127I2 abs, comb + doppler): offset + slope, 3-sigma clipped, weighted."""
    ref = defaultdict(list)
    for p in precision:
        if p["kind"] == "abs" and p["isotopologue"] == "127I2" and p["category"] in training.PRECISION:
            ref[(p["branch"], p["J_lower"], p["v_upper"], p["v_lower"])].append((p["value"], p["sigma"]))
    d, s, nu = [], [], []
    for b, J, vu, vl, val, sig in atlas:
        if (b, J, vu, vl) in ref:
            vals = np.array(ref[(b, J, vu, vl)])
            w = 1/vals[:, 1]**2
            d.append(val-np.sum(w*vals[:, 0])/np.sum(w))
            s.append(np.hypot(sig, 1/np.sqrt(np.sum(w))))
            nu.append(val/MHZ_PER_CM)
    d, s, nu = map(np.array, (d, s, nu))
    w = 1/s
    x = (nu-17000)/1000
    keep = np.ones(len(d), bool)
    for _ in range(10):
        A = np.c_[np.ones_like(x), x][keep]*w[keep, None]
        coef, *_ = np.linalg.lstsq(A, d[keep]*w[keep], rcond=None)
        keep = np.abs(d-coef[0]-coef[1]*x) <= 3*s
    A = np.c_[np.ones_like(x), x][keep]*w[keep, None]
    coef, *_ = np.linalg.lstsq(A, d[keep]*w[keep], rcond=None)
    z = (d[keep]-coef[0]-coef[1]*x[keep])/s[keep]
    cov = np.linalg.inv(A.T@A)*max(float(np.mean(z**2))*keep.sum()/max(keep.sum()-2, 1), 1.)
    return dict(offset_MHz=float(coef[0]), offset_unc_MHz=float(np.sqrt(cov[0, 0])),
                slope_MHz_per_1000cm=float(coef[1]), slope_unc=float(np.sqrt(cov[1, 1])),
                n_shared=int(len(d)), n_used=int(keep.sum()), rms_z=float(np.sqrt(np.mean(z**2))),
                nu_range_cm=[float(nu.min()), float(nu.max())])


def calibrated(val, cal):
    x = (val/MHZ_PER_CM-17000)/1000
    return val-cal["offset_MHz"]-cal["slope_MHz_per_1000cm"]*x


def uncalibrated(val, cal):
    """The atlas value before `calibrated` subtracted cal."""
    a, b = cal["offset_MHz"], cal["slope_MHz_per_1000cm"]
    return (val+a-17*b)/(1-b/(1000*MHZ_PER_CM))


def orsay_rows(name):
    fname = {"orsay1": "orsay1982_part1_assigned.csv", "orsay4": "orsay1983_part4_assigned.csv"}[name]
    scale = 1+ORSAY_SCALE_PPB[name]*1e-9
    out = []
    with open(ATLAS/fname, newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            if r.get("n_assignments", "1") != "1":
                continue
            eps = float(r["eps_mk"])*1e-3*MHZ_PER_CM if r["eps_mk"] else 150.
            out.append((r["branch"], int(r["J_lower"]), int(r["v_upper"]), int(r["v_lower"]),
                        float(r["sigma_cm1"])*scale*MHZ_PER_CM, float(np.hypot(eps, ORSAY_FLOOR[name])), eps))
    return out


def martin1986_rows():
    out = []
    with open(training.ROOT/"data/x_levels/martin1986/data.csv") as handle:
        for r in csv.DictReader(handle):
            f = lambda k, s=1.: float(r[k])*s if r.get(k) else 0.      # noqa: E731
            v = int(r["v"])
            G, B, D, H, L, M = (f("G"), f("B_x10", .1), f("D_x1e8", 1e-8), -f("negH_x1e14", 1e-14),
                                -f("negL_x1e20", 1e-20), -f("negM_x1e23", 1e-23))
            Jmax = 120 if v <= 85 else 60 if v <= 96 else 20
            for J in range(0, Jmax+1, 20):
                x = J*(J+1.)
                E = G+B*x-D*x**2+H*x**3+L*x**4+M*x**5
                s = .005 if J <= 80 else .01 if J <= 100 else .02
                out.append((v, J, E*MHZ_PER_CM, s*MHZ_PER_CM))
    return out


def _table(recs):
    count = defaultdict(int)
    for r in recs:
        r["id"] = f"{r['source']}:{count[r['source']]:05d}"
        count[r["source"]] += 1
    cols = {k: np.array([r[k] for r in recs]) for k in training.COLUMNS}
    cols["id"] = cols["id"].astype(str)
    return Rows(cols, [r["levels"] for r in recs])


def build(offset, center, datasets=None, raw=None):
    """(rows, calibrations, skipped) from the data files: precision rows, the calibrated atlases, Orsay, Martin 1986."""
    recs, skipped = precision_rows(offset, center, datasets, raw)
    calib = {}
    blank = dict(isotopologue2="", branch2="", J_lower2=-1, v_upper2=-1, v_lower2=-1, state="", v=-1, J=-1, new=False,
                 n_components=1, isotopologue="127I2")
    for name, cat in (("salami_ross_2005", "atlas_sr"), ("apo_nist_2009", "atlas_apo")):
        raw_rows = atlas_raw(name)
        cal = calibrate(raw_rows, recs)
        calib[cat] = cal
        for b, J, vu, vl, val, sig in raw_rows:
            v = calibrated(val, cal)
            recs.append(dict(blank, source=name, kind="abs", category=cat, branch=b, J_lower=J, v_upper=vu, v_lower=vl,
                             value=v, sigma=sig, sigma_stated=sig, nu=v/MHZ_PER_CM,
                             levels=merge(line_levels("127I2", b, J, vu, vl))))
    for cat, name in (("orsay1", "orsay1982_part1"), ("orsay4", "orsay1983_part4")):
        for b, J, vu, vl, val, sig, eps in orsay_rows(cat):
            recs.append(dict(blank, source=name, kind="abs", category=cat, branch=b, J_lower=J, v_upper=vu, v_lower=vl,
                             value=val, sigma=sig, sigma_stated=eps, nu=val/MHZ_PER_CM,
                             levels=merge(line_levels("127I2", b, J, vu, vl))))
    for v, J, E, s in martin1986_rows():
        recs.append(dict(blank, source="martin1986", kind="level",
                         category="martin_dissoc" if v >= MARTIN_DISSOC_V else "martin", branch="", J_lower=-1,
                         v_upper=-1, v_lower=-1, state="X", v=v, J=J, value=E, sigma=s, sigma_stated=s, nu=np.nan,
                         levels=merge([(("127I2", "X", v, J), 1.), (("127I2", "X", 0, 0), -1.)])))
    return _table(recs), calib, skipped


def martin1983_rows(rows):
    """Martin 1983 Table II levels that Martin 1986 leaves out, as Martin 1986 rows (J = 0-60 by 20 up to v''=96, else
    0 and 20; sigma 0.005 cm-1, or the printed 1-sigma T if larger), on the 1986 origin."""
    template = rows.record(int(np.flatnonzero(rows.source == "martin1986")[0]))
    out = []
    with open(training.ROOT/"data/x_levels/martin1983/data.csv") as f:
        for r in csv.DictReader(f):
            v = int(r["v"])
            if v not in MARTIN1983_LEVELS:
                continue
            T, B = float(r["T"])+MARTIN1983_ORIGIN, float(r["B"])
            D, H = float(r["D_x1e7"] or 0)*1e-7, float(r["H_x1e12"] or 0)*1e-12
            s = max(.005, float(r["T_unc"] or 0)/2)*MHZ_PER_CM
            for J in range(0, (60 if v <= 96 else 20)+1, 20):
                x = J*(J+1.)
                out.append(dict(template, id=f"martin1983:{len(out):05d}", source="martin1983",
                                category="martin_dissoc", v=v, J=J, value=(T+B*x-D*x*x+H*x**3)*MHZ_PER_CM, sigma=s,
                                sigma_stated=s, new=False,
                                levels=merge([(("127I2", "X", v, J), 1.), (("127I2", "X", 0, 0), -1.)])))
    cols = {k: np.concatenate([a, np.array([o[k] for o in out], dtype=a.dtype)]) for k, a in rows.cols.items()}
    return Rows(cols, list(rows.levels)+[o["levels"] for o in out])


def excluded(rows, text=EXCLUSION):
    """Mask of rows either of whose lines is excluded."""
    line = Line.parse(text)
    mask = np.zeros(len(rows), bool)
    for suffix in ("", "2"):
        mask |= ((rows["isotopologue"+suffix] == line.isotopologue) & (rows["branch"+suffix] == line.branch)
                 & (rows["J_lower"+suffix] == line.J_lower) & (rows["v_upper"+suffix] == line.v_upper)
                 & (rows["v_lower"+suffix] == line.v_lower))
    return mask


def recalibrate(rows, calibration, precision_rows_for_fit):
    """The atlas rows of `rows` (built on `calibration`) moved to a calibration fitted on `precision_rows_for_fit`
    (the training rows' precision records); returns (rows, the new calibrations)."""
    rows, new = rows.copy(), {}
    for category in ("atlas_sr", "atlas_apo"):
        m = rows.category == category
        if not m.any():
            continue
        raw = uncalibrated(rows.value[m], calibration[category])
        new[category] = calibrate([(rows.branch[i], rows.J_lower[i], rows.v_upper[i], rows.v_lower[i], z, rows.sigma[i])
                                   for i, z in zip(np.flatnonzero(m), raw)], precision_rows_for_fit)
        rows.cols["value"][m] = calibrated(raw, new[category])
        rows.cols["nu"][m] = rows.value[m]/MHZ_PER_CM
    return rows, new


def stages(rows, calibration):
    """(core train, core benchmark, continuation train, continuation benchmark, atlas calibrations) as the fits take
    them: 127I2 only, P(13) 43-0 out, the newest sets held out; the core within X v''<=17, B v'<=43; the atlases on a
    calibration from the core's training precision rows; the continuation with Martin 1983's levels added."""
    rows = rows[np.array([all(k[0] == "127I2" for k, _ in lv) for lv in rows.levels]) & ~excluded(rows)]
    train_mask = ~rows.new.astype(bool)
    core = np.array([all(k[2] <= CORE_V[k[1]] for k, _ in lv) for lv in rows.levels])
    precision = [rows.record(i) for i in np.flatnonzero(core & train_mask & np.isin(rows.category, training.PRECISION))]
    rows, cal = recalibrate(rows, calibration, precision)
    cont = martin1983_rows(rows)
    top = np.array([all(k[2] <= DATA_V[k[1]] for k, _ in lv) for lv in cont.levels])
    return (rows[core & train_mask], rows[core & ~train_mask], cont[top & ~cont.new.astype(bool)],
            cont[cont.new.astype(bool)], cal)


def identity(rows, i):
    return tuple(str(rows.cols[k][i]) for k in IDENTITY)


def compare(a, b, tolerance=1e-9):
    """Row-by-row differences of two tables, matched by identity (not id): {only_a, only_b, changed}."""
    ia = {identity(a, i): i for i in range(len(a))}
    ib = {identity(b, i): i for i in range(len(b))}
    changed = []
    for key in sorted(set(ia) & set(ib)):
        i, j = ia[key], ib[key]
        dv, ds = float(b.value[j]-a.value[i]), float(b.sigma[j]-a.sigma[i])
        if abs(dv) > tolerance or abs(ds) > tolerance:
            changed.append(dict(identity=key, id_a=str(a.id[i]), id_b=str(b.id[j]), value_a=float(a.value[i]),
                                change_MHz=dv, sigma_change_MHz=ds))
    return dict(only_a=sorted(set(ia)-set(ib)), only_b=sorted(set(ib)-set(ia)), changed=changed)


# ------------------------------------------------------------------------------------------- the model's offsets

def far(line):
    return line.isotopologue == '127I2' and (line.v_upper > CORE_V['B'] or line.v_lower > CORE_V['X'])


class OwnHyperfine:
    def __init__(self, constants, b, theta, xmodel, xlv=None, levels=None, intervals=None, cache=None, orders=None):
        """constants, b, theta: the B fit (hyperfine_beyond; orders: its background orders where they are not the
        core fit's); xmodel: x_hyperfine.fit's model on xlv."""
        self.levels = levels or hb.Levels()
        self.intervals = intervals if intervals is not None else hb.load_data()[0]
        self.constants = {s: dict(c) for s, c in constants.items()}
        self.b, self.theta = np.asarray(b, float), dict(theta)
        self.xlv = xlv or xmodel['rad'].lv
        self.xmodel = xmodel
        self.d_x, self.delta_x = hp.D_X, hp.DELTA_X          # kHz, the core's X spin-spin constants
        self.cache_dir = cache
        self.orders = orders

    def model(self, levels_needed):
        return hb.Model(self.levels, hb.fill_terms(levels_needed, self.constants, cache=self.cache_dir),
                        orders=self.orders)

    def x_structure(self, v, J, bkt02=False):
        if v <= CORE_V['X']:
            return hb.x_structure(v, J)
        E0 = self.xlv.energy_x00(v, J)

        def params(j):
            q, c = xh.constants(self.xlv, self.xmodel, v, j)
            if bkt02:                     # eqQ_X from BKT02, the alternative; its C_X pole is not an option
                q = self.xlv.bkt02(v, j)[0]
            return HyperfineParameters(eqQ=q, C=c, d=self.d_x, delta=self.delta_x)
        return level_structure(J, params, lambda j: (self.xlv.energy_x00(v, j)-E0)*MHZ_PER_CM, symmetry='g', dJ=2)

    def offsets(self, line, constants_by_J, bkt02=False):
        Jup = hb.upper_J(line)
        E0 = self.levels.energy(line.v_upper, Jup)

        def params(j):
            c = constants_by_J[j]
            return HyperfineParameters(eqQ=c[0], C=c[1], d=c[2], delta=c[3])
        s = level_structure(Jup, params, lambda j: (self.levels.energy(line.v_upper, j)-E0)*MHZ_PER_CM,
                            symmetry='u', dJ=2)
        comps = line_components(s, self.x_structure(line.v_lower, line.J_lower, bkt02), Jup, line.J_lower)
        return {component_rank(c.label): c.offset for c in comps if c.label}


def per_line_fit(own, model, line, rows):
    """B constants of `line` (array per J of J', J'+-2) fitted to its own intervals: one shift of the four
    constants shared by the three J values, by linearized least squares (repeated until the step is below 1e-4,
    at most three times)."""
    fn = lambda v, J: model.constants(own.b, own.theta, v, J)       # noqa: E731
    shift = np.zeros(4)
    for _ in range(3):
        Js, c0, ranks, base, D = hb.linearize(own.levels, line, lambda v, J: fn(v, J)+shift)
        index = {r: i for i, r in enumerate(ranks)}
        A = D.reshape(len(ranks), len(Js), 4).sum(axis=1)
        M, y = [], []
        for r in rows:
            w = 1/np.hypot(r['sigma'], .002)
            M.append((A[index[r['rank']]]-A[index[r['ref_rank']]])*w)
            y.append((r['observed']-(base[index[r['rank']]]-base[index[r['ref_rank']]]))*w)
        M, y = np.array(M), np.array(y)
        step = np.linalg.lstsq(M, y, rcond=None)[0]
        shift = shift+step
        if np.max(np.abs(step)) < 1e-4:
            break
    resid = M@step-y
    return {j: fn(line.v_upper, j)+shift for j in Js}, shift, float(np.sqrt(np.mean(resid**2)))


class ModelOffsets:
    """offset(line, component) of the model's own hyperfine for any 127I2 line: the core fit's hyperfine in the core
    (core_offsets: {line: {rank: MHz}}, core_hyperfine.Problem.evaluate's third output), beyond it `own`
    (OwnHyperfine) with per-line fits where the line's intervals were measured."""

    def __init__(self, core_offsets, own, beyond_lines=()):
        self.core_offsets, self.own = core_offsets, own
        self.table, self.extra, self.report = beyond_offsets(own, beyond_lines) if own else ({}, {}, {})

    def model_sigma(self, rows):
        """rows with the beyond-core model's error added (in quadrature) to the sigma of the precision rows whose
        first line has model offsets only."""
        rows = rows.copy()
        for i in np.flatnonzero(np.isin(rows.category, ("comb", "doppler", "emission"))):
            ln = Line(str(rows.isotopologue[i]), str(rows.branch[i]), int(rows.J_lower[i]), int(rows.v_upper[i]),
                      int(rows.v_lower[i]))
            if ln in self.extra:
                rows.cols["sigma"][i] = float(np.hypot(rows.sigma[i], self.extra[ln]))
        return rows

    def __call__(self, line, component):
        if line.isotopologue != "127I2":
            return None
        rank = component_rank(component)
        if line in self.core_offsets:
            return self.core_offsets[line][rank]
        if line in self.table:
            return self.table[line][0][rank]
        return None


def beyond_offsets(own, lines):
    """{line: (offsets by rank, spread by rank, method)}, sigma additions {line: MHz} and a report, for 127I2 lines
    beyond the core."""
    lines = [ln for ln in lines if far(ln)]
    needed = sorted({(ln.v_upper, j) for ln in lines for j in hb.neighbours(hb.upper_J(ln))})
    model = own.model(needed)
    by_line = {}
    for r in own.intervals:
        by_line.setdefault(r["line"], []).append(r)
    table, report, model_error = {}, {}, {}
    fn = lambda v, J: model.constants(own.b, own.theta, v, J)    # noqa: E731
    for ln in sorted(set(lines), key=str):
        start = {j: fn(ln.v_upper, j) for j in hb.neighbours(hb.upper_J(ln))}
        if ln in by_line and ln.v_lower <= CORE_V["X"]:
            fitted, shift, rms = per_line_fit(own, model, ln, by_line[ln])
            o_fit, o_model = own.offsets(ln, fitted), own.offsets(ln, start)
            table[ln] = (o_fit, 0., "per-line fit to its own intervals")
            model_error[ln] = [o_model[r]-o_fit[r] for r in o_fit]
            report[str(ln)] = dict(method="per-line", n_intervals=len(by_line[ln]), fit_rms_sigma=rms,
                                   constant_shift=shift.tolist())
        else:
            o = own.offsets(ln, start)
            ob = own.offsets(ln, start, bkt02=True)
            table[ln] = (o, {r: abs(o[r]-ob[r]) for r in o}, "model offsets")
            report[str(ln)] = dict(method="model", x_eqq_sensitivity_MHz=max(table[ln][1].values()))
    extra = {}
    for ln, (o, spread, method) in table.items():
        if method != "model offsets":
            continue
        e = max(spread.values())                  # conservative: the largest component shift
        if ln.v_upper > CORE_V["B"]:
            # the per-line fits of B v' within 3; with none there, those of the nearest v' that has any
            if not model_error:
                raise ValueError(f"{ln}: no per-line fit sets the B model's error")
            gap = max(3, min(abs(k.v_upper-ln.v_upper) for k in model_error))
            near = [d for k, err in model_error.items() if abs(k.v_upper-ln.v_upper) <= gap for d in err]
            rms = float(np.sqrt(np.mean(np.square(near))))
            e = float(np.hypot(e, rms))
            report[str(ln)]["b_model_rms_MHz"] = rms
        extra[ln] = e
    return table, extra, report


def observed_lines(datasets):
    lines = set()
    for ds in datasets:
        for o in ds.observations:
            lines.update(ln for ln in (o.line, o.ref_line) if ln is not None)
    return lines


def core_offsets(lines, spec=None, theta=None):
    """{line: {rank: MHz}} of the 127I2 core lines among `lines`, from the core fit's potentials and hyperfine
    coefficients (default: the release's)."""
    from .core import release_spec
    from .curves import Config, Curve
    spec = spec or release_spec()
    config = Config.from_spec(spec["config"])
    curves = {s: Curve(s, config, spec["values"][s]) for s in ("X", "B")}
    inside = sorted({ln for ln in lines if ln.isotopologue == "127I2" and ln.v_upper <= CORE_V["B"]
                     and ln.v_lower <= CORE_V["X"]}, key=str)
    records = [r for r in hyperfine_rows.release()[0] if r["inside"]]
    energies, means = core_hyperfine.moments(curves, spec["values"], [r["line"] for r in records]+inside)
    problem = core_hyperfine.Problem(energies, means, records, extra_lines=inside)
    theta = core_hyperfine.release_coefficients(spec, problem.names) if theta is None else theta
    return problem.evaluate(np.asarray(theta, float))[2]


def release_offsets(lines, beyond=True, cache=None):
    """ModelOffsets of the released hyperfine for `lines`: core lines from the core fit, and with `beyond` the
    127I2 lines beyond the core from the released B and X models (hyperfine.json)."""
    from ..research import _entry, _read
    table = _read(_entry("balanced_expanded"), "hyperfine.json")
    core = core_offsets(lines)
    if not beyond:
        return ModelOffsets(core, None, ())
    lv = xh.XLevels()
    model = table["B"]["model"]
    own = OwnHyperfine(model["curve_constants"], model["background"], model["theta_deg"],
                       xh.model_from(table["X"]["model"], lv), xlv=lv, cache=cache)
    return ModelOffsets(core, own, [ln for ln in lines if ln.isotopologue == "127I2"])


def model_center():
    """center(line) for the `nu` labels: the released model's hyperfine-free position (nan outside its domain)."""
    from ..research import ResearchModel
    model = ResearchModel("i2spec")

    def center(line):
        try:
            return model.transition(line.v_upper, line.v_lower, line.J_lower, line.branch)*MHZ_PER_CM
        except ValueError:
            return np.nan
    return center


def release_rows(beyond=True, datasets=None, raw=None, cache=None):
    """The training stages rebuilt from the data files with the released hyperfine: (core train, core benchmark,
    continuation train, continuation benchmark, atlas calibrations, offsets). Without `beyond` the rows that need
    hyperfine beyond the core are left out (only the core stage is then complete)."""
    datasets = datasets if datasets is not None else fit_datasets()
    raw = raw if raw is not None else fit_datasets(raw=True)
    offsets = release_offsets(observed_lines(datasets), beyond=beyond, cache=cache)
    rows, calibration, _ = build(offsets, model_center(), datasets, raw)
    return (*stages(offsets.model_sigma(rows), calibration), offsets)
