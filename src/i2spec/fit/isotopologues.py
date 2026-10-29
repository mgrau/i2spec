"""129I2 and 127I129I on the 127I2 potentials: the isotope shift and the conditioning of the covariance on it.

Levels: the same potentials solved with the isotopologue's reduced mass. Born-Oppenheimer breakdown in Le Roy's
convention, with rho = mu(127I2)/mu: the rotational function q(R), fitted to 127I2, scales by rho, and an adiabatic
correction (1 - rho) dV_B(R) on the B state, relative to X, dV_B = (1 - y)^p (a0 + a1 y), y = (R - Re_B)/(R + Re_B),
is fitted to the isotope shifts the CIPM compilation ties to 127I2 lines (bipm2003a/d), p chosen among POWERS by the
leave-one-line-out rms (a training objective).
Hyperfine, per nucleus (eqQ x Q ratio, C x g ratio, d and delta x g_a g_b; Salumbides et al. 2006 Table 1 ratios):
  * B core (127I2-equivalent v' <= 43): the core fit's radial functions averaged over the isotopologue's own
    wavefunctions;
  * X core: the fixed formulae at the 127I2-equivalent u = v + 1/2 (same J = 0 energy) and x = rho J(J+1);
  * beyond the core: the 127I2 constants interpolated at the same J = 0 energy and rho J(J+1) (the B v' 77/78
    1g resonance is 127I2's own and is not carried over).
127I129I has distinguishable nuclei: every J, both parities of I, nucleus-specific eqQ and C.
"""
from __future__ import annotations

from math import sqrt

import numpy as np
from numpy.polynomial.polynomial import polyvander
from scipy.interpolate import CubicSpline

from ..constants import ISOTOPOLOGUES, MHZ_PER_CM, NUCLEAR_SPIN, reduced_mass
from ..hyperfine import HyperfineParameters, level_structure, line_components
from ..observations import component_rank, load_all
from .. import research_expanded as rx
from ..research import ResearchModel

POWERS = (0, 3, 5, 8)
FLOOR_MHZ = .1                     # added in quadrature to each line's sigma in the shift fit (hyperfine, transfer)
J_SPLIT, J_MAX = 200, 445

GAMMA_Q = {127: 1.0, 129: 0.701213}          # eqQ(129I)/eqQ(127I) (S06 Table 1)
GAMMA_MU = {127: 1.0, 129: 0.6655}           # g(129I)/g(127I)
ISOS = ("129I2", "127I129I")
ALL = ("127I2",)+ISOS


def ratio(iso):
    return reduced_mass("127I2")/reduced_mass(iso)


class IsoModel:
    """Line centers and hyperfine components of one isotopologue."""

    def __init__(self, base, iso, shift=None):
        self.base, self.iso = base, iso
        # the runtime's isotopologue levels (research_expanded.ExpandedLevels; shift: an IsotopeShift or None)
        self.levels = rx.ExpandedLevels(base._base.coefficients, base._spec["structure"], base._spec["values"], iso, shift)
        self.rho = self.levels.rho
        a, b = ISOTOPOLOGUES[iso]
        self.i1, self.i2 = NUCLEAR_SPIN[a], NUCLEAR_SPIN[b]
        self.symmetry = (lambda s: "g" if s == "X" else "u") if a == b else (lambda s: None)
        self.scale = dict(eqQ=GAMMA_Q[a], C=GAMMA_MU[a], d=GAMMA_MU[a]*GAMMA_MU[b], delta=GAMMA_MU[a]*GAMMA_MU[b])
        self.eqQ_ratio, self.C_ratio = GAMMA_Q[b]/GAMMA_Q[a], GAMMA_MU[b]/GAMMA_MU[a]
        self.core = base._spec["core"]
        self._vref, self._moments, self._ref = {}, {}, {}

    def energy(self, state, v, J):
        return self.levels.energy(state, v, J)

    # --- the 127I2-equivalent level -----------------------------------------------------------------
    def v_ref(self, state, v):
        """The fractional 127I2 v with this level's J = 0 energy."""
        if state not in self._vref:
            ref = self.base._levels_x
            n = self.levels.grids[state][-1][0]
            E = []
            for w in range(n+1):
                try:
                    E.append(ref.energy(state, w, 0))
                except ValueError:
                    break
            self._vref[state] = CubicSpline(np.array(E), np.arange(len(E)))
        return float(self._vref[state](self.energy(state, v, 0)))

    def j_ref(self, J):
        return sqrt(self.rho*J*(J+1)+.25)-.5

    def _ref_constants(self, state, w, j):
        """127I2 constants (eqQ, C, d, delta) of integer level (w, j): core couplings or the beyond-core tables."""
        key = state, w, j
        if key not in self._ref:
            p = self.base._couplings(state, w, j)
            self._ref[key] = np.array([p.eqQ, p.C, p.d, p.delta])
        return self._ref[key]

    def _interpolated(self, state, v, J):
        vr, jr = self.v_ref(state, v), self.j_ref(J)
        w0, j0 = int(np.floor(vr)), int(np.floor(jr))
        fv, fj = vr-w0, jr-j0
        out = 0.
        for dw, cw in ((0, 1-fv), (1, fv)):
            for dj, cj in ((0, 1-fj), (1, fj)):
                if cw*cj:
                    out = out+cw*cj*self._ref_constants(state, w0+dw, j0+dj)
        return out

    def _b_core(self, v, J):
        """B core: the radial functions averaged over this isotopologue's wavefunction."""
        key = v, J
        if key not in self._moments:
            R, density, _ = self.levels.densities("B", v, J)
            spec = self.base._base._spec
            re = np.exp(spec["values"]["B"]["logRe"])
            z = (R**2-re**2)/(R**2+re**2)/spec["hyperfine"]["radial_scale"]["B"]
            mean = density[:, v]@polyvander(z, 7)
            coefficients = spec["hyperfine"]["coefficients"]["B"]
            self._moments[key] = np.array([float(np.dot(coefficients[p], mean[:len(coefficients[p])]))
                                           for p in ("eqQ", "C", "d", "delta")])
        return self._moments[key]

    def _x_core(self, v, J):
        hf = self.base._base._spec["hyperfine"]
        u, x = self.v_ref("X", v)+.5, self.rho*J*(J+1)
        q, c = hf["fixed_x_constants"]["eqQ"], hf["fixed_x_constants"]["C"]
        e = self.energy("X", v, 0)-self.base._levels_x.energy("X", 0, 0)
        return np.array([q[0]+q[1]*u+q[2]*u*u+q[3]*x+q[4]*u*x, c[0]+c[1]*u+c[2]/(e-c[3]),
                         hf["fixed_x_constants"]["d"], hf["fixed_x_constants"]["delta"]])

    def where(self, state, v):
        """'core' or 'beyond': by the 127I2-equivalent v."""
        return "core" if self.v_ref(state, v) <= self.core[state]+1e-9 else "beyond"

    def couplings(self, state, v, J):
        if self.where(state, v) == "core":
            c = self._b_core(v, J) if state == "B" else self._x_core(v, J)
        else:
            c = self._interpolated(state, v, J)
        s = self.scale
        return HyperfineParameters(eqQ=c[0]*s["eqQ"], C=c[1]*s["C"], d=c[2]*s["d"], delta=c[3]*s["delta"])

    def structure(self, state, v, J):
        E = self.energy(state, v, J)
        return level_structure(J, lambda j: self.couplings(state, v, j),
                               lambda j: (self.energy(state, v, j)-E)*MHZ_PER_CM, i1=self.i1, i2=self.i2,
                               symmetry=self.symmetry(state), eqQ_ratio=self.eqQ_ratio, C_ratio=self.C_ratio, dJ=2)

    def line(self, branch, J, vu, vl):
        Ju = J+(1 if branch == "R" else -1)
        center = (self.energy("B", vu, Ju)-self.energy("X", vl, J))*MHZ_PER_CM
        return center, line_components(self.structure("B", vu, Ju), self.structure("X", vl, J), Ju, J)

    def mean_y(self, v, J, re):
        """<(R - Re)/(R + Re)> of a B level (its first-order sensitivity to the shift's slope)."""
        R, density, _ = self.levels.densities("B", v, J)
        return float(density[:, v]@((R-re)/(R+re)))


# --- the CIPM isotopologue rows -----------------------------------------------------------------------

def _label(line):
    return f"{line.isotopologue} {line.branch}({line.J_lower}) {line.v_upper}-{line.v_lower}"


def cipm_rows():
    """Isotopologue component frequencies, each tied through measured intervals to a measured 127I2 frequency."""
    obs = [o for ds in load_all() if ds.id in ("bipm2003a", "bipm2003d") for o in ds.observations]
    absolute = {(_label(o.line), o.component): (o.value, o.uncertainty) for o in obs if o.kind == "frequency"}
    links = {(_label(o.line), o.component): o for o in obs if o.ref_line is not None and o.ref_line != o.line}

    def resolve(key):
        if key in absolute:
            return absolute[key]
        o = links[key]
        f, u = resolve((_label(o.ref_line), o.ref_component))
        return f+o.value, float(np.hypot(u, o.uncertainty))

    rows = []
    for key, o in links.items():
        if o.line.isotopologue == "127I2":
            continue
        f, u = resolve(key)
        rows.append(dict(iso=o.line.isotopologue, line=key[0], branch=o.line.branch, J=o.line.J_lower,
                         vu=o.line.v_upper, vl=o.line.v_lower, component=o.component, rank=component_rank(o.component),
                         value=f, sigma=u, own=o.uncertainty, note=o.note or ""))
    return rows


def evaluate(models, rows, s06=None):
    """Each row's model component frequency (and the 2008-formula one, if s06 holds RovibronicModels)."""
    cache = {}
    for r in rows:
        key = r["iso"], r["line"]
        if key not in cache:
            center, comps = models[r["iso"]].line(r["branch"], r["J"], r["vu"], r["vl"])
            offsets = {c.label: c.offset for c in comps if c.label}
            other = None
            if s06 is not None:
                nu0, oc = s06[r["iso"]].hyperfine_components(r["vu"], r["vl"], r["J"], r["branch"])
                other = {c.label: c.offset for c in oc if c.label}
            cache[key] = center, offsets, other
        center, offsets, other = cache[key]
        name = f"a{r['rank']}"
        r["center"], r["offset"] = center, offsets[name]
        r["offset_s06"] = other[name] if other else None
        r["n_main"] = len(offsets)
    return rows


def intraline(rows, which="offset"):
    """Component residuals about each line's mean (MHz): the hyperfine pattern's error alone."""
    out = []
    by = {}
    for r in rows:
        by.setdefault(r["line"], []).append(r)
    for line, rs in by.items():
        if len(rs) < 2:
            continue
        d = np.array([r["value"]-r[which] for r in rs])
        w = 1/np.array([max(r["own"], .03) for r in rs])**2
        mean = float(np.sum(w*d)/np.sum(w))
        for r, x in zip(rs, d):
            out.append(dict(line=line, component=r["component"], resid=float(x-mean), own=r["own"]))
    return out


def fit_shift(rows=None):
    """The isotope shift: every form p of POWERS fitted, the one with the lowest leave-one-line-out rms kept; and
    the hyperfine floor of the isotopologue components. Compared with the 2008 model's formulae on the same rows."""
    from ..model import RovibronicModel
    base = ResearchModel("balanced_expanded")
    rows = rows if rows is not None else cipm_rows()
    models = {iso: IsoModel(base, iso) for iso in ISOS}
    evaluate(models, rows, {iso: RovibronicModel(iso) for iso in ISOS})
    ok = [r for r in rows if r["own"] <= 1.]
    re_b = float(np.exp(base._base._spec["values"]["B"]["logRe"]))
    lines = {}
    for r in ok:
        lines.setdefault((r["iso"], r["line"]), []).append(r)
    data = []
    for (iso, line), rs in sorted(lines.items()):
        c = np.array([r["value"]-r["offset"] for r in rs])
        w = 1/np.array([max(r["sigma"], .03) for r in rs])**2
        r0 = rs[0]
        Ju = r0["J"]+(1 if r0["branch"] == "R" else -1)
        R, density, _ = models[iso].levels.densities("B", r0["vu"], Ju)
        y = (R-re_b)/(R+re_b)
        means = {p: [float(density[:, r0["vu"]]@(1-y)**p), float(density[:, r0["vu"]]@((1-y)**p*y))] for p in POWERS}
        data.append(dict(isotopologue=iso, line=line, n_components=len(rs), measured_center_MHz=float(np.sum(w*c)/np.sum(w)),
                         sigma_MHz=float(min(r["sigma"] for r in rs)), unshifted_model_MHz=r0["center"],
                         one_minus_rho=1-models[iso].rho, v_upper=r0["vu"], J_upper=Ju,
                         v_ref=models[iso].v_ref("B", r0["vu"]), means=means))
    d = np.array([x["unshifted_model_MHz"]-x["measured_center_MHz"] for x in data])
    s = np.sqrt(np.array([x["sigma_MHz"] for x in data])**2+FLOOR_MHZ**2)
    forms = {}
    for p in POWERS:
        A = np.array([x["one_minus_rho"]*MHZ_PER_CM*np.array(x["means"][p]) for x in data])
        c = np.linalg.lstsq(A/s[:, None], -d/s, rcond=None)[0]
        res = d+A@c
        chi2 = float(np.sum((res/s)**2))
        cov = np.linalg.inv((A/s[:, None]).T@(A/s[:, None]))*max(1., chi2/(len(d)-2))
        loo = []
        for k in range(len(d)):
            keep = np.arange(len(d)) != k
            ck = np.linalg.lstsq(A[keep]/s[keep, None], -d[keep]/s[keep], rcond=None)[0]
            loo.append(float(d[k]+A[k]@ck))
        forms[p] = dict(power=p, coefficients=c.tolist(), covariance=cov.tolist(), residuals_MHz=res.tolist(),
                        chi2=chi2, leave_one_line_out_MHz=loo, loo_rms_MHz=float(np.sqrt(np.mean(np.square(loo)))))
    best = min(forms.values(), key=lambda f: f["loo_rms_MHz"])
    others = sorted((f for f in forms.values() if f["power"] != best["power"] and f["power"] > 0),
                    key=lambda f: abs(f["power"]-best["power"]))[:2]
    # the hyperfine floor of a component: the excess of the intraline residuals over the stated uncertainties
    res = intraline(ok, "offset")
    r2, own2 = np.array([x["resid"]**2 for x in res]), np.array([x["own"]**2 for x in res])
    lo, hi = 0., 1.
    for _ in range(60):
        mid = (lo+hi)/2
        lo, hi = (mid, hi) if np.mean(r2/(own2+mid**2)) > 1 else (lo, mid)
    s06 = intraline(ok, "offset_s06")
    hf = dict(floor_MHz=round(hi, 4), n_intervals=len(res),
              rms_kHz=1e3*float(np.sqrt(np.mean(r2))), rms_kHz_2006_formulae_fitted_to_them=1e3*float(np.sqrt(np.mean([x["resid"]**2 for x in s06]))),
              median_stated_kHz=1e3*float(np.median([x["own"] for x in res])),
              status="127I2 functions scaled per nucleus, tested on the CIPM 129I2 and 127I129I patterns (v' 6-14)")
    out = dict(data=data, forms={str(k): v for k, v in forms.items()}, chosen=best["power"],
               alternatives=[o["power"] for o in others], re_B=re_b, hyperfine=hf,
               rho={iso: reduced_mass("127I2")/reduced_mass(iso) for iso in ALL})
    return out


def core_rows(data):
    """fit_shift's line centers (its "data") as core-fit rows (category "isotopologue"): with
    core.Settings.isotope_shift the core fit sees the other isotopologues, which fix the split between the potential's
    shape and q(R) that the 127I2 lines leave open."""
    from . import rows as training
    cols = {k: [] for k in training.COLUMNS}
    levels = []
    for x in data:
        iso, label = x["isotopologue"], x["line"].split(" ", 1)[1]
        branch, J = label[0], int(label[2:label.index(")")])
        vu, vl = map(int, label.split()[1].split("-"))
        row = dict(id=f"cipm:{x['line']}", source="cipm", kind="abs", category="isotopologue", isotopologue=iso,
                   branch=branch, J_lower=J, v_upper=vu, v_lower=vl, isotopologue2="", branch2="", J_lower2=-1,
                   v_upper2=-1, v_lower2=-1, state="", v=-1, J=-1, value=x["measured_center_MHz"],
                   sigma=x["sigma_MHz"], sigma_stated=x["sigma_MHz"], n_components=x["n_components"],
                   nu=x["measured_center_MHz"]/MHZ_PER_CM, new=False)
        for k in training.COLUMNS:
            cols[k].append(row[k])
        levels.append(tuple(training.line_levels(iso, branch, J, vu, vl)))
    return training.Rows({k: np.array(v) for k, v in cols.items()}, levels)


def domains(extension):
    """Each isotopologue's v range: the levels whose 127I2 equivalent lies in 127I2's (X 108, B 80)."""
    base = ResearchModel("balanced_expanded")
    out = {"127I2": dict(X_v_max=108, B_v_max=80, J_max=J_MAX)}
    for iso in ("129I2", "127I129I"):
        levels = rx.ExpandedLevels(base._base.coefficients, base._spec["structure"], base._spec["values"], iso,
                                   rx.IsotopeShift(extension["isotope_shift"]), J_SPLIT)
        m = rx.IsotopologueMap(levels, base._levels_x)
        top = {}
        for state, vmax in (("X", 108), ("B", 80)):
            v = vmax
            while True:
                try:
                    if m.v_ref(state, v+1) > vmax+.5:
                        break
                except ValueError:
                    break
                v += 1
            top[state] = v
        out[iso] = dict(X_v_max=top["X"], B_v_max=top["B"], J_max=J_MAX)
    return out


def extension_record(fitted, model_id="research-20261021-i2spec", base_id="research-20261020-balanced-expanded"):
    """isotopologues.json from fit_shift's result (before `condition` adds the conditioned covariance)."""
    best = fitted["forms"][str(fitted["chosen"])]
    rho129 = fitted["rho"]["129I2"]
    shift = dict(
        form="dV_B(R) = (1 - y)^p (a0 + a1 y), y = (R - Re_B)/(R + Re_B), cm-1 per unit (1 - rho), rho = mu(127I2)/mu;"
             " added to the B potential (X is the reference); vanishes towards the separated atoms",
        power=best["power"], re_B=fitted["re_B"], coefficients=best["coefficients"], covariance=best["covariance"],
        alternatives=[dict(power=p, coefficients=fitted["forms"][str(p)]["coefficients"]) for p in fitted["alternatives"]],
        floor_MHz=best["loo_rms_MHz"], floor_one_minus_rho=1-rho129,
        data_v_ref=[min(x["v_ref"] for x in fitted["data"]), max(x["v_ref"] for x in fitted["data"])],
        selection={str(p): dict(loo_rms_MHz=f["loo_rms_MHz"], chi2=f["chi2"]) for p, f in fitted["forms"].items()},
        selected_by="leave-one-line-out rms over the seven training lines (a training objective)",
        data=[{k: v for k, v in x.items() if k != "means"} for x in fitted["data"]],
        source="CIPM isotope-shift intervals of bipm2003a (633 nm) and bipm2003d (532 nm), each component tied to a"
               " measured 127I2 frequency; the pending sets of the same lines stay blind")
    extension = dict(
        schema_version=1, model=model_id, extends=base_id, isotopologues=list(ALL),
        rotation="q(R) of each state scaled by rho (Le Roy; as the 2008 model's alpha)",
        isotope_shift=shift,
        nuclear_ratios=dict(eqQ={"127": 1.0, "129": 0.701213}, g={"127": 1.0, "129": 0.6655},
                            source="Salumbides et al., Mol. Phys. 104, 2641 (2006), Table 1"),
        hyperfine=dict(fitted["hyperfine"], rule=(
            "per nucleus: eqQ x Q ratio, C x g ratio, d and delta x g_1 g_2; B in the core: the core fit's radial "
            "functions over the isotopologue's wavefunction; X in the core: the fixed formulae at the 127I2-equivalent "
            "u, x = rho J(J+1) and energy; beyond the core: the 127I2 constants interpolated at the equivalent (v, J); "
            "127I2's B v' 77/78 1g resonance is not carried over")),
        j_split=J_SPLIT,
        j_rule=("levels above J = j_split are solved on the long grid with the barrier rule; constants beyond the "
                "hyperfine tables' last J are held at it; a level's model sigma grows as J(J+1)/(J_r(J_r+1)) beyond "
                "J_r, the J the training rows reach at its v or within 2 in v, plus 15 (a rule, not fitted)"),
        equivalents="an isotopologue level's 127I2 equivalent: the fractional v of the same J = 0 energy, and J with "
                    "J(J+1) scaled by rho; the frozen uncertainty tables are read at the nearest integer level",
        pending_outcomes_used=False)
    extension["domains"] = domains(extension)
    return extension


def condition(ext, models):
    """Condition the fit's covariance on the seven CIPM isotopologue lines (research_expanded.IsotopeShiftUncertainty):
    the joint precision of the whitened fit parameters (prior identity) and the two shift coefficients (flat).
    ext: the extension record (updated in place and returned); models: {iso: the research model with that extension}."""
    from .. import research
    ms = models
    shift = rx.IsotopeShift(ext["isotope_shift"])
    A, Bm, w = [], [], []
    for x in ext["isotope_shift"]["data"]:
        m = ms[x["isotopologue"]]
        unc = m._uncertainty_model()
        branch, J, vu, vl = research.parse_label(x["line"].split(" ", 1)[1])
        Ju = J+(1 if branch == "R" else -1)
        A.append(unc.level_factor(m._levels_x, "B", vu, Ju)-unc.level_factor(m._levels_x, "X", vl, J))
        R, density, _ = m._levels_x.densities("B", vu, Ju)
        Bm.append((1-m._levels_x.rho)*MHZ_PER_CM*(density[:, vu]@shift.basis(R)))
        w.append(1/(x["sigma_MHz"]**2+FLOOR_MHZ**2))
    A, Bm, W = np.array(A), np.array(Bm), np.diag(w)
    K = A.shape[1]
    P = np.block([[np.eye(K)+A.T@W@A, A.T@W@Bm], [Bm.T@W@A, Bm.T@W@Bm]])
    S = np.linalg.inv(P)
    S = (S+S.T)/2
    L = np.linalg.cholesky(S)
    s = ext["isotope_shift"]
    s["joint_factor"] = L.tolist()
    s["joint_note"] = ("Cholesky factor of the covariance of (the fit's %d whitened parameters, a0, a1) given the "
                       "seven CIPM isotopologue lines, each with sigma (+) %.1f MHz; the shift's marginal covariance "
                       "is its last 2 x 2 block" % (K, FLOOR_MHZ))
    s["covariance_marginal"] = S[K:, K:].tolist()
    line = np.hstack([A, Bm])
    s["conditioned_line_sigma_MHz"] = {x["line"]: float(np.sqrt(v@S@v)) for x, v in zip(s["data"], line)}
    return ext


REFERENCE_LINES = {
    "127I2": ["R(56) 32-0", "P(53) 32-0", "R(127) 11-5", "R(47) 9-2", "R(98) 58-1", "R(20) 70-0", "R(10) 2-90",
              "P(210) 20-3", "R(260) 10-2", "P(330) 15-4", "R(401) 22-6"],
    "129I2": ["P(54) 8-4", "P(33) 6-3", "P(69) 12-6", "R(60) 8-4", "P(110) 10-2", "R(113) 14-4", "R(56) 32-0",
              "P(120) 45-1", "R(20) 70-0", "R(300) 12-3"],
    "127I129I": ["P(33) 6-3", "R(56) 32-0", "R(41) 51-0", "P(87) 20-1", "R(250) 15-2"],
}


def references(model_id="research-20261021-i2spec"):
    """Reference predictions of every isotopologue from the installed model `model_id`."""
    from .. import research
    rows = []
    for iso, labels in REFERENCE_LINES.items():
        m = research.ResearchModel(model_id, iso)
        for label in labels:
            try:
                center, comps = m._line(research.parse_label(label))
            except ValueError:          # outside this isotopologue's domain
                continue
            rows.append(dict(isotopologue=iso, line=label, component=None, ref_component=None, predicted_MHz=center))
            names = [c.label for c in comps if c.label]
            for name in names:
                rows.append(dict(isotopologue=iso, line=label, component=name, ref_component=None,
                                 predicted_MHz=m.frequency(label, name)))
            rows.append(dict(isotopologue=iso, line=label, component=names[-1], ref_component=names[0],
                             predicted_MHz=m.frequency(label, names[-1])-m.frequency(label, names[0])))
    return rows
