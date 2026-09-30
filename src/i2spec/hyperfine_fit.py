"""Stage 2 of the refit: the hyperfine parameters, with the potentials held fixed.

Why this is a separate stage: `bodermann1998b` measures an absolute frequency and a hyperfine
splitting of the same lines at v'' = 16 and 17. The splitting is reproduced to 36 kHz; the positions
are out by 2.5 and 5.8 MHz. The error is in the term values, not in the hyperfine Hamiltonian, so a
joint fit would let the hyperfine constants absorb potential error.

What it fits to: only intra-line splittings. Within one line the rovibronic centre cancels exactly, so
these carry no information about the potentials, and Stage 1 discards them. The two stages use
disjoint data.

What it fits: additive corrections to the published BKT02 and S06 formulae (`hfs_params.Corrections`),
low-order polynomials in (v + 1/2) and J(J + 1). `free` chooses which terms vary.

The result, in docs/design/hyperfine-fit.md: no set of correction terms predicts held-out lines better
than the published formulae in every data set, so the model keeps the published formulae. What the fit
did establish:

- The Hamiltonian is right. Freeing eqQ, C, d and delta of the B state line by line brings 25 of the
  29 lines measured to <= 25 kHz within 1.5 sigma (`per_line_fit`); the published formulae manage 1.
- The error is mostly in C_B: freed alone, it takes R(145) 37-0 from 1235 kHz to 31 kHz.
- A correction fitted to five of the six v' groups of the BIPM 532 nm lines predicts the sixth better
  in three cases and worse in three (`cross_validate`), and one fitted to everything damages lines at
  v' = 6-21. Within those lines v' and J' rise together, so the data cannot fix the correction's form.
"""

from __future__ import annotations

import time
from dataclasses import dataclass, fields

import numpy as np
from scipy.optimize import least_squares

from .constants import MHZ_PER_CM, DEFAULT_PARAMETERS
from .hfs_params import Corrections
from .model import RovibronicModel
from .observations import component_rank, load_all

#: 1 sigma priors, by parameter: the 2 sigma figures of BKT02 section 7 (eqQ 50 kHz, C 2 kHz, delta
#: 3 kHz, d 4 kHz), halved. Every term of a parameter gets the same width, which the scaled variables
#: of Corrections make a comparable constraint on each.
FAMILY_PRIOR = {"eqQ": 0.025, "C": 1.0, "d": 2.0, "delta": 1.5}

NAMES = tuple(f.name for f in fields(Corrections))
PRIORS = {n: FAMILY_PRIOR[n.split("_")[0]] for n in NAMES}

#: Model error of the published formulae, MHz, added to every measurement uncertainty in quadrature.
#:
#: The splittings here are measured to 1-2 kHz, but the formulae are interpolations whose error is a
#: smooth function of (v', J') rather than noise: the best line in the whole set, R(56) 32-0, sits at
#: 24 kHz rms, and BKT02 section 7 quotes that figure itself. Weighting by the measurement
#: uncertainty alone therefore treats that structured error as if it were random, and hands the fit
#: to whichever line the formulae happen to fit worst. Without this floor one line, R(145) 37-0,
#: carried 91% of chi-squared -- 20 of 642 points -- and drove eqQ_X 53 sigma from its published
#: value while every other data set got worse. See docs/design/hyperfine-fit.md.
MODEL_FLOOR = 0.025

#: The B-state parameters freed in a single-line fit. The X-state ones are left out on purpose: within
#: one line eqQ_X is perfectly correlated with eqQ_B, and C_X with C_B, so only B is identifiable.
PER_LINE_TERMS = ("eqQ_b", "C_b", "d_b", "delta_b")


def rms(a):
    """Root mean square, which is what matters here -- not the scatter about the mean."""
    a = np.asarray(a, float)
    return float(np.sqrt(np.mean(a ** 2))) if a.size else 0.0


@dataclass
class Splitting:
    """One measured interval between two hyperfine components of a single line."""

    isotopologue: str
    branch: str
    J_lower: int
    v_upper: int
    v_lower: int
    rank: int
    ref_rank: int
    value: float          # MHz
    uncertainty: float    # MHz
    dataset: str
    key: str
    weight: float = 1.0   # 1/sqrt(sigma_meas^2 + sigma_model^2), MHz^-1


class HyperfineFit:
    #: Highest v' the published hyperfine formulae claim: BKT02 stops at 43, S06 extends to 53 and
    #: freezes above it. R(98) 58-1 of bipm2005a sits beyond both, and the formulae miss its
    #: splittings by up to 81 MHz against 5 kHz measurements -- an error no offset of the kind fitted
    #: here can absorb, and one that would otherwise dominate chi-squared completely.
    V_UPPER_VALID = 53

    def __init__(self, datasets=None, parameters=DEFAULT_PARAMETERS, max_uncertainty=None,
                 max_v_upper=V_UPPER_VALID, model_floor=MODEL_FLOOR, loss="soft_l1", f_scale=3.0,
                 free=None):
        self.parameters = parameters
        self.free = tuple(free) if free is not None else NAMES
        unknown = set(self.free) - set(NAMES)
        if unknown:
            raise ValueError(f"not a correction term: {sorted(unknown)}")
        self.model_floor = float(model_floor)
        self.loss, self.f_scale = loss, float(f_scale)
        self._models: dict[str, RovibronicModel] = {}
        self.data: list[Splitting] = []
        self.skipped: list[tuple[str, str]] = []
        for ds in (datasets if datasets is not None else load_all()):
            scale = MHZ_PER_CM if ds.meta["unit"] == "cm-1" else 1.0
            for o in ds.observations:
                ref = o.ref_line or o.line
                if o.kind != "interval" or ref != o.line or not o.component or not o.ref_component:
                    continue
                if max_uncertainty and o.uncertainty * scale > max_uncertainty:
                    continue
                if max_v_upper is not None and o.line.v_upper > max_v_upper:
                    self.skipped.append((str(o.line), f"v' = {o.line.v_upper} is beyond the formulae"))
                    continue
                try:
                    rank, ref_rank = component_rank(o.component), component_rank(o.ref_component)
                except ValueError as exc:
                    self.skipped.append((str(o.line), str(exc)))
                    continue
                L = o.line
                sigma = o.uncertainty * scale
                self.data.append(Splitting(L.isotopologue, L.branch, L.J_lower, L.v_upper, L.v_lower,
                                           rank, ref_rank, o.value * scale, sigma,
                                           ds.id, f"{L} {o.component}-{o.ref_component}",
                                           1.0 / np.hypot(sigma, self.model_floor)))
        self.lines = sorted({(d.isotopologue, d.branch, d.J_lower, d.v_upper, d.v_lower)
                             for d in self.data})
        self.x0 = np.zeros(len(self.free))
        self.prior = np.array([PRIORS[n] for n in self.free])
        self.weight = np.array([d.weight for d in self.data])
        self.sigma = np.array([d.uncertainty for d in self.data])

    def _model(self, isotopologue) -> RovibronicModel:
        if isotopologue not in self._models:
            self._models[isotopologue] = RovibronicModel(isotopologue, self.parameters)
        return self._models[isotopologue]

    def corrections(self, x) -> Corrections:
        return Corrections(**dict(zip(self.free, map(float, x))))

    def offsets(self, x):
        """{line: array of main-component offsets in MHz} at the given corrections."""
        c = self.corrections(x)
        out = {}
        for key in self.lines:
            iso, br, J, vu, vl = key
            _, comps = self._model(iso).hyperfine_components(vu, vl, J, br, corrections=c, table=None)
            main = sorted([q for q in comps if q.label], key=lambda q: q.offset)
            out[key] = np.array([q.offset for q in main])
        return out

    def deviations(self, x):
        """model - observed, in MHz, for every splitting."""
        off = self.offsets(x)
        out = np.empty(len(self.data))
        for i, d in enumerate(self.data):
            a = off[(d.isotopologue, d.branch, d.J_lower, d.v_upper, d.v_lower)]
            if max(d.rank, d.ref_rank) > len(a):
                out[i] = 0.0                       # component beyond those the model resolves
                continue
            out[i] = (a[d.rank - 1] - a[d.ref_rank - 1]) - d.value
        return out

    def residuals(self, x):
        """(model - observed), weighted by measurement and model error together."""
        return self.deviations(x) * self.weight

    def scaled_residuals(self, u):
        """Residuals plus the prior block, as a function of parameters in units of their priors."""
        return np.concatenate([self.residuals(u * self.prior), np.asarray(u, float)])

    def report(self, x, label=""):
        khz = self.deviations(x) * 1e3
        chi2 = np.mean((khz * 1e-3 / self.sigma) ** 2)
        print(f"{label:<10} n={len(khz)}  rms {rms(khz):8.1f} kHz  max |{np.abs(khz).max():9.1f}|  "
              f"chi2/n {chi2:9.2f}")
        return khz

    def rms_by_dataset(self, x):
        khz = self.deviations(x) * 1e3
        return {name: float(rms(khz[[d.dataset == name for d in self.data]]))
                for name in sorted({d.dataset for d in self.data})}

    def rms_by_line(self, x):
        khz = self.deviations(x) * 1e3
        out = {}
        for key in self.lines:
            sel = [(d.isotopologue, d.branch, d.J_lower, d.v_upper, d.v_lower) == key for d in self.data]
            out[key] = float(rms(khz[sel]))
        return out

    def line_of(self, d: Splitting):
        return (d.isotopologue, d.branch, d.J_lower, d.v_upper, d.v_lower)

    def groups(self, kind="line"):
        """A label per splitting for cross-validation: its line, or its (data set, v') group."""
        if kind == "line":
            keys = [self.line_of(d) for d in self.data]
        elif kind == "v_upper":
            keys = [(d.dataset, d.v_upper) for d in self.data]
        else:
            raise ValueError(kind)
        index = {k: i for i, k in enumerate(dict.fromkeys(keys))}
        return np.array([index[k] for k in keys])

    def per_line_fit(self, key, terms=PER_LINE_TERMS, max_nfev=60):
        """Fit the given terms to one line alone, weighted by measurement uncertainty only.

        Returns the rms before and after, in kHz and normalised (rms of deviation/sigma, which is 1 at
        the measurement noise), the fitted corrections and their 1 sigma errors. The errors are scaled
        up by sqrt(chi2/dof) when that exceeds one, never down. Judge a line by the normalised figure:
        the tables often mix 5 kHz and 1 MHz components within one line, and the kHz rms then reports
        only the loose ones.
        """
        rows = [i for i, d in enumerate(self.data) if self.line_of(d) == key]
        D = [self.data[i] for i in rows]
        iso, br, J, vu, vl = key
        model = self._model(iso)
        step = np.array([PRIORS[t] for t in terms])
        sig = np.array([d.uncertainty for d in D])

        def dev(u):
            c = Corrections(**dict(zip(terms, map(float, u * step))))
            _, comps = model.hyperfine_components(vu, vl, J, br, corrections=c, table=None)
            a = np.array(sorted(q.offset for q in comps if q.label))
            return np.array([(a[d.rank - 1] - a[d.ref_rank - 1]) - d.value
                             if max(d.rank, d.ref_rank) <= len(a) else 0.0 for d in D])

        res = least_squares(lambda u: dev(u) / sig, np.zeros(len(terms)), diff_step=1e-3, max_nfev=max_nfev)
        scale = max(np.sum(res.fun ** 2) / max(len(D) - len(terms), 1), 1.0)
        cov = np.linalg.pinv(res.jac.T @ res.jac) * scale
        before, after = dev(np.zeros(len(terms))), dev(res.x)
        return {"n": len(D), "rms_before": rms(before) * 1e3, "rms_after": rms(after) * 1e3,
                "normalised_before": rms(before / sig), "normalised_after": rms(after / sig),
                "x": dict(zip(terms, res.x * step)), "error": dict(zip(terms, np.sqrt(np.diag(cov)) * step)),
                "median_uncertainty": float(np.median(sig)) * 1e3}

    def jacobian(self, fraction=0.01):
        """d(deviation)/d(term), MHz per unit of each free term, by central differences.

        The step is `fraction` of each term's prior. The splittings are linear in the corrections to
        better than 0.001 kHz over a full prior step, but a step of several priors can reorder
        near-degenerate components, which the rank matching then reads as a jump.
        """
        n = len(self.free)
        out = np.empty((len(self.data), n))
        for k in range(n):
            e = np.zeros(n)
            e[k] = fraction * self.prior[k]
            out[:, k] = (self.deviations(e) - self.deviations(-e)) / (2 * e[k])
        return out

    def cross_validate(self, terms, groups, jacobian=None, deviations=None):
        """Held-out predictions (MHz) of every splitting, linearised about the published formulae.

        For each group, the terms are fitted to every other group -- ridge-regularised by their priors,
        weighted by `weight` -- and the result predicts the held-out group. Returns deviations the same
        shape as `deviations(x0)`, so rms() of any subset is a held-out error.
        """
        J = self.jacobian() if jacobian is None else jacobian
        dev0 = self.deviations(self.x0) if deviations is None else deviations
        cols = [self.free.index(t) for t in terms]
        pred = dev0.copy()
        if not cols:
            return pred
        for g in np.unique(groups):
            test = groups == g
            pred[test] = dev0[test] + J[np.ix_(test, cols)] @ self._ridge(J, dev0, cols, ~test)
        return pred

    def linear_fit(self, terms, jacobian=None, deviations=None):
        """The terms fitted to all splittings on the linearised problem: (coefficients, deviations)."""
        J = self.jacobian() if jacobian is None else jacobian
        dev0 = self.deviations(self.x0) if deviations is None else deviations
        cols = [self.free.index(t) for t in terms]
        coef = self._ridge(J, dev0, cols, np.ones(len(dev0), bool))
        return dict(zip(terms, coef)), dev0 + J[:, cols] @ coef

    def _ridge(self, J, dev0, cols, rows):
        prior = self.prior[cols]
        A = np.vstack([J[np.ix_(rows, cols)] * prior * self.weight[rows, None], np.eye(len(cols))])
        b = np.concatenate([-dev0[rows] * self.weight[rows], np.zeros(len(cols))])
        u, *_ = np.linalg.lstsq(A, b, rcond=None)
        return u * prior

    def run(self, max_nfev=40, verbose=2):
        t0 = time.perf_counter()
        before = self.rms_by_dataset(self.x0)
        res = least_squares(self.scaled_residuals, np.zeros(len(self.free)), method="trf",
                            diff_step=1e-2, max_nfev=max_nfev, verbose=verbose,
                            loss=self.loss, f_scale=self.f_scale)
        x = res.x * self.prior
        return {"x": x, "corrections": self.corrections(x),
                "rms_before": before, "rms_after": self.rms_by_dataset(x),
                "sigma": dict(zip(self.free, res.x)), "message": str(res.message),
                "seconds": time.perf_counter() - t0}
