"""The core fit's potentials: one Morse/long-range (MLR) curve and one rotational function per state.

V(R) = Te + De (1 - u(R)/u(Re) exp(-beta(R) y_p(R)))^2, u(R) = sum_n C_n/R^n (inherited dispersion series),
beta(R) = y_p ln(2 De/u(Re)) + (1 - y_p) sum_i c_i T_i(z), z = y_q in the reference radius ref Re (X: mapped to
(y_q - 0.2)/0.6), y_k(R) = (R^k - ref^k)/(R^k + ref^k). The rotational function f(R) = 1 + (2Re/(R+Re)) sum_i
alpha_i T_i(w), w = (R-Re)/(R+Re) mapped to (w - 0.3)/0.4, multiplies the centrifugal term hbar^2 J(J+1)/(2 mu R^2).

Parameters per state: logDe, logRe, (B) Te, c0..c(n-1), alpha0..alpha(m-1). Derivatives with respect to c and
alpha are analytic; logDe, logRe and Te by central differences.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.polynomial.chebyshev import chebval, chebvander


LONG_RANGE = {"X": {6: 1.48e6, 8: 3.86e7, 10: 1.0e8},
              "B": {5: 3.161e5, 6: 1.506e6, 8: 2.48e7, 10: 4.20e8}}


@dataclass(frozen=True)
class Config:
    """Structure of the two curves (the release's model.json "config"; only the MLR family is supported)."""
    family: str = "mlr"
    nx: int = 14
    nb: int = 20
    alpha_x: int = 2
    alpha_b: int = 6
    px: int = 5
    pb: int = 6
    qx: int = 3
    qb: int = 4
    ref_x: float = 1.
    ref_b: float = 1.5
    mapped_x: bool = True
    mapped_alpha: bool = True
    #: power p of the B isotope shift (1 - rho)(1 - y)^p (iso0 + iso1 y), y = (R - Re)/(R + Re), seen only by the
    #: other isotopologues (isotopologues.py's form); None: no shift terms (the release)
    iso_shift: int | None = None

    @classmethod
    def from_spec(cls, spec: dict):
        """From a release's config, refusing options this implementation does not carry."""
        unsupported = dict(adiabatic_b=0, damping=False, free_lr=False, fit_reference=False, correct_b_c10=False,
                           positive_rotation=False, family_x=None, family_b=None, denominator_x=0, denominator_b=0)
        for key, default in unsupported.items():
            if spec.get(key, default) != default:
                raise ValueError(f"config option {key}={spec[key]} is not supported")
        if spec.get("family", "mlr") != "mlr":
            raise ValueError(f"family {spec['family']} is not supported")
        return cls(**{k: spec[k] for k in cls.__dataclass_fields__ if k in spec})

    def order(self, state):
        return self.nx if state == "X" else self.nb

    def powers(self, state):
        return (self.px, self.qx) if state == "X" else (self.pb, self.qb)

    def ref(self, state):
        return self.ref_x if state == "X" else self.ref_b

    def nalpha(self, state):
        return self.alpha_x if state == "X" else self.alpha_b


def long_range(state, values=None):
    """The inherited dispersion coefficients C_n (cm-1 A^n), each scaled by exp(logC_n) where `values` has one."""
    values = values or {}
    return {n: c*np.exp(values.get(f"logC{n}", 0.)) for n, c in LONG_RANGE[state].items()}


def initial(config: Config):
    """Generic starting values, bounds and scales: (values per state, [(state, name)], lower, upper, scales)."""
    values, names, lower, upper, scales = {}, [], [], [], []
    for state in ("X", "B"):
        v = dict(logDe=np.log(12547.335 if state == "X" else 4381.245), logRe=np.log(2.6664 if state == "X" else 3.0263))
        lo = dict(logDe=np.log(9000. if state == "X" else 3000.), logRe=np.log(2.5 if state == "X" else 2.8))
        hi = dict(logDe=np.log(15000. if state == "X" else 6500.), logRe=np.log(2.85 if state == "X" else 3.25))
        sc = dict(logDe=.01, logRe=.005)
        if state == "B":
            v["Te"], lo["Te"], hi["Te"], sc["Te"] = 15769.06, 15300., 16200., 50.
        for i in range(config.order(state)):
            v[f"c{i}"], lo[f"c{i}"], hi[f"c{i}"], sc[f"c{i}"] = ((-1.25 if state == "X" else -.45) if i == 0 else 0.,
                                                                  -100., 100., .1)
        for i in range(config.nalpha(state)):
            v[f"alpha{i}"], lo[f"alpha{i}"], hi[f"alpha{i}"], sc[f"alpha{i}"] = 0., -.1, .1, 1e-4
        if state == "B" and config.iso_shift is not None:
            for i in range(2):
                v[f"iso{i}"], lo[f"iso{i}"], hi[f"iso{i}"], sc[f"iso{i}"] = 0., -100., 100., .1
        values[state] = v
        for key in v:
            names.append((state, key)); lower.append(lo[key]); upper.append(hi[key]); scales.append(sc[key])
    return values, names, np.array(lower), np.array(upper), np.array(scales)


class Curve:
    """The MLR potential and rotational function of one state at given parameter values."""

    def __init__(self, state, config: Config, values: dict):
        self.state, self.config, self.values = state, config, dict(values)

    def _coordinates(self, values, R):
        cfg, s = self.config, self.state
        re = np.exp(values["logRe"])
        p, q = cfg.powers(s)
        ref = re*cfg.ref(s)
        y = (R**p-ref**p)/(R**p+ref**p)
        z = (R**q-ref**q)/(R**q+ref**q)
        if cfg.mapped_x and s == "X":
            z = (z-.2)/.6
        yp = (R**p-re**p)/(R**p+re**p)
        return re, y, z, yp

    def potential_values(self, values, R):
        R = np.asarray(R, float)
        re, y, z, yp = self._coordinates(values, R)
        de = np.exp(values["logDe"])
        C = long_range(self.state, values)
        u = sum(c/R**n for n, c in C.items())
        ue = float(sum(c/re**n for n, c in C.items()))
        exponent = chebval(z, [values[f"c{i}"] for i in range(self.config.order(self.state))])
        beta = y*np.log(2*de/ue)+(1-y)*exponent
        shape = u/ue*np.exp(np.clip(-beta*yp, -100, 100))
        return values.get("Te", 0.)+de*(1-shape)**2

    def asymptote(self, values):
        return float(self.potential_values(values, np.array([400.]))[0])

    def _rotation(self, values, R):
        re = np.exp(values["logRe"])
        w = (R-re)/(R+re)
        if self.config.mapped_alpha:
            w = (w-.3)/.4
        return 2*re/(R+re), w

    def curves(self, values, R, ratio):
        """V(R) and f(R); `ratio` = mu(127I2)/mu scales the rotational correction for other isotopologues."""
        R = np.asarray(R, float)
        f = np.ones_like(R)
        m = self.config.nalpha(self.state)
        if m:
            envelope, w = self._rotation(values, R)
            f = 1+ratio*envelope*chebval(w, [values[f"alpha{i}"] for i in range(m)])
        V = self.potential_values(values, R)
        if ratio != 1. and self._shifted():
            V = V+(1-ratio)*self._iso_basis(values, R)@[values["iso0"], values["iso1"]]
        return V, f

    def _shifted(self):
        return self.state == "B" and self.config.iso_shift is not None

    def _iso_basis(self, values, R):
        re = np.exp(values["logRe"])
        y = (R-re)/(R+re)
        return np.stack([(1-y)**self.config.iso_shift, (1-y)**self.config.iso_shift*y], axis=-1)

    def derivatives(self, values, R, ratio, names):
        """dV/dp and df/dp on R for the parameters `names` (columns)."""
        R = np.asarray(R, float)
        re, y, z, yp = self._coordinates(values, R)
        de = np.exp(values["logDe"])
        n = self.config.order(self.state)
        C = long_range(self.state, values)
        u = sum(c/R**k for k, c in C.items())
        ue = float(sum(c/re**k for k, c in C.items()))
        beta = y*np.log(2*de/ue)+(1-y)*chebval(z, [values[f"c{i}"] for i in range(n)])
        argument = -beta*yp
        shape = u/ue*np.exp(np.clip(argument, -100, 100))
        factor = 2*de*(1-shape)*shape*yp*(abs(argument) < 100)
        coefficient = (factor*(1-y))[:, None]*chebvander(z, n-1)
        de_derivative = de*(1-shape)**2+factor*y
        envelope, w = self._rotation(values, R)
        dV, df = np.zeros((len(R), len(names))), np.zeros((len(R), len(names)))
        for col, name in enumerate(names):
            if name.startswith("c") and name[1:].isdigit():
                dV[:, col] = coefficient[:, int(name[1:])]
            elif name == "logDe":
                dV[:, col] = de_derivative
            elif name == "Te":
                dV[:, col] = 1.
            elif name.startswith("alpha"):
                df[:, col] = ratio*envelope*chebval(w, [0.]*int(name[5:])+[1.])
            elif name.startswith("iso"):
                if ratio != 1.:
                    dV[:, col] = (1-ratio)*self._iso_basis(values, R)[:, int(name[3:])]
            else:
                h = max(abs(values[name])*1e-7, 1e-6)
                plus, minus = dict(values), dict(values)
                plus[name] += h
                minus[name] -= h
                vp, fp = self.curves(plus, R, ratio)
                vm, fm = self.curves(minus, R, ratio)
                dV[:, col], df[:, col] = (vp-vm)/(2*h), (fp-fm)/(2*h)
        return dV, df

    def log_shape(self, values, R):
        """-ln(Morse factor) = base + basis @ c, linear in the exponent coefficients: (base, basis)."""
        R = np.asarray(R, float)
        re, y, z, yp = self._coordinates(values, R)
        de = np.exp(values["logDe"])
        C = long_range(self.state, values)
        u = sum(c/R**k for k, c in C.items())
        ue = float(sum(c/re**k for k, c in C.items()))
        return yp*y*np.log(2*de/ue)-np.log(u/ue), (yp*(1-y))[:, None]*chebvander(z, self.config.order(self.state)-1)

    def log_slope(self, values, R, names=None):
        """d/dR of -ln(Morse factor); positive everywhere for a single well, even where V itself is too close to
        its asymptote for differences of V to show it. With `names`, also its parameter derivatives."""
        R = np.asarray(R, float)
        step = 1e-5
        n = self.config.order(self.state)
        bp, ap = self.log_shape(values, R+step)
        bm, am = self.log_shape(values, R-step)
        basis = (ap-am)/(2*step)
        slope = (bp-bm)/(2*step)+basis@np.array([values[f"c{i}"] for i in range(n)])
        if names is None:
            return slope
        jac = np.zeros((len(R), len(names)))
        for i, name in enumerate(names):
            if name.startswith("c") and name[1:].isdigit():
                jac[:, i] = basis[:, int(name[1:])]
            elif name == "Te" or name.startswith("alpha"):
                continue
            else:
                h = max(abs(values[name])*1e-6, 1e-5)
                plus, minus = dict(values), dict(values)
                plus[name] += h
                minus[name] -= h
                jac[:, i] = (self.log_slope(plus, R)-self.log_slope(minus, R))/(2*h)
        return slope, jac

