"""The model's named parameters (decision 9): Be(X), Be(B), T0 and D0(X) stored in place of logRe(X), logRe(B), Te
and logDe(X), and in a continued model D0(X) of the whole model in place of the asymptote shift dDe. An exact change
of variables: a model file in named coordinates loads to the same coefficients, and predicts the same.

  Be     = hbar^2/(2 mu Re^2) (1 + q(Re))       q the rotational function, whose argument is 0 at Re: q(Re) does not
                                                 depend on Re, and Re follows in closed form
  T0     = E_B(0, 0) - E_X(0, 0)
  D0(X)  = De(X) - E_X(0, 0)                     (core); + dDe de_unit in a continued model, whose shift moves both
                                                 asymptotes (the levels v = 0 do not reach the tails)
  De(B)  = De(X) + the atomic splitting - Te     the B limit is D0(X) + the atomic splitting, exactly ("asymptote":
                                                 "tied"); a tied model stores no De(B)

E(v, J) are the runtime's own levels (research._Potential on the core file's grids), so a file written by `to_named`
and read back by `to_native` returns its coefficients to the solver's precision (tests/test_named.py)."""
from __future__ import annotations

import copy

import numpy as np
from numpy.polynomial.chebyshev import chebval

from .constants import ATOMIC_SPLITTING, HBAR2_2U, reduced_mass

NAMES = {"X": ("Be", "D0"), "B": ("Be", "T0")}
TOL = 1e-10                     # cm-1 (3 Hz), the inversions' convergence


def q_re(spec, state):
    """The rotational function at Re."""
    cfg, values = spec["config"], spec["values"][state]
    count = cfg["alpha_x" if state == "X" else "alpha_b"]
    if not count:
        return 0.
    return float(chebval(-.75 if cfg["mapped_alpha"] else 0., [values[f"alpha{i}"] for i in range(count)]))


def be(spec, state):
    re = np.exp(spec["values"][state]["logRe"])
    return HBAR2_2U/reduced_mass("127I2")/re**2*(1+q_re(spec, state))


def log_re(spec, state, be_value):
    return .5*np.log(HBAR2_2U*(1+q_re(spec, state))/(reduced_mass("127I2")*be_value))


def ground(spec, state):
    """E(v=0, J=0) of the core potential (cm-1), on the file's grids."""
    from .bspline import BSplineSolver
    from .research import _Potential
    solver = BSplineSolver(_Potential(state, spec), reduced_mass("127I2"), **spec["grids"][state])
    return float(solver.levels(0)[0])


def tie(spec):
    """De(B) from the tie: logDe(B) = log(De(X) + atomic splitting - Te)."""
    v = spec["values"]
    v["B"]["logDe"] = float(np.log(np.exp(v["X"]["logDe"])+ATOMIC_SPLITTING-v["B"]["Te"]))


def named(spec):
    """{state: {name: value}} of a core spec in native coordinates."""
    ex = ground(spec, "X")
    return {"X": dict(Be=be(spec, "X"), D0=float(np.exp(spec["values"]["X"]["logDe"])-ex)),
            "B": dict(Be=be(spec, "B"), T0=ground(spec, "B")-ex)}


def to_named(spec):
    """A tied core spec in named coordinates (a copy): Be, D0 for X; Be, T0 for B; no logRe, logDe(X), Te, logDe(B)."""
    if spec.get("asymptote") != "tied":
        raise ValueError("named coordinates need the B limit tied to D0(X) + the atomic splitting")
    out = copy.deepcopy(spec)
    for state, values in named(spec).items():
        for k in ("logRe", "logDe", "Te"):
            out["values"][state].pop(k, None)
        out["values"][state].update(values)
    out["coordinates"] = "named"
    return out


def to_native(spec):
    """A core spec in native coordinates (a copy; one in native coordinates is returned unchanged), its named
    values kept under "named"."""
    if spec.get("coordinates") != "named":
        return spec
    out = copy.deepcopy(spec)
    v = out["values"]
    keep = {s: {k: v[s].pop(k) for k in NAMES[s]} for s in ("X", "B")}
    for s in ("X", "B"):
        v[s]["logRe"] = float(log_re(out, s, keep[s]["Be"]))
    # De(X): D0 = De - E_X(0,0)(De), by the secant method from E_X(0,0) at the D0 + 107 cm-1 guess
    def miss(de):
        v["X"]["logDe"] = float(np.log(de))
        return de-ground(out, "X")-keep["X"]["D0"]
    a, b = keep["X"]["D0"]+100., keep["X"]["D0"]+110.
    fa, fb = miss(a), miss(b)
    for _ in range(50):
        if abs(fb) < TOL or fb == fa:
            break
        a, b, fa = b, b-fb*(b-a)/(fb-fa), fb
        fb = miss(b)
    else:
        raise ValueError("D0(X) did not invert")
    miss(b)
    ex = ground(out, "X")
    # Te: T0 = E_B(0,0)(Te, De(B)(Te)) - E_X(0,0); E_B(0,0) - Te depends on Te only through De(B), weakly
    v["B"]["Te"] = keep["B"]["T0"]+ex
    for _ in range(50):
        tie(out)
        step = keep["B"]["T0"]-(ground(out, "B")-ex)
        v["B"]["Te"] += step
        if abs(step) < TOL:
            break
    else:
        raise ValueError("T0 did not invert")
    tie(out)
    out["named"] = keep
    out["coordinates"] = "native"
    return out


def continued_to_named(spec, core_named):
    """A continued model's spec (a copy) with dDe replaced by D0(X) of the whole model: D0(core) + dDe de_unit."""
    out = copy.deepcopy(spec)
    dde = out["values"]["X"].pop("dDe")
    out["values"]["B"].pop("dDe", None)
    out["values"]["X"]["D0"] = core_named["X"]["D0"]+dde*out["structure"]["de_unit"]
    out["coordinates"] = "named"
    return out


def continued_to_native(spec, core_named):
    """The inverse of continued_to_named (a copy; native specs unchanged)."""
    if spec.get("coordinates") != "named":
        return spec
    out = copy.deepcopy(spec)
    d0 = out["values"]["X"].pop("D0")
    dde = (d0-core_named["X"]["D0"])/out["structure"]["de_unit"]
    out["values"]["X"]["dDe"] = out["values"]["B"]["dDe"] = float(dde)
    out["named"] = {"X": dict(D0=d0)}
    out["coordinates"] = "native"
    return out


def covariance(spec, names, factor, de_unit=None, step=1e-6):
    """Values, 1-sigma and correlations of the named parameters of a native core spec, given the fit's covariance
    Cov = F F^T over `names` ("state:name"; the core's values, its inherited "logC<n>", the shared "X:dDe" with
    `de_unit` for a continued model's D0): dN/dp by central differences of `named`, through the tie."""
    keys = [(s, k) for s in ("X", "B") for k in NAMES[s]]
    base = named(spec)
    value = np.array([base[s][k] for s, k in keys])
    if de_unit is not None:
        value[keys.index(("X", "D0"))] += spec["values"]["X"]["dDe"]*de_unit
    J = np.zeros((len(keys), len(names)))
    for j, name in enumerate(names):
        s, n = name.split(":", 1)
        if n == "dDe" and de_unit is not None:
            J[keys.index(("X", "D0")), j] = de_unit
            continue
        if s not in ("X", "B") or not (n in spec["values"][s] or n.startswith("logC")):
            continue
        cols = []
        for sign in (1, -1):
            p = copy.deepcopy(spec)
            if n.startswith("logC"):
                p["inherited_potential_coefficients"][s][n[4:]] *= np.exp(sign*step)
                h = step
            else:
                h = step*max(1., abs(p["values"][s][n]))
                p["values"][s][n] += sign*h
            if p.get("asymptote") == "tied":
                tie(p)
            got = named(p)
            cols.append(np.array([got[a][b] for a, b in keys]))
        J[:, j] = (cols[0]-cols[1])/(2*h)
    C = (J@np.asarray(factor))@(J@np.asarray(factor)).T
    sigma = np.sqrt(np.diag(C))
    label = [f"{k}({s})" if k != "T0" else "T0" for s, k in keys]
    return dict(names=label, values=dict(zip(label, value.tolist())), sigma=dict(zip(label, sigma.tolist())),
                correlation=(C/np.outer(sigma, sigma)).tolist(),
                note="linear propagation of the fit's statistical covariance (uncertainty.json) through the exact map")
