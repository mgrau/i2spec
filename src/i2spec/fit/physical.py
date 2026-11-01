"""The model in physical constants: an exact change of coordinates for the parameters spectroscopy names.

Per state, four coefficients of the core curves are exchanged for mechanical constants of the potential at Re
(Dunham, Phys. Rev. 41, 721 (1932)); the predictions do not change, because the map is inverted exactly:

  logRe, c0, c1, c2   <->   Be, omega_e, alpha_e, omega_e x_e

  V(R) - Te = A2 xi^2 (1 + a1 xi + a2 xi^2 + ...),   xi = (R - Re)/Re
  Be0 = hbar^2/(2 mu Re^2),   Be = Be0 (1 + q(Re))            (the rotational function q = f - 1 included)
  omega_e = 2 sqrt(A2 Be0),   omega_e x_e = -Y20 = -(3 Be0/2)(a2 - 5 a1^2/4)
  alpha_e = -Y11 = -Be (2 Be0/omega_e) (b2 - 3 a1 b1/2),   b1 = -2 + Re q'/(1+q),   b2 = 3 - 2 Re q'/(1+q) + Re^2 q''/(2(1+q))

(for q = 0, alpha_e = -(6 Be0^2/omega_e)(1 + a1)). The expansion is exact (jet arithmetic on the MLR formula), so
from_physical inverts to machine precision and the levels do not move (tests/test_fit_physical.py). The other
coefficients (c3..., the rotational alpha_i, De, Te) stay as they are.

What the data determine: Be to 3e-7 relative, but not the curvature at the very bottom of the well. In the core's
14- and 20-term exponents, omega_e, omega_e x_e and alpha_e at Re move along directions the levels barely see to first
order (the levels' own Dunham constants are fixed to 2e-5, 7e-6 and 1e-9 cm-1), so their linear sigma is large and
not meaningful; Re trades with q(R).
"""
from __future__ import annotations

import numpy as np
from numpy.polynomial import Chebyshev, Polynomial

from ..constants import HBAR2_2U, reduced_mass
from .curves import Curve, long_range

NATIVE = ("logRe", "c0", "c1", "c2")
PHYSICAL = ("Be", "omega_e", "alpha_e", "omega_e_x_e")
H, DEGREE, NODES = .03, 14, 64          # the Taylor expansion about Re: Chebyshev fit on Re (1 +- H)


def _taylor(f, re):
    """Power-series coefficients of f(Re (1 + xi)) in xi, by a Chebyshev fit (a check of the exact `expansion`)."""
    t = np.cos(np.pi*(np.arange(NODES)+.5)/NODES)
    p = Chebyshev.fit(t, f(re*(1+H*t)), DEGREE, domain=[-1, 1]).convert(kind=Polynomial).coef
    return p/H**np.arange(len(p))


class Jet:
    """A truncated power series in xi (coefficients c[0..N-1]): exact Taylor coefficients of the MLR formula."""

    def __init__(self, c):
        self.c = np.asarray(c, float)

    @staticmethod
    def _j(x, n):
        if isinstance(x, Jet):
            return x
        c = np.zeros(n)
        c[0] = x
        return Jet(c)

    def __add__(self, o):
        return Jet(self.c+self._j(o, len(self.c)).c)

    __radd__ = __add__

    def __sub__(self, o):
        return Jet(self.c-self._j(o, len(self.c)).c)

    def __rsub__(self, o):
        return Jet(self._j(o, len(self.c)).c-self.c)

    def __neg__(self):
        return Jet(-self.c)

    def __mul__(self, o):
        o = self._j(o, len(self.c))
        return Jet(np.convolve(self.c, o.c)[:len(self.c)])

    __rmul__ = __mul__

    def __truediv__(self, o):
        o = self._j(o, len(self.c))
        q = np.zeros(len(self.c))
        for k in range(len(q)):
            q[k] = (self.c[k]-np.dot(q[:k], o.c[k:0:-1]))/o.c[0]
        return Jet(q)

    def __rtruediv__(self, o):
        return self._j(o, len(self.c))/self

    def exp(self):
        b = np.zeros(len(self.c))
        b[0] = np.exp(self.c[0])
        for k in range(1, len(b)):
            b[k] = sum(j*self.c[j]*b[k-j] for j in range(1, k+1))/k
        return Jet(b)

    def log(self):
        b = np.zeros(len(self.c))
        b[0] = np.log(self.c[0])
        for k in range(1, len(b)):
            b[k] = (self.c[k]-sum(j*b[j]*self.c[k-j] for j in range(1, k))/k)/self.c[0]
        return Jet(b)

    def __pow__(self, n):
        return (self.log()*float(n)).exp()


def _chebval(z, coefficients):
    """Clenshaw's recurrence on jets."""
    b1, b2 = 0., 0.
    for c in coefficients[:0:-1]:
        b1, b2 = z*b1*2.-b2+c, b1
    return z*b1-b2+coefficients[0]


def expansion(curve: Curve, values, order=6):
    """Exact Taylor coefficients in xi = R/Re - 1 of V(R) - Te and of q(R) = f(R) - 1 (the MLR curve and rotational
    function of curves.Curve), by jet arithmetic."""
    cfg, s = curve.config, curve.state
    re, de = float(np.exp(values["logRe"])), float(np.exp(values["logDe"]))
    R = Jet(np.r_[re, re, np.zeros(order-1)])
    p, q = cfg.powers(s)
    ref = re*cfg.ref(s)
    y = (R**p-ref**p)/(R**p+ref**p)
    z = (R**q-ref**q)/(R**q+ref**q)
    if cfg.mapped_x and s == "X":
        z = (z-.2)/.6
    yp = (R**p-re**p)/(R**p+re**p)
    C = long_range(s, values)
    u = sum((R**n).__rtruediv__(c) for n, c in C.items())       # c is a numpy scalar: keep numpy out of the jets
    ue = float(sum(c/re**n for n, c in C.items()))
    exponent = _chebval(z, [values[f"c{i}"] for i in range(cfg.order(s))])
    beta = y*float(np.log(2*de/ue))+(1-y)*exponent
    shape = u/ue*(-(beta*yp)).exp()
    V = (1-shape)*(1-shape)*de
    m = cfg.nalpha(s)
    if m:
        w = (R-re)/(R+re)
        if cfg.mapped_alpha:
            w = (w-.3)/.4
        f = (2*re/(R+re))*_chebval(w, [values[f"alpha{i}"] for i in range(m)])
        qc = f.c
    else:
        qc = np.zeros(order+1)
    return V.c, qc


def mechanical(curve: Curve, values, isotopologue="127I2"):
    """The mechanical constants of one state's curves at Re (cm-1), and the expansion they come from."""
    re = float(np.exp(values["logRe"]))
    mu = reduced_mass(isotopologue)
    A, q = expansion(curve, values)
    a1, a2 = A[3]/A[2], A[4]/A[2]
    be0 = HBAR2_2U/mu/re**2
    be = be0*(1+q[0])
    we = 2*np.sqrt(A[2]*be0)
    b1 = -2+q[1]/(1+q[0])
    b2 = 3-2*q[1]/(1+q[0])+q[2]/(1+q[0])
    return dict(Re=re, Be0=be0, Be=be, omega_e=we, omega_e_x_e=-1.5*be0*(a2-1.25*a1*a1),
                alpha_e=-be*(2*be0/we)*(b2-1.5*a1*b1), a1=a1, a2=a2, A2=A[2], q_Re=q[0], slope=A[1])


def to_physical(curve, values):
    """values with logRe, c0, c1, c2 replaced by Be, omega_e, alpha_e, omega_e_x_e."""
    m = mechanical(curve, values)
    out = {k: z for k, z in values.items() if k not in NATIVE}
    out.update({k: m[k] for k in PHYSICAL})
    return out


def from_physical(curve, physical, start, tol=1e-14, iterations=30):
    """Native values whose mechanical constants are `physical`'s (Newton from `start`, the native values to begin
    with; columns by central differences), to `tol` relative."""
    x = np.array([start[k] for k in NATIVE], float)
    target = np.array([physical[k] for k in PHYSICAL])
    base = {k: z for k, z in physical.items() if k not in PHYSICAL}

    def phi(x):
        m = mechanical(curve, dict(base, **dict(zip(NATIVE, x))))
        return np.array([m[k] for k in PHYSICAL])

    for _ in range(iterations):
        f = phi(x)
        r = (f-target)/np.abs(target)
        if np.max(np.abs(r)) < tol:
            for _ in range(2):           # two more steps to the noise floor
                step = 1e-7*np.maximum(np.abs(x), 1e-3)
                J = np.column_stack([(phi(x+e*s)-phi(x-e*s))/(2*s) for e, s in zip(np.eye(4), step)])
                x = x-np.linalg.solve(J/np.abs(target)[:, None], (phi(x)-target)/np.abs(target))
            break
        step = 1e-7*np.maximum(np.abs(x), 1e-3)
        J = np.column_stack([(phi(x+e*s)-phi(x-e*s))/(2*s) for e, s in zip(np.eye(4), step)])/np.abs(target)[:, None]
        x = x-np.linalg.solve(J, r)
    else:
        raise ValueError("the physical constants could not be inverted")
    return dict(base, **dict(zip(NATIVE, x)))


def jacobian(curve, values, names):
    """d(physical constants)/d(native parameter) for the native parameters `names` (central differences)."""
    out = np.zeros((len(PHYSICAL), len(names)))
    for j, n in enumerate(names):
        h = 1e-6*max(1., abs(values[n]))
        plus, minus = dict(values), dict(values)
        plus[n] += h
        minus[n] -= h
        a, b = mechanical(curve, plus), mechanical(curve, minus)
        out[:, j] = [(a[k]-b[k])/(2*h) for k in PHYSICAL]
    return out
