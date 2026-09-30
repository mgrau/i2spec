"""Feasibility check: rebuild the Knöckel, Bodermann & Tiemann (EPJD 28, 199, 2004) B/X
potentials from their Table 4, then compare the rovibronic line positions against
(1) their local NIR Dunham model (Table 1, hyperfine-free, <200 kHz 1σ), and
(2) known absolute hyperfine-component frequencies (coarse: hyperfine offsets ~±500 MHz).
Solver: sinc-DVR (Colbert–Miller), checked for grid convergence.
"""
import os
import sys
import numpy as np
from numpy.polynomial import polynomial as P
from scipy.linalg import eigh
from scipy import constants as K

C_MHZ = K.c / 1e4  # MHz per cm^-1 (29979.2458)
# hbar^2/(2 u Å^2) in cm^-1
HB2 = K.hbar**2 / (2 * K.atomic_mass * 1e-20) / (K.h * K.c * 100)
M127 = float(sys.argv[1]) if len(sys.argv) > 1 else 126.9044719  # AME2020 atomic mass, u
MU = M127 / 2
SIGN_AO = +1  # sign of the A_O exponential in the long-range branch (checked below)

B = dict(
    Rm=3.02669183, b=-0.75,
    a=[15769.0678546, 0.888929155287492279e1, 0.850235519686038970e4, 0.695917336737764890e4,
       0.994920182775927628e3, -0.556131449320786032e4, -0.107555928627667199e5,
       -0.169188457283073476e5, -0.153304207986605579e5, 0.433363227519699431e5,
       0.625822290457034105e5, -0.389266247379213921e6, -0.432218520110222744e6,
       0.272340536816511489e7, 0.321856089579221839e7, -0.127338456382425241e8,
       -0.172753342751943097e8, 0.425256255365765169e8, 0.712696682021539360e8,
       -0.913057098990224898e8, -0.210868067425342679e9, 0.926719291015004069e8,
       0.413233196395634711e9, 0.608534295451399982e8, -0.487970270485778809e9,
       -0.305239338118491709e9, 0.269185593565317631e9, 0.349461001962127566e9,
       0.277101477602953948e8, -0.133371646903035790e9, -0.765550455925225466e8,
       -0.136649829213649314e8],
    RI=2.647, AI=0.19603113e5, BI=0.13702730e1,
    RO=4.9, AO=0.44452892e2, BO=0.71764654, De=20150.317,
    C={5: 0.3161e6, 6: 0.1506e7, 8: 0.2480e8, 10: 0.0420e10},
    alpha=[-0.257024561963860657e-3, 0.304354657912661525e-4, -0.824779020512065772e-3,
           -0.321720674688543535e-3, -0.484148886990049172e-3],
)
X = dict(
    Rm=2.66638233, b=-0.60,
    a=[0.0, 0.947051006049874633e1, 0.492602562595701165e5, 0.213682631498861410e5,
       -0.181879788791815008e5, -0.499045801469034486e5, -0.702742764701535925e5,
       -0.882312608143186953e5, -0.102331582920833287e5, 0.339528118110142939e6,
       0.628206506769069820e5, -0.717230249674846511e7, 0.855012363130142097e6,
       0.40925020305406115e8],
    RI=2.40, AI=0.45807577e4, BI=0.90114269e1,
    RO=3.30, AO=0.10523133e4, BO=0.12830221e1, De=12547.340,
    C={6: 0.148e7, 8: 0.386e8, 10: 0.100e9},
    alpha=None,
)
assert len(B["a"]) == 32 and len(X["a"]) == 14

# Salumbides, Eikema, Ubachs, Hollenstein, Knöckel, Tiemann, EPJD 47, 171 (2008), Table 1.
# Long-range branch there: V = De - sum Cn/R^n - A_O exp(-B_O (R - R_O))  -> sign = -1.
B08 = dict(
    Rm=3.02669183, b=-0.75,
    a=[15769.067837, 0.813930e1, 0.850290979e4, 0.695677594e4, 0.995264129e3, -0.556338644e4,
       -0.1074389661e5, -0.1685419562e5, -0.1533164369e5, 0.4315045697e5, 0.6251096110e5,
       -0.38914630400e6, -0.43219553147e6, 0.27235279905e7, 0.32186808611e7,
       -0.127338590547e8, -0.172752514874e8, 0.425255265406e8, 0.712695731657e8,
       -0.9130578539938e8, -0.2108682403457e9, 0.9267186544061e8, 0.4132330557779e9,
       0.6085337864116e8, -0.4879703383624e9, -0.3052394100029e9, 0.2691855865112e9,
       0.3494609136897e9, 0.2771015696362e8, -0.1333716748717e9, -0.7655507971445e8,
       -0.1366495310913e8],
    RI=2.647, AI=0.19605269e5, BI=0.13736027e1,
    RO=4.9, AO=-0.8105289769e1, BO=0.9347221716e1, De=20150.317, sign=-1,
    C={5: 0.3161e6, 6: 0.1506e7, 8: 0.2480e8, 10: 0.4200e9},
    alpha=[-0.23586e-3, 0.23445e-4, -0.43703e-3, -0.39630e-3, -0.10191e-2, -0.16462e-3],
)
X08 = dict(
    Rm=2.66638233, b=-0.60,
    a=[0.0, 0.155389e0, 0.492601477e5, 0.213629963e5, -0.18173865e5, -0.49468474e5,
       -0.71353991e5, -0.10110925e6, -0.7156108e5, 0.4710538e6, -0.2904654e6,
       -0.73300648e7, 0.22578115e7, 0.38250822e8],
    RI=2.40, AI=0.45809329e4, BI=0.90168320e1,
    RO=3.30, AO=0.1053614048e4, BO=0.1259030368e1, De=12547.340, sign=-1,
    C={6: 0.148e7, 8: 0.386e8, 10: 0.10e9},
    alpha=None,
)
assert len(B08["a"]) == 32 and len(X08["a"]) == 14
if os.environ.get("MODEL") == "2008":
    B, X = B08, X08


def xvar(R, p):
    return (R - p["Rm"]) / (R + p["b"] * p["Rm"])


def v_series(R, p):
    return P.polyval(xvar(R, p), p["a"])


def v_long(R, p, sign=None):
    sign = p.get("sign", SIGN_AO) if sign is None else sign
    return p["De"] - sum(c / R**n for n, c in p["C"].items()) + sign * p["AO"] * np.exp(-p["BO"] * (R - p["RO"]))


def V(R, p):
    R = np.asarray(R, float)
    with np.errstate(all="ignore"):
        out = v_series(R, p)
        out = np.where(R < p["RI"], p["AI"] * np.exp(-p["BI"] * (R - p["RI"])), out)
        out = np.where(R > p["RO"], v_long(R, p), out)
    return out


def alpha(R, p):
    """Effective nonadiabatic BOC for the reference isotopologue (mu = mu_ref)."""
    if p["alpha"] is None:
        return 0.0
    return 2 * p["Rm"] / (R + p["Rm"]) * P.polyval(xvar(R, p), p["alpha"])


def dv_series(R, p):
    Rm, b = p["Rm"], p["b"]
    return P.polyval(xvar(R, p), P.polyder(p["a"])) * Rm * (1 + b) / (R + b * Rm) ** 2


def enforce_continuity(p):
    """Re-derive A_I, B_I, A_O, B_O from value+slope continuity of the series (as the
    paper says was done after each modification of the central part)."""
    p["tab"] = {k: p[k] for k in ("AI", "BI", "AO", "BO")}
    RI, RO = p["RI"], p["RO"]
    p["AI"] = v_series(RI, p)
    p["BI"] = -dv_series(RI, p) / p["AI"]
    disp = p["De"] - sum(c / RO**n for n, c in p["C"].items())
    ddisp = sum(n * c / RO ** (n + 1) for n, c in p["C"].items())
    p["AO"] = v_series(RO, p) - disp
    p["BO"] = -(dv_series(RO, p) - ddisp) / p["AO"]
    p["sign"] = +1  # A_O derived with the "+" convention


CONTINUITY = os.environ.get("TAB_EXT") is None  # set TAB_EXT=1 to use the printed table values
if CONTINUITY:
    enforce_continuity(B)
    enforce_continuity(X)


def continuity_report():
    for name, p in (("B", B), ("X", X)):
        tab = p.get("tab", {k: p[k] for k in ("AI", "BI", "AO", "BO")})
        print(f"{name}: used  A_I={p['AI']:.6f} B_I={p['BI']:.6f} A_O={p['AO']:.6f} B_O={p['BO']:.6f}")
        print(f"{name}: table A_I={tab['AI']:.6f} B_I={tab['BI']:.6f} A_O={tab['AO']:.6f} B_O={tab['BO']:.6f}")
        for R in (p["RI"], p["RO"]):
            print(f"   R={R}: V(R-)={V(R - 1e-9, p):.6f} V(R+)={V(R + 1e-9, p):.6f}")


class Solver:
    def __init__(self, p, rmin, rmax, dr, nlev):
        self.p, self.nlev = p, nlev
        self.R = np.arange(rmin, rmax + dr / 2, dr)
        n = len(self.R)
        i = np.arange(n)
        d = i[:, None] - i[None, :]
        with np.errstate(divide="ignore"):
            T = 2.0 * (-1.0) ** d / (d.astype(float) ** 2)
        T[i, i] = np.pi**2 / 3
        self.T = HB2 / MU * T / dr**2
        self.V = V(self.R, p)
        self.cent = HB2 / MU * (1 + alpha(self.R, p)) / self.R**2
        self.cache = {}

    def levels(self, J):
        if J not in self.cache:
            H = self.T + np.diag(self.V + J * (J + 1) * self.cent)
            self.cache[J] = eigh(H, eigvals_only=True, subset_by_index=[0, self.nlev - 1])
        return self.cache[J]


def build(scale=1.0):
    return (Solver(B, 2.35, 7.0, 0.005 * scale, 50), Solver(X, 2.10, 4.0, 0.004 * scale, 20))


def line(sB, sX, vp, vpp, branch, Jpp):
    Jp = Jpp + 1 if branch == "R" else Jpp - 1
    return sB.levels(Jp)[vp] - sX.levels(Jpp)[vpp]


def dunham_nir(vpp, branch, Jpp):
    """Knöckel 2004 Table 1 local Dunham model, bands 0-v'' (v''=12..17)."""
    Jp = Jpp + 1 if branch == "R" else Jpp - 1
    xp, xpp = Jp * (Jp + 1), Jpp * (Jpp + 1)
    S = [15832.7358035, 0.289260041202e-1, -0.62353693595e-8, -0.22689095e-14, -0.4269941e-20]
    up = S[0] + sum(S[k] * xp**k for k in range(1, 5))
    SX = {(1, 0): 0.21479058554e3, (2, 0): -0.6349014655, (4, 0): -0.31552791e-4,
          (0, 1): 0.373816471278e-1, (1, 1): -0.1166181096942e-3, (2, 1): -0.996620979e-7,
          (3, 1): -0.895649619e-8, (0, 2): -0.4608461804e-8, (1, 2): -0.1196947480e-10,
          (2, 2): -0.9829624898e-12, (0, 3): -0.36003776e-15, (1, 3): -0.33513423e-16,
          (0, 4): 0.6778317e-21, (1, 4): -0.1441870e-21}
    lo = sum(c * (vpp + 0.5) ** l * xpp**k for (l, k), c in SX.items())
    return up - lo


# Absolute frequencies of single hyperfine components (kHz). Coarse check only:
# component offsets from the rovibronic centre are up to several hundred MHz.
ABS = [
    ("R(56)32-0 a10, 532 nm (BIPM)", 32, 0, "R", 56, 563260223513),
    ("R(127)11-5 a16 (f), 633 nm (BIPM)", 11, 5, "R", 127, 473612353604),
    ("R(37)16-1 a1, 578 nm (Hong 2009)", 16, 1, "R", 37, 518304551833),
]

if __name__ == "__main__":
    print(f"hbar^2/2u = {HB2:.9f} cm^-1 Å^2 ; m(127I) = {M127} u ; mu = {MU} u ; c = {C_MHZ} MHz/cm^-1")
    continuity_report()
    sB, sX = build(1.0)
    sB2, sX2 = build(0.8)
    e0x, e0b = sX.levels(0)[0], sB.levels(0)[0]
    print(f"\nE_X(0,0) = {e0x:.6f} cm^-1  (GLL 1991 ZPE: De-D0 = 107.101)")
    print(f"T0 = E_B(0,0)-E_X(0,0) = {e0b - e0x:.6f} cm^-1  (GLL 1991: 15724.587)")
    print(f"grid check: dE_X(0,0) = {(sX2.levels(0)[0] - e0x) * C_MHZ * 1e3:.3f} kHz, "
          f"dE_B(0,0) = {(sB2.levels(0)[0] - e0b) * C_MHZ * 1e3:.3f} kHz, "
          f"dE_B(v=45,J=0) = {(sB2.levels(0)[45] - sB.levels(0)[45]) * C_MHZ * 1e3:.3f} kHz")

    print("\nNIR check: potential model minus local Dunham model (MHz), R(J'') lines")
    Js = [20, 50, 80, 100, 120, 150, 180, 210]
    print("band  " + "".join(f"{J:>9d}" for J in Js))
    for vpp in range(12, 18):
        row = [(line(sB, sX, 0, vpp, "R", J) - dunham_nir(vpp, "R", J)) * C_MHZ for J in Js]
        print(f"0-{vpp:<3d}" + "".join(f"{d:9.2f}" for d in row))
    J = 100
    dconv = (line(sB2, sX2, 0, 15, "R", J) - line(sB, sX, 0, 15, "R", J)) * C_MHZ * 1e3
    print(f"grid convergence of R(100) 0-15: {dconv:.3f} kHz")

    print("\nAbsolute components (coarse check): model rovibronic minus measured component")
    for name, vp, vpp, br, Jpp, fk in ABS:
        nu = line(sB, sX, vp, vpp, br, Jpp) * C_MHZ
        print(f"  {name:36s} model {nu:18.3f} MHz   meas {fk / 1e3:18.3f} MHz   diff {nu - fk / 1e3:9.2f} MHz")
