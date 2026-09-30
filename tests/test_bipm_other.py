"""Hyperfine structure of the other iodine standards against the BIPM mise en pratique tables.

Model: Hannover 2008 potentials; Bodermann et al. (2002) hyperfine formulae for ¹²⁷I₂ with v' <= 43,
Salumbides et al. (2006) otherwise (docs/research/isotopologue-hyperfine.md).
Reference: BIPM, "Recommended values of standard frequencies", iodine at 515 nm (MEP 2005) and at 543,
576, 612, 633 and 640 nm (MEP 2003), retrieved 2026-09-14. URLs are in data/catalog/precision.yaml;
transcription notes and per-line residuals in docs/research/bipm-hyperfine-tables.md. The tables
number the ΔF = ΔJ components of each line by increasing frequency (a1, a2, ... or b1, ..., m1, ...),
as i2spec labels them; the letter only names the line, so every label below is written aN.
"""

import numpy as np
import pytest

from i2spec import RovibronicModel

# Each line: isotopologue, (v', v''), J'', branch; values [f(aN) - f(zero)] / MHz with u_c / MHz,
# where the zero is the table's reference component (sometimes on another line). s: the model scatter
# (kHz) achieved on the well-measured components, rounded up (see test_component_positions).
LINES = {
    "633 127I2 R(127) 11-5": dict(iso="127I2", v=(11, 5), J=127, branch="R", s=30, data={  # Table 1
        "a2": (-721.8, 0.5), "a3": (-697.8, 0.5), "a4": (-459.62, 0.01), "a5": (-431.58, 0.05),
        "a6": (-429.18, 0.05), "a7": (-402.09, 0.01), "a8": (-301.706, 0.005), "a9": (-292.693, 0.005),
        "a10": (-276.886, 0.005), "a11": (-268.842, 0.005), "a12": (-160.457, 0.005), "a13": (-138.892, 0.005),
        "a14": (-116.953, 0.005), "a15": (-13.198, 0.005), "a16": (0.0, 0.0), "a17": (13.363, 0.005),
        "a18": (26.224, 0.005), "a19": (144.114, 0.005), "a20": (152.208, 0.005), "a21": (161.039, 0.005)}),
    "633 127I2 P(33) 6-3": dict(iso="127I2", v=(6, 3), J=33, branch="P", s=40, data={  # Table 2, bN
        "a1": (-922.571, 0.008), "a2": (-895.064, 0.008), "a3": (-869.67, 0.01), "a4": (-660.50, 0.02),
        "a5": (-610.697, 0.008), "a6": (-593.996, 0.008), "a7": (-547.40, 0.02), "a8": (-487.074, 0.009),
        "a9": (-461.30, 0.03), "a10": (-453.21, 0.03), "a11": (-439.01, 0.01), "a12": (-347.354, 0.007),
        "a13": (-310.30, 0.01), "a14": (-263.588, 0.009), "a15": (-214.53, 0.02), "a16": (-179.312, 0.005),
        "a17": (-153.942, 0.005), "a18": (-118.228, 0.007), "a19": (-36.73, 0.01), "a20": (-21.980, 0.007),
        "a21": (0.0, 0.0)}),
    "515 127I2 P(13) 43-0": dict(iso="127I2", v=(43, 0), J=13, branch="P", s=3, data={  # Table 1 (kHz)
        "a1": (-131.770, 0.001), "a2": (-59.905, 0.001), "a3": (0.0, 0.0), "a4": (76.049, 0.001),
        "a5": (203.229, 0.005), "a6": (240.774, 0.005), "a7": (255.005, 0.001), "a8": (338.699, 0.005),
        "a9": (349.717, 0.005), "a10": (369.0, 1.0), "a11": (393.962, 0.002), "a12": (435.599, 0.003),
        "a13": (499.712, 0.005), "a14": (518.0, 1.0), "a15": (587.396, 0.002), "a16": (616.756, 0.005),
        "a17": (660.932, 0.005), "a18": (740.0, 1.0), "a19": (742.0, 1.0), "a20": (757.631, 0.010),
        "a21": (817.337, 0.005)}),
    "515 127I2 R(15) 43-0": dict(iso="127I2", v=(43, 0), J=15, branch="R", s=3, data={  # Table 2, bN (kHz)
        "a1": (0.0, 0.0), "a2": (69.739, 0.005), "a3": (129.155, 0.005), "a4": (217.0, 1.0),
        "a5": (335.828, 0.005), "a6": (368.0, 1.0), "a7": (396.442, 0.005), "a8": (471.0, 1.0),
        "a9": (472.0, 1.0), "a10": (500.627, 0.005), "a11": (525.207, 0.005), "a12": (566.287, 0.005),
        "a13": (630.782, 0.005), "a14": (658.178, 0.005), "a15": (725.166, 0.005), "a16": (739.394, 0.005),
        "a17": (791.673, 0.005), "a18": (865.523, 0.005), "a19": (874.840, 0.005), "a20": (892.895, 0.010),
        "a21": (947.278, 0.010)}),
    "515 127I2 R(98) 58-1": dict(iso="127I2", v=(58, 1), J=98, branch="R", s=5, data={  # Table 3, dN (kHz)
        "a1": (-413.488, 0.005), "a2": (-359.553, 0.005), "a3": (-194.521, 0.005), "a4": (-159.158, 0.005),
        "a5": (-105.769, 0.005), "a6": (0.0, 0.0), "a7": (172.200, 0.005), "a8": (200.478, 0.005),
        "a9": (225.980, 0.005), "a10": (253.0, 1.0), "a11": (254.0, 1.0), "a12": (314.131, 0.005),
        "a13": (426.691, 0.005), "a14": (481.574, 0.005), "a15": (510.246, 0.005)}),
    "543 127I2 R(12) 26-0": dict(iso="127I2", v=(26, 0), J=12, branch="R", s=25, data={  # Table 1, zero b10 of R(106)
        "a1": (-1162.24, 0.02), "a2": (-909.87, 0.02), "a3": (-900.11, 0.03), "a4": (-853.336, 0.005),
        "a5": (-848.131, 0.005), "a6": (-795.92, 0.01), "a7": (-752.382, 0.005), "a8": (-733.134, 0.005),
        "a9": (-679.420, 0.005), "a10": (-596.134, 0.005), "a11": (-485.61, 0.01), "a12": (-476.35, 0.01),
        "a13": (-423.23, 0.01), "a14": (-410.01, 0.01), "a15": (-305.910, 0.005)}),
    "543 127I2 R(106) 28-0": dict(iso="127I2", v=(28, 0), J=106, branch="R", s=45, data={  # Table 2, bN
        "a1": (-573.765, 0.005), "a2": (-320.462, 0.005), "a3": (-291.59, 0.01), "a4": (-282.143, 0.005),
        "a5": (-253.675, 0.005), "a6": (-172.693, 0.005), "a7": (-159.428, 0.005), "a8": (-127.760, 0.005),
        "a9": (-114.575, 0.005), "a10": (0.0, 0.0), "a11": (124.83, 0.01), "a12": (132.31, 0.01),
        "a13": (154.51, 0.01), "a14": (162.65, 0.01), "a15": (287.24, 0.01)}),
    "576 127I2 P(62) 17-1": dict(iso="127I2", v=(17, 1), J=62, branch="P", s=100, data={  # Table 1
        "a1": (0.0, 0.0), "a2": (275.03, 0.02), "a3": (287.05, 0.02), "a4": (292.57, 0.02), "a5": (304.26, 0.02),
        "a6": (416.67, 0.02), "a7": (428.51, 0.02), "a8": (440.17, 0.02), "a9": (452.30, 0.02),
        "a10": (579.43, 0.03), "a15": (869.53, 0.03)}),
    "612 127I2 R(47) 9-2": dict(iso="127I2", v=(9, 2), J=47, branch="R", s=110, data={  # Table 1
        "a1": (-357.16, 0.02), "a2": (-333.97, 0.01), "a3": (-312.46, 0.02), "a4": (-86.168, 0.007),
        "a5": (-47.274, 0.004), "a6": (-36.773, 0.003), "a7": (0.0, 0.0), "a8": (81.452, 0.003),
        "a9": (99.103, 0.003), "a10": (107.463, 0.005), "a11": (119.045, 0.006), "a12": (219.602, 0.006),
        "a13": (249.60, 0.01), "a14": (284.30, 0.01), "a15": (358.37, 0.03), "a16": (384.66, 0.01),
        "a17": (403.76, 0.02), "a18": (429.99, 0.02), "a19": (527.16, 0.02), "a20": (539.22, 0.02),
        "a21": (555.09, 0.02)}),
    "612 127I2 P(48) 11-3": dict(iso="127I2", v=(11, 3), J=48, branch="P", s=75, data={  # Table 2, bN; zero a7 of R(47)
        "a1": (-1034.75, 0.07), "a2": (-755.86, 0.05), "a3": (-748.28, 0.03), "a4": (-738.35, 0.04),
        "a5": (-731.396, 0.006), "a6": (-616.01, 0.03), "a7": (-602.42, 0.03), "a8": (-593.98, 0.01),
        "a9": (-579.91, 0.01), "a10": (-452.163, 0.005), "a11": (-316.6, 0.4), "a12": (-315.8, 0.4),
        "a13": (-297.42, 0.03), "a14": (-294.72, 0.03), "a15": (-160.318, 0.003)}),
    "612 127I2 R(48) 15-5": dict(iso="127I2", v=(15, 5), J=48, branch="R", s=25, data={  # Table 3, cN; zero a7 of R(47)
        "a1": (-513.83, 0.03), "a2": (-237.40, 0.03), "a3": (-228.08, 0.03), "a4": (-218.78, 0.03),
        "a5": (-209.96, 0.03), "a6": (-97.74, 0.03), "a8": (-73.92, 0.03), "a9": (-59.30, 0.03)}),
    "640 127I2 P(10) 8-5": dict(iso="127I2", v=(8, 5), J=10, branch="P", s=90, data={  # Table 1
        "a1": (-495.4, 0.4), "a2": (-241.5, 0.7), "a3": (-233.0, 0.4), "a4": (-177.8, 1.3), "a5": (-175.2, 0.6),
        "a6": (-130.8, 0.1), "a7": (-82.45, 0.03), "a8": (-61.85, 0.14), "a9": (0.0, 0.0), "a10": (77.84, 0.03),
        "a11": (186.22, 0.07), "a12": (199.51, 0.07), "a13": (256.6, 0.2), "a14": (272.75, 0.07),
        "a15": (374.0, 0.2)}),
    "640 127I2 R(16) 8-5": dict(iso="127I2", v=(8, 5), J=16, branch="R", s=30, data={  # Table 2, bN; zero a9 of P(10)
        "a1": (62.834, 0.01), "a2": (329.8, 0.2), "a3": (335.99, 0.02)}),
    "633 129I2 P(54) 8-4": dict(iso="129I2", v=(8, 4), J=54, branch="P", s=45, data={  # Table 3
        "a2": (-449.0, 2.0), "a3": (-443.0, 2.0), "a4": (-434.0, 2.0), "a5": (-429.0, 2.0), "a6": (-360.9, 1.0),
        "a7": (-345.1, 1.0), "a8": (-340.8, 1.0), "a9": (-325.4, 1.0), "a10": (-307.0, 1.0), "a11": (-298.2, 1.0),
        "a12": (-293.1, 1.0), "a13": (-289.7, 1.0), "a14": (-282.7, 1.0), "a15": (-206.1, 0.2),
        "a16": (-197.73, 0.08), "a17": (-193.23, 0.08), "a18": (-182.74, 0.03), "a19": (-162.61, 0.05),
        "a20": (-155.72, 0.05), "a21": (-138.66, 0.05), "a22": (-130.46, 0.05), "a23": (-98.22, 0.03),
        "a24": (-55.6, 0.5), "a25": (-55.6, 0.5), "a26": (-43.08, 0.03), "a27": (-41.24, 0.05), "a28": (0.0, 0.0)}),
    "633 129I2 P(69) 12-6": dict(iso="129I2", v=(12, 6), J=69, branch="P", s=170, data={  # Table 4, bN; zero a28 of P(54)
        "a1": (99.12, 0.05), "a2": (116.08, 0.05), "a3": (132.05, 0.05), "a4": (234.54, 0.05), "a5": (256.90, 0.05),
        "a6": (264.84, 0.05), "a7": (288.06, 0.05), "a8": (337.75, 0.1), "a9": (358.8, 0.5), "a10": (358.8, 0.5),
        "a11": (373.80, 0.05), "a12": (387.24, 0.05), "a13": (395.3, 0.2), "a14": (402.45, 0.05), "a15": (407.0, 4.0),
        "a16": (412.37, 0.05), "a17": (417.0, 4.0), "a21": (507.66, 0.10), "a22": (532.65, 0.10),
        "a23": (536.59, 0.10), "a24": (545.06, 0.05), "a25": (560.94, 0.05), "a26": (566.19, 0.05),
        "a27": (586.27, 0.03), "a28": (601.78, 0.03), "a29": (620.85, 0.03), "a30": (632.42, 0.03),
        "a31": (644.09, 0.03), "a32": (655.47, 0.03), "a33": (666.81, 0.10), "a34": (692.45, 0.10),
        "a35": (697.96, 0.10), "a36": (705.43, 0.10)}),
    "633 129I2 R(60) 8-4": dict(iso="129I2", v=(8, 4), J=60, branch="R", s=250, data={  # Table 5, dN; zero a28 of P(54)
        "a23": (-555.0, 5.0), "a24": (-511.0, 2.0), "a25": (-511.0, 2.0), "a26": (-499.0, 2.0), "a27": (-499.0, 2.0),
        "a28": (-456.0, 2.0)}),
    "612 129I2 P(110) 10-2": dict(iso="129I2", v=(10, 2), J=110, branch="P", s=90, data={  # Table 4; zero a7 of 127I2 R(47)
        "a1": (-376.29, 0.05), "a2": (-244.76, 0.10), "a3": (-230.79, 0.20), "a4": (-229.40, 0.20),
        "a5": (-216.10, 0.05), "a6": (-149.37, 0.10), "a7": (-134.68, 0.10), "a8": (-130.98, 0.10),
        "a9": (-116.67, 0.05), "a10": (-96.26, 0.20), "a11": (-90.70, 0.20), "a12": (-84.12, 0.20),
        "a13": (-77.79, 0.20), "a14": (-72.70, 0.20), "a15": (1.61, 0.20), "a16": (10.63, 0.15), "a17": (15.82, 0.20),
        "a18": (25.32, 0.10), "a19": (49.44, 0.15), "a20": (54.66, 0.20), "a21": (69.02, 0.10), "a22": (74.47, 0.15),
        "a23": (110.60, 0.10), "a24": (153.09, 0.20), "a25": (154.70, 0.20), "a26": (163.98, 0.20),
        "a27": (166.22, 0.20), "a28": (208.29, 0.10)}),
    "612 129I2 R(113) 14-4": dict(iso="129I2", v=(14, 4), J=113, branch="R", s=250, data={  # Table 5, bN; zero a7 of 127I2 R(47)
        "a19": (-410.4, 0.3), "a20": (-390.0, 0.3), "a21": (-383.9, 0.5), "a22": (-362.8, 0.3), "a23": (-352.9, 0.3),
        "a24": (-346.4, 0.3), "a25": (-330.0, 0.3), "a26": (-324.9, 0.3), "a27": (-304.7, 0.3), "a28": (-289.4, 0.5),
        "a29": (-273.1, 0.3), "a30": (-255.7, 0.5), "a31": (-247.0, 5.0), "a32": (-237.0, 5.0), "a33": (-223.0, 5.0),
        "a34": (-198.6, 0.3), "a35": (-193.1, 0.3), "a36": (-187.0, 0.3)}),
    "633 129I2 P(33) 6-3": dict(iso="129I2", v=(6, 3), J=33, branch="P", s=40, data={  # Table 6, eN; zero e2
        "a1": (-19.82, 0.05), "a2": (0.0, 0.0), "a3": (17.83, 0.03), "a4": (102.58, 0.05), "a5": (141.0, 2.0),
        "a6": (157.0, 2.0), "a7": (191.0, 2.0), "a8": (208.0, 2.0), "a9": (239.0, 2.0), "a10": (249.0, 2.0),
        "a11": (260.0, 2.0), "a12": (269.0, 3.0), "a13": (273.0, 4.0), "a14": (287.0, 4.0), "a15": (293.0, 5.0),
        "a16": (295.0, 5.0), "a17": (306.0, 6.0)}),
    "633 127I129I P(33) 6-3": dict(iso="127I129I", v=(6, 3), J=33, branch="P", s=110, data={  # Table 7, mN; zero a28 of 129I2 P(54)
        "a1": (-254.0, 3.0), "a2": (-233.71, 0.10), "a3": (-226.14, 0.10), "a4": (-207.0, 2.0), "a5": (-117.79, 0.10),
        "a6": (-87.83, 0.15), "a7": (-78.2, 0.5), "a8": (-56.0, 1.0), "a9": (-17.55, 0.05), "a10": (12.04, 0.03),
        "a11": (15.60, 0.03), "a12": (33.16, 0.03), "a13": (39.9, 0.2), "a14": (41.3, 0.2), "a15": (50.72, 0.03),
        "a16": (54.06, 0.10), "a17": (69.33, 0.03), "a18": (75.06, 0.03), "a19": (80.00, 0.03), "a20": (95.00, 0.03),
        "a21": (160.74, 0.03), "a22": (199.52, 0.03), "a23": (205.06, 0.05), "a24": (207.9, 0.5), "a25": (207.9, 0.5),
        "a26": (212.80, 0.05), "a27": (219.43, 0.05), "a28": (256.90, 0.10), "a29": (264.84, 0.05),
        "a30": (299.22, 0.05), "a31": (312.43, 0.05), "a32": (324.52, 0.03), "a33": (333.14, 0.03),
        "a34": (337.7, 0.5), "a35": (337.7, 0.5), "a36": (345.05, 0.05), "a37": (362.18, 0.10), "a38": (369.78, 0.03),
        "a39": (380.37, 0.03), "a40": (385.0, 4.0), "a41": (431.0, 4.0), "a42": (445.0, 4.0), "a43": (456.7, 0.5),
        "a44": (477.17, 0.05), "a45": (486.43, 0.05), "a46": (495.16, 0.05), "a47": (503.55, 0.05),
        "a48": (515.11, 0.05)}),
}
# 515 nm tables are in kHz; the values above are converted to MHz.

#: Features the tables assign to both a 129I2 and a 127I129I component (the same letter in the mixed-cell
#: spectrum, or an explicit "see" note). Their positions are not independent measurements, so they are skipped.
BLENDS = {"633 129I2 P(54) 8-4": {"a24", "a25"},  # n2, n1 = 127I129I m8 (n)
          "633 129I2 P(69) 12-6": {"a5", "a6", "a8", "a12"},  # r'', q'', k'', d''
          "633 127I129I P(33) 6-3": {"a8", "a28", "a29", "a34", "a35", "a40"}}
#: Components outside 3 sigma, checked instead against the 1 MHz that S06 quotes for isotopologue predictions.
OUTLIERS = {("633 127I129I P(33) 6-3", "a47"): "m47 (r') is 0.42 MHz (8 u_c) below the model; restoring the "
            "heteronuclear spin-spin terms that S06 omits does not move it"}
XFAIL = {"515 127I2 R(98) 58-1": "v' = 58 is beyond the validated range: the S06 formulae stop at v' = 53 (frozen "
         "above), and the Hannover potentials put this line 41 GHz from its measured position"}
UC_GOOD = (0.1, 0.5, np.inf)  # MHz: the components that fix the common offset (tightest cut leaving >= 3)


@pytest.fixture(scope="module")
def models():
    return {iso: RovibronicModel(iso) for iso in ("127I2", "129I2", "127I129I")}


def components(models, name):
    L = LINES[name]
    nu0, comps = models[L["iso"]].hyperfine_components(*L["v"], L["J"], L["branch"])  # dJ = 2
    return nu0, {c.label: c.offset for c in comps if c.label}


def residuals(models, name):
    """Model - table (kHz) per unblended component, less their mean over the well-measured ones."""
    L = LINES[name]
    _, o = components(models, name)
    keys = [k for k in L["data"] if k not in BLENDS.get(name, ())]
    r = np.array([o[k] - L["data"][k][0] for k in keys]) * 1e3
    uc = np.array([L["data"][k][1] for k in keys]) * 1e3
    for cut in UC_GOOD:
        good = uc <= cut * 1e3
        if good.sum() >= 3:
            break
    return keys, r - r[good].mean(), uc


@pytest.mark.parametrize("name", [pytest.param(n, marks=pytest.mark.xfail(reason=XFAIL[n], strict=True))
                                  if n in XFAIL else n for n in LINES])
def test_component_positions(models, name):
    """Every component within 3 sqrt(s^2 + u_c^2) of the table, s being the line's achieved scatter."""
    keys, r, uc = residuals(models, name)
    s = LINES[name]["s"]
    misses = {k: round(x) for k, x, u in zip(keys, r, uc) if (name, k) not in OUTLIERS and abs(x) > 3 * np.hypot(s, u)}
    assert not misses
    for line, k in OUTLIERS:
        if line == name:
            assert abs(r[keys.index(k)]) < 1000.0


#: CIPM recommended frequencies (MHz) of the reference components.
ABSOLUTE = {("633 127I2 R(127) 11-5", "a16"): 473_612_353.604, ("515 127I2 P(13) 43-0", "a3"): 582_490_603.442,
            ("543 127I2 R(106) 28-0", "a10"): 551_580_162.400, ("576 127I2 P(62) 17-1", "a1"): 520_206_808.4,
            ("612 127I2 R(47) 9-2", "a7"): 489_880_354.93, ("640 127I2 P(10) 8-5", "a9"): 468_218_332.4}
#: [f(line 1, comp 1) - f(line 2, comp 2)] / MHz within one isotopologue, from the table notes and zeros.
INTERVALS = [("633 127I2 P(33) 6-3", "a21", "633 127I2 R(127) 11-5", "a16", -532.42),
             ("515 127I2 R(15) 43-0", "a1", "515 127I2 P(13) 43-0", "a1", 283.835),
             ("543 127I2 R(12) 26-0", "a4", "543 127I2 R(106) 28-0", "a10", -853.336),
             ("612 127I2 P(48) 11-3", "a15", "612 127I2 R(47) 9-2", "a7", -160.318),
             ("612 127I2 R(48) 15-5", "a1", "612 127I2 R(47) 9-2", "a7", -513.83),
             ("640 127I2 R(16) 8-5", "a1", "640 127I2 P(10) 8-5", "a9", 62.834),
             ("633 129I2 P(69) 12-6", "a27", "633 129I2 P(54) 8-4", "a28", 586.27),
             ("633 129I2 R(60) 8-4", "a28", "633 129I2 P(54) 8-4", "a28", -456.0)]
#: The same between isotopologues.
ISOTOPE_SHIFTS = [("633 129I2 P(54) 8-4", "a28", "633 127I2 R(127) 11-5", "a16", -42.99),
                  ("633 129I2 P(33) 6-3", "a2", "633 127I2 R(127) 11-5", "a16", 849.4),
                  ("612 129I2 P(110) 10-2", "a1", "612 127I2 R(47) 9-2", "a7", -376.29),
                  ("612 129I2 R(113) 14-4", "a36", "612 127I2 R(47) 9-2", "a7", -187.0),
                  ("633 127I129I P(33) 6-3", "a10", "633 129I2 P(54) 8-4", "a28", 12.04)]


def frequency(models, name, label):
    nu0, o = components(models, name)
    return nu0 + o[label]


def test_absolute_frequencies(models):
    """Achieved -1.9 to +1.4 MHz (rms 1.3 MHz); the rovibronic model is good to about 1.5 MHz."""
    diffs = [frequency(models, *key) - f for key, f in ABSOLUTE.items()]
    assert np.abs(diffs).max() < 4.0


def test_intervals_between_lines(models):
    """Achieved -4.2 to +0.5 MHz (rms 1.6 MHz)."""
    diffs = [frequency(models, a, ka) - frequency(models, b, kb) - f for a, ka, b, kb, f in INTERVALS]
    assert np.abs(diffs).max() < 5.0


def test_isotope_shifts(models):
    """With hannover2008 as printed, 129I2 lines came out 7-14 MHz low against 127I2 and 127I129I about
    5 MHz high against 129I2 (Salumbides 2006 Fig. 2 shows the same). i2spec2026a adds 0.02154 cm-1 to
    the constant term of the B-state V_ad, which fixes that to 2.4 MHz rms; R(113) 14-4 stays at +5.1."""
    diffs = [frequency(models, a, ka) - frequency(models, b, kb) - f for a, ka, b, kb, f in ISOTOPE_SHIFTS]
    assert np.abs(diffs).max() < 6.0
    assert np.sqrt(np.mean(np.square(diffs))) < 3.5


def test_hannover2008_as_printed_misses_the_isotope_shifts():
    from i2spec.model import RovibronicModel
    ms = {iso: RovibronicModel(iso, "hannover2008") for iso in ("127I2", "129I2", "127I129I")}
    diffs = [frequency(ms, a, ka) - frequency(ms, b, kb) - f for a, ka, b, kb, f in ISOTOPE_SHIFTS]
    assert np.abs(diffs).max() > 7.0
