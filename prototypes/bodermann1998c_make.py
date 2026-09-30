"""Write data/observations/bodermann1998c/data.csv from Bodermann's thesis, Tables 6.12 and 7.1 (verified on page images).

Run from the repository root: python prototypes/bodermann1998c_make.py. The transcription and every
exclusion and relabel are documented in data/observations/bodermann1998c/meta.toml.

The table prints unsigned beat frequencies between Uebergang 1 and 2. The sign is fixed by the model:
line = the higher-frequency component, ref = the lower, so every interval value is positive. Each row
is also checked against the model magnitude - a transcription slip or misread label would show as MHz.
"""
import csv
from i2spec.observations import Line, Predictor, Observation

T71 = """R(180) 0-16|a1|P(34) 0-17|a10|15270134.3|30.9
R(180) 0-16|a1|P(34) 0-17|a15|14976019.3|30.9
R(180) 0-16|a1|P(239) 0-15|a2|7759664.6|45.0
R(180) 0-16|a1|P(239) 0-15|a13|7172466.0|40.0
P(34) 0-17|a10|R(44) 0-17|a10|13233567.7|34.3
P(34) 0-17|a15|R(44) 0-17|a15|13233686.5|31.2
R(44) 0-17|a10|P(172) 0-16|a10|20997223.5|41.3
R(44) 0-17|a15|P(172) 0-16|a15|20997488.2|30.3
P(171) 0-16|a14|R(180) 0-16|a1|19594882.0|30.0
P(36) 0-16|c1|R(45) 0-16|b20|114991.7|15.6
P(36) 0-16|c10|R(45) 0-16|b20|701958.7|15.3
P(36) 0-16|c1|R(45) 0-16|b14|370890.6|16.0
P(36) 0-16|c1|P(171) 0-15|a20|420458.0|19.7
P(36) 0-16|c1|P(171) 0-15|a11|847510.4|8.3
R(45) 0-16|b20|P(171) 0-15|a11|732528.0|14.6
R(45) 0-16|b14|P(171) 0-15|a2|930989.3|22.4
P(171) 0-15|a2|R(46) 0-16|a10|15217083.8|26.4
P(171) 0-15|a2|R(46) 0-16|a15|14922948.9|26.8
R(96) 0-15|a10|P(88) 0-15|a15|26006452.0|15.9
R(96) 0-15|a1|P(88) 0-15|a15|25418492.5|11.7
P(54) 0-15|a10|R(63) 0-15|a13|5766740.3|17.6
P(54) 0-15|a1|R(63) 0-15|a13|5178962.7|24.3
P(54) 0-15|a1|R(63) 0-15|a2|5766820.5|19.0
R(62) 0-15|a10|R(63) 0-15|a13|23206957.5|25.9
R(62) 0-15|a1|R(63) 0-15|a13|22619114.7|25.5
R(63) 0-15|a13|R(64) 0-15|a10|23604195.8|25.1
R(63) 0-15|a13|R(64) 0-15|a1|24192057.7|32.9
P(55) 0-15|a13|R(64) 0-15|a10|5884611.0|38.7
P(55) 0-15|a2|R(64) 0-15|a10|5296814.8|41.6
R(64) 0-15|a10|P(56) 0-15|a10|17997970.6|20.3
R(64) 0-15|a10|P(56) 0-15|a1|18585789.3|14.5
R(64) 0-15|a10|R(65) 0-15|a13|24001486.5|15.6
R(64) 0-15|a10|R(65) 0-15|a2|24589360.3|16.1
R(10) 0-15|a9|R(12) 0-15|a9|5630438.4|30.0
R(12) 0-15|a9|R(17) 0-15|a8|21069912.8|24.8
R(17) 0-15|a8|P(12) 0-15|a9|22288102.7|30.5
R(174) 0-14|a1|R(16) 0-15|a9|11537920.6|20.2
R(16) 0-15|a9|R(19) 0-15|a15|15955960.1|35.4
R(16) 0-15|a9|R(19) 0-15|a2|16622728.8|29.7
R(16) 0-15|a9|P(166) 0-14|a15|15900116.9|27.4
R(16) 0-15|a9|P(166) 0-14|a1|16782292.7|28.0
P(87) 0-14|a2|R(96) 0-14|a10|15976628.4|4.5
R(138) 1-14|a10|R(16) 0-14|a1|20628625.3|11.9
R(138) 1-14|a15|R(16) 0-14|a15|20043421.1|13.0
R(240) 0-12|a1|R(138) 1-14|a10|24257255.8|30.0
P(164) 0-13|a10|R(240) 0-12|a1|16912959.8|30.7
P(164) 0-13|a10|R(240) 0-12|a15|16030721.1|30.7
P(129) 1-14|b3|P(164) 0-13|a1|504154.6|23.5
P(129) 1-14|b18|P(164) 0-13|a1|1237286.5|19.3"""
DUPLICATES = "R(138) 1-14 a10 - R(16) 0-14 a10, R(240) 0-12 a10 - R(138) 1-14 a10, P(164) 0-13 a10 - R(240) 0-12 a10"

P = Predictor("i2spec2026d")
L = lambda s: Line.parse("127I2 " + s)   # noqa: E731
f = lambda s, c: P(Observation(L(s), c, "frequency", 0.0, 1.0))   # noqa: E731
rows = []
# Table 6.12: the three primary wavelength-comparison lines as absolute frequencies (kHz -> MHz) ...
for s, c, v, u, note in (("R(96) 0-15", "a10", 377796992463, 48, "Table 6.12, wavelength comparison"),
                         ("P(166) 0-14", "a1", 379431381762, 48, "Table 6.12, wavelength comparison"),
                         ("R(96) 0-14", "a10", 383617264515, 49, "Table 6.12, wavelength comparison")):
    rows.append(["127I2 " + s, c, "frequency", f"{v/1e3:.3f}", f"{u/1e3:.3f}", "", "", "bodermann1998c", note])
# ... and the two it connects to them by beat measurements, as intervals to the line they were beaten against
for s, c, v, u, rs, rc, rv, ru in (("P(228) 1-14", "b10", 377796446422, 80, "R(96) 0-15", "a10", 377796992463, 48),
                                   ("R(19) 0-15", "b2", 379431541326, 62, "P(166) 0-14", "a1", 379431381762, 48)):
    rows.append(["127I2 " + s, c, "interval", f"{(v - rv)/1e3:.3f}", f"{(u*u - ru*ru)**0.5/1e3:.3f}",
                 "127I2 " + rs, rc, "bodermann1998c",
                 "Table 6.12 absolute values differenced: the thesis calibrated this line by a beat against the line of ref_line; "
                 "uncertainty sqrt(u^2 - u_ref^2)"])
worst = 0
for r in T71.splitlines():
    l1, c1, l2, c2, v, sd = r.split("|")
    d = f(l2, c2) - f(l1, c1)
    line, comp, ref, rcomp = (l2, c2, l1, c1) if d > 0 else (l1, c1, l2, c2)
    mism = abs(abs(d) - float(v) / 1e3)
    worst = max(worst, mism)
    flag = "" if mism < 20 else f"  <-- model differs by {mism:.1f} MHz"
    print(f"{l1:12s} {c1:4s} {l2:12s} {c2:4s} {float(v)/1e3:12.3f}  model |d| {abs(d):12.3f}  diff {abs(d) - float(v)/1e3:+8.3f}{flag}")
    note = "Table 7.1 beat frequency; sign from the order of the two lines"
    if (l1, c1, l2, c2) == ("P(54) 0-15", "a1", "R(63) 0-15", "a2"):
        note += "; RELABELLED: printed R(63) 0-15 a01, entered as a02 (see meta.toml)"
    rows.append(["127I2 " + line, comp, "interval", f"{float(v)/1e3:.4f}", f"{float(sd)/1e3:.4f}", "127I2 " + ref, rcomp,
                 "bodermann1998c", note])
print("worst |model - printed|:", round(worst, 3), "MHz;", len(rows), "rows")
with open("data/observations/bodermann1998c/data.csv", "w", newline="") as fh:
    w = csv.writer(fh); w.writerow("line,component,kind,value,uncertainty,ref_line,ref_component,group,note".split(","))
    w.writerows(rows)
