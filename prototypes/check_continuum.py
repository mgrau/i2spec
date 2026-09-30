"""Continuum and band-averaged absorption against Tellinghuisen (2011) and Spietz et al. (2006).

Usage:  uv run python prototypes/check_continuum.py

1. Continuum (A + C + B bound-free) against Tellinghuisen, JCP 135, 054301 (2011) [T11C] Table II at
   390-495 and 650-900 nm, 0 and 35 °C, with the table's wavelengths read as vacuum and as air.
2. A + C against the ε_c of Tellinghuisen, JCP 134, 084301 (2011) [T11B] at 530-635 nm (cell ~309 K).
3. Lines + continuum, averaged over the instrument function (weak-absorption limit):
   ε(436 nm), ε(500.2 nm, 35 °C) = 587.4(8) [T11C], σ(500 nm air) = 2.186(21)e-18 cm² at 0.58 nm
   [Spietz 2006], and [T11C] Table II at 500-650 nm (their pseudocontinuum; both sides smoothed
   over 4 nm, so only the band envelope is compared).
Values: docs/research/continuum-model.md.
"""

import numpy as np

from i2spec.continuum import continuum_model
from i2spec.intensity import intensity_model, master_line_list
from i2spec.spectrum import air_to_vacuum

EPS = 3.82353e-21  # cm² per (L mol⁻¹ cm⁻¹)
T0, T35 = 273.15, 308.15
TABLE2 = {  # [T11C] Table II: λ (nm) -> ε at 0 °C, 35 °C
    390: (0.18, 0.28), 395: (0.37, 0.56), 400: (0.73, 1.05), 405: (1.38, 1.89), 410: (2.47, 3.29),
    415: (4.27, 5.49), 420: (7.08, 8.83), 425: (11.34, 13.73), 430: (17.52, 20.69), 435: (26.22, 30.24),
    440: (38.07, 43.02), 445: (53.78, 59.66), 450: (74.08, 80.85), 455: (99.74, 107.28), 460: (131.49, 139.57),
    465: (170.00, 178.24), 470: (215.71, 223.55), 475: (268.73, 275.37), 480: (328.57, 332.99),
    485: (394.05, 395.08), 490: (463.20, 459.63), 495: (533.16, 523.97), 500: (600.45, 584.98),
    505: (661.21, 639.36), 510: (711.99, 684.29), 515: (749.65, 717.30), 520: (772.28, 737.03),
    525: (778.42, 742.43), 530: (768.32, 733.85), 535: (742.51, 711.76), 540: (703.36, 678.15),
    545: (653.11, 634.77), 550: (595.24, 584.47), 555: (532.61, 529.50), 560: (468.56, 472.66),
    565: (405.54, 415.96), 570: (345.89, 361.42), 575: (291.04, 310.35), 580: (242.13, 263.84),
    585: (199.56, 222.44), 590: (163.45, 186.41), 595: (133.52, 155.69), 600: (109.30, 130.01),
    605: (90.14, 108.95), 610: (75.36, 92.02), 615: (64.28, 78.68), 620: (56.23, 68.39), 625: (50.60, 60.62),
    630: (46.84, 54.87), 635: (44.51, 50.74), 640: (43.21, 47.86), 645: (42.62, 45.89), 650: (42.47, 44.59),
    660: (42.75, 43.08), 670: (42.92, 42.13), 680: (42.37, 41.00), 690: (40.91, 39.33), 700: (38.55, 37.07),
    720: (31.81, 30.98), 740: (24.04, 24.06), 760: (16.80, 17.45), 780: (10.99, 11.95), 800: (6.79, 7.79),
    850: (1.69, 2.26), 900: (0.35, 0.57),
}
T11B_EC = {530: (84.7, 85.3), 534: (76.9, 78.0), 540: (63.8, 58.6), 545: (57.4, 53.8), 555: (38.1, 41.2),
           570: (28.1, 25.8), 580: (21.2, 21.9), 590: (19.5, 23.1), 600: (21.1, 21.9),
           610: (26.3, 27.0, 23.0, 23.9), 625: (31.5, 31.1), 635: (35.1, 32.8)}


def smoothed(lines, cont, T, lam_vac, fwhm_nm, bin_cm=0.05):
    """Lines + continuum (cm²) averaged over a Gaussian of fwhm_nm centred on each vacuum wavelength."""
    lam_vac = np.atleast_1d(np.asarray(lam_vac, dtype=float))
    nu0 = 1e7 / lam_vac
    half = 3 * 1e7 * fwhm_nm / lam_vac.min() ** 2
    edges = np.arange(nu0.min() - half, nu0.max() + half, bin_cm)
    centres = (edges[1:] + edges[:-1]) / 2
    sigma = np.histogram(lines.nu, edges, weights=lines.S)[0] / bin_cm + cont.cross_section(centres, T)
    out = []
    for n, lam in zip(nu0, lam_vac):
        g = np.exp(-4 * np.log(2) * ((centres - n) / (1e7 * fwhm_nm / lam**2)) ** 2)
        out.append((g * sigma).sum() / g.sum())
    return np.array(out)


def main():
    master = master_line_list(intensity_model("127I2"), 11000.0, 20100.0, T_range=(200.0, 600.0), S_min=1e-27)
    cont = continuum_model("127I2", master.upper_cut, T_max=600.0)

    print("1. Continuum / [T11C] Table II (table λ read as vacuum | as air)")
    for band in ((390, 495), (700, 900)):
        lam = np.array([l for l in TABLE2 if band[0] <= l <= band[1]], dtype=float)
        for T, col in ((T0, 0), (T35, 1)):
            ref = np.array([TABLE2[int(l)][col] for l in lam])
            rv = cont.cross_section(1e7 / lam, T) / EPS / ref
            ra = cont.cross_section(1e7 / air_to_vacuum(lam), T) / EPS / ref
            print(f"   {band[0]}-{band[1]} nm, {T - 273.15:.0f} °C: vacuum {rv.min():.3f}-{rv.max():.3f} "
                  f"(mean {rv.mean():.3f}) | air {ra.min():.3f}-{ra.max():.3f} (mean {ra.mean():.3f})")
        if band[0] == 390:
            print("   per λ at 0 °C (vacuum | air): " + ", ".join(
                f"{l:.0f} {float(cont.cross_section(1e7 / l, T0)) / EPS / TABLE2[int(l)][0]:.3f}|"
                f"{float(cont.cross_section(1e7 / air_to_vacuum(l), T0)) / EPS / TABLE2[int(l)][0]:.3f}"
                for l in lam[::3]))

    print("\n2. A + C at 309 K / [T11B] ε_c (their spread in brackets)")
    for l, vals in T11B_EC.items():
        ac = cont.cross_section(np.array([1e7 / l]), 309.0, states=("A", "C"))[0] / EPS
        print(f"   {l} nm: model {ac:6.1f}  data {np.mean(vals):6.1f} [{min(vals)}-{max(vals)}]  ratio {ac / np.mean(vals):.3f}")

    print("\n3. Lines + continuum, weak-absorption limit")
    at = lambda T, lam, fw: smoothed(master.at(T, S_min=0.0), cont, T, lam, fw)[0]  # noqa: E731
    for T in (296.0, 298.0):
        s = at(T, float(air_to_vacuum(500.0)), 0.58)
        print(f"   σ(500.0 nm air, 0.58 nm, {T:.0f} K) = {s:.4e} cm²; Spietz 2.186(21)e-18 -> ratio {s / 2.186e-18:.3f}")
    for fw in (0.5, 1.0, 2.0):
        sv, sa = at(T35, 500.2, fw), at(T35, float(air_to_vacuum(500.2)), fw)
        print(f"   ε(500.2 nm, 35 °C, {fw} nm): vacuum {sv / EPS:6.1f}, air {sa / EPS:6.1f}; [T11C] 587.4(8)")
    e436 = cont.cross_section(np.array([1e7 / 436.0, 1e7 / air_to_vacuum(436.0)]), 298.0) / EPS
    print(f"   ε(436 nm, 298 K): vacuum {e436[0]:.2f}, air {e436[1]:.2f}; [T11C] 31.0(4)")

    print("\n   [T11C] Table II at 500-650 nm vs lines + continuum, both smoothed over 4 nm (λ as vacuum)")
    lam = np.arange(505.0, 646.0, 5.0)
    fine = np.arange(495.0, 655.0, 0.25)
    for T, col in ((T0, 0), (T35, 1)):
        ref_fine = np.interp(fine, sorted(TABLE2), [TABLE2[k][col] for k in sorted(TABLE2)])
        ref = [np.sum(np.exp(-4 * np.log(2) * ((fine - l) / 4.0) ** 2) * ref_fine)
               / np.sum(np.exp(-4 * np.log(2) * ((fine - l) / 4.0) ** 2)) for l in lam]
        model = smoothed(master.at(T, S_min=0.0), cont, T, lam, 4.0) / EPS
        r = model / np.array(ref)
        print(f"   {T - 273.15:.0f} °C: " + ", ".join(f"{l:.0f} {x:.3f}" for l, x in zip(lam, r)))


if __name__ == "__main__":
    main()
