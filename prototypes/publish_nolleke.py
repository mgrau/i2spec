"""Publish the Nölleke et al. (2018) 915-985 nm line list and scans (CC BY 4.0) in the repository.

The source (data/external/nolleke_2018, from Mendeley Data doi:10.17632/cvr6vv2zps.1) is redistributable with
attribution and a note of the changes. Changes made here: the columns are renamed, the vacuum wavenumber is added
(1e7 / vacuum wavelength), and the scans are gzipped with their decimal separator and delimiter unchanged. The
lines carry no assignments (docs/research/atlas-line-positions.md says why).

usage: publish_nolleke.py   -> data/atlas_lines/nolleke_2018.{csv,toml}, data/spectra/nolleke_2018_scans.csv.gz
"""
from __future__ import annotations

import csv
import gzip
import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "data" / "external" / "nolleke_2018"
ATLAS = ROOT / "data" / "atlas_lines"
SPECTRA = ROOT / "data" / "spectra"

CITATION = ("C. Nölleke, C. Raab, R. Neuhaus, S. Falke, Absolute frequency atlas from 915 nm to 985 nm based on laser "
            "absorption spectroscopy of iodine, J. Mol. Spectrosc. 346, 19-22 (2018), doi:10.1016/j.jms.2017.12.013. "
            "Data: C. Nölleke et al., Laser absorption spectroscopy of iodine between 915 and 985 nm, Mendeley Data, "
            "V1 (2017), doi:10.17632/cvr6vv2zps.1")


def main():
    rows = list(csv.DictReader(open(SRC / "iodine_atlas.csv", newline="", encoding="utf-8")))
    with open(ATLAS / "nolleke_2018.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["nu_cm1", "lambda_vac_nm", "absorption", "water_line_lambda_vac_nm", "water_line_absorption",
                    "reference_lambda_vac_nm"])
        for r in rows:
            lam = float(r["Vacuum-wavelength identified line (nm)"])
            w.writerow([f"{1e7 / lam:.4f}", r["Vacuum-wavelength identified line (nm)"], r["Absorption"],
                        r["Wavelength near water line (nm)"], r["Absorption water line (a.u.)"],
                        r["Vacuum-wavelength in [1] (nm)"]])
    nus = [1e7 / float(r["Vacuum-wavelength identified line (nm)"]) for r in rows]
    n_water = sum(1 for r in rows if r["Wavelength near water line (nm)"])
    n_ref = sum(1 for r in rows if r["Vacuum-wavelength in [1] (nm)"])
    SPECTRA.mkdir(exist_ok=True)
    with open(SRC / "scan_data.csv", "rb") as fin, gzip.open(SPECTRA / "nolleke_2018_scans.csv.gz", "wb", compresslevel=9) as fout:
        shutil.copyfileobj(fin, fout)
    (ATLAS / "nolleke_2018.toml").write_text(f'''id = "nolleke_2018"
unit = "cm-1"
isotopologue = "127I2"
kind = "unassigned atlas line positions and absorptions"
citation = """{CITATION}"""
doi = "10.1016/j.jms.2017.12.013"
url = "https://doi.org/10.17632/cvr6vv2zps.1"
licence = "CC BY 4.0 (https://creativecommons.org/licenses/by/4.0/): cite the paper and the data set"
instrument = """Doppler-limited absorption of a mode-hop-free tunable diode laser, wavemeter referenced; 225 cm path, cell body at
300 °C, cold finger 39 °C (127 Pa)"""
calibration = """Absolute, stated to 50 MHz by the authors"""
extraction = """The authors' line list (iodine_atlas.csv of the data set) with the columns renamed and nu_cm1 = 1e7 / lambda_vac_nm
added (prototypes/publish_nolleke.py). Nothing else is changed. The {n_water} lines with a water_line entry lie near an
atmospheric H2O line in the open beam and may be blended with it; reference_lambda_vac_nm is the authors' comparison line
from their reference [1] where they give one ({n_ref} lines). The raw scans, transmission against vacuum wavelength, are in
data/spectra/nolleke_2018_scans.csv.gz (the authors' scan_data.csv, gzipped). The lines are not assigned: two attempts
against the model failed (docs/research/atlas-line-positions.md), because they come from X v'' = 24-45, the levels the
assignment would measure."""
range_cm-1 = [{min(nus):.1f}, {max(nus):.1f}]
n_lines = {len(rows)}
''')
    print(f"{len(rows)} lines, {min(nus):.1f}-{max(nus):.1f} cm-1 -> {ATLAS / 'nolleke_2018.csv'}; scans -> "
          f"{SPECTRA / 'nolleke_2018_scans.csv.gz'} ({(SPECTRA / 'nolleke_2018_scans.csv.gz').stat().st_size / 1e6:.1f} MB)")


if __name__ == "__main__":
    main()
