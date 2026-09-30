"""The local near-infrared Dunham models, against the measurements they claim to describe.

Data: ``data/observations/liao2010a`` (Liao et al., JOSA B 27, 1208 (2010): 31 comb-referenced
components, 755-810 nm) and ``data/observations/bodermann2000a`` (Bodermann et al., EPJD 11, 213
(2000): 29 absolute frequencies, 778-795 nm).

Every tolerance below is set just outside what the models actually achieve today, so that a change
in either direction fails the test. The numbers are tabulated in ``docs/research/nir-model-2010.md``.
"""

import numpy as np
import pytest

from i2spec.local_nir import load_local_nir
from i2spec.observations import DATA_DIR, Predictor, load_dataset

#: Liao et al. 2010 section 4: their own fit adds this pressure shift to the published frequencies.
PRESSURE_SHIFT_MHz = 0.114


@pytest.fixture(scope="module")
def residuals():
    """(model − observed) in MHz for every absolute frequency of both NIR data sets."""
    predictor = Predictor("i2spec2026a")     # the bare potentials: i2spec2026b carries corrections fitted to these lines
    k04 = load_local_nir("knockel2004")
    l10 = load_local_nir("liao2010", trust_printed_digits=True)
    rows = []
    for ds_id, shift in (("liao2010a", PRESSURE_SHIFT_MHz), ("bodermann2000a", 0.0)):
        for o in load_dataset(DATA_DIR / ds_id).observations:
            if o.kind != "frequency":
                continue
            line, target = o.line, o.value + shift
            centre = predictor.position(line, None)
            offset = 0.0 if o.component is None else predictor.position(line, o.component) - centre
            covered = k04.covers(line.v_upper, line.v_lower, line.J_lower)
            local = (lambda m: m.transition_MHz(line.v_upper, line.v_lower, line.J_lower, line.branch)
                     + offset - target)
            rows.append(dict(dataset=ds_id, line=str(line), band=f"{line.v_upper}-{line.v_lower}",
                             v_upper=line.v_upper, wavelength=1e7 / (centre / 29979.2458),
                             potentials=centre + offset - target,
                             knockel2004=local(k04) if covered else None,
                             liao2010=local(l10) if covered else None))
    return rows


def rms(rows, key):
    values = [r[key] for r in rows if r[key] is not None]
    assert values, "no rows selected"
    return float(np.sqrt(np.mean(np.square(values))))


def worst(rows, key):
    """Largest |residual| over the rows the model actually covers (a local model declines some)."""
    values = [abs(r[key]) for r in rows if r[key] is not None]
    assert values, "no rows selected"
    return max(values)


def best(rows, key):
    values = [abs(r[key]) for r in rows if r[key] is not None]
    assert values, "no rows selected"
    return min(values)


def test_liao2010_is_refused_unless_the_caller_insists():
    """Its printed digits cannot reproduce its own fit, so loading it needs an explicit opt-in."""
    with pytest.raises(ValueError, match="printed"):
        load_local_nir("liao2010")
    assert load_local_nir("liao2010", trust_printed_digits=True).name == "liao2010"
    with pytest.raises(ValueError, match="unknown local NIR model"):
        load_local_nir("no such model")


def test_domain():
    model = load_local_nir("knockel2004")
    assert model.covers(0, 14, 100) and model.covers(0, 12, 242) and model.covers(0, 17, 0)
    assert not model.covers(1, 14, 100)  # v′ > 0 is outside both published local models
    assert not model.covers(0, 11, 100) and not model.covers(0, 18, 100)
    assert not model.covers(0, 12, 243)  # the stated J″ <= 242
    assert model.domain.covers_wavelength(800.0) and not model.domain.covers_wavelength(760.0)
    assert load_local_nir("liao2010", trust_printed_digits=True).domain.covers_wavelength(760.0)
    with pytest.raises(ValueError, match="only v' = 0"):
        model.transition(1, 14, 50, "R")
    with pytest.raises(ValueError, match="branch"):
        model.transition(0, 14, 50, "Q")


def test_knockel2004_beats_the_potentials_inside_its_own_window(residuals):
    """0.23 MHz against 2.31 MHz on the 34 anchors between 775 and 815 nm."""
    window = [r for r in residuals if r["v_upper"] == 0 and r["wavelength"] >= 775]
    assert len(window) >= 34
    assert rms(window, "knockel2004") < 0.30
    assert rms(window, "potentials") > 2.0
    assert worst(window, "knockel2004") < 0.80
    # The one anchor the local model declines is P(243) 0-12, whose J″ = 243 is past the stated
    # J″ <= 242 — a self-consistent refusal (docs/research/nir-anchor-data-2000.md §5.4).
    assert [r["line"] for r in window if r["knockel2004"] is None] == ["127I2 P(243) 0-12"]


def test_knockel2004_fails_where_liao_measured(residuals):
    """Below 775 nm the 2004 model extrapolates and breaks down; the potentials do not."""
    head = [r for r in residuals if r["v_upper"] == 0 and r["wavelength"] < 775]
    assert len(head) >= 21
    assert rms(head, "knockel2004") > 20.0
    assert rms(head, "potentials") < 3.0


def test_the_755_nm_band_head(residuals):
    """R(48) 0-12 at 755.6 nm settles the open question of nir-model-2010.md §6.1.

    The ≈49 MHz there is the 2004 local model extrapolating, not an error in the potentials.
    """
    head = [r for r in residuals if r["line"].startswith("127I2 R(48) 0-12")]
    assert len(head) == 3
    assert worst(head, "potentials") < 0.5
    assert best(head, "knockel2004") > 40.0


def test_potentials_against_the_liao_measurements(residuals):
    """Phase A accuracy at 755-810 nm: a few MHz, not the 30-60 MHz the literature assumes."""
    liao = [r for r in residuals if r["dataset"] == "liao2010a"]
    assert len(liao) == 31
    assert rms(liao, "potentials") < 2.6
    assert worst(liao, "potentials") < 4.5


def test_liao2010_as_printed_does_not_reproduce_its_own_measurements(residuals):
    """Documented defect: Table 3 is rounded too hard (docs/research/nir-model-2010.md §3.3).

    If a future parameter source fixes this, the assertion below fails and should be revisited.
    """
    liao = [r for r in residuals if r["dataset"] == "liao2010a" and r["liao2010"] is not None]
    assert rms(liao, "liao2010") > 20.0
    calibration = [r for r in liao if r["band"] in ("0-14", "0-16")]
    assert worst(calibration, "knockel2004") < 0.1  # the 2004 model gets these right
    assert best(calibration, "liao2010") > 20.0  # its own successor does not


def test_model_metadata():
    k04, l10 = load_local_nir("knockel2004"), load_local_nir("liao2010", trust_printed_digits=True)
    assert k04.uncertainty_MHz == 0.2 and l10.uncertainty_MHz == 0.2
    assert k04.domain.wavelength_nm == (775.0, 815.0) and l10.domain.wavelength_nm == (755.0, 815.0)
    assert len(k04.upper) == len(l10.upper) == 5
    assert len(k04.lower) == len(l10.lower) == 14  # both papers use 19 parameters in total
    assert (30, 32) not in k04.lower and (3, 2) in l10.lower  # the 2010 refit moved S22 to S32
    assert l10.reproduction_MHz > 20.0 and k04.reproduction_MHz < 0.3
