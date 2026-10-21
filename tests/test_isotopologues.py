"""Release 20261021, the i2spec model: balanced-expanded unchanged for 127I2, carried to 129I2 and 127I129I (isotope
shift fitted to the CIPM isotopologue lines, hyperfine scaled per nucleus) and to J'' = 445."""
import numpy as np
import pytest

from i2spec.research import ResearchModel, _entry, _read


@pytest.fixture(scope="module")
def models():
    return {iso: ResearchModel("i2spec", iso) for iso in ("127I2", "129I2", "127I129I")}


def test_the_release_extends_the_frozen_model_unchanged(models):
    entry, base = _entry("i2spec"), _entry("balanced_expanded")
    assert entry["extends"]["id"] == base["id"]
    assert all(base["files_sha256"][k] == v for k, v in entry["extends"]["files_sha256"].items())
    old, new = ResearchModel("balanced_expanded"), models["127I2"]
    for line in ("R(56) 32-0", "R(98) 58-1", "R(20) 70-0", "P(13) 43-0"):
        assert new.frequency(line) == old.frequency(line)
        assert new.frequency(line, "a1") == old.frequency(line, "a1")
        assert new.uncertainty(line).sigma_MHz == old.uncertainty(line).sigma_MHz


@pytest.mark.parametrize("iso", ["127I2", "129I2", "127I129I"])
def test_frozen_references_replay(models, iso):
    assert models[iso].verify()["passed"]


def test_the_isotope_shift_reproduces_the_cipm_lines():
    ext = _read(_entry("i2spec"), "isotopologues.json")
    shift = ext["isotope_shift"]
    assert shift["power"] == 5 and len(shift["coefficients"]) == 2
    # chosen by the training objective: the lowest leave-one-line-out rms of the forms tried
    assert shift["selection"][str(shift["power"])]["loo_rms_MHz"] == min(f["loo_rms_MHz"] for f in shift["selection"].values())
    assert shift["floor_MHz"] < 0.1
    # the 633 nm isotope shift: 129I2 P(54) 8-4 a28 is 42.99(4) MHz below 127I2 R(127) 11-5 a16 (bipm2003a)
    i127, i129 = ResearchModel("i2spec"), ResearchModel("i2spec", "129I2")
    d = i129.frequency("P(54) 8-4", "a28")-i127.frequency("R(127) 11-5", "a16")
    sigma = np.hypot(i129.uncertainty("P(54) 8-4", "a28").sigma_MHz, i127.uncertainty("R(127) 11-5", "a16").sigma_MHz)
    assert abs(d-(-42.99)) < 3*sigma


def test_isotopologue_hyperfine_and_uncertainty(models):
    m129, m127129 = models["129I2"], models["127I129I"]
    # 129I2: I = 7/2 nuclei, 28 or 36 main components; 127I129I: distinguishable nuclei, 48
    assert len([c for c in m129.hyperfine_components(8, 4, 54, "P")[1] if c.label]) == 28
    assert len([c for c in m129.hyperfine_components(6, 3, 33, "P")[1] if c.label]) == 36
    assert len([c for c in m127129.hyperfine_components(6, 3, 33, "P")[1] if c.label]) == 48
    # near the CIPM lines the uncertainty is small and interpolated; at 532 nm the shift is extrapolated
    near, far = m129.uncertainty("P(54) 8-4"), m129.uncertainty("R(56) 32-0")
    assert near.sigma_MHz < 0.5 and near.status == "interpolated"
    assert 1 < far.sigma_MHz < 20 and far.status == "extrapolated"
    # the conditioned covariance: two 129I2 lines that share their levels' shift cancel it in an interval
    both = m129.uncertainty("P(54) 8-4", ref_line="R(60) 8-4")
    assert both.sigma_MHz < np.hypot(near.sigma_MHz, m129.uncertainty("R(60) 8-4").sigma_MHz)


def test_domains_and_high_J(models):
    assert models["129I2"].info.domain.isotopologue == "129I2"
    assert models["129I2"].info.domain.B_v_max == 81 and models["127I2"].info.domain.J_lower_max == 445
    u = models["127I2"].uncertainty("R(260) 10-2")
    assert u.status == "extrapolated" and np.isfinite(models["127I2"].frequency("R(260) 10-2", "a1"))
    with pytest.raises(ValueError, match="dissociated"):
        models["127I2"].frequency("R(401) 22-6")
