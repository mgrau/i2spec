"""balanced-expanded (release 20261020): selection, frozen base, domain, port and replay contracts."""
import numpy as np
import pytest

from i2spec import ResearchModel, research_models
from i2spec.research import _entry, _read
from i2spec.research_expanded import _ExpandedPotential


@pytest.fixture(scope="module")
def model():
    return ResearchModel("balanced_expanded")


@pytest.fixture(scope="module")
def references():
    return _read(_entry("balanced_expanded"), "reference_predictions.json")


def test_selection_ledger_and_domain():
    info = next(m for m in research_models() if m.role == "balanced_expanded")
    assert info.id == "research-20261020-balanced-expanded"
    assert (info.domain.X_v_max, info.domain.B_v_max, info.domain.J_lower_max) == (108, 80, 200)
    assert info.inventory.line_specific == info.inventory.level_specific == 0
    assert info.inventory.fitted_calibration == 5
    assert ResearchModel("research-20261020-balanced-expanded").info == info


def test_base_coefficients_are_the_frozen_balanced_release(model):
    spec = model.coefficients
    assert spec["base"]["model_sha256"] == _entry("balanced")["files_sha256"]["model.json"]
    assert spec["base_coefficients"] == ResearchModel("balanced").coefficients
    assert spec["values"]["B"]["dDe"] == spec["values"]["X"]["dDe"]       # one shared dissociation-energy change


def test_tails_leave_the_base_potentials_unchanged_inside_the_knots(model):
    from i2spec.research import _Potential
    spec = model.coefficients
    base = spec["base_coefficients"]
    for state, (lo, hi) in (("X", (2.43, 3.07)), ("B", (2.66, 4.44))):
        R = np.linspace(lo, hi, 400)
        expanded = _ExpandedPotential(state, base, spec["structure"], spec["values"])
        frozen = _Potential(state, base)
        np.testing.assert_allclose(expanded(R), frozen(R), atol=1e-9, rtol=0)
        np.testing.assert_allclose(expanded.nonadiabatic(R), frozen.nonadiabatic(R), atol=1e-15, rtol=0)


@pytest.mark.parametrize("line", ["R(20) 81-0", "R(20) 32-109", "R(201) 32-0", "P(0) 50-0", ("Q", 3, 50, 0)])
def test_outside_domain_is_rejected(model, line):
    with pytest.raises(ValueError):
        model.frequency(line)


def test_dissociated_levels_and_missing_hyperfine_do_not_fall_back(model):
    with pytest.raises(ValueError, match="dissociated"):
        model.frequency("R(60) 80-0")            # B v'80 holds no J'=61 level
    center = model.frequency("R(190) 44-0")      # quasi-bound B v'44 J'=191: a center, no tabulated hyperfine
    assert np.isfinite(center)
    with pytest.raises(ValueError, match="hyperfine"):
        model.frequency("R(190) 44-0", "a1")


def test_core_lines_keep_balanced_hyperfine(model):
    balanced = ResearchModel("balanced")
    for label in ("R(56) 32-0", "P(53) 32-0"):
        _, mine = model.hyperfine_components(*_tuple(label))
        _, theirs = balanced.hyperfine_components(*_tuple(label))
        offsets = {c.label: c.offset for c in mine if c.label}
        for c in theirs:
            if c.label:
                assert abs(offsets[c.label]-c.offset) < 1e-6                   # MHz
        assert abs(model.frequency(label)-balanced.frequency(label)) < 1.     # MHz: the tails barely reach v'32


def _tuple(label):
    branch, rest = label[0], label[2:]
    J, bands = rest.split(") ")
    vu, vl = bands.split("-")
    return int(vu), int(vl), int(J), branch


@pytest.mark.parametrize("kind", ["core", "B beyond", "X beyond", "resonance"])
def test_replays_one_reference_line_of_each_kind(model, references, kind):
    def category(line):
        b, J, vu, vl = line
        return ("resonance" if vu in (77, 78) else "core" if vu <= 43 and vl <= 17 else
                "B beyond" if vl <= 17 else "X beyond" if vu <= 43 else "both")
    line = next(r["line"] for r in references if category(r["line"]) == kind and r.get("ref_component"))
    for row in (r for r in references if r["line"] == line):
        value = model.frequency(row["line"], row["component"])
        if row.get("ref_component"):
            value -= model.frequency(row["line"], row["ref_component"])
        assert abs(value-row["predicted_MHz"]) < 100e-6


@pytest.fixture(scope="module")
def uncertainty_references():
    return _read(_entry("balanced_expanded"), "uncertainty_references.json")


def _replay(model, row):
    if row.get("ref_line"):
        return model.uncertainty(row["line"], row["component"], row["ref_component"], row["ref_line"])
    return model.uncertainty(row["line"], row["component"], row["ref_component"])


@pytest.mark.parametrize("pick", ["R(56) 32-0", "R(41) 51-0", "R(10) 77-0", "R(56) 32-53", "R(3) 1-108", "interline"])
def test_uncertainties_replay_their_frozen_references(model, uncertainty_references, pick):
    rows = [r for r in uncertainty_references if "error" not in r and
            (r.get("ref_line") if pick == "interline" else r["line"] == pick and not r.get("ref_line"))]
    assert rows
    for row in rows:
        u = _replay(model, row)
        # the statistical part rests on finite-difference level derivatives, which turn last-bit differences
        # between platforms' eigensolves into about 1e-6 of sigma (Linux against the macOS references)
        for key in ("sigma_MHz", "statistical_MHz", "model_MHz", "hyperfine_MHz"):
            assert u.as_dict()[key] == pytest.approx(row[key], rel=1e-4, abs=1e-9)
        assert (u.status, u.line_class, u.hyperfine_status) == (row["status"], row["line_class"], row["hyperfine_status"])


def test_uncertainty_parts_add_in_quadrature_and_intervals_drop_the_center(model):
    center = model.uncertainty("R(56) 32-0")
    component = model.uncertainty("R(56) 32-0", "a10")
    interval = model.uncertainty("R(56) 32-0", "a15", "a1")
    assert center.hyperfine_MHz == 0. and center.hyperfine_status is None
    assert center.sigma_MHz == pytest.approx(np.hypot(center.statistical_MHz, center.model_MHz))
    assert component.statistical_MHz == center.statistical_MHz and component.hyperfine_MHz > 0
    assert interval.statistical_MHz == interval.model_MHz == 0. and interval.sigma_MHz == interval.hyperfine_MHz
    assert interval.sigma_MHz < .1 < center.sigma_MHz < 1.        # core: tens of kHz within a line, ~0.4 MHz absolute
    far = model.uncertainty("R(3) 1-108")
    assert far.status == "extrapolated" and far.sigma_MHz > 100.
    assert model.uncertainty("P(13) 43-0").sigma_MHz > 0        # the excluded line still has a prediction and a sigma
    with pytest.raises(ValueError, match="unavailable"):
        model.uncertainty("R(56) 32-0", "a40")
    with pytest.raises(ValueError):
        model.uncertainty("R(20) 81-0")


def test_uncertainty_model_is_frozen_without_pending_outcomes(model):
    table = _read(_entry("balanced_expanded"), "uncertainty.json")
    assert table["pending_outcomes_used"] is False and table["model"] == model.info.id
    assert set(table) >= {"recipe", "statistical", "model_error", "hyperfine", "coverage", "checks", "inputs_sha256"}
    entry = dict(_entry("balanced_expanded"))
    entry["files_sha256"] = dict(entry["files_sha256"], **{"uncertainty.json": "0"*64})
    with pytest.raises(ValueError, match="SHA256"):
        _read(entry, "uncertainty.json")


def test_tampered_tables_are_rejected(tmp_path, monkeypatch):
    import i2spec.research as research
    entry = dict(_entry("balanced_expanded"))
    entry["files_sha256"] = dict(entry["files_sha256"], **{"hyperfine.json": "0"*64})
    with pytest.raises(ValueError, match="SHA256"):
        _read(entry, "hyperfine.json")
    assert research.ResearchModel("balanced").info.id == "research-20261018-balanced"
