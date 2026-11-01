"""The fit's packaging (i2spec.fit.package, release) rebuilds the release's tables and references."""
import json

import numpy as np
import pytest

from i2spec.fit import b1g_resonance as rz
from i2spec.research import _entry, _read


def test_resonance_overlaps_carry_one_phase_convention(table):
    """The B~1g coupling mixes B J and J +- 2, so the overlap of the 1g level with B (v, J) must keep one phase across
    J: its sign at fixed 1g rotation is the same for every B J. The magnitudes are the release's."""
    m = table["B"]["model"]
    res = rz.Resonance(m["curve_constants"], m["background"], m["theta_deg"])
    v = 78
    n = res.resonant_index(v, rz.PIQUE[v]["dE"])
    released = table["resonance"]["levels"][str(v)]["overlaps"]
    for Jg in (1, 2, 3, 6):
        s = [res.overlap(n, Jg, v, jb) for jb in range(max(0, Jg-4), Jg+5)]
        assert len(set(np.sign(s))) == 1
        assert np.allclose(np.abs(s), [abs(released[f"{Jg},{jb}"]) for jb in range(max(0, Jg-4), Jg+5)], rtol=1e-8)


def test_a_level_without_neighbours_at_its_J_is_not_matched():
    from i2spec.fit import package
    E = {("B", 50, 10): 100., ("B", 51, 10): 110.}
    assert package._match(E, "B", 50, 10, [99.9, 109.8]) == 0
    assert package._match({("B", 50, 10): 100.}, "B", 50, 10, [99.9]) is None


@pytest.mark.slow
def test_rebuild_reproduces_the_release():
    from i2spec.fit import release
    built = release.rebuild()
    files = release.records(built)
    entry = _entry("balanced_expanded")
    model, ref = files["model.json"], _read(entry, "model.json")
    assert json.loads(json.dumps(model["structure"])) == ref["structure"]
    for s in ref["values"]:
        for n, z in ref["values"][s].items():
            assert model["values"][s][n] == pytest.approx(z, rel=1e-6, abs=1e-6), (s, n)
    assert model["inventory"] == ref["inventory"]
    hf, href = files["hyperfine.json"], _read(entry, "hyperfine.json")
    assert len(hf["B"]["rows"]) == len(href["B"]["rows"])
    # near a perturber level a second-order term reaches 1e5 kHz, and a fresh solve rounds it 5e-8 apart
    assert np.max(abs(np.array(hf["B"]["rows"])-np.array(href["B"]["rows"]))) < 1e-2
    assert np.max(abs(np.array(hf["X"]["rows"])-np.array(href["X"]["rows"]))[:, 2]) < 1e-3     # X eqQ, MHz
    refs, rref = files["reference_predictions.json"], _read(entry, "reference_predictions.json")
    assert [(r["line"], r["component"]) for r in refs] == [(r["line"], r["component"]) for r in rref]
    far = [abs(a["predicted_MHz"]-b["predicted_MHz"]) for a, b in zip(refs, rref) if a["line"][2] not in (77, 78)]
    assert max(far) < 1e-3                  # MHz; B v' 77, 78 now carry consistent wavefunction phases
    unc, uref = files["uncertainty.json"], _read(entry, "uncertainty.json")
    assert unc["statistical"]["parameters"] == uref["statistical"]["parameters"]
    for key in ("all", "precise", "margin_J", "neighbor_v"):
        assert json.loads(json.dumps(unc["coverage"][key])) == uref["coverage"][key]
