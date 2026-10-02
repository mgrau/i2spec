"""The explorer's crossovers: the weak components the export carries, and the level identities rebuilt from
their indices (webapp.weak_links, app.js withLevels), find the same crossovers as the eigenstate labels."""
from types import SimpleNamespace

import pytest

from i2spec import saturation, webapp
from i2spec.model import RovibronicModel


@pytest.fixture(scope="module")
def model():
    return RovibronicModel("127I2")


def comp(offset, strength, upper, lower):
    return SimpleNamespace(label=None, offset=offset, strength=strength, upper_level=upper, lower_level=lower)


def rebuilt(main, total, weak):
    """app.js withLevels(): main component k has levels (k, k); a weak one the numbers the export gives it."""
    comps = [comp(c.offset, c.strength / total, k, k) for k, c in enumerate(main)]
    comps += [comp(o, s, iu, il) for o, s, iu, il in weak]
    return comps


def crossovers(components):
    return sorted((r.offset, r.amplitude, r.sharing) for r in saturation.resonances(components, 400.0)
                  if r.kind == "crossover")


@pytest.mark.parametrize("vu, vl, J, branch", [(43, 0, 10, "P"), (32, 0, 30, "R"), (5, 12, 3, "R")])
def test_crossovers_from_indices_match_the_eigenstates(model, vu, vl, J, branch):
    _, comps = model.hyperfine_components(vu, vl, J, branch, dJ=0)
    main = sorted((c for c in comps if c.label), key=lambda c: c.offset)
    total = sum(c.strength for c in main)
    weak = webapp.weak_links(comps, main, total)
    assert weak, "a low-J line has weak components above the threshold"
    # the same components with their own eigenstate labels, normalised as the export is
    chosen = [c for c in comps if not c.label and c.strength >= webapp.WEAK_FRACTION * total]
    assert len(chosen) == len(weak)
    true = [comp(c.offset, c.strength / total, c.upper_level, c.lower_level) for c in main + chosen]
    want, got = crossovers(true), crossovers(rebuilt(main, total, weak))
    assert want and len(got) == len(want)
    # offsets are exported to 0.1 MHz and weak strengths to 1e-6 of the line (1% at the threshold)
    for o1, a1, s1 in want:
        o2, a2, s2 = min(got, key=lambda r: (r[2] != s1, abs(r[0] - o1) + abs(r[1] - a1) / max(a1, 1e-12)))
        assert s1 == s2 and abs(o1 - o2) < 0.06 and abs(a1 - a2) <= 0.02 * a1 + 1e-9


def test_main_components_share_no_level(model):
    # what lets the export name a weak component's levels by the main components alone
    for vu, vl, J, branch in [(43, 0, 10, "P"), (32, 0, 56, "R"), (5, 12, 3, "R")]:
        _, comps = model.hyperfine_components(vu, vl, J, branch, dJ=0)
        main = [c for c in comps if c.label]
        assert len({c.upper_level for c in main}) == len({c.lower_level for c in main}) == len(main)


def test_high_J_lines_carry_no_weak_components(model):
    _, comps = model.hyperfine_components(32, 0, 56, "R", dJ=0)
    main = [c for c in comps if c.label]
    assert webapp.weak_links(comps, main, sum(c.strength for c in main)) == []
