"""The documentation's numbers are the package's: every {{ key }} of a page resolves, the worked example's values agree
with the public interface, and its Python snippets return what their comments say."""
import importlib.util
import math
import re
from pathlib import Path

import pytest

DOCS = Path(__file__).resolve().parents[1] / "documentation"


def _load(name):
    spec = importlib.util.spec_from_file_location(f"i2spec_docs_{name}", DOCS / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def hooks():
    return _load("hooks")


@pytest.fixture(scope="module")
def numbers(hooks):
    return hooks.model_numbers()


def test_every_placeholder_resolves(hooks):
    for page in sorted(DOCS.glob("*.md")):
        filled = hooks.fill(page.read_text(encoding="utf-8"), page.name)
        assert not re.search(r"\{\{.*?\}\}", filled), page.name


def test_unknown_placeholder_stops_the_build(hooks):
    with pytest.raises(KeyError, match="no generated number"):
        hooks.fill("{{ we.no_such_number }}")


def test_worked_example_is_the_package(numbers):
    from i2spec import MHZ_PER_CM, ResearchModel, RovibronicModel

    w = numbers.worked_example()
    m, levels = ResearchModel("i2spec"), RovibronicModel("127I2")
    assert w["model_id"] == m.info.id
    nu = m.transition(32, 0, 56, "R")
    assert w["nu"] == nu and w["nu_MHz"] == nu * MHZ_PER_CM
    assert w["ex"] == levels.energy("X", 0, 56) and w["eb"] == levels.energy("B", 32, 57)
    assert abs((w["eb"] - w["ex"]) - nu) < 1e-8                            # the line is the difference of its levels
    assert w["comp_MHz"] == m.frequency("R(56) 32-0", "a10")
    center, components = m.hyperfine_components(32, 0, 56, "R")
    a10 = next(c for c in components if c.label == "a10")
    assert abs(center + a10.offset - w["comp_MHz"]) < 1e-6 and abs(center - nu * MHZ_PER_CM) < 1e-6
    assert w["n_main"] == 15 and w["b_eqq"] == m.level_hyperfine("B", 32, 57).eqQ
    u = m.uncertainty("R(56) 32-0", "a10")
    assert (w["sigma_comp"], w["sigma_hf"]) == (u.sigma_MHz, u.hyperfine_MHz)
    assert math.isclose(w["sigma"] ** 2, w["sigma_stat"] ** 2 + w["sigma_model"] ** 2)
    assert math.isclose(w["sigma_comp"] ** 2, w["sigma"] ** 2 + w["sigma_hf"] ** 2)
    assert math.isclose(w["sigma_interval"], math.sqrt(2 * w["hf_a"] ** 2 + (w["hf_kappa"] * w["span"]) ** 2))
    assert w["e_lower"] == pytest.approx(w["E_strength"], abs=1e-6)       # line list and level model agree on E''


def test_worked_strength_is_the_line_lists(numbers):
    from i2spec.lookup import Catalog

    w = numbers.worked_example()
    line = Catalog().line("R(56) 32-0", temperature=numbers.TEMPERATURE)
    assert line.strength == pytest.approx(w["S"], rel=1e-9, abs=0)      # ~1e-18: the default abs would pass anything
    assert line.E_lower == pytest.approx(w["E_strength"], abs=1e-9)
    assert line.doppler_fwhm_MHz == pytest.approx(w["doppler_MHz"], rel=1e-12)


def test_worked_measurement_was_held_out(numbers):
    w = numbers.worked_example()
    assert w["held_out_precision"] >= 1 and w["measured_set"] in w["held_out_sources"]
    assert abs(w["z"]) < 2                                                 # the model's 1 sigma covers it


def _snippets(text):
    return re.findall(r"```python\n(.*?)```", text, flags=re.S)


@pytest.mark.parametrize("page", ["calculation", "levels", "lines", "hyperfine", "uncertainty", "getting-started"])
def test_snippets_return_what_their_comments_say(hooks, page):
    """Each expression line of a page's Python snippets is evaluated; a comment that starts with a number must be
    its value, to the digits shown."""
    filled = hooks.fill((DOCS / f"{page}.md").read_text(encoding="utf-8"), page)
    if page == "getting-started":                           # its CLI-oriented first snippet needs no check
        filled = filled[filled.index("## Python interface"):]
    checked = 0
    for block in _snippets(filled):
        namespace = {}
        for raw in block.splitlines():
            code, _, comment = raw.strip().partition("#")
            code = code.strip()
            if not code or code.startswith(("print(", "Catalog().line")):
                continue
            try:
                value = eval(compile(code, page, "eval"), namespace)
            except SyntaxError:
                exec(compile(code, page, "exec"), namespace)
                continue
            number = re.match(r"\s*([-+]?[\d ]+\.?\d*(?:e[-+]?\d+)?)", comment.replace("−", "-"))
            if isinstance(value, float) and number:
                text = number[1].replace(" ", "")
                decimals = len(text.split(".")[1].split("e")[0]) if "." in text else 0
                tolerance = abs(float(text)) * 10 ** -3 if "e" in text else 0.5 * 10 ** -decimals
                assert abs(value - float(text)) <= tolerance * (1 + 1e-9), (page, code, value, text)
                checked += 1
    if page != "getting-started":
        assert checked >= 1, page


def test_docs_speak_of_one_model():
    for page in DOCS.glob("*.md"):
        text = page.read_text(encoding="utf-8")
        assert not re.search(r"2026[a-z]\b|balanced|research-20", text), page.name
