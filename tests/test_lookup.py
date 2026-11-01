"""Line lookup, the command line and the TUI (src/i2spec/lookup.py, cli.py, tui.py)."""

import asyncio
import numpy as np
import pytest

from i2spec import cli
from i2spec.intensity import MasterLineList
from i2spec.lookup import (Catalog, Line, from_wavenumber, parse_label, parse_quantity, parse_range, to_wavenumber,
                           uncertainty)


@pytest.fixture
def catalog(hand_made):
    """A catalog holding four hand-made lines near 532 nm, so no line list has to be built."""
    nu = np.array([18788.33560, 18788.0, 18790.0, 18795.0])
    master = MasterLineList(nu=nu, strength0=np.array([4e-17, 1e-17, 2e-19, 3e-17]),
                            v_upper=np.array([32, 58, 32, 44]), v_lower=np.array([0, 5, 0, 1]),
                            J_lower=np.array([56, 40, 60, 12]), branch=np.array([1, 1, -1, 1]),
                            E_lower=np.array([119.052, 1114.38, 135.0, 250.0]),
                            x_levels=np.array([[0.0, 213.3], [0.075, 213.4]]), upper_cut=np.full(3, 19949.5),
                            isotopologue="127I2")
    return Catalog(masters={"127I2": master})


def test_units_round_trip():
    for unit in ("nm", "nm-air", "cm-1", "MHz", "GHz", "THz"):
        nu = to_wavenumber(from_wavenumber(18788.3356, unit), unit)
        assert nu == pytest.approx(18788.3356, rel=1e-10)
    # the same number read as an air wavelength means a longer vacuum wavelength, so a smaller wavenumber
    assert to_wavenumber(532.245, "nm-air") < to_wavenumber(532.245, "nm")
    assert from_wavenumber(18788.3356, "nm-air") < from_wavenumber(18788.3356, "nm")


def test_parse():
    assert parse_quantity("532.245 nm") == (532.245, "nm")
    assert parse_quantity("18788.3cm-1") == (18788.3, "cm-1")
    assert parse_quantity("563260.2 GHz") == (563260.2, "GHz")
    assert parse_quantity("532.1", "nm-air") == (532.1, "nm-air")
    assert parse_label("R(56) 32-0") == ("R", 56, 32, 0)
    assert parse_label("p12 5-3") == ("P", 12, 5, 3)
    for bad in ("532 furlongs", "R56"):
        with pytest.raises(ValueError):
            parse_quantity(bad) if " " in bad else parse_label(bad)
    for bad in ("0", "-532 nm", "0 MHz"):  # no wavelength or frequency is zero or negative
        with pytest.raises(ValueError, match="positive"):
            parse_quantity(bad)
    with pytest.raises(ValueError, match="positive"):
        to_wavenumber(0.0, "nm")


def test_parse_range():
    # a unit written once applies to both ends, and the "-" in "cm-1" is not a separator
    assert parse_range("532.2-532.3 nm") == (532.2, 532.3, "nm")
    assert parse_range("18788-18789 cm-1") == (18788.0, 18789.0, "cm-1")
    assert parse_range("18788cm-1-18789cm-1") == (18788.0, 18789.0, "cm-1")
    assert parse_range("532.2nm - 532.3nm") == (532.2, 532.3, "nm")
    assert parse_range("532.2 to 532.3 nm-air") == (532.2, 532.3, "nm-air")
    assert parse_range("532.2..532.3") == (532.2, 532.3, "nm")
    assert parse_range("18788 18789", "cm-1") == (18788.0, 18789.0, "cm-1")
    assert parse_range("563260-563261 GHz") == (563260.0, 563261.0, "GHz")
    for bad in ("532.2", "532.2-532.3-532.4", "532 nm - 18790 cm-1", "532.2-532.3 furlongs"):
        with pytest.raises(ValueError):
            parse_range(bad)


def test_search_range_sort_and_limit(catalog):
    result = catalog.search(532.2, 532.3, "nm")
    assert [line.label for line in result] == ["R(40) 58-5", "R(56) 32-0"]  # sorted by wavenumber
    assert result.total == 2 and not result.truncated
    assert catalog.search(532.2, 532.3, "nm", limit=1).truncated
    strongest = catalog.search(18780.0, 18800.0, "cm-1", sort="strength").lines[0]
    assert strongest.label == "R(56) 32-0"  # v″ = 0 beats the hot bands at 293 K
    assert catalog.search(532.3, 532.2, "nm").total == 2  # the range may be given either way round


def test_line_fields(catalog):
    line, = catalog.search(18788.3, 18788.4, "cm-1").lines
    assert line.position("nm") == pytest.approx(532.24512, abs=1e-5)
    assert line.position("MHz") == pytest.approx(563260131.1, abs=0.5)
    assert line.position("nm-air") < line.position("nm")
    assert line.J_upper == 57 and line.label == "R(56) 32-0"
    assert line.doppler_fwhm_MHz == pytest.approx(433.6, rel=0.01)
    assert 1e-19 < line.strength < 1e-17


def test_uncertainty_estimates():
    """Every line the level tables hold carries the model's 1σ and data status; the rest a rough bound."""
    fitted = Line("127I2", "R", 56, 32, 0, 18788.3356, 1e-18, 119.0, 293.15)
    value, why = uncertainty(fitted)
    assert 0.1 < value < 1.0 and fitted.flags == ("levels measured",) and fitted.source == "i2spec"
    assert why.startswith("The fit's statistical uncertainty")
    iso = Line("129I2", "R", 56, 32, 0, 18788.0, 1e-18, 119.0, 293.15)
    assert uncertainty(iso)[0] > value and "isotope shift" in uncertainty(iso)[1]
    hot = Line("127I2", "R", 230, 10, 2, 18000.0, 1e-20, 100.0, 300.0)       # J'' = 230: beyond the fit's J
    assert hot.flags == ("levels extrapolated",) and uncertainty(hot)[0] > 1.0
    edge = Line("127I2", "R", 20, 81, 0, 19500.0, 1e-22, 0.0, 293.15)         # B v' = 81: beyond the fitted domain,
    assert 1.0 < uncertainty(edge)[0] < 1e5 and edge.flags[0] == "outside the fitted domain"   # with the edge rule
    assert edge.source == "i2spec, outside the fitted domain" and "Danyluk" in uncertainty(edge)[1]
    outside = Line("127I2", "R", 20, 87, 0, 19500.0, 1e-22, 0.0, 293.15)      # B v' = 87: not in the tables
    assert uncertainty(outside)[0] == 1000.0 and outside.flags == ("outside the fitted domain",)
    assert outside.source == "i2spec, outside the fitted domain"


def test_components_match_the_bipm_reference():
    """The CIPM 532 nm reference: a10 of R(56) 32-0."""
    catalog = Catalog()
    line = Line("127I2", "R", 56, 32, 0, 0.0, 0.0, 0.0, 293.15)
    comps = catalog.components(line, main_only=True)
    assert [c.label for c in comps] == [f"a{n}" for n in range(1, 16)]
    a10 = next(c for c in comps if c.label == "a10")
    assert a10.frequency_MHz == pytest.approx(563260223.513, abs=0.5)
    assert sum(c.strength for c in comps) == pytest.approx(1.0, abs=0.02)
    assert len(catalog.components(line)) > len(comps)  # the weak ΔF ≠ ΔJ components


def test_cli_table_and_parser(catalog, capsys):
    result = catalog.search(532.2, 532.3, "nm")
    text = cli._table(result.lines)
    assert "λ vac (nm)" in text and "R(56) 32-0" in text and "outside the fitted domain" in text
    assert len(text.splitlines()) == 4
    args = cli.build_parser().parse_args(["-i", "129I2", "lines", "532", "533", "--sort", "strength"])
    assert (args.isotopologue, args.low, args.high, args.sort) == ("129I2", "532", "533", "strength")
    after = cli.build_parser().parse_args(["lines", "532", "533", "-i", "129I2"])  # options work after the subcommand
    assert after.isotopologue == "129I2"
    bare = cli.build_parser().parse_args([])  # no arguments opens the browser
    assert bare.func is cli.cmd_tui and (bare.low, bare.high, bare.unit) == (None, None, "nm")
    assert cli.main(["line", "Q(1) 0-0"]) == 2  # a Q line cannot be parsed as a label
    assert "i2spec:" in capsys.readouterr().err


def test_cli_range_units_and_errors(catalog, monkeypatch, capsys):
    import json

    monkeypatch.setattr(cli, "_catalog", lambda args: catalog)
    # a bare upper bound takes the lower bound's unit: 563 255-563 262 GHz, not 563 255 GHz down to 563 262 MHz
    for argv in (["563255GHz", "563262"], ["563255", "563262", "-u", "GHz"], ["563255 GHz", "563262000 MHz"]):
        assert cli.main(["lines", *argv, "--json"]) == 0, argv
        assert [d["line"] for d in json.loads(capsys.readouterr().out)] == ["R(56) 32-0"], argv
    assert cli.main(["lines", "0", "500"]) == 2  # was a ZeroDivisionError traceback
    assert "positive" in capsys.readouterr().err

    def missing(args):
        raise FileNotFoundError("no such file: data/observations")
    monkeypatch.setattr(cli, "_catalog", missing)
    assert cli.main(["lines", "532", "533"]) == 2
    assert "i2spec: no such file" in capsys.readouterr().err


def test_tui_starts_and_shows_lines(catalog):
    pytest.importorskip("textual")
    from i2spec.tui import PLOT_MODES, LineBrowser

    async def go():
        app = LineBrowser(catalog=catalog)  # the default range, "532.2-532.3 nm", read from the search box
        async with app.run_test() as pilot:
            await pilot.pause()
            for _ in range(50):  # the search and the hyperfine patterns run in worker threads
                if app.query_one("#components").row_count:
                    break
                await pilot.pause(0.1)
            assert app.query_one("#results").row_count == 2
            assert app.query_one("#components").row_count >= 15
            assert "2 of 2 lines" in app.status_text
            assert app.current.label == "R(56) 32-0"  # opens on the strongest line, not the first row
            for mode in PLOT_MODES:  # every plot mode produces a curve and draws it
                app.plot_mode = mode
                curve, title, left, right = app.plot_data()
                assert curve.max() > 0 and title and left and right, mode
                app.draw_plot()
                await pilot.pause()

    asyncio.run(go())


def test_braille_plot():
    from i2spec.tui import braille_plot

    rows = braille_plot(np.array([0.0, 1.0, 0.0]), width=10, height=3)
    assert len(rows) == 3 and all(len(r) == 10 for r in rows)
    assert rows[0].strip() and not rows[0].startswith("⣿")  # the peak reaches the top row, the edges do not
    assert braille_plot(np.array([]), 10, 3) == [" " * 10] * 3
    assert braille_plot(np.zeros(5), 10, 3) == [" " * 10] * 3  # a flat zero curve draws nothing
