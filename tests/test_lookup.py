"""Line lookup, the command line and the TUI (src/i2spec/lookup.py, cli.py, tui.py)."""

import asyncio
import math

import numpy as np
import pytest

from i2spec import cli
from i2spec.intensity import MasterLineList
from i2spec.lookup import Catalog, Line, from_wavenumber, parse_label, parse_quantity, to_wavenumber, uncertainty


@pytest.fixture(scope="module")
def catalog():
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
    for unit in ("nm", "nm-air", "cm-1", "MHz", "THz"):
        nu = to_wavenumber(from_wavenumber(18788.3356, unit), unit)
        assert nu == pytest.approx(18788.3356, rel=1e-10)
    # the same number read as an air wavelength means a longer vacuum wavelength, so a smaller wavenumber
    assert to_wavenumber(532.245, "nm-air") < to_wavenumber(532.245, "nm")
    assert from_wavenumber(18788.3356, "nm-air") < from_wavenumber(18788.3356, "nm")


def test_parse():
    assert parse_quantity("532.245 nm") == (532.245, "nm")
    assert parse_quantity("18788.3cm-1") == (18788.3, "cm-1")
    assert parse_quantity("563260.2 GHz") == (563260200.0, "MHz")
    assert parse_quantity("532.1", "nm-air") == (532.1, "nm-air")
    assert parse_label("R(56) 32-0") == ("R", 56, 32, 0)
    assert parse_label("p12 5-3") == ("P", 12, 5, 3)
    for bad in ("532 furlongs", "R56"):
        with pytest.raises(ValueError):
            parse_quantity(bad) if " " in bad else parse_label(bad)


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


def test_uncertainty_estimates(catalog):
    fitted = Line("127I2", "R", 56, 32, 0, 18788.3356, 1e-18, 119.0, 293.15)
    edge = Line("127I2", "R", 56, 44, 0, 19500.0, 1e-18, 119.0, 293.15)
    beyond = Line("127I2", "R", 40, 58, 5, 18787.85, 1e-18, 1114.0, 293.15)
    red = Line("127I2", "P", 60, 12, 6, 14000.0, 1e-18, 500.0, 293.15)  # 714 nm, near the 671 nm anchor
    nir = Line("127I2", "P", 60, 0, 13, 12800.0, 1e-18, 500.0, 293.15)  # 781 nm, inside the measured 0-13 band
    nir_hot = Line("127I2", "R", 139, 1, 14, 12850.0, 1e-18, 500.0, 293.15)  # 778 nm, v′ = 1: inside the 1-14 band's lines
    nir_far = Line("127I2", "R", 60, 3, 16, 12750.0, 1e-18, 500.0, 293.15)  # v′ = 3: one line at J′ = 42; J′ = 61 is outside
    far = Line("127I2", "P", 60, 0, 18, 12000.0, 1e-18, 500.0, 293.15)  # 833 nm: X v'' = 18, measured by the Orsay atlas
    farther = Line("127I2", "P", 60, 0, 27, 10900.0, 1e-18, 500.0, 293.15)  # v'' = 27: beyond the atlas
    iso = Line("129I2", "R", 56, 32, 0, 18788.0, 1e-18, 119.0, 293.15)
    # v' = 32 at 532 nm: the B level is corrected from the BIPM lines, 0.68 MHz held out (i2spec2026c)
    assert uncertainty(fitted)[0] == 0.74 and fitted.flags == ("corrected levels",)
    low = Line("127I2", "R", 56, 20, 0, 17500.0, 1e-18, 119.0, 293.15)
    assert uncertainty(low)[0] == 3.0 and low.flags == ()
    # v' = 44 at J' = 57 is outside the yoshiki lines' J (33-39)
    assert uncertainty(edge)[0] == 15.0 and edge.flags == ("v′ > 43: extended model",)
    # v' = 58 at J' = 41 is inside the Orsay Partie IV lines of that level (J' 2-99, i2spec2026k): the level's
    # held-out figure; far beyond them in J the refitted potential is only an extrapolation
    assert 15.0 <= uncertainty(beyond)[0] <= 30.0
    assert beyond.flags == ("v′ 51-79: atlas-measured", "corrected levels")
    outside = Line("127I2", "R", 150, 58, 0, 19500.0, 1e-22, 1114.0, 293.15)
    assert uncertainty(outside)[0] == 1000.0 and outside.flags == ("v′ > 50: extrapolated",)
    assert uncertainty(red)[0] == 5.0 and "few data" in red.flags
    # 755-815 nm: i2spec2026c carries level corrections fitted to Liao 2010, Bodermann 2000 and the other
    # comb-referenced NIR lines, trusted within 15 in J of the measured range; the value is the levels' own
    # held-out rms (X v'' = 13: 1.2 MHz; B v' = 1 and X v'' = 14 together: several MHz, the 1-14 band's spread).
    assert 0.5 < uncertainty(nir)[0] < 2.0 and nir.flags == ("corrected levels",)
    # i2spec2026f: Bodermann's thesis adds beats between 1-14 and 0-12..0-14 lines, and B v' = 1 / X v'' = 14 drop
    # from several MHz held out to under 2 (the 1-14 band is no longer the worst-known part of the NIR)
    assert 1.0 < uncertainty(nir_hot)[0] < 3.0 and nir_hot.flags == ("corrected levels",)
    assert uncertainty(nir_far)[0] == 8.0 and nir_far.flags == ("few data",)
    # Reaching 833 nm requires v'' > 17: no B level lies low enough to get there from v'' <= 17.
    # i2spec2026g: the Orsay atlas part I measures X v'' = 18-25, so inside its J range such a line carries its
    # level's measured uncertainty (15-50 MHz); beyond v'' = 25 the v'' rule still applies.
    assert 10.0 < uncertainty(far)[0] < 30.0
    assert far.flags == ("v″ 18-25: atlas-measured", "corrected levels")
    assert uncertainty(farther)[0] == 300.0
    assert farther.flags == ("v″ > 17: unmeasured", "beyond 815 nm")
    assert "671 nm" in uncertainty(red)[1] and "level corrections" in uncertainty(nir)[1]
    assert "Orsay atlas" in uncertainty(far)[1] and "Orsay atlas" in uncertainty(farther)[1] and "v′ > 0" in uncertainty(nir_far)[1]
    assert uncertainty(iso)[0] == pytest.approx(math.hypot(5.0, 3.0), abs=0.01) and "isotope shift" in iso.flags
    assert "level corrections" in uncertainty(fitted)[1] and "Partie IV" in uncertainty(beyond)[1]


def test_components_match_the_bipm_reference(catalog):
    """The CIPM 532 nm reference: a10 of R(56) 32-0, which the model puts 2.0 MHz high."""
    line, = catalog.search(18788.3, 18788.4, "cm-1").lines
    comps = catalog.components(line, main_only=True)
    assert [c.label for c in comps] == [f"a{n}" for n in range(1, 16)]
    a10 = next(c for c in comps if c.label == "a10")
    assert a10.frequency_MHz == pytest.approx(563260223.513, abs=3.0)
    assert sum(c.strength for c in comps) == pytest.approx(1.0, abs=0.02)
    assert len(catalog.components(line)) > len(comps)  # the weak ΔF ≠ ΔJ components


def test_cli_table_and_parser(catalog, capsys):
    result = catalog.search(532.2, 532.3, "nm")
    text = cli._table(result.lines)
    assert "λ vac (nm)" in text and "R(56) 32-0" in text and "v′ 51-79: atlas-measured" in text
    assert len(text.splitlines()) == 4
    args = cli.build_parser().parse_args(["-i", "129I2", "lines", "532", "533", "--sort", "strength"])
    assert (args.isotopologue, args.low, args.high, args.sort) == ("129I2", "532", "533", "strength")
    after = cli.build_parser().parse_args(["lines", "532", "533", "-i", "129I2"])  # options work after the subcommand
    assert after.isotopologue == "129I2"
    bare = cli.build_parser().parse_args([])  # no arguments opens the browser
    assert bare.func is cli.cmd_tui and (bare.low, bare.high, bare.unit) == (None, None, "nm")
    assert cli.main(["line", "Q(1) 0-0"]) == 2  # a Q line cannot be parsed as a label
    assert "i2spec:" in capsys.readouterr().err


def test_tui_starts_and_shows_lines(catalog):
    pytest.importorskip("textual")
    from i2spec.tui import LineBrowser

    async def go():
        app = LineBrowser(catalog=catalog, low="532.2", high="532.3")
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
            for mode in ("line", "sub-Doppler", "window"):  # every plot mode produces a curve and draws it
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


def test_high_v_lower_is_flagged_as_unreliable():
    """Above v'' = 17 the published potentials extrapolate, and the measured error is enormous.

    The rule used to branch on v' alone, so a line like R(86) 33-48 was reported at +/-3 MHz when
    the model is 7.76 cm-1 (233 GHz) out there. See docs/design/fitting.md.
    """
    from i2spec.lookup import V_LOWER_FITTED, Line, uncertainty

    def line(v_upper, v_lower, nu):
        return Line("127I2", "R", 60, v_upper, v_lower, nu, 1e-20, 100.0, 300.0)

    assert V_LOWER_FITTED == 17
    ok = line(20, 0, 17500.0)
    assert uncertainty(ok)[0] == 3.0                       # unchanged where the data are (v' = 20, v'' = 0)
    assert not any("v″" in f for f in ok.flags)

    # i2spec2026d: the extended-range MLR X takes over above v'' = 17 -- 0.3 GHz where Martin 1986's levels
    # check it (to v'' = 28), 1 GHz beyond, 20 MHz at v'' = 48-54 where emission lines were measured
    for v_lower, expect in ((18, 300.0), (28, 300.0), (29, 1000.0), (47, 1000.0), (48, 20.0), (54, 20.0)):
        bad = line(33, v_lower, 10188.0)
        value, why = uncertainty(bad)
        assert value == expect, (v_lower, value)
        assert any("v″" in f for f in bad.flags)
        assert "v″" in why and len(why) > 60

    # a line whose levels are both measured (X v'' = 48 by emission, B v' = 58 by the Orsay atlas Partie IV,
    # i2spec2026k) carries their combined held-out figures
    assert uncertainty(line(58, 48, 19429.0))[0] == pytest.approx(21.0, abs=0.5)
