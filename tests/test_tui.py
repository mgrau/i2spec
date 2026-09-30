"""The TUI, headless: a search fills the spectrum panel, and zoom, pan and log scale redraw it."""
import asyncio

import numpy as np
import pytest

textual = pytest.importorskip("textual")

from i2spec.intensity import MasterLineList  # noqa: E402
from i2spec.lookup import Catalog  # noqa: E402


class _Model:
    """Stands in for RovibronicModel: the hyperfine pattern is one component, no solve needed."""

    def hyperfine_components(self, v_upper, v_lower, J_lower, branch, dJ=2):
        from types import SimpleNamespace
        return 0.0, [SimpleNamespace(label="a1", offset=0.0, strength=1.0, F_upper=1, F_lower=1, I=None)]


def _catalog():
    nu = np.linspace(18787.0, 18790.0, 40)
    n = nu.size
    master = MasterLineList(nu=nu, strength0=np.full(n, 1e-20), v_upper=np.full(n, 32), v_lower=np.zeros(n, int),
                            J_lower=np.arange(n), branch=np.where(np.arange(n) % 2, 1, -1), E_lower=np.zeros(n),
                            x_levels=np.zeros((5, 3)), upper_cut=np.full(6, np.nan), isotopologue="127I2")
    return Catalog(masters={"127I2": master}, models={"127I2": _Model()})


def test_spectrum_panel_draws_and_zooms():
    from i2spec.tui import LineBrowser

    async def run():
        app = LineBrowser(_catalog(), low="532.25", high="532.30", unit="nm")
        async with app.run_test(size=(140, 45)) as pilot:
            for _ in range(50):
                await pilot.pause(0.05)
                if app.spectrum is not None:
                    break
            assert app.spectrum is not None and app.view is not None
            grid, sigma, n = app.spectrum
            assert n > 0 and sigma.max() > 0
            width0 = app.view[1] - app.view[0]
            text = str(app.query_one("#spectrum").render())
            assert "lines" in text and "nm" in text
            await pilot.press("plus")
            await pilot.pause(0.3)
            assert app.view[1] - app.view[0] == pytest.approx(width0 / 2)
            before = app.view
            await pilot.press("right_square_bracket")
            await pilot.pause(0.3)
            assert app.view[0] > before[0]
            await pilot.press("l")
            await pilot.pause(0.1)
            assert app.log_scale

    asyncio.run(run())
