"""Terminal browser for iodine lines, built on i2spec.lookup and Textual."""

from __future__ import annotations

from math import log

import numpy as np
from textual.app import App, ComposeResult
from textual.containers import Horizontal
from textual.widgets import DataTable, Footer, Header, Input, Select, Static

from . import saturation
from .constants import ISOTOPOLOGUES
from .lookup import Catalog, Component, parse_range, to_wavenumber, uncertainty
from .spectrum import cross_section, doppler_fwhm

RESULT_COLUMNS = ("λ vac (nm)", "ν (cm⁻¹)", "f (THz)", "line", "S (cm)", "u (MHz)", "notes")
COMPONENT_COLUMNS = ("component", "offset (MHz)", "frequency (MHz)", "strength")
PLOT_MODES = ("line", "sub-Doppler")
#: Spectrum samples: at least this many, and at least four per Doppler width, capped for speed.
SPECTRUM_MIN_POINTS, SPECTRUM_MAX_POINTS = 4000, 400_000
#: Homogeneous width (MHz) drawn for sub-Doppler resonances: transit time plus pressure broadening
#: in a typical cell. saturation.transit_time_width() gives about 60 kHz for a 2 mm beam alone.
SUB_DOPPLER_FWHM = 2.0
#: Braille dot bits: [dot column][dot row within the cell].
DOTS = ((0x01, 0x02, 0x04, 0x40), (0x08, 0x10, 0x20, 0x80))


def braille_plot(y, width, height):
    """A filled curve drawn in braille: ``height`` strings of ``width`` characters.

    Each cell holds 2 x 4 dots, so the canvas is 2·width by 4·height. Samples are reduced to one
    value per dot column by taking the maximum, which keeps narrow lines from vanishing.
    """
    y = np.asarray(y, dtype=float)
    blank = [" " * width] * height
    if y.size == 0 or not np.isfinite(y).any():
        return blank
    columns, rows = 2 * width, 4 * height
    index = np.minimum(np.arange(y.size) * columns // y.size, columns - 1)
    peaks = np.zeros(columns)
    np.maximum.at(peaks, index, np.nan_to_num(y, nan=0.0))
    top = peaks.max()
    if not top > 0:
        return blank
    cells = [[0] * width for _ in range(height)]
    for c, value in enumerate(peaks):
        filled = int(round(max(value, 0.0) / top * (rows - 1)))
        for r in range(rows - 1 - filled, rows):
            cells[r // 4][c // 2] |= DOTS[c % 2][r % 4]
    return ["".join(chr(0x2800 + value) for value in row) for row in cells]


class LineBrowser(App):
    """Search a wavelength range, read the hyperfine structure of a line, and plot its spectrum."""

    CSS = """
    Screen { layout: vertical; }
    #query { height: 3; padding: 0 1; }
    #query Input { width: 1fr; }
    #query Select { width: 24; }
    #status { height: 1; padding: 0 2; color: $text-muted; }
    #spectrum { height: 14; padding: 0 1; border-bottom: solid $panel; }
    #results { height: 1fr; }
    #lower { height: 12; }
    #detail { width: 36; padding: 0 2; overflow-y: auto; }
    #components { width: 56; border-left: solid $panel; }
    #plot { width: 1fr; height: 100%; padding: 0 1; border-left: solid $panel; }
    """
    BINDINGS = [("q", "quit", "Quit"), ("slash", "focus_range", "Search"), ("s", "toggle_sort", "Sort"),
                ("a", "toggle_components", "All components"), ("p", "cycle_plot", "Line plot"), ("r", "search", "Refresh"),
                ("plus,equals_sign", "zoom(0.5)", "Zoom in"), ("minus", "zoom(2)", "Zoom out"),
                ("left_square_bracket", "pan(-0.25)", "◀"), ("right_square_bracket", "pan(0.25)", "▶"),
                ("l", "toggle_log", "Log"), ("c", "center", "Centre on line")]

    def __init__(self, catalog=None, isotopologue="127I2", temperature=293.15, low=None, high=None, unit="nm",
                 limit=200):
        super().__init__()
        self.catalog = catalog or Catalog(temperature=temperature)
        self.isotopologue = isotopologue
        self.temperature = float(temperature)
        self.unit = unit
        self.limit = limit
        self.sort = "nu"
        self.lines = []
        self.status_text = ""
        self.all_components = False
        self.current = None
        self.plot_mode = "line"
        self.components = []  # i2spec.lookup.Component, for the table
        self.raw = []  # i2spec.hyperfine.Component, which carry the nuclear spins the crossovers need
        self.window = None  # (nu_lo, nu_hi) of the last search, cm-1
        self.view = None    # (nu_lo, nu_hi) the spectrum panel shows
        self.spectrum = None  # (grid, sigma) for self.view
        self.log_scale = False
        self.range_text = f"{low}-{high}" if low and high else "532.2-532.3 nm"

    def compose(self) -> ComposeResult:
        yield Header()
        with Horizontal(id="query"):
            yield Input(value=self.range_text, placeholder="532.2-532.3 nm, or 18788-18789 cm-1", id="range")
            yield Select([(name, name) for name in sorted(ISOTOPOLOGUES)], value=self.isotopologue, id="iso",
                         allow_blank=False)
            yield Input(value=f"{self.temperature:g}", placeholder="T (K)", id="temperature")
        yield Static("", id="status")
        yield Static("", id="spectrum", markup=False)
        yield DataTable(id="results", cursor_type="row", zebra_stripes=True)
        with Horizontal(id="lower"):
            yield Static("", id="detail", markup=True)
            yield DataTable(id="components", cursor_type="row")
            yield Static("", id="plot", markup=False)
        yield Footer()

    def on_mount(self):
        self.title = "i2spec"
        self.sub_title = "iodine line lookup"
        self.query_one("#results", DataTable).add_columns(*RESULT_COLUMNS)
        self.query_one("#components", DataTable).add_columns(*COMPONENT_COLUMNS)
        self.action_search()

    # --- searching -------------------------------------------------------------------------------
    def action_search(self):
        text = self.query_one("#range", Input).value
        try:
            low, high, unit = parse_range(text, self.unit)
            self.unit = unit
            self.temperature = float(self.query_one("#temperature", Input).value)
        except ValueError as e:
            self.set_status(f"[red]{e}[/red]")
            return
        self.window = tuple(sorted(to_wavenumber(v, unit) for v in (low, high)))
        self.set_status(f"searching {self.isotopologue} …  (the first search after a model change builds the line "
                        f"list, about a minute)")
        self.run_worker(lambda: self._search(low, high, unit), thread=True, exclusive=True, group="search")

    def _search(self, low, high, unit):
        try:
            result = self.catalog.search(low, high, unit, self.isotopologue, temperature=self.temperature,
                                         limit=self.limit, sort=self.sort)
        except Exception as e:  # surfaced in the status line rather than crashing the app
            self.call_from_thread(self.set_status, f"[red]{e}[/red]")
            return
        self.call_from_thread(self.show_results, result)

    def show_results(self, result):
        self.lines = list(result)
        self.view = self.window
        self.refresh_spectrum()
        table = self.query_one("#results", DataTable)
        table.clear()
        for line in self.lines:
            table.add_row(f"{line.position('nm'):.5f}", f"{line.nu:.5f}", f"{line.position('THz'):.7f}", line.label,
                          f"{line.strength:.2e}", f"{line.uncertainty_MHz:,.0f}", ", ".join(line.flags))
        more = f", {result.total - len(self.lines)} more not shown" if result.truncated else ""
        self.set_status(f"{len(self.lines)} of {result.total} lines · {self.isotopologue} · "
                        f"{result.temperature:.2f} K · sorted by {'position' if self.sort == 'nu' else 'strength'}{more}")
        if self.lines:
            # open on the strongest line: in a wavelength range most rows are weak or badly placed hot bands
            row = max(range(len(self.lines)), key=lambda k: self.lines[k].strength)
            table.move_cursor(row=row)
            self.show_line(self.lines[row])

    # --- one line --------------------------------------------------------------------------------
    def on_data_table_row_highlighted(self, event):
        if event.data_table.id == "results" and 0 <= event.cursor_row < len(self.lines):
            self.show_line(self.lines[event.cursor_row])

    def show_line(self, line):
        self.current = line
        value, basis = uncertainty(line)
        self.query_one("#detail", Static).update(
            f"[b]{line.isotopologue} {line.label}[/b]\n"
            f"{line.position('nm'):.5f} nm vacuum · {line.position('nm-air'):.5f} nm air\n"
            f"{line.nu:.5f} cm⁻¹ · {line.position('MHz'):,.3f} MHz\n\n"
            f"S = {line.strength:.3e} cm at {line.temperature:.2f} K\n"
            f"E″ = {line.E_lower:.3f} cm⁻¹ · v′ = {line.v_upper}, v″ = {line.v_lower}, J″ = {line.J_lower}\n"
            f"Doppler FWHM {line.doppler_fwhm_MHz:.1f} MHz\n\n"
            f"[b]position uncertainty ≈ {value:,.0f} MHz[/b] (1σ estimate)\n{basis}")
        self.query_one("#components", DataTable).clear()
        self.run_worker(lambda: self._components(line), thread=True, exclusive=True, group="components")
        self.draw_spectrum()

    def _components(self, line):
        try:
            model = self.catalog.model(line.isotopologue)
            nu0, raw = model.hyperfine_components(line.v_upper, line.v_lower, line.J_lower, line.branch, dJ=2)
        except Exception as e:
            self.call_from_thread(self.set_status, f"[red]hyperfine: {e}[/red]")
            return
        self.call_from_thread(self.show_components, nu0, raw)

    def show_components(self, nu0, raw):
        self.raw = sorted(raw, key=lambda c: c.offset)
        self.components = [Component(c.label, c.offset, c.strength, c.F_upper, c.F_lower, nu0 + c.offset)
                           for c in self.raw if c.label or self.all_components]
        table = self.query_one("#components", DataTable)
        table.clear()
        for c in self.components:
            table.add_row(c.name, f"{c.offset_MHz:.3f}", f"{c.frequency_MHz:,.3f}", f"{c.strength:.4f}")
        self.call_after_refresh(self.draw_plot)  # the pane's height is only known once layout has settled

    # --- the spectrum panel ------------------------------------------------------------------------
    def refresh_spectrum(self):
        """Recompute the cross section over self.view in a worker; the panel redraws when it arrives."""
        if self.view is None:
            return
        view, T, iso = self.view, self.temperature, self.isotopologue
        self.run_worker(lambda: self._compute_spectrum(view, T, iso), thread=True, exclusive=True, group="spectrum")

    def _compute_spectrum(self, view, T, iso):
        lo, hi = view
        width = hi - lo
        step = min(width / SPECTRUM_MIN_POINTS, doppler_fwhm(0.5 * (lo + hi), T) / 4)
        step = max(step, width / SPECTRUM_MAX_POINTS)
        pad = 0.01 * width + 0.05
        lines = self.catalog.master(iso).at(T, lo - pad, hi + pad)
        grid = np.arange(lo, hi, step)
        sigma = cross_section(grid, lines, T)
        self.call_from_thread(self._spectrum_ready, view, grid, sigma, len(lines))

    def _spectrum_ready(self, view, grid, sigma, n):
        if view != self.view:
            return  # superseded by a later zoom or pan
        self.spectrum = (grid, sigma, n)
        self.draw_spectrum()

    def draw_spectrum(self):
        pane = self.query_one("#spectrum", Static)
        if self.spectrum is None or self.view is None:
            pane.update("")
            return
        grid, sigma, n = self.spectrum
        region = pane.content_size if pane.content_size.height > 3 else pane.size
        width = max(20, region.width)
        height = max(3, region.height - 3)  # title, marker row, axis
        y = sigma
        if self.log_scale:
            top = float(sigma.max()) if sigma.size else 0.0
            y = np.log10(np.maximum(sigma, top * 1e-6) / (top * 1e-6)) if top > 0 else sigma
        lo, hi = grid[0], grid[-1]
        marker = [" "] * width
        if self.current is not None and lo <= self.current.nu <= hi:
            marker[min(width - 1, int((self.current.nu - lo) / (hi - lo) * width))] = "▲"
        scale = "log σ, 6 decades" if self.log_scale else f"peak σ {sigma.max():.2e} cm²"
        title = (f"{self.isotopologue} {self.temperature:g} K · {1e7 / hi:.4f}–{1e7 / lo:.4f} nm "
                 f"({hi - lo:.3f} cm⁻¹) · {n} lines · {scale}")
        mid = 0.5 * (lo + hi)
        axis = f"{1e7 / lo:<.4f} nm"
        centre = f"{1e7 / mid:.4f}"
        right = f"{1e7 / hi:.4f} nm"
        gap = max(1, width - len(axis) - len(centre) - len(right))
        axis = axis + " " * (gap // 2) + centre + " " * (gap - gap // 2) + right
        pane.update("\n".join([title[:width], *braille_plot(y, width, height), "".join(marker), axis[:width]]))

    def action_zoom(self, factor):
        if self.view is None:
            return
        lo, hi = self.view
        centre = self.current.nu if self.current is not None and lo <= self.current.nu <= hi else 0.5 * (lo + hi)
        half = max(0.5 * (hi - lo) * float(factor), 0.002)
        self.view = (centre - half, centre + half)
        self.refresh_spectrum()

    def action_pan(self, fraction):
        if self.view is None:
            return
        lo, hi = self.view
        d = (hi - lo) * float(fraction)
        self.view = (lo + d, hi + d)
        self.refresh_spectrum()

    def action_center(self):
        if self.view is None or self.current is None:
            return
        lo, hi = self.view
        half = 0.5 * (hi - lo)
        self.view = (self.current.nu - half, self.current.nu + half)
        self.refresh_spectrum()

    def action_toggle_log(self):
        self.log_scale = not self.log_scale
        self.draw_spectrum()

    # --- one line's profile (lower pane) -----------------------------------------------------------
    def plot_data(self):
        """(y, title, left label, right label) for the current mode, or None when there is nothing to draw."""
        line = self.current
        if line is None or not self.raw:
            return None
        offsets = np.array([c.offset for c in self.raw])
        strengths = np.array([c.strength for c in self.raw])
        width = line.doppler_fwhm_MHz
        span = float(np.abs(offsets).max()) + 3 * width
        grid = np.linspace(-span, span, 1500)
        if self.plot_mode == "line":
            profile = strengths * np.exp(-4 * log(2) * ((grid[:, None] - offsets) / width) ** 2)
            return (profile.sum(axis=1), f"{line.label} · Doppler profile · {width:.0f} MHz FWHM",
                    f"−{span / 1000:.2f} GHz", f"+{span / 1000:.2f} GHz")
        resonances = saturation.resonances(self.raw, width)
        signal = saturation.signal(grid, resonances, SUB_DOPPLER_FWHM)
        crossovers = sum(r.kind == "crossover" for r in resonances)
        return (signal, f"{line.label} · sub-Doppler · {len(resonances) - crossovers} dips, {crossovers} crossovers",
                f"−{span / 1000:.2f} GHz", f"+{span / 1000:.2f} GHz")

    def draw_plot(self):
        pane = self.query_one("#plot", Static)
        data = self.plot_data()
        if data is None:
            pane.update("")
            return
        y, title, left, right = data
        region = pane.content_size if pane.content_size.height > 2 else pane.size
        width = max(20, region.width)
        height = max(3, region.height - 2)  # one row for the title, one for the axis
        axis = f"{left:<{width // 2}}{right:>{width - width // 2}}"
        pane.update("\n".join([title[:width], *braille_plot(y, width, height), axis]))

    def action_cycle_plot(self):
        self.plot_mode = PLOT_MODES[(PLOT_MODES.index(self.plot_mode) + 1) % len(PLOT_MODES)]
        self.draw_plot()

    def on_resize(self, event):
        self.draw_plot()
        self.draw_spectrum()

    def on_ready(self):
        self.draw_plot()  # the first layout pass finishes after on_mount, so size the canvases again here
        self.draw_spectrum()

    # --- small actions ---------------------------------------------------------------------------
    def set_status(self, text):
        self.status_text = text
        self.query_one("#status", Static).update(text)

    def action_focus_range(self):
        self.query_one("#range", Input).focus()

    def action_toggle_sort(self):
        self.sort = "strength" if self.sort == "nu" else "nu"
        self.action_search()

    def action_toggle_components(self):
        """Show every component, including the weak ΔF ≠ ΔJ ones, or only the main a1, a2, … series."""
        self.all_components = not self.all_components
        if self.current is not None:
            self.show_line(self.current)

    def on_input_submitted(self, event):
        # a new temperature changes the spectrum even when the range does not
        self.action_search()

    def on_select_changed(self, event):
        if event.value and event.value != self.isotopologue:
            self.isotopologue = str(event.value)
            self.action_search()


def run(catalog=None, isotopologue="127I2", temperature=293.15, low=None, high=None, unit="nm"):
    LineBrowser(catalog, isotopologue, temperature, low, high, unit).run()
