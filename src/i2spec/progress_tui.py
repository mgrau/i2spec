"""Read-only terminal dashboard for the physical model search."""
from __future__ import annotations

import math

from rich.console import Group
from rich.panel import Panel
from rich.table import Table
from rich.text import Text
from textual.app import App, ComposeResult
from textual.binding import Binding
from textual.containers import VerticalScroll
from textual.widgets import DataTable, Footer, Header, Static, TabbedContent, TabPane

from .research_progress import OUT, ProgressMonitor, _date

GATE_LABELS = {"complete": "Complete", "awaiting_measurements": "Awaiting new measurements",
               "requires_primary_source": "Requires primary source", "not_established": "Not established"}


def completion_panel(completion):
    text = Text()
    ready = completion["computational_release_ready"]
    text.append("Computational release: " + ("ready" if ready else "pending"),
                style="bold green" if ready else "bold yellow")
    fresh = completion["scientific_validation_complete"]
    text.append("\nFresh experimental validation: " + ("complete" if fresh else "still open"),
                style="bold green" if fresh else "bold yellow")
    text.append("\nNumerical checks and frozen prediction replay are separate from fresh experimental evidence.", style="dim")
    return Panel(text, title="Completion gates", title_align="left", border_style="green" if fresh else "yellow")


def eta_text(eta):
    if eta is None:
        return "Awaiting timing data"
    if not eta["days_high"]:
        return "Finished"
    hours = eta["days_high"] < 1
    scale = 24 if hours else 1
    low = max(1, math.floor(eta["days_low"] * scale))
    high = max(low, math.ceil(eta["days_high"] * scale))
    value = str(low) if low == high else f"{low}–{high}"
    unit = ("hour" if hours else "day") + ("s" if high != 1 else "")
    return f"{value} {unit}"


def fold_label(value):
    if value.startswith("physical_level:"):
        _, state, level = value.split(":")
        return f"{state}(v={level})"
    return value.removeprefix("physical_source:")


def elapsed_text(hours):
    if hours is None:
        return "—"
    if hours < 1:
        return f"{int(hours * 60)} min"
    if hours < 24:
        return f"{hours:.1f} h"
    return f"{hours / 24:.1f} d"


def log_age(value, now):
    stamp = _date(value)
    if stamp is None:
        return "No log"
    seconds = max(0, (_date(now) - stamp).total_seconds())
    return "Just now" if seconds < 60 else elapsed_text(seconds / 3600) + " ago"


def _progress_bar(campaign, width=30):
    total = campaign["total"]
    text = Text()
    # Allocate every cell by cumulative boundaries so rounding cannot overflow.
    end = 0
    done = 0
    for key, style in (("accepted", "green"), ("failed", "yellow"), ("active", "cyan")):
        done += campaign[key]
        boundary = min(width, round(width * done / total)) if total else 0
        text.append("━" * max(0, boundary-end), style=style)
        end = boundary
    text.append("━" * (width-end), style="bright_black")
    return text


def queue_panel(c):
    if not c["available"]:
        return Panel(Text("Queue files are unavailable. Check the artifact directory.", style="yellow"), title=c["label"])
    lines = Text()
    lines.append(f"{c['finished']:,} / {c['total']:,} finished ({c['progress_pct']:.1f}%)", style="bold")
    lines.append(f"   ·   {c['remaining']:,} remaining\n")
    lines.append_text(_progress_bar(c))
    lines.append(f"   ETA {eta_text(c['eta'])}\n", style="bold cyan")
    for key, label, style in (("accepted", "accepted", "green"), ("failed", "failed checks", "yellow"),
                              ("active", "launched", "cyan"), ("queued", "pending", "default")):
        lines.append(f"{c[key]:,} {label}   ", style=style)
    rate = c["rates"][-1]
    rate_text = f"{rate['per_day']:.0f} finished/day over {rate['observed_hours']:.1f} h" if rate["per_day"] is not None else "Collecting completion rates"
    lines.append(f"\n{rate_text}   ·   {c['configured_workers']} workers configured", style="dim")
    return Panel(lines, title=c["label"], title_align="left", border_style="cyan")


def leader_panel(leader):
    if leader is None:
        return Panel("No accepted model comparison is available.", title="Absolute-accuracy leader")
    text = Text()
    text.append(leader["model"] + "\n", style="bold")
    text.append(f"{leader['absolute_rms_MHz'] * 1000:.1f} kHz absolute RMS", style="bold cyan")
    if leader["intraline_rms_MHz"] is not None:
        text.append(f"   ·   {leader['intraline_rms_MHz'] * 1000:.1f} kHz hyperfine-interval RMS")
    text.append(f"\n{leader['free']} fitted + {leader['inherited']} inherited molecular coefficients; "
                f"{leader['calibration']} fitted calibration terms additional.")
    validation = leader["validation"]
    if validation:
        text.append(f"\n{validation['accepted']} / {validation['total']} validation refits accepted", style="green")
    text.append("\nExploratory benchmark: 60 absolute components and 592 intervals. Accuracy measures can trade off.", style="dim")
    return Panel(text, title="Best absolute-line RMS so far", title_align="left", border_style="cyan")


def finalization_panel(data):
    table = Table(expand=True, box=None, padding=(0, 1))
    table.add_column("Study")
    table.add_column("Current stage")
    table.add_column("Ended", justify="right")
    table.add_column("Active", justify="right")
    for track in data["tracks"]:
        phase = {"combine": "Combination screens", "exact_screen": "Exact-J finalist refits", "validate": "Validation", "stability": "Training starts"}.get(track["phase"], track["phase"])
        table.add_row(Text(track["label"]), Text(phase + " · " + track["status"].replace("_", " ")),
                      f"{track['finished']} / {track['total']}", str(track["active"]))
    note = data["notes"]
    if data.get("review"):
        review = data["review"]
        numerical = "ready" if review["numerical_review_ready"] else "pending"
        note += f"\nFinal review: {review['status'].replace('_', ' ')} · numerical review: {numerical} · snapshot replay: {review['snapshot_verification']}."
    return Panel(Group(table, Text("\n"+note, style="dim")),
                 title="Model finalization · live follow-up studies", title_align="left", border_style="cyan")


class ProgressDashboard(App):
    """Display saved fit evidence; keyboard controls change only this display."""

    TITLE = "i2spec · Model optimization"
    CSS = """
    Screen { layout: vertical; }
    #status { height: auto; min-height: 1; max-height: 4; padding: 0 1; color: $text-muted; }
    #summary { height: auto; min-height: 5; padding: 0 1; border: round $primary; }
    #views { height: 1fr; }
    TabPane { padding: 0 1; }
    VerticalScroll { height: 1fr; }
    .report { height: auto; margin-bottom: 1; }
    #worker-filter { height: 2; }
    #workers-table { height: 1fr; min-height: 4; }
    #worker-detail { height: auto; min-height: 4; max-height: 8; padding: 1; border-top: solid $primary; }
    #fold-results { height: 11; }
    #recent-results { height: 15; }
    #daily-results { height: 10; }
    .section-label { height: 2; padding-top: 1; text-style: bold; color: $accent; }
    """
    BINDINGS = [
        ("q", "quit", "Quit"), ("r", "refresh", "Refresh"), ("space", "toggle_pause", "Pause updates"),
        Binding("1", "view('overview')", show=False), Binding("2", "view('workers')", show=False),
        Binding("3", "view('results')", show=False), Binding("4", "view('next')", show=False),
        ("f", "filter_workers", "Filter workers"),
    ]

    def __init__(self, out=OUT, interval=30, *, monitor=None):
        super().__init__()
        if not math.isfinite(interval) or interval < 5:
            raise ValueError("Use a refresh interval of at least 5 seconds")
        self.monitor = monitor or ProgressMonitor(out)
        self.interval = interval
        self.snapshot_data = None
        self.display_paused = False
        self._loading = False
        self._last_error = None
        self.worker_filter = 0
        self.worker_rows = []

    def compose(self) -> ComposeResult:
        yield Header(show_clock=True)
        yield Static("Reading local fitting results…", id="status", markup=False)
        yield Static("", id="summary")
        with TabbedContent(initial="overview", id="views"):
            with TabPane("1 Overview", id="overview"):
                with VerticalScroll():
                    yield Static("", id="completion", classes="report")
                    yield Static("", id="finalization", classes="report")
                    yield Static("", id="potential", classes="report")
                    yield Static("", id="hyperfine", classes="report")
                    yield Static("", id="leader-overview", classes="report")
                    yield Static("", id="estimate-note", classes="report", markup=False)
            with TabPane("2 Workers", id="workers"):
                yield Static("", id="worker-filter", markup=False)
                yield DataTable(id="workers-table", cursor_type="row", zebra_stripes=True)
                yield Static("Select a fit to see its details.", id="worker-detail", markup=False)
            with TabPane("3 Results", id="results"):
                with VerticalScroll():
                    yield Static("", id="leader-results", classes="report")
                    yield Static("Withheld bands and sources · independent starts, RMS in kHz", classes="section-label")
                    yield DataTable(id="fold-results", cursor_type="row", zebra_stripes=True)
                    yield Static("Recent validation completions", classes="section-label")
                    yield DataTable(id="recent-results", cursor_type="row", zebra_stripes=True)
                    yield Static("Daily finished attempts · UTC, today is partial", classes="section-label")
                    yield DataTable(id="daily-results", cursor_type="row", zebra_stripes=True)
            with TabPane("4 Next steps", id="next"):
                with VerticalScroll():
                    yield Static("", id="milestones", classes="report")
                    yield Static("", id="follow-up", classes="report")
                    yield Static("", id="counting-notes", classes="report", markup=False)
        yield Footer()

    def on_mount(self):
        self.sub_title = "read-only research progress"
        table = self.query_one("#workers-table", DataTable)
        for label, width in (("Model", 30), ("Withheld data", 20), ("Elapsed", 8), ("Last log", 12),
                             ("Stage", 30), ("Start", 5), ("Evaluation", 10), ("Campaign", 21)):
            table.add_column(label, width=width)
        self.query_one("#fold-results", DataTable).add_columns("Withheld data", "Start 0 RMS", "Start 2 RMS")
        self.query_one("#recent-results", DataTable).add_columns("Finished (local)", "Result", "Campaign", "Model", "Withheld data", "Start")
        self.query_one("#daily-results", DataTable).add_columns("Date (UTC)", "Finished", "Accepted", "Failed checks")
        self.set_interval(self.interval, self._tick)
        self.action_refresh()

    def _tick(self):
        if not self.display_paused:
            self.action_refresh()

    def action_refresh(self):
        if self._loading:
            return
        self._loading = True
        self.run_worker(self._load, thread=True, group="progress", exit_on_error=False)

    def _load(self):
        try:
            try:
                data = self.monitor.snapshot(force=True)
            except Exception as error:
                self.call_from_thread(self._show_error, str(error))
            else:
                self.call_from_thread(self.show_snapshot, data)
        except RuntimeError:
            # The user may have closed the app while a filesystem read finished.
            pass
        finally:
            self._loading = False

    def _show_error(self, message):
        self._last_error = message
        self._update_status()

    def _update_status(self):
        if self._last_error:
            message = ("Refresh failed; showing the last successful snapshot. " if self.snapshot_data
                       else "Could not read fitting results. ")
            self.query_one("#status", Static).update(message + self._last_error)
            return
        if not self.snapshot_data:
            return
        data = self.snapshot_data
        updated = _date(data["generated_utc"]).astimezone().strftime("%b %d %H:%M:%S %Z")
        state = "DISPLAY PAUSED — fits keep running" if self.display_paused else f"Auto-refresh every {self.interval:g}s"
        text = f"Updated {updated} · {state} · read-only"
        if not data["available"]:
            text += "\nSome queue files are unavailable; displayed totals are partial. Check --out."
        if data["warnings"]:
            text += "\n" + "  ".join(data["warnings"])
        completion = data.get("completion")
        if completion and not completion["scientific_validation_complete"]:
            text += "\nFresh experimental validation remains open; finished fit queues are computational progress."
        self.query_one("#status", Static).update(text)

    def show_snapshot(self, data):
        self.snapshot_data = data
        self._last_error = None
        totals = data["totals"]
        grid = Table.grid(expand=True, padding=(0, 2))
        for _ in range(4):
            grid.add_column(ratio=1)
        grid.add_row(Text("REFITS LEFT", style="dim"), Text("FINISHED", style="dim"),
                     Text("LOGS / WORKERS", style="dim"), Text("QUEUE ETA", style="dim"))
        pct = f"{totals['progress_pct']:.1f}%" if totals["progress_pct"] is not None else "—"
        grid.add_row(Text(f"{totals['remaining']:,}" if data["available"] else "Partial data", style="bold cyan"),
                     Text(pct, style="bold"), Text(f"{totals['recently_active']} / {totals['configured_workers']}", style="bold green"),
                     Text(eta_text(totals["eta"]) if data["available"] else "Unavailable", style="bold cyan"))
        grid.add_row(Text(f"{totals['queued']:,} pending", style="dim"),
                     Text(f"{totals['finished']:,} / {totals['total']:,}", style="dim"),
                     Text("Recent / set", style="dim"), Text("Current queues", style="dim"))
        self.query_one("#summary", Static).update(grid)
        completion = data.get("completion")
        completion_widget = self.query_one("#completion", Static)
        completion_widget.display = bool(completion)
        completion_widget.update(completion_panel(completion) if completion else "")
        finalization = data.get("finalization")
        self.query_one("#finalization", Static).update(finalization_panel(finalization) if finalization else "")
        if finalization:
            grid = Table.grid(expand=True, padding=(0, 2))
            for _ in range(4):
                grid.add_column(ratio=1)
            grid.add_row("FOLLOW-UPS LEFT", "STEPS ENDED", "LOGS / LIMIT", "STUDIES ENDED")
            grid.add_row(Text(str(finalization["remaining"]), style="bold cyan"),
                         Text(f"{finalization['finished']} / {finalization['total']}", style="bold"),
                         Text(f"{finalization['recently_active']} / {finalization['worker_budget']}", style="bold green"),
                         Text(f"{finalization['completed_tracks']} / {len(finalization['tracks'])}", style="bold"))
            grid.add_row("Current stages", "Includes pilot steps", "Recent / CPU cap",
                         "Numerical release ready" if completion and completion["computational_release_ready"] else "Review still required")
            self.query_one("#summary", Static).update(grid)
        for widget, campaign in zip(("potential", "hyperfine"), data["campaigns"]):
            self.query_one("#"+widget, Static).update(queue_panel(campaign))
        for widget in ("leader-overview", "leader-results"):
            self.query_one("#"+widget, Static).update(leader_panel(data["leader"]))
        self.query_one("#estimate-note", Static).update(
            ("Finalization counts cover each study's current stage; later validation can be conditional. "
             "Ended pilot steps do not mean final-model acceptance.\n" if finalization else "") +
            "Green: accepted · yellow: failed checks · cyan: launched · grey: pending.\n"
            "Failed attempts count as finished. Log activity is not a process inspection.\n"
            "ETA spans recent 6/12/24-hour rates at the current allocation, with 20% extra time on the slowest rate. "
            "It assumes continuous running and excludes recovery fits and unqueued studies.")
        self._render_workers()
        folds = {}
        for row in ((data["leader"] or {}).get("validation") or {}).get("folds", []):
            folds.setdefault(row["fold"], {})[row["seed"]] = row
        table = self.query_one("#fold-results", DataTable)
        table.clear()
        for fold, seeds in folds.items():
            scores = []
            for seed in (0, 2):
                row = seeds.get(seed)
                if row is None:
                    scores.append("Pending")
                elif not row["accepted"]:
                    scores.append(Text("Failed checks", style="yellow"))
                elif row["absolute_rms_MHz"] is None:
                    scores.append("Unavailable")
                else:
                    scores.append(f"{row['absolute_rms_MHz']*1000:,.0f}")
            table.add_row(Text(fold_label(fold)), *scores)
        recent = self.query_one("#recent-results", DataTable)
        recent.clear()
        for row in data["recent_completions"]:
            recent.add_row(_date(row["finished_utc"]).astimezone().strftime("%b %d %H:%M"),
                           Text("Accepted" if row["accepted"] else "Failed checks", style="green" if row["accepted"] else "yellow"),
                           Text(row["campaign"]), Text(row["model"]), Text(fold_label(row["fold"])), str(row["seed"]))
        daily = self.query_one("#daily-results", DataTable)
        daily.clear()
        for row in data["daily_completions"]:
            daily.add_row(row["date"], str(row["total"]), str(row["accepted"]), str(row["failed"]))
        milestones = Table(title="Earlier stages", expand=True, box=None)
        milestones.add_column("Stage")
        milestones.add_column("Finished / planned", justify="right")
        for row in data["milestones"]:
            milestones.add_row(Text(row["label"]), f"{row['finished']:,} / {row['total']:,}" if row["total"] else "Unavailable")
        self.query_one("#milestones", Static).update(Panel(milestones, border_style="green"))
        panels = []
        if finalization and not completion:
            for track in finalization["tracks"]:
                panels.append(Panel(Text(track["detail"]), title=Text(track["label"]),
                                    subtitle=Text(track["status"]), title_align="left", border_style="cyan"))
        for row in data["follow_up"]:
            status = GATE_LABELS.get(row["status"], row["status"]) if completion else row["status"]
            color = "green" if completion and row["status"] == "complete" else "yellow"
            panels.append(Panel(Text(row["detail"]), title=Text(row["title"]),
                                subtitle=Text(status), title_align="left", border_style=color))
        heading = ("Completion gates · computational release and fresh validation" if completion else
                   "Model-finalization studies" if finalization else "Beyond these queues — outside the progress percentage and ETA")
        self.query_one("#follow-up", Static).update(Group(Text(heading, style="bold"), *panels))
        self.query_one("#counting-notes", Static).update("\n".join(data["notes"]))
        self._update_status()

    def _render_workers(self):
        if not self.snapshot_data:
            return
        extra = (self.snapshot_data.get("finalization") or {}).get("tracks", [])
        filters = [None, "Potential models", "Hyperfine allocations"] + [track["label"] for track in extra]
        campaign = filters[self.worker_filter % len(filters)]
        table = self.query_one("#workers-table", DataTable)
        selected = self.worker_rows[table.cursor_row]["tag"] if 0 <= table.cursor_row < len(self.worker_rows) else None
        self.worker_rows = [row for c in [*self.snapshot_data["campaigns"], *extra] for row in c["active_jobs"]
                            if campaign is None or row["campaign"] == campaign]
        table.clear()
        for row in self.worker_rows:
            evaluation = row["evaluation"].get("evaluation")
            table.add_row(Text(row["model"]), Text(fold_label(row["fold"])), elapsed_text(row["elapsed_hours"]),
                          Text(log_age(row["last_activity_utc"], self.snapshot_data["generated_utc"]),
                               style="green" if row["activity"] == "recent" else "yellow"),
                          Text(row["stage"]), str(row["seed"]), str(evaluation) if evaluation is not None else "—",
                          Text(row["campaign"]), key=row["tag"])
        if self.worker_rows:
            index = next((i for i, row in enumerate(self.worker_rows) if row["tag"] == selected), 0)
            table.move_cursor(row=index)
            self._show_worker(index)
        else:
            self.query_one("#worker-detail", Static).update("No unfinished launches in this view.")
        self.query_one("#worker-filter", Static).update(
            f"Filter: {campaign or 'All campaigns'} · {len(self.worker_rows)} unfinished launches · f cycles filter\n"
            "Arrow keys select a fit; scroll horizontally to see all columns. Quiet logs do not prove a worker stopped.")

    def _show_worker(self, index):
        if not 0 <= index < len(self.worker_rows):
            return
        row = self.worker_rows[index]
        detail = (f"{row['model']}\n{row['campaign']} · withheld {fold_label(row['fold'])} · start {row['seed']}\n"
                  f"{row['stage']} · elapsed {elapsed_text(row['elapsed_hours'])} · "
                  f"last log {log_age(row['last_activity_utc'], self.snapshot_data['generated_utc'])}")
        evaluation = row["evaluation"]
        if evaluation.get("evaluation") is not None:
            detail += f"\nEvaluation {evaluation['evaluation']}"
            if evaluation.get("jacobian_evaluations") is not None:
                detail += f" · {evaluation['jacobian_evaluations']} Jacobian evaluations"
        self.query_one("#worker-detail", Static).update(detail)

    def on_data_table_row_highlighted(self, event):
        if event.data_table.id == "workers-table":
            self._show_worker(event.cursor_row)

    def action_toggle_pause(self):
        self.display_paused = not self.display_paused
        self._update_status()
        if not self.display_paused:
            self.action_refresh()

    def action_filter_workers(self):
        count = 3 + len(((self.snapshot_data or {}).get("finalization") or {}).get("tracks", []))
        self.worker_filter = (self.worker_filter+1) % count
        self._render_workers()

    def action_view(self, view):
        self.query_one("#views", TabbedContent).active = view


def run(out=OUT, interval=30):
    ProgressDashboard(out=out, interval=interval).run()
