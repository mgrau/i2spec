"""Exercise the dashboard's controls and live-file refresh without running fits."""
import asyncio
from datetime import datetime, timedelta, timezone
import json

import pytest
from rich.console import Console
from textual.widgets import DataTable, Static, TabbedContent

from i2spec import cli, research_progress as progress
from i2spec.progress_tui import ProgressDashboard


NOW = datetime(2026, 10, 11, 20, tzinfo=timezone.utc)
POTENTIAL, HYPERFINE = (item[0] for item in progress.CAMPAIGNS)
HELPER = progress.CAMPAIGNS[1][2]


def write(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data), encoding="utf-8")


def job(seed):
    return dict(phase="validate", cell=dict(id="candidate"), fold="physical_level:B:39",
                seed=seed, suffix=f"validation_s{seed}")


def finish(out, campaign, seed, accepted=True):
    row = job(seed)
    last = dict(free=61, inherited=11, score=dict(unique_absolute_comb=dict(rms_MHz=.389)),
                hyperfine_held_out=dict(rms_MHz=.108))
    write(out/campaign/"attempts"/(progress._tag(row, campaign)+".json"),
          dict(row, tag=progress._tag(row, campaign), accepted=accepted,
               started_utc=(NOW-timedelta(hours=30)).isoformat(),
               finished_utc=(NOW-timedelta(hours=1)).isoformat(), last=last))


class Monitor:
    def __init__(self, out):
        self.reader = progress.ProgressMonitor(out)
        self.calls = 0
        self.fail = False

    def snapshot(self, *, force):
        self.calls += 1
        if self.fail:
            raise OSError("Artifact directory temporarily unavailable")
        return self.reader.snapshot(now=NOW, force=force)


@pytest.fixture
def monitor(tmp_path):
    for campaign, seeds, workers in ((POTENTIAL, (0, 2), 6), (HYPERFINE, (0, 2, 4), 2)):
        write(tmp_path/campaign/"jobs_validate.json", [job(seed) for seed in seeds])
        write(tmp_path/campaign/"state.json", dict(status="running", workers=workers))
    write(tmp_path/HELPER/"state.json", dict(status="running", workers=4))
    write(tmp_path/HELPER/"launch.json", dict(started_utc=(NOW-timedelta(hours=18)).isoformat()))
    finish(tmp_path, POTENTIAL, 0)
    finish(tmp_path, HYPERFINE, 0, accepted=False)
    finish(tmp_path, HYPERFINE, 2)
    write(tmp_path/HYPERFINE/"comparison.json", dict(selected=[dict(cell=dict(id="candidate"),
          last=dict(free=61, inherited=11, score=dict(unique_absolute_comb=dict(rms_MHz=.156)),
                    hyperfine_held_out=dict(rms_MHz=.108)))]))
    for owner, seed in ((POTENTIAL, 2), (HELPER, 4)):
        row = job(seed)
        tag = progress._tag(row, owner)
        write(tmp_path/owner/"launches"/(tag+".json"),
              dict(row, tag=tag, started_utc=(NOW-timedelta(hours=2)).isoformat()))
        write(tmp_path/(tag+"_potential.log"), dict(evaluation=10, jacobian_evaluations=3))
    return Monitor(tmp_path)


async def settle(app, pilot):
    for _ in range(100):
        await pilot.pause(.02)
        if not app._loading:
            await pilot.pause(.05)
            return
    pytest.fail("Dashboard did not finish its background read")


def screen_text(app):
    return "\n".join(strip.text for strip in app.screen._compositor.render_strips())


def test_counts_navigation_filters_and_result_states(monitor):
    async def run():
        app = ProgressDashboard(monitor=monitor)
        async with app.run_test(size=(120, 40)) as pilot:
            await settle(app, pilot)
            assert app.snapshot_data["totals"]["remaining"] == 2
            assert "2 pending" not in screen_text(app)  # Both remaining fits have launched.
            assert "60.0%" in screen_text(app)
            assert "1 failed checks" in screen_text(app)
            views = app.query_one("#views", TabbedContent)
            await pilot.press("2")
            assert views.active == "workers"
            workers = app.query_one("#workers-table", DataTable)
            assert workers.row_count == 2
            await pilot.press("f")
            assert [row["campaign"] for row in app.worker_rows] == ["Potential models"]
            await pilot.press("f")
            assert [row["campaign"] for row in app.worker_rows] == ["Hyperfine allocations"]
            assert "start 4" in str(app.query_one("#worker-detail", Static).render())
            await pilot.press("f", "3")
            assert workers.row_count == 2 and views.active == "results"
            fold = app.query_one("#fold-results", DataTable).get_row_at(0)
            assert list(map(str, fold)) == ["B(v=39)", "Failed checks", "389"]
            assert app.query_one("#recent-results", DataTable).row_count == 3
            await pilot.press("4")
            assert views.active == "next"
            assert "Not queued" in screen_text(app)
            await pilot.press("1")
            assert views.active == "overview"
            await pilot.press("q")
        assert not app.is_running

    asyncio.run(run())


def test_pause_manual_refresh_selection_and_read_failure(monitor):
    async def run():
        app = ProgressDashboard(monitor=monitor)
        async with app.run_test(size=(100, 35)) as pilot:
            await settle(app, pilot)
            await pilot.press("2")
            workers = app.query_one("#workers-table", DataTable)
            workers.focus()
            await pilot.press("down")
            selected = app.worker_rows[workers.cursor_row]["tag"]
            await pilot.press("space")
            assert app.display_paused
            before = monitor.calls
            app._tick()
            await pilot.pause(.05)
            assert monitor.calls == before
            finish(monitor.reader.out, POTENTIAL, 2)
            await pilot.press("r")
            await settle(app, pilot)
            assert app.display_paused and monitor.calls == before+1
            assert app.snapshot_data["totals"]["remaining"] == 1
            assert app.worker_rows[workers.cursor_row]["tag"] == selected
            good_data = app.snapshot_data
            monitor.fail = True
            await pilot.press("r")
            await settle(app, pilot)
            assert app.snapshot_data is good_data
            assert "last successful snapshot" in str(app.query_one("#status", Static).render())
            monitor.fail = False
            await pilot.press("space")
            await settle(app, pilot)
            assert not app.display_paused and app._last_error is None
            assert "Auto-refresh" in str(app.query_one("#status", Static).render())

    asyncio.run(run())


@pytest.mark.parametrize("size", [(80, 24), (60, 20)])
def test_small_terminal_missing_data_and_resize(tmp_path, size):
    async def run():
        app = ProgressDashboard(out=tmp_path)
        async with app.run_test(size=size) as pilot:
            await settle(app, pilot)
            assert app.snapshot_data and not app.snapshot_data["available"]
            assert "Partial data" in screen_text(app)
            assert "unavailable" in str(app.query_one("#status", Static).render())
            await pilot.press("2")
            assert app.query_one("#workers-table", DataTable).row_count == 0
            assert "No unfinished launches" in str(app.query_one("#worker-detail", Static).render())
            await pilot.resize_terminal(140, 45)
            await pilot.press("3", "4", "1")
            assert "Queue files are unavailable" in screen_text(app)

    asyncio.run(run())


def test_cli_snapshot_and_interval_validation(monitor, capsys):
    out = str(monitor.reader.out)
    assert cli.main(["progress", "--out", out, "--json"]) == 0
    data = json.loads(capsys.readouterr().out)
    assert data["totals"]["finished"] == 3
    assert data["totals"]["failed"] == 1
    assert cli.main(["--json", "progress", "--out", out]) == 0
    assert json.loads(capsys.readouterr().out)["totals"]["total"] == 5
    for interval in ("0", "4", "nan", "inf"):
        assert cli.main(["progress", "--interval", interval, "--json"]) == 2
        assert "at least 5 seconds" in capsys.readouterr().err


def test_finalization_studies_are_visible_in_the_dashboard(monitor):
    out = monitor.reader.out
    write(out/"finalization_registry.json", dict(worker_budget=12, notes="Pilots require review", tracks=[
        dict(id="joint", label="Joint raw-data pilots", root="joint", kind="pilot", detail="Bounded pilots")]))
    write(out/"joint/state.json", dict(status="running", total=6, finished=2, active=1, workers=2,
        jobs=[dict(status="finished"), dict(status="queued")], active_jobs=[
            dict(id="pilot:fit", tag="pilot", stage="fit", fold="future")]))

    async def run():
        app = ProgressDashboard(monitor=monitor)
        async with app.run_test(size=(120, 40)) as pilot:
            await settle(app, pilot)
            assert "Joint raw-data pilots" in screen_text(app)
            assert "FOLLOW-UPS LEFT" in screen_text(app)
            assert "STUDIES ENDED" in screen_text(app)
            await pilot.press("2", "f", "f", "f")
            assert len(app.worker_rows) == 1
            assert app.worker_rows[0]["campaign"] == "Joint raw-data pilots"
            await pilot.press("4")
            assert "Bounded pilots" in screen_text(app)

    asyncio.run(run())


@pytest.mark.parametrize("size", [(120, 45), (80, 24)])
def test_completion_gates_separate_computation_from_fresh_measurements(monitor, size):
    out = monitor.reader.out
    # Complete all frozen queue attempts, including the existing failed check.
    finish(out, POTENTIAL, 2)
    finish(out, HYPERFINE, 4)
    write(out/"completion_state.json", dict(schema_version=1,
        status="computational_release_ready_awaiting_new_measurements",
        computational_release_ready=True, scientific_validation_complete=False, gates=[
            dict(id="numerical", label="Numerical and replay checks", status="complete", detail="All snapshots reproduce"),
            dict(id="blind", label="Fresh experimental validation", status="awaiting_measurements", detail="Acquire new frequencies"),
            dict(id="primary", label="Source evidence", status="requires_primary_source", detail="Obtain source documentation"),
            dict(id="coverage", label="Probability coverage", status="not_established", detail="95% coverage remains unestablished")]))

    async def run():
        app = ProgressDashboard(monitor=monitor)
        async with app.run_test(size=size) as pilot:
            await settle(app, pilot)
            assert app.snapshot_data["totals"]["progress_pct"] == 100
            assert "Computational release: ready" in screen_text(app)
            assert "Fresh experimental validation: still open" in screen_text(app)
            assert "Fresh experimental validation remains open" in str(app.query_one("#status", Static).render())
            await pilot.press("4")
            # All gate text is available in Next steps, including statuses that
            # do not mean a failed fit or a pending numerical worker.
            assert {r["status"] for r in app.snapshot_data["follow_up"]} == {
                "complete", "awaiting_measurements", "requires_primary_source", "not_established"}
            assert app.query_one("#views", TabbedContent).active == "next"
            console = Console(width=size[0])
            with console.capture() as capture:
                console.print(app.query_one("#follow-up", Static).content)
            gates = capture.get()
            for status in ("Complete", "Awaiting new measurements", "Requires primary source", "Not established"):
                assert status in gates
            assert "95% coverage remains unestablished" in gates
            assert "Not queued" not in gates

    asyncio.run(run())
