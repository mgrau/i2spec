"""Command line for i2spec: line lookup, terminal tools and the web app."""

from __future__ import annotations

import argparse
import json
import sys

from .constants import ISOTOPOLOGUES
from .lookup import UNITS, Catalog, format_mhz, from_wavenumber, parse_quantity, to_wavenumber, uncertainty

COLUMNS = (("λ vac (nm)", 12, "{:.5f}"), ("λ air (nm)", 12, "{:.5f}"), ("ν (cm⁻¹)", 12, "{:.5f}"),
           ("f (THz)", 11, "{:.7f}"), ("line", 14, "{}"), ("S (cm)", 10, "{:.2e}"),
           ("E″ (cm⁻¹)", 10, "{:.2f}"), ("u (MHz)", 9, "{}"), ("notes", 0, "{}"))


def _row(line):
    return [line.position("nm"), line.position("nm-air"), line.nu, line.position("THz"), line.label, line.strength,
            line.E_lower, format_mhz(line.uncertainty_MHz), ", ".join(line.flags)]


def _table(lines):
    head = "  ".join(f"{name:>{width}}" if width else name for name, width, _ in COLUMNS)
    rows = ["  ".join(f"{fmt.format(v):>{width}}" if width else fmt.format(v)
                      for v, (_, width, fmt) in zip(_row(line), COLUMNS)) for line in lines]
    return "\n".join([head, "-" * min(len(head), 120), *rows])


def _as_dict(line):
    value, basis = uncertainty(line)
    return dict(isotopologue=line.isotopologue, line=line.label, branch=line.branch, J_lower=line.J_lower,
                v_upper=line.v_upper, v_lower=line.v_lower, wavelength_vacuum_nm=line.position("nm"),
                wavelength_air_nm=line.position("nm-air"), wavenumber_cm1=line.nu, frequency_MHz=line.position("MHz"),
                strength_cm=line.strength, E_lower_cm1=line.E_lower, temperature_K=line.temperature,
                doppler_fwhm_MHz=line.doppler_fwhm_MHz, uncertainty_MHz=value, uncertainty_basis=basis,
                model=line.source, flags=list(line.flags))


def _catalog(args):
    return Catalog(temperature=args.temperature, s_min=args.cache_strength)


def cmd_lines(args):
    low, unit = parse_quantity(args.low, args.unit)
    high, high_unit = parse_quantity(args.high, unit)  # a bare upper bound takes the lower bound's unit
    if high_unit != unit:
        high = from_wavenumber(to_wavenumber(high, high_unit), unit)
    result = _catalog(args).search(low, high, unit, args.isotopologue, min_strength=args.min_strength,
                                   limit=args.limit, sort=args.sort)
    if args.json:
        print(json.dumps([_as_dict(line) for line in result], indent=1))
        return
    print(f"{len(result)} of {result.total} lines of {args.isotopologue} at {result.temperature:.2f} K"
          f"{', strongest first' if args.sort == 'strength' else ''}")
    print(_table(result))
    if result.truncated:
        print(f"... {result.total - len(result)} more; raise --limit or --min-strength")


def cmd_line(args):
    catalog = _catalog(args)
    line = catalog.line(args.label, args.isotopologue)
    value, basis = uncertainty(line)
    try:
        comps = catalog.components(line, dJ=args.dJ, main_only=not args.all)
    except ValueError as error:                       # the model's reason it gives the center only
        comps, note = [], str(error)
    else:
        note = catalog.center_only(line, comps)
    if args.json:
        print(json.dumps(dict(_as_dict(line), components=[
            dict(label=c.label, offset_MHz=c.offset_MHz, frequency_MHz=c.frequency_MHz, strength=c.strength,
                 F_upper=c.F_upper, F_lower=c.F_lower, uncertainty_MHz=c.uncertainty_MHz) for c in comps],
            **({"center_only": note} if note else {})), indent=1))
        return
    print(f"{line.isotopologue} {line.label}   {line.position('nm'):.5f} nm vacuum   {line.nu:.5f} cm⁻¹   "
          f"{line.position('MHz'):,.3f} MHz")
    print(f"  S = {line.strength:.3e} cm at {line.temperature:.2f} K; E″ = {line.E_lower:.3f} cm⁻¹; "
          f"Doppler FWHM {line.doppler_fwhm_MHz:.1f} MHz")
    print(f"  model: {line.source}")
    print(f"  position uncertainty ≈ {format_mhz(value)} MHz (1σ): {basis}")
    if line.flags:
        print(f"  flags: {', '.join(line.flags)}")
    if note:
        print(f"\n  no hyperfine components: {note}")
        return
    with_u = any(c.uncertainty_MHz is not None for c in comps)
    print(f"\n  {'component':<12}{'offset (MHz)':>14}{'frequency (MHz)':>20}{'u (MHz)':>10}{'strength':>10}" if with_u
          else f"\n  {'component':<12}{'offset (MHz)':>14}{'frequency (MHz)':>20}{'strength':>10}")
    for c in comps:
        u = f"{format_mhz(c.uncertainty_MHz):>10}" if with_u else ""
        print(f"  {c.name:<12}{c.offset_MHz:>14.4f}{c.frequency_MHz:>20.4f}{u}{c.strength:>10.4f}")


def cmd_tui(args):
    from .tui import run

    run(isotopologue=args.isotopologue, temperature=args.temperature, low=args.low, high=args.high, unit=args.unit)


def cmd_model(args):
    from .preferred import model

    m = model(args.isotopologue)
    if args.action == "verify":
        print(json.dumps(m.verify(), indent=2))
        return
    info = m.info
    if args.json:
        print(json.dumps(info.as_dict(), indent=2))
        return
    ledger, domain = info.inventory, info.domain
    print(f"i2spec model, release {info.id.split('-')[1]}: {domain.isotopologue}, X v″ ≤ {domain.X_v_max}, "
          f"B v′ ≤ {domain.B_v_max}, J″ ≤ {domain.J_lower_max}")
    print(f"  {ledger.fitted_molecular} fitted and {ledger.inherited_molecular} inherited molecular coefficients, "
          f"{ledger.fitted_calibration} calibration coefficients; none specific to a line or level")
    print(f"  held-out comb data: {info.absolute_rms_MHz * 1000:.0f} kHz rms (absolute components), "
          f"{info.intraline_rms_MHz * 1000:.0f} kHz (hyperfine intervals)")


def cmd_web(args):
    from . import webapp

    if args.action == "inputs":                     # the digest of what an export reads (the workflow's cache key)
        print(webapp.export_inputs())
        return
    # before any export, which takes minutes: the files the command reads exist only in a source checkout
    webapp.check_checkout(site=args.action != "export", docs=args.action != "export")
    if args.action == "export":
        webapp.export(s_min=args.s_min)
        return
    # the site is built from data exported from the current model and data; serving locally, older data will do
    stale = webapp.has_data() and webapp.data_inputs() != webapp.export_inputs()
    if args.rebuild or not webapp.has_data() or (stale and args.action == "build"):
        if stale and not args.rebuild:
            print("the model or its data changed since the last export: exporting again ...", flush=True)
        elif args.action == "serve":
            print("exporting the model for the web app first (a few minutes the first time) ...", flush=True)
        webapp.export(s_min=args.s_min)
    else:
        if stale:
            print("i2spec: note: the model or its data changed since the last export; `i2spec web export` writes "
                  "the data again", file=sys.stderr, flush=True)
        cutoff = webapp.data_cutoff()
        if cutoff is not None and cutoff != args.s_min:
            print(f"i2spec: note: the exported data has --s-min {cutoff:g}, not {args.s_min:g}; "
                  f"add --rebuild to export again", file=sys.stderr, flush=True)
    if args.action == "build":
        print(f"site written to {webapp.build(args.dest)}")
    else:
        site = webapp.build(webapp.ROOT / "site")        # the explorer and the documentation, as published
        webapp.serve(port=args.port, open_browser=not args.no_browser, directory=site,
                     host="0.0.0.0" if args.lan else args.host)


def _common(parser, suppress=False):
    """Options accepted both before and after the subcommand; after it, they win."""
    default = (lambda value: argparse.SUPPRESS) if suppress else (lambda value: value)
    parser.add_argument("--isotopologue", "-i", default=default("127I2"), choices=sorted(ISOTOPOLOGUES))
    parser.add_argument("--temperature", "-T", type=float, default=default(293.15),
                        help="cell temperature in K (default 293.15)")
    parser.add_argument("--cache-strength", type=float, default=default(1e-27),
                        help="weakest line to keep in the cached list (cm); a smaller value means a slower first run")
    parser.add_argument("--json", action="store_true", default=default(False), help="print JSON instead of a table")
    return parser


def build_parser():
    parser = _common(argparse.ArgumentParser(
        prog="i2spec", description="Look up molecular iodine B-X lines. With no arguments, opens the terminal browser; `i2spec web` opens the web app."))
    # no subcommand opens the TUI, so that plain `i2spec` is the interactive tool
    parser.set_defaults(func=cmd_tui, command="tui", low=None, high=None, unit="nm")
    sub = parser.add_subparsers(dest="command")

    p = _common(sub.add_parser("lines", help="list the lines in a wavelength or frequency range"), suppress=True)
    p.add_argument("low")
    p.add_argument("high")
    p.add_argument("--unit", "-u", default="nm", choices=UNITS, help="unit of low/high when they carry none")
    p.add_argument("--min-strength", type=float, default=0.0, help="ignore lines weaker than this (cm)")
    p.add_argument("--limit", type=int, default=50)
    p.add_argument("--sort", choices=("nu", "strength"), default="nu")
    p.set_defaults(func=cmd_lines)

    p = _common(sub.add_parser("line", help="one line with its hyperfine components"), suppress=True)
    p.add_argument("label", help="for example 'R(56) 32-0'")
    p.add_argument("--all", action="store_true", help="include the weak ΔF ≠ ΔJ components")
    p.add_argument("--dJ", type=int, default=2, choices=(0, 2), help="rotational mixing in the hyperfine matrix")
    p.set_defaults(func=cmd_line)

    p = _common(sub.add_parser("tui", help="interactive line browser (needs Textual)"), suppress=True)
    p.add_argument("--low", default=None)
    p.add_argument("--high", default=None)
    p.add_argument("--unit", "-u", default="nm", choices=UNITS)
    p.set_defaults(func=cmd_tui)

    p = _common(sub.add_parser("model", help="the model's domain and coefficients, or replay its frozen predictions"),
                suppress=True)
    p.add_argument("action", nargs="?", default="info", choices=("info", "verify"))
    p.set_defaults(func=cmd_model)

    p = sub.add_parser("web", help="the explorer and its documentation: serve locally, export data, build the site")
    p.add_argument("action", nargs="?", default="serve", choices=("serve", "export", "build", "inputs"),
                   help="serve (default); export the data; build the site, exporting first if the model or its data "
                        "changed; inputs: print the digest of what an export reads")
    p.add_argument("dest", nargs="?", default="site", help="build: the directory to write the site to")
    p.add_argument("--port", type=int, default=8777)
    p.add_argument("--no-browser", action="store_true", help="serve without opening a browser")
    p.add_argument("--host", default="127.0.0.1", help="address to serve on (default: this machine only)")
    p.add_argument("--lan", action="store_true", help="serve on the local network (same as --host 0.0.0.0)")
    p.add_argument("--rebuild", action="store_true", help="re-export the data even if it exists")
    p.add_argument("--s-min", type=float, default=1e-24,
                   help="weakest unmeasured line to export, cm at 300 K (measured lines are always kept)")
    p.set_defaults(func=cmd_web)
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    try:
        args.func(args)
    except (ValueError, LookupError, OSError) as e:  # OSError: a missing or unreadable data file
        print(f"i2spec: {e}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
