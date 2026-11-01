"""The core and continuation the fit stages build on, and which release's models take the roles.

By default the stages build on the releases the pipeline was first written for: the core of release 20261018
(core.release_spec) and the continuation of release 20261020 (continuation.release_spec), with the comb floor 0.5
MHz, and the models resolve by role to their first release. A rebuild of a later release names its own (`based_on`):

  core, continuation   specs in native coordinates (model.json records) that replace the released files
  comb_floor           the comb rows' weight floor of the core and continuation fits (core.FLOORS)
  s2_frozen            (state, curve constants) -> True where that perturber state's second-order terms are
                       computed on the B potential the release 20261021 hyperfine stage was fitted on
                       (hyperfine_beyond.potential) rather than on the current continuation
  s2_fixed             (state, curve constants) -> {(v, J): terms} that replace the computed second-order terms
                       of those levels (a release's own, where no committed input reproduces them), or None
  candidate            a release folder (its own manifest.json) whose models take over their ids and roles in
                       i2spec.research, so the stages that load the models by role (uncertainty, isotopologues, the
                       references) read the rebuild's

The worker processes of the parallel stages receive the same settings (`current`, passed to their initializer and
applied there with `use`).
"""
from __future__ import annotations

import copy
from contextlib import contextmanager

_STATE = dict(core=None, continuation=None, comb_floor=None, s2_frozen=None, s2_fixed=None)


def spec(name):
    """A copy of the named base spec ("core" or "continuation"), or None for the released file's."""
    value = _STATE[name]
    return copy.deepcopy(value) if value is not None else None


def s2_frozen(state, constants):
    """Whether this state's second-order terms are computed on the frozen 20261021 potential (None: the default
    rule of hyperfine_beyond.potential)."""
    rule = _STATE["s2_frozen"]
    return None if rule is None else bool(rule(state, constants))


def s2_fixed(state, constants):
    """{(v, J): terms} replacing the computed second-order terms of this state, or None."""
    rule = _STATE["s2_fixed"]
    return None if rule is None else rule(state, constants)


def current():
    """The settings a worker process needs (picklable: the rules must be module-level functions)."""
    return dict(_STATE)


def use(core=None, continuation=None, comb_floor=None, s2_frozen=None, s2_fixed=None):
    """Set every base (None restores that default)."""
    from . import core as core_module
    _STATE.update(core=core, continuation=continuation, comb_floor=comb_floor, s2_frozen=s2_frozen, s2_fixed=s2_fixed)
    core_module.FLOORS["comb"] = core_module.COMB_FLOOR if comb_floor is None else float(comb_floor)


@contextmanager
def based_on(core=None, continuation=None, comb_floor=None, s2_frozen=None, s2_fixed=None, candidate=None):
    """The stages inside build on these bases (and read the candidate release's models); the previous settings
    return on exit."""
    from .. import research
    before, folder = current(), research.candidate()
    use(core, continuation, comb_floor, s2_frozen, s2_fixed)
    research.use_candidate(candidate)
    try:
        yield
    finally:
        use(**before)
        research.use_candidate(folder)
