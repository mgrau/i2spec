"""The preferred model: balanced-expanded where it applies, i2spec2026o elsewhere.

Matt, 2026-10-21: balanced-expanded (research-20261020) predicts the lines it was not fitted to better than
i2spec2026o, so the line list, the terminal browser, the command line and the explorer prefer it. For 127I2 lines
inside its domain (X v'' <= 108, B v' <= 80, J'' <= 200, both levels bound or held behind their centrifugal
barrier) they give its position, hyperfine components and uncertainty (``ResearchModel.uncertainty``).
i2spec2026o still gives every line strength, the other isotopologues and the 127I2 lines outside that domain, and
each line says which model gave it (``source``).

Every level of the domain is solved once, in parallel (80 s on 12 cores), and kept in ``cache_dir()`` with its
statistical factor and its data status, keyed by the release files' SHA256: positions and uncertainties of whole
line lists are then table lookups.
"""
from __future__ import annotations

from functools import lru_cache
import hashlib
import json
import math
import os

import numpy as np


MODEL = "balanced_expanded"
LABEL = "balanced-expanded"
FALLBACK = "i2spec2026o"
#: how a line's data status reads in the apps, from the worse of its two levels (ResearchModel.uncertainty)
STATUS = ("measured", "interpolated", "extrapolated")
_TABLE_VERSION = 1


def active():
    """Do the apps prefer balanced-expanded? Yes, unless I2SPEC_MODEL names i2spec2026o."""
    return os.environ.get("I2SPEC_MODEL", MODEL).strip().lower() not in (FALLBACK, "2026o")


@lru_cache(maxsize=1)
def model():
    """The preferred research model (ResearchModel), loaded once."""
    from .research import ResearchModel
    return ResearchModel(MODEL)


def _key():
    from .research import _entry
    entry = _entry(MODEL)
    files = entry["files_sha256"]
    text = json.dumps([_TABLE_VERSION, entry["id"], files["model.json"], files["hyperfine.json"],
                       files["uncertainty.json"]])
    return hashlib.sha256(text.encode()).hexdigest()[:16]


def _task(args):
    """Energies (cm-1), statistical factors (MHz) and data status of every level of one state at one J."""
    state, J = args
    m = model()
    unc = m._uncertainty_model()
    core = m._spec["core"]
    vmax = m.info.domain.X_v_max if state == "X" else m.info.domain.B_v_max
    rows = []
    for v in range(vmax+1):
        try:
            E = m._levels_x.energy(state, v, J)
        except ValueError:                           # dissociated at this J
            continue
        a = unc.level_factor(m._levels_x, state, v, J)
        status = STATUS.index(unc.status(state, v, J))
        near = unc._near("precise", state, v, J)
        sig = [unc.model_sigma(core, state, v, J, cls) for cls in ("precise", "atlas", "beyond")]
        rows.append((v, E, a, status, near, sig))
    return state, J, rows


def _build(workers):
    import multiprocessing as mp
    from .intensity import _THREAD_VARS
    m = model()
    tasks = [("X", J) for J in range(m.info.domain.J_lower_max+1)]
    tasks += [("B", J) for J in range(m.info.domain.J_lower_max+2)]
    n = workers or max(1, min(len(tasks), (os.cpu_count() or 2)-1))
    main = __import__("__main__")
    if n > 1 and getattr(main, "__file__", None):
        saved = {k: os.environ.get(k) for k in _THREAD_VARS}
        os.environ.update({k: "1" for k in _THREAD_VARS})
        try:
            with mp.get_context("spawn").Pool(n) as pool:
                results = pool.map(_task, tasks, chunksize=4)
        finally:
            for k, v in saved.items():
                os.environ.pop(k, None) if v is None else os.environ.__setitem__(k, v)
    else:
        results = [_task(t) for t in tasks]
    K = m._uncertainty_model().F.shape[1]
    out = {}
    for state in ("X", "B"):
        vmax = m.info.domain.X_v_max if state == "X" else m.info.domain.B_v_max
        Jmax = m.info.domain.J_lower_max+(1 if state == "B" else 0)
        E = np.full((vmax+1, Jmax+1), np.nan)
        A = np.zeros((vmax+1, Jmax+1, K), np.float32)
        status = np.full((vmax+1, Jmax+1), -1, np.int8)
        near = np.zeros((vmax+1, Jmax+1), bool)
        sigma = np.zeros((vmax+1, Jmax+1, 3))
        for s, J, rows in results:
            if s != state:
                continue
            for v, e, a, st, nr, sg in rows:
                E[v, J], A[v, J], status[v, J], near[v, J], sigma[v, J] = e, a, st, nr, sg
        out.update({f"E_{state}": E, f"A_{state}": A, f"status_{state}": status, f"near_{state}": near,
                    f"sigma_{state}": sigma})
    return out


class Preferred:
    """Whole-list positions and uncertainties of the preferred model, from its cached level table."""

    def __init__(self, workers=None, cache=True):
        from .intensity import cache_dir
        path = cache_dir()/f"preferred_{MODEL}_{_key()}.npz"
        if cache and path.exists():
            z = np.load(path)
            table = {k: z[k] for k in z.files}
        else:
            print(f"Solving every level of {LABEL} once (about 80 s on 12 cores; cached in {path.parent}) ...", flush=True)
            table = _build(workers)
            if cache:
                path.parent.mkdir(parents=True, exist_ok=True)
                np.savez_compressed(path, **table)
        self.t = table
        self.core = model()._spec["core"]

    # --- lookups -------------------------------------------------------------------------------------------

    def _index(self, branch, J, vu, vl):
        """Upper/lower array indices of lines (numpy), and the mask of lines whose levels the table holds."""
        branch, J, vu, vl = (np.asarray(x) for x in (branch, J, vu, vl))
        sign = np.where((branch == "R") | (branch == 1), 1, -1) if branch.dtype.kind in "UOS" else np.sign(branch)
        Ju = J+sign
        EB, EX = self.t["E_B"], self.t["E_X"]
        ok = ((vu >= 0) & (vu < EB.shape[0]) & (Ju >= 0) & (Ju < EB.shape[1]) & (vl >= 0) & (vl < EX.shape[0])
              & (J >= 0) & (J < EX.shape[1]) & ~((sign < 0) & (J == 0)))
        iu = (np.where(ok, vu, 0), np.where(ok, Ju, 0))
        il = (np.where(ok, vl, 0), np.where(ok, J, 0))
        ok &= np.isfinite(EB[iu]) & np.isfinite(EX[il])
        return iu, il, ok

    def positions(self, branch, J, vu, vl):
        """Wavenumbers (cm-1) of the lines, NaN where the preferred model has no line."""
        iu, il, ok = self._index(branch, J, vu, vl)
        return np.where(ok, self.t["E_B"][iu]-self.t["E_X"][il], np.nan)

    def uncertainties(self, branch, J, vu, vl, chunk=200_000):
        """(sigma_center MHz, status index 0/1/2, class index 0 precise/1 atlas/2 beyond), NaN/-1 where not held."""
        iu, il, ok = self._index(branch, J, vu, vl)
        n = len(ok)
        stat = np.full(n, np.nan)
        for s in range(0, n, chunk):
            sl = slice(s, s+chunk)
            a = self.t["A_B"][iu[0][sl], iu[1][sl]].astype(float)-self.t["A_X"][il[0][sl], il[1][sl]]
            stat[sl] = np.sqrt(np.sum(a*a, axis=1))
        beyond = (np.asarray(vu) > self.core["B"]) | (np.asarray(vl) > self.core["X"])
        precise = self.t["near_B"][iu] & self.t["near_X"][il] & ~beyond
        cls = np.where(beyond, 2, np.where(precise, 0, 1))
        sm_b = np.take_along_axis(self.t["sigma_B"][iu], cls[:, None], axis=1)[:, 0]
        sm_x = np.take_along_axis(self.t["sigma_X"][il], cls[:, None], axis=1)[:, 0]
        sigma = np.sqrt(stat**2+sm_b**2+sm_x**2)
        status = np.maximum(self.t["status_B"][iu], self.t["status_X"][il])
        return np.where(ok, sigma, np.nan), np.where(ok, status, -1), np.where(ok, cls, -1)

    def covers(self, isotopologue, branch, J, vu, vl):
        if isotopologue != "127I2" or J > model().info.domain.J_lower_max:
            return False
        return bool(self._index([branch], [J], [vu], [vl])[2][0])


@lru_cache(maxsize=1)
def preferred() -> Preferred:
    """The process-wide Preferred table."""
    return Preferred()


def basis(status, cls):
    """The sentence the apps show for a preferred line's uncertainty (status, cls: Preferred.uncertainties codes)."""
    where = {0: "both levels are within 15 in J of comb or Doppler-free measurements at their v",
             1: "its levels are known from atlases and other lines (no precision data within 15 in J)",
             2: "a level lies beyond the 2008 fit's range (X v″ > 17 or B v′ > 43)"}[int(cls)]
    data = {0: "the data constrain both levels directly",
            1: "a level is interpolated between data in J or v",
            2: "a level is extrapolated beyond the data in J or v: treat the value with care"}[int(status)]
    return (f"{LABEL} (research-20261020): the fit's statistical uncertainty through both levels, with the error "
            f"of predicting a level from the data around it fitted to held-out measurements; {where}; {data}")


def flags(status, cls):
    return (f"{LABEL}", f"levels {STATUS[int(status)]}")


def sigma_component(center_sigma, branch, J, vu, vl, offset):
    """1-sigma (MHz) of a hyperfine component's frequency: the center's with the hyperfine floor and scale error."""
    m = model()
    a, kappa = m._uncertainty_model().hyperfine_terms(m._spec["core"], vu, vl)
    return math.sqrt(center_sigma**2+a**2+(kappa*offset)**2)


def components(branch, J, vu, vl):
    """(center MHz, hyperfine components) of a preferred line, with the J +- 2 couplings."""
    return model().hyperfine_components(vu, vl, J, branch)


__all__ = ["MODEL", "LABEL", "FALLBACK", "Preferred", "preferred", "active", "model", "basis", "flags",
           "components", "sigma_component"]
