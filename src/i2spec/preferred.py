"""The preferred model: the i2spec model (release 20261021) for every isotopologue, i2spec2026o elsewhere.

Matt, 2026-10-21: balanced-expanded (research-20261020) predicts the lines it was not fitted to better than
i2spec2026o, so the line list, the terminal browser, the command line and the explorer prefer it; and "our goal is
to produce one single model". Release 20261021 (``research-20261021-i2spec``, role ``i2spec``) is balanced-expanded
unchanged for 127I2 up to J'' = 200, carried to 129I2 and 127I129I (fitted isotope shift, scaled hyperfine) and to
J'' = 445. For every line whose levels it holds (bound, or held behind their centrifugal barrier) the apps give its
position, hyperfine components and uncertainty (``ResearchModel.uncertainty``). i2spec2026o now gives only lines
the model does not hold, and every line strength until the line lists are rebuilt on the model's wavefunctions.

Every level of the domain is solved once per isotopologue, in parallel (a few minutes on 12 cores), and kept in
``cache_dir()`` with its statistical vector, data status and model sigma, keyed by the release files' SHA256:
positions and uncertainties of whole line lists are then table lookups.
"""
from __future__ import annotations

from functools import lru_cache
import hashlib
import json
import math
import os

import numpy as np


MODEL = "i2spec"
#: Internal names, for research notes and I2SPEC_MODEL only. The apps and the documentation present one model,
#: i2spec, and name neither (Matt, 2026-10-21): a line it does not hold is flagged OUTSIDE.
LABEL = "balanced-expanded"
FALLBACK = "i2spec2026o"
PUBLIC = "i2spec"
OUTSIDE = "outside the fitted domain"
ISOTOPOLOGUES = ("127I2", "129I2", "127I129I")
#: how a line's data status reads in the apps, from the worse of its two levels (ResearchModel.uncertainty)
STATUS = ("measured", "interpolated", "extrapolated")
_TABLE_VERSION = 2


def active():
    """Do the apps prefer the model? Yes, unless I2SPEC_MODEL names i2spec2026o."""
    return os.environ.get("I2SPEC_MODEL", MODEL).strip().lower() not in (FALLBACK, "2026o")


def parameters():
    """The parameter set the line lists (strengths, lower-level energies) are built on: the model's own release,
    or i2spec2026o under I2SPEC_MODEL."""
    if not active():
        from .constants import DEFAULT_PARAMETERS
        return DEFAULT_PARAMETERS
    from .research import _entry
    return _entry(MODEL)["id"]


@lru_cache(maxsize=None)
def model(isotopologue="127I2"):
    """The preferred research model (ResearchModel) of one isotopologue, loaded once."""
    from .research import ResearchModel
    return ResearchModel(MODEL, isotopologue)


def _key(isotopologue):
    from .research import _entry
    entry = _entry(MODEL)
    core = _entry(entry["extends"]["id"])["files_sha256"]
    text = json.dumps([_TABLE_VERSION, entry["id"], isotopologue, entry["files_sha256"]["isotopologues.json"],
                       core["model.json"], core["hyperfine.json"], core["uncertainty.json"]])
    return hashlib.sha256(text.encode()).hexdigest()[:16]


def _task(args):
    """Energies (cm-1) and uncertainty terms (ResearchModel.level_terms) of every level of one state at one J."""
    isotopologue, state, J = args
    m = model(isotopologue)
    vmax = m.info.domain.X_v_max if state == "X" else m.info.domain.B_v_max
    rows = []
    for v in range(vmax+1):
        try:
            E = m._levels_x.energy(state, v, J)
        except ValueError:                           # dissociated at this J
            continue
        a, status, near, sigma, beyond = m.level_terms(state, v, J)
        rows.append((v, E, a, STATUS.index(status), near, sigma, beyond))
    return state, J, rows


def _build(isotopologue, workers):
    import multiprocessing as mp
    from .intensity import _THREAD_VARS
    domain = model(isotopologue).info.domain
    tasks = [(isotopologue, "X", J) for J in range(domain.J_lower_max+1)]
    tasks += [(isotopologue, "B", J) for J in range(domain.J_lower_max+2)]
    n = workers or max(1, min(len(tasks), (os.cpu_count() or 2)-1))
    main = __import__("__main__")
    if n > 1 and getattr(main, "__file__", None):
        saved = {k: os.environ.get(k) for k in _THREAD_VARS}
        os.environ.update({k: "1" for k in _THREAD_VARS})
        try:
            with mp.get_context("spawn").Pool(n) as pool:
                results = pool.map(_task, tasks, chunksize=2)
        finally:
            for k, v in saved.items():
                os.environ.pop(k, None) if v is None else os.environ.__setitem__(k, v)
    else:
        results = [_task(t) for t in tasks]
    K = max(len(r[2]) for _, _, rows in results for r in rows)
    out = {}
    for state in ("X", "B"):
        vmax = domain.X_v_max if state == "X" else domain.B_v_max
        Jmax = domain.J_lower_max+(1 if state == "B" else 0)
        E = np.full((vmax+1, Jmax+1), np.nan)
        A = np.zeros((vmax+1, Jmax+1, K), np.float32)
        status = np.full((vmax+1, Jmax+1), -1, np.int8)
        near = np.zeros((vmax+1, Jmax+1), bool)
        sigma = np.zeros((vmax+1, Jmax+1, 3))
        beyond = np.zeros((vmax+1, Jmax+1), bool)
        for s, J, rows in results:
            if s != state:
                continue
            for v, e, a, st, nr, sg, by in rows:
                E[v, J], A[v, J], status[v, J], near[v, J], sigma[v, J], beyond[v, J] = e, a, st, nr, sg, by
        out.update({f"E_{state}": E, f"A_{state}": A, f"status_{state}": status, f"near_{state}": near,
                    f"sigma_{state}": sigma, f"beyond_{state}": beyond})
    return out


class Preferred:
    """Whole-list positions and uncertainties of the preferred model for one isotopologue, from its cached
    level table."""

    def __init__(self, isotopologue="127I2", workers=None, cache=True):
        from .intensity import cache_dir
        path = cache_dir()/f"preferred_{MODEL}_{isotopologue}_{_key(isotopologue)}.npz"
        if cache and path.exists():
            z = np.load(path)
            table = {k: z[k] for k in z.files}
        else:
            print(f"Solving every {isotopologue} level of the model once (a few minutes on 12 cores; cached in "
                  f"{path.parent}) ...", flush=True)
            table = _build(isotopologue, workers)
            if cache:
                path.parent.mkdir(parents=True, exist_ok=True)
                np.savez_compressed(path, **table)
        self.t = table
        self.isotopologue = isotopologue

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

    def lower_energies(self, J, vl):
        """Lower-level term values above X(0, 0) (cm-1), NaN where the table has no level."""
        J, vl = np.asarray(J), np.asarray(vl)
        EX = self.t["E_X"]
        ok = (vl >= 0) & (vl < EX.shape[0]) & (J >= 0) & (J < EX.shape[1])
        return np.where(ok, EX[np.where(ok, vl, 0), np.where(ok, J, 0)]-EX[0, 0], np.nan)

    def uncertainties(self, branch, J, vu, vl, chunk=200_000):
        """(sigma_center MHz, status index 0/1/2, class index 0 precise/1 atlas/2 beyond), NaN/-1 where not held."""
        iu, il, ok = self._index(branch, J, vu, vl)
        n = len(ok)
        stat = np.full(n, np.nan)
        for s in range(0, n, chunk):
            sl = slice(s, s+chunk)
            a = self.t["A_B"][iu[0][sl], iu[1][sl]].astype(float)-self.t["A_X"][il[0][sl], il[1][sl]]
            stat[sl] = np.sqrt(np.sum(a*a, axis=1))
        beyond = self.t["beyond_B"][iu] | self.t["beyond_X"][il]
        precise = self.t["near_B"][iu] & self.t["near_X"][il] & ~beyond
        cls = np.where(beyond, 2, np.where(precise, 0, 1))
        sm_b = np.take_along_axis(self.t["sigma_B"][iu], cls[:, None], axis=1)[:, 0]
        sm_x = np.take_along_axis(self.t["sigma_X"][il], cls[:, None], axis=1)[:, 0]
        sigma = np.sqrt(stat**2+sm_b**2+sm_x**2)
        status = np.maximum(self.t["status_B"][iu], self.t["status_X"][il])
        return np.where(ok, sigma, np.nan), np.where(ok, status, -1), np.where(ok, cls, -1)

    def covers(self, branch, J, vu, vl):
        if J > model(self.isotopologue).info.domain.J_lower_max:
            return False
        return bool(self._index([branch], [J], [vu], [vl])[2][0])


@lru_cache(maxsize=None)
def preferred(isotopologue="127I2") -> Preferred:
    """The process-wide Preferred table of one isotopologue."""
    return Preferred(isotopologue)


def basis(status, cls, isotopologue="127I2"):
    """The sentence the apps show for a preferred line's uncertainty (status, cls: Preferred.uncertainties codes)."""
    where = {0: "both levels are within 15 in J of comb or Doppler-free measurements at their v",
             1: "its levels are known from atlases and other lines (no precision data within 15 in J)",
             2: "a level lies beyond the core fit (X v″ > 17 or B v′ > 43)"}[int(cls)]
    data = {0: "the data constrain both levels directly",
            1: "a level is interpolated between data in J or v",
            2: "a level is extrapolated beyond the data in J or v: treat the value with care"}[int(status)]
    iso = "" if isotopologue == "127I2" else (
        "; the potentials are the ¹²⁷I₂ ones at this isotopologue's mass, with an isotope shift fitted to the seven "
        "CIPM isotopologue lines, and the fit's uncertainty is conditioned on those lines")
    return (f"The fit's statistical uncertainty through both levels, with the error "
            f"of predicting a level from the data around it fitted to held-out measurements; {where}; {data}{iso}")


def flags(status, cls):
    return (f"levels {STATUS[int(status)]}",)


def sigma_component(center_sigma, branch, J, vu, vl, offset, isotopologue="127I2"):
    """1-sigma (MHz) of a hyperfine component's frequency: the center's with the hyperfine floor and scale error."""
    Ju = J+(1 if branch == "R" else -1)
    a, kappa = model(isotopologue).hyperfine_terms(vu, Ju, vl, J)
    return math.sqrt(center_sigma**2+a**2+(kappa*offset)**2)


def components(branch, J, vu, vl, isotopologue="127I2"):
    """(center MHz, hyperfine components) of a preferred line, with the J +- 2 couplings."""
    return model(isotopologue).hyperfine_components(vu, vl, J, branch)


__all__ = ["MODEL", "LABEL", "FALLBACK", "PUBLIC", "OUTSIDE", "ISOTOPOLOGUES", "Preferred", "preferred", "active",
           "model", "parameters", "basis", "flags", "components", "sigma_component"]
