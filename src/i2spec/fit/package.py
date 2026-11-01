"""The release files from the fits: per-level energies, the hyperfine tables beyond the core, the B~1g resonance
tables and the reference predictions the packaged runtime must reproduce.

  levels      every (state, v, J <= 203) level of the final potentials, each on the fit engine's grid for its own v
              (box + 0.6 A, refined X tiers, the curves' joins), quasi-bound levels kept by the engine's barrier rule;
              a level a short tier puts at or above the limit is solved on the long tier (fit.levels)
  hyperfine   the constants per level beyond the core: B v' 44-80 (hyperfine_beyond: anchored background +
              second-order sums) and X v'' 18-108 (x_hyperfine), each level found on the hyperfine model's own grid by
              its energy (so no box or quasi-bound state can be mislabeled)
  resonance   B v' 77, 78: the constants without the resonant 1g(1Pig) level, the level's vibrational overlaps and
              couplings (b1g_resonance), per J
  references  line predictions: centers from `levels`; components from the core fit's constants on the core
              levels, the out-of-core offsets (reduction.OwnHyperfine) beyond, and the B~1g diagonalization at B v' 77,
              78, rotational spacings from `levels`
"""
from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor
from multiprocessing import get_context

import numpy as np

from ..bspline import BSplineSolver
from ..constants import HBAR2_2U, MHZ_PER_CM, reduced_mass
from ..hyperfine import level_structure, line_components
from ..observations import Line, component_rank
from . import b1g_resonance as rz, continuation, hyperfine_beyond as hb, x_hyperfine as xh
from .continuation import SolverPotential
from . import levels as levels_module
from .levels import TIERS
from .rows import CORE_V, EXCLUSION

J_MAX = 203                     # upper J <= 201 and the J + 2 neighbours of hyperfine mixing
V_MAX = dict(X=108, B=80)
BOX_PADDING = .6                # A: the continuation Problem's box padding


def tiers(joins):
    """The fit engine's tiers as the final fit used them: box + 0.6 A, X tiers >= 33 halved, the curves' joins."""
    return levels_module.tiers({s: tuple(j) for s, j in joins.items()}, BOX_PADDING, refine_x=True)


_W = {}


def _init_levels(values=None):
    curves = continuation.build_curves()
    values = values or continuation.release_spec()['values']
    vals = {s: dict(c.values, **values[s]) for s, c in curves.items()}
    _W['curves'], _W['values'] = curves, vals
    _W['tiers'] = tiers({s: c.joins for s, c in curves.items()})
    _W['solvers'] = {}


def kept_levels(solver, pot, grid, J, asym):
    """The engine's levels at J: the lowest nlev, and on a mesh (long box) only bound levels and quasi-bound
    resonances (>= 90% of the density inside the centrifugal barrier)."""
    E, psi = solver.wavefunctions(J)
    if 'mesh' not in grid:
        return E
    x = J*(J+1)
    R = solver.R
    w = solver.W[:, None]*psi**2
    veff = pot(R)+x*HBAR2_2U/reduced_mass('127I2')/R**2*(1+pot.nonadiabatic(R))
    outer = R > 4.5
    if veff[outer].max() > asym:
        Rb = R[outer][int(np.argmax(veff[outer]))]
        frac = (w*(R[:, None] < Rb)).sum(0)/w.sum(0)
        keep = (E < asym) | (frac >= .9)
    else:
        keep = E < asym
    return E[keep]


def _levels_task(args):
    state, tier, J = args
    curves, values = _W['curves'], _W['values']
    limit, grid = _W['tiers'][state][tier]
    key = state, tier
    pot = SolverPotential(curves[state], values[state])
    if key not in _W['solvers']:
        g = dict(grid)
        nlev, joins = g.pop('nlev'), g.pop('joins')
        _W['solvers'][key] = BSplineSolver(pot, reduced_mass('127I2'), joins=joins, nlev=nlev, **g)
    asym = float(pot(np.array([400.]))[0])
    E = kept_levels(_W['solvers'][key], pot, grid, J, asym)
    lo = 0 if tier == 0 else _W['tiers'][state][tier-1][0]+1
    if not levels_module.ABOVE_LIMIT_ON_LONG_GRID:
        return [(state, v, J, float(E[v])) for v in range(lo, min(limit, V_MAX[state])+1) if v < len(E)]
    # a level a short tier puts at or above the limit is the long tier's (its barrier rule): fit.levels' rule
    if 'mesh' in grid:
        vs = [v for v in range(min(limit, V_MAX[state])+1) if v < len(E) and (v >= lo or E[v] >= asym)]
    else:
        vs = [v for v in range(lo, min(limit, V_MAX[state])+1) if v < len(E) and E[v] < asym]
    return [(state, v, J, float(E[v])) for v in vs]


def levels(values=None, workers=12):
    """{(state, v, J): energy cm-1} of every level the potentials bind (or hold behind their barrier), J <= J_MAX,
    each on the fit engine's grid for its own v: the reference the packaged solver must reproduce."""
    jobs = [(s, t, J) for s, grids in TIERS.items() for t, (lim, _) in enumerate(grids)
            if (0 if t == 0 else grids[t-1][0]+1) <= V_MAX[s] for J in range(J_MAX+1)]
    jobs.sort(key=lambda j: -(j[1] if j[0] == 'X' else j[1]-1))     # the slow long-box tiers first
    with ProcessPoolExecutor(workers, mp_context=get_context('spawn'), initializer=_init_levels,
                             initargs=(values,)) as pool:
        rows = [r for part in pool.map(_levels_task, jobs, chunksize=2) for r in part]
    return {(s, v, J): E for s, v, J, E in sorted(rows)}


def structure_spec(parameters, structure=None):
    """The tails' structure (with the fit engine's grids and the curves' joins) and the fitted values, in
    research_expanded's form, from a continuation fit's parameters {"X:to1": ...}."""
    structure = dict(structure or continuation.release_spec()['structure'])
    curves = continuation.build_curves(structure)
    structure['grids'] = {s: [[limit, {k: v for k, v in g.items() if k != 'joins'}] for limit, g in t]
                          for s, t in tiers({s: c.joins for s, c in curves.items()}).items()}
    structure['joins'] = {s: list(c.joins) for s, c in curves.items()}
    values = {s: {k.split(':')[1]: v for k, v in parameters.items() if k.startswith(f'{s}:')} for s in 'XB'}
    return structure, values


def _match(E, state, v, J, E_grid):
    """Index of the hyperfine grid's level that is the model's (state, v, J) level: the nearest in energy, within a
    quarter of the local vibrational spacing (the hyperfine grids are coarser: up to 0.05 cm-1 off at high v). None
    when no level matches, or when neither v - 1 nor v + 1 exists at J (no spacing to judge the match by)."""
    e = E[(state, v, J)]
    gaps = [abs(E[(state, w, J)]-e) for w in (v-1, v+1) if (state, w, J) in E]
    if not gaps:
        return None
    spacing = min(gaps)
    d = np.abs(np.asarray(E_grid)-e)
    i = int(np.argmin(d))
    others = np.delete(d, i)
    if d[i] > .25*spacing or (len(others) and others.min() < .5*spacing):
        return None
    return i


def hyperfine_tables(E, b, theta, constants, xmodel, xlv, blv=None, workers=12, cache=None, orders=None):
    """B v' 44-80 and X v'' 18-108 constants per level (J <= J_MAX), on the hyperfine models' own grids (orders: the
    B background's where they are not the core fit's)."""
    report = {}
    blv = blv or hb.Levels()
    pairs, worst, skipped = [], 0., []
    for J in range(J_MAX+1):
        grid = blv.level(J)[0]
        for v in range(CORE_V["B"]+1, V_MAX["B"]+1):
            if ("B", v, J) not in E:
                continue
            i = _match(E, "B", v, J, grid)
            if i is None:
                skipped.append((v, J))
                continue
            worst = max(worst, abs(grid[i]-E[("B", v, J)]))
            pairs.append((v, J, i))
    terms = hb.fill_terms(sorted({(i, J) for _, J, i in pairs}), constants, workers=workers, cache=cache)
    model = hb.Model(blv, terms, orders=orders)
    rows_b = [[v, J, *map(float, model.constants(np.asarray(b), theta, i, J))] for v, J, i in pairs]
    report["B"] = dict(levels=len(rows_b), skipped=len(skipped), skipped_levels=skipped[:50],
                       index_shifted=sum(1 for v, _, i in pairs if i != v), max_energy_difference_cm=worst)
    rows_x, worst, skipped = [], 0., []
    for J in range(J_MAX+1):
        grid = xlv.level(J)[0]
        for v in range(CORE_V["X"]+1, V_MAX["X"]+1):
            if ("X", v, J) not in E:
                continue
            i = _match(E, "X", v, J, grid)
            if i is None or i <= CORE_V["X"]:
                skipped.append((v, J))
                continue
            worst = max(worst, abs(grid[i]-E[("X", v, J)]))
            rows_x.append([v, J, *map(float, xh.constants(xlv, xmodel, i, J))])
    report["X"] = dict(levels=len(rows_x), skipped=len(skipped), skipped_levels=skipped[:50],
                       max_energy_difference_cm=worst)
    return dict(B=rows_b, X=rows_x, report=report)


def resonance_table(E, res):
    """B v' 77, 78: constants without the resonant 1g(1Pig) level's term, its overlaps and couplings, per J."""
    ok = {v: [J for J in range(J_MAX+1) if ("B", v, J) in E and _match(E, "B", v, J, res.levels.level(J)[0]) == v]
          for v in rz.PIQUE}
    res.prepare({(v, J) for v, Js in ok.items() for J in Js})
    omega, parity, v0, v1, v2 = res.coup[rz.STATE]
    out = dict(state=rz.STATE, eqQ_1g_MHz=rz.G_EQQ, C_1g_MHz_JJ1=rz.G_C,
               source="Pique et al., J. Physique 47, 1917 (1986), sect. 3.3 (P86 II)", levels={})
    for v, Js in ok.items():
        p = rz.PIQUE[v]
        n = res.resonant_index(v, p["dE"])
        constants = {str(J): [float(x) for x in res.b_constants(v, J, n)] for J in Js}
        overlaps = {f"{Jg},{jb}": float(res.overlap(n, Jg, v, jb)) for jb in Js for Jg in range(max(1, jb-4), jb+5)}
        out["levels"][str(v)] = dict(model_level=n, dE_MHz=p["dE"], dE_sigma_MHz=p["dE_sigma"], B_MHz=p["B"],
                                     B_sigma_MHz=p["B_sigma"],
                                     coupling=[int(omega), parity, float(v0), float(v1), float(v2)],
                                     J=Js, constants=constants, overlaps=overlaps)
    return out


REFERENCE_LINES = dict(
    b_out=dict(vu=(44, 45, 47, 50, 53, 56, 58, 62, 66, 70, 74, 76, 79, 80), vl=(0, 2), J=(1, 9, 25, 60, 100, 150, 199)),
    x_out=dict(vu=(0, 10, 32, 43), vl=(18, 24, 33, 48, 54, 69, 85, 95, 100, 104, 108), J=(1, 8, 20, 40, 80, 150)),
    both=dict(vu=(50, 62, 70), vl=(20, 30, 60), J=(5, 30)),
    resonance=dict(vu=(77, 78), vl=(0, 1), J=(0, 1, 2, 5, 10, 26, 27, 40)))


def references(E, tables, resonance, own, res, base, core_references):
    """Predictions of the replay lines: centers from `levels`; components from the core fit's constants on the core
    levels (base: the core release's ResearchModel), `own` beyond, and the B~1g diagonalization at B v' 77, 78, all
    rotational spacings from `levels`. core_references: the core release's own reference predictions, recomputed."""
    have = {s: {(int(r[0]), int(r[1])) for r in tables[s]} for s in "XB"}
    rJ = {int(v): set(d["J"]) for v, d in resonance["levels"].items()}
    own.levels.energy = lambda v, J: E[("B", v, J)]
    own.xlv.energy_x00 = lambda v, J: E[("X", v, J)]-E[("X", 0, 0)]
    res.levels.energy = lambda v, J: E[("B", v, J)]

    def upper(branch, J):
        return J+(1 if branch == "R" else -1)

    def available(branch, J, vu, vl):
        Ju = upper(branch, J)
        for state, v, j0 in (("B", vu, Ju), ("X", vl, J)):
            for j in (j0-2, j0, j0+2):
                if j < 0:
                    continue
                if (state, v, j) not in E:
                    return False
                if state == "B" and v in rJ:
                    if j not in rJ[v]:
                        return False
                elif v > CORE_V[state]:
                    if (v, j) not in have[state]:
                        return False
                    grid = own.levels.level(j)[0] if state == "B" else own.xlv.level(j)[0]
                    if _match(E, state, v, j, grid) != v:
                        return False
        return True

    def core_structure(state, v, j0):
        e0 = E[(state, v, j0)]
        return level_structure(j0, lambda j: base._couplings(state, v, j), lambda j: (E[(state, v, j)]-e0)*MHZ_PER_CM,
                               symmetry="g" if state == "X" else "u", dJ=2)

    core_cache = {}

    def core_offsets(branch, J, vu, vl):
        """{label: MHz} of a core line from the core fit's constants (once per line)."""
        key = branch, J, vu, vl
        if key not in core_cache:
            Ju = upper(branch, J)
            comps = line_components(core_structure("B", vu, Ju), core_structure("X", vl, J), Ju, J)
            core_cache[key] = {c.label: c.offset for c in comps if c.label}
        return core_cache[key]

    lines = {tuple(ref["line"]) for ref in core_references}
    for spec in REFERENCE_LINES.values():
        for vu in spec["vu"]:
            for vl in spec["vl"]:
                for J in spec["J"]:
                    for branch in "PR":
                        if not (branch == "P" and J == 0) and available(branch, J, vu, vl):
                            lines.add((branch, J, vu, vl))
    needed = sorted({(vu, j) for b, J, vu, vl in lines if vu > CORE_V["B"] and vu not in rJ
                     for j in hb.neighbours(upper(b, J))})
    model = own.model(needed)
    res.prepare({(vu, j) for b, J, vu, vl in lines if vu in rJ for j in range(max(upper(b, J)-4, 0), upper(b, J)+5)
                 if j in rJ[vu]})
    out = []
    for branch, J, vu, vl in sorted(lines, key=lambda x: (x[2], x[3], x[1], x[0])):
        Ju = upper(branch, J)
        line = Line("127I2", branch, J, vu, vl)
        center = (E[("B", vu, Ju)]-E[("X", vl, J)])*MHZ_PER_CM
        if vu in rJ:
            s_up, _ = res.structure(vu, Ju)
            comps = line_components(s_up, own.x_structure(vl, J), Ju, J)
            offsets = {c.label: c.offset for c in comps if c.label}
        elif vu <= CORE_V["B"] and vl <= CORE_V["X"]:
            offsets = core_offsets(branch, J, vu, vl)
        else:
            if vu <= CORE_V["B"]:
                cb = {j: np.array([getattr(base._couplings("B", vu, j), p) for p in hb.PARAMS])
                      for j in hb.neighbours(Ju)}
            else:
                cb = {j: model.constants(own.b, own.theta, vu, j) for j in hb.neighbours(Ju)}
            offsets = {f"a{r}": o for r, o in own.offsets(line, cb).items()}
        labels = sorted(offsets, key=lambda c: int(c[1:]))
        out.append(dict(line=[branch, J, vu, vl], component=None, predicted_MHz=center))
        for c in dict.fromkeys([labels[0], labels[len(labels)//2], labels[-1]]):
            out.append(dict(line=[branch, J, vu, vl], component=c, predicted_MHz=center+offsets[c]))
        out.append(dict(line=[branch, J, vu, vl], component=labels[-1], ref_component=labels[0],
                        predicted_MHz=offsets[labels[-1]]-offsets[labels[0]]))
    centers = {tuple(r["line"]): r["predicted_MHz"] for r in out if r["component"] is None}
    rank = lambda c: f"a{component_rank(c)}"          # noqa: E731
    for ref in core_references:
        offsets = core_offsets(*ref["line"])
        row = dict(line=ref["line"], component=ref["component"],
                   predicted_MHz=centers[tuple(ref["line"])]+offsets[rank(ref["component"])], source_id=ref["id"])
        if "ref_component" in ref:
            row.update(ref_component=ref["ref_component"],
                       predicted_MHz=offsets[rank(ref["component"])]-offsets[rank(ref["ref_component"])])
        out.append(row)
    return out


def rounded(rows, digits=7):
    return [[int(r[0]), int(r[1]), *(round(float(x), digits) for x in r[2:])] for r in rows]


HF_FITTED = ["B theta_u (1u mixing)", "B gamma (eqQ continuation curvature)",
             "1g(1Pig) C5 scale (selected by training cost)", "X eqQ lambda", "X eqQ p", "X C slope"]


def inventory(values, curve_constants, base_fitted=61):
    """(inventory, detail) of the continued model: the core's fitted (61 in the release)/18/4, the tails, the
    hyperfine beyond the core."""
    new_terms = sorted({f"{s}:{n}" for s in "XB" for n in values[s]}-{"B:dDe"})
    n_perturber = sum(len([k for k in c if k in ("C5", "C6", "De", "re", "we")]) for c in curve_constants.values())
    n_inherited = 18+3+n_perturber+1+6
    inv = dict(fitted_molecular=base_fitted+len(new_terms)+len(HF_FITTED), inherited_molecular=n_inherited,
               fitted_calibration=5, line_specific=0, level_specific=0)
    return inv, new_terms, n_perturber


def model_record(structure, values, orsay4_offset, fit, curve_constants, counts, base_entry, derived=None):
    """model.json of the continued model (the release's schema 3). derived: the counts of the hyperfine
    coefficients set by inherited or core numbers (B background, X eqQ and C polynomials)."""
    base_fitted = base_entry.get("inventory", {}).get("fitted_molecular", 61)
    inv, new_terms, n_perturber = inventory(values, curve_constants, base_fitted)
    detail = dict(
        base=dict(model=base_entry["id"], fitted_molecular=base_fitted, inherited_molecular=18, fitted_calibration=4),
        potential_tails=dict(fitted=len(new_terms), names=new_terms,
                             note="B:dDe is X:dDe (one shared dissociation-energy change keeps the atomic splitting)"),
        hyperfine_fitted=dict(count=len(HF_FITTED), names=HF_FITTED,
                              note="fitted to measured B intra-line intervals out of core (B) and Bacis et al. 1980 "
                                   "X eqQ (X); the X C slope rests on one value (59 kHz at v''83) and is unvalidated"),
        hyperfine_derived=dict(derived or {}, note="anchored to the base model's own core predictions (B) and to BKT02 "
                                                   "in the core (X): set by inherited/base numbers, not by new "
                                                   "measurements"),
        inherited_new=dict(atomic=["A(2P3/2)", "A(2P1/2)", "B(2P3/2)"], perturber_curves=n_perturber,
                           theta_g="30 deg, the separated-atom quadrupole mixing (Saute 1982)",
                           resonance=["dE(77)", "B_1g(77)", "dE(78)", "B_1g(78)", "eqQ(1g)", "C(1g)"]),
        inherited_total=inv["inherited_molecular"],
        calibration=dict(base=4, orsay4_offset_MHz=orsay4_offset),
        structural=dict(knots=structure["knots"], orders=structure["orders"], window_knots=structure["window_knots"],
                        grids="the fit engine's tiers + 0.6 A box, X tiers >= 33 halved", J_tables=J_MAX))
    return dict(schema_version=3, id=base_entry.get("continued_id", "research-20261020-balanced-expanded"),
                kind="expanded", isotopologue="127I2",
                base=dict(release=base_entry["release"], model=base_entry["id"],
                          model_sha256=base_entry["files_sha256"]["model.json"]),
                core=dict(CORE_V), domain=dict(X_v_max=V_MAX["X"], B_v_max=V_MAX["B"], J_max=200,
                                                levels="a line needs its levels (and their J +- 2 partners) bound or "
                                                       "held behind the centrifugal barrier, and hyperfine components "
                                                       "need tabulated constants"),
                inventory=inv, inventory_detail=detail, structure=structure, values=values,
                calibration=dict(orsay4_offset_MHz=orsay4_offset,
                                 note="Orsay IV atlas offset; the base atlas calibrations are in the base model"),
                hyperfine=dict(core=f"{base_entry['id']}: B radial functions and fixed X constants",
                               beyond="hyperfine.json: B v'44-80 and X v''18-108 per level; B v'77, 78 with the "
                                      "1g(1Pig) level diagonalized", dJ=2),
                fit=fit, exclusion=EXCLUSION, training_fold="future", reference_counts=counts)


def hyperfine_record(tables, resonance, b_model, x_model, spin_spin):
    """hyperfine.json: the per-level tables (rounded to 7 digits), their models and the resonance."""
    return dict(
        B=dict(columns=["v", "J", "eqQ_MHz", "C_kHz", "d_kHz", "delta_kHz"], rows=rounded(tables["B"]), model=b_model),
        X=dict(columns=["v", "J", "eqQ_MHz", "C_kHz"], rows=rounded(tables["X"]), spin_spin=spin_spin, model=x_model),
        resonance=resonance, report=tables["report"], units=dict(eqQ="MHz", C="kHz", d="kHz", delta="kHz"))


def reference_counts(references):
    return dict(center=sum(1 for r in references if r["component"] is None),
                absolute=sum(1 for r in references if r["component"] and not r.get("ref_component")),
                intraline=sum(1 for r in references if r.get("ref_component")))
