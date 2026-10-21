"""Explicitly selected, frozen 127I2 research models.

The October 18, 2026 release holds four models with shared MLR potentials and
radial hyperfine functions, without line/level corrections. The October 20
release continues its balanced model to X v<=108 and B v<=80
(research_expanded). Their data are installed with i2spec; the fitting
workspace and observation archive are not runtime dependencies. Production
continues to use i2spec2026o. Benchmark errors describe reused exploratory
validation, not per-prediction uncertainties.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass, replace
from functools import lru_cache
import hashlib
from importlib import resources
import json
import math
from numbers import Integral
from threading import RLock

import numpy as np
from numpy.polynomial.chebyshev import chebval
from numpy.polynomial.polynomial import polyvander

from .bspline import BSplineSolver
from .constants import MHZ_PER_CM, reduced_mass
from .hyperfine import HyperfineParameters, level_structure, line_components
from .lookup import parse_label
from .observations import component_rank

_RELEASE = "20261018"
_RELEASES = (_RELEASE, "20261020", "20261021")


def _data(release=_RELEASE):
    return resources.files("i2spec").joinpath("data", "research", release)


@lru_cache(maxsize=None)
def _manifest(release=_RELEASE):
    return json.loads(_data(release).joinpath("manifest.json").read_text())


def _entries():
    """Every release's models in release order, each entry tagged with its release."""
    return [dict(value, release=release) for release in _RELEASES if _data(release).joinpath("manifest.json").is_file()
            for value in _manifest(release)["models"]]


def _entry(name):
    for value in _entries():
        if name in (value["id"], value["role"]):
            return value
    choices = ", ".join(value["id"] for value in _entries())
    raise ValueError(f"unknown research model {name!r}; choose from {choices}")


def _read(entry, filename):
    folder = _data(entry["release"]) if entry.get("release", _RELEASE) != _RELEASE else _data()
    raw = folder.joinpath(entry["directory"], filename).read_bytes()
    if hashlib.sha256(raw).hexdigest() != entry["files_sha256"][filename]:
        raise ValueError(f"frozen research artifact failed its SHA256 check: {entry['id']}/{filename}")
    return json.loads(raw)


@dataclass(frozen=True)
class CoefficientInventory:
    fitted_molecular: int
    inherited_molecular: int
    fitted_calibration: int
    line_specific: int = 0
    level_specific: int = 0


@dataclass(frozen=True)
class ResearchDomain:
    isotopologue: str = "127I2"
    X_v_max: int = 17
    B_v_max: int = 43
    J_lower_max: int = 200

    def validate(self, branch, J, vu, vl):
        if branch not in ("P", "R"):
            raise ValueError("research models support P and R lines")
        if any(isinstance(n, bool) or not isinstance(n, Integral) or n < 0 for n in (J, vu, vl)):
            raise ValueError("quantum numbers must be nonnegative integers")
        if (branch == "P" and J == 0) or J > self.J_lower_max or vu > self.B_v_max or vl > self.X_v_max:
            raise ValueError(f"research domain: {self.isotopologue}, X v<={self.X_v_max}, B v<={self.B_v_max}, "
                             f"J_lower<={self.J_lower_max}, P J>=1")


@dataclass(frozen=True)
class ResearchModelInfo:
    id: str
    role: str
    title: str
    architecture: str
    inventory: CoefficientInventory
    domain: ResearchDomain
    exclusion: str
    absolute_rms_MHz: float
    intraline_rms_MHz: float
    validated_folds: int
    verification_max_absolute_Hz: float
    verification_max_intraline_Hz: float
    training_fold: str = "future"
    validation_scope: str = "Reused exploratory benchmark and seven overlapping physical holdouts"

    def as_dict(self):
        return asdict(self)


def _info(entry, isotopologue=None):
    isotopologue = isotopologue or entry["isotopologue"]
    domain = entry["domains"][isotopologue] if "domains" in entry else entry["domain"]
    errors = entry["verification"]["maximum_difference_Hz"]
    return ResearchModelInfo(entry["id"], entry["role"], entry["title"], entry["architecture"],
        CoefficientInventory(**entry["inventory"]),
        ResearchDomain(isotopologue, domain["X_v_max"], domain["B_v_max"], domain["J_max"]),
        entry["exclusion"], entry["absolute_rms_MHz"], entry["intraline_rms_MHz"],
        entry["validated_folds"], errors["absolute"], errors["intraline"],
        **{key: entry[key] for key in ("training_fold", "validation_scope") if key in entry})


@lru_cache(maxsize=None)
def _reference_model(model_id):
    """The 127I2 model of a release, which the other isotopologues' hyperfine and equivalents are read from."""
    return ResearchModel(model_id, "127I2")


def research_models():
    """Metadata for every selectable research model (all releases), without radial solves."""
    return tuple(_info(entry) for entry in _entries())


@dataclass(frozen=True)
class ResearchPrediction:
    model: str
    isotopologue: str
    line: str
    component: str | None
    frequency_MHz: float
    wavenumber_cm1: float
    excluded_from_fit: bool

    def as_dict(self):
        return asdict(self)


@dataclass(frozen=True)
class PredictionUncertainty:
    """1-sigma of a balanced-expanded prediction and where it comes from (MHz).

    ``sigma_MHz`` is the total: statistical and model terms in quadrature for a line center, plus the hyperfine
    term for a component; an intraline interval carries the hyperfine term alone. ``status`` is the worse of the
    two levels': measured, interpolated or extrapolated (see ``ResearchModel.uncertainty``)."""
    model: str
    line: str
    component: str | None
    ref_component: str | None
    sigma_MHz: float
    statistical_MHz: float
    model_MHz: float
    hyperfine_MHz: float
    status: str
    upper_status: str
    lower_status: str
    line_class: str
    hyperfine_status: str | None

    def as_dict(self):
        return asdict(self)


class _Potential:
    """The selected frozen MLR equations, including their affine coordinates."""
    def __init__(self, state, spec):
        self.state = state
        self.values = spec["values"][state]
        self.config = spec["config"]
        self.C = {int(n): c for n, c in spec["inherited_potential_coefficients"][state].items()}

    def __call__(self, R):
        R = np.asarray(R, float)
        cfg, values, state = self.config, self.values, self.state
        re, de = np.exp(values["logRe"]), np.exp(values["logDe"])
        n = cfg["nx" if state == "X" else "nb"]
        p, q = (cfg["px"], cfg["qx"]) if state == "X" else (cfg["pb"], cfg["qb"])
        ref = re * cfg["ref_x" if state == "X" else "ref_b"]
        y = (R**p-ref**p)/(R**p+ref**p)
        yq = (R**q-ref**q)/(R**q+ref**q)
        yp = (R**p-re**p)/(R**p+re**p)
        coordinate = (yq-.2)/.6 if state == "X" and cfg["mapped_x"] else yq
        exponent = chebval(coordinate, [values[f"c{i}"] for i in range(n)])
        u = sum(c/R**power for power, c in self.C.items())
        ue = float(sum(c/re**power for power, c in self.C.items()))
        beta = y*np.log(2*de/ue)+(1-y)*exponent
        shape = u/ue*np.exp(np.clip(-beta*yp, -100, 100))
        return values.get("Te", 0.)+de*(1-shape)**2

    def adiabatic(self, R, mass_ratio=1.):
        if mass_ratio != 1.:
            raise ValueError("research models are restricted to reference-isotope 127I2")
        return np.zeros_like(R)

    def nonadiabatic(self, R, mass_ratio=1.):
        R = np.asarray(R, float)
        values, cfg = self.values, self.config
        re = np.exp(values["logRe"])
        y = (R-re)/(R+re)
        coordinate = (y-.3)/.4 if cfg["mapped_alpha"] else y
        count = cfg["alpha_x" if self.state == "X" else "alpha_b"]
        if not count:
            return np.zeros_like(R)
        return mass_ratio*2*re/(R+re)*chebval(coordinate, [values[f"alpha{i}"] for i in range(count)])


class ResearchModel:
    """One immutable coefficient snapshot selected by canonical ID or role.

    ``frequency('P(53) 32-0', 'a1')`` returns MHz. ``transition`` and
    ``hyperfine_components`` follow RovibronicModel's argument order and units.
    Fitted radial functions include J±2 hyperfine mixing; no measured correction
    tables are applied. This API does not provide a line-intensity catalog or
    attach production uncertainty estimates to research predictions.
    """
    def __new__(cls, model, isotopologue="127I2"):
        if cls is ResearchModel and isinstance(model, str) and _entry(model).get("kind") == "expanded":
            cls = _ExpandedModel
        return super().__new__(cls)

    def __init__(self, model, isotopologue="127I2"):
        if isotopologue != "127I2":
            raise ValueError("research models support only 127I2")
        self._entry = _entry(model)
        self.info = _info(self._entry)
        self._spec = _read(self._entry, "model.json")
        self._validate()
        self.isotopologue = isotopologue
        self.parameters = self.info.id
        self._solvers, self._radial, self._structures, self._lines = {}, {}, {}, {}
        self._lock = RLock()

    def _validate(self):
        spec, cfg = self._spec, self._spec["config"]
        unsupported = [key for key in ("damping", "free_lr", "fit_reference", "positive_rotation", "correct_b_c10", "adiabatic_b") if cfg[key]]
        if cfg["family"] != "mlr" or cfg["family_x"] or cfg["family_b"] or unsupported:
            raise ValueError("unsupported change to the frozen research MLR configuration")
        inventory = self.info.inventory
        actual_free = sum(len(v) for v in spec["optimizer_parameters"].values())
        if spec["schema_version"] != 2 or spec["inventory"] != asdict(inventory) or actual_free != inventory.fitted_molecular:
            raise ValueError("frozen research coefficient inventory is inconsistent")
        if inventory.line_specific or inventory.level_specific or spec["isotopologue"] != "127I2":
            raise ValueError("frozen research model violates its shared-function domain")
        if spec["exclusion"] != self.info.exclusion or not self._entry["verification"]["passed"]:
            raise ValueError("research model lacks its exclusion or successful frozen verification")

    @property
    def coefficients(self):
        """A copy of all frozen coefficients/configuration, safe for inspection."""
        return json.loads(json.dumps(self._spec))

    def _moments(self, state, J):
        key = state, int(J)
        if key not in self._radial:
            if state not in self._solvers:
                self._solvers[state] = BSplineSolver(_Potential(state, self._spec), reduced_mass("127I2"),
                                                     **self._spec["grids"][state])
            solver = self._solvers[state]
            energy, psi = solver.wavefunctions(J)
            re = np.exp(self._spec["values"][state]["logRe"])
            z = (solver.R**2-re**2)/(solver.R**2+re**2)/self._spec["hyperfine"]["radial_scale"][state]
            mean = (solver.W[:, None]*psi**2).T@polyvander(z, 7)
            if np.max(abs(mean[:, 0]-1.)) > 1e-10:
                raise ValueError("research radial wavefunction normalization failed")
            self._radial[key] = energy, mean
        return self._radial[key]

    def _energy(self, state, v, J):
        return float(self._moments(state, J)[0][v])

    def _couplings(self, state, v, J):
        hf = self._spec["hyperfine"]
        if state == "X" and hf["fixed_x"]:
            u, x = v+.5, J*(J+1)
            q, c = hf["fixed_x_constants"]["eqQ"], hf["fixed_x_constants"]["C"]
            e = self._energy("X", v, 0)-self._energy("X", 0, 0)
            return HyperfineParameters(q[0]+q[1]*u+q[2]*u*u+q[3]*x+q[4]*u*x,
                c[0]+c[1]*u+c[2]/(e-c[3]), hf["fixed_x_constants"]["d"], hf["fixed_x_constants"]["delta"])
        moments = self._moments(state, J)[1][v]
        coefficients = hf["coefficients"][state]
        values = {p: float(np.dot(c, moments[:len(c)])) for p, c in coefficients.items()}
        if state == "X":
            if hf["anchor_x"]:
                anchor = self._moments("X", 0)[1][0]
                for p in ("eqQ", "C"):
                    c = coefficients[p]
                    values[p] += hf["x_ground_anchors"][p]-float(np.dot(c, anchor[:len(c)]))
            values.update(hf["x_spin_spin"])
        return HyperfineParameters(**values)

    def _structure(self, state, v, J):
        key = state, v, J
        if key not in self._structures:
            energy = self._energy(state, v, J)
            self._structures[key] = level_structure(J, lambda j: self._couplings(state, v, j),
                lambda j: (self._energy(state, v, j)-energy)*MHZ_PER_CM,
                symmetry="g" if state == "X" else "u", dJ=2)
        return self._structures[key]

    def _line(self, line):
        branch, J, vu, vl = line
        self.info.domain.validate(branch, J, vu, vl)
        key = branch, int(J), int(vu), int(vl)
        if key not in self._lines:
            upper = J+(1 if branch == "R" else -1)
            center = (self._energy("B", vu, upper)-self._energy("X", vl, J))*MHZ_PER_CM
            components = line_components(self._structure("B", vu, upper), self._structure("X", vl, J), upper, J)
            self._lines[key] = center, tuple(components)
        return self._lines[key]

    def hyperfine_components(self, v_upper, v_lower, J_lower, branch, dJ=2):
        """Return center MHz and all components, including weak unlabeled ones."""
        if dJ != 2:
            raise ValueError("frozen research predictions require dJ=2")
        with self._lock:
            center, components = self._line((branch, J_lower, v_upper, v_lower))
            return center, [replace(component) for component in components]

    def transition(self, v_upper, v_lower, J_lower, branch):
        """Hyperfine-free transition wavenumber in cm^-1."""
        with self._lock:
            self.info.domain.validate(branch, J_lower, v_upper, v_lower)
            upper = J_lower+(1 if branch == "R" else -1)
            return self._energy("B", v_upper, upper)-self._energy("X", v_lower, J_lower)

    def predict(self, line, component=None):
        """Predict a labeled line (or branch/J/v_upper/v_lower tuple), in MHz."""
        try:
            parsed = parse_label(line) if isinstance(line, str) else tuple(line)
        except TypeError as error:
            raise ValueError("line must be a label or (branch, J_lower, v_upper, v_lower) tuple") from error
        if len(parsed) != 4:
            raise ValueError("line tuple must be (branch, J_lower, v_upper, v_lower)")
        branch, J, vu, vl = parsed
        self.info.domain.validate(branch, J, vu, vl)
        canonical = None
        if component is not None:
            if not isinstance(component, str):
                raise ValueError("component must be a label such as 'a1'")
            canonical = f"a{component_rank(component)}"
        with self._lock:
            if canonical is None:
                value = self.transition(vu, vl, J, branch)*MHZ_PER_CM
            else:
                center, components = self._line(parsed)
                offsets = {c.label: c.offset for c in components if c.label}
                if canonical not in offsets:
                    raise ValueError(f"component {component!r} is unavailable for this line; choose from {', '.join(offsets)}")
                value = center+offsets[canonical]
        label = f"{branch}({J}) {vu}-{vl}"
        return ResearchPrediction(self.info.id, self.isotopologue, label, canonical,
            float(value), float(value/MHZ_PER_CM), f"127I2 {label}" == self.info.exclusion)

    def frequency(self, line, component=None):
        """Frequency in MHz; omitted component requests the hyperfine-free center."""
        return self.predict(line, component).frequency_MHz

    def verify(self, tolerance_Hz=100.):
        """Read-only replay of 60 absolute and 592 intraline frozen references.

        Checks numerical reproduction, not predictive accuracy. Coefficients,
        source assets and installed package files are never modified.
        """
        if not math.isfinite(tolerance_Hz) or tolerance_Hz <= 0:
            raise ValueError("verification tolerance must be finite and positive")
        rows = _read(self._entry, "reference_predictions.json")
        counts, maximum = dict(absolute=0, intraline=0), dict(absolute=0., intraline=0.)
        for row in rows:
            value = self.frequency(row["line"], row["component"])
            kind = "absolute"
            if "ref_component" in row:
                value -= self.frequency(row["line"], row["ref_component"])
                kind = "intraline"
            counts[kind] += 1
            maximum[kind] = max(maximum[kind], abs(value-row["predicted_MHz"])*1e6)
        result = dict(model=self.info.id, counts=counts, maximum_difference_Hz=maximum,
            tolerance_Hz=tolerance_Hz, passed=counts == self._spec["reference_counts"]
            and all(error <= tolerance_Hz for error in maximum.values()),
            scope="Numerical reproduction of frozen predictions; not independent predictive accuracy")
        if not result["passed"]:
            raise ValueError("frozen research prediction replay failed: "+json.dumps(result))
        return result


class _ExpandedModel(ResearchModel):
    """balanced-expanded (release 20261020): research-balanced's frozen coefficients with tails to X v<=108,
    B v<=80, balanced's hyperfine in its own domain and tabulated second-order hyperfine beyond
    (i2spec.research_expanded). Selected as ResearchModel("balanced_expanded")."""

    def __init__(self, model, isotopologue="127I2"):
        from . import research_expanded
        self._entry = _entry(model)
        # release 20261021 extends the frozen 20261020 model: its files are read, SHA256-checked, from there
        extends = self._entry.get("extends")
        self._core_entry = _entry(extends["id"]) if extends else self._entry
        if extends and any(self._core_entry["files_sha256"][k] != v for k, v in extends["files_sha256"].items()):
            raise ValueError(f"{self._entry['id']} extends a frozen model whose files have changed")
        supported = self._entry.get("isotopologues", [self._entry["isotopologue"]])
        if isotopologue not in supported:
            raise ValueError(f"{self._entry['id']} supports {', '.join(supported)} only")
        self.info = _info(self._entry, isotopologue)
        self._spec = _read(self._core_entry, "model.json")
        tables = _read(self._core_entry, "hyperfine.json")
        self._extension = _read(self._entry, "isotopologues.json") if extends else None
        base = self._spec["base"]
        if _entry(base["model"])["files_sha256"]["model.json"] != base["model_sha256"]:
            raise ValueError("balanced-expanded's base coefficients differ from the frozen balanced release")
        self._validate()
        self.isotopologue = isotopologue
        self.parameters = self.info.id
        self._solvers, self._radial, self._structures, self._lines = {}, {}, {}, {}
        self._lock = RLock()
        reference = _reference_model(self._entry["id"]) if isotopologue != "127I2" else None
        research_expanded.install(self, self._spec, tables, ResearchModel(base["model"]), self._extension,
                                  isotopologue, reference)
        self._unc = None

    def _uncertainty_model(self):
        if self._unc is None:
            from .research_expanded import ExpandedUncertainty
            table = _read(self._core_entry, "uncertainty.json")
            if table["model"] != self._core_entry["id"] or table["pending_outcomes_used"]:
                raise ValueError("balanced-expanded's uncertainty model does not belong to this frozen model")
            self._unc = ExpandedUncertainty(table, self._base.coefficients, self._spec["structure"],
                                            self._spec["values"])
        return self._unc

    def _uncertainty_terms(self, line, component):
        """(statistical factor, {level: (sign, model sigma)}, class, (upper, lower) status, hyperfine sigma and status)
        of one line or component."""
        parsed = parse_label(line) if isinstance(line, str) else tuple(line)
        branch, J, vu, vl = parsed
        self.info.domain.validate(branch, J, vu, vl)
        unc, core = self._uncertainty_model(), self._spec["core"]
        Ju = J+(1 if branch == "R" else -1)
        (vu_e, Ju_e), (vl_e, Jl_e) = self._equivalent("B", vu, Ju), self._equivalent("X", vl, J)
        if self._iso is None:
            a = unc.level_factor(self._levels_x, "B", vu, Ju)-unc.level_factor(self._levels_x, "X", vl, J)
        else:                       # conditioned on the isotopologue lines, with the isotope shift's terms
            shift = self._shift_unc()
            a = (shift.level_vector("B", unc.level_factor(self._levels_x, "B", vu, Ju), vu, Ju)
                 - shift.level_vector("X", unc.level_factor(self._levels_x, "X", vl, J), vl, J))
        cls = unc.line_class(core, (vu_e, Ju_e), (vl_e, Jl_e))
        levels = {("B", vu, Ju): (1., self._level_sigma(unc, core, "B", vu, Ju, cls)),
                  ("X", vl, J): (-1., self._level_sigma(unc, core, "X", vl, J, cls))}
        status = unc.status("B", vu_e, Ju_e), unc.status("X", vl_e, Jl_e)
        if self._iso is not None:
            order = ["measured", "interpolated", "extrapolated"]
            shifted = self._shift_unc().status(self._iso.map.v_ref("B", vu))
            status = max(status[0], shifted, key=order.index), status[1]
        hf, hf_status, name, offset = 0., None, None, 0.
        if component is not None:
            _, components = self._line((branch, J, vu, vl))
            offsets = {c.label: c.offset for c in components if c.label}
            name = f"a{component_rank(component)}"
            if name not in offsets:
                raise ValueError(f"component {component!r} is unavailable for this line; choose from {', '.join(offsets)}")
            a_hf, kappa = unc.hyperfine_terms(core, vu_e, vl_e)
            if self._iso is not None:
                a_hf = math.hypot(a_hf, self._extension["hyperfine"]["floor_MHz"])
            offset = offsets[name]
            hf = math.hypot(a_hf, kappa*offset)
            order = ["core", "validated", "interpolated", "resonance", "extrapolated"]
            hf_status = max((unc.hyperfine_status("B", vu_e), unc.hyperfine_status("X", vl_e)),
                            key=lambda text: order.index(text.split(":")[0].split(" ")[0]))
            if self._iso is not None and hf_status.startswith("core"):
                hf_status = "validated: " + self._extension["hyperfine"]["status"]
        return dict(label=f"{branch}({J}) {vu}-{vl}", a=a, levels=levels, cls=cls, status=status, hf=hf,
                    hf_status=hf_status, component=name, offset=offset, kappa_floor=(unc, core, vu_e, vl_e))

    def level_terms(self, state, v, J):
        """What the line uncertainties need of one level (preferred.py tabulates it): its statistical vector,
        data status, whether its 127I2 equivalent is within 15 in J of precision rows, its model sigma for each
        line class (precise, atlas, beyond) and whether it lies beyond the core."""
        unc, core = self._uncertainty_model(), self._spec["core"]
        a = unc.level_factor(self._levels_x, state, v, J)
        ve, Je = self._equivalent(state, v, J)
        status = unc.status(state, ve, Je)
        if self._iso is not None:
            shift = self._shift_unc()
            a = shift.level_vector(state, a, v, J)
            if state == "B":
                order = ["measured", "interpolated", "extrapolated"]
                status = max(status, shift.status(self._iso.map.v_ref("B", v)), key=order.index)
        sigma = [self._level_sigma(unc, core, state, v, J, cls) for cls in ("precise", "atlas", "beyond")]
        return a, status, unc._near("precise", state, ve, Je), sigma, ve > core[state]

    def hyperfine_terms(self, vu, Ju, vl, J):
        """(floor a MHz, scale error kappa) of a line's components: the frozen recipe at the 127I2 equivalents,
        with the isotopologue floor."""
        unc, core = self._uncertainty_model(), self._spec["core"]
        a, kappa = unc.hyperfine_terms(core, self._equivalent("B", vu, Ju)[0], self._equivalent("X", vl, J)[0])
        if self._iso is not None:
            a = math.hypot(a, self._extension["hyperfine"]["floor_MHz"])
        return a, kappa

    def _equivalent(self, state, v, J):
        """The 127I2 level the frozen uncertainty tables are read at: the level itself, or an isotopologue's
        equivalent (the same J = 0 energy, rho J(J+1))."""
        return (v, J) if self._iso is None else self._iso.map.equivalent(state, v, J)

    def _level_sigma(self, unc, core, state, v, J, cls):
        """One level's model sigma (MHz): the frozen recipe at the 127I2 equivalent; with release 20261021's
        extension, grown as J(J+1) beyond the J the data reach, and for an isotopologue's B level with the isotope
        shift's floor."""
        ve, Je = self._equivalent(state, v, J)
        sigma = unc.model_sigma(core, state, ve, Je, cls)
        if self._extension is None:
            return sigma
        reach = unc.reach_J(state, ve)
        if reach is not None and Je > reach:
            sigma *= Je*(Je+1)/(reach*(reach+1))
        if self._iso is not None and state == "B":
            sigma = math.hypot(sigma, self._shift_unc().floor)
        return sigma

    def _shift_unc(self):
        if getattr(self, "_shift_terms", None) is None:
            from .research_expanded import IsotopeShiftUncertainty
            self._shift_terms = IsotopeShiftUncertainty(self._extension, self._levels_x)
        return self._shift_terms

    def uncertainty(self, line, component=None, ref_component=None, ref_line=None):
        """1-sigma (PredictionUncertainty, MHz) of ``frequency(line, component)``; with ``ref_component``, of the
        difference ``frequency(line, component) - frequency(ref_line or line, ref_component)``.

        Statistical: the fit's linearized covariance through the levels; model: the error of predicting a level
        from the data around it, fitted to held-out residuals (a level the two lines share cancels); hyperfine: a
        floor and a pattern-scale error, fitted to held-out intervals (within one line the center cancels and the
        interval carries sqrt(2 a^2 + (kappa Delta)^2)). Frozen recipe and checks: release 20261020 uncertainty.json.
        """
        if (ref_component is not None or ref_line is not None) and component is None and ref_line is None:
            raise ValueError("an interval needs both components")
        with self._lock:
            one = self._uncertainty_terms(line, component)
            order = ["measured", "interpolated", "extrapolated"]
            statuses = list(one["status"])
            if ref_line is None and ref_component is not None:        # intraline: hyperfine only
                two = self._uncertainty_terms(line, ref_component)
                unc, core, vu, vl = one["kappa_floor"]
                a_hf, kappa = unc.hyperfine_terms(core, vu, vl)
                if self._iso is not None:
                    a_hf = math.hypot(a_hf, self._extension["hyperfine"]["floor_MHz"])
                stat = model = 0.
                hf = math.sqrt(2*a_hf**2+(kappa*(one["offset"]-two["offset"]))**2)
                ref_name = two["component"]
            elif ref_line is not None:                                 # one line against another
                two = self._uncertainty_terms(ref_line, ref_component)
                stat = float(np.linalg.norm(one["a"]-two["a"]))
                net = {}
                for terms, sign in ((one, 1.), (two, -1.)):
                    for key, (s, sigma) in terms["levels"].items():
                        c, _ = net.get(key, (0., sigma))
                        net[key] = (c+sign*s, sigma)
                model = math.sqrt(sum((c*sigma)**2 for c, sigma in net.values()))
                hf = math.hypot(one["hf"], two["hf"])
                statuses += list(two["status"])
                ref_name = two["component"]
            else:
                stat = float(np.linalg.norm(one["a"]))
                model = math.sqrt(sum(sigma**2 for _, sigma in one["levels"].values()))
                hf, ref_name = one["hf"], None
        label = one["label"] if ref_line is None else f"{one['label']} - {two['label']}"
        return PredictionUncertainty(
            self.info.id, label, one["component"], ref_name, math.sqrt(stat**2+model**2+hf**2), stat, model, hf,
            max(statuses, key=order.index), one["status"][0], one["status"][1], one["cls"], one["hf_status"])

    def _validate(self):
        spec, inventory = self._spec, self.info.inventory
        if spec["schema_version"] != 3 or spec["inventory"] != asdict(_info(self._core_entry).inventory):
            raise ValueError("frozen research coefficient inventory is inconsistent")
        if inventory.line_specific or inventory.level_specific or spec["isotopologue"] != "127I2":
            raise ValueError("frozen research model violates its shared-function domain")
        if spec["exclusion"] != self.info.exclusion or not self._entry["verification"]["passed"]:
            raise ValueError("research model lacks its exclusion or successful frozen verification")

    @property
    def coefficients(self):
        """A copy of the frozen tails and settings, the base model's coefficients and the hyperfine tables'
        SHA256, safe for inspection."""
        return json.loads(json.dumps(dict(self._spec, base_coefficients=self._base.coefficients,
                                          hyperfine_sha256=self._core_entry["files_sha256"]["hyperfine.json"],
                                          extension=self._extension)))

    def _moments(self, state, J):
        raise NotImplementedError("balanced-expanded keeps no radial moments: see research_expanded")

    def _energy(self, state, v, J):
        return self._levels_x.energy(state, v, J)

    def _couplings(self, state, v, J):
        if self._iso is not None:
            return self._iso.couplings(state, v, J)
        split = self._levels_x.j_split
        if v <= self._spec["core"][state]:
            if split is None or J <= split:
                return self._base._couplings(state, v, J)
            # above the fit's J the core fit's own grid does not hold the level: its functions on the model's
            from .research_expanded import core_couplings
            e = self._energy("X", v, 0)-self._energy("X", 0, 0) if state == "X" else None
            return core_couplings(self._base._spec, self._levels_x, state, v, J, u=v+.5, x=J*(J+1), e=e)
        return self._hf.constants(state, v, J, clamp=split is not None)

    def _structure(self, state, v, J):
        from .research_expanded import structure
        key = state, v, J
        if key not in self._structures:
            self._structures[key] = structure(self, state, v, J)
        return self._structures[key]

    def verify(self, tolerance_Hz=100.):
        """Read-only replay of the frozen centers, components and intraline intervals (see ResearchModel.verify)."""
        if not math.isfinite(tolerance_Hz) or tolerance_Hz <= 0:
            raise ValueError("verification tolerance must be finite and positive")
        rows = _read(self._entry, "reference_predictions.json")
        rows = [row for row in rows if row.get("isotopologue", "127I2") == self.isotopologue]
        expected = self._entry["reference_counts"][self.isotopologue] if self._extension else self._spec["reference_counts"]
        counts, maximum = dict(center=0, absolute=0, intraline=0), dict(center=0., absolute=0., intraline=0.)
        for row in rows:
            value = self.frequency(row["line"], row["component"])
            kind = "absolute" if row["component"] else "center"
            if row.get("ref_component"):
                value -= self.frequency(row["line"], row["ref_component"])
                kind = "intraline"
            counts[kind] += 1
            maximum[kind] = max(maximum[kind], abs(value-row["predicted_MHz"])*1e6)
        result = dict(model=self.info.id, counts=counts, maximum_difference_Hz=maximum,
            tolerance_Hz=tolerance_Hz, passed=counts == expected
            and all(error <= tolerance_Hz for error in maximum.values()),
            scope="Numerical reproduction of frozen predictions; not independent predictive accuracy")
        if not result["passed"]:
            raise ValueError("frozen research prediction replay failed: "+json.dumps(result))
        return result
