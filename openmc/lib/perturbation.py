"""Constrained runtime sensitivity perturbations of native nuclear data.

These shape coefficients are user-selected sensitivity parameters, not an
independent sampling engine or evaluated covariance uncertainties.
"""

from contextlib import contextmanager
from copy import deepcopy
from ctypes import POINTER, byref, c_bool, c_double, c_int, c_uint64

import numpy as np

from . import _dll
from .error import _error_handler
from .secondary import _index, _nuclide_index

__all__ = ["perturb_xs", "perturb_tabular"]


def perturb_xs(values, coefficients=(), basis=None, *, scale=1.0):
    """Return nonnegative XS times scale*exp(coefficients @ basis).

    Parameters
    ----------
    values : array-like
        Baseline microscopic cross sections; zeros remain zero.
    coefficients : array-like, optional
        Log-shape coefficients.
    basis : array-like, optional
        Basis evaluated on the native grid, shape (n_coefficients, n_values).
    scale : float, optional
        Positive constant multiplicative factor.

    Returns
    -------
    numpy.ndarray
        Independent perturbed array.
    """
    values = np.asarray(values, dtype=float)
    coefficients = np.asarray(coefficients, dtype=float)
    if (
        values.ndim != 1
        or not np.all(np.isfinite(values))
        or np.any(values < 0)
        or not np.isfinite(scale)
        or scale <= 0
    ):
        raise ValueError("Invalid baseline cross sections or scale")
    if coefficients.ndim != 1 or not np.all(np.isfinite(coefficients)):
        raise ValueError("Invalid coefficients")
    exponent = np.zeros_like(values)
    if coefficients.size:
        basis = np.asarray(basis, dtype=float)
        if basis.shape != (len(coefficients), len(values)) or not np.all(
            np.isfinite(basis)
        ):
            raise ValueError(
                "Basis must be finite with shape (coefficients, grid)"
            )
        exponent = coefficients @ basis
    with np.errstate(over="raise", invalid="raise"):
        result = values * scale * np.exp(exponent)
    if not np.all(np.isfinite(result)):
        raise ValueError("Nonfinite perturbed cross sections")
    return result


def _areas(x, p, interpolation):
    if interpolation == 1:
        return np.diff(x) * p[:-1]
    if interpolation == 2:
        return np.diff(x) * (p[:-1] + p[1:]) / 2
    raise ValueError(
        "Only native histogram/linear-linear interpolation supported"
    )


def perturb_tabular(law, coefficients, basis):
    """Exponentially tilt continuous PDFs while retaining native grids.

    Parameters
    ----------
    law : dict
        Continuous-tabular energy or tabular-angle runtime snapshot.
    coefficients : array-like
        Shape coefficients.
    basis : callable
        Called as basis(x, incident_energy), returning shape (n_coefficients,
        n_points). x is outgoing energy in eV or scattering cosine.

    Returns
    -------
    dict
        Candidate normalized snapshot with rebuilt CDFs. Discrete ENERGY line
        probabilities are kept unchanged; angular tables contain no atoms.

    Notes
    -----
    Exponentials are evaluated at native grid points. The resulting native
    interpolant is normalized using its original interpolation rule; this is
    not a continuously evaluated exponential between grid points.
    """
    if law["type"] not in ("continuous-tabular", "tabular-angle"):
        raise ValueError("Only existing independent tabular laws supported")
    result = deepcopy(law)
    coefficients = np.asarray(coefficients, dtype=float)
    if coefficients.ndim != 1 or not np.all(np.isfinite(coefficients)):
        raise ValueError("Invalid coefficients")
    key = "energy" if law["type"] == "continuous-tabular" else "x"
    for incident, table in zip(result["incident_energy"], result["tables"]):
        x, p = table[key], table["pdf"]
        n = table.get("n_discrete", 0)
        functions = np.asarray(basis(x, incident), dtype=float)
        if functions.shape != (len(coefficients), len(x)) or not np.all(
            np.isfinite(functions)
        ):
            raise ValueError("Invalid shape basis")
        discrete_mass = float(np.sum(p[:n]))
        if (
            discrete_mass < 0
            or discrete_mass > 1
            or not np.all(np.isfinite(p))
            or np.any(p < 0)
        ):
            raise ValueError("Invalid native probabilities")
        if n == len(x):
            continue
        exponent = coefficients @ functions
        exponent -= np.max(
            exponent[n:]
        )  # stable normalization; no shape change
        p[n:] *= np.exp(exponent[n:])
        areas = _areas(x[n:], p[n:], table["interpolation"])
        mass = np.sum(areas)
        if not np.isfinite(mass) or mass <= 0:
            raise ValueError("Invalid continuous mass")
        p[n:] *= (1.0 - discrete_mass) / mass
        cdf = np.empty_like(p)
        cdf[:n] = np.cumsum(p[:n])
        cdf[n:] = (
            discrete_mass
            + np.r_[
                0.0, np.cumsum(_areas(x[n:], p[n:], table["interpolation"]))
            ]
        )
        table["cdf"] = cdf
    return result


__all__ += ["nuclear_data_perturbation", "get_reaction_info"]

_dll.openmc_nuclear_data_backup.argtypes = [c_int] * 5 + [POINTER(c_uint64)]
_dll.openmc_nuclear_data_backup.restype = c_int
_dll.openmc_nuclear_data_backup.errcheck = _error_handler
_dll.openmc_nuclear_data_validate_xs.argtypes = [c_int]
_dll.openmc_nuclear_data_validate_xs.restype = c_int
_dll.openmc_nuclear_data_validate_xs.errcheck = _error_handler
_dll.openmc_nuclear_data_restore.argtypes = [c_uint64]
_dll.openmc_nuclear_data_restore.restype = c_int
_dll.openmc_nuclear_data_restore.errcheck = _error_handler
_dll.openmc_nuclear_data_backup_free.argtypes = [c_uint64]
_dll.openmc_nuclear_data_backup_free.restype = None


def _validate_fixed_tables(candidate, original, key):
    if (
        candidate["type"] != original["type"]
        or not np.array_equal(
            candidate["incident_energy"], original["incident_energy"]
        )
        or len(candidate["tables"]) != len(original["tables"])
    ):
        raise ValueError(
            "Incident grid and native representation must be preserved"
        )
    for c, o in zip(candidate["tables"], original["tables"]):
        if (
            not np.array_equal(c[key], o[key])
            or c["interpolation"] != o["interpolation"]
            or c.get("n_discrete", 0) != o.get("n_discrete", 0)
        ):
            raise ValueError(
                "Native coordinates/interpolation/discrete counts must be preserved"
            )
        p = np.asarray(c["pdf"], dtype=float)
        n = c.get("n_discrete", 0)
        if (
            p.shape != o["pdf"].shape
            or not np.all(np.isfinite(p))
            or np.any(p < 0)
        ):
            raise ValueError("Invalid probabilities")
        if n and not np.array_equal(p[:n], o["pdf"][:n]):
            raise ValueError("Discrete line probabilities must be preserved")
        mass = np.sum(p[:n]) + np.sum(
            _areas(c[key][n:], p[n:], c["interpolation"])
        )
        if not np.isclose(mass, 1.0, atol=1.0e-6, rtol=0.0):
            raise ValueError("Candidate PDF is not normalized")


@contextmanager
def nuclear_data_perturbation(*, xs=(), energy=(), angle=()):
    """Apply validated changes and restore exact native data on exit.

    Parameters
    ----------
    xs : sequence of dict, optional
        Entries: nuclide, mt, values (full native grid), temperature (default 0).
        Native threshold is preserved; zero below it is required. Redundant,
        URR, multipole, resonance and fission nuclides are rejected by backup.
    energy, angle : sequence of dict, optional
        Entries: nuclide, mt, law; product/distribution default to zero. Only
        existing independent tabular laws on unchanged grids are supported.

    Notes
    -----
    Uses owned native deep copies, not setter-based restoration. Every candidate
    is checked and every backup acquired before applying changes. A native
    setter/rebuild failure rolls back all changes before propagating the error.
    This is rollback atomicity under serialized library access, not concurrent
    transaction isolation. Do not finalize the library inside the context.
    Reads/tallies should be extracted inside the context after run(), before
    another simulation initialization resets tally accumulation.
    """
    from . import core, secondary

    entries, seen = [], set()
    for kind, items in enumerate((xs, energy, angle)):
        for item in items:
            item = deepcopy(item)
            index = _nuclide_index(item["nuclide"])
            mt = _index(item["mt"], "mt")
            pi = _index(item.get("product", 0), "product")
            di = _index(item.get("distribution", 0), "distribution")
            identity = (kind, index, mt, pi, di)
            if identity in seen:
                raise ValueError("Duplicate mutation target")
            seen.add(identity)
            if kind == 0:
                temp = _index(item.get("temperature", 0), "temperature")
                grid = core.nuclide_energy_grid(index, temp)
                threshold = core.get_reaction_threshold_energy(index, mt, temp)
                values = np.ascontiguousarray(item["values"], dtype=float)
                if (
                    values.shape != grid.shape
                    or not np.all(np.isfinite(values))
                    or np.any(values < 0)
                    or np.any(values[grid < threshold] != 0)
                ):
                    raise ValueError(
                        "Invalid XS values or nonzero values below threshold"
                    )
                item.update(values=values, grid=grid, threshold=threshold)
            else:
                original = secondary.get_emission_distribution(
                    index, mt, pi, di
                )
                if original["type"] != "uncorrelated":
                    raise ValueError("Conditional laws remain read-only")
                component = original["energy" if kind == 1 else "angle"]
                if component is None or component["type"] not in (
                    "continuous-tabular",
                    "tabular-angle",
                ):
                    raise ValueError(
                        "Only existing independent tabular components are mutable"
                    )
                _validate_fixed_tables(
                    item["law"], component, "energy" if kind == 1 else "x"
                )
            entries.append((kind, index, mt, pi, di, item))
    backups = []
    try:
        for kind, index, mt, pi, di, item in entries:
            handle = c_uint64()
            _dll.openmc_nuclear_data_backup(
                index, mt, pi, di, kind, byref(handle)
            )
            backups.append(handle.value)
        core.simulation_init()
        for kind, index, mt, pi, di, item in entries:
            if kind == 0:
                core.new_nuclide_xs_with_threshold(
                    index,
                    mt,
                    item["grid"],
                    item["values"],
                    item["threshold"],
                    item.get("temperature", 0),
                )
            elif kind == 1:
                core.set_emission_spectrum(index, mt, item["law"], pi, di)
            else:
                secondary.set_angular_distribution(
                    index, mt, item["law"], pi, di
                )
        for index in sorted({i for k, i, *_ in entries if k == 0}):
            core.nuclide_rebuild_derived_xs(index)
            _dll.openmc_nuclear_data_validate_xs(index)
        yield
    finally:
        # Every backup is consumed by restore. Continue restoring other targets
        # if one fails, then report the first failure instead of hiding it.
        errors = []
        for handle in reversed(backups):
            try:
                _dll.openmc_nuclear_data_restore(handle)
            except Exception as error:
                errors.append(error)
            finally:
                _dll.openmc_nuclear_data_backup_free(handle)
        if errors:
            raise RuntimeError("Native restoration failed") from errors[0]


_dll.openmc_nuclear_data_reaction_info.argtypes = [
    c_int,
    c_int,
    POINTER(c_double),
    POINTER(c_bool),
    POINTER(c_bool),
]
_dll.openmc_nuclear_data_reaction_info.restype = c_int
_dll.openmc_nuclear_data_reaction_info.errcheck = _error_handler


def get_reaction_info(nuclide, mt):
    """Inspect native reaction metadata even when no secondary product is stored.

    Parameters
    ----------
    nuclide : str or int
        Loaded nuclide name or index.
    mt : int
        Reaction MT.

    Returns
    -------
    dict
        MT, reaction name, native Q value in eV, redundant flag and frame.
    """
    from openmc.data import REACTION_NAME

    q, redundant, cm = c_double(), c_bool(), c_bool()
    mt = _index(mt, "mt")
    _dll.openmc_nuclear_data_reaction_info(
        _nuclide_index(nuclide), mt, byref(q), byref(redundant), byref(cm)
    )
    return dict(
        mt=mt,
        name=REACTION_NAME.get(mt, str(mt)),
        q_value=q.value,
        redundant=redundant.value,
        reference_frame="center-of-mass" if cm.value else "laboratory",
    )
