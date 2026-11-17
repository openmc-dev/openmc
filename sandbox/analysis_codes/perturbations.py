"""Representation-aware native-grid candidate generation, never Python sampling."""

from copy import deepcopy
import numpy as np
import openmc.lib as lib
from openmc.lib.perturbation import _areas
from .registry import Parameter


def equal(a, b) -> bool:
    """Compare snapshots exactly, including PDF/CDF arrays."""
    if isinstance(a, dict):
        return a.keys() == b.keys() and all(equal(a[k], b[k]) for k in a)
    if isinstance(a, (list, tuple)):
        return len(a) == len(b) and all(equal(x, y) for x, y in zip(a, b))
    return bool(np.array_equal(a, b))


def shape(law: dict, parameters: list[Parameter], values: np.ndarray) -> dict:
    """Preserve native grids/discrete mass and rebuild normalized PDF/CDF.

    Modes sharing a native target share a total fractional budget. Strict mode
    centers each basis with the native PDF integration rule and bounds it by
    one, making normalization independent of the coefficient. Normalized mode
    uses raw bounded basis functions, then normalizes; its final ratio can
    exceed 1±limit. Coefficients are independent uniforms in [-limit,+limit].
    """
    result = deepcopy(law)
    key = "energy" if law["type"] == "continuous-tabular" else "x"
    if len({p.shape for p in parameters}) != 1:
        raise ValueError("Do not mix shape semantics on the same native law")
    for incident, table in zip(result["incident_energy"], result["tables"]):
        active = [
            (p, v)
            for p, v in zip(parameters, values)
            if p.interval[0] <= incident <= p.interval[1]
        ]
        if not active or all(v == 0 for _, v in active):
            continue
        n = table.get("n_discrete", 0)
        x, pdf = table[key][n:], table["pdf"][n:].copy()
        if len(x) < 2:
            continue
        mass = _areas(x, pdf, table["interpolation"]).sum()
        target = 1 - table["pdf"][:n].sum()
        if mass <= 0 or target <= 0:
            raise ValueError("No mutable continuous probability mass")
        change = np.zeros_like(pdf)
        # The two declared basis modes always share the budget, including
        # when one is screened out. This preserves parameter identity.
        for p, value in active:
            z = x / max(float(np.max(x)), 1.0) if key == "energy" else x
            f = (
                z
                if p.mode == 1
                else (z * z if key == "energy" else (3 * z * z - 1) / 2)
            )
            if p.shape == "strict":
                f = f - _areas(x, pdf * f, table["interpolation"]).sum() / mass
            f = f / max(float(np.max(abs(f))), 1e-300)
            change += value * f / 2
        new = pdf * (1 + change)
        new *= target / _areas(x, new, table["interpolation"]).sum()
        if not np.all(np.isfinite(new)) or np.any(new < 0):
            raise ValueError("Invalid candidate PDF")
        if parameters[0].shape == "strict":
            positive = pdf > 0
            ratio = new[positive] / pdf[positive]
            limit = max(p.limit for p, _ in active)
            if np.any(abs(ratio - 1) > limit + 1e-6):
                raise ValueError(
                    "Native baseline normalization prevents strict bound"
                )
        table["pdf"][n:] = new
        table["cdf"][:n] = np.cumsum(table["pdf"][:n])
        table["cdf"][n:] = (
            table["pdf"][:n].sum()
            + np.r_[0.0, np.cumsum(_areas(x, new, table["interpolation"]))]
        )
    return result


class NativeState:
    """Owned Python baseline snapshots; native context owns rollback copies."""

    def __init__(self, parameters: list[Parameter]):
        self.parameters = parameters
        self.baseline = {}
        self.grids = {}
        for p in parameters:
            if p.target in self.baseline:
                continue
            if p.family == "xs":
                i = lib.nuclides[p.nuclide]._index
                grid = lib.nuclide_energy_grid(i)
                xs = lib.get_nuclide_xs(i, p.mt)
                full = np.zeros_like(grid)
                full[len(grid) - len(xs) :] = xs
                self.baseline[p.target] = full
                self.grids[p.target] = grid
            else:
                self.baseline[p.target] = self.read(p)

    def read(self, p: Parameter):
        if p.family == "xs":
            i = lib.nuclides[p.nuclide]._index
            xs = lib.get_nuclide_xs(i, p.mt)
            full = np.zeros_like(self.grids[p.target])
            full[len(full) - len(xs) :] = xs
            return full
        return lib.get_emission_distribution(
            p.nuclide, p.mt, p.product, p.distribution
        )[p.family]

    def verify(self) -> None:
        for p in self.parameters:
            if not equal(self.read(p), self.baseline[p.target]):
                raise RuntimeError(f"Native baseline changed: {p.identity}")

    def candidates(self, values: np.ndarray) -> dict:
        values = np.asarray(values, float)
        if values.shape != (len(self.parameters),) or not np.all(
            np.isfinite(values)
        ):
            raise ValueError("Invalid complete parameter vector")
        grouped = {}
        for p, v in zip(self.parameters, values):
            if abs(v) > p.limit + 1e-14:
                raise ValueError("Parameter outside declared prior")
            grouped.setdefault(p.target, []).append((p, v))
        result = dict(xs=[], energy=[], angle=[])
        for target, entries in grouped.items():
            if all(v == 0 for p, v in entries):
                continue
            p = entries[0][0]
            original = self.baseline[target]
            if p.family == "xs":
                from .cross_sections import perturb

                array = perturb(
                    original,
                    self.grids[target],
                    [e[0] for e in entries],
                    [e[1] for e in entries],
                )
                result["xs"].append(
                    dict(nuclide=p.nuclide, mt=p.mt, values=array)
                )
            else:
                result[p.family].append(
                    dict(
                        nuclide=p.nuclide,
                        mt=p.mt,
                        product=p.product,
                        distribution=p.distribution,
                        law=shape(
                            original,
                            [e[0] for e in entries],
                            np.array([e[1] for e in entries]),
                        ),
                    )
                )
        return result

    def verify_candidates(self, candidates: dict) -> None:
        for family, entries in candidates.items():
            for item in entries:
                p = next(
                    p
                    for p in self.parameters
                    if p.family == family
                    and p.nuclide == item["nuclide"]
                    and p.mt == item["mt"]
                    and p.product == item.get("product", 0)
                    and p.distribution == item.get("distribution", 0)
                )
                actual = self.read(p)
                if family == "xs":
                    np.testing.assert_array_equal(actual, item["values"])
                else:
                    for a, b in zip(actual["tables"], item["law"]["tables"]):
                        np.testing.assert_allclose(
                            a["pdf"], b["pdf"], rtol=1e-12, atol=1e-14
                        )
                        np.testing.assert_allclose(
                            a["cdf"], b["cdf"], rtol=1e-10, atol=1e-12
                        )
