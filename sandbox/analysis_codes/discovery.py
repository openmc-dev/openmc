"""Discover loaded native objects; report unsupported targets explicitly."""

from ctypes import byref, c_uint64
from dataclasses import asdict
import numpy as np
import openmc.lib as lib
from .registry import Parameter


def discover(
    nuclides: tuple[str, ...],
    intervals=((0.0, 20e6),),
    limit=0.2,
    mts=None,
    families=("xs", "energy", "angle"),
) -> tuple[list[Parameter], list[dict]]:
    """Inspect all loaded reactions, probing native backup eligibility.

    No mutation is performed. Backup handles are immediately freed. Secondary
    candidates are neutron products with existing independent tabular laws;
    other products and conditional/analytic representations remain reported.
    """
    parameters, report = [], []
    for name in nuclides:
        index = lib.nuclides[name]._index
        grid = lib.nuclide_energy_grid(index)
        try:
            lib.nuclide_energy_grid(index, 1)
            multiple_temperatures = True
        except RuntimeError:
            multiple_temperatures = False
        for mt in lib.get_nuclide_mt_numbers(index):
            mt = int(mt)
            if mts is not None and mt not in mts:
                continue
            info = lib.get_reaction_info(name, mt)
            row = dict(nuclide=name, family="xs", **info)
            try:
                if multiple_temperatures:
                    raise ValueError(
                        "Framework XS adapter requires a single loaded temperature; transactional multi-temperature parameterization is not implemented"
                    )
                handle = c_uint64()
                lib._dll.openmc_nuclear_data_backup(
                    index, mt, 0, 0, 0, byref(handle)
                )
                lib._dll.openmc_nuclear_data_backup_free(handle.value)
                xs = lib.get_nuclide_xs(index, mt)
                threshold = lib.get_reaction_threshold_energy(index, mt)
                row.update(
                    supported=True,
                    representation="native pointwise XS",
                    threshold_eV=threshold,
                )
                if "xs" in families:
                    for interval in intervals:
                        mask = (grid >= max(interval[0], threshold)) & (
                            grid <= interval[1]
                        )
                        if np.any(mask) and np.any(xs > 0):
                            parameters.append(
                                Parameter(
                                    name,
                                    mt,
                                    "xs",
                                    interval=tuple(interval),
                                    limit=limit,
                                )
                            )
            except Exception as error:
                row.update(supported=False, reason=str(error))
            report.append(row)
        for product in lib.get_reaction_products(name):
            mt = product["mt"]
            if mts is not None and mt not in mts:
                continue
            for descriptor in product["distributions"]:
                pi, di = (
                    product["product_index"],
                    descriptor["distribution_index"],
                )
                law = lib.get_emission_distribution(name, mt, pi, di)
                for family in ("energy", "angle"):
                    component = law.get(family)
                    representation = (component or {}).get("type", law["type"])
                    supported = (
                        product["particle"] == "neutron"
                        and law["type"] == "uncorrelated"
                        and representation
                        == (
                            "continuous-tabular"
                            if family == "energy"
                            else "tabular-angle"
                        )
                    )
                    row = dict(
                        nuclide=name,
                        mt=mt,
                        family=family,
                        product=pi,
                        distribution=di,
                        particle=product["particle"],
                        reference_frame=product["reference_frame"],
                        representation=representation,
                        supported=supported,
                        reason=(
                            ""
                            if supported
                            else "No independent transactional tabular adapter; no representation conversion"
                        ),
                    )
                    if supported:
                        incident = component["incident_energy"]
                        key = "energy" if family == "energy" else "x"
                        row.update(
                            incident_range_eV=[
                                float(incident[0]),
                                float(incident[-1]),
                            ],
                            outgoing_range=[
                                float(
                                    min(
                                        t[key].min()
                                        for t in component["tables"]
                                    )
                                ),
                                float(
                                    max(
                                        t[key].max()
                                        for t in component["tables"]
                                    )
                                ),
                            ],
                        )
                        if family in families:
                            for interval in intervals:
                                if np.any(
                                    (incident >= interval[0])
                                    & (incident <= interval[1])
                                ):
                                    for mode in (1, 2):
                                        parameters.append(
                                            Parameter(
                                                name,
                                                mt,
                                                family,
                                                mode,
                                                pi,
                                                di,
                                                tuple(interval),
                                                limit,
                                            )
                                        )
                    report.append(row)
    unique = {p.identity: p for p in parameters}
    return list(unique.values()), report
