from pathlib import Path
from typing import Sequence, TypedDict

import numpy as np

import openmc
import openmc.checkvalue as cv
import openmc.deplete
import openmc.lib
from openmc.deplete import GROUP_STRUCTURES, REACTION_MT, Chain
from openmc.deplete.microxs import (_find_cross_sections,
                                    _get_nuclides_with_data,
                                    _resolve_chain_file_path)
from openmc.mgxs import GROUP_STRUCTURES
from openmc.utility_funcs import change_directory


class ActivationDict(TypedDict):
    """Contains necessary data for material activation calculations.

    This TypedDict defines the required parameters for performing activation
    calculations on materials using the deplete_materials function.
    """

    material: openmc.Material
    """Material to be depleted during activation calculations"""

    multigroup_flux: Sequence[float]
    """Energy-dependent multigroup flux values, length must match energy groups"""

    energy_group_structure: Sequence[float] | str
    """Energy group boundaries in [eV] or name of predefined group structure"""

    source_rate: Sequence[float]
    """Source rates for each timestep in the depletion calculation"""

    timesteps: Sequence[float]
    """Time intervals for the depletion steps in timestep_units"""


def deplete_materials(
    activation_data: Sequence[ActivationDict],
    timestep_units: str = "s",
    chain_file: cv.PathLike | None = None,
    nuclides: Sequence[str] | None = None,
    reactions: Sequence[str] | None = None,
    **init_kwargs: dict,
) -> Sequence[openmc.Materials]:
    """Deplete materials using multigroup flux and microscopic cross sections.

    This function computes the microscopic cross sections for the specified
    materials based on the provided multigroup flux and energy group boundaries.
    It then sets up a depletion simulation using OpenMC's deplete module and
    depletes the materials over the specified timesteps.

    .. versionadded:: 0.15.3

    Parameters
    ----------
    activation_data : Sequence[ActivationDict]
        Sequence of dictionaries containing the necessary data for each material
        to be depleted. Each dictionary should contain:
        - 'material': openmc.Material
            Material to be depleted, should have a temperature and volume set.
        - 'multigroup_flux': Sequence[float]
            Energy-dependent multigroup flux values, where each sublist corresponds
            to a specific material. Will be normalized so that it sums to 1.
        - 'energy_group_structure': Sequence[float] | str
            Energy group boundaries in [eV] or the name of the group structure.
        - 'source_rate': Sequence[float]
            Source rates for each timestep.
        - 'timesteps': Sequence[float]
            Timesteps for the depletion simulation.
    timestep_units : {'s', 'min', 'h', 'd', 'a', 'MWd/kg'}
        Units for values specified in the `timesteps` argument. 's' means
        seconds, 'min' means minutes, 'h' means hours, 'a' means Julian years
        and 'MWd/kg' indicates that the values are given in burnup (MW-d of
        energy deposited per kilogram of initial heavy metal).
    chain_file : str, optional
        Path to the depletion chain XML file that will be used in depletion
        simulation. Defaults to ``openmc.config['chain_file']``.
    nuclides : list of str, optional
        Nuclides to get cross sections for. If not specified, all burnable
        nuclides from the depletion chain file are used.
    reactions : list of str, optional
        Reactions to get cross sections for. If not specified, all neutron
        reactions listed in the depletion chain file are used.
    **init_kwargs : dict
        Keyword arguments passed to :func:`openmc.lib.init` defaults to
        {'output': False}

    Returns
    -------
    Sequence[openmc.Materials]

    """

    # default to not print terminal output
    if init_kwargs == {}:
        init_kwargs = {"output": False}

    # Validate activation_data entries
    for i, entry in enumerate(activation_data):
        # Check that entry is a dict with all required keys
        if not isinstance(entry, dict):
            raise TypeError(f"Entry {i} in activation_data is not a dictionary")

        # Check for required keys (using the ActivationDict definition)
        required_keys = {
            "material",
            "multigroup_flux",
            "energy_group_structure",
            "source_rate",
            "timesteps",
        }
        if not all(key in entry for key in required_keys):
            missing = required_keys - set(entry.keys())
            raise KeyError(f"Entry {i} is missing required keys: {missing}")

        # Check material field
        if not isinstance(entry["material"], openmc.Material):
            raise TypeError(f"Entry {i}: 'material' must be an openmc.Material object")

        # Check that every value in source_rate is a float
        for key in ["source_rate", "timesteps"]:
            if not all(isinstance(rate, (int, float)) for rate in entry[key]):
                raise TypeError(
                    f"Entry {i}: All values in 'source_rate' must be numeric (int or float)"
                )

        # Check temperature and volume are set
        if entry["material"].temperature is None:
            raise ValueError(
                f"Entry {i}: Material temperature must be set before depletion"
            )
        if entry["material"].volume is None:
            raise ValueError(f"Entry {i}: Material volume must be set before depletion")

    chain_file_path = _resolve_chain_file_path(Path(chain_file)).resolve()
    chain = Chain.from_xml(chain_file_path)

    cross_sections = _find_cross_sections(model=None)
    nuclides_with_data = _get_nuclides_with_data(cross_sections)

    # If no nuclides were specified, default to all nuclides from the chain
    if not nuclides:
        nuclides = chain.nuclides
        nuclides = [nuc.name for nuc in nuclides]

    # Get reaction MT values. If no reactions specified, default to the
    # reactions available in the chain file
    if reactions is None:
        reactions = chain.reactions
    mts = [REACTION_MT[name] for name in reactions]

    # Create 3D array for microscopic cross sections
    microxs_arr = np.zeros((len(nuclides), len(mts), 1))

    # Create a material with all nuclides
    mat_all_nucs = openmc.Material()
    for nuc in nuclides:
        if nuc in nuclides_with_data:
            mat_all_nucs.add_nuclide(nuc, 1.0)
    mat_all_nucs.set_density("atom/b-cm", 1.0)

    # Create simple model containing the above material
    surf1 = openmc.Sphere(boundary_type="vacuum")
    surf1_cell = openmc.Cell(fill=mat_all_nucs, region=-surf1)
    model = openmc.Model()
    model.geometry = openmc.Geometry([surf1_cell])
    model.settings = openmc.Settings(particles=1, batches=1, output={"summary": False})

    with change_directory(tmpdir=True):
        # Export model within temporary directory
        model.export_to_model_xml()

        with openmc.lib.run_in_memory(**init_kwargs):
            # For each material, energy and multigroup flux compute the flux-averaged cross section for the nuclides and reactions

            all_depleted_materials = []
            for entry in activation_data:
                material = entry["material"]
                multigroup_flux = entry["multigroup_flux"]
                energy = entry["energy_group_structure"]
                source_rates = entry["source_rate"]
                timesteps = entry["timesteps"]

                # Normalize multigroup flux
                multigroup_flux = np.array(multigroup_flux)
                multigroup_flux /= multigroup_flux.sum()

                # check_type("temperature", temperature, (int, float))
                # if energy is string then use group structure of that name
                if isinstance(energy, str):
                    energy = GROUP_STRUCTURES[energy]
                else:
                    # if user inputs energy check they are ascending (low to high) as
                    # some depletion codes use high energy to low energy.
                    if not np.all(np.diff(energy) > 0):
                        raise ValueError(
                            "Energy group boundaries must be in ascending order"
                        )

                # check dimension consistency
                if len(multigroup_flux) != len(energy) - 1:
                    msg = (
                        "Length of flux array should be len(energy)-1, but "
                        f"got {len(multigroup_flux)} multigroup_flux entries "
                        f"and {len(energy)} energy group boundaries"
                    )
                    raise ValueError(msg)

                for nuc_index, nuc in enumerate(nuclides):
                    if nuc not in nuclides_with_data:
                        continue
                    lib_nuc = openmc.lib.nuclides[nuc]
                    for mt_index, mt in enumerate(mts):
                        xs = lib_nuc.collapse_rate(
                            mt, material.temperature, energy, multigroup_flux
                        )
                        microxs_arr[nuc_index, mt_index, 0] = xs

                micro_xs = openmc.deplete.MicroXS(microxs_arr, nuclides, reactions)

                operator = openmc.deplete.IndependentOperator(
                    materials=openmc.Materials([material]),
                    fluxes=[material.volume],
                    micros=[micro_xs],
                    normalization_mode="source-rate",
                    reduce_chain_level=5,
                    chain_file=chain_file_path,
                )

                integrator = openmc.deplete.PredictorIntegrator(
                    operator=operator,
                    timesteps=timesteps,
                    source_rates=source_rates,
                    timestep_units=timestep_units,
                )

                integrator.integrate(path="depletion_results.h5")

                results = openmc.deplete.ResultsList.from_hdf5(
                    filename="depletion_results.h5"
                )
                depleted_material = openmc.Materials()
                for result in results:
                    mat = result.get_material(str(material.id))
                    depleted_material.append(mat)
                all_depleted_materials.append(depleted_material)

    return all_depleted_materials
