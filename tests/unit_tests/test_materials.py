from pathlib import Path

import openmc
from openmc.deplete import Chain


def test_materials_deplete():
    pristine_material_1 = openmc.Material()
    pristine_material_1.add_nuclide("Ni58", 1.)
    pristine_material_1.set_density("g/cm3", 7.87)
    pristine_material_1.depletable = True
    pristine_material_1.temperature = 293.6
    pristine_material_1.volume = 1.

    pristine_material_2 = openmc.Material()
    pristine_material_2.add_nuclide("Ni60", 1.)
    pristine_material_2.set_density("g/cm3", 7.87)
    pristine_material_2.depletable = True
    pristine_material_2.temperature = 293.6
    pristine_material_2.volume = 1.

    pristine_materials = openmc.Materials([pristine_material_1, pristine_material_2])

    mg_flux = [0.5e11] * 42

    chain = Chain.from_xml(
        Path(__file__).parents[1] / "chain_ni.xml"
    )

    depleted_material = pristine_materials.deplete(
        multigroup_fluxes=[mg_flux, mg_flux],
        energy_group_structures=["VITAMIN-J-42", "VITAMIN-J-42"],
        timesteps=[100, 100],
        source_rates=[1e19, 0.0],
        timestep_units="d",
        chain_file=chain,
    )

    assert list(depleted_material.keys()) == [pristine_material_1.id, pristine_material_2.id]
    for mat_id, materials in depleted_material.items():
        for material in materials:
            assert isinstance(material, openmc.Material)
            assert len(material.get_nuclides()) > 1
            assert mat_id == material.id

    mats = depleted_material[pristine_material_1.id]
    Co58_mat_1_step_0 = mats[0].get_nuclide_atom_densities("Co58")["Co58"]
    Co58_mat_1_step_1 = mats[1].get_nuclide_atom_densities("Co58")["Co58"]
    Co58_mat_1_step_2 = mats[2].get_nuclide_atom_densities("Co58")["Co58"]

    assert Co58_mat_1_step_0 == 0.0
    # Co58 is the main activation product of Ni58 in the first irradiation step.
    # It then decays in the second cooling step (flux = 0)
    assert Co58_mat_1_step_1 > 0.0 and Co58_mat_1_step_1 > Co58_mat_1_step_2

    Ni59_mat_1_step_0 = mats[0].get_nuclide_atom_densities("Ni59")["Ni59"]
    Ni59_mat_1_step_1 = mats[1].get_nuclide_atom_densities("Ni59")["Ni59"]
    Ni59_mat_1_step_2 = mats[2].get_nuclide_atom_densities("Ni59")["Ni59"]

    assert Ni59_mat_1_step_0 == 0.0
    # Ni59 is one of the main activation product of Ni60 in the first irradiation
    # step. It then decays in the second cooling step (flux = 0)
    assert Ni59_mat_1_step_1 > 0.0 and Ni59_mat_1_step_1 > Ni59_mat_1_step_2
