import openmc
from pathlib import Path


def test_materials_deplete():

    from openmc.deplete import Chain

    pristine_material_1 = openmc.Material(material_id=1)
    pristine_material_1.add_nuclide("Ni58", 1.)
    pristine_material_1.set_density("g/cm3", 7.87)
    pristine_material_1.depletable = True
    pristine_material_1.temperature = 293.6
    pristine_material_1.volume = 1.

    pristine_material_2 = openmc.Material(material_id=1)
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
        timesteps=[10,10],
        source_rates=[1e19, 0.0],
        timestep_units='s',
        chain_file= chain,
    )

    for mat_id, materials in depleted_material.items():
        for material in materials:
            assert isinstance(material, openmc.Material)
            assert len(material.get_nuclides()) > 1
            assert mat_id in [1, 2]
