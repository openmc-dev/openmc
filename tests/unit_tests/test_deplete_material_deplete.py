from pathlib import Path

import numpy as np
import pytest

import openmc
from openmc.deplete import deplete_materials
from openmc.mgxs import GROUP_STRUCTURES


def test_deplete_materials():

    # Ensure the OpenMC data directory is set
    openmc.data.data_directory = "data"
    mat_ni58 = openmc.Material()
    mat_ni58.add_nuclide("Ni58", 0.95)
    mat_ni58.set_density("g/cm3", 7.87)
    mat_ni58.depletable = True
    mat_ni58.volume = 1.0
    mat_ni58.temperature = 293.6

    mat_ni60 = openmc.Material()
    mat_ni60.add_nuclide("Ni60", 1.0)
    mat_ni60.set_density("g/cm3", 11.34)
    mat_ni60.depletable = True
    mat_ni60.volume = 1.0
    mat_ni60.temperature = 293.6

    mg_flux = np.array([0.5e11] * 42)
    multigroup_flux = list(mg_flux / sum(mg_flux))

    depleted_materials = deplete_materials(
        activation_data=[
            {
                "material": mat_ni58,
                "multigroup_flux": multigroup_flux,
                "energy_group_structure": GROUP_STRUCTURES["VITAMIN-J-42"],
                "source_rate": [1e24, 0, 0],
                "timesteps": [10, 10, 10],
            },
            {
                "material": mat_ni60,
                "multigroup_flux": multigroup_flux,
                "energy_group_structure": "VITAMIN-J-42",
                "source_rate": [1e10, 1e10, 1e10],
                "timesteps": [1000, 1000, 200],
            },
        ],
        timestep_units="s",
        chain_file=Path(__file__).parents[1] / "chain_ni.xml",
        # nuclides=["Ni56", "Ni60"],  # limit to two nuclides to speed up the test
        reactions=None,
    )

    for material in depleted_materials[0]:
        assert "Ni57" in material.get_nuclides()

    mat_2_fe55 = depleted_materials[0][1].get_nuclide_atom_densities("Ni57")["Ni57"]
    mat_3_fe55 = depleted_materials[0][2].get_nuclide_atom_densities("Ni57")["Ni57"]
    assert mat_3_fe55 < mat_2_fe55

    for material in depleted_materials[1]:
        assert "Ni59" in material.get_nuclides()

    mat_2_cu66 = depleted_materials[1][1].get_nuclide_atom_densities("Ni59")["Ni59"]
    mat_3_cu66 = depleted_materials[1][2].get_nuclide_atom_densities("Ni59")["Ni59"]
    assert mat_3_cu66 > mat_2_cu66


def test_deplete_materials_error_handling():

    mat_ni = openmc.Material()
    mat_ni.add_nuclide("Ni58", 1.0)
    mat_ni.set_density("g/cm3", 19.32)

    mg_flux = np.array([0.5e11] * 42)
    multigroup_flux = list(mg_flux / sum(mg_flux))

    activation_data = [
        {
            "material": mat_ni,
            "multigroup_flux": multigroup_flux,
            "energy_group_structure": GROUP_STRUCTURES["VITAMIN-J-42"],
            "source_rate": [1e19, 0, 0],
            "timesteps": [10, 10, 10],
        }
    ]

    with pytest.raises(
        ValueError, match="Material temperature must be set before depletion"
    ):
        deplete_materials(
            activation_data=activation_data,
            timestep_units="s",
            chain_file=Path(__file__).parents[1] / "chain_ni.xml",
            reactions=None,
        )
    mat_ni.temperature = 293.6

    with pytest.raises(
        ValueError, match="Material volume must be set before depletion"
    ):
        deplete_materials(
            activation_data=activation_data,
            timestep_units="s",
            chain_file=Path(__file__).parents[1] / "chain_ni.xml",
            reactions=None,
        )

    mat_ni.volume = 1.0
    with pytest.raises(
        RuntimeError, match="No depletable materials were found in the model"
    ):
        deplete_materials(
            activation_data=activation_data,
            timestep_units="s",
            chain_file=Path(__file__).parents[1] / "chain_ni.xml",
            reactions=None,
        )
    mat_ni.depletable = True
    with pytest.raises(ValueError, match="Number of time steps"):
        deplete_materials(
            activation_data=[
                {
                    "material": mat_ni,
                    "multigroup_flux": multigroup_flux,
                    "energy_group_structure": GROUP_STRUCTURES["VITAMIN-J-42"],
                    "source_rate": [0, 0],  # wrong length
                    "timesteps": [10, 10, 10],
                }
            ],
            timestep_units="s",
            chain_file=Path(__file__).parents[1] / "chain_ni.xml",
            reactions=None,
        )
    with pytest.raises(ValueError, match="Length of flux array should be"):
        deplete_materials(
            activation_data=[
                {
                    "material": mat_ni,
                    "multigroup_flux": multigroup_flux,
                    "energy_group_structure": [1.0, 2.0, 3.0],  # wrong length
                    "source_rate": [1e19, 0, 0],
                    "timesteps": [10, 10, 10],
                }
            ],
            timestep_units="s",
            chain_file=Path(__file__).parents[1] / "chain_ni.xml",
            reactions=None,
        )
