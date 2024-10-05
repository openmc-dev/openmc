"""Basic unit tests for openmc.deplete.CoupledOperator instantiation

"""

from pathlib import Path

import pytest
from openmc.deplete import CoupledOperator
import openmc
import numpy as np

CHAIN_PATH = Path(__file__).parents[1] / "chain_simple.xml"


@pytest.fixture(scope="module")
def model():
    fuel = openmc.Material(name="uo2")
    fuel.add_element("U", 1, percent_type="ao", enrichment=4.25)
    fuel.add_element("O", 2)
    fuel.set_density("g/cc", 10.4)

    clad = openmc.Material(name="clad")
    clad.add_element("Zr", 1)
    clad.set_density("g/cc", 6)

    water = openmc.Material(name="water")
    water.add_element("O", 1)
    water.add_element("H", 2)
    water.set_density("g/cc", 1.0)
    water.add_s_alpha_beta("c_H_in_H2O")

    radii = [0.42, 0.45]
    fuel.volume = np.pi * radii[0] ** 2
    clad.volume = np.pi * (radii[1]**2 - radii[0]**2)
    water.volume = 1.24**2 - (np.pi * radii[1]**2)

    materials = openmc.Materials([fuel, clad, water])

    pin_surfaces = [openmc.ZCylinder(r=r) for r in radii]
    pin_univ = openmc.model.pin(pin_surfaces, materials)
    bound_box = openmc.model.RectangularPrism(
        1.24, 1.24, boundary_type="reflective")
    root_cell = openmc.Cell(fill=pin_univ, region=-bound_box)
    geometry = openmc.Geometry([root_cell])

    settings = openmc.Settings()
    settings.particles = 1000
    settings.inactive = 10
    settings.batches = 50

    return openmc.Model(geometry, materials, settings)


@pytest.fixture()
def model_with_volumes():
    mat1 = openmc.Material()
    mat1.add_element("Ag", 1, percent_type="ao")
    mat1.set_density("g/cm3", 10.49)
    mat1.depletable = True
    mat1.volume = 102

    mat2 = openmc.Material()
    mat2.add_element("Ag", 1, percent_type="ao")
    mat2.set_density("g/cm3", 10.49)

    sph1 = openmc.Sphere(r=1.0)
    sph2 = openmc.Sphere(r=2.0, x0=3)
    sph3 = openmc.Sphere(r=5.0, boundary_type="vacuum")

    cell1 = openmc.Cell(region=-sph1, fill=mat1)
    cell1.volume = 4.19
    cell2 = openmc.Cell(region=-sph2, fill=mat1)
    cell2.volume = 33.51
    cell3 = openmc.Cell(region=-sph3 & +sph1 & +sph2, fill=mat2)
    cell3.volume = 485.9

    geometry = openmc.Geometry([cell1, cell2, cell3])

    return openmc.Model(geometry)


def test_operator_init(model):
    """The test uses a temporary dummy chain. This file will be removed
    at the end of the test, and only contains a depletion_chain node."""

    CoupledOperator(model, CHAIN_PATH)


def test_diff_volume_method_match_cell(model_with_volumes):
    """Tests the volumes assigned to the materials match the cell volumes"""

    operator = openmc.deplete.CoupledOperator(
        model=model_with_volumes,
        diff_burnable_mats=True,
        diff_volume_method='match cell',
        chain_file=CHAIN_PATH
    )

    all_cells = list(operator.model.geometry.get_all_cells().values())
    assert all_cells[0].fill.volume == 4.19
    assert all_cells[1].fill.volume == 33.51
    # mat2 is not depletable
    assert all_cells[2].fill.volume is None


def test_diff_volume_method_divide_equally(model_with_volumes):
    """Tests the volumes assigned to the materials are divided equally"""

    operator = openmc.deplete.CoupledOperator(
        model=model_with_volumes,
        diff_burnable_mats=True,
        diff_volume_method='divide equally',
        chain_file=CHAIN_PATH
    )

    all_cells = list(operator.model.geometry.get_all_cells().values())
    assert all_cells[0].fill[0].volume == 51
    assert all_cells[1].fill[0].volume == 51
    # mat2 is not depletable
    assert all_cells[2].fill.volume is None
