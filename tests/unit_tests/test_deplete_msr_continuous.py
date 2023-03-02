""" Tests for Msr continuos class """

from pathlib import Path
from math import exp

import pytest
import openmc
from openmc.deplete import CoupledOperator
from openmc.deplete import MsrContinuous
import numpy as np

CHAIN_PATH = Path(__file__).parents[1] / "chain_simple.xml"

@pytest.fixture
def model():
    f = openmc.Material(name="f")
    f.add_element("U", 1, percent_type="ao", enrichment=4.25)
    f.add_element("O", 2)
    f.set_density("g/cc", 10.4)

    w = openmc.Material(name="w")
    w.add_element("O", 1)
    w.add_element("H", 2)
    w.set_density("g/cc", 1.0)
    w.depletable = True

    radii = [0.42, 0.45]
    f.volume = np.pi * radii[0] ** 2
    w.volume = np.pi * (radii[1]**2 - radii[0]**2)

    materials = openmc.Materials([f, w])

    surf_f = openmc.Sphere(r=radii[0])
    surf_w = openmc.Sphere(r=radii[1], boundary_type='vacuum')
    cell_f = openmc.Cell(fill=f, region=-surf_f)
    cell_w = openmc.Cell(fill=w, region=+surf_f & -surf_w)
    geometry = openmc.Geometry([cell_f,cell_w])

    settings = openmc.Settings()
    settings.particles = 1000
    settings.inactive = 10
    settings.batches = 50

    return openmc.Model(geometry, materials, settings)

def test_get_set(model):
    """Tests the get/set methods"""

    # create removal rates for U and Xe
    removal_rates = {'U': 0.01, 'Xe': 0.1}
    op = CoupledOperator(model, CHAIN_PATH)
    msr = MsrContinuous(op, model)

    # Test by Openmc material, material name and material id
    material, dest_material = [m for m in model.materials if m.depletable]
    material_inputs = [material, material.name, material.id]:
    dest_material_inputs = [dest_material,
                            dest_material.name,
                            dest_material.id]:
    material_inputs = np.array(np.meshgrid(material_inputs
                                     dest_material_inputs)).T.reshape(-1, 2)
    for material_input, dest_material_input in material_inputs:
        for element, removal_rate in removal_rates.items():
            msr.set_removal_rate(material_input, [element], removal_rate,
                                 destination_material=dest_material_input)
            assert msr.get_removal_rate(material_input, element) == removal_rate
            assert msr.get_destination_material(material_input, element) == str(dest_mat.id)
        assert msr.get_elements(material_input) == removal_rates.keys()
