"""Test one-group cross section generation"""
from pathlib import Path

import numpy as np
import pytest
import openmc
from openmc.deplete import MicroXS, get_microxs_and_flux

from tests.regression_tests import config

CHAIN_FILE = Path(__file__).parents[2] / "chain_simple.xml"

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

    materials = openmc.Materials([fuel, clad, water])

    pin_surfaces = [openmc.ZCylinder(r=r) for r in radii]
    pin_univ = openmc.model.pin(pin_surfaces, materials)
    bound_box = openmc.rectangular_prism(1.24, 1.24, boundary_type="reflective")
    root_cell = openmc.Cell(fill=pin_univ, region=bound_box)
    geometry = openmc.Geometry([root_cell])

    settings = openmc.Settings()
    settings.particles = 1000
    settings.inactive = 5
    settings.batches = 10

    return openmc.Model(geometry, materials, settings)


@pytest.mark.parametrize("domain_type", ["materials", "mesh"])
def test_from_model(model, domain_type):
    if domain_type == 'materials':
        domains = model.materials[:1]
    elif domain_type == 'mesh':
        mesh = openmc.RegularMesh()
        mesh.lower_left = (-0.62, -0.62)
        mesh.upper_right = (0.62, 0.62)
        mesh.dimension = (3, 3)
        domains = mesh
    nuclides = ['U234', 'U235', 'U238', 'U236', 'O16', 'O17', 'I135', 'Xe135',
                'Xe136', 'Cs135', 'Gd157', 'Gd156']
    _, test_xs = get_microxs_and_flux(model, domains, nuclides, chain_file=CHAIN_FILE)
    if config['update']:
        test_xs[0].to_csv(f'test_reference_{domain_type}.csv')

    ref_xs = MicroXS.from_csv(f'test_reference_{domain_type}.csv')

    np.testing.assert_allclose(test_xs[0].data, ref_xs.data, rtol=1e-11)
