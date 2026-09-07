""" Tests for KeffSearchControl class """

from pathlib import Path

import pytest
import numpy as np

import openmc
import openmc.lib
from openmc.deplete import CoupledOperator

CHAIN_PATH = Path(__file__).parents[1] / "chain_simple.xml"


def make_model():
    f = openmc.Material(name="fuel")
    f.add_element("U", 1, percent_type="ao", enrichment=4.25)
    f.add_element("O", 2)
    f.set_density("g/cc", 10.4)
    f.temperature = 293.15

    w = openmc.Material(name="water")
    w.add_element("O", 1)
    w.add_element("H", 2)
    w.set_density("g/cc", 1.0)
    w.temperature = 293.15
    w.depletable = True

    h = openmc.Material(name='helium')
    h.add_element('He', 1)
    h.set_density('g/cm3', 0.001598)

    radii = [0.42, 0.45]
    height = 0.5

    f.volume = np.pi * radii[0] ** 2 * height
    w.volume = np.pi * (radii[1]**2 - radii[0]**2) * height/2

    materials = openmc.Materials([f, w, h])

    surf_interface = openmc.ZPlane(z0=0)
    surf_top = openmc.ZPlane(z0=height/2)
    surf_bot = openmc.ZPlane(z0=-height/2)
    surf_in = openmc.Sphere(r=radii[0])
    surf_out = openmc.Sphere(r=radii[1], boundary_type='vacuum')

    cell_water = openmc.Cell(fill=w, region=-surf_interface)
    cell_helium = openmc.Cell(fill=h, region=+surf_interface)
    universe = openmc.Universe(cells=(cell_water, cell_helium))
    cell_fuel = openmc.Cell(name='fuel_cell', fill=f,
                            region=-surf_in & -surf_top & +surf_bot)
    cell_universe = openmc.Cell(name='universe_cell',fill=universe,
                            region=+surf_in & -surf_out & -surf_top & +surf_bot)
    geometry = openmc.Geometry([cell_fuel, cell_universe])

    settings = openmc.Settings()
    settings.particles = 1000
    settings.inactive = 10
    settings.batches = 50

    return openmc.Model(geometry, materials, settings)


def translate_cell(position):
    """Helper function to translate a cell"""
    cell = [c for c in openmc.lib.cells.values() if c.name == 'universe_cell'][0]
    openmc.lib.cells[cell.id].translation = [0, 0, position]
    return position


def rotate_cell(angle):
    """Helper function to rotate a cell"""
    cell = [c for c in openmc.lib.cells.values() if c.name == 'universe_cell'][0]
    openmc.lib.cells[cell.id].rotation = [0, 0, angle]
    return angle


def set_u235_density(u235_density):
    """Helper function to set the U235 density directly"""
    fuel = [m for m in openmc.lib.materials.values() if m.name == 'fuel'][0]
    nuclides = openmc.lib.materials[fuel.id].nuclides
    densities = openmc.lib.materials[fuel.id].densities
    u235_idx = nuclides.index('U235')
    densities[u235_idx] = u235_density
    openmc.lib.materials[fuel.id].set_densities(nuclides, densities)
    return u235_density


@pytest.mark.parametrize("function, x0, x1, bracket", [
    (translate_cell, -1.0, 1.0, (-5.0, 5.0)),
    (rotate_cell, -45.0, 45.0, (-90.0, 90.0)),
    (set_u235_density, 0.8, 1.2, (0.5, 1.5))
])
def test_integrator_add_keff_search_control(run_in_tmpdir, function, x0, x1, bracket):
    """Test adding add_keff_search_control to integrator"""
    model = make_model()
    operator = CoupledOperator(model, CHAIN_PATH)
    integrator = openmc.deplete.PredictorIntegrator(
        operator, [1, 1], 0.0, timestep_units='d')

    integrator.add_keff_search_control(
        function=function,
        x0=x0,
        x1=x1,
        bracket=bracket,
        k_tol=0.1,
        output=False,
    )

    assert integrator._keff_search_control.x0 == x0
    assert integrator._keff_search_control.x1 == x1
    assert integrator._keff_search_control.function == function
    assert integrator._keff_search_control.search_kwargs['x_min'] == bracket[0]
    assert integrator._keff_search_control.search_kwargs['x_max'] == bracket[1]
    assert integrator._keff_search_control.search_kwargs['k_tol'] == 0.1
    assert not integrator._keff_search_control.search_kwargs['output']
