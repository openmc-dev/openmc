"""Test migration area calculation."""

import openmc
import pytest

from tests.testing_harness import PyAPITestHarness

@pytest.fixture()
def model():
    model = openmc.Model()

    # Material
    material = openmc.Material(name="Hydrogen")
    material.add_nuclide("H1", 1.0)
    material.set_density('g/cm3', 1.0)
    material.add_s_alpha_beta('c_H_in_H2O')

    # Geometry
    radius = 10.0
    sphere = openmc.Sphere(r=radius, boundary_type="reflective")
    cell = openmc.Cell(region=-sphere, fill=material)
    model.geometry = openmc.Geometry([cell])

    # Settings
    model.settings.particles = 10000
    model.settings.batches = 20
    model.settings.run_mode = 'fixed source'
    model.settings.source = openmc.IndependentSource()
    

    # Tally 
    groups = openmc.mgxs.EnergyGroups('CASMO-70')
    energy_filter = openmc.EnergyFilter(groups.group_edges)
    tally = openmc.Tally()
    tally.scores = ["migration-area", "flux", "total"]
    tally.filters = [energy_filter]
    model.tallies = [tally]

    return model


def test_migration_area(model):
    harness = PyAPITestHarness("statepoint.20.h5", model=model)
    harness.main()
