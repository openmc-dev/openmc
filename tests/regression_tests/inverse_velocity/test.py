"""Test tallying inverse-velocity using tracklength and collision estimators"""

import openmc
import pytest

from tests.testing_harness import PyAPITestHarness

@pytest.fixture()
def simple_model():
    model = openmc.Model()

    # Material
    material = openmc.Material(name="core")
    material.add_nuclide("U235", 1.0)
    material.set_density('g/cm3', 16.0)

    # Geometry
    radius = 10.0
    sphere = openmc.Sphere(r=radius, boundary_type="vacuum")
    cell = openmc.Cell(region=-sphere, fill=material)
    model.geometry = openmc.Geometry([cell])

    # Settings
    model.settings.particles = 1000
    model.settings.batches = 20
    model.settings.inactive = 5

    space = openmc.stats.Box(*cell.bounding_box)
    model.settings.source = openmc.IndependentSource(
        space=space, constraints={'fissionable': True})

    # Tally inverse-velocity score
    tally1 = openmc.Tally(name="tracklength")
    tally1.scores = ["inverse-velocity"]
    
    tally2 = openmc.Tally(name="collision")
    tally2.scores = ["inverse-velocity"]
    tally2.estimator = "collision"
    
    tally3 = openmc.Tally(name="analog")
    tally3.scores = ["inverse-velocity"]
    tally3.estimator = "analog"
    
    model.tallies = [tally1, tally2, tally3]

    return model


def test_inverse_velocity(simple_model):
    harness = PyAPITestHarness("statepoint.20.h5", model=simple_model)
    harness.main()
