"""Test the Iterated Fission Probability (IFP) method to compute adjoint-weighted
kinetics parameters using dedicated tallies."""

import openmc
import pytest

from tests.testing_harness import PyAPITestHarness

@pytest.fixture()
def ifp_model():
    # Material
    material = openmc.Material(name="core")
    material.add_nuclide("U235", 1.0)
    material.add_nuclide("Pu239", 1.0)
    material.set_density('g/cm3', 16.0)

    # Geometry
    radius = 10.0
    sphere = openmc.Sphere(r=radius, boundary_type="vacuum")
    cell = openmc.Cell(region=-sphere, fill=material)
    geometry = openmc.Geometry([cell])

    # Settings
    settings = openmc.Settings()
    settings.particles = 100000
    settings.batches = 20
    settings.inactive = 5
    settings.ifp_n_generation = 5

    model = openmc.Model(settings=settings, geometry=geometry)

    space = openmc.stats.Box(*cell.bounding_box)
    model.settings.source = openmc.IndependentSource(
        space=space, constraints={'fissionable': True})
    model.add_kinetics_parameters_tallies(num_groups=6, nuclides = ["U235", "Pu239"])
    return model


def test_iterated_fission_probability(ifp_model):
    harness = PyAPITestHarness("statepoint.20.h5", model=ifp_model)
    harness.main()
