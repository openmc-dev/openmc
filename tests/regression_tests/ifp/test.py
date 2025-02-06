"""Test the Iterated Fission Probability (IFP) method to compute adjoint-weighted
kinetics parameters using dedicated tallies."""

import openmc
import pytest

from tests.testing_harness import PyAPITestHarness


@pytest.fixture()
def ifp_model():
    openmc.reset_auto_ids()
    model = openmc.model.Model()

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
    settings = openmc.Settings()
    settings.run_mode = 'eigenvalue'
    settings.particles = 1000
    settings.batches = 20
    settings.inactive = 5
    settings.photon_transport = False
    settings.ifp_n_generation = 5

    bounds = [
        -radius, -radius, -radius,
        radius, radius, radius
    ]
    space = openmc.stats.Box(bounds[:3], bounds[3:])
    settings.source = openmc.IndependentSource(space=space, constraints={'fissionable': True})

    model.settings = settings

    # Tally IFP scores
    tally = openmc.Tally(name="ifp-scores")
    tally.scores = ["ifp-time-numerator", "ifp-beta-numerator", "ifp-denominator"]

    tallies = openmc.Tallies([tally])
    model.tallies = tallies

    return model


def test_iterated_fission_probability(ifp_model):
    harness = PyAPITestHarness("statepoint.20.h5", model=ifp_model)
    harness.main()
