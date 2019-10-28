"""The purpose of this test is to provide coverage of energy distributions that
are not covered in other tests. It has a single material with the following
nuclides:

U233: Only nuclide that has a Watt fission spectrum

Am244: One of a few nuclides that has a Maxwell fission spectrum

H2: Only nuclide that has an N-body phase space distribution, in this case for
(n,2n)

Na23: Has an evaporation spectrum and also has reactions that have multiple
angle-energy distributions, so it provides coverage for both of those
situations.

Ta181: One of a few nuclides that has reactions with Kalbach-Mann distributions
that use linear-linear interpolation.

"""

import openmc
import pytest

from tests.testing_harness import PyAPITestHarness


@pytest.fixture
def model():
    model = openmc.model.Model()

    m = openmc.Material()
    m.set_density('g/cm3', 20.0)
    m.add_nuclide('U233', 1.0)
    m.add_nuclide('Am244', 1.0)
    m.add_nuclide('H2', 1.0)
    m.add_nuclide('Na23', 1.0)
    m.add_nuclide('Ta181', 1.0)

    s = openmc.Sphere(r=100.0, boundary_type='reflective')
    c = openmc.Cell(fill=m, region=-s)
    model.geometry = openmc.Geometry([c])

    model.settings.batches = 10
    model.settings.inactive = 5
    model.settings.particles = 1000

    return model


def test_energy_laws(model):
    harness = PyAPITestHarness('statepoint.10.h5', model)
    harness.main()
