import numpy as np
from math import pi

import openmc
import pytest

from tests.testing_harness import PyAPITestHarness


@pytest.fixture
def model():
    model = openmc.Model()
    fuel = openmc.Material()
    fuel.set_density("g/cm3", 10.0)
    fuel.add_nuclide("U235", 1.0)
    zr = openmc.Material()
    zr.set_density("g/cm3", 1.0)
    zr.add_nuclide("Zr90", 1.0)

    cyl1 = openmc.ZCylinder(r=1.0)
    cyl2 = openmc.ZCylinder(r=3.0, boundary_type="vacuum")
    cell1 = openmc.Cell(fill=fuel, region=-cyl1)
    cell2 = openmc.Cell(fill=zr, region=+cyl1 & -cyl2)
    model.geometry = openmc.Geometry([cell1, cell2])

    model.settings.batches = 5
    model.settings.inactive = 0
    model.settings.particles = 1000

    # Create a tally for current through the first surface binned by mu
    surf_filter = openmc.SurfaceFilter([cyl1])
    mu_filter = openmc.MuSurfaceFilter([-1.0, -0.5, 0.0, 0.5, 1.0])
    tally = openmc.Tally()
    tally.filters = [surf_filter, mu_filter]
    tally.scores = ["current"]
    model.tallies.append(tally)

    return model


def test_filter_musurface(model):
    harness = PyAPITestHarness("statepoint.5.h5", model)
    harness.main()
