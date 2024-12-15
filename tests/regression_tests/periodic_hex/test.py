import openmc
import pytest

from tests.testing_harness import PyAPITestHarness


@pytest.fixture
def hex_model():
    model = openmc.model.Model()

    fuel = openmc.Material()
    fuel.add_nuclide("U235", 1.0)
    fuel.set_density("g/cc", 4.5)

    hex_prism = openmc.model.HexagonalPrism(10.0, boundary_type="periodic")
    cell = openmc.Cell(fill=fuel, region=-hex_prism)
    model.geometry = openmc.Geometry([cell])

    # Define settings
    model.settings.particles = 1000
    model.settings.batches = 5
    model.settings.inactive = 0
    return model


def test_periodic_hex(hex_model):
    harness = PyAPITestHarness("statepoint.5.h5", hex_model)
    harness.main()
