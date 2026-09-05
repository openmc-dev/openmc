import openmc
import pytest

from tests.testing_harness import PyAPITestHarness


@pytest.fixture
def temperature_field_model():
    """Create a temperature field from a regular mesh over a box with
    different temperature for each cell. This box is surrounded by a larger
    box where the temperature field is not defined and has reflective
    boundary conditions.

    """
    model = openmc.Model()

    # Material
    mat = openmc.Material()
    mat.add_nuclide("U235", 0.2)
    mat.add_nuclide("U238", 0.8)
    mat.add_element("O", 2.0)
    mat.add_element("H", 4.0)
    mat.set_density("g/cm3", 5.0)
    mat.add_s_alpha_beta('c_H_in_H2O')
    model.materials = openmc.Materials([mat])

    # Create mesh
    dim = 2
    lower_left = (0.0, 0.0, 0.0)
    upper_right = (5.0, 5.0, 5.0)
    mesh = openmc.RegularMesh()
    mesh.lower_left = lower_left
    mesh.upper_right = upper_right
    mesh.dimension = (dim, dim, dim)

    # Temperature values
    temperature_values = [294.0 + i * 100 for i in range(dim**3)]

    # Create the temperature field
    temperature_field = openmc.TemperatureField(mesh, temperature_values)

    # Create geometry
    box = openmc.model.RectangularParallelepiped(
        xmin=lower_left[0], xmax=upper_right[0],
        ymin=lower_left[1], ymax=upper_right[1],
        zmin=lower_left[2], zmax=upper_right[2]
    )
    cell_1 = openmc.Cell(fill=mat, region=-box)

    surrounding_box = openmc.model.RectangularParallelepiped(
        xmin=-5.0, xmax=10.0,
        ymin=-5.0, ymax=10.0,
        zmin=-5.0, zmax=10.0,
        boundary_type="reflective"
    )
    cell_2 = openmc.Cell(fill=mat, region=(-surrounding_box & +box))
    cell_2.temperature = 250.0

    model.geometry = openmc.Geometry([cell_1, cell_2])

    # Settings
    settings = openmc.Settings()
    settings.batches = 20
    settings.particles = 200
    settings.temperature_field = temperature_field
    spatial_dist = openmc.stats.Box((-5.0, -5.0, -5.0), (10.0, 10.0, 10.0))
    settings.source = openmc.IndependentSource(
        space=spatial_dist, constraints={"fissionable": True})
    settings.temperature = {'tolerance': 1000, 'multipole': True}
    model.settings = settings

    # Add tallies
    mesh_filter = openmc.MeshFilter(mesh)
    mesh_tally = openmc.Tally(name="total reaction rate")
    mesh_tally.filters = [mesh_filter]
    mesh_tally.scores = ["total"]
    tallies = openmc.Tallies([mesh_tally])
    model.tallies = tallies

    return model


def test_temperature_field(temperature_field_model):
    harness = PyAPITestHarness('statepoint.20.h5', temperature_field_model)
    harness.main()
