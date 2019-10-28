import openmc
import pytest

from tests.testing_harness import PyAPITestHarness


@pytest.fixture
def model():
    model = openmc.model.Model()

    fuel = openmc.Material()
    fuel.set_density('g/cm3', 10.0)
    fuel.add_nuclide('U235', 1.0)
    zr = openmc.Material()
    zr.set_density('g/cm3', 1.0)
    zr.add_nuclide('Zr90', 1.0)
    model.materials.extend([fuel, zr])

    box1 = openmc.model.rectangular_prism(10.0, 10.0)
    box2 = openmc.model.rectangular_prism(20.0, 20.0, boundary_type='reflective')
    top = openmc.ZPlane(z0=10.0, boundary_type='vacuum')
    bottom = openmc.ZPlane(z0=-10.0, boundary_type='vacuum')
    cell1 = openmc.Cell(fill=fuel, region=box1 & +bottom & -top)
    cell2 = openmc.Cell(fill=zr, region=~box1 & box2 & +bottom & -top)
    model.geometry = openmc.Geometry([cell1, cell2])

    model.settings.batches = 5
    model.settings.inactive = 0
    model.settings.particles = 1000


    mesh = openmc.Mesh()
    mesh.lower_left = (-10.0, -10.0, -10.0)
    mesh.upper_right = (10.0, 10.0, 10.0)
    mesh.dimension = (3, 3, 3)

    mesh_surface_filter = openmc.MeshSurfaceFilter(mesh)
    energy_filter = openmc.EnergyFilter([0.0, 0.253, 20.0e6])

    tally1 = openmc.Tally()
    tally1.filters = [mesh_surface_filter]
    tally1.scores = ['current']
    tally2 = openmc.Tally()
    tally2.filters = [mesh_surface_filter, energy_filter]
    tally2.scores = ['current']
    model.tallies.extend([tally1, tally2])

    return model


def test_score_current(model):
    harness = PyAPITestHarness('statepoint.5.h5', model)
    harness.main()
