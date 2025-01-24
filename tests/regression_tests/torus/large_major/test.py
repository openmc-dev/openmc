import openmc
import pytest

from tests.testing_harness import PyAPITestHarness


@pytest.fixture
def model():
    model = openmc.Model()
    tungsten = openmc.Material()
    tungsten.set_density('g/cm3', 1.0)
    tungsten.add_nuclide('W184', 1.0)
    ss = openmc.Material()
    ss.set_density('g/cm3', 5.0)
    ss.add_nuclide('Fe56', 1.0)
    model.materials.extend([tungsten, ss])

    # Create nested torii with very large major radii
    R = 1000.0
    vacuum = openmc.ZTorus(a=R, b=30.0, c=30.0)
    first_wall = openmc.ZTorus(a=R, b=35.0, c=35.0)
    vessel = openmc.ZTorus(a=R, b=40.0, c=40.0, boundary_type='vacuum')
    cell1 = openmc.Cell(region=-vacuum)
    cell2 = openmc.Cell(fill=tungsten, region=+vacuum & -first_wall)
    cell3 = openmc.Cell(fill=ss, region=+first_wall & -vessel)
    model.geometry = openmc.Geometry([cell1, cell2, cell3])

    model.settings.run_mode ='fixed source'
    model.settings.particles = 1000
    model.settings.batches = 10
    model.settings.source = openmc.IndependentSource(space=openmc.stats.Point((-R, 0, 0,)))

    tally = openmc.Tally()
    tally.scores = ['flux']
    model.tallies.append(tally)
    return model


def test_torus_large_major(model):
    harness = PyAPITestHarness('statepoint.10.h5', model)
    harness.main()
