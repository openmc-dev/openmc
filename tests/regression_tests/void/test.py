import numpy as np
import openmc
import pytest

from tests.testing_harness import PyAPITestHarness

@pytest.fixture
def model():
    model = openmc.model.Model()

    zn = openmc.Material()
    zn.set_density('g/cm3', 7.14)
    zn.add_nuclide('Zn64', 1.0)
    model.materials.append(zn)

    radii = np.linspace(1.0, 100.0)
    surfs = [openmc.Sphere(r=r) for r in radii]
    surfs[-1].boundary_type = 'vacuum'
    cells = [openmc.Cell(fill=(zn if i % 2 == 0 else None), region=region)
             for i, region in enumerate(openmc.model.subdivide(surfs))]
    model.geometry = openmc.Geometry(cells)

    model.settings.run_mode = 'fixed source'
    model.settings.batches = 3
    model.settings.particles = 1000
    model.settings.source = openmc.Source(space=openmc.stats.Point())

    cell_filter = openmc.CellFilter(cells)
    tally = openmc.Tally()
    tally.filters = [cell_filter]
    tally.scores = ['total']
    model.tallies.append(tally)

    return model


def test_void(model):
    harness = PyAPITestHarness('statepoint.3.h5', model)
    harness.main()
