import openmc
import openmc.model
import pytest

from tests.testing_harness import PyAPITestHarness


@pytest.fixture
def model():
    model = openmc.model.Model()

    # Materials
    m1 = openmc.Material()
    m1.set_density('g/cc', 4.5)
    m1.add_nuclide('U235', 1.0)
    m2 = openmc.Material()
    m2.set_density('g/cc', 1.0)
    m2.add_nuclide('H1', 1.0)
    model.materials += [m1, m2]

    # Geometry
    cyl1 = openmc.ZCylinder(r=0.7)
    c1 = openmc.Cell(fill=m1, region=-cyl1)
    c2 = openmc.Cell(fill=m2, region=+cyl1)
    u1 = openmc.Universe(cells=[c1, c2])

    cyl2 = openmc.ZCylinder(r=0.5)
    c3 = openmc.Cell(fill=m1, region=-cyl2)
    c4 = openmc.Cell(fill=m2, region=+cyl2)
    u2 = openmc.Universe(cells=[c3, c4])

    lat = openmc.RectLattice()
    lat.lower_left = (-4, -4)
    lat.pitch = (2, 2)
    lat.universes = [
        [u1, u2, u2, u2],
        [u2, u1, u2, u2],
        [u2, u2, u1, u2],
        [u2, u2, u2, u1]
    ]
    box = openmc.model.rectangular_prism(8.0, 8.0, boundary_type='reflective')
    main_cell = openmc.Cell(fill=lat, region=box)
    model.geometry.root_universe = openmc.Universe(cells=[main_cell])
    model.geometry.determine_paths()

    # Settings
    model.settings.batches = 5
    model.settings.inactive = 0
    model.settings.particles = 1000
    model.settings.source = openmc.Source(space=openmc.stats.Point())

    instances = ([(c3, i) for i in range(c3.num_instances)] +
                 [(c2, i) for i in range(c2.num_instances)])
    f1 = openmc.CellInstanceFilter(instances)
    f2 = openmc.CellInstanceFilter(instances[::-1])
    t1 = openmc.Tally()
    t1.filters = [f1]
    t1.scores = ['total']
    t2 = openmc.Tally()
    t2.filters = [f2]
    t2.scores = ['total']
    model.tallies += [t1, t2]

    return model


def test_cell_instance(model):
    harness = PyAPITestHarness('statepoint.5.h5', model)
    harness.main()
