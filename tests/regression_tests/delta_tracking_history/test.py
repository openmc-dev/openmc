import openmc
import pytest
from openmc.utility_funcs import change_directory

from tests.testing_harness import PyAPITestHarness


@pytest.fixture
def model():
    openmc.reset_auto_ids()
    model = openmc.Model()

    uo2 = openmc.Material(name='UO2')
    uo2.set_density('g/cm3', 10.0)
    uo2.add_nuclide('U235', 1.0)
    uo2.add_nuclide('O16', 2.0)
    water = openmc.Material(name='light water')
    water.add_nuclide('H1', 2.0)
    water.add_nuclide('O16', 1.0)
    water.set_density('g/cm3', 1.0)
    water.add_s_alpha_beta('c_H_in_H2O')
    model.materials.extend([uo2, water])

    cyl = openmc.ZCylinder(r=0.4)
    pin = openmc.model.pin([cyl], [uo2, water])
    d = 1.0

    lattice = openmc.RectLattice()
    lattice.lower_left = (-d, -d)
    lattice.pitch = (d, d)
    lattice.universes = [[pin, pin],
                         [pin, pin]]
    box = openmc.model.RectangularPrism(
        2.0 * d, 2.0 * d,
        origin=(0.0, 0.0),
        boundary_type='reflective'
    )

    model.geometry = openmc.Geometry([openmc.Cell(fill=lattice, region=-box)])
    model.geometry.merge_surfaces = True

    msh = openmc.RegularMesh(mesh_id=0)
    msh.lower_left = (-d, -d)
    msh.upper_right = (d, d)
    msh.dimension = (2, 2)
    t = openmc.Tally()
    t.filters = [openmc.ParticleFilter(bins='neutron'), openmc.MeshFilter(mesh=msh)]
    t.scores = ['total']
    t.estimator = 'collision'
    model.tallies.append(t)

    model.settings.batches = 10
    model.settings.inactive = 5
    model.settings.particles = 1000
    model.settings.delta_tracking = True

    return model


@pytest.mark.parametrize("photon", [False, True])
def test_lattice_checkerboard(model, photon):
    with change_directory(str(photon)):
        model.settings.photon_transport = photon
        if photon:
            msh = openmc.RegularMesh(mesh_id=1)
            msh.lower_left = (-1.0, -1.0)
            msh.upper_right = (1.0, 1.0)
            msh.dimension = (2, 2)
            t = openmc.Tally()
            t.filters = [openmc.ParticleFilter(bins='photon'), openmc.MeshFilter(mesh=msh)]
            t.scores = ['total']
            t.estimator = 'collision'
            model.tallies.append(t)
        harness = PyAPITestHarness('statepoint.10.h5', model)
        harness.main()
