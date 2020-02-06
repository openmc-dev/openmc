import openmc
import pytest

from tests.testing_harness import PyAPITestHarness


@pytest.fixture
def model():
    model = openmc.model.Model()
    mat = openmc.Material()
    mat.set_density('g/cm3', 2.6989)
    mat.add_nuclide('Al27', 1.0)
    model.materials.append(mat)

    cyl = openmc.XCylinder(r=1.0, boundary_type='vacuum')
    x_plane_left = openmc.XPlane(-1.0, 'vacuum')
    x_plane_center = openmc.XPlane(1.0)
    x_plane_right = openmc.XPlane(1.0e9, 'vacuum')

    inner_cyl_left = openmc.Cell()
    inner_cyl_right = openmc.Cell()
    outer_cyl = openmc.Cell()

    inner_cyl_left.region = -cyl & +x_plane_left & -x_plane_center
    inner_cyl_right.region = -cyl & +x_plane_center & -x_plane_right
    outer_cyl.region = ~(-cyl & +x_plane_left & -x_plane_right)
    inner_cyl_right.fill = mat
    model.geometry = openmc.Geometry([inner_cyl_left, inner_cyl_right, outer_cyl])

    source = openmc.Source()
    source.space = openmc.stats.Point((0, 0, 0))
    source.angle = openmc.stats.Monodirectional()
    source.energy = openmc.stats.Discrete([14.0e6], [1.0])
    source.particle = 'neutron'

    model.settings.particles = 10000
    model.settings.run_mode = 'fixed source'
    model.settings.batches = 1
    model.settings.photon_transport = True
    model.settings.electron_treatment = 'ttb'
    model.settings.cutoff = {'energy_photon' : 1000.0}
    model.settings.source = source

    surface_filter = openmc.SurfaceFilter(cyl)
    particle_filter = openmc.ParticleFilter(['neutron', 'photon', 'electron', 'positron'])
    current_tally = openmc.Tally()
    current_tally.filters = [surface_filter, particle_filter]
    current_tally.scores = ['current']
    tally_tracklength = openmc.Tally()
    tally_tracklength.filters = [particle_filter]
    tally_tracklength.scores = ['total']  # heating doesn't work with tracklength
    tally_tracklength.nuclides = ['Al27', 'total']
    tally_tracklength.estimator = 'tracklength'
    tally_collision = openmc.Tally()
    tally_collision.filters = [particle_filter]
    tally_collision.scores = ['total', 'heating']
    tally_collision.nuclides = ['Al27', 'total']
    tally_collision.estimator = 'collision'
    tally_analog = openmc.Tally()
    tally_analog.filters = [particle_filter]
    tally_analog.scores = ['total', 'heating']
    tally_analog.nuclides = ['Al27', 'total']
    tally_analog.estimator = 'analog'
    model.tallies.extend([current_tally, tally_tracklength,
                          tally_collision, tally_analog])

    return model


def test_photon_production(model):
    harness = PyAPITestHarness('statepoint.1.h5', model)
    harness.main()
