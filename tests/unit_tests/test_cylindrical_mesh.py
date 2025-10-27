from itertools import product, permutations

import openmc
import numpy as np

import pytest

geom_size = 5

@pytest.fixture()
def model():
    openmc.reset_auto_ids()

    water = openmc.Material(name='water')
    water.add_element('H', 2.0)
    water.add_element('O', 1.0)
    water.set_density('g/cc', 1.0)

    rpp = openmc.model.RectangularParallelepiped(*([-geom_size, geom_size] * 3),
                                                 boundary_type='vacuum')

    cell = openmc.Cell(region=-rpp, fill=water)

    geom = openmc.Geometry([cell])

    source = openmc.IndependentSource()
    source.space = openmc.stats.Point()
    source.energy = openmc.stats.Discrete([10000], [1.0])

    settings = openmc.Settings()
    settings.particles = 2000
    settings.batches = 10
    settings.run_mode = 'fixed source'

    # build
    mesh = openmc.CylindricalMesh(
        phi_grid=np.linspace(0, 2*np.pi, 21),
        z_grid=np.linspace(-geom_size, geom_size, 11),
        r_grid=np.linspace(0, geom_size, geom_size)
    )
    tally = openmc.Tally()

    mesh_filter = openmc.MeshFilter(mesh)
    tally.filters.append(mesh_filter)

    tally.scores.append("flux")

    tallies = openmc.Tallies([tally])

    return openmc.Model(geometry=geom, settings=settings, tallies=tallies)


def test_origin_read_write_to_xml(run_in_tmpdir, model):
    """Tests that the origin attribute can be written and read back to XML
    """
    mesh = model.tallies[0].filters[0].mesh
    mesh.origin = [0.1, 0.2, 0.3]
    model.tallies.export_to_xml()

    # read back
    new_tallies = openmc.Tallies.from_xml()
    new_tally = new_tallies[0]
    new_mesh = new_tally.filters[0].mesh
    np.testing.assert_equal(new_mesh.origin, mesh.origin)

estimators = ('tracklength', 'collision')
offset = geom_size + 0.001
origins = set(permutations((-offset , 0, 0)))
origins |= set(permutations((offset, 0, 0)))

test_cases = product(estimators, origins)

def label(p):
    if isinstance(p, tuple):
        return f'origin:{p}'
    if isinstance(p, str):
        return f'estimator:{p}'


@pytest.mark.parametrize('estimator,origin', test_cases, ids=label)
def test_offset_mesh(model, estimator, origin):
    """Tests that the mesh has been moved based on tally results
    """
    mesh = model.tallies[0].filters[0].mesh
    model.tallies[0].estimator = estimator
    # move the center of the cylinder mesh upwards
    mesh.origin = origin

    sp_filename = model.run()

    with openmc.StatePoint(sp_filename) as sp:
        tally = sp.tallies[1]

        # we've translated half of the cylinder mesh above the model,
        # so ensure that half of the bins are populated
        assert np.count_nonzero(tally.mean) == tally.mean.size / 2

        # check that the half of the mesh that is outside of the geometry
        # contains the zero values
        mean = tally.get_reshaped_data('mean', expand_dims=True)
        centroids = mesh.centroids
        for ijk in mesh.indices:
            i, j, k = np.array(ijk) - 1
            if model.geometry.find(centroids[i, j, k]):
                mean[i, j, k] == 0.0
            else:
                mean[i, j, k] != 0.0


@pytest.fixture()
def void_coincident_geom_model():
    """A model with many geometric boundaries coincident with mesh boundaries
       across many scales
    """
    openmc.reset_auto_ids()
    model = openmc.model.Model()

    model.materials = openmc.Materials()
    radii = [0.1,1, 5, 50, 100, 150, 250]
    cylinders = [openmc.ZCylinder(r=ri) for ri in radii]
    cylinders[-1].boundary_type = 'vacuum'

    regions = openmc.model.subdivide(cylinders)[:-1]
    cells = [openmc.Cell(region=r, fill=None) for r in regions]
    geom = openmc.Geometry(cells)

    model.geometry = geom

    settings = openmc.Settings(run_mode='fixed source')
    settings.batches = 2
    settings.particles = 1000
    model.settings = settings

    mesh = openmc.CylindricalMesh(
        r_grid=np.linspace(0, 250, 501),
        z_grid=[-250, 250],
        phi_grid=np.linspace(0, 2*np.pi, 2),
    )
    mesh_filter = openmc.MeshFilter(mesh)

    tally = openmc.Tally()
    tally.scores = ['flux']
    tally.filters = [mesh_filter]

    model.tallies = openmc.Tallies([tally])

    return model


# convenience function for checking tally results
# in the following tests
def _check_void_cylindrical_tally(statepoint_filename):
    with openmc.StatePoint(statepoint_filename) as sp:
        flux_tally = sp.tallies[1]
        mesh = flux_tally.find_filter(openmc.MeshFilter).mesh
        neutron_flux = flux_tally.get_reshaped_data().squeeze()
        # we expect the tally results to be the same as the mesh grid width
        # for these cases
        d_r = mesh.r_grid[1] - mesh.r_grid[0]
        assert neutron_flux == pytest.approx(d_r)


def test_void_geom_pnt_src(run_in_tmpdir, void_coincident_geom_model):
    src = openmc.IndependentSource()
    src.space = openmc.stats.Point()
    src.angle = openmc.stats.PolarAzimuthal(mu=openmc.stats.Discrete([0.0], [1.0]))
    src.energy = openmc.stats.Discrete([14.06e6], [1])
    void_coincident_geom_model.settings.source = src

    sp_filename = void_coincident_geom_model.run()
    _check_void_cylindrical_tally(sp_filename)


def test_void_geom_boundary_src(run_in_tmpdir, void_coincident_geom_model):
    # update source to a number of points on the outside of the cylinder
    # with directions pointing toward the origin
    bbox = void_coincident_geom_model.geometry.bounding_box

    # can't source particle directly on the geometry boundary
    outer_r = bbox[1][0] - 1e-08

    n_sources = 100
    radial_vals = np.linspace(0.0, 2.0*np.pi, n_sources)

    sources = []

    energy = openmc.stats.Discrete([14.06e6], [1])
    for val in radial_vals:
        src = openmc.IndependentSource()
        src.energy = energy

        pnt = np.array([np.cos(val), np.sin(val), 0.0])
        u = -pnt
        src.space = openmc.stats.Point(outer_r*pnt)
        src.angle = openmc.stats.Monodirectional(u)
        src.strength = 0.5/n_sources
        sources.append(src)

    void_coincident_geom_model.settings.source = sources
    sp_filename = void_coincident_geom_model.run()

    _check_void_cylindrical_tally(sp_filename)