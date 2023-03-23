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

    source = openmc.Source()
    source.space = openmc.stats.Point()
    source.energy = openmc.stats.Discrete([10000], [1.0])

    settings = openmc.Settings()
    settings.particles = 2000
    settings.batches = 10
    settings.run_mode = 'fixed source'

    # build
    mesh = openmc.SphericalMesh()
    mesh.phi_grid = np.linspace(0, 2*np.pi, 13)
    mesh.theta_grid = np.linspace(0, np.pi, 7)
    mesh.r_grid = np.linspace(0, geom_size, geom_size)

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
# TODO: determine why this is needed for spherical mesh
# but not cylindrical mesh
offset = geom_size + 0.001

origins = set(permutations((-offset, 0, 0)))
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
    # move the center of the spherical mesh
    mesh.origin = origin

    sp_filename = model.run()

    with openmc.StatePoint(sp_filename) as sp:
        tally = sp.tallies[1]

        # we've translated half of the spherical mesh above the model,
        # so ensure that half of the bins are populated
        assert np.count_nonzero(tally.mean) == tally.mean.size / 2

        # check that the half of the mesh that is outside of the geometry
        # contains the zero values
        mean = tally.get_reshaped_data('mean', expand_dims=True)
        centroids = mesh.centroids
        for ijk in mesh.indices:
            i, j, k = np.array(ijk) - 1
            print(centroids[:, i, j, k])
            if model.geometry.find(centroids[:, i, j, k]):
                mean[i, j, k] == 0.0
            else:
                mean[i, j, k] != 0.0
