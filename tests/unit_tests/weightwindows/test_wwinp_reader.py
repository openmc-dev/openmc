from pathlib import Path

import numpy as np
import openmc
import pytest

# check that we can successfully read wwinp files with the following contents:
#
#  - neutrons on a rectilinear mesh
#  - neutrons and photons on a rectilinear mesh

# check that the following raises the correct exceptions (for now):
#
#  - wwinp file with multiple time steps
#  - wwinp file with cylindrical or spherical mesh


# expected retults - neutron data only
n_mesh = openmc.RectilinearMesh()
n_mesh.x_grid = np.array([-100.0,
                            -99.0,
                            -97.0,
                            -79.3636,
                            -61.7273,
                            -44.0909,
                            -26.4546,
                            -8.81818,
                            8.81818,
                            26.4546,
                            44.0909,
                            61.7273,
                            79.3636,
                            97.0,
                            99.0,
                            100])
n_mesh.y_grid = np.array([-100.0,
                            -50.0,
                            -13.3333,
                            23.3333,
                            60.0,
                            70.0,
                            80.0,
                            90.0,
                            100.0])
n_mesh.z_grid = np.array([-100.0,
                            -66.6667,
                            -33.3333,
                            0.0,
                            33.3333,
                            66.6667,
                            100.0])
n_e_bounds = (np.array([0.0,
                         100000.0,
                         146780.0]),)
n_particles = ('neutron',)

# expected results - neutron and photon data
np_mesh = openmc.RectilinearMesh()
np_mesh.x_grid = np.array([-100.0, 100.0])
# y grid and z grid are the same as the previous mesh
np_mesh.y_grid = n_mesh.y_grid
np_mesh.z_grid = n_mesh.z_grid

np_e_bounds = (np.array([0.0, 100000.0, 146780.0, 215440.0]),
               np.array([0.0, 1.0E8]))
np_particles = ('neutron', 'photon')

# expected results - photon data only
p_mesh = openmc.RectilinearMesh()
# adopts z grid from previous meshes as its x grid
p_mesh.x_grid = np_mesh.z_grid
# uses the same y grid
p_mesh.y_grid = np_mesh.y_grid
p_mesh.z_grid = np.array([-50.0, 50.0])

p_e_bounds = (np.array([0.0, 100000.0, 146780.0, 215440.0, 316230.0]),)
p_particles = ('photon',)

expected_results = [('wwinp_n', n_mesh, n_particles, n_e_bounds),
                    ('wwinp_np', np_mesh, np_particles, np_e_bounds),
                    ('wwinp_p', p_mesh, p_particles, p_e_bounds)]


# function for printing readable test labels
def id_fn(params):
    suffix = params[0].split('_')[-1]
    if suffix == 'n':
        return 'neutron-only'
    elif suffix == 'np':
        return 'neutron-photon'
    elif suffix == 'p':
        return 'photon-only'


@pytest.mark.parametrize('wwinp_data', expected_results, ids=id_fn)
def test_wwinp_reader(wwinp_data, request):
    wwinp_file, mesh, particle_types, energy_bounds = wwinp_data

    wws = openmc.wwinp_to_wws(request.node.path.parent / wwinp_file)

    for i, ww in enumerate(wws):
        e_bounds = energy_bounds[i]
        particle_type = particle_types[i]

        assert ww.particle_type == particle_type

        # check the mesh grid
        # there will be some very small changes due to the number of digits
        # provided in the wwinp format and the use of np.linspace to compute
        # boundaries of the fine mesh intervals
        np.testing.assert_allclose(mesh.x_grid, ww.mesh.x_grid, rtol=1e-6)
        np.testing.assert_allclose(mesh.y_grid, ww.mesh.y_grid, rtol=1e-6)
        np.testing.assert_allclose(mesh.z_grid, ww.mesh.z_grid, rtol=1e-6)

        # check the energy bounds
        np.testing.assert_array_equal(e_bounds, ww.energy_bounds)

        # check the expected weight window values mocked in the file --
        # a reversed array of the flat index into the numpy array
        n_wws = np.prod((*mesh.dimension, e_bounds.size - 1))
        exp_ww_lb = np.linspace(1, n_wws, n_wws)[::-1]
        np.testing.assert_array_equal(exp_ww_lb, ww.lower_ww_bounds.flatten())


# check expected failures
def fail_id_fn(params):
    suffix = params[0].split('_')[-1]
    if suffix == 't':
        return 'time-steps-failure'
    elif suffix == 'cyl':
        return 'cyl-mesh-failure'


expected_failure_data = (('wwinp_t', ValueError),
                         ('wwinp_cyl', NotImplementedError))


@pytest.mark.parametrize('wwinp_data', expected_failure_data, ids=fail_id_fn)
def test_wwinp_reader_failures(wwinp_data, request):
    filename, expected_failure = wwinp_data

    with pytest.raises(expected_failure):
        _ = openmc.wwinp_to_wws(request.node.path.parent / filename)
