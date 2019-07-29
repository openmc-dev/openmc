import os
import warnings
import itertools

from operator import itemgetter
from nose.tools import assert_equal, with_setup, assert_almost_equal, assert_raises
from random import uniform, seed

import numpy as np
from numpy.testing import assert_array_equal, assert_array_almost_equal

try:
    from pyne.mesh import Mesh
    # see if the source sampling module exists but do not import it
    import imp
    pyne_info = imp.find_module('pyne')
    pyne_mod = imp.load_module('pyne', *pyne_info)
    imp.find_module('source_sampling', pyne_mod.__path__)
except ImportError:
    from nose.plugins.skip import SkipTest
    raise SkipTest

from pyne.source_sampling import Sampler, AliasTable
from pyne.mesh import Mesh, NativeMeshTag
from pymoab import core as mb_core, types
from pyne.utils import QAWarning
warnings.simplefilter("ignore", QAWarning)


# Define modes
DEFAULT_ANALOG = 0
DEFAULT_UNIFORM = 1
DEFAULT_USER = 2
SUBVOXEL_ANALOG = 3
SUBVOXEL_UNIFORM = 4
SUBVOXEL_USER = 5


def try_rm_file(filename):
    return lambda: os.remove(filename) if os.path.exists(filename) else None


@with_setup(None, try_rm_file('sampling_mesh.h5m'))
def test_single_hex_tag_names_map():
    """This test tests uniform sampling within a single hex volume element.
    This is done by dividing the volume element in 4 smaller hex and ensuring
    that each sub-hex is sampled equally.
    """
    seed(1953)
    m = Mesh(structured=True,
             structured_coords=[[0, 3, 3.5], [0, 1], [0, 1]],
             mats=None)
    m.src = NativeMeshTag(2, float)
    m.src[:] = [[2.0, 1.0], [9.0, 3.0]]
    e_bounds = np.array([0, 0.5, 1.0])
    m.bias = NativeMeshTag(2, float)
    cell_fracs = np.zeros(2, dtype=[('idx', np.int64),
                                    ('cell', np.int64),
                                    ('vol_frac', np.float64),
                                    ('rel_error', np.float64)])
    cell_fracs[:] = [(0, 11, 1.0, 0.0), (1, 11, 1.0, 0.0)]
    m.tag_cell_fracs(cell_fracs)
    m.bias[:] = [[1.0, 2.0], [3.0, 3.0]]
    filename = "sampling_mesh.h5m"
    m.write_hdf5(filename)

    # right condition
    tag_names = {"src_tag_name": "src",
                 "cell_number_tag_name": "cell_number",
                 "cell_fracs_tag_name": "cell_fracs"}
    e_bounds = np.array([0, 1])
    sampler = Sampler(filename, tag_names, e_bounds, DEFAULT_ANALOG)

    # src_tag_name not given
    tag_names = {"cell_number_tag_name": "cell_number",
                 "cell_fracs_tag_name": "cell_fracs"}
    assert_raises(ValueError, Sampler, filename,
                  tag_names, e_bounds, DEFAULT_ANALOG)
    assert_raises(ValueError, Sampler, filename,
                  tag_names, e_bounds, DEFAULT_UNIFORM)
    assert_raises(ValueError, Sampler, filename,
                  tag_names, e_bounds, DEFAULT_USER)
    assert_raises(ValueError, Sampler, filename,
                  tag_names, e_bounds, SUBVOXEL_ANALOG)
    assert_raises(ValueError, Sampler, filename,
                  tag_names, e_bounds, SUBVOXEL_UNIFORM)
    assert_raises(ValueError, Sampler, filename,
                  tag_names, e_bounds, SUBVOXEL_USER)

    # bias_tag_name not given
    tag_names = {"src_tag_name": "src",
                 "cell_number_tag_name": "cell_number",
                 "cell_fracs_tag_name": "cell_fracs"}
    assert_raises(ValueError, Sampler, filename,
                  tag_names, e_bounds, DEFAULT_USER)
    assert_raises(ValueError, Sampler, filename,
                  tag_names, e_bounds, SUBVOXEL_USER)

    # cell_number_tag_name not given
    tag_names = {"src_tag_name": "src",
                 "cell_fracs_tag_name": "cell_fracs"}
    assert_raises(ValueError, Sampler, filename,
                  tag_names, e_bounds, DEFAULT_ANALOG)
    assert_raises(ValueError, Sampler, filename,
                  tag_names, e_bounds, DEFAULT_UNIFORM)
    assert_raises(ValueError, Sampler, filename,
                  tag_names, e_bounds, DEFAULT_USER)
    assert_raises(ValueError, Sampler, filename,
                  tag_names, e_bounds, SUBVOXEL_ANALOG)
    assert_raises(ValueError, Sampler, filename,
                  tag_names, e_bounds, SUBVOXEL_UNIFORM)
    assert_raises(ValueError, Sampler, filename,
                  tag_names, e_bounds, SUBVOXEL_USER)

    # cell_fracs_tag_name not given
    tag_names = {"src_tag_name": "src",
                 "cell_number_tag_name": "cell_number"}
    assert_raises(ValueError, Sampler, filename,
                  tag_names, e_bounds, DEFAULT_ANALOG)
    assert_raises(ValueError, Sampler, filename,
                  tag_names, e_bounds, DEFAULT_UNIFORM)
    assert_raises(ValueError, Sampler, filename,
                  tag_names, e_bounds, DEFAULT_USER)
    assert_raises(ValueError, Sampler, filename,
                  tag_names, e_bounds, SUBVOXEL_ANALOG)
    assert_raises(ValueError, Sampler, filename,
                  tag_names, e_bounds, SUBVOXEL_UNIFORM)
    assert_raises(ValueError, Sampler, filename,
                  tag_names, e_bounds, SUBVOXEL_USER)

    # wrong bias_tag data (non-zero source_density biased to zero -> NAN weight)
    m.src = NativeMeshTag(2, float)
    m.src[:] = [[1.0, 1.0]]
    m.bias = NativeMeshTag(2, float)
    m.bias[:] = [[0.0, 0.0]]
    m.write_hdf5(filename)
    tag_names = {"src_tag_name": "src",
                 "cell_number_tag_name": "cell_number",
                 "cell_fracs_tag_name": "cell_fracs",
                 "bias_tag_name": "bias"}
    assert_raises(RuntimeError, Sampler, filename,
                  tag_names, e_bounds, SUBVOXEL_USER)


@with_setup(None, try_rm_file('sampling_mesh.h5m'))
def test_analog_single_hex():
    """This test tests that particles of sampled evenly within the phase-space
    of a single mesh volume element with one energy group in an analog sampling
    scheme. This done by dividing each dimension (x, y, z, E) in half, then
    sampling particles and tallying on the basis of which of the 2^4 = 8 regions
    of phase space the particle is born into.
    """
    seed(1953)
    m = Mesh(structured=True, structured_coords=[[0, 1], [0, 1], [0, 1]],
             mats=None)
    m.src = NativeMeshTag(1, float)
    m.src[0] = 1.0
    cell_fracs = np.zeros(1, dtype=[('idx', np.int64),
                                    ('cell', np.int64),
                                    ('vol_frac', np.float64),
                                    ('rel_error', np.float64)])
    cell_fracs[:] = [(0, 11, 1.0, 0.0)]
    m.tag_cell_fracs(cell_fracs)

    filename = "sampling_mesh.h5m"
    m.write_hdf5(filename)
    tag_names = {"src_tag_name": "src",
                 "cell_number_tag_name": "cell_number",
                 "cell_fracs_tag_name": "cell_fracs"}
    sampler = Sampler(filename, tag_names, np.array([0, 1]), DEFAULT_ANALOG)

    num_samples = 5000
    score = 1.0/num_samples
    num_divs = 2
    tally = np.zeros(shape=(num_divs, num_divs, num_divs, num_divs))

    for i in range(num_samples):
        s = sampler.particle_birth(np.array([uniform(0, 1) for x in range(6)]))
        assert_equal(s.w, 1.0)  # analog: all weights must be one
        tally[int(s.x*num_divs), int(s.y*num_divs), int(s.z*num_divs),
              int(s.e*num_divs)] += score

    # Test that each half-space of phase space (e.g. x > 0.5) is sampled about
    # half the time.
    for i in range(0, 4):
        for j in range(0, 2):
            assert(abs(np.sum(np.rollaxis(tally, i)[j, :, :, :]) - 0.5) < 0.05)


@with_setup(None, try_rm_file('sampling_mesh.h5m'))
def test_analog_multiple_hex():
    """This test tests that particle are sampled uniformly from a uniform source
    defined on eight mesh volume elements in two energy groups. This is done
    using the exact same method ass test_analog_multiple_hex.
    """
    seed(1953)
    m = Mesh(structured=True,
             structured_coords=[[0, 0.5, 1], [0, 0.5, 1], [0, 0.5, 1]],
             mats=None)
    m.src = NativeMeshTag(2, float)
    m.src[:] = np.ones(shape=(8, 2))
    cell_fracs = np.zeros(8, dtype=[('idx', np.int64),
                                    ('cell', np.int64),
                                    ('vol_frac', np.float64),
                                    ('rel_error', np.float64)])
    cell_fracs[:] = [(0, 11, 1.0, 0.0),
                     (1, 11, 1.0, 0.0),
                     (2, 11, 1.0, 0.0),
                     (3, 11, 1.0, 0.0),
                     (4, 11, 1.0, 0.0),
                     (5, 11, 1.0, 0.0),
                     (6, 11, 1.0, 0.0),
                     (7, 11, 1.0, 0.0)]
    m.tag_cell_fracs(cell_fracs)
    filename = "sampling_mesh.h5m"
    m.write_hdf5(filename)
    tag_names = {"src_tag_name": "src",
                 "cell_number_tag_name": "cell_number",
                 "cell_fracs_tag_name": "cell_fracs"}
    sampler = Sampler(filename, tag_names, np.array(
        [0, 0.5, 1]), DEFAULT_ANALOG)

    num_samples = 5000
    score = 1.0/num_samples
    num_divs = 2
    tally = np.zeros(shape=(num_divs, num_divs, num_divs, num_divs))
    for i in range(num_samples):
        s = sampler.particle_birth([uniform(0, 1) for x in range(6)])
        assert_equal(s.w, 1.0)
        tally[int(s.x*num_divs), int(s.y*num_divs), int(s.z*num_divs),
              int(s.e*num_divs)] += score

    for i in range(0, 4):
        for j in range(0, 2):
            halfspace_sum = np.sum(np.rollaxis(tally, i)[j, :, :, :])
            assert(abs(halfspace_sum - 0.5)/0.5 < 0.1)


@with_setup(None, try_rm_file('tet.h5m'))
def test_analog_single_tet():
    """This test tests uniform sampling within a single tetrahedron. This is
    done by dividing the tetrahedron in 4 smaller tetrahedrons and ensuring
    that each sub-tet is sampled equally.
    """
    seed(1953)
    mesh = mb_core.Core()
    v1 = [0., 0., 0.]
    v2 = [1., 0., 0.]
    v3 = [0., 1., 0.]
    v4 = [0., 0., 1.]
    verts = mesh.create_vertices([v1, v2, v3, v4])
    mesh.create_element(types.MBTET, verts)
    m = Mesh(structured=False, mesh=mesh)
    m.src = NativeMeshTag(1, float)
    m.src[:] = np.array([1])
    filename = "tet.h5m"
    m.write_hdf5(filename)
    center = m.ve_center(list(m.iter_ve())[0])

    subtets = [[center, v1, v2, v3],
               [center, v1, v2, v4],
               [center, v1, v3, v4],
               [center, v2, v3, v4]]
    tag_names = {"src_tag_name": "src"}
    sampler = Sampler(filename, tag_names, np.array([0., 1.]), DEFAULT_ANALOG)
    num_samples = 5000
    score = 1.0/num_samples
    tally = np.zeros(shape=(4))
    for i in range(num_samples):
        s = sampler.particle_birth([uniform(0., 1.) for x in range(6)])
        assert_equal(s.w, 1.0)
        for i, tet in enumerate(subtets):
            if point_in_tet(tet, [s.x, s.y, s.z]):
                tally[i] += score
                break

    for t in tally:
        assert(abs(t - 0.25)/0.25 < 0.2)


@with_setup(None, try_rm_file('sampling_mesh.h5m'))
def test_uniform():
    """This test tests that the uniform biasing scheme:
    1. Samples space uniformly. This is checked using the same method
       described in test_analog_single_hex().
    2. Adjusts weights accordingly. Sample calculations are provided in Case 1
       in the Theory Manual.
    """
    seed(1953)
    m = Mesh(structured=True,
             structured_coords=[[0, 3, 3.5], [0., 1.], [0., 1.]],
             mats=None)
    m.src = NativeMeshTag(2, float)
    m.src[:] = [[2.0, 1.0], [9.0, 3.0]]
    e_bounds = np.array([0., 0.5, 1.0])
    filename = "sampling_mesh.h5m"
    cell_fracs = np.zeros(2, dtype=[('idx', np.int64),
                                    ('cell', np.int64),
                                    ('vol_frac', np.float64),
                                    ('rel_error', np.float64)])
    cell_fracs[:] = [(0, 11, 1.0, 0.0),
                     (1, 11, 1.0, 0.0)]
    m.tag_cell_fracs(cell_fracs)
    m.write_hdf5(filename)
    tag_names = {"src_tag_name": "src",
                 "cell_number_tag_name": "cell_number",
                 "cell_fracs_tag_name": "cell_fracs"}
    sampler = Sampler(filename, tag_names, e_bounds, DEFAULT_UNIFORM)

    num_samples = 10000
    score = 1.0/num_samples
    num_divs = 2
    num_e = 2
    spatial_tally = np.zeros(shape=(num_divs, num_divs, num_divs))
    e_tally = np.zeros(shape=(4))  # number of phase space groups
    for i in range(num_samples):
        s = sampler.particle_birth(np.array([uniform(0, 1) for x in range(6)]))
        if s.x < 3.0:
            assert_almost_equal(s.w, 0.7)  # hand calcs
        else:
            assert_almost_equal(s.w, 2.8)  # hand calcs

        spatial_tally[int(s.x*num_divs/3.5),
                      int(s.y*num_divs/1.0),
                      int(s.z*num_divs/1.0)] += score

        if s.x < 3 and s.e < 0.5:
            e_tally[0] += score
        elif s.x < 3 and s.e > 0.5:
            e_tally[1] += score
        if s.x > 3 and s.e < 0.5:
            e_tally[2] += score
        if s.x > 3 and s.e > 0.5:
            e_tally[3] += score

    for i in range(0, 3):
        for j in range(0, 2):
            halfspace_sum = np.sum(np.rollaxis(spatial_tally, i)[j, :, :])
            assert(abs(halfspace_sum - 0.5)/0.5 < 0.1)

    expected_e_tally = [4./7, 2./7, 3./28, 1./28]  # hand calcs
    for i in range(4):
        assert(abs(e_tally[i] - expected_e_tally[i])
               / expected_e_tally[i] < 0.1)


@with_setup(None, try_rm_file('sampling_mesh.h5m'))
def test_single_hex_single_subvoxel_analog():
    """This test tests that particles of sampled evenly within the phase-space
    of a single mesh volume element (also a sub-voxel) with one energy group
    in an analog sampling scheme. This done by dividing each dimension
    (x, y, z, E) in half, then sampling particles and tallying on the basis of
    which of the 2^4 = 16 regions of phase space the particle is born into.
    """
    seed(1953)
    m = Mesh(structured=True, structured_coords=[[0., 1.], [0., 1.], [0., 1.]],
             mats=None)
    m.src = NativeMeshTag(1, float)
    m.src[0] = 1.0
    cell_fracs = np.zeros(1, dtype=[('idx', np.int64),
                                    ('cell', np.int64),
                                    ('vol_frac', np.float64),
                                    ('rel_error', np.float64)])
    cell_fracs[:] = [(0, 11, 1.0, 0.0)]
    m.tag_cell_fracs(cell_fracs)
    filename = "sampling_mesh.h5m"
    m.write_hdf5(filename)
    tag_names = {"src_tag_name": "src",
                 "cell_number_tag_name": "cell_number",
                 "cell_fracs_tag_name": "cell_fracs"}
    sampler = Sampler(filename, tag_names, np.array([0., 1.]), SUBVOXEL_ANALOG)

    num_samples = 5000
    score = 1.0/num_samples
    num_divs = 2
    tally = np.zeros(shape=(num_divs, num_divs, num_divs, num_divs))

    for i in range(num_samples):
        s = sampler.particle_birth(np.array([uniform(0., 1.) for x in range(6)]))
        assert_equal(s.w, 1.0)  # analog: all weights must be one
        assert_equal(s.cell_list[0], 11)  # analog: the cell number
        tally[int(s.x*num_divs), int(s.y*num_divs), int(s.z*num_divs),
              int(s.e*num_divs)] += score

    # Test that each half-space of phase space (e.g. x > 0.5) is sampled about
    # half the time.
    for i in range(0, 4):
        for j in range(0, 2):
            assert(abs(np.sum(np.rollaxis(tally, i)[j, :, :, :]) - 0.5) < 0.05)


@with_setup(None, try_rm_file('sampling_mesh.h5m'))
def test_single_hex_multiple_subvoxel_analog():
    """This test tests that particles of sampled analog within the phase-space
    of a single mesh volume element but multiple sub-voxels with one energy
    group in an analog sampling scheme. Then sampling particles and tallying
    the particles and check the probability of particles born in each
    sub-voxel and the cell_number.
    """
    seed(1953)
    m = Mesh(structured=True, structured_coords=[[0., 1.], [0., 1.], [0., 1.]],
             mats=None)
    m.src = NativeMeshTag(3, float)
    m.src[:] = np.empty(shape=(1, 3))
    m.src[0] = [0, 0.2, 0.8]
    cell_fracs = np.zeros(3, dtype=[('idx', np.int64),
                                    ('cell', np.int64),
                                    ('vol_frac', np.float64),
                                    ('rel_error', np.float64)])
    cell_fracs[:] = [(0, 11, 0.4, 0.0), (0, 12, 0.3, 0.0), (0, 13, 0.3, 0.0)]
    m.tag_cell_fracs(cell_fracs) # cell_fracs will be sorted
    filename = "sampling_mesh.h5m"
    m.write_hdf5(filename)
    tag_names = {"src_tag_name": "src",
                 "cell_number_tag_name": "cell_number",
                 "cell_fracs_tag_name": "cell_fracs"}
    sampler = Sampler(filename, tag_names, np.array([0, 1]), SUBVOXEL_ANALOG)
    num_samples = 50000
    score = 1.0/num_samples
    num_divs = 2
    tally = [0.0] * 3
    for i in range(num_samples):
        s = sampler.particle_birth(np.array([uniform(0, 1) for x in range(6)]))
        assert_equal(s.w, 1.0)  # analog: all weights must be one
        if s.cell_list[0] == 11:
            tally[0] += score
        elif s.cell_list[0] == 12:
            tally[1] += score
        elif s.cell_list[0] == 13:
            tally[2] += score

    # Test that each source particle in each cell has right frequency
    assert_equal(tally[0], 0.0)
    assert(abs(tally[1] - 0.2)/0.2 < 0.05)
    assert(abs(tally[2] - 0.8)/0.8 < 0.05)


@with_setup(None, try_rm_file('sampling_mesh.h5m'))
def test_multiple_hex_multiple_subvoxel_analog():
    """This test tests that particle are sampled analog from a uniform source
    defined on eight mesh volume elements in two energy groups.
    """
    seed(1953)
    m = Mesh(structured=True,
             structured_coords=[[0, 0.5, 1], [0, 0.5, 1], [0, 0.5, 1]],
             mats=None)
    m.src = NativeMeshTag(2, float)
    m.src[:] = np.ones(shape=(8, 2))
    cell_fracs = np.zeros(8, dtype=[('idx', np.int64),
                                    ('cell', np.int64),
                                    ('vol_frac', np.float64),
                                    ('rel_error', np.float64)])
    cell_fracs[:] = [(0, 1, 1.0, 0.0), (1, 2, 1.0, 0.0), (2, 3, 1.0, 0.0),
                     (3, 4, 1.0, 0.0), (4, 5, 1.0, 0.0), (5, 6, 1.0, 0.0),
                     (6, 7, 1.0, 0.0), (7, 8, 1.0, 0.0)]
    m.tag_cell_fracs(cell_fracs)
    filename = "sampling_mesh.h5m"
    m.write_hdf5(filename)
    tag_names = {"src_tag_name": "src",
                 "cell_number_tag_name": "cell_number",
                 "cell_fracs_tag_name": "cell_fracs"}
    sampler = Sampler(filename, tag_names, np.array(
        [0, 0.5, 1]), SUBVOXEL_ANALOG)
    num_samples = 5000
    score = 1.0/num_samples
    num_divs = 2
    tally = np.zeros(shape=(num_divs, num_divs, num_divs, num_divs))
    for i in range(num_samples):
        s = sampler.particle_birth([uniform(0, 1) for x in range(6)])
        assert_equal(s.w, 1.0)
        assert_equal(s.cell_list[0], 4*int(s.x*num_divs) + 2*int(s.y*num_divs)
                     + int(s.z*num_divs) + 1)
        tally[int(s.x*num_divs), int(s.y*num_divs), int(s.z*num_divs),
              int(s.e*num_divs)] += score

    for i in range(0, 4):
        for j in range(0, 2):
            halfspace_sum = np.sum(np.rollaxis(tally, i)[j, :, :, :])
            assert(abs(halfspace_sum - 0.5)/0.5 < 0.1)


@with_setup(None, try_rm_file('sampling_mesh.h5m'))
def test_single_hex_subvoxel_uniform():
    """This test tests that particles of sampled evenly within the phase-space
    of a single mesh volume element with one energy group in an uniform sampling
    scheme. This done by dividing each dimension (x, y, z, E) in half, then
    sampling particles and tallying on the basis of which of the 2^4 = 8 regions
    of phase space the particle is born into.
    """
    seed(1953)
    m = Mesh(structured=True, structured_coords=[[0., 1.], [0., 1.], [0., 1.]],
             mats=None)
    m.src = NativeMeshTag(1, float)
    m.src[0] = 1.0
    cell_fracs = np.zeros(1, dtype=[('idx', np.int64),
                                    ('cell', np.int64),
                                    ('vol_frac', np.float64),
                                    ('rel_error', np.float64)])
    cell_fracs[:] = [(0, 11, 1.0, 0.0)]
    m.tag_cell_fracs(cell_fracs)
    filename = "sampling_mesh.h5m"
    m.write_hdf5(filename)
    tag_names = {"src_tag_name": "src",
                 "cell_number_tag_name": "cell_number",
                 "cell_fracs_tag_name": "cell_fracs"}
    sampler = Sampler(filename, tag_names, np.array([0., 1.]), SUBVOXEL_UNIFORM)

    num_samples = 5000
    score = 1.0/num_samples
    num_divs = 2
    tally = np.zeros(shape=(num_divs, num_divs, num_divs, num_divs))

    for i in range(num_samples):
        s = sampler.particle_birth(np.array([uniform(0., 1.) for x in range(6)]))
        assert_equal(s.w, 1.0)  # analog: all weights must be one
        assert_equal(s.cell_list[0], 11)  # analog: the cell number
        tally[int(s.x*num_divs), int(s.y*num_divs), int(s.z*num_divs),
              int(s.e*num_divs)] += score

     # Test that each half-space of phase space (e.g. x > 0.5) is sampled about
     # half the time.
    for i in range(0, 4):
        for j in range(0, 2):
            assert(abs(np.sum(np.rollaxis(tally, i)[j, :, :, :]) - 0.5) < 0.05)


@with_setup(None, try_rm_file('sampling_mesh.h5m'))
def test_single_hex_multiple_subvoxel_uniform():
    """This test tests that particles of sampled evenly within the phase-space
    of a single mesh volume element with one energy group in an uniform sampling
    scheme. This done by dividing each dimension (x, y, z, E) in half, then
    sampling particles and tallying on the basis of which of the 2^4 = 8 regions
    of phase space the particle is born into.
    """
    seed(1953)
    m = Mesh(structured=True, structured_coords=[[0., 1.], [0., 1.], [0., 1.]],
             mats=None)
    m.src = NativeMeshTag(3, float)
    m.src[:] = np.empty(shape=(1, 3))
    m.src[0] = [0, 0.2, 0.8]
    cell_fracs = np.zeros(3, dtype=[('idx', np.int64),
                                    ('cell', np.int64),
                                    ('vol_frac', np.float64),
                                    ('rel_error', np.float64)])
    cell_fracs[:] = [(0, 11, 0.4, 0.0), (0, 12, 0.3, 0.0), (0, 13, 0.3, 0.0)]
    m.tag_cell_fracs(cell_fracs)
    filename = "sampling_mesh.h5m"
    m.write_hdf5(filename)
    tag_names = {"src_tag_name": "src",
                 "cell_number_tag_name": "cell_number",
                 "cell_fracs_tag_name": "cell_fracs"}
    sampler = Sampler(filename, tag_names, np.array([0., 1.]), SUBVOXEL_UNIFORM)
    num_samples = 5000
    score = 1.0/num_samples
    num_divs = 2
    tally = [0.0] * 3
    for i in range(num_samples):
        s = sampler.particle_birth(np.array([uniform(0.0, 1.0) for x in range(6)]))
        if s.cell_list[0] == 11:
            tally[0] += score
        if s.cell_list[0] == 12:
            tally[1] += score
            # analog: all weights must be one
            assert(abs(s.w - 0.4)/0.4 < 0.05)
        if s.cell_list[0] == 13:
            tally[2] += score
            assert(abs(s.w - 1.6)/1.6 < 0.05)

    # Test that each source particle in each cell has right frequency
    assert_equal(tally[0], 0.0)
    assert(abs(tally[1] - 0.5) < 0.05)
    assert(abs(tally[2] - 0.5) < 0.05)


@with_setup(None, try_rm_file('sampling_mesh.h5m'))
def test_multiple_hex_multiple_subvoxel_uniform():
    """This test tests that particle are sampled uniformly from a uniform source
    defined on eight mesh volume elements in two energy groups.
    """
    seed(1953)
    m = Mesh(structured=True,
             structured_coords=[[0, 0.5, 1], [0, 0.5, 1], [0, 0.5, 1]],
             mats=None)
    m.src = NativeMeshTag(2, float)
    m.src[:] = np.empty(shape=(8, 2), dtype=float)
    m.src[:] = [[0, 0], [1, 0], [0, 0], [2, 0],
                [0, 0], [3, 0], [0, 0], [4, 0]]
    cell_fracs = np.zeros(8, dtype=[('idx', np.int64),
                                    ('cell', np.int64),
                                    ('vol_frac', np.float64),
                                    ('rel_error', np.float64)])
    cell_fracs[:] = [(0, 0, 1.0, 0.0), (1, 1, 1.0, 0.0), (2, 2, 1.0, 0.0),
                     (3, 3, 1.0, 0.0), (4, 4, 1.0, 0.0), (5, 5, 1.0, 0.0),
                     (6, 6, 1.0, 0.0), (7, 7, 1.0, 0.0)]
    empty_cells = [0, 2, 4, 6]
    m.tag_cell_fracs(cell_fracs)
    filename = "sampling_mesh.h5m"
    m.write_hdf5(filename)
    tag_names = {"src_tag_name": "src",
                 "cell_number_tag_name": "cell_number",
                 "cell_fracs_tag_name": "cell_fracs"}
    sampler = Sampler(filename, tag_names, np.array(
        [0, 0.5, 1]), SUBVOXEL_UNIFORM)
    num_samples = 50000
    score = 1.0/num_samples
    num_divs = 2
    tally = [0.0] * 8
    for i in range(num_samples):
        s = sampler.particle_birth([uniform(0, 1) for x in range(6)])
        # check the cell_number
        assert_equal(s.cell_list[0], 4*int(s.x*num_divs) + 2*int(s.y*num_divs)
                     + int(s.z*num_divs))
        # check the weight of each subvoxel
        if s.cell_list[0] not in empty_cells:
            # weight for cell 1, 3, 5, 7 should be: 0.4, 0.8, 1.2, 1.6
            exp_w = (s.cell_list[0] + 1) / 2 * 0.4
            out_w = s.w
            assert(abs(out_w - exp_w)/exp_w < 0.05)  # hand calculate
        # count the tally
        tally[s.cell_list[0]] += score

    # check the real sample rate
    for i, item in enumerate(tally):
        if i not in empty_cells:
            assert(abs(item - 0.25)/0.25 < 0.05)


@with_setup(None, try_rm_file('sampling_mesh.h5m'))
def test_bias():
    """This test tests that a user-specified biasing scheme:
    1. Samples space uniformly according to the scheme.
    2. Adjusts weights accordingly. Sample calculations are provided in Case 2
       in the Theory Manual.
    """
    seed(1953)
    m = Mesh(structured=True,
             structured_coords=[[0, 3, 3.5], [0, 1], [0, 1]],
             mats=None)
    m.src = NativeMeshTag(2, float)
    m.src[:] = [[2.0, 1.0], [9.0, 3.0]]
    e_bounds = np.array([0, 0.5, 1.0])
    m.bias = NativeMeshTag(2, float)
    m.bias[:] = [[1.0, 2.0], [3.0, 3.0]]
    cell_fracs = np.zeros(2, dtype=[('idx', np.int64),
                                    ('cell', np.int64),
                                    ('vol_frac', np.float64),
                                    ('rel_error', np.float64)])
    cell_fracs[:] = [(0, 11, 1.0, 0.0),
                     (1, 11, 1.0, 0.0)]
    m.tag_cell_fracs(cell_fracs)
    tag_names = {"src_tag_name": "src",
                 "cell_number_tag_name": "cell_number",
                 "cell_fracs_tag_name": "cell_fracs",
                 "bias_tag_name": "bias"}
    filename = "sampling_mesh.h5m"
    m.write_hdf5(filename)
    sampler = Sampler(filename, tag_names, e_bounds, DEFAULT_USER)

    num_samples = 10000
    score = 1.0/num_samples
    num_divs = 2
    tally = np.zeros(shape=(4))
    for i in range(num_samples):
        s = sampler.particle_birth(np.array([uniform(0, 1) for x in range(6)]))
        if s.x < 3:
            if s.e < 0.5:
                assert_almost_equal(s.w, 1.6)  # hand calcs
                tally[0] += score
            else:
                assert_almost_equal(s.w, 0.4)  # hand calcs
                tally[1] += score
        else:
            if s.e < 0.5:
                assert_almost_equal(s.w, 2.4)  # hand calcs
                tally[2] += score
            else:
                assert_almost_equal(s.w, 0.8)  # hand calcs
                tally[3] += score

    expected_tally = [0.25, 0.5, 0.125, 0.125]  # hand calcs
    for a, b in zip(tally, expected_tally):
        assert(abs(a-b)/b < 0.25)


@with_setup(None, try_rm_file('sampling_mesh.h5m'))
def test_bias_spatial():
    """This test tests a user-specified biasing scheme for which the only 1
    bias group is supplied for a source distribution containing two energy
    groups. This bias group is applied to both energy groups. In this test,
    the user-supplied bias distribution that was choosen, correspondes to
    uniform sampling, so that results can be checked against Case 1 in the
    theory manual.
    """
    seed(1953)
    m = Mesh(structured=True,
             structured_coords=[[0, 3, 3.5], [0, 1], [0, 1]],
             mats=None)
    m.src = NativeMeshTag(2, float)
    m.src[:] = [[2.0, 1.0], [9.0, 3.0]]
    m.bias = NativeMeshTag(1, float)
    m.bias[:] = [1, 1]
    e_bounds = np.array([0, 0.5, 1.0])
    filename = "sampling_mesh.h5m"
    cell_fracs = np.zeros(2, dtype=[('idx', np.int64),
                                    ('cell', np.int64),
                                    ('vol_frac', np.float64),
                                    ('rel_error', np.float64)])
    cell_fracs[:] = [(0, 11, 1.0, 0.0),
                     (1, 11, 1.0, 0.0)]
    m.tag_cell_fracs(cell_fracs)
    tag_names = {"src_tag_name": "src",
                 "cell_number_tag_name": "cell_number",
                 "cell_fracs_tag_name": "cell_fracs",
                 "bias_tag_name": "bias"}

    m.write_hdf5(filename)
    sampler = Sampler(filename, tag_names, e_bounds, DEFAULT_USER)

    num_samples = 10000
    score = 1.0/num_samples
    num_divs = 2
    num_e = 2
    spatial_tally = np.zeros(shape=(num_divs, num_divs, num_divs))
    e_tally = np.zeros(shape=(4))  # number of phase space groups
    for i in range(num_samples):
        s = sampler.particle_birth(np.array([uniform(0, 1) for x in range(6)]))
        if s.x < 3.0:
            assert_almost_equal(s.w, 0.7)  # hand calcs
        else:
            assert_almost_equal(s.w, 2.8)  # hand calcs

        spatial_tally[int(s.x*num_divs/3.5),
                      int(s.y*num_divs/1.0),
                      int(s.z*num_divs/1.0)] += score

        if s.x < 3 and s.e < 0.5:
            e_tally[0] += score
        elif s.x < 3 and s.e > 0.5:
            e_tally[1] += score
        if s.x > 3 and s.e < 0.5:
            e_tally[2] += score
        if s.x > 3 and s.e > 0.5:
            e_tally[3] += score

    for i in range(0, 3):
        for j in range(0, 2):
            halfspace_sum = np.sum(np.rollaxis(spatial_tally, i)[j, :, :])
            assert(abs(halfspace_sum - 0.5)/0.5 < 0.1)

    expected_e_tally = [4./7, 2./7, 3./28, 1./28]  # hand calcs
    for i in range(4):
        assert(abs(e_tally[i] - expected_e_tally[i])
               / expected_e_tally[i] < 0.1)


@with_setup(None, try_rm_file('sampling_mesh.h5m'))
def test_subvoxel_multiple_hex_bias_1():
    """This test tests that particle are sampled from a biased source
    defined on two voxels (2*2 = 4 sub-voxels) with the biased tag length of 1.
    """
    seed(1953)
    # mesh contains two voxels. 2 * 1 * 1 = 2
    m = Mesh(structured=True,
             structured_coords=[[0, 0.5, 1], [0, 1], [0, 1]],
             mats=None)

    # max_num_cells = 2. 4 sub-voxels
    cell_fracs = np.zeros(4, dtype=[('idx', np.int64),
                                    ('cell', np.int64),
                                    ('vol_frac', np.float64),
                                    ('rel_error', np.float64)])
    cell_fracs[:] = [(0, 11, 0.5, 0.0), (0, 12, 0.5, 0.0),
                     (1, 21, 0.5, 0.0), (1, 22, 0.5, 0.0)]
    m.tag_cell_fracs(cell_fracs)

    # the photon emitting rate of 4 sub-voxels is 0.1, 0.2, 0.3, 0.4
    m.src = NativeMeshTag(4, float)
    m.src[:] = np.empty(shape=(2, 4), dtype=float)
    m.src[:] = [[0.05, 0.05, 0.10, 0.10],
                [0.15, 0.15, 0.20, 0.20]]
    e_bounds = np.array([0, 0.5, 1.0])
    # bias, tag size = 1
    m.bias = NativeMeshTag(1, float)
    m.bias[:] = [[0.4], [0.6]]

    filename = "sampling_mesh.h5m"
    m.write_hdf5(filename)
    tag_names = {"src_tag_name": "src",
                 "cell_number_tag_name": "cell_number",
                 "cell_fracs_tag_name": "cell_fracs",
                 "bias_tag_name": "bias"}
    sampler = Sampler(filename, tag_names, e_bounds, SUBVOXEL_USER)

    num_samples = 50000
    score = 1.0/num_samples
    num_divs = 2
    # tally shape (v, c, e)
    tally = np.zeros(shape=(num_divs, num_divs, num_divs))
    for i in range(num_samples):
        s = sampler.particle_birth([uniform(0, 1) for x in range(6)])
        vid = s.cell_list[0]//10 - 1
        cid = s.cell_list[0] % 10 - 1
        eid = 0 if s.e < 0.5 else 1
        # check the cell_number
        if s.x < 0.5:
            assert(s.cell_list[0] in [11, 12])
        if s.x > 0.5:
            assert(s.cell_list[0] in [21, 22])
        # check the weight of each subvoxel
        if vid == 0:
            assert(abs(s.w - 0.746) / 0.746 < 0.05)
        if vid == 1:
            assert(abs(s.w - 1.163) / 1.163 < 0.05)
        # count the tally
        tally[vid, cid, eid] += score

    # check the real sample rate
    # exp_tally calculated by hand
    exp_tally = np.zeros(shape=(2, 2, 2))
    exp_tally[:] = [[[0.067, 0.067],
                     [0.133, 0.133]],
                    [[0.129, 0.129],
                     [0.171, 0.171]]]
    for v in range(2):
        for c in range(2):
            for e in range(2):
                assert(abs(tally[v, c, e] - exp_tally[v, c, e]
                           ) / exp_tally[v, c, e] < 0.05)


@with_setup(None, try_rm_file('sampling_mesh.h5m'))
def test_subvoxel_multiple_hex_bias_max_num_cells_num_e_groups():
    """This test tests that particle are sampled from a biased source
    defined on two voxels (2*2 = 4 sub-voxels) with the biased tag length
    of max_num_cells*num_e_group.
    """
    seed(1953)
    # mesh contains two voxels. 2 * 1 * 1 = 2
    m = Mesh(structured=True,
             structured_coords=[[0, 0.5, 1], [0, 1], [0, 1]],
             mats=None)

    # max_num_cells = 2. 4 sub-voxels
    cell_fracs = np.zeros(4, dtype=[('idx', np.int64),
                                    ('cell', np.int64),
                                    ('vol_frac', np.float64),
                                    ('rel_error', np.float64)])
    cell_fracs[:] = [(0, 11, 0.5, 0.0), (0, 12, 0.5, 0.0),
                     (1, 21, 0.5, 0.0), (1, 22, 0.5, 0.0)]
    m.tag_cell_fracs(cell_fracs)

    # the photon emitting rate of 4 sub-voxels is 0.1, 0.2, 0.3, 0.4
    m.src = NativeMeshTag(4, float)
    m.src[:] = np.empty(shape=(2, 4), dtype=float)
    m.src[:] = [[0.125, 0.125, 0.125, 0.125],
                [0.125, 0.125, 0.125, 0.125]]
    e_bounds = np.array([0, 0.5, 1.0])
    # bias, tag size = 1
    m.bias = NativeMeshTag(4, float)
    m.bias[:] = [[0.125, 0.125, 0.1, 0.15], [0.1, 0.1, 0.15, 0.15]]

    filename = "sampling_mesh.h5m"
    m.write_hdf5(filename)
    tag_names = {"src_tag_name": "src",
                 "cell_number_tag_name": "cell_number",
                 "cell_fracs_tag_name": "cell_fracs",
                 "bias_tag_name": "bias"}
    sampler = Sampler(filename, tag_names, e_bounds, SUBVOXEL_USER)

    num_samples = 50000
    score = 1.0/num_samples
    num_divs = 2
    # tally shape (v, c, e)
    tally = np.zeros(shape=(num_divs, num_divs, num_divs))
    exp_wgt = np.zeros(shape=(num_divs, num_divs, num_divs))
    exp_wgt[:] = [[[1.0, 1.0], [1.25, 0.83]], [[1.25, 1.25], [0.83, 0.83]]]
    for i in range(num_samples):
        s = sampler.particle_birth([uniform(0, 1) for x in range(6)])
        vid = s.cell_list[0]//10 - 1
        cid = s.cell_list[0] % 10 - 1
        eid = 0 if s.e < 0.5 else 1
        # check the cell_number
        if s.x < 0.5:
            assert(s.cell_list[0] in [11, 12])
        if s.x > 0.5:
            assert(s.cell_list[0] in [21, 22])
        # check the weight of each subvoxel
        assert(abs(s.w - exp_wgt[vid, cid, eid]) /
               exp_wgt[vid, cid, eid] < 0.05)
        # count the tally
        tally[vid, cid, eid] += score

    # check the real sample rate
    exp_tally = np.zeros(shape=(2, 2, 2))
    exp_tally[:] = [[[0.125, 0.125], [0.100, 0.150]],
                    [[0.100, 0.100], [0.150, 0.150]]]
    for v in range(2):
        for c in range(2):
            for e in range(2):
                assert(abs(tally[v, c, e] - exp_tally[v, c, e]
                           ) / exp_tally[v, c, e] < 0.05)


@with_setup(None, try_rm_file('sampling_mesh.h5m'))
def test_subvoxel_multiple_hex_bias_e_groups():
    """This test tests that particle are sampled from a biased source
    defined on two voxels (2*2 = 4 sub-voxels) with the biased tag length
    of energy groups.
    """
    seed(1953)
    # mesh contains two voxels. 2 * 1 * 1 = 2
    m = Mesh(structured=True,
             structured_coords=[[0, 0.5, 1], [0, 1], [0, 1]],
             mats=None)

    # max_num_cells = 2. 4 sub-voxels
    cell_fracs = np.zeros(4, dtype=[('idx', np.int64),
                                    ('cell', np.int64),
                                    ('vol_frac', np.float64),
                                    ('rel_error', np.float64)])
    cell_fracs[:] = [(0, 11, 0.5, 0.0), (0, 12, 0.5, 0.0),
                     (1, 21, 0.5, 0.0), (1, 22, 0.5, 0.0)]
    m.tag_cell_fracs(cell_fracs)

    # the photon emitting rate of 4 sub-voxels is 0.1, 0.2, 0.3, 0.4
    m.src = NativeMeshTag(4, float)
    m.src[:] = np.empty(shape=(2, 4), dtype=float)
    m.src[:] = [[0.05, 0.05, 0.10, 0.10],
                [0.15, 0.15, 0.20, 0.20]]
    e_bounds = np.array([0, 0.5, 1.0])
    # bias, tag size = 1
    m.bias = NativeMeshTag(2, float)
    m.bias[:] = [[0.1, 0.3], [0.2, 0.4]]

    filename = "sampling_mesh.h5m"
    m.write_hdf5(filename)
    tag_names = {"src_tag_name": "src",
                 "cell_number_tag_name": "cell_number",
                 "cell_fracs_tag_name": "cell_fracs",
                 "bias_tag_name": "bias"}
    sampler = Sampler(filename, tag_names, e_bounds, SUBVOXEL_USER)

    num_samples = 50000
    score = 1.0/num_samples
    num_divs = 2
    # tally shape (v, c, e)
    tally = np.zeros(shape=(num_divs, num_divs, num_divs))
    for i in range(num_samples):
        s = sampler.particle_birth([uniform(0, 1) for x in range(6)])
        vid = s.cell_list[0] // 10 - 1
        cid = s.cell_list[0] % 10 - 1
        eid = 0 if s.e < 0.5 else 1
        # check the cell_number
        if s.x < 0.5:
            assert(s.cell_list[0] in [11, 12])
        if s.x > 0.5:
            assert(s.cell_list[0] in [21, 22])
        # check the weight of each subvoxel
        if vid == 0 and eid == 0:
            assert(abs(s.w - 1.5)/1.5 < 0.05)
        if vid == 0 and eid == 1:
            assert(abs(s.w - 0.5)/0.5 < 0.05)
        if vid == 1 and eid == 0:
            assert(abs(s.w - 1.75)/1.75 < 0.05)
        if vid == 1 and eid == 1:
            assert(abs(s.w - 0.875)/0.875 < 0.05)
        # count the tally
        tally[vid, cid, eid] += score

    # check the real sample rate
    exp_tally = np.zeros(shape=(2, 2, 2))
    exp_tally[:] = [[[0.0333, 0.1000],
                     [0.0667, 0.2000]],
                    [[0.0857, 0.1714],
                     [0.1143, 0.2286]]]
    for v in range(2):
        for c in range(2):
            for e in range(2):
                assert(abs(tally[v, c, e] - exp_tally[v, c, e]
                           ) / exp_tally[v, c, e] < 0.05)


def test_alias_table():
    """This tests that the AliasTable class produces samples in the ratios
    consistant with the supplied PDF.
    """
    seed(1953)
    pdf = np.array([0.1, 0.2, 0.7])
    at = AliasTable(pdf)
    num_samples = 50000
    score = 1.0/num_samples
    tally = np.zeros(shape=(3))

    for i in range(num_samples):
        s = at.sample_pdf(uniform(0, 1), uniform(0, 1))
        tally[s] += score

    for i in range(0, 3):
        assert(abs(tally[i] - pdf[i])/pdf[i] < 0.05)


def point_in_tet(t, p):
    """ This function determines if some point <p> lies within some tetrahedron
    <t> using the method described here:
    http://steve.hollasch.net/cgindex/geometry/ptintet.html
    """
    matricies = [
        np.array([[t[0][0], t[0][1], t[0][2], 1],
                  [t[1][0], t[1][1], t[1][2], 1],
                  [t[2][0], t[2][1], t[2][2], 1],
                  [t[3][0], t[3][1], t[3][2], 1]]),
        np.array([[p[0], p[1], p[2], 1],
                  [t[1][0], t[1][1], t[1][2], 1],
                  [t[2][0], t[2][1], t[2][2], 1],
                  [t[3][0], t[3][1], t[3][2], 1]]),
        np.array([[t[0][0], t[0][1], t[0][2], 1],
                  [p[0], p[1], p[2], 1],
                  [t[2][0], t[2][1], t[2][2], 1],
                  [t[3][0], t[3][1], t[3][2], 1]]),
        np.array([[t[0][0], t[0][1], t[0][2], 1],
                  [t[1][0], t[1][1], t[1][2], 1],
                  [p[0], p[1], p[2], 1],
                  [t[3][0], t[3][1], t[3][2], 1]]),
        np.array([[t[0][0], t[0][1], t[0][2], 1],
                  [t[1][0], t[1][1], t[1][2], 1],
                  [t[2][0], t[2][1], t[2][2], 1],
                  [p[0], p[1], p[2], 1]])]

    determinates = [np.linalg.det(x) for x in matricies]
    return all(x >= 0 for x in determinates) or all(x < 0 for x in determinates)


def test_template_examples():
    """
    An example of using source_sampling test template to do the test
    """
    # DEFAULT and SUBVOXEL
    for mode in (DEFAULT_ANALOG, DEFAULT_UNIFORM, DEFAULT_USER,
                 SUBVOXEL_ANALOG, SUBVOXEL_UNIFORM, SUBVOXEL_USER):
        for num_e_groups in (1, 2):
            # num_bias_groups could be:
            # 1, num_e_groups, and max_num_cells*num_e_groups

            # test case: 1 voxel, 1 subvoxel
            cell_fracs_list = [(0, 1, 1.0, 0.0)]
            src_tag = [[1.0]*num_e_groups]
            if mode == DEFAULT_USER or mode == SUBVOXEL_USER:
                for num_bias_groups in (1, num_e_groups):
                    bias_tag = [[1.0]*num_bias_groups]
                    _source_sampling_test_template(
                        mode, cell_fracs_list, src_tag, bias_tag)
            else:
                _source_sampling_test_template(mode, cell_fracs_list, src_tag)

            # test case: 1 voxel, 2 subvoxels
            # create src and cell_fracs tag data
            if mode in (0, 1, 2):
                src_tag = [[1.0]*num_e_groups]
                cell_fracs_list = [(0, 1, 1.0, 0.0)]
            elif mode in (3, 4, 5):
                src_tag = [[1.0, 1.0]*num_e_groups]
                cell_fracs_list = [(0, 1, 0.5, 0.0), (0, 2, 0.5, 0.0)]

            if mode == DEFAULT_USER:
                for num_bias_groups in (1, num_e_groups):
                    bias_tag = [[1.0]*num_bias_groups]
                    _source_sampling_test_template(
                        mode, cell_fracs_list, src_tag, bias_tag)
            elif mode == SUBVOXEL_USER:
                for num_bias_groups in (1, num_e_groups, 2*num_e_groups):
                    bias_tag = [[1.0]*num_bias_groups]
                    _source_sampling_test_template(
                        mode, cell_fracs_list, src_tag, bias_tag)
            else:
                _source_sampling_test_template(mode, cell_fracs_list, src_tag)

            # test case: 2 voxel, 2 subvoxels
            cell_fracs_list = [(0, 1, 1.0, 0.0), (1, 2, 1.0, 0.0)]
            src_tag = [[1.0]*num_e_groups, [1.0]*num_e_groups]
            if mode == DEFAULT_USER or mode == SUBVOXEL_USER:
                for num_bias_groups in (1, num_e_groups):
                    bias_tag = [[1.0]*num_bias_groups, [1.0]*num_bias_groups]
                    _source_sampling_test_template(
                        mode, cell_fracs_list, src_tag, bias_tag)
            else:
                _source_sampling_test_template(mode, cell_fracs_list, src_tag)

            # test case: 2 voxel, 4 subvoxels
            # create src and cell_fracs tag data
            if mode in (0, 1, 2):
                src_tag = [[1.0]*num_e_groups, [1.0]*num_e_groups]
                cell_fracs_list = [(0, 1, 1.0, 0.0),
                                   (1, 2, 1.0, 0.0)]
            elif mode in (3, 4, 5):
                src_tag = [[1.0, 1.0]*num_e_groups, [1.0, 1.0]*num_e_groups]
                cell_fracs_list = [(0, 1, 0.5, 0.0), (0, 2, 0.5, 0.0),
                                   (1, 3, 0.5, 0.0), (1, 4, 0.5, 0.0)]

            if mode == DEFAULT_USER:
                for num_bias_groups in (1, num_e_groups):
                    bias_tag = [[1.0]*num_bias_groups, [1.0]*num_bias_groups]
                    _source_sampling_test_template(
                        mode, cell_fracs_list, src_tag, bias_tag)
            elif mode == SUBVOXEL_USER:
                for num_bias_groups in (1, num_e_groups, 2*num_e_groups):
                    bias_tag = [[1.0]*num_bias_groups, [1.0]*num_bias_groups]
                    _source_sampling_test_template(
                        mode, cell_fracs_list, src_tag, bias_tag)
            else:
                _source_sampling_test_template(mode, cell_fracs_list, src_tag)


def _get_num_ve_sve_and_max_num_cells(cell_fracs):
    """
    Calculate the num_ve, num_sve and max_num_cells

    Parameters
    ----------
    cell_fracs : structured array, optional
        A sorted, one dimensional array,
        each entry containing the following fields:

            :idx: int
                The volume element index.
            :cell: int
                The geometry cell number.
            :vol_frac: float
                The volume fraction of the cell withing the mesh ve.
            :rel_error: float
                The relative error associated with the volume fraction.

    Returns
    -------
    num_ve : int
        Number of the total voxels
    num_sve : int
        Number of the total subvoxels, eqaul to or greater than num_ve
    max_num_cells : int
        Max number of cells (subvoxels) in a voxel
    """
    num_sve = len(cell_fracs)
    num_ve = len(set(cell_fracs['idx']))
    max_num_cells = -1
    for i in range(num_sve):
        max_num_cells = max(max_num_cells,
                            len(cell_fracs[cell_fracs['idx'] == i]))
    return num_ve, num_sve, max_num_cells


def _create_mesh_via_num_ve(num_ve):
    """
    This function creates mesh from number of voxels

    Parameters
    ----------
    num_ve : int
        Number of voxels

    Returns
    -------
    mesh. MOAB mesh.
    """
    x_bounds = [v*1.0/(num_ve) for v in range(num_ve+1)]
    mesh = Mesh(structured=True,
                structured_coords=[x_bounds, [0, 1], [0, 1]],
                mats=None)
    return mesh


def _cal_pdf_and_biased_pdf(cell_fracs, src_tag, bias_tag=None):
    """
    This function calcualtes the normalized pdf of source.

    Parameters
    ----------
    cell_fracs : structured array
        A sorted, one dimensional array,
        each entry containing the following fields:

            :idx: int
                The volume element index.
            :cell: int
                The geometry cell number.
            :vol_frac: float
                The volume fraction of the cell withing the mesh ve.
            :rel_error: float
                The relative error associated with the volume fraction.
    src_tag : numpy array
        An one or two dimentional array contains data of the source tag.
    bias_tag : numpy array, optional
        An one or two dimentional array contains data of bias tag

    Returns
    -------
    pdf : numpy array
        A three dimentional numpy array, shape=(num_ve, num_sve, num_e_groups)
    biased_pdf : numpy array
        A three dimentional numpy array, shape=(num_ve, num_sve, num_e_groups)
    """
    num_ve, num_sve, max_num_cells = _get_num_ve_sve_and_max_num_cells(
        cell_fracs)
    num_e_groups = len(src_tag[0])//max_num_cells
    pdf = np.empty(shape=(num_ve, max_num_cells, num_e_groups),
                   dtype=np.float64)
    pdf.fill(0.0)
    for vid in range(num_ve):
        for svid in range(max_num_cells):
            for eid in range(num_e_groups):
                pdf[vid, svid, eid] = src_tag[vid][svid*num_e_groups + eid] \
                    * cell_fracs[vid*max_num_cells + svid]['vol_frac']
    # normalize
    pdf = pdf/pdf.sum()

    # calculate biased_pdf
    biased_pdf = np.empty(shape=(num_ve, max_num_cells, num_e_groups),
                          dtype=np.float64)
    biased_pdf.fill(0.0)
    # set up bias_array to proper value
    if bias_tag == None:
        # UNIFORM mode, set default bias_group and bias_array
        num_bias_groups = 1
        bias_array = np.empty(shape=(num_ve, max_num_cells, num_e_groups),
                              dtype=np.float64)
        for vid in range(num_ve):
            for svid in range(max_num_cells):
                for eid in range(num_e_groups):
                    bias_array[vid, svid, eid] = src_tag[vid][svid*num_e_groups+eid] /\
                        np.array(src_tag[vid]).sum()
    else:
        # USER mode, set bias_array according to bias_tag
        num_bias_groups = len(bias_tag[0])
        bias_array = np.empty(shape=(num_ve, max_num_cells, num_e_groups),
                              dtype=np.float64)
        bias_array.fill(0.0)
        for vid in range(num_ve):
            for svid in range(max_num_cells):
                for eid in range(num_e_groups):
                    if num_bias_groups == 1:
                        bias_array[vid, svid, eid] = bias_tag[vid][0]
                    elif num_bias_groups == num_e_groups:
                        bias_array[vid, svid, eid] = bias_tag[vid][eid]
                    elif num_bias_groups == max_num_cells*num_e_groups:
                        bias_array[vid, svid,
                                   eid] = bias_tag[vid][svid*num_e_groups + eid]
                    else:
                        raise ValueError("Wrong bias_tag length")
    # calculate biased_pdf
    if num_bias_groups == 1:
        for vid in range(num_ve):
            for svid in range(max_num_cells):
                current_ve = cell_fracs[cell_fracs['idx'] == vid]
                biased_pdf[vid, svid, :] = bias_array[vid, svid, :] \
                    * current_ve[svid]['vol_frac']
    elif num_bias_groups == num_e_groups:
        for vid in range(num_ve):
            for eid in range(num_e_groups):
                for svid in range(max_num_cells):
                    current_ve = cell_fracs[cell_fracs['idx'] == vid]
                    biased_pdf[vid, svid, eid] = bias_array[vid, svid, eid] \
                        * current_ve[svid]['vol_frac']
    elif num_bias_groups == max_num_cells*num_e_groups:
        for vid in range(num_ve):
            for svid in range(max_num_cells):
                for eid in range(num_e_groups):
                    biased_pdf[vid, svid, eid] = bias_array[vid, svid, eid] \
                        * cell_fracs[vid]['vol_frac']
    # normalize biased_pdf
    biased_pdf = np.divide(biased_pdf, biased_pdf.sum())
    return pdf, biased_pdf


def _cal_exp_w_c(s, mode, cell_fracs, src_tag, bias_tag):
    """
    This function calcualtes the exptected weight and cell_number
    for a given particle (according to it's x coordinate)

    Parameters
    ----------
    s : SourceParticle
        The given particle
    mode : int
        Mode of the source_sampling
    cell_fracs : structured array
        A sorted, one dimensional array,
        each entry containing the following fields:

            :idx: int
                The volume element index.
            :cell: int
                The geometry cell number.
            :vol_frac: float
                The volume fraction of the cell withing the mesh ve.
            :rel_error: float
                The relative error associated with the volume fraction.
    src_tag : numpy array
        An one or two dimentional array contains data of the source tag.
    bias_tag : numpy array, optional
        An one or two dimentional array contains data of bias tag

    Returns
    -------
    exp_w : float
        Expected weight of the source particle
    exp_c : set of available cell numbers
        Expected cell number of the source particle
    """
    num_ve, num_sve, max_num_cells = \
        _get_num_ve_sve_and_max_num_cells(cell_fracs)
    # calculate vid
    x_bounds = [v*1.0/(num_ve) for v in range(num_ve+1)]
    vid = -1
    for i in range(num_ve):
        if x_bounds[i] <= s.x <= x_bounds[i+1]:
            vid = i
            break
    if vid == -1:
        raise ValueError("x coordinate of particle not in (0, 1), s.x = {0}"
                         .format(str(s.x)))
    # calculate svid
    # get number of cells/subvoxels of current voxel
    current_cell_fracs = cell_fracs[cell_fracs['idx'] == vid]
    num_cells = len(current_cell_fracs)
    x_bounds = np.array([0.0]*(num_cells+1))
    # the x_bounds of the vid start from 1.0/num_ve*vid
    x_bounds[0] = 1.0/num_ve*vid
    for svid in range(num_cells):
        x_bounds[svid+1] = x_bounds[svid] + 1.0/num_ve *\
            current_cell_fracs[svid]['vol_frac']
    svid = -1
    for i in range(num_cells):
        if x_bounds[i] <= s.x <= x_bounds[i+1]:
            svid = i
            break
    if svid == -1:
        raise ValueError("x coordinate not in the voxel, s.x = {0}"
                         .format(str(s.x)))
    # get the cell_number
    exp_c = set(list(current_cell_fracs['cell']))

    # calculate eid
    if mode in (0, 1, 2):
        num_e_groups = len(src_tag[0])
    elif mode in (3, 4, 5):
        num_e_groups = len(src_tag[0])//max_num_cells
    e_bounds = np.array([i*1.0/num_e_groups for i in range(num_e_groups+1)])
    eid = -1
    for i in range(num_e_groups):
        if e_bounds[i] <= s.e <= e_bounds[i+1]:
            eid = i
            break
    if eid == -1:
        raise ValueError("energy not in (0, 1), s.e = ".format(str(s.e)))
    # calculate exp_w, weight is determined by mode, vid, svid and energy
    if mode in (0, 3):
        # ANALOG
        exp_w = 1.0
    elif mode in (1, 4):
        # UNIFORM
        pdf, biased_pdf = _cal_pdf_and_biased_pdf(
            cell_fracs, src_tag, bias_tag)
        exp_w = pdf[vid, svid, eid] / biased_pdf[vid, svid, eid]
    else:
        # USER
        pdf, biased_pdf = _cal_pdf_and_biased_pdf(
            cell_fracs, src_tag, bias_tag)
        exp_w = pdf[vid, svid, eid] / biased_pdf[vid, svid, eid]
    return exp_w, exp_c


def _get_p_y_z_halfspace(particles):
    """
    This function calcualtes the probabilities of y and z half space
    for a given set of particles

    Parameters
    ----------
    particles : list
        List of SourceParticle

    Returns
    -------
    p_y_halfspace : float
        The probability of y half space
    p_z_halfspace : float
        The probability of z half space
    """
    y_count, z_count = 0, 0
    for s in particles:
        if s.y < 0.5:
            y_count = y_count + 1
        if s.z < 0.5:
            z_count = z_count + 1
    p_y_halfspace = float(y_count)/len(particles)
    p_z_halfspace = float(z_count)/len(particles)
    return p_y_halfspace, p_z_halfspace


def _get_x_dis(particles, num_ve):
    """
    This function calcualtes the particle distribution along x direction
    for a given set of particles

    Parameters
    ----------
    particles : list
        List of SourceParticle
    num_ve : int
        Number of voxels

    Returns
    -------
    x_dis : one dimentional numpy array
        The particle direction along x direction
    """
    x_bounds = [v*1.0/(num_ve) for v in range(num_ve+1)]
    x_dis = np.array([0.0]*num_ve)
    for i in range(num_ve):
        for s in particles:
            if x_bounds[i] <= s.x <= x_bounds[i+1]:
                x_dis[i] = x_dis[i] + 1
    x_dis = np.divide(x_dis, len(particles))
    return x_dis


def _get_x_dis_exp(mode, cell_fracs, src_tag, bias_tag=None):
    """
    This function calcualtes the exptected particle distribution along x
    direction.

    Parameters
    ----------
    mode : int
        Mode of the source_sampling
    cell_fracs : structured array
        A sorted, one dimensional array,
        each entry containing the following fields:

            :idx: int
                The volume element index.
            :cell: int
                The geometry cell number.
            :vol_frac: float
                The volume fraction of the cell withing the mesh ve.
            :rel_error: float
                The relative error associated with the volume fraction.
    src_tag : numpy array
        An one or two dimentional array contains data of the source tag.
    bias_tag : numpy array, optional
        An one or two dimentional array contains data of bias tag

    Returns
    -------
    x_dis_exp : one dimentional numpy array
        The expected particle direction along x direction
    """
    num_ve, num_sve, max_num_cells = \
        _get_num_ve_sve_and_max_num_cells(cell_fracs)
    if mode in (0, 1, 2):
        num_e_groups = len(src_tag[0])
    elif mode in (3, 4, 5):
        num_e_groups = len(src_tag[0])//max_num_cells
    x_bounds = [v*1.0/(num_ve) for v in range(num_ve+1)]
    x_dis_exp = np.array([0.0]*num_ve)
    if mode in (0, 3):
        # ANALOG, particles distribution according to the src_tag
        for vid in range(num_ve):
            current_ve = cell_fracs[cell_fracs['idx'] == vid]
            for svid in range(len(current_ve)):
                x_dis_exp[vid] += current_ve[svid]['vol_frac'] *\
                    np.array(
                        src_tag[vid][svid*num_e_groups:(svid+1)*num_e_groups]).sum()
    elif mode in (1, 4):
        # UNIFORM, particles distribution uniformly in x direction
        x_dis_exp = np.array([1.0/num_ve]*num_ve)
    elif mode in (2, 5):
        if bias_tag == None:
            raise ValueError("bias_tag must be provided when mode is {0}"
                             .format(str(mode)))
        # USER, particles distribute accroding to the bias_tag
        for vid in range(num_ve):
            current_ve = cell_fracs[cell_fracs['idx'] == vid]
            for svid in range(len(current_ve)):
                x_dis_exp[vid] += current_ve[svid]['vol_frac'] *\
                    np.array(
                        bias_tag[vid][svid*num_e_groups:(svid+1)*num_e_groups]).sum()
    # normalize x_dis_exp
    x_dis_exp = np.divide(x_dis_exp, x_dis_exp.sum())
    return x_dis_exp


def _get_e_dis(particles, num_e_groups):
    """
    This function calcualtes the particle distribution along energy
    for a given set of particles

    Parameters
    ----------
    particles : list
        List of SourceParticle
    num_e_groups : int
        Number of energy groups

    Returns
    -------
    e_dis : one dimentional numpy array
        The particle direction along energy
    """
    e_bounds = [e*1.0/(num_e_groups) for e in range(num_e_groups+1)]
    e_dis = np.array([0.0]*num_e_groups)
    for i in range(num_e_groups):
        for s in particles:
            if e_bounds[i] <= s.e <= e_bounds[i+1]:
                e_dis[i] = e_dis[i] + 1
    e_dis = np.divide(e_dis, len(particles))
    return e_dis


def _get_e_dis_exp(mode, cell_fracs, src_tag, bias_tag=None):
    """
    This function calcualtes the exptected particle distribution along energy

    Parameters
    ----------
    mode : int
        Mode of the source_sampling
    cell_fracs : structured array
        A sorted, one dimensional array,
        each entry containing the following fields:

            :idx: int
                The volume element index.
            :cell: int
                The geometry cell number.
            :vol_frac: float
                The volume fraction of the cell withing the mesh ve.
            :rel_error: float
                The relative error associated with the volume fraction.
    src_tag : numpy array
        An one or two dimentional array contains data of the source tag.
    bias_tag : numpy array, optional
        An one or two dimentional array contains data of bias tag

    Returns
    -------
    e_dis_exp : one dimentional numpy array
        The expected particle direction along energy
    """
    # input check
    if mode in (2, 5) and bias_tag == None:
        raise ValueError("bias_tag must be provided when mode is {0}"
                         .format(str(mode)))
    num_ve, num_sve, max_num_cells = \
        _get_num_ve_sve_and_max_num_cells(cell_fracs)
    if mode in (0, 1, 2):
        num_e_groups = len(src_tag[0])
    elif mode in (3, 4, 5):
        num_e_groups = len(src_tag[0])//max_num_cells
    e_bounds = [e*1.0/(num_e_groups) for e in range(num_e_groups+1)]
    e_dis_exp = np.array([0.0]*num_e_groups)
    if mode in (0, 1, 3, 4) or (mode in (2, 5) and len(bias_tag[0]) == 1):
        # when mode is ANALOG and UNIFORM, or mode is USER but num_bias_groups is 1
        # particles distribution according to the src_tag
        for vid in range(num_ve):
            current_ve = cell_fracs[cell_fracs['idx'] == vid]
            for svid in range(len(current_ve)):
                for eid in range(num_e_groups):
                    e_dis_exp[eid] += current_ve[svid]['vol_frac'] *\
                        src_tag[vid][svid*num_e_groups+eid]
    elif mode == 2 or (mode == 5 and len(bias_tag[0]) == num_e_groups):
        # Energy is biased according to the bias_tag
        for vid in range(num_ve):
            current_ve = cell_fracs[cell_fracs['idx'] == vid]
            for svid in range(len(current_ve)):
                for eid in range(num_e_groups):
                    e_dis_exp[eid] += current_ve[svid]['vol_frac'] *\
                        bias_tag[vid][eid]
    else:
        for vid in range(num_ve):
            current_ve = cell_fracs[cell_fracs['idx'] == vid]
            for svid in range(len(current_ve)):
                for eid in range(num_e_groups):
                    e_dis_exp[eid] += current_ve[svid]['vol_frac'] *\
                        bias_tag[vid][svid*num_e_groups+eid]
    # normalize x_dis_exp
    e_dis_exp = np.divide(e_dis_exp, e_dis_exp.sum())
    return e_dis_exp


@with_setup(None, try_rm_file('sampling_mesh.h5m'))
def _source_sampling_test_template(mode, cell_fracs_list, src_tag,
                                   bias_tag=None):
    """
    This function serve as a template for all source_sampling test cases.
    It constrcut Sampler from input parameters.
    And then perform a standardized sampling and tally,
    Finally, it compares tallied results with exp_answers.

    Assumptions:
        * Use unit cube for all the meshes
        * Use structured meshes for all the tests
        * filename will always be: "sampling_mesh.h5m"
        * distribution changes only on X direction
        * uniform distribution in Y and Z directions
        * cell_number always equal to the index of sve + 1, no void cell
        * voxels have the same volume
        For example:
        cell_fracs = [(0, 1, 0.4, 0.0),
                      (0, 2, 0.6, 0.0),
                      (1, 3, 1.0, 0.0),
                      (2, 4, 1.0, 0.0), ...]

        voxel idx           v0           v1           v2
                      |------------|------------|------------|---       y
                      |     |      |            |            |          ^  z
        subvoxel      | sve0| sve1 |    sve2    |    sve3    | . . .    | /
                      |     |      |            |            |          |/
                      |------------|------------|------------|---       ----> x
        cell_number      c1    c2        c3           c4


        * Energy have only two options:
            - [0.0, 1.0]
            - [0.0, 0.5, 1.0]
        * Voxel number of meshes could be:
            - 1 voxel 1 subvoxel -> Single voxel single subvoxel
            - 1 voxel 2 subvoxel -> Single voxel multiple subvoxel
            - 2 voxel 2 subvoxel -> Multiple voxel multiple subvoxel
            - 2 voxel 4 subvoxel -> Multiple voxel multiple subvoxel

    Under these assumptions:
        * Mesh could be derived from cell_fracs
        * e_bounds could be derived from src_tag
        * construct_paras contain:
            - mode
            - cell_fracs
            - src_tag
            - bias_tag (optional, required for bias_mode == USER)

    Check items:
        * weight for each particle
        * cell_number for each particle
        * position distribution
        * energy distribution

    Parameters
    ----------
    mode : int
        Mode of the source sampling, could be 0, 1, 2, 3, 4 or 5
    cell_fracs_list : numpy array
        A one dimentional numpy array used to construct cell_fracs,
        Element: (idx, cell, vol_frac, rel_error)
    src_tag : numpy array
        An one or two dimentional array contains data of the source tag.
    bias_tag : numpy array, optional
        An one or two dimentional array contains data of bias tag

    Returns
    -------
    None
    """
    sub_mode_r2s = (0, 1, 2)
    sub_mode_subvoxel = (3, 4, 5)
    avail_mode = (0, 1, 2, 3, 4, 5)
    # input check
    # check mode
    if mode not in avail_mode:
        raise ValueError("mode must be in (0, 1, 2, 3, 4, 5)")
    # set cell_fracs
    cell_fracs = np.zeros(len(cell_fracs_list),
                          dtype=[('idx', np.int64),
                                 ('cell', np.int64),
                                 ('vol_frac', np.float64),
                                 ('rel_error', np.float64)])
    cell_fracs[:] = cell_fracs_list
    # check bias_tag
    if mode in (2, 5) and bias_tag == None:  # bias_mode == USER
        raise ValueError(
            "bias_tag must be given when mode is {0}".format(str(mode)))

    # get number of voxel, max_num_cells
    num_ve, num_sve, max_num_cells = _get_num_ve_sve_and_max_num_cells(
        cell_fracs)
    # set up e_bounds
    if mode in (0, 1, 2):
        num_e_groups = len(src_tag[0])
    elif mode in (3, 4, 5):
        num_e_groups = len(src_tag[0])//max_num_cells
    e_bounds = [i*1.0/num_e_groups for i in range(num_e_groups+1)]
    e_bounds = np.array(e_bounds)
    # set up mesh
    m = _create_mesh_via_num_ve(num_ve)
    # set up src tag
    if mode in (0, 1, 2):
        m.src = NativeMeshTag(num_e_groups, float)
    elif mode in (3, 4, 5):
        m.src = NativeMeshTag(max_num_cells*num_e_groups, float)
    m.src[:] = src_tag
    # set up cell_number and cell_fracs tag
    m.tag_cell_fracs(cell_fracs)
    # set up bias tag
    if mode in (2, 5):
        bias_tag_lenght = len(bias_tag[0])
        m.bias = NativeMeshTag(bias_tag_lenght, float)
        m.bias[:] = bias_tag
    # set up tag_names
    tag_names = {"src_tag_name": "src",
                 "cell_number_tag_name": "cell_number",
                 "cell_fracs_tag_name": "cell_fracs"}
    if mode in (2, 5):
        tag_names["bias_tag_name"] = "bias"
    # save the mesh into h5m file
    filename = "sampling_mesh.h5m"
    m.write_hdf5(filename)

    # construct Sampler
    sampler = Sampler(filename, tag_names, e_bounds, mode)
    # remove the temporary file
    os.remove(filename)

    # sampling and tally, tally should be defined by the mesh cell_fracs
    num_samples = 5000
    particles = []

    seed(1953)
    for i in range(num_samples):
        rands = np.array([uniform(0, 1) for x in range(6)])
        s = sampler.particle_birth(rands)
        # check w, and c for each particle
        # calculate the expected weight and cell_number
        exp_w, exp_c = _cal_exp_w_c(s, mode, cell_fracs, src_tag, bias_tag)
        assert_equal(s.w, exp_w)
        # when mode in (0, 1, 2), the set exp_c is (-1), otherwise it contains
        # several available cell number
        if mode in (0, 1, 2):
            assert(set(s.cell_list) == exp_c)
        elif mode in (3, 4, 5):
            assert(set(s.cell_list).issubset(exp_c))
        # store all the particles for the convinent of distribution check
        particles.append(s)

    # check position distribution
    # X direction follow specified distribution
    x_dis = _get_x_dis(particles, num_ve)
    x_dis_exp = _get_x_dis_exp(mode, cell_fracs, src_tag, bias_tag)
    for i in range(len(x_dis)):
        assert(abs(x_dis[i] - x_dis_exp[i]) / x_dis_exp[i] < 0.05)
    # uniform in Y and Z directions
    p_y_halfspace, p_z_halfspace = _get_p_y_z_halfspace(particles)
    assert(abs(p_y_halfspace - 0.5) / 0.5 < 0.05)
    assert(abs(p_z_halfspace - 0.5) / 0.5 < 0.05)

    # check energy distribution
    e_dis = _get_e_dis(particles, num_e_groups)
    e_dis_exp = _get_e_dis_exp(mode, cell_fracs, src_tag, bias_tag)
    for i in range(len(e_dis)):
        if e_dis_exp[i] > 0:
            assert(abs(e_dis[i] - e_dis_exp[i]) / e_dis_exp[i] < 0.05)
        else:
            assert_equal(e_dis[i], 0.0)
