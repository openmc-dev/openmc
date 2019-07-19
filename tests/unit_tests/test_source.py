import openmc
import openmc.stats
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

def test_source():
    space = openmc.stats.Point()
    energy = openmc.stats.Discrete([1.0e6], [1.0])
    angle = openmc.stats.Isotropic()

    src = openmc.Source(space=space, angle=angle, energy=energy)
    assert src.space == space
    assert src.angle == angle
    assert src.energy == energy

    elem = src.to_xml_element()
    assert 'strength' in elem.attrib
    assert elem.find('space') is not None
    assert elem.find('angle') is not None
    assert elem.find('energy') is not None

    src = openmc.Source.from_xml_element(elem)
    assert isinstance(src.angle, openmc.stats.Isotropic)
    assert src.space.xyz == [0.0, 0.0, 0.0]
    assert src.energy.x == [1.0e6]
    assert src.energy.p == [1.0]
    assert src.strength == 1.0


def test_source_file():
    filename = 'source.h5'
    src = openmc.Source(filename=filename)
    assert src.file == filename

    elem = src.to_xml_element()
    assert 'strength' in elem.attrib
    assert 'file' in elem.attrib

# The following tests are designed for source_sampling,
# copied from pyne/tests/test_source_sampling.py
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
