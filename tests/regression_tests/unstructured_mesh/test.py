from itertools import product
from pathlib import Path

import openmc
import openmc.lib
import numpy as np

import pytest
from tests.regression_tests import config


def run_and_check(model, tmp_path, mesh_filename, holes=None,
                  scale_factor=10.0, elements_per_voxel=12):
    """Compare unstructured tallies with a regular mesh in the same run."""
    for surface in model.geometry.get_all_surfaces().values():
        if surface.boundary_type == 'vacuum':
            for coeff in surface._coefficients:
                surface._coefficients[coeff] *= scale_factor

    if config['build_inputs']:
        model.export_to_model_xml(tmp_path / 'model.xml')
        return

    mpi_args = ([config['mpiexec'], '-n', config['mpi_np']]
                if config['mpi'] else None)
    statepoint = model.run(cwd=tmp_path, openmc_exec=config['exe'],
                           mpi_args=mpi_args, event_based=config['event'])
    with openmc.StatePoint(statepoint) as sp:
        # check some properties of the unstructured mesh
        umesh = None
        for m in sp.meshes.values():
            if isinstance(m, openmc.UnstructuredMesh):
                umesh = m
        assert umesh is not None
        assert Path(umesh.filename).name == mesh_filename

        # check that the first element centroid is correct
        # this will depend on whether the tet mesh or hex mesh
        # file is being used in this test
        if umesh.element_types[0] == umesh._LINEAR_TET:
            exp_vertex = (-10.0, -10.0, -10.0)
            exp_centroid = (-8.75, -9.75, -9.25)
        else:
            exp_vertex = (-10.0, -10.0, 10.0)
            exp_centroid = (-9.0, -9.0, 9.0)

        np.testing.assert_array_equal(umesh.vertices[0], exp_vertex)
        np.testing.assert_array_equal(umesh.centroid(0), exp_centroid)

        # loop over the tallies and get data
        for tally in sp.tallies.values():
            # find the regular and unstructured meshes
            if tally.contains_filter(openmc.MeshFilter):
                flt = tally.find_filter(openmc.MeshFilter)

                if isinstance(flt.mesh, openmc.RegularMesh):
                    reg_mesh_data = get_mesh_tally_data(tally)
                    if holes:
                        reg_mesh_data = np.delete(reg_mesh_data, holes)
                else:
                    umesh_tally = tally
                    unstructured_data = get_mesh_tally_data(
                        tally, elements_per_voxel)

    # Collision tallies agree to ten decimal places; tracklength to six.
    decimals = 10 if umesh_tally.estimator == 'collision' else 6
    np.testing.assert_array_almost_equal(np.sort(unstructured_data),
                                        np.sort(reg_mesh_data),
                                        decimals)


def get_mesh_tally_data(tally, elements_per_voxel=1):
    data = tally.get_reshaped_data(value='mean')
    return data.reshape(-1, elements_per_voxel).sum(axis=1)


@pytest.fixture
def model():
    openmc.reset_auto_ids()

    model = openmc.Model()

    ### Materials ###
    materials = openmc.Materials()

    fuel_mat = openmc.Material(name="fuel")
    fuel_mat.add_nuclide("U235", 1.0)
    fuel_mat.set_density('g/cc', 4.5)
    materials.append(fuel_mat)

    zirc_mat = openmc.Material(name="zircaloy")
    zirc_mat.add_element("Zr", 1.0)
    zirc_mat.set_density("g/cc", 5.77)
    materials.append(zirc_mat)

    water_mat = openmc.Material(name="water")
    water_mat.add_nuclide("H1", 2.0)
    water_mat.add_nuclide("O16", 1.0)
    water_mat.set_density("atom/b-cm", 0.07416)
    materials.append(water_mat)

    model.materials = materials

    ### Geometry ###
    fuel_box = openmc.model.RectangularParallelepiped(-5.0, 5.0, -5.0, 5.0, -5.0, 5.0)
    fuel_cell = openmc.Cell(name="fuel", region=-fuel_box, fill=fuel_mat)

    clad_box = openmc.model.RectangularParallelepiped(-6.0, 6.0, -6.0, 6.0, -6.0, 6.0)
    clad_cell = openmc.Cell(name="clad", region=-clad_box & +fuel_box, fill=zirc_mat)

    # set bounding cell dimension to one
    # this will be updated later according to the test case parameters
    water_box = openmc.model.RectangularParallelepiped(-1.0, 1.0, -1.0, 1.0, -1.0, 1.0, boundary_type='vacuum')
    water_cell = openmc.Cell(name="water", region=-water_box & +clad_box, fill=water_mat)

    # create a containing universe
    model.geometry = openmc.Geometry([fuel_cell, clad_cell, water_cell])

    ### Reference Tally ###

    # create meshes and mesh filters
    regular_mesh = openmc.RegularMesh()
    regular_mesh.dimension = (10, 10, 10)
    regular_mesh.lower_left = (-10.0, -10.0, -10.0)
    regular_mesh.upper_right = (10.0, 10.0, 10.0)

    regular_mesh_filter = openmc.MeshFilter(mesh=regular_mesh)
    regular_mesh_tally = openmc.Tally(name="regular mesh tally")
    regular_mesh_tally.filters = [regular_mesh_filter]
    regular_mesh_tally.scores = ['flux']

    model.tallies = openmc.Tallies([regular_mesh_tally])

    ### Settings ###
    model.settings.run_mode = 'fixed source'
    model.settings.particles = 1000
    model.settings.batches = 10

    # source setup
    space = openmc.stats.spherical_uniform(r_outer=9.0)
    energy = openmc.stats.delta_function(15.e6)
    source = openmc.IndependentSource(space=space, energy=energy)
    model.settings.source = source

    return model


param_values = (['libmesh', 'moab'], # mesh libraries
                ['native', 'xdg'], # mesh interfaces
                ['collision', 'tracklength'], # estimators
                [True, False], # geometry outside of the mesh
                [(333, 90, 77), None]) # location of holes in the mesh
test_cases = []
for lib, interface, estimator, ext_geom, holes in product(*param_values):
    if lib == 'libmesh' and interface == 'native' and estimator == 'tracklength':
        continue
    test_cases.append({'library' : lib,
                       'interface': interface,
                       'estimator' : estimator,
                       'external_geom' : ext_geom,
                       'holes' : holes})

# Retain the XDG collision tests with source sites along the positive z-axis.
# The .exo cases also exercise libMesh's alternate filename extension.
for external_geom, holes, extension in product(
        (False, True), (None, (333, 90, 77)), ('e', 'exo')):
    stem = 'test_mesh_tets_w_holes' if holes else 'test_mesh_tets'
    libraries = ('moab', 'libmesh') if extension == 'e' else ('libmesh',)
    for library in libraries:
        test_cases.append({'library': library,
                           'interface': 'xdg',
                           'estimator': 'collision',
                           'external_geom': external_geom,
                           'holes': holes,
                           'mesh_filename': f'{stem}.{extension}',
                           'source_kind': 'axis'})

def param_ids(test_case):
    case_id = (
        f"{test_case['library']}_{test_case['interface']}_{test_case['estimator']}"
        f"_holes_{test_case['holes']}_external_geom_{test_case['external_geom']}"
    )
    if test_case.get('source_kind') == 'axis':
        case_id += f"_source_axis_file_{Path(test_case['mesh_filename']).suffix[1:]}"
    return case_id

@pytest.mark.parametrize("test_opts", test_cases, ids=param_ids)
def test_unstructured_mesh_tets(model, test_opts, tmp_path):
    # skip the test if appropriate libraries or interfaces are not enabled
    if test_opts['interface'] == 'xdg' and not openmc.lib.feature_enabled('xdg'):
        pytest.skip("XDG interface is not enabled in this build.")
    elif test_opts['interface'] == 'native':
        if test_opts['library'] == 'moab' and not openmc.lib.feature_enabled('dagmc'):
            pytest.skip("DAGMC (and MOAB) mesh not enabled in this build.")

        if test_opts['library'] == 'libmesh' and not openmc.lib.feature_enabled('libmesh'):
            pytest.skip("LibMesh is not enabled in this build.")

    mesh_filename = test_opts.get('mesh_filename')
    if mesh_filename is None:
        mesh_filename = ("test_mesh_tets_w_holes.e" if test_opts['holes']
                         else "test_mesh_tets.e")

    if test_opts.get('source_kind') == 'axis':
        r = openmc.stats.Uniform(a=0.0, b=9.0)
        cos_theta = openmc.stats.delta_function(1.0)
        phi = openmc.stats.delta_function(0.0)
        space = openmc.stats.SphericalIndependent(r, cos_theta, phi)
        energy = openmc.stats.delta_function(15e6)
        model.settings.source = openmc.IndependentSource(space=space, energy=energy)

    interface = test_opts['interface']

    # add reference mesh tally
    regular_mesh_tally = model.tallies[0]
    regular_mesh_tally.estimator = test_opts['estimator']

    # add analagous unstructured mesh tally
    uscd_mesh = openmc.UnstructuredMesh(
        Path(__file__).with_name(mesh_filename), test_opts['library'])
    if test_opts['library'] == 'moab':
        uscd_mesh.options = 'MAX_DEPTH=15;PLANE_SET=2'
    uscd_filter = openmc.MeshFilter(mesh=uscd_mesh)

    uscd_mesh.interface = interface

    # create tallies
    uscd_tally = openmc.Tally(name="unstructured mesh tally")
    uscd_tally.filters = [uscd_filter]
    uscd_tally.scores = ['flux']
    uscd_tally.estimator = test_opts['estimator']
    model.tallies.append(uscd_tally)

    # modify model geometry according to test opts
    if test_opts['external_geom']:
        scale_factor = 15.0
    else:
        scale_factor = 10.0

    run_and_check(model, tmp_path, mesh_filename, test_opts['holes'], scale_factor)


param_values = (['libmesh', 'moab'], # mesh libraries
                ['native', 'xdg'], # mesh interfaces
                ['collision', 'tracklength']) # estimators
test_cases = []
for lib, interface, estimator in product(*param_values):
    if lib == 'moab' and interface != 'xdg':
        continue
    if lib == 'libmesh' and interface == 'native' and estimator == 'tracklength':
        continue
    test_cases.append((lib, interface, estimator))

@pytest.mark.parametrize("test_opts", test_cases, ids=lambda x: f"{x[0]}_{x[1]}_{x[2]}")
def test_unstructured_mesh_hexes(model, test_opts, tmp_path):

    library, interface, estimator = test_opts

    if library == 'libmesh' and interface == 'native' and not openmc.lib.feature_enabled('libmesh'):
        pytest.skip("LibMesh is not enabled in this build.")
    if library == 'moab' and interface == 'native' and not openmc.lib.feature_enabled('dagmc'):
        pytest.skip("DAGMC (and MOAB) mesh not enabled in this build.")
    if interface == 'xdg' and not openmc.lib.feature_enabled('xdg'):
        pytest.skip("XDG interface is not enabled in this build.")

    regular_mesh_tally = model.tallies[0]
    regular_mesh_tally.estimator = estimator

    # add analagous unstructured mesh tally
    filename = "test_mesh_hexes.e" if library == 'libmesh' else "test_mesh_hexes.exo"
    uscd_mesh = openmc.UnstructuredMesh(Path(__file__).with_name(filename), library)
    uscd_mesh.interface = interface
    uscd_filter = openmc.MeshFilter(mesh=uscd_mesh)

    # create tallies
    uscd_tally = openmc.Tally(name="unstructured mesh tally")
    uscd_tally.filters = [uscd_filter]
    uscd_tally.scores = ['flux']
    uscd_tally.estimator = estimator
    model.tallies.append(uscd_tally)

    run_and_check(model, tmp_path, filename, elements_per_voxel=1)
