import filecmp
import glob
from itertools import product
import os
import warnings

import openmc
import openmc.lib
import numpy as np

import pytest
from tests.testing_harness import PyAPITestHarness


class UnstructuredMeshTest(PyAPITestHarness):

    ELEM_PER_VOXEL = 12

    def __init__(self,
                 statepoint_name,
                 model,
                 inputs_true='inputs_true.dat',
                 holes=False,
                 scale_factor=10.0):

        super().__init__(statepoint_name, model, inputs_true)
        self.holes = holes # holes in the test mesh
        self.scale_bounding_cell(scale_factor)

    def scale_bounding_cell(self, scale_factor):
        geometry = self._model.geometry
        for surface in geometry.get_all_surfaces().values():
            if surface.boundary_type != 'vacuum':
                continue
            for coeff in surface._coefficients:
                surface._coefficients[coeff] *= scale_factor

    def _compare_results(self):
        with openmc.StatePoint(self._sp_name) as sp:
            # check some properties of the unstructured mesh
            umesh = None
            for m in sp.meshes.values():
                if isinstance(m, openmc.UnstructuredMesh):
                    umesh = m
            assert umesh is not None

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
                        reg_mesh_data = self.get_mesh_tally_data(tally)
                        if self.holes:
                            reg_mesh_data = np.delete(reg_mesh_data, self.holes)
                    else:
                        umesh_tally = tally
                        unstructured_data = self.get_mesh_tally_data(tally, True)

        # we expect these results to be the same to within at least ten
        # decimal places
        decimals = 10 if umesh_tally.estimator == 'collision' else 6
        np.testing.assert_array_almost_equal(np.sort(unstructured_data),
                                            np.sort(reg_mesh_data),
                                            decimals)

    def get_mesh_tally_data(self, tally, structured=False):
        data = tally.get_reshaped_data(value='mean')
        if structured:
           data = data.reshape((-1, self.ELEM_PER_VOXEL))
        else:
            data.shape = (data.size, 1)
        return np.sum(data, axis=1)

    def update_results(self):
        """Update results_true.dat and inputs_true.dat"""
        try:
            self._build_inputs()
            inputs = self._get_inputs()
            self._write_inputs(inputs)
            self._overwrite_inputs()
            self._run_openmc()
            self._test_output_created()
        finally:
            self._cleanup()

    def _cleanup(self):
        super()._cleanup()
        output = glob.glob('tally*.vtk')
        output += glob.glob('tally*.e')
        for f in output:
            if os.path.exists(f):
                os.remove(f)


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
    fuel_cell = openmc.Cell(name="fuel", region=-fuel_box)
    fuel_cell.fill = fuel_mat

    clad_box = openmc.model.RectangularParallelepiped(-6.0, 6.0, -6.0, 6.0, -6.0, 6.0)
    clad_cell = openmc.Cell(name="clad", region=-clad_box & +fuel_box)
    clad_cell.fill = zirc_mat

    # set bounding cell dimension to one
    # this will be updated later according to the test case parameters
    water_box = openmc.model.RectangularParallelepiped(-1.0, 1.0, -1.0, 1.0, -1.0, 1.0, boundary_type='vacuum')
    water_cell = openmc.Cell(name="water", region=-water_box & +clad_box)
    water_cell.fill = water_mat

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
    settings = openmc.Settings()
    settings.run_mode = 'fixed source'
    settings.particles = 1000
    settings.batches = 10

    # source setup
    space = openmc.stats.spherical_uniform(r_outer = 9.0)
    energy = openmc.stats.Discrete(x=[15.e+06], p=[1.0])
    source = openmc.IndependentSource(space=space, energy=energy)
    settings.source = source

    model.settings = settings

    return model


param_values = (['libmesh', 'moab'], # mesh libraries
                ['native', 'xdg'], # mesh interfaces
                ['collision', 'tracklength'], # estimators
                [True, False], # geometry outside of the mesh
                [(333, 90, 77), None]) # location of holes in the mesh
test_cases = []
for i, (lib, interface, estimator, ext_geom, holes) in enumerate(product(*param_values)):
    if lib == 'libmesh' and estimator == 'tracklength':
        continue
    test_cases.append({'library' : lib,
                       'interface': interface,
                       'estimator' : estimator,
                       'external_geom' : ext_geom,
                       'holes' : holes,
                       'inputs_true' : f'inputs_tets_true{i}.dat'})

def param_ids(test_case):
    return f"{test_case['library']}_{test_case['interface']}_{test_case['estimator']}_holes_{test_case['holes']}_external_geom_{test_case['external_geom']}"

@pytest.mark.parametrize("test_opts", test_cases, ids=param_ids)
def test_unstructured_mesh_tets(model, test_opts):
    # skip the test if appropriate libraries or interfaces are not enabled
    if test_opts['interface'] == 'xdg' and not openmc.lib.feature_enabled('xdg'):
        pytest.skip("XDG interface is not enabled in this build.")
    elif test_opts['interface'] == 'native':
        if test_opts['library'] == 'moab' and not openmc.lib.feature_enabled('dagmc'):
            pytest.skip("DAGMC (and MOAB) mesh not enabled in this build.")

        if test_opts['library'] == 'libmesh' and not openmc.lib.feature_enabled('libmesh'):
            pytest.skip("LibMesh is not enabled in this build.")

    # skip the tracklength test for libmesh
    if test_opts['library'] == 'libmesh' and \
       test_opts['estimator'] == 'tracklength' and \
       test_opts['interface'] != 'xdg':
       pytest.skip("Tracklength tallies are not supported using libmesh.")

    if test_opts['holes']:
        mesh_filename = "test_mesh_tets_w_holes.e"
    else:
        mesh_filename = "test_mesh_tets.e"

    interface = test_opts['interface']

    # add reference mesh tally
    regular_mesh_tally = model.tallies[0]
    regular_mesh_tally.estimator = test_opts['estimator']

    # add analagous unstructured mesh tally
    uscd_mesh = openmc.UnstructuredMesh(mesh_filename, test_opts['library'])
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

    harness = UnstructuredMeshTest('statepoint.10.h5',
                                   model,
                                   test_opts['inputs_true'],
                                   test_opts['holes'],
                                   scale_factor)
    harness.main()


param_values = (['libmesh', 'moab'], # mesh libraries
                ['native', 'xdg'], # mesh interfaces
                ['collision', 'tracklength']) # estimators
test_cases = []
for i, (lib, interface, estimator) in enumerate(product(*param_values)):
    if lib == 'moab' and interface != 'xdg':
        continue
    if lib == 'libmesh' and estimator == 'tracklength':
        continue
    test_cases.append((lib, interface, estimator, f'inputs_hexes_true{i}.dat'))

@pytest.mark.parametrize("test_opts", test_cases, ids=lambda x: f"{x[0]}_{x[1]}_{x[2]}")
def test_unstructured_mesh_hexes(model, test_opts):

    library, interface, estimator, inputs_true = test_opts

    if library == 'libmesh' and not openmc.lib.feature_enabled('libmesh'):
        pytest.skip("LibMesh is not enabled in this build.")
    if library == 'moab' and not openmc.lib.feature_enabled('dagmc'):
        pytest.skip("DAGMC (and MOAB) mesh not enabled in this build.")
    if interface == 'xdg' and not openmc.lib.feature_enabled('xdg'):
        pytest.skip("XDG interface is not enabled in this build.")

    regular_mesh_tally = model.tallies[0]
    regular_mesh_tally.estimator = estimator

    # add analagous unstructured mesh tally
    filename = "test_mesh_hexes.e" if library == 'libmesh' else "test_mesh_hexes.exo"
    uscd_mesh = openmc.UnstructuredMesh(filename, library)
    uscd_mesh.interface = interface
    uscd_filter = openmc.MeshFilter(mesh=uscd_mesh)

    # create tallies
    uscd_tally = openmc.Tally(name="unstructured mesh tally")
    uscd_tally.filters = [uscd_filter]
    uscd_tally.scores = ['flux']
    uscd_tally.estimator = estimator
    model.tallies.append(uscd_tally)

    harness = UnstructuredMeshTest('statepoint.10.h5',
                                   model,
                                   inputs_true)
    harness.ELEM_PER_VOXEL = 1

    harness.main()
