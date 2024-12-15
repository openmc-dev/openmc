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

    def __init__(
        self,
        statepoint_name,
        model,
        inputs_true="inputs_true.dat",
        holes=False,
        scale_factor=10.0,
    ):

        super().__init__(statepoint_name, model, inputs_true)
        self.holes = holes  # holes in the test mesh
        self.scale_bounding_cell(scale_factor)

    def scale_bounding_cell(self, scale_factor):
        geometry = self._model.geometry
        for surface in geometry.get_all_surfaces().values():
            if surface.boundary_type != "vacuum":
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
        decimals = 10 if umesh_tally.estimator == "collision" else 8
        np.testing.assert_array_almost_equal(
            np.sort(unstructured_data), np.sort(reg_mesh_data), decimals
        )

    def get_mesh_tally_data(self, tally, structured=False):
        data = tally.get_reshaped_data(value="mean")
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
        output = glob.glob("tally*.vtk")
        output += glob.glob("tally*.e")
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
    fuel_mat.set_density("g/cc", 4.5)
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
    fuel_min_x = openmc.XPlane(-5.0, name="minimum x")
    fuel_max_x = openmc.XPlane(5.0, name="maximum x")

    fuel_min_y = openmc.YPlane(-5.0, name="minimum y")
    fuel_max_y = openmc.YPlane(5.0, name="maximum y")

    fuel_min_z = openmc.ZPlane(-5.0, name="minimum z")
    fuel_max_z = openmc.ZPlane(5.0, name="maximum z")

    fuel_cell = openmc.Cell(name="fuel")
    fuel_cell.region = (
        +fuel_min_x
        & -fuel_max_x
        & +fuel_min_y
        & -fuel_max_y
        & +fuel_min_z
        & -fuel_max_z
    )
    fuel_cell.fill = fuel_mat

    clad_min_x = openmc.XPlane(-6.0, name="minimum x")
    clad_max_x = openmc.XPlane(6.0, name="maximum x")

    clad_min_y = openmc.YPlane(-6.0, name="minimum y")
    clad_max_y = openmc.YPlane(6.0, name="maximum y")

    clad_min_z = openmc.ZPlane(-6.0, name="minimum z")
    clad_max_z = openmc.ZPlane(6.0, name="maximum z")

    clad_cell = openmc.Cell(name="clad")
    clad_cell.region = (
        -fuel_min_x
        | +fuel_max_x
        | -fuel_min_y
        | +fuel_max_y
        | -fuel_min_z
        | +fuel_max_z
    ) & (
        +clad_min_x
        & -clad_max_x
        & +clad_min_y
        & -clad_max_y
        & +clad_min_z
        & -clad_max_z
    )
    clad_cell.fill = zirc_mat

    # set bounding cell dimension to one
    # this will be updated later according to the test case parameters
    water_min_x = openmc.XPlane(x0=-1.0, name="minimum x", boundary_type="vacuum")
    water_max_x = openmc.XPlane(x0=1.0, name="maximum x", boundary_type="vacuum")

    water_min_y = openmc.YPlane(y0=-1.0, name="minimum y", boundary_type="vacuum")
    water_max_y = openmc.YPlane(y0=1.0, name="maximum y", boundary_type="vacuum")

    water_min_z = openmc.ZPlane(z0=-1.0, name="minimum z", boundary_type="vacuum")
    water_max_z = openmc.ZPlane(z0=1.0, name="maximum z", boundary_type="vacuum")

    water_cell = openmc.Cell(name="water")
    water_cell.region = (
        -clad_min_x
        | +clad_max_x
        | -clad_min_y
        | +clad_max_y
        | -clad_min_z
        | +clad_max_z
    ) & (
        +water_min_x
        & -water_max_x
        & +water_min_y
        & -water_max_y
        & +water_min_z
        & -water_max_z
    )
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
    regular_mesh_tally.scores = ["flux"]

    model.tallies = openmc.Tallies([regular_mesh_tally])

    ### Settings ###
    settings = openmc.Settings()
    settings.run_mode = "fixed source"
    settings.particles = 1000
    settings.batches = 10

    # source setup
    r = openmc.stats.Uniform(a=0.0, b=0.0)
    cos_theta = openmc.stats.Discrete(x=[1.0], p=[1.0])
    phi = openmc.stats.Discrete(x=[0.0], p=[1.0])

    space = openmc.stats.SphericalIndependent(r, cos_theta, phi)
    energy = openmc.stats.Discrete(x=[15.0e06], p=[1.0])
    source = openmc.IndependentSource(space=space, energy=energy)
    settings.source = source

    model.settings = settings

    return model


param_values = (
    ["libmesh", "moab"],  # mesh libraries
    ["collision", "tracklength"],  # estimators
    [True, False],  # geometry outside of the mesh
    [(333, 90, 77), None],
)  # location of holes in the mesh
test_cases = []
for i, (lib, estimator, ext_geom, holes) in enumerate(product(*param_values)):
    test_cases.append(
        {
            "library": lib,
            "estimator": estimator,
            "external_geom": ext_geom,
            "holes": holes,
            "inputs_true": "inputs_true{}.dat".format(i),
        }
    )


@pytest.mark.parametrize("test_opts", test_cases)
def test_unstructured_mesh_tets(model, test_opts):
    # skip the test if the library is not enabled
    if test_opts["library"] == "moab" and not openmc.lib._dagmc_enabled():
        pytest.skip("DAGMC (and MOAB) mesh not enbaled in this build.")

    if test_opts["library"] == "libmesh" and not openmc.lib._libmesh_enabled():
        pytest.skip("LibMesh is not enabled in this build.")

    # skip the tracklength test for libmesh
    if test_opts["library"] == "libmesh" and test_opts["estimator"] == "tracklength":
        pytest.skip("Tracklength tallies are not supported using libmesh.")

    if test_opts["holes"]:
        mesh_filename = "test_mesh_tets_w_holes.e"
    else:
        mesh_filename = "test_mesh_tets.e"

    # add reference mesh tally
    regular_mesh_tally = model.tallies[0]
    regular_mesh_tally.estimator = test_opts["estimator"]

    # add analagous unstructured mesh tally
    uscd_mesh = openmc.UnstructuredMesh(mesh_filename, test_opts["library"])
    if test_opts["library"] == "moab":
        uscd_mesh.options = "MAX_DEPTH=15;PLANE_SET=2"
    uscd_filter = openmc.MeshFilter(mesh=uscd_mesh)

    # create tallies
    uscd_tally = openmc.Tally(name="unstructured mesh tally")
    uscd_tally.filters = [uscd_filter]
    uscd_tally.scores = ["flux"]
    uscd_tally.estimator = test_opts["estimator"]
    model.tallies.append(uscd_tally)

    # modify model geometry according to test opts
    if test_opts["external_geom"]:
        scale_factor = 15.0
    else:
        scale_factor = 10.0

    harness = UnstructuredMeshTest(
        "statepoint.10.h5",
        model,
        test_opts["inputs_true"],
        test_opts["holes"],
        scale_factor,
    )
    harness.main()


@pytest.mark.skipif(
    not openmc.lib._libmesh_enabled(), reason="LibMesh is not enabled in this build."
)
def test_unstructured_mesh_hexes(model):
    regular_mesh_tally = model.tallies[0]
    regular_mesh_tally.estimator = "collision"

    # add analagous unstructured mesh tally
    uscd_mesh = openmc.UnstructuredMesh("test_mesh_hexes.e", "libmesh")
    uscd_filter = openmc.MeshFilter(mesh=uscd_mesh)

    # create tallies
    uscd_tally = openmc.Tally(name="unstructured mesh tally")
    uscd_tally.filters = [uscd_filter]
    uscd_tally.scores = ["flux"]
    uscd_tally.estimator = "collision"
    model.tallies.append(uscd_tally)

    harness = UnstructuredMeshTest("statepoint.10.h5", model)
    harness.ELEM_PER_VOXEL = 1

    harness.main()
