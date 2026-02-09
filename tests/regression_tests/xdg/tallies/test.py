import glob
import os

import numpy as np
import pytest

import openmc
import openmc.lib
from openmc.xdg import XDGMesh

from tests.testing_harness import PyAPITestHarness


class XDGMeshTallyTest(PyAPITestHarness):

    ELEM_PER_VOXEL = 12

    def __init__(
        self,
        statepoint_name,
        model,
        inputs_true,
        mesh_filename,
        mesh_kind,
        holes=None,
        scale_factor=10.0,
    ):
        super().__init__(statepoint_name, model, inputs_true)
        print(f"Running test with mesh file: {mesh_filename} and inputs file: {inputs_true}")
        self.mesh_filename = mesh_filename
        self.mesh_kind = mesh_kind
        self.holes = holes
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
            xdg_mesh = None
            for mesh in sp.meshes.values():
                if isinstance(mesh, openmc.UnstructuredMesh) and mesh.library == "xdg":
                    xdg_mesh = mesh
                    break

            assert xdg_mesh is not None
            assert os.path.basename(xdg_mesh.filename) == os.path.basename(
                self.mesh_filename
            )

            if self.mesh_kind == "tet":
                exp_vertex = (-10.0, -10.0, -10.0)
                exp_centroid = (-8.75, -9.75, -9.25)
            else:
                exp_vertex = (-10.0, -10.0, 10.0)
                exp_centroid = (-9.0, -9.0, 9.0)

            np.testing.assert_array_equal(xdg_mesh.vertices[0], exp_vertex)
            np.testing.assert_array_equal(xdg_mesh.centroid(0), exp_centroid)

            reg_mesh_data = None
            xdg_mesh_data = None
            xdg_tally = None
            for tally in sp.tallies.values():
                if tally.contains_filter(openmc.MeshFilter):
                    flt = tally.find_filter(openmc.MeshFilter)
                    if isinstance(flt.mesh, openmc.RegularMesh):
                        reg_mesh_data = self.get_mesh_tally_data(tally)
                    else:
                        xdg_tally = tally
                        xdg_mesh_data = self.get_mesh_tally_data(tally, True)

        assert reg_mesh_data is not None
        assert xdg_mesh_data is not None

        if self.holes:
            reg_mesh_data = np.delete(reg_mesh_data, self.holes)

        decimals = 10 if xdg_tally.estimator == 'collision' else 8
        np.testing.assert_array_almost_equal(
            np.sort(xdg_mesh_data),
            np.sort(reg_mesh_data),
            decimals
        )

    def get_mesh_tally_data(self, tally, structured=False):
        data = tally.get_reshaped_data(value='mean')
        if structured:
            data = data.reshape((-1, self.ELEM_PER_VOXEL))
        else:
            data.shape = (data.size, 1)
        return np.sum(data, axis=1)

    def update_results(self):
        """Update inputs_true.dat without storing results_true.dat."""
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
    fuel_min_x = openmc.XPlane(-5.0, name="minimum x")
    fuel_max_x = openmc.XPlane(5.0, name="maximum x")

    fuel_min_y = openmc.YPlane(-5.0, name="minimum y")
    fuel_max_y = openmc.YPlane(5.0, name="maximum y")

    fuel_min_z = openmc.ZPlane(-5.0, name="minimum z")
    fuel_max_z = openmc.ZPlane(5.0, name="maximum z")

    fuel_cell = openmc.Cell(name="fuel")
    fuel_cell.region = +fuel_min_x & -fuel_max_x & \
                       +fuel_min_y & -fuel_max_y & \
                       +fuel_min_z & -fuel_max_z
    fuel_cell.fill = fuel_mat

    clad_min_x = openmc.XPlane(-6.0, name="minimum x")
    clad_max_x = openmc.XPlane(6.0, name="maximum x")

    clad_min_y = openmc.YPlane(-6.0, name="minimum y")
    clad_max_y = openmc.YPlane(6.0, name="maximum y")

    clad_min_z = openmc.ZPlane(-6.0, name="minimum z")
    clad_max_z = openmc.ZPlane(6.0, name="maximum z")

    clad_cell = openmc.Cell(name="clad")
    clad_cell.region = (-fuel_min_x | +fuel_max_x |
                        -fuel_min_y | +fuel_max_y |
                        -fuel_min_z | +fuel_max_z) & \
                        (+clad_min_x & -clad_max_x &
                         +clad_min_y & -clad_max_y &
                         +clad_min_z & -clad_max_z)
    clad_cell.fill = zirc_mat

    # set bounding cell dimension to one
    # this will be updated later according to the test case parameters
    water_min_x = openmc.XPlane(x0=-1.0,
                                name="minimum x",
                                boundary_type='vacuum')
    water_max_x = openmc.XPlane(x0=1.0,
                                name="maximum x",
                                boundary_type='vacuum')

    water_min_y = openmc.YPlane(y0=-1.0,
                                name="minimum y",
                                boundary_type='vacuum')
    water_max_y = openmc.YPlane(y0=1.0,
                                name="maximum y",
                                boundary_type='vacuum')

    water_min_z = openmc.ZPlane(z0=-1.0,
                                name="minimum z",
                                boundary_type='vacuum')
    water_max_z = openmc.ZPlane(z0=1.0,
                                name="maximum z",
                                boundary_type='vacuum')

    water_cell = openmc.Cell(name="water")
    water_cell.region = (-clad_min_x | +clad_max_x |
                         -clad_min_y | +clad_max_y |
                         -clad_min_z | +clad_max_z) & \
                         (+water_min_x & -water_max_x &
                          +water_min_y & -water_max_y &
                          +water_min_z & -water_max_z)
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
    r = openmc.stats.Uniform(a=0.0, b=0.0)
    cos_theta = openmc.stats.Discrete(x=[1.0], p=[1.0])
    phi = openmc.stats.Discrete(x=[0.0], p=[1.0])

    space = openmc.stats.SphericalIndependent(r, cos_theta, phi)
    energy = openmc.stats.Discrete(x=[15.e+06], p=[1.0])
    source = openmc.IndependentSource(space=space, energy=energy)
    settings.source = source

    model.settings = settings

    return model


MESH_CASES = (
    {
        "mesh_filename": "test_mesh_tets.e",
        "mesh_kind": "tet",
        "holes": None,
        "elem_per_voxel": 12,
        "libraries": ("moab", "libmesh"),
    },
    {
        "mesh_filename": "test_mesh_tets.exo",
        "mesh_kind": "tet",
        "holes": None,
        "elem_per_voxel": 12,
        "libraries": ("libmesh",),
    },
    {
        "mesh_filename": "test_mesh_tets_w_holes.e",
        "mesh_kind": "tet",
        "holes": (333, 90, 77),
        "elem_per_voxel": 12,
        "libraries": ("moab", "libmesh"),
    },
    {
        "mesh_filename": "test_mesh_tets_w_holes.exo",
        "mesh_kind": "tet",
        "holes": (333, 90, 77),
        "elem_per_voxel": 12,
        "libraries": ("libmesh",),
    },
    # {
    #     "mesh_filename": "test_mesh_hexes.e",
    #     "mesh_kind": "hex",
    #     "holes": None,
    #     "elem_per_voxel": 1,
    #     "libraries": ("libmesh",),
    # },
    # {
    #     "mesh_filename": "test_mesh_hexes.exo",
    #     "mesh_kind": "hex",
    #     "holes": None,
    #     "elem_per_voxel": 1,
    #     "libraries": ("libmesh",),
    # },
)

test_cases = []
for i, mesh_case in enumerate(MESH_CASES):
    for library in mesh_case["libraries"]:
        test_cases.append({
            "inputs_true": f"inputs_true{i}_{library}.dat",
            "library": library,
            **{k: v for k, v in mesh_case.items() if k != "libraries"},
        })


@pytest.mark.parametrize("test_opts", test_cases)
def test_xdg_mesh_tallies(model, test_opts):
    if not openmc.lib._xdg_enabled():
        pytest.skip("XDG is not enabled in this build.")

    if test_opts["library"] == "moab" and not openmc.lib._dagmc_enabled():
        pytest.skip("DAGMC (and MOAB) mesh not enabled in this build.")

    if test_opts["library"] == "libmesh" and not openmc.lib._libmesh_enabled():
        pytest.skip("LibMesh is not enabled in this build.")

    # reference mesh tally
    regular_mesh_tally = model.tallies[0]
    regular_mesh_tally.estimator = 'collision'

    # add analogous XDG mesh tally
    xdg_mesh = XDGMesh(test_opts["mesh_filename"], test_opts["library"])
    xdg_filter = openmc.MeshFilter(mesh=xdg_mesh)

    xdg_tally = openmc.Tally(name="xdg mesh tally")
    xdg_tally.filters = [xdg_filter]
    xdg_tally.scores = ['flux']
    xdg_tally.estimator = 'collision'
    model.tallies.append(xdg_tally)

    harness = XDGMeshTallyTest(
        'statepoint.10.h5',
        model,
        test_opts["inputs_true"],
        test_opts["mesh_filename"],
        test_opts["mesh_kind"],
        test_opts["holes"],
        scale_factor=10.0,
    )
    harness.ELEM_PER_VOXEL = test_opts["elem_per_voxel"]
    harness.main()
