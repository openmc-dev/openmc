import glob
from pathlib import Path

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
                if isinstance(mesh, XDGMesh):
                    xdg_mesh = mesh
                    break

            assert xdg_mesh is not None
            assert Path(xdg_mesh.filename).name == Path(self.mesh_filename).name

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
                        xdg_mesh_data = self.get_mesh_tally_data(tally, structured=True)

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
            Path(f).unlink(missing_ok=True)


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
    fuel_boundary = openmc.model.RectangularParallelepiped(*[-5.0, 5.0]*3)
    clad_boundary = openmc.model.RectangularParallelepiped(*[-6.0, 6.0]*3)

    fuel_cell = openmc.Cell(name="fuel")
    fuel_cell.region = -fuel_boundary
    fuel_cell.fill = fuel_mat

    clad_cell = openmc.Cell(name="clad")
    clad_cell.region = +fuel_boundary & -clad_boundary
    clad_cell.fill = zirc_mat

    # set bounding cell dimension to one
    # this will be updated later according to the test case parameters
    boundary = openmc.model.RectangularParallelepiped(*[-1.0, 1.0]*3,
                                                      boundary_type='vacuum')

    water_cell = openmc.Cell(name="water")
    water_cell.region = +clad_boundary & -boundary
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
    model.settings.run_mode = 'fixed source'
    model.settings.particles = 1000
    model.settings.batches = 10

    # source setup
    r = openmc.stats.Uniform(a=0.0, b=9.0)
    cos_theta = openmc.stats.delta_function(1.0)
    phi = openmc.stats.delta_function(0.0)

    space = openmc.stats.SphericalIndependent(r, cos_theta, phi)
    energy = openmc.stats.delta_function(15e6)
    source = openmc.IndependentSource(space=space, energy=energy)
    model.settings.source = source

    return model


MESH_CASES = (
    {
        "mesh_filename": "test_mesh_tets.e",
        "mesh_kind": "tet",
        "holes": None,
        "elem_per_voxel": 12,
        "external_geom": False,
        "libraries": ("moab", "libmesh"),
    },
    {
        "mesh_filename": "test_mesh_tets.exo",
        "mesh_kind": "tet",
        "holes": None,
        "elem_per_voxel": 12,
        "external_geom": False,
        "libraries": ("libmesh",),
    },
    {
        "mesh_filename": "test_mesh_tets_w_holes.e",
        "mesh_kind": "tet",
        "holes": (333, 90, 77),
        "elem_per_voxel": 12,
        "external_geom": False,
        "libraries": ("moab", "libmesh"),
    },
    {
        "mesh_filename": "test_mesh_tets_w_holes.exo",
        "mesh_kind": "tet",
        "holes": (333, 90, 77),
        "elem_per_voxel": 12,
        "external_geom": False,
        "libraries": ("libmesh",),
    },
    {
        "mesh_filename": "test_mesh_tets.e",
        "mesh_kind": "tet",
        "holes": None,
        "elem_per_voxel": 12,
        "external_geom": True,
        "libraries": ("moab", "libmesh"),
    },
    {
        "mesh_filename": "test_mesh_tets.exo",
        "mesh_kind": "tet",
        "holes": None,
        "elem_per_voxel": 12,
        "external_geom": True,
        "libraries": ("libmesh",),
    },
    {
        "mesh_filename": "test_mesh_tets_w_holes.e",
        "mesh_kind": "tet",
        "holes": (333, 90, 77),
        "elem_per_voxel": 12,
        "external_geom": True,
        "libraries": ("moab", "libmesh"),
    },
    {
        "mesh_filename": "test_mesh_tets_w_holes.exo",
        "mesh_kind": "tet",
        "holes": (333, 90, 77),
        "elem_per_voxel": 12,
        "external_geom": True,
        "libraries": ("libmesh",),
    }
)

test_cases = []
test_case_ids = []
for i, mesh_case in enumerate(MESH_CASES):
    for library in mesh_case["libraries"]:
        test_cases.append({
            "inputs_true": f"inputs_true{i}_{library}.dat",
            "library": library,
            **{k: v for k, v in mesh_case.items() if k != "libraries"},
        })
        mesh_stem = Path(mesh_case["mesh_filename"]).stem
        geom_mode = "external-geom" if mesh_case["external_geom"] else "bounded-geom"
        holes = "holes" if mesh_case["holes"] else "solid"
        test_case_ids.append(f"{library}-{mesh_stem}-{holes}-{geom_mode}")


@pytest.mark.parametrize("test_opts", test_cases, ids=test_case_ids)
def test_xdg_mesh_tallies(model, test_opts):
    if not openmc.lib._xdg_enabled():
        pytest.skip("XDG is not enabled in this build.")

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
        scale_factor=15.0 if test_opts["external_geom"] else 10.0,
    )
    harness.ELEM_PER_VOXEL = test_opts["elem_per_voxel"]
    harness.main()
