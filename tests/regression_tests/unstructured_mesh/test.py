from collections import defaultdict
import glob
from itertools import product
import os

import openmc
import openmc.lib
import numpy as np

import pytest
from tests.testing_harness import PyAPITestHarness

pytestmark = pytest.mark.skipif(
    not openmc.lib._dagmc_enabled(),
    reason="Mesh library is not available.")

class UnstructuredMeshTest(PyAPITestHarness):

    def __init__(self, statepoint_name, **kwargs):
        super().__init__(statepoint_name)

        # defaults
        self.estimator = "collision" # tally estimator type
        self.external_geom = False # geometry size matches mesh
        self.mesh_has_holes = False # holes in the mesh
        self.holes = ()
        self.mesh_filename = "test_mesh_tets.h5m" # mesh file to use

        # set parameters for the test
        self.estimator = kwargs.get('estimator', self.estimator)
        self.external_geom = kwargs.get('external_geom', self.external_geom)
        self.holes = kwargs.get('holes', self.holes)
        if self.holes:
            self.mesh_filename = "test_mesh_tets_w_holes.h5m"

    def _build_inputs(self):
        ### Materials ###
        materials = openmc.Materials()

        fuel_mat = openmc.Material(name="fuel")
        fuel_mat.add_nuclide("U235", 1.0)
        fuel_mat.set_density('g/cc', 4.5)
        materials.append(fuel_mat)

        zirc_mat = openmc.Material(name="zircaloy")
        zirc_mat.add_nuclide("Zr90", 0.5145)
        zirc_mat.add_nuclide("Zr91", 0.1122)
        zirc_mat.add_nuclide("Zr92", 0.1715)
        zirc_mat.add_nuclide("Zr94", 0.1738)
        zirc_mat.add_nuclide("Zr96", 0.028)
        zirc_mat.set_density("g/cc", 5.77)
        materials.append(zirc_mat)

        water_mat = openmc.Material(name="water")
        water_mat.add_nuclide("H1", 2.0)
        water_mat.add_nuclide("O16", 1.0)
        water_mat.set_density("atom/b-cm", 0.07416)
        materials.append(water_mat)

        materials.export_to_xml()

        ### Geometry ###
        fuel_min_x = openmc.XPlane(x0=-5.0, name="minimum x")
        fuel_max_x = openmc.XPlane(x0=5.0, name="maximum x")

        fuel_min_y = openmc.YPlane(y0=-5.0, name="minimum y")
        fuel_max_y = openmc.YPlane(y0=5.0, name="maximum y")

        fuel_min_z = openmc.ZPlane(z0=-5.0, name="minimum z")
        fuel_max_z = openmc.ZPlane(z0=5.0, name="maximum z")

        fuel_cell = openmc.Cell(name="fuel")
        fuel_cell.region = +fuel_min_x & -fuel_max_x & \
            +fuel_min_y & -fuel_max_y & \
            +fuel_min_z & -fuel_max_z
        fuel_cell.fill = fuel_mat

        clad_min_x = openmc.XPlane(x0=-6.0, name="minimum x")
        clad_max_x = openmc.XPlane(x0=6.0, name="maximum x")

        clad_min_y = openmc.YPlane(y0=-6.0, name="minimum y")
        clad_max_y = openmc.YPlane(y0=6.0, name="maximum y")

        clad_min_z = openmc.ZPlane(z0=-6.0, name="minimum z")
        clad_max_z = openmc.ZPlane(z0=6.0, name="maximum z")

        clad_cell = openmc.Cell(name="clad")
        clad_cell.region = (-fuel_min_x | +fuel_max_x | \
                            -fuel_min_y | +fuel_max_y | \
                            -fuel_min_z | +fuel_max_z) & \
                            (+clad_min_x & -clad_max_x & \
                             +clad_min_y & -clad_max_y & \
                             +clad_min_z & -clad_max_z)
        clad_cell.fill = zirc_mat

        if self.external_geom:
            bounds = (15, 15, 15)
        else:
            bounds = (10, 10, 10)

        water_min_x = openmc.XPlane(x0=-bounds[0],
                                    name="minimum x",
                                    boundary_type='vacuum')
        water_max_x = openmc.XPlane(x0=bounds[0],
                                    name="maximum x",
                                    boundary_type='vacuum')

        water_min_y = openmc.YPlane(y0=-bounds[1],
                                    name="minimum y",
                                    boundary_type='vacuum')
        water_max_y = openmc.YPlane(y0=bounds[1],
                                    name="maximum y",
                                    boundary_type='vacuum')

        water_min_z = openmc.ZPlane(z0=-bounds[2],
                                    name="minimum z",
                                    boundary_type='vacuum')
        water_max_z = openmc.ZPlane(z0=bounds[2],
                                    name="maximum z",
                                    boundary_type='vacuum')

        water_cell = openmc.Cell(name="water")
        water_cell.region = (-clad_min_x | +clad_max_x | \
                             -clad_min_y | +clad_max_y | \
                             -clad_min_z | +clad_max_z) & \
                             (+water_min_x & -water_max_x & \
                              +water_min_y & -water_max_y & \
                              +water_min_z & -water_max_z)
        water_cell.fill = water_mat

        # create a containing universe
        root_univ = openmc.Universe()
        root_univ.add_cells([fuel_cell, clad_cell, water_cell])

        geom = openmc.Geometry(root=root_univ)

        geom.export_to_xml()

        ### Tallies ###

        # create meshes
        regular_mesh = openmc.RegularMesh()
        regular_mesh.dimension = (10, 10, 10)
        regular_mesh.lower_left = (-10.0, -10.0, -10.0)
        regular_mesh.upper_right = (10.0, 10.0, 10.0)

        regular_mesh_filter = openmc.MeshFilter(mesh=regular_mesh)

        uscd_mesh = openmc.UnstructuredMesh()
        uscd_mesh.filename = self.mesh_filename
        uscd_mesh.mesh_lib = 'moab'
        uscd_filter = openmc.MeshFilter(mesh=uscd_mesh)

        # create tallies
        tallies = openmc.Tallies()

        regular_mesh_tally = openmc.Tally(name="regular mesh tally")
        regular_mesh_tally.filters = [regular_mesh_filter]
        regular_mesh_tally.scores = ['flux']
        regular_mesh_tally.estimator = self.estimator
        tallies.append(regular_mesh_tally)

        uscd_tally = openmc.Tally(name="unstructured mesh tally")
        uscd_tally.filters = [uscd_filter]
        uscd_tally.scores = ['flux']
        uscd_tally.estimator = self.estimator
        tallies.append(uscd_tally)

        tallies.export_to_xml()

        ### Settings ###
        settings = openmc.Settings()
        settings.run_mode = 'fixed source'
        settings.particles = 100
        settings.batches = 10

        # source setup
        r = openmc.stats.Uniform(a=0.0, b=0.0)
        theta = openmc.stats.Discrete(x=[0.0], p=[1.0])
        phi = openmc.stats.Discrete(x=[0.0], p=[1.0])
        origin = (0.0, 0.0, 0.0)

        space = openmc.stats.SphericalIndependent(r=r,
                                                  theta=theta,
                                                  phi=phi,
                                                  origin=origin)

        angle = openmc.stats.Monodirectional((-1.0, 0.0, 0.0))

        energy = openmc.stats.Discrete(x=[15.0E6], p=[1.0])

        source = openmc.Source(space=space, energy=energy, angle=angle)

        settings.source = source

        settings.export_to_xml()

    def _compare_inputs(self):
        pass

    def _compare_results(self):
        with openmc.StatePoint(self._sp_name) as sp:
            # loop over the tallies
            regular_data = None
            regular_std_dev = None
            unstructured_data = None
            unstructured_std_dev = None
            for tally in sp.tallies.values():
                # find the regular and unstructured meshes
                if tally.contains_filter(openmc.MeshFilter):
                    flt = tally.find_filter(openmc.MeshFilter)

                    if isinstance(flt.mesh, openmc.RegularMesh):
                        regular_data, regular_std_dev = self.get_mesh_tally_data(tally)
                        if self.holes:
                            regular_data = np.delete(regular_data, holes)
                            regular_std_dev = np.delete(regular_std_dev, holes)
                    else:
                        unstructured_data, unstructured_std_dev = self.get_mesh_tally_data(tally, True)

            # successively check how many decimals the results are equal to
            decimals = 1
            while True:
                try:
                    np.testing.assert_array_almost_equal(unstructured_data,
                                                         regular_data,
                                                         decimals)
                except AssertionError as ae:
                    print(ae)
                    print()
                    break
                # increment decimals
                decimals += 1

                print("Results equal to within {} decimal places.\n".format(decimals))

            assert decimals >= 10

    @staticmethod
    def get_mesh_tally_data(tally, structured=False):
        TETS_PER_VOXEL = 12
        data = tally.get_reshaped_data(value='mean')
        std_dev = tally.get_reshaped_data(value='std_dev')
        if structured:
           data.shape = (data.size // TETS_PER_VOXEL, TETS_PER_VOXEL)
           std_dev.shape = (std_dev.size // TETS_PER_VOXEL, TETS_PER_VOXEL)
        else:
            data.shape = (data.size, 1)
            std_dev.shape = (std_dev.size, 1)
        return np.sum(data, axis=1), np.sum(std_dev, axis=1)

    def _cleanup(self):
        super()._cleanup()

        output = glob.glob('tally*.h5m')
        for f in output:
            if os.path.exists(f):
                os.remove(f)


hole_indexes = (333, 90, 777)
param_values = (('collision', 'tracklength'), (True, False), (hole_indexes, ()))
test_cases = []
for estimator, holes, ext_geom in product(*param_values):
    test_cases.append({'estimator' : estimator,
                       'holes' : holes,
                       'external_geom' : ext_geom})
@pytest.mark.parametrize("opts", test_cases)
def test_unstructured_mesh(opts):
    harness = UnstructuredMeshTest('statepoint.10.h5', kwargs=opts)
    harness.main()
