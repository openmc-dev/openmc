from pathlib import Path
import os
import shutil
import subprocess
from subprocess import CalledProcessError
import textwrap
import glob
from itertools import product

import openmc
import openmc.lib
import numpy as np
import pytest

from tests.regression_tests import config
from tests.testing_harness import PyAPITestHarness

pytestmark = pytest.mark.skipif(
    not openmc.lib._dagmc_enabled(), reason="DAGMC is not enabled."
)

TETS_PER_VOXEL = 12

# Test that an external moab instance can be passed in through the C API


@pytest.fixture
def cpp_driver(request):
    """Compile the external source"""

    # Get build directory and write CMakeLists.txt file
    openmc_dir = Path(str(request.config.rootdir)) / "build"
    with open("CMakeLists.txt", "w") as f:
        f.write(
            textwrap.dedent(
                """
            cmake_minimum_required(VERSION 3.10 FATAL_ERROR)
            project(openmc_cpp_driver CXX)
            add_executable(main main.cpp)
            find_package(OpenMC REQUIRED HINTS {})
            target_link_libraries(main OpenMC::libopenmc)
            target_compile_features(main PUBLIC cxx_std_14)
            set(CMAKE_CXX_FLAGS "-pedantic-errors")
            add_compile_definitions(DAGMC=1)
            """.format(
                    openmc_dir
                )
            )
        )

    # Create temporary build directory and change to there
    local_builddir = Path("build")
    local_builddir.mkdir(exist_ok=True)
    os.chdir(str(local_builddir))

    if config["mpi"]:
        mpi_arg = "On"
    else:
        mpi_arg = "Off"

    try:
        print("Building driver")
        # Run cmake/make to build the shared libary
        subprocess.run(
            ["cmake", os.path.pardir, f"-DOPENMC_USE_MPI={mpi_arg}"], check=True
        )
        subprocess.run(["make"], check=True)
        os.chdir(os.path.pardir)

        yield "./build/main"

    finally:
        # Remove local build directory when test is complete
        shutil.rmtree("build")
        os.remove("CMakeLists.txt")


class ExternalMoabTest(PyAPITestHarness):
    def __init__(self, executable, statepoint_name, model):
        super().__init__(statepoint_name, model)
        self.executable = executable

    def _run_openmc(self):
        if config["mpi"]:
            mpi_args = [config["mpiexec"], "-n", config["mpi_np"]]
            openmc.run(
                openmc_exec=self.executable,
                mpi_args=mpi_args,
                event_based=config["event"],
            )
        else:
            openmc.run(openmc_exec=self.executable, event_based=config["event"])

    # Override some methods to do nothing
    def _get_results(self):
        pass

    def _write_results(self, results_string):
        pass

    def _overwrite_results(self):
        pass

    def _test_output_created(self):
        pass

    # Directly compare results of unstructured  mesh with internal and
    # external moab
    def _compare_results(self):

        with openmc.StatePoint(self._sp_name) as sp:
            # loop over the tallies and get data

            ext_data = []
            unstr_data = []

            for tally in sp.tallies.values():

                # Safety check that mesh filter is correct
                if tally.contains_filter(openmc.MeshFilter):
                    flt = tally.find_filter(openmc.MeshFilter)

                    if isinstance(flt.mesh, openmc.UnstructuredMesh):

                        if tally.name == "external mesh tally":
                            ext_data = tally.get_reshaped_data(value="mean")

                        elif tally.name == "unstructured mesh tally":
                            unstr_data = tally.get_reshaped_data(value="mean")

            # we expect these results to be the same to within at 8
            # decimal places
            decimals = 8
            np.testing.assert_array_almost_equal(unstr_data, ext_data, decimals)

    @staticmethod
    def get_mesh_tally_data(tally):
        data = tally.get_reshaped_data(value="mean")
        std_dev = tally.get_reshaped_data(value="std_dev")
        data.shape = (data.size, 1)
        std_dev.shape = (std_dev.size, 1)
        return np.sum(data, axis=1), np.sum(std_dev, axis=1)

    def _cleanup(self):
        super()._cleanup()
        output = glob.glob("tally*.vtk")
        for f in output:
            if os.path.exists(f):
                os.remove(f)


def test_external_mesh(cpp_driver):

    # Materials
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

    # Geometry
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

    bounds = (10, 10, 10)

    water_min_x = openmc.XPlane(x0=-bounds[0], name="minimum x", boundary_type="vacuum")
    water_max_x = openmc.XPlane(x0=bounds[0], name="maximum x", boundary_type="vacuum")

    water_min_y = openmc.YPlane(y0=-bounds[1], name="minimum y", boundary_type="vacuum")
    water_max_y = openmc.YPlane(y0=bounds[1], name="maximum y", boundary_type="vacuum")

    water_min_z = openmc.ZPlane(z0=-bounds[2], name="minimum z", boundary_type="vacuum")
    water_max_z = openmc.ZPlane(z0=bounds[2], name="maximum z", boundary_type="vacuum")

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
    geometry = openmc.Geometry([fuel_cell, clad_cell, water_cell])

    # Meshes
    mesh_filename = "test_mesh_tets.h5m"

    # Create a normal unstructured mesh to compare to
    uscd_mesh = openmc.UnstructuredMesh(mesh_filename, "moab")

    # Create filters
    uscd_filter = openmc.MeshFilter(mesh=uscd_mesh)

    # Tallies
    tallies = openmc.Tallies()
    uscd_tally = openmc.Tally(name="unstructured mesh tally")
    uscd_tally.filters = [uscd_filter]
    uscd_tally.scores = ["flux"]
    uscd_tally.estimator = "tracklength"
    tallies.append(uscd_tally)

    # Settings
    settings = openmc.Settings()
    settings.run_mode = "fixed source"
    settings.particles = 100
    settings.batches = 10

    # Source setup
    space = openmc.stats.Point()
    angle = openmc.stats.Monodirectional((-1.0, 0.0, 0.0))
    energy = openmc.stats.Discrete(x=[15.0e06], p=[1.0])
    source = openmc.IndependentSource(space=space, energy=energy, angle=angle)
    settings.source = source

    model = openmc.model.Model(
        geometry=geometry, materials=materials, tallies=tallies, settings=settings
    )

    harness = ExternalMoabTest(cpp_driver, "statepoint.10.h5", model)

    # Run open MC and check results
    harness.main()
