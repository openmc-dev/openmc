from pathlib import Path
import os
import shutil
import subprocess
import textwrap

from numpy.testing import assert_allclose
import openmc
import pytest

from tests.regression_tests import config
from tests.testing_harness import PyAPITestHarness


@pytest.fixture
def cpp_driver(request):
    """Compile the external source"""

    # Get build directory and write CMakeLists.txt file
    openmc_dir = Path(str(request.config.rootdir)) / 'build'
    with open('CMakeLists.txt', 'w') as f:
        f.write(textwrap.dedent("""
            cmake_minimum_required(VERSION 3.3 FATAL_ERROR)
            project(openmc_cpp_driver CXX)
            add_executable(cpp_driver driver.cpp)
            find_package(OpenMC REQUIRED HINTS {})
            target_link_libraries(cpp_driver OpenMC::libopenmc)
            """.format(openmc_dir)))

    # Create temporary build directory and change to there
    local_builddir = Path('build')
    local_builddir.mkdir(exist_ok=True)
    os.chdir(str(local_builddir))

    if config['mpi']:
        os.environ['CXX'] = 'mpicxx'

    try:
        print("Building driver")
        # Run cmake/make to build the shared libary
        subprocess.run(['cmake', os.path.pardir], check=True)
        subprocess.run(['make'], check=True)
        os.chdir(os.path.pardir)

        yield "./build/cpp_driver"

    finally:
        # Remove local build directory when test is complete
        shutil.rmtree('build')
        os.remove('CMakeLists.txt')


@pytest.fixture
def model():
    model = openmc.model.Model()

    # materials
    u235 = openmc.Material(name="fuel")
    u235.add_nuclide('U235', 1.0, 'ao')
    u235.set_density('g/cc', 11)

    zirc = openmc.Material(name='cladding')
    zirc.add_nuclide('Zr90', 1.0)
    zirc.set_density('g/cc', 6.44)

    water = openmc.Material(name="water")
    water.add_nuclide('H1', 2.0, 'ao')
    water.add_nuclide('O16', 1.0, 'ao')
    water.set_density('g/cc', 1.0)

    mats = openmc.Materials([u235, zirc, water])
    model.materials = mats

    # geometry
    fuel_or = openmc.ZCylinder(r=1.5)
    cladding_or = openmc.ZCylinder(r=1.7)

    fuel = openmc.Cell(fill=u235, region=-fuel_or)
    cladding = openmc.Cell(fill=zirc, region=+fuel_or & -cladding_or)
    moderator = openmc.Cell(fill=water, region=+cladding_or)

    pincell_univ = openmc.Universe(cells=[fuel, cladding, moderator])

    # lattice
    lattice = openmc.RectLattice()
    lattice.pitch = (4.0, 4.0)
    lattice.lower_left = (-4.0, -4.0)
    lattice.universes = [[pincell_univ, pincell_univ], [pincell_univ, pincell_univ]]
    lattice_region = openmc.model.rectangular_prism(8.0,
                                                    8.0,
                                                    boundary_type='reflective')
    lattice_cell = openmc.Cell(fill=lattice, region=lattice_region)

    model.geometry = openmc.Geometry([lattice_cell])

    model.settings.particles = 100
    model.settings.batches = 10
    model.settings.inactive = 1

    return model


class ExternalDriverTestHarness(PyAPITestHarness):

    def __init__(self, executable, statepoint_name, model=None):
        super().__init__(statepoint_name, model)
        self.executable = executable

    def _run_openmc(self):
        if config['mpi']:
            mpi_args = [config['mpiexec'], '-n', config['mpi_np']]
            openmc.run(openmc_exec=self.executable,
                       mpi_args=mpi_args,
                       event_based=config['event'])
        else:
            openmc.run(openmc_exec=self.executable,
                       event_based=config['event'])

    def _compare_results(self):
        super()._compare_results()

        # load the summary file
        summary = openmc.Summary('summary.h5')

        # get the summary cells
        cells = summary.geometry.get_all_cells()

        # for the 2 by 2 lattice, each cell should have 4
        # temperature values set to 300 K
        for cell in cells.values():
            if isinstance(cell.fill, openmc.Material):
                assert len(cell.temperature) == 4
                assert_allclose(cell.temperature, 300.0)


def test_cpp_driver(cpp_driver, model):
    harness = ExternalDriverTestHarness(cpp_driver, 'statepoint.10.h5', model)
    harness.main()
