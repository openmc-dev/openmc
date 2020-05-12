from pathlib import Path
import os
import shutil
import subprocess
import textwrap

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

    # Run cmake/make to build the shared libary
    subprocess.run(['cmake', os.path.pardir], check=True)
    subprocess.run(['make', 'VERBOSE=1'], check=True)
    os.chdir(os.path.pardir)

    yield "./build/cpp_driver"

    # Remove local build directory when test is complete
    shutil.rmtree('build')

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

def test_cpp_driver(cpp_driver):
    model = openmc.examples.pwr_pin_cell()
    model.settings.particles = 100
    model.settings.batches = 10
    model.settings.inactive = 1

    harness = ExternalDriverTestHarness(cpp_driver, 'statepoint.10.h5', model)
    harness.main()
