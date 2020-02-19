from pathlib import Path
import os
import shutil
import subprocess
import textwrap

import openmc
import pytest

from tests.testing_harness import PyAPITestHarness


@pytest.fixture
def compile_source(request):
    """Compile the external source"""

    # Get build directory and write CMakeLists.txt file
    openmc_dir = Path(str(request.config.rootdir)) / 'build'
    with open('CMakeLists.txt', 'w') as f:
        f.write(textwrap.dedent("""
            cmake_minimum_required(VERSION 3.3 FATAL_ERROR)
            project(openmc_sources CXX)
            add_library(source SHARED source_sampling.cpp)
            find_package(OpenMC REQUIRED HINTS {})
            target_link_libraries(source OpenMC::libopenmc)
            """.format(openmc_dir)))

    # Create temporary build directory and change to there
    local_builddir = Path('build')
    local_builddir.mkdir(exist_ok=True)
    os.chdir(str(local_builddir))

    # Run cmake/make to build the shared libary
    subprocess.run(['cmake', os.path.pardir], check=True)
    subprocess.run(['make'], check=True)
    os.chdir(os.path.pardir)

    yield

    # Remove local build directory when test is complete
    shutil.rmtree('build')


@pytest.fixture
def model():
    model = openmc.model.Model()
    natural_lead = openmc.Material(name="natural_lead")
    natural_lead.add_element('Pb', 1.0)
    natural_lead.set_density('g/cm3', 11.34)
    model.materials.append(natural_lead)

    # geometry
    surface_sph1 = openmc.Sphere(r=100, boundary_type='vacuum')
    cell_1 = openmc.Cell(fill=natural_lead, region=-surface_sph1)
    model.geometry = openmc.Geometry([cell_1])

    # settings
    model.settings.batches = 10
    model.settings.inactive = 0
    model.settings.particles = 1000
    model.settings.run_mode = 'fixed source'

    # custom source from shared library
    source = openmc.Source()
    source.library = 'build/libsource.so'
    model.settings.source = source

    return model


def test_dlopen_source(compile_source, model):
    harness = PyAPITestHarness('statepoint.10.h5', model)
    harness.main()
