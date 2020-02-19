import openmc
import pytest
import os
import glob

from tests.testing_harness import PyAPITestHarness

def __write_cmake_file(openmc_dir, build_dir):
    cmake_string = "cmake_minimum_required(VERSION 3.3 FATAL_ERROR)\n" +
                   "project(openmc_sources CXX)\n" +
                   "add_library(source SHARED source_ring.cpp)\n" +
                   "find_package(OpenMC REQUIRED HINTS {})\n" +
                   "target_link_libraries(source OpenMC::libopenmc)"

    f = open('CMakeLists.txt','w')
    f.write(cmake_string.format(openmc_dir))
    f.close()

    return

# compile the external source
def compile_source():

    openmc_home_dir = os.environ['OPENMC_INSTALL_DIR']
    __write_cmake_file(openmc_home_dir,'build')

    # copy the cmakefile
    status = os.system('sh clear_and_build.sh')
    assert os.WEXITSTATUS(status) == 0

    return 0

class SourceTestHarness(PyAPITestHarness):
    # build the test geometry
    def _build_inputs(self):
        mats = openmc.Materials()

        natural_lead = openmc.Material(1, "natural_lead")
        natural_lead.add_element('Pb', 1,'ao')
        natural_lead.set_density('g/cm3', 11.34)
        mats.append(natural_lead)

        # surfaces
        surface_sph1 = openmc.Sphere(r=100,boundary_type='vacuum') 
        volume_sph1 = -surface_sph1 
        
        # cell
        cell_1 = openmc.Cell(region=volume_sph1)
        cell_1.fill = natural_lead #assigning a material to a cell
        universe = openmc.Universe(cells=[cell_1]) #hint, this list will need to include the new cell
        geom = openmc.Geometry(universe)

        # settings
        sett = openmc.Settings()
        batches = 10
        sett.batches = batches
        sett.inactive = 0
        sett.particles = 1000
        sett.particle = "neutron"
        sett.run_mode = 'fixed source'

        #source
        source = openmc.Source()
        source.library = 'build/libsource.so'
        sett.source = source

        # run
        model = openmc.model.Model(geom,mats,sett)
        model.export_to_xml()
        return 

def test_dlopen_source():
    compile_source()    
    harness = SourceTestHarness('statepoint.10.h5')    
    harness.main()
