#!/usr/bin/env python

from __future__ import print_function

import os
import shutil
import re
from subprocess import call 
from collections import OrderedDict
from optparse import OptionParser

# Compiler paths
FC_DEFAULT='gfortran'
MPI_DIR='/opt/mpich/3.0.4-gnu'
HDF5_DIR='/opt/hdf5/1.8.12-gnu'
PHDF5_DIR='/opt/phdf5/1.8.12-gnu'
PETSC_DIR='/opt/petsc/3.4.3-gnu'

# Command line parsing
parser = OptionParser()
parser.add_option('-j', '--parallel', dest='n_procs', default='1',
                  help="Number of parallel jobs.")
parser.add_option('-R', '--tests-regex', dest='regex_tests',
                  help="Run tests matching regular expression. \
                  Test names are the directories present in tests folder.\
                  This uses standard regex syntax to select tests.")
parser.add_option('-C', '--build-config', dest='build_config',
                  help="Build configurations matching regular expression. \
                        Specific build configurations can be printed out with \
                        optional argument -p, --print. This uses standard \
                        regex syntax to select build configurations.")
parser.add_option('-p', '--print', action="store_true", 
                  dest="print_build_configs", default=False,
                  help="Print out build configurations.")
parser.add_option("-b", "--branch", dest="branch", default="",
                  help="branch name for build")
parser.add_option("-D", "--dashboard", dest="dash", default="Experimental",
                  help="Dash name -- Experimental, Nightly, Continuous")
(options, args) = parser.parse_args()

# CTest script template
ctest_str = """set (CTEST_SOURCE_DIRECTORY "{source_dir}")
set (CTEST_BINARY_DIRECTORY "{build_dir}")

set(CTEST_SITE "{host_name}")
set (CTEST_BUILD_NAME "{build_name}")
set (CTEST_CMAKE_GENERATOR "Unix Makefiles")
set (CTEST_BUILD_OPTIONS "{build_opts}")

set(CTEST_UPDATE_COMMAND "git")

set(CTEST_CONFIGURE_COMMAND "${{CMAKE_COMMAND}} -H${{CTEST_SOURCE_DIRECTORY}} -B${{CTEST_BINARY_DIRECTORY}} ${{CTEST_BUILD_OPTIONS}}")
set(CTEST_MEMORYCHECK_COMMAND "/usr/bin/valgrind")
set(CTEST_MEMORYCHECK_COMMAND_OPTIONS "--tool=memcheck --leak-check=yes --show-reachable=yes --num-callers=20 --track-fds=yes")
set(MEM_CHECK {mem_check})
set(ENV{{MEM_CHECK}} ${{MEM_CHECK}})

set(CTEST_COVERAGE_COMMAND "/usr/bin/gcov")
set(COVERAGE {coverage})
set(ENV{{COVERAGE}} ${{COVERAGE}})

ctest_start("{dashboard}")
ctest_configure()
ctest_update()
ctest_build()
if(MEM_CHECK)
ctest_test({tests} PARALLEL_LEVEL {n_procs})
ctest_memcheck({tests})
else(MEM_CHECK)
ctest_test({tests} PARALLEL_LEVEL {n_procs})
endif(MEM_CHECK)
if(COVERAGE)
ctest_coverage()
endif(COVERAGE)
ctest_submit()"""

# Define test data structure
tests = OrderedDict()

class Test(object):
    def __init__(self, name, debug=False, optimize=False, mpi=False, openmp=False,
                 hdf5=False, petsc=False, valgrind=False, coverage=False):
        self.name = name
        self.debug = debug
        self.optimize = optimize
        self.mpi = mpi
        self.openmp = openmp
        self.hdf5 = hdf5
        self.petsc = petsc
        self.valgrind = valgrind
        self.coverage = coverage

        # Check for MPI/HDF5
        if self.mpi and not self.hdf5:
            self.fc = MPI_DIR+'/bin/mpif90'
        elif not self.mpi and self.hdf5:
            self.fc = HDF5_DIR+'/bin/h5fc'
        elif self.mpi and self.hdf5:
            self.fc = PHDF5_DIR+'/bin/h5pfc'
        else:
            self.fc = FC_DEFAULT

    def get_build_name(self):
        self.build_name =  options.branch + '_' + self.name
        return self.build_name

    def get_build_opts(self):
        build_str = ""
        if self.debug:
            build_str += "-Ddebug=ON "
        if self.optimize:
            build_str += "-Doptimize=ON "
        if self.openmp:
            build_str += "-Dopenmp=ON "
        if self.petsc:
            build_str += "-Dpetsc=ON "
        if self.coverage:
            build_str += "-Dcoverage=ON "
        self.build_opts = build_str
        return self.build_opts

    def create_ctest_script(self, ctest_vars):
        with open('ctestscript.run', 'w') as fh:
            fh.write(ctest_str.format(**ctest_vars))

    def run_ctest(self):
        os.environ['FC'] = self.fc
        if self.petsc:
            os.environ['PETSC_DIR'] = PETSC_DIR
        if self.mpi:
            os.environ['MPI_DIR'] = MPI_DIR
        call(['ctest', '-S', 'ctestscript.run','-V'])

def add_test(name, debug=False, optimize=False, mpi=False, openmp=False,\
             hdf5=False, petsc=False, valgrind=False, coverage=False):
    tests.update({name:Test(name, debug, optimize, mpi, openmp, hdf5, petsc, valgrind, coverage)})

# List of tests
add_test('basic-normal')
add_test('basic-debug', debug=True)
add_test('basic-optimize', optimize=True)
add_test('omp-normal', openmp=True)
add_test('omp-debug', openmp=True, debug=True)
add_test('omp-optimize', openmp=True, optimize=True)
add_test('hdf5-normal', hdf5=True)
add_test('hdf5-debug', hdf5=True, debug=True)
add_test('hdf5-optimize', hdf5=True, optimize=True)
add_test('omp-hdf5-normal', openmp=True, hdf5=True)
add_test('omp-hdf5-debug', openmp=True, hdf5=True, debug=True)
add_test('omp-hdf5-optimize', openmp=True, hdf5=True, optimize=True)
add_test('mpi-normal', mpi=True)
add_test('mpi-debug', mpi=True, debug=True)
add_test('mpi-optimize', mpi=True, optimize=True)
add_test('mpi-omp-normal', mpi=True, openmp=True)
add_test('mpi-omp-debug', mpi=True, openmp=True, debug=True)
add_test('mpi-omp-optimize', mpi=True, openmp=True, optimize=True)
add_test('phdf5-normal', mpi=True, hdf5=True)
add_test('phdf5-debug', mpi=True, hdf5=True, debug=True)
add_test('phdf5-optimize', mpi=True, hdf5=True, optimize=True)
add_test('phdf5-omp-normal', mpi=True, hdf5=True, openmp=True)
add_test('phdf5-omp-debug', mpi=True, hdf5=True, openmp=True, debug=True)
add_test('phdf5-omp-optimize', mpi=True, hdf5=True, openmp=True, optimize=True)
add_test('petsc-normal', petsc=True, mpi=True)
add_test('petsc-debug', petsc=True, mpi=True, debug=True)
add_test('petsc-optimize', petsc=True, mpi=True, optimize=True)
add_test('phdf5-petsc-normal', mpi=True, hdf5=True, petsc=True)
add_test('phdf5-petsc-debug', mpi=True, hdf5=True, petsc=True, debug=True)
add_test('phdf5-petsc-optimize', mpi=True, hdf5=True, petsc=True, optimize=True)
add_test('omp-phdf5-petsc-normal', openmp=True, mpi=True, hdf5=True, petsc=True)
add_test('omp-phdf5-petsc-debug', openmp=True, mpi=True, hdf5=True, petsc=True, debug=True)
add_test('omp-phdf5-petsc-optimize', openmp=True, mpi=True, hdf5=True, petsc=True, optimize=True)
add_test('basic-debug_valgrind', debug=True, valgrind=True)
add_test('hdf5-debug_valgrind', hdf5=True, debug=True, valgrind=True)
add_test('petsc-debug_valgrind', petsc=True, mpi=True, debug=True, valgrind=True)
add_test('basic-debug_coverage', debug=True, coverage=True)
add_test('hdf5-debug_coverage', debug=True, hdf5=True, coverage=True)
add_test('mpi-debug_coverage', debug=True, mpi=True, coverage=True)
add_test('petsc-debug_coverage', debug=True, petsc=True, mpi=True, coverage=True)

# Check to see if we are to just print build configuratinos
if options.print_build_configs:
    for key in tests:
        print('Configuration Name: {0}'.format(key))
        print('  Debug Flags:..........{0}'.format(tests[key].debug))
        print('  Optimization Flags:...{0}'.format(tests[key].optimize))
        print('  HDF5 Active:..........{0}'.format(tests[key].hdf5))
        print('  MPI Active:...........{0}'.format(tests[key].mpi))
        print('  OpenMP Active:........{0}'.format(tests[key].openmp))
        print('  PETSc Active:.........{0}'.format(tests[key].petsc))
        print('  Valgrind Test:........{0}'.format(tests[key].valgrind))
        print('  Coverage Test:........{0}\n'.format(tests[key].coverage))
    exit()

# Delete items of dictionary that don't match regular expression
if options.build_config is not None:
    for key in tests:
        if not re.search(options.build_config, key):
            del tests[key]

# Setup CTest vars
pwd = os.environ['PWD']
ctest_vars = {
'source_dir' : pwd + '/../src',
'build_dir' :  pwd + '/build',
'host_name' : 'neutronbalance',
'dashboard' : options.dash,
'n_procs'   : options.n_procs
}

# Set up default valgrind tests
# Currently takes too long to run all the tests with valgrind
valgrind_default_tests = "basic|cmfd_feed|confidence_intervals| \
    density_atombcm|eigenvalue_genperbatch|energy_grid|entropy| \
    filter_cell|lattice_multiple|output|plot_background|reflective_plane| \
    rotation|salphabeta_multiple|score_absorption|seed|source_energy_mono| \
    sourcepoint_batch|statepoint_interval|survival_biasing| \
    tally_assumesep|translation|uniform_fs|universe|void"

# Begin testing
shutils.rmtree('build', ignore_errors=True)
os.remove('ctestscript.run')
call(['./cleanup'])
for key in iter(tests):
    test = tests[key]

    # Set test specific CTest vars
    ctest_vars.update({'build_name' : test.get_build_name()})
    ctest_vars.update({'build_opts' : test.get_build_opts()})
    ctest_vars.update({'mem_check'  : test.valgrind})
    ctest_vars.update({'coverage'  : test.coverage})

    # Check for user custom tests
    # INCLUDE is a CTest command that allows for a subset
    # of tests to be executed.
    if options.regex_tests is None:
        ctest_vars.update({'tests' : ''})

        # No user tests, use default valgrind tests
        if test.valgrind:
            ctest_vars.update({'tests' : 'INCLUDE {0}'.
                              format(valgrind_default_tests)})
    else:
        ctest_vars.update({'tests' : 'INCLUDE {0}'.
                          format(options.regex_tests)})

    # Create ctest script
    test.create_ctest_script(ctest_vars)

    # Run test
    test.run_ctest()

    # Clear build directory
    shutils.rmtree('build', ignore_errors=True)
    os.remove('ctestscript.run')
    call(['./cleanup'])
