#!/usr/bin/env python

from __future__ import print_function

import os
import sys
import shutil
import re
import glob
import socket
from subprocess import call 
from collections import OrderedDict
from optparse import OptionParser

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
parser.add_option('-l', '--list', action="store_true", 
                  dest="list_build_configs", default=False,
                  help="List out build configurations.")
parser.add_option("-p", "--project", dest="project", default="",
                  help="project name for build")
parser.add_option("-D", "--dashboard", dest="dash",
                  help="Dash name -- Experimental, Nightly, Continuous")
parser.add_option("-u", "--update", action="store_true", dest="update",
                  help="Allow CTest to update repo. (WARNING: may overwrite\
                        changes that were not pushed.")
parser.add_option("-s", "--script", action="store_true", dest="script",
                  help="Activate CTest scripting mode for coverage, valgrind\
                        and dashboard capability.")
(options, args) = parser.parse_args()

# Default compiler paths
FC='gfortran'
MPI_DIR='/opt/mpich/3.1-gnu'
HDF5_DIR='/opt/hdf5/1.8.12-gnu'
PHDF5_DIR='/opt/phdf5/1.8.12-gnu'
PETSC_DIR='/opt/petsc/3.4.4-gnu'

# Script mode for extra capability
script_mode = False

# Override default compiler paths if environmental vars are found
if 'FC' in os.environ:
    FC = os.environ['FC']
if 'MPI_DIR' in os.environ:
    MPI_DIR = os.environ['MPI_DIR']
if 'HDF5_DIR' in os.environ:
    HDF5_DIR = os.environ['HDF5_DIR']
if 'PHDF5_DIR' in os.environ:
    PHDF5_DIR = os.environ['PHDF5_DIR']
if 'PETSC_DIR' in os.environ:
    PETSC_DIR = os.environ['PETSC_DIR']

# CTest script template
ctest_str = """set (CTEST_SOURCE_DIRECTORY "{source_dir}")
set (CTEST_BINARY_DIRECTORY "{build_dir}")

set(CTEST_SITE "{host_name}")
set (CTEST_BUILD_NAME "{build_name}")
set (CTEST_CMAKE_GENERATOR "Unix Makefiles")
set (CTEST_BUILD_OPTIONS "{build_opts}")

set(CTEST_UPDATE_COMMAND "git")

set(CTEST_CONFIGURE_COMMAND "${{CMAKE_COMMAND}} -H${{CTEST_SOURCE_DIRECTORY}} -B${{CTEST_BINARY_DIRECTORY}} ${{CTEST_BUILD_OPTIONS}}")
set(CTEST_MEMORYCHECK_COMMAND "{valgrind_cmd}")
set(CTEST_MEMORYCHECK_COMMAND_OPTIONS "--tool=memcheck --leak-check=yes --show-reachable=yes --num-callers=20 --track-fds=yes")
set(CTEST_MEMORYCHECK_SUPPRESSIONS_FILE ${{CTEST_SOURCE_DIRECTORY}}/../tests/valgrind.supp)
set(MEM_CHECK {mem_check})
set(ENV{{MEM_CHECK}} ${{MEM_CHECK}})

set(CTEST_COVERAGE_COMMAND "{gcov_cmd}")
set(COVERAGE {coverage})
set(ENV{{COVERAGE}} ${{COVERAGE}})

{subproject}

ctest_start("{dashboard}")
ctest_configure()
{update}
ctest_build()
ctest_test({tests} PARALLEL_LEVEL {n_procs})
if(MEM_CHECK)
ctest_memcheck({tests})
endif(MEM_CHECK)
if(COVERAGE)
ctest_coverage()
endif(COVERAGE)
{submit}"""

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
        self.success = True
        self.msg = None
        self.skipped = False
        self.valgrind_cmd = "" 
        self.gcov_cmd = ""
        self.cmake = ['cmake', '-H../src', '-Bbuild']

        # Check for MPI/HDF5
        if self.mpi and not self.hdf5:
            self.fc = MPI_DIR+'/bin/mpif90'
        elif not self.mpi and self.hdf5:
            self.fc = HDF5_DIR+'/bin/h5fc'
        elif self.mpi and self.hdf5:
            self.fc = PHDF5_DIR+'/bin/h5pfc'
        else:
            self.fc = FC

    # Sets the build name that will show up on the CDash
    def get_build_name(self):
        self.build_name =  options.project + '_' + self.name
        return self.build_name

    # Sets up build options for various tests. It is used both
    # in script and non-script modes
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

    # Write out the ctest script to tests directory
    def create_ctest_script(self, ctest_vars):
        with open('ctestscript.run', 'w') as fh:
            fh.write(ctest_str.format(**ctest_vars))

    # Runs the ctest script which performs all the cmake/ctest/cdash
    def run_ctest_script(self):
        os.environ['FC'] = self.fc
        if self.petsc:
            os.environ['PETSC_DIR'] = PETSC_DIR
        if self.mpi:
            os.environ['MPI_DIR'] = MPI_DIR
        rc = call(['ctest', '-S', 'ctestscript.run','-V'])
        if rc != 0:
            self.success = False
            self.msg = 'Failed on ctest script.'

    # Runs cmake when in non-script mode
    def run_cmake(self):
        os.environ['FC'] = self.fc
        build_opts = self.build_opts.split()
        self.cmake += build_opts
        rc = call(self.cmake)
        if rc != 0:
            self.success = False
            self.msg = 'Failed on cmake.'

    # Runs make when in non-script mode
    def run_make(self):
        if not self.success:
            return

        # Default make string
        make_list = ['make','-s']

        # Check for parallel
        if options.n_procs is not None:
            make_list.append('-j')
            make_list.append(options.n_procs)

        # Run make
        rc = call(make_list)
        if rc != 0:
            self.success = False
            self.msg = 'Failed on make.'

    # Runs ctest when in non-script mode
    def run_ctests(self):
        if not self.success:
            return

        # Default ctest string
        ctest_list = ['ctest']

        # Check for parallel
        if options.n_procs is not None:
            ctest_list.append('-j')
            ctest_list.append(options.n_procs)

        # Check for subset of tests
        if options.regex_tests is not None:
            ctest_list.append('-R')
            ctest_list.append(options.regex_tests)

        # Run ctests
        rc = call(ctest_list)
        if rc != 0:
            self.success = False
            self.msg = 'Failed on testing.'

    # Checks to see if file exists in PWD or PATH
    def check_compiler(self):
        result = False
        if os.path.isfile(self.fc):
            result = True
        for path in os.environ["PATH"].split(":"):
            if os.path.isfile(os.path.join(path, self.fc)):
                result = True
        if not result:
            self.msg = 'Compiler not found: {0}'.\
                       format((os.path.join(path, self.fc)))
            self.success = False

    # Get valgrind command from user's environment
    def find_valgrind(self):
        result = False
        for path in os.environ["PATH"].split(":"):
            if os.path.isfile(os.path.join(path, 'valgrind')):
                self.valgrind_cmd = os.path.join(path, 'valgrind')
                result = True
                break
        if not result:
            self.msg = 'valgrind not found.'
            self.success = False

    # Get coverage command from user's environment
    def find_coverage(self):
        result = False
        for path in os.environ["PATH"].split(":"):
            if os.path.isfile(os.path.join(path, 'gcov')):
                self.gcov_cmd = os.path.join(path, 'gcov')
                result = True
                break
        if not result:
            self.msg = 'gcov not found.'
            self.success = False

# Simple function to add a test to the global tests dictionary
def add_test(name, debug=False, optimize=False, mpi=False, openmp=False,\
             hdf5=False, petsc=False, valgrind=False, coverage=False):
    tests.update({name:Test(name, debug, optimize, mpi, openmp, hdf5, petsc, 
    valgrind, coverage)})

# List of all tests that may be run. User can add -C to command line to specify
# a subset of these configurations
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

# Check to see if we should just print build configuration information to user
if options.list_build_configs:
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

# Check for dashboard and determine whether to push results to server
# Note that there are only 3 basic dashboards: 
# Experimental, Nightly, Continuous. On the CDash end, these can be
# reorganized into groups when a hostname, dashboard and build name
# are matched.
if options.dash is None:
    dash = 'Experimental'
    submit = ''
else:
    dash = options.dash
    submit = 'ctest_submit()'

# Check for update command, which will run git fetch/merge and will delete
# any changes to repo that were not pushed to remote origin
if options.update:
    update = 'ctest_update()'
else:
    update = ''

# Check for CTest scipts mode
# Sets up whether we should use just the basic ctest command or use 
# CTest scripting to perform tests.
if not options.dash is None or options.script:
    script_mode = True
else:
    script_mode = False

# Setup CTest script vars. Not used in non-script mode
pwd = os.environ['PWD']
ctest_vars = {
'source_dir' : pwd + '/../src',
'build_dir' :  pwd + '/build',
'host_name' : socket.gethostname(),
'dashboard' : dash,
'submit'    : submit,
'update'    : update,
'n_procs'   : options.n_procs
}

# Check project name
subprop = """set_property(GLOBAL PROPERTY SubProject {0})"""
if options.project == "" :
    ctest_vars.update({'subproject':''})
elif options.project == 'develop':
    ctest_vars.update({'subproject':''})
else:
    ctest_vars.update({'subproject':subprop.format(options.project)})

# Set up default valgrind tests (subset of all tests)
# Currently takes too long to run all the tests with valgrind
# Only used in script mode
valgrind_default_tests = "basic|cmfd_feed|confidence_intervals|\
density_atombcm|eigenvalue_genperbatch|energy_grid|entropy|\
filter_cell|lattice_multiple|output|plot_background|reflective_plane|\
rotation|salphabeta_multiple|score_absorption|seed|source_energy_mono|\
sourcepoint_batch|statepoint_interval|survival_biasing|\
tally_assumesep|translation|uniform_fs|universe|void"

# Delete items of dictionary if valgrind or coverage and not in script mode
if not script_mode:
    for key in tests:
        if re.search('valgrind|coverage', key):
            del tests[key]

# Check if tests empty
if len(tests.keys()) == 0:
    print('No tests to run.')
    exit()

# Begin testing
shutil.rmtree('build', ignore_errors=True)
call(['./cleanup']) # removes all binary and hdf5 output files from tests
for key in iter(tests):
    test = tests[key]

    # Extra display if not in script mode
    if not script_mode:
        print('-'*(len(key) + 6))
        print(key + ' tests')
        print('-'*(len(key) + 6))
        sys.stdout.flush()

    # Verify fortran compiler exists
    test.check_compiler()
    if not test.success:
        continue

    # Get valgrind command 
    if test.valgrind:
        test.find_valgrind()
    if not test.success:
        continue

    # Get coverage command
    if test.coverage:
        test.find_coverage()
    if not test.success:
        continue
 
    # Set test specific CTest script vars. Not used in non-script mode
    ctest_vars.update({'build_name' : test.get_build_name()})
    ctest_vars.update({'build_opts' : test.get_build_opts()})
    ctest_vars.update({'mem_check'  : test.valgrind})
    ctest_vars.update({'coverage'  : test.coverage})
    ctest_vars.update({'valgrind_cmd'  : test.valgrind_cmd})
    ctest_vars.update({'gcov_cmd'  : test.gcov_cmd})

    # Check for user custom tests
    # INCLUDE is a CTest command that allows for a subset
    # of tests to be executed. Only used in script mode.
    if options.regex_tests is None:
        ctest_vars.update({'tests' : ''})

        # No user tests, use default valgrind tests
        if test.valgrind:
            ctest_vars.update({'tests' : 'INCLUDE {0}'.
                              format(valgrind_default_tests)})
    else:
        ctest_vars.update({'tests' : 'INCLUDE {0}'.
                          format(options.regex_tests)})

    # Main part of code that does the ctest execution.
    # It is broken up by two modes, script and non-script
    if script_mode:

        # Create ctest script
        test.create_ctest_script(ctest_vars)

        # Run test
        test.run_ctest_script()

    else:

        # Run CMAKE to configure build
        test.run_cmake()

        # Go into build directory
        os.chdir('build')

        # Build OpenMC
        test.run_make()

        # Run tests
        test.run_ctests()

        # Leave build directory
        os.chdir('..')

    # Copy over log file
    if script_mode:
        logfile = glob.glob('build/Testing/Temporary/LastTest_*.log')
    else:
        logfile = glob.glob('build/Testing/Temporary/LastTest.log')
    if len(logfile) > 0:
        logfilename = os.path.split(logfile[0])[1]
        logfilename = os.path.splitext(logfilename)[0]
        logfilename = logfilename + '_{0}.log'.format(test.name)
        shutil.copy(logfile[0], logfilename)

    # Clear build directory and remove binary and hdf5 files
    shutil.rmtree('build', ignore_errors=True)
    if script_mode:
        os.remove('ctestscript.run')
    call(['./cleanup'])

# Print out summary of results
print('\n' + '='*54)
print('Summary of Compilation Option Testing:\n')

if sys.stdout.isatty():
    OK = '\033[92m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
else:
    OK = ''
    FAIL = ''
    ENDC = ''
    BOLD = ''

for test in tests:
    print(test + '.'*(50 - len(test)), end='')
    if tests[test].success:
        print(BOLD + OK + '[OK]' + ENDC)
    else:
        print(BOLD + FAIL + '[FAILED]' + ENDC)
        print(' '*len(test)+tests[test].msg)
