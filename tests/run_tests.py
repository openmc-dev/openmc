#!/usr/bin/env python

from __future__ import print_function

import os
import sys
import shutil
import re
import glob
import socket
from subprocess import call, check_output
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
CC='gcc'
CXX='g++'
MPI_DIR='/opt/mpich/3.2-gnu'
HDF5_DIR='/opt/hdf5/1.8.16-gnu'
PHDF5_DIR='/opt/phdf5/1.8.16-gnu'

# Script mode for extra capability
script_mode = False

# Override default compiler paths if environmental vars are found
if 'FC' in os.environ:
    FC = os.environ['FC']
if 'CC' in os.environ:
    CC = os.environ['CC']
if 'CXX' in os.environ:
    CXX = os.environ['CXX']
if 'MPI_DIR' in os.environ:
    MPI_DIR = os.environ['MPI_DIR']
if 'HDF5_DIR' in os.environ:
    HDF5_DIR = os.environ['HDF5_DIR']
if 'PHDF5_DIR' in os.environ:
    PHDF5_DIR = os.environ['PHDF5_DIR']

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
#set(CTEST_MEMORYCHECK_SUPPRESSIONS_FILE ${{CTEST_SOURCE_DIRECTORY}}/../tests/valgrind.supp)
set(MEM_CHECK {mem_check})
if(MEM_CHECK)
set(ENV{{MEM_CHECK}} ${{MEM_CHECK}})
endif()

set(CTEST_COVERAGE_COMMAND "gcov")
set(COVERAGE {coverage})
set(ENV{{COVERAGE}} ${{COVERAGE}})

{subproject}

ctest_start("{dashboard}")
ctest_configure(RETURN_VALUE res)
{update}
ctest_build(RETURN_VALUE res)
if(NOT MEM_CHECK)
ctest_test({tests} PARALLEL_LEVEL {n_procs}, RETURN_VALUE res)
endif()
if(MEM_CHECK)
ctest_memcheck({tests} RETURN_VALUE res)
endif(MEM_CHECK)
if(COVERAGE)
ctest_coverage(RETURN_VALUE res)
endif(COVERAGE)
{submit}

if (res EQUAL 0)
else()
message(FATAL_ERROR "")
endif()
"""

# Define test data structure
tests = OrderedDict()

def cleanup(path):
    """Remove generated output files."""
    for dirpath, dirnames, filenames in os.walk(path):
        for fname in filenames:
            for ext in ['.h5', '.ppm', '.voxel']:
                if fname.endswith(ext) and fname != '1d_mgxs.h5':
                    os.remove(os.path.join(dirpath, fname))


def which(program):
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file
    return None


class Test(object):
    def __init__(self, name, debug=False, optimize=False, mpi=False, openmp=False,
                 phdf5=False, valgrind=False, coverage=False):
        self.name = name
        self.debug = debug
        self.optimize = optimize
        self.mpi = mpi
        self.openmp = openmp
        self.phdf5 = phdf5
        self.valgrind = valgrind
        self.coverage = coverage
        self.success = True
        self.msg = None
        self.skipped = False
        self.cmake = ['cmake', '-H..', '-Bbuild',
                      '-DPYTHON_EXECUTABLE=' + sys.executable]

        # Check for MPI
        if self.mpi:
            if os.path.exists(os.path.join(MPI_DIR, 'bin', 'mpifort')):
                self.fc = os.path.join(MPI_DIR, 'bin', 'mpifort')
            else:
                self.fc = os.path.join(MPI_DIR, 'bin', 'mpif90')
            self.cc = os.path.join(MPI_DIR, 'bin', 'mpicc')
            self.cxx = os.path.join(MPI_DIR, 'bin', 'mpicxx')
        else:
            self.fc = FC
            self.cc = CC
            self.cxx = CXX

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
        if not self.openmp:
            build_str += "-Dopenmp=OFF "
        if self.coverage:
            build_str += "-Dcoverage=ON "
        if self.phdf5:
            build_str += "-DHDF5_PREFER_PARALLEL=ON "
        else:
            build_str += "-DHDF5_PREFER_PARALLEL=OFF "
        self.build_opts = build_str
        return self.build_opts

    # Write out the ctest script to tests directory
    def create_ctest_script(self, ctest_vars):
        with open('ctestscript.run', 'w') as fh:
            fh.write(ctest_str.format(**ctest_vars))

    # Runs the ctest script which performs all the cmake/ctest/cdash
    def run_ctest_script(self):
        os.environ['FC'] = self.fc
        os.environ['CC'] = self.cc
        os.environ['CXX'] = self.cxx
        if self.mpi:
            os.environ['MPI_DIR'] = MPI_DIR
        if self.phdf5:
            os.environ['HDF5_ROOT'] = PHDF5_DIR
        else:
            os.environ['HDF5_ROOT'] = HDF5_DIR
        rc = call(['ctest', '-S', 'ctestscript.run','-V'])
        if rc != 0:
            self.success = False
            self.msg = 'Failed on ctest script.'

    # Runs cmake when in non-script mode
    def run_cmake(self):
        build_opts = self.build_opts.split()
        self.cmake += build_opts

        os.environ['FC'] = self.fc
        os.environ['CC'] = self.cc
        os.environ['CXX'] = self.cxx
        if self.mpi:
            os.environ['MPI_DIR'] = MPI_DIR
        if self.phdf5:
            os.environ['HDF5_ROOT'] = PHDF5_DIR
            self.cmake.append('-DHDF5_PREFER_PARALLEL=ON')
        else:
            os.environ['HDF5_ROOT'] = HDF5_DIR
            self.cmake.append('-DHDF5_PREFER_PARALLEL=OFF')
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


# Simple function to add a test to the global tests dictionary
def add_test(name, debug=False, optimize=False, mpi=False, openmp=False,\
             phdf5=False, valgrind=False, coverage=False):
    tests.update({name: Test(name, debug, optimize, mpi, openmp, phdf5,
                             valgrind, coverage)})

# List of all tests that may be run. User can add -C to command line to specify
# a subset of these configurations
add_test('hdf5-normal')
add_test('hdf5-debug', debug=True)
add_test('hdf5-optimize', optimize=True)
add_test('omp-hdf5-normal', openmp=True)
add_test('omp-hdf5-debug', openmp=True, debug=True)
add_test('omp-hdf5-optimize', openmp=True, optimize=True)
add_test('mpi-hdf5-normal', mpi=True)
add_test('mpi-hdf5-debug', mpi=True, debug=True)
add_test('mpi-hdf5-optimize', mpi=True, optimize=True)
add_test('phdf5-normal', mpi=True, phdf5=True)
add_test('phdf5-debug', mpi=True, phdf5=True, debug=True)
add_test('phdf5-optimize', mpi=True, phdf5=True, optimize=True)
add_test('phdf5-omp-normal', mpi=True, phdf5=True, openmp=True)
add_test('phdf5-omp-debug', mpi=True, phdf5=True, openmp=True, debug=True)
add_test('phdf5-omp-optimize', mpi=True, phdf5=True, openmp=True, optimize=True)
add_test('hdf5-debug_valgrind', debug=True, valgrind=True)
add_test('hdf5-debug_coverage', debug=True, coverage=True)

# Check to see if we should just print build configuration information to user
if options.list_build_configs:
    for key in tests:
        print('Configuration Name: {0}'.format(key))
        print('  Debug Flags:..........{0}'.format(tests[key].debug))
        print('  Optimization Flags:...{0}'.format(tests[key].optimize))
        print('  MPI Active:...........{0}'.format(tests[key].mpi))
        print('  OpenMP Active:........{0}'.format(tests[key].openmp))
        print('  Valgrind Test:........{0}'.format(tests[key].valgrind))
        print('  Coverage Test:........{0}\n'.format(tests[key].coverage))
    exit()

# Delete items of dictionary that don't match regular expression
if options.build_config is not None:
    to_delete = []
    for key in tests:
        if not re.search(options.build_config, key):
            to_delete.append(key)
    for key in to_delete:
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
pwd = os.getcwd()
ctest_vars = {
    'source_dir': os.path.join(pwd, os.pardir),
    'build_dir': os.path.join(pwd, 'build'),
    'host_name': socket.gethostname(),
    'dashboard': dash,
    'submit': submit,
    'update': update,
    'n_procs': options.n_procs
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
valgrind_default_tests = "cmfd_feed|confidence_intervals|\
density|eigenvalue_genperbatch|energy_grid|entropy|\
lattice_multiple|output|plotreflective_plane|\
rotation|salphabetascore_absorption|seed|source_energy_mono|\
sourcepoint_batch|statepoint_interval|survival_biasing|\
tally_assumesep|translation|uniform_fs|universe|void"

# Delete items of dictionary if valgrind or coverage and not in script mode
to_delete = []
if not script_mode:
    for key in tests:
        if re.search('valgrind|coverage', key):
            to_delete.append(key)

for key in to_delete:
    del tests[key]

# Check if tests empty
if len(list(tests.keys())) == 0:
    print('No tests to run.')
    exit()

# Begin testing
shutil.rmtree('build', ignore_errors=True)
cleanup('.')
for key in iter(tests):
    test = tests[key]

    # Extra display if not in script mode
    if not script_mode:
        print('-'*(len(key) + 6))
        print(key + ' tests')
        print('-'*(len(key) + 6))
        sys.stdout.flush()

    # Verify fortran compiler exists
    if which(test.fc) is None:
        self.msg = 'Compiler not found: {0}'.format(test.fc)
        self.success = False
        continue

    # Verify valgrind command exists
    if test.valgrind:
        valgrind_cmd = which('valgrind')
        if valgrind_cmd is None:
            self.msg = 'No valgrind executable found.'
            self.success = False
            continue
    else:
        valgrind_cmd = ''

    # Verify gcov/lcov exist
    if test.coverage:
        if which('gcov') is None:
            self.msg = 'No {} executable found.'.format(exe)
            self.success = False
            continue

    # Set test specific CTest script vars. Not used in non-script mode
    ctest_vars.update({'build_name': test.get_build_name()})
    ctest_vars.update({'build_opts': test.get_build_opts()})
    ctest_vars.update({'mem_check': test.valgrind})
    ctest_vars.update({'coverage': test.coverage})
    ctest_vars.update({'valgrind_cmd': valgrind_cmd})

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
        os.chdir(os.pardir)

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

    # For coverage builds, use lcov to generate HTML output
    if test.coverage:
        if which('lcov') is None or which('genhtml') is None:
            print('No lcov/genhtml command found. '
                  'Could not generate coverage report.')
        else:
            shutil.rmtree('coverage', ignore_errors=True)
            call(['lcov', '--directory', '.', '--capture',
                  '--output-file', 'coverage.info'])
            call(['genhtml', '--output-directory', 'coverage', 'coverage.info'])
            os.remove('coverage.info')

    if test.valgrind:
        # Copy memcheck output to memcheck directory
        shutil.rmtree('memcheck', ignore_errors=True)
        os.mkdir('memcheck')
        memcheck_out = glob.glob('build/Testing/Temporary/MemoryChecker.*.log')
        for fname in memcheck_out:
            shutil.copy(fname, 'memcheck/')

        # Remove generated XML files
        xml_files = check_output(['git', 'ls-files', '.',  '--exclude-standard',
                                  '--others']).split()
        for f in xml_files:
            os.remove(f)

    # Clear build directory and remove binary and hdf5 files
    shutil.rmtree('build', ignore_errors=True)
    if script_mode:
        os.remove('ctestscript.run')
    cleanup('.')

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

return_code = 0

for test in tests:
    print(test + '.'*(50 - len(test)), end='')
    if tests[test].success:
        print(BOLD + OK + '[OK]' + ENDC)
    else:
        print(BOLD + FAIL + '[FAILED]' + ENDC)
        print(' '*len(test)+tests[test].msg)
        return_code = 1

sys.exit(return_code)
