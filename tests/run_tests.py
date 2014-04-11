#!/usr/bin/env python

from __future__ import print_function

import os
import shutil
import re
import sys
from subprocess import call 
from collections import OrderedDict
from optparse import OptionParser

parser = OptionParser()
parser.add_option('-j', '--parallel', dest='n_procs',
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
(opts, args) = parser.parse_args()

# Default compiler paths
FC='gfortran'
MPI_DIR='/opt/mpich/3.1-gnu'
HDF5_DIR='/opt/hdf5/1.8.12-gnu'
PHDF5_DIR='/opt/phdf5/1.8.12-gnu'
PETSC_DIR='/opt/petsc/3.4.4-gnu'

# Override default compiler paths if environmental vars are found
if 'FC' in os.environ:
    FC = os.environ['FC']
    if FC is not 'gfortran':
        print('NOTE: Test suite only verifed for gfortran compiler.')
if 'MPI_DIR' in os.environ:
    MPI_DIR = os.environ['MPI_DIR']
if 'HDF5_DIR' in os.environ:
    HDF5_DIR = os.environ['HDF5_DIR']
if 'PHDF5_DIR' in os.environ:
    PHDF5_DIR = os.environ['PHDF5_DIR']
if 'PETSC_DIR' in os.environ:
    PETSC_DIR = os.environ['PETSC_DIR']

# Define test data structure
tests = OrderedDict()

class Test(object):
    def __init__(self, debug=False, optimize=False, mpi=False, openmp=False,
                 hdf5=False, petsc=False):
        self.debug = debug
        self.optimize = optimize
        self.mpi = mpi
        self.openmp = openmp
        self.hdf5 = hdf5
        self.petsc = petsc
        self.success = True
        self.msg = None
        self.setup_cmake()

    def setup_cmake(self):
        # Default cmake
        self.cmake = ['cmake','-H../src','-Bbuild']

        # Check for MPI/HDF5
        if self.mpi and not self.hdf5:
            self.fc = MPI_DIR+'/bin/mpif90'
        elif not self.mpi and self.hdf5:
            self.fc = HDF5_DIR+'/bin/h5fc'
        elif self.mpi and self.hdf5:
            self.fc = PHDF5_DIR+'/bin/h5pfc'
        else:
            self.fc = FC

        # Set rest of options
        if self.debug:
            self.cmake.append('-Ddebug=on')
        if self.optimize:
            self.cmake.append('-Doptimize=on')
        if self.openmp:
            self.cmake.append('-Dopenmp=on')
        if self.petsc:
            self.cmake.append('-Dpetsc=on')
            os.environ['PETSC_DIR'] = PETSC_DIR
        if self.mpi:
            os.environ['MPI_DIR'] = MPI_DIR

    def run_cmake(self):
        os.environ['FC'] = self.fc
        rc = call(self.cmake)
        if rc != 0:
            self.success = False
            self.msg = 'Failed on cmake.'

    def run_make(self):
        if not self.success:
            return

        # Default make string
        make_list = ['make','-s']

        # Check for parallel
        if opts.n_procs is not None:
            make_list.append('-j')
            make_list.append(opts.n_procs)

        # Run make
        rc = call(make_list)
        if rc != 0:
            self.success = False
            self.msg = 'Failed on make.'

    def run_ctests(self):
        if not self.success:
            return

        # Default ctest string
        ctest_list = ['ctest']

        # Check for parallel
        if opts.n_procs is not None:
            ctest_list.append('-j')
            ctest_list.append(opts.n_procs)

        # Check for subset of tests
        if opts.regex_tests is not None:
            ctest_list.append('-R')
            ctest_list.append(opts.regex_tests)

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
            raise Exception("Compiler path '{0}' does not exist."
                           .format(self.fc)+
                           "Please set appropriate environmental variable(s).")

def add_test(name, debug=False, optimize=False, mpi=False, openmp=False,\
             hdf5=False, petsc=False):
    tests.update({name:Test(debug, optimize, mpi, openmp, hdf5, petsc)})

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
add_test('omp-phdf5-petsc-debug', openmp=True, mpi=True, hdf5=True, petsc=True,
                                  debug=True)
add_test('omp-phdf5-petsc-optimize', openmp=True, mpi=True, hdf5=True, petsc=True,
                                     optimize=True)

# Check to see if we are to just print build configuratinos
if opts.print_build_configs:
    for key in tests:
        print('Configuration Name: {0}'.format(key))
        print('  Debug Flags:..........{0}'.format(tests[key].debug))
        print('  Optimization Flags:...{0}'.format(tests[key].optimize))
        print('  HDF5 Active:..........{0}'.format(tests[key].hdf5))
        print('  MPI Active:...........{0}'.format(tests[key].mpi))
        print('  OpenMP Active:........{0}'.format(tests[key].openmp))
        print('  PETSc Active:.........{0}\n'.format(tests[key].petsc))
    exit()

# Delete items of dictionary that don't match regular expression
if opts.build_config is not None:
    for key in tests:
        if not re.search(opts.build_config, key):
            del tests[key]

# Begin testing
shutil.rmtree('build', ignore_errors=True)
for test in tests:
    print('-'*(len(test) + 6))
    print(test + ' tests')
    print('-'*(len(test) + 6))
    sys.stdout.flush()

    # Verify fortran compiler exists
    tests[test].check_compiler()

    # Run CMAKE to configure build
    tests[test].run_cmake()

    # Go into build directory
    os.chdir('build')

    # Build OpenMC
    tests[test].run_make()

    # Run tests
    tests[test].run_ctests()

    # Leave build directory
    os.chdir('..')

    # Copy test log file if failed
    if tests[test].msg == 'Failed on testing.':
        shutil.copy('build/Testing/Temporary/LastTest.log',
                    'LastTest_{0}.log'.format(test))

    # Clean up build
    shutil.rmtree('build', ignore_errors=True)

# Print out summary of results
print('\n' + '='*54)
print('Summary of Compilation Option Testing:\n')

OK = '\033[92m'
FAIL = '\033[91m'
ENDC = '\033[0m'
BOLD = '\033[1m'

for test in tests:
    print(test + '.'*(50 - len(test)), end='')
    if tests[test].success:
        print(BOLD + OK + '[OK]' + ENDC)
    else:
        print(BOLD + FAIL + '[FAILED]' + ENDC)
        print(' '*len(test)+tests[test].msg)
