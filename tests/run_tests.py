#!/usr/bin/env python

from __future__ import print_function

import os
import sys
from subprocess import call 
from collections import OrderedDict

# Compiler paths
FC_DEFAULT='gfortran'
MPI_DIR='/opt/mpich/3.0.4-gnu'
HDF5_DIR='/opt/hdf5/1.8.12-gnu'
PHDF5_DIR='/opt/phdf5/1.8.12-gnu'
PETSC_DIR='/opt/petsc/3.4.3-gnu'

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
        self.setup_cmake()

    def setup_cmake(self):
        # Default cmake
        self.cmake = ['cmake','-H.','-Bbuild']

        # Check for MPI/HDF5
        if self.mpi and not self.hdf5:
            self.fc = MPI_DIR+'/bin/mpif90'
        elif not self.mpi and self.hdf5:
            self.fc = HDF5_DIR+'/bin/h5fc'
        elif self.mpi and self.hdf5:
            self.fc = PHDF5_DIR+'/bin/h5pfc'
        else:
            self.fc = FC_DEFAULT

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

    def run_cmake(self):
        os.environ['FC'] = self.fc
        rc = call(self.cmake)
        if rc != 0:
            self.success = False

    def run_make(self):
        if not self.success:
            return
        rc = call(['make','-j', '-s','-C','build'])
        if rc != 0:
            self.success = False

    def run_ctests(self):
        if not self.success:
            return
        rc = call(['make','test','-C','build'])
        if rc != 0:
            self.success = False

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
add_test('omp-phdf5-petsc-debug', openmp=True, mpi=True, hdf5=True, petsc=True, debug=True)
add_test('omp-phdf5-petsc-optimize', openmp=True, mpi=True, hdf5=True, petsc=True, optimize=True)

# Change to source directory and clean it up
os.chdir('../src')
call(['make','distclean'])

# Process command line arguments
if len(sys.argv) > 1:
    flags = [i for i in sys.argv[1:] if i.startswith('-')]
    tests_ = [i for i in sys.argv[1:] if not i.startswith('-')]

    # Check for special subsets of tests
    tests__ = []
    for i in tests_:

        # This checks for any subsets of tests. The string after
        # all-XXXX will be used to search through all tests.
        if i.startswith('all-'):
            suffix = i.split('all-')[1]
            for j in tests:
                if j.rfind(suffix) != -1:
                    try:
                        tests__.index(j)
                    except ValueError:
                        tests__.append(j)

        # Test name specified on command line
        else:
            try:
                tests__.index(i)
            except ValueError:
                tests__.append(i)  # append specific test (e.g., mpi-debug)

    # Delete items of dictionary that are not in tests__
    for key in iter(tests):
        try:
            tests__.index(key)
        except ValueError:
            del tests[key]

# Begin testing
for test in iter(tests):
    print('-'*(len(test) + 6))
    print(test + ' tests')
    print('-'*(len(test) + 6))


    # Run CMAKE to configure build
    tests[test].run_cmake()

    # Build OpenMC
    tests[test].run_make()

    # Run tests
    tests[test].run_ctests()

    # Clean up build
    call(['make','distclean'])

# Print out summary of results
print('\n' + '='*54)
print('Summary of Compilation Option Testing:\n')

OK = '\033[92m'
FAIL = '\033[91m'
ENDC = '\033[0m'
BOLD = '\033[1m'

for test in iter(tests):
    print(test + '.'*(50 - len(test)), end='')
    if tests[test].success:
        print(BOLD + OK + '[OK]' + ENDC)
    else:
        print(BOLD + FAIL + '[FAILED]' + ENDC)
