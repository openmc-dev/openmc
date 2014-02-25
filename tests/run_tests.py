#!/usr/bin/env python

from __future__ import print_function

import os
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
            self.sucess = False

def add_test(name, debug=False, optimize=False, mpi=False, openmp=False,\
             hdf5=False, petsc=False):
    tests.update({name:Test(debug, optimize, mpi, openmp, hdf5, petsc)})

# List of tests
add_test('normal')
add_test('debug', debug=True)
add_test('optimize', optimize=True)
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

# Call tests
os.chdir('../src')
call(['make','distclean'])
for test in iter(tests):
    tests[test].run_cmake()
    call(['make','distclean'])
