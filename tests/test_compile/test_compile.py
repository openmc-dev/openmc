#!/usr/bin/env python

"""Compilation tests

This set of tests makes sure that OpenMC can compile for many different
combinations of compilation flags and options. By default, only the gfortran and
ifort compilers are tested. This requires that the MPI, HDF5, and PETSC
libraries are already set up appropriately in the Makefile.

"""

import os
from subprocess import Popen, STDOUT, PIPE
import shutil

pwd = os.path.dirname(__file__)

def setup():
    # Change to source directory
    os.chdir(pwd + '/../../src')

def test_gfortran():
    returncode = run(['make','distclean'])
    returncode = run(['make'])
    assert returncode == 0
    shutil.move('openmc', 'openmc-gfortran')

def test_gfortran_debug():
    returncode = run(['make','distclean'])
    returncode = run(['make', 'DEBUG=yes'])
    assert returncode == 0
    shutil.move('openmc', 'openmc-gfortran-dbg')

def test_gfortran_profile():
    returncode = run(['make','distclean'])
    returncode = run(['make', 'PROFILE=yes'])
    assert returncode == 0

def test_gfortran_optimize():
    returncode = run(['make','distclean'])
    returncode = run(['make', 'OPTIMIZE=yes'])
    assert returncode == 0
    shutil.move('openmc', 'openmc-gfortran-opt')

def test_gfortran_mpi():
    returncode = run(['make','distclean'])
    returncode = run(['make', 'MPI=yes'])
    assert returncode == 0
    shutil.move('openmc', 'openmc-gfortran-mpi')

def test_gfortran_hdf5():
    returncode = run(['make','distclean'])
    returncode = run(['make', 'HDF5=yes'])
    assert returncode == 0
    shutil.move('openmc', 'openmc-gfortran-hdf5')

def test_gfortran_petsc():
    returncode = run(['make','distclean'])
    returncode = run(['make', 'MPI=yes', 'PETSC=yes'])
    assert returncode == 0
    shutil.move('openmc', 'openmc-gfortran-petsc')

def test_gfortran_mpi_hdf5():
    returncode = run(['make','distclean'])
    returncode = run(['make', 'MPI=yes', 'HDF5=yes'])
    assert returncode == 0
    shutil.move('openmc', 'openmc-gfortran-phdf5')

def test_gfortran_mpi_hdf5_petsc():
    returncode = run(['make','distclean'])
    returncode = run(['make', 'MPI=yes', 'HDF5=yes', 'PETSC=yes'])
    assert returncode == 0
    shutil.move('openmc', 'openmc-gfortran-phdf5-petsc')

def test_gfortran_mpi_hdf5_petsc_optimize():
    returncode = run(['make','distclean'])
    returncode = run(['make', 'MPI=yes', 'HDF5=yes', 'PETSC=yes', 'OPTIMIZE=yes'])
    assert returncode == 0
    shutil.move('openmc', 'openmc-gfortran-phdf5-petsc-opt')

def run(commands):
    proc = Popen(commands, stderr=STDOUT, stdout=PIPE)
    returncode = proc.wait()
    print(proc.communicate()[0])
    return returncode

def teardown(commands):
    returncode = run(['make','distclean'])
    shutil.copy('openmc-gfortran','openmc')
