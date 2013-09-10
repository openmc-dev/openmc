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

if os.environ.has_key('COMPILER'):
    compiler = 'COMPILER=' + os.environ['COMPILER']
else:
    compiler = 'COMPILER=gnu'

def setup():
    # Change to source directory
    os.chdir(pwd + '/../../src')

def test_normal():
    returncode = run(['make','distclean'])
    returncode = run(['make', compiler])
    assert returncode == 0
    shutil.move('openmc', 'openmc-normal')

def test_debug():
    returncode = run(['make','distclean'])
    returncode = run(['make', compiler, 'DEBUG=yes'])
    assert returncode == 0
    shutil.move('openmc', 'openmc-debug')

def test_profile():
    returncode = run(['make','distclean'])
    returncode = run(['make', compiler, 'PROFILE=yes'])
    assert returncode == 0

def test_optimize():
    returncode = run(['make','distclean'])
    returncode = run(['make', compiler, 'OPTIMIZE=yes'])
    assert returncode == 0
    shutil.move('openmc', 'openmc-optimize')

def test_mpi():
    returncode = run(['make','distclean'])
    returncode = run(['make', compiler, 'MPI=yes'])
    assert returncode == 0
    shutil.move('openmc', 'openmc-mpi')

def test_mpi_debug():
    returncode = run(['make','distclean'])
    returncode = run(['make', compiler, 'MPI=yes', 'DEBUG=yes'])
    assert returncode == 0
    shutil.move('openmc', 'openmc-mpi-debug')

def test_mpi_optimize():
    returncode = run(['make','distclean'])
    returncode = run(['make', compiler, 'MPI=yes', 'OPTIMIZE=yes'])
    assert returncode == 0
    shutil.move('openmc', 'openmc-mpi-optimize')

def test_hdf5():
    returncode = run(['make','distclean'])
    returncode = run(['make', compiler, 'HDF5=yes'])
    assert returncode == 0
    shutil.move('openmc', 'openmc-hdf5')

def test_hdf5_debug():
    returncode = run(['make','distclean'])
    returncode = run(['make', compiler, 'HDF5=yes', 'DEBUG=yes'])
    assert returncode == 0
    shutil.move('openmc', 'openmc-hdf5-debug')

def test_hdf5_optimize():
    returncode = run(['make','distclean'])
    returncode = run(['make', compiler, 'HDF5=yes', 'OPTIMIZE=yes'])
    assert returncode == 0
    shutil.move('openmc', 'openmc-hdf5-optimize')

def test_petsc():
    returncode = run(['make','distclean'])
    returncode = run(['make', compiler, 'MPI=yes', 'PETSC=yes'])
    assert returncode == 0
    shutil.move('openmc', 'openmc-petsc')

def test_petsc_debug():
    returncode = run(['make','distclean'])
    returncode = run(['make', compiler, 'MPI=yes', 'PETSC=yes', 'DEBUG=yes'])
    assert returncode == 0
    shutil.move('openmc', 'openmc-petsc-debug')

def test_petsc_optimize():
    returncode = run(['make','distclean'])
    returncode = run(['make', compiler, 'MPI=yes', 'PETSC=yes', 'OPTIMIZE=yes'])
    assert returncode == 0
    shutil.move('openmc', 'openmc-petsc-optimize')

def test_mpi_hdf5():
    returncode = run(['make','distclean'])
    returncode = run(['make', compiler, 'MPI=yes', 'HDF5=yes'])
    assert returncode == 0
    shutil.move('openmc', 'openmc-phdf5')

def test_mpi_hdf5_debug():
    returncode = run(['make','distclean'])
    returncode = run(['make', compiler, 'MPI=yes', 'HDF5=yes', 'DEBUG=yes'])
    assert returncode == 0
    shutil.move('openmc', 'openmc-phdf5-debug')

def test_mpi_hdf5_optimize():
    returncode = run(['make','distclean'])
    returncode = run(['make', compiler, 'MPI=yes', 'HDF5=yes', 'OPTIMIZE=yes'])
    assert returncode == 0
    shutil.move('openmc', 'openmc-phdf5-optimize')

def test_mpi_hdf5_petsc():
    returncode = run(['make','distclean'])
    returncode = run(['make', compiler, 'MPI=yes', 'HDF5=yes', 'PETSC=yes'])
    assert returncode == 0
    shutil.move('openmc', 'openmc-phdf5-petsc')

def test_mpi_hdf5_petsc_debug():
    returncode = run(['make','distclean'])
    returncode = run(['make', compiler, 'MPI=yes', 'HDF5=yes', 'PETSC=yes', 'DEBUG=yes'])
    assert returncode == 0
    shutil.move('openmc', 'openmc-phdf5-petsc-debug')

def test_mpi_hdf5_petsc_optimize():
    returncode = run(['make','distclean'])
    returncode = run(['make', compiler, 'MPI=yes', 'HDF5=yes', 'PETSC=yes', 'OPTIMIZE=yes'])
    assert returncode == 0
    shutil.move('openmc', 'openmc-phdf5-petsc-optimize')

def run(commands):
    proc = Popen(commands, stderr=STDOUT, stdout=PIPE)
    returncode = proc.wait()
    print(proc.communicate()[0])
    return returncode

def teardown(commands):
    returncode = run(['make','distclean'])
    shutil.copy('openmc-normal','openmc')
