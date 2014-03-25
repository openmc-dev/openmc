#!/usr/bin/env python

"""Compilation tests

This set of tests makes sure that OpenMC can compile for many different
combinations of compilation flags and options. By default, only the gfortran
and ifort compilers are tested. This requires that the MPI, HDF5, and PETSC
libraries are already set up appropriately in the Makefile.

"""

import os
from subprocess import Popen, STDOUT, PIPE
import shutil

pwd = os.path.dirname(__file__)

# Set default copmilers
fc = 'gfortran'
h5fc = 'h5fc'
h5pfc = 'h5pfc'
mpifc = 'mpif90'


def setup():
    # Change to source directory
    os.chdir(pwd + '/../../src')


def test_normal():
    shutil.rmtree('build', ignore_errors=True)
    returncode = run('FC={0} cmake -H. -Bbuild'.format(fc))
    assert returncode == 0
    returncode = run('make -s -C build')
    assert returncode == 0
    shutil.move('build/bin/openmc', 'openmc-normal')


def test_debug():
    shutil.rmtree('build', ignore_errors=True)
    returncode = run('FC={0} cmake -Ddebug=on -H. -Bbuild'.format(fc))
    assert returncode == 0
    returncode = run('make -s -C build')
    assert returncode == 0
    shutil.move('build/bin/openmc', 'openmc-debug')


def test_profile():
    shutil.rmtree('build', ignore_errors=True)
    returncode = run('FC={0} cmake -Dprofile=on -H. -Bbuild'.format(fc))
    assert returncode == 0
    returncode = run('make -s -C build')
    assert returncode == 0
    shutil.move('build/bin/openmc', 'openmc-profile')


def test_optimize():
    shutil.rmtree('build', ignore_errors=True)
    returncode = run('FC={0} cmake -Doptimize=on -H. -Bbuild'.format(fc))
    assert returncode == 0
    returncode = run('make -s -C build')
    assert returncode == 0
    shutil.move('build/bin/openmc', 'openmc-optimize')


def test_mpi():
    shutil.rmtree('build', ignore_errors=True)
    returncode = run('FC={0} cmake -H. -Bbuild'.format(mpifc))
    assert returncode == 0
    returncode = run('make -s -C build')
    assert returncode == 0
    shutil.move('build/bin/openmc', 'openmc-mpi')


def test_mpi_debug():
    shutil.rmtree('build', ignore_errors=True)
    returncode = run('FC={0} cmake -Ddebug=on -H. -Bbuild'.format(mpifc))
    assert returncode == 0
    returncode = run('make -s -C build')
    assert returncode == 0
    shutil.move('build/bin/openmc', 'openmc-mpi-debug')


def test_mpi_optimize():
    shutil.rmtree('build', ignore_errors=True)
    returncode = run('FC={0} cmake -Doptimize=on -H. -Bbuild'.format(mpifc))
    assert returncode == 0
    returncode = run('make -s -C build')
    assert returncode == 0
    shutil.move('build/bin/openmc', 'openmc-mpi-optimize')


def test_omp():
    shutil.rmtree('build', ignore_errors=True)
    returncode = run('FC={0} cmake -Dopenmp=on -H. -Bbuild'.format(fc))
    assert returncode == 0
    returncode = run('make -s -C build')
    assert returncode == 0
    shutil.move('build/bin/openmc', 'openmc-omp')


def test_omp_debug():
    shutil.rmtree('build', ignore_errors=True)
    returncode = run('FC={0} cmake -Dopenmp=on -Ddebug=on '.format(fc) +
                     '-H. -Bbuild')
    assert returncode == 0
    returncode = run('make -s -C build')
    assert returncode == 0
    shutil.move('build/bin/openmc', 'openmc-omp-debug')


def test_omp_optimize():
    shutil.rmtree('build', ignore_errors=True)
    returncode = run('FC={0} cmake -Dopenmp=on -Doptimize=on '.format(fc) +
                     '-H. -Bbuild')
    assert returncode == 0
    returncode = run('make -s -C build')
    assert returncode == 0
    shutil.move('build/bin/openmc', 'openmc-omp-optimize')


def test_hdf5():
    shutil.rmtree('build', ignore_errors=True)
    returncode = run('FC={0} cmake -H. -Bbuild'.format(h5fc))
    assert returncode == 0
    returncode = run('make -s -C build')
    assert returncode == 0
    shutil.move('build/bin/openmc', 'openmc-hdf5')


def test_hdf5_debug():
    shutil.rmtree('build', ignore_errors=True)
    returncode = run('FC={0} cmake -Ddebug=on -H. -Bbuild'.format(h5fc))
    assert returncode == 0
    returncode = run('make -s -C build')
    assert returncode == 0
    shutil.move('build/bin/openmc', 'openmc-hdf5-debug')


def test_hdf5_optimize():
    shutil.rmtree('build', ignore_errors=True)
    returncode = run('FC={0} cmake -Doptimize=on -H. -Bbuild'.format(h5fc))
    assert returncode == 0
    returncode = run('make -s -C build')
    assert returncode == 0
    shutil.move('build/bin/openmc', 'openmc-hdf5-optimize')


def test_omp_hdf5():
    shutil.rmtree('build', ignore_errors=True)
    returncode = run('FC={0} cmake -Dopenmp=on -H. -Bbuild'.format(h5fc))
    assert returncode == 0
    returncode = run('make -s -C build')
    assert returncode == 0
    shutil.move('build/bin/openmc', 'openmc-omp-hdf5')


def test_omp_hdf5_debug():
    shutil.rmtree('build', ignore_errors=True)
    returncode = run('FC={0} cmake -Dopenmp=on -Ddebug=on '.format(h5fc) +
                     '-H. -Bbuild')
    assert returncode == 0
    returncode = run('make -s -C build')
    assert returncode == 0
    shutil.move('build/bin/openmc', 'openmc-omp-hdf5-debug')


def test_omp_hdf5_optimize():
    shutil.rmtree('build', ignore_errors=True)
    returncode = run('FC={0} cmake -Dopenmp=on -Doptimize=on '.format(h5fc) +
                     '-H. -Bbuild')
    assert returncode == 0
    returncode = run('make -s -C build')
    assert returncode == 0
    shutil.move('build/bin/openmc', 'openmc-omp-hdf5-optimize')


def test_petsc():
    shutil.rmtree('build', ignore_errors=True)
    returncode = run('FC={0} cmake -Dpetsc=on -H. -Bbuild'.format(mpifc))
    assert returncode == 0
    returncode = run('make -s -C build')
    assert returncode == 0
    shutil.move('build/bin/openmc', 'openmc-petsc')


def test_petsc_debug():
    shutil.rmtree('build', ignore_errors=True)
    returncode = run('FC={0} cmake -Dpetsc=on -Ddebug=on '.format(mpifc) +
                     '-H. -Bbuild')
    assert returncode == 0
    returncode = run('make -s -C build')
    assert returncode == 0
    shutil.move('build/bin/openmc', 'openmc-petsc-debug')


def test_petsc_optimize():
    shutil.rmtree('build', ignore_errors=True)
    returncode = run('FC={0} cmake -Dpetsc=on -Doptimize=on '.format(mpifc) +
                     '-H. -Bbuild')
    assert returncode == 0
    returncode = run('make -s -C build')
    assert returncode == 0
    shutil.move('build/bin/openmc', 'openmc-petsc-optimize')


def test_mpi_omp():
    shutil.rmtree('build', ignore_errors=True)
    returncode = run('FC={0} cmake -Dopenmp=on -H. -Bbuild'.format(mpifc))
    assert returncode == 0
    returncode = run('make -s -C build')
    assert returncode == 0
    shutil.move('build/bin/openmc', 'openmc-mpi-omp')


def test_mpi_omp_debug():
    shutil.rmtree('build', ignore_errors=True)
    returncode = run('FC={0} cmake -Dopenmp=on -Ddebug=on '.format(mpifc) +
                     '-H. -Bbuild')
    assert returncode == 0
    returncode = run('make -s -C build')
    assert returncode == 0
    shutil.move('build/bin/openmc', 'openmc-mpi-omp-debug')


def test_mpi_omp_optimize():
    shutil.rmtree('build', ignore_errors=True)
    returncode = run('FC={0} cmake -Dopenmp=on -Doptimize=on '.format(mpifc) +
                     '-H. -Bbuild')
    assert returncode == 0
    returncode = run('make -s -C build')
    assert returncode == 0
    shutil.move('build/bin/openmc', 'openmc-mpi-omp-optimize')


def test_mpi_hdf5():
    shutil.rmtree('build', ignore_errors=True)
    returncode = run('FC={0} cmake -H. -Bbuild'.format(h5pfc))
    assert returncode == 0
    returncode = run('make -s -C build')
    assert returncode == 0
    shutil.move('build/bin/openmc', 'openmc-phdf5')


def test_mpi_hdf5_debug():
    shutil.rmtree('build', ignore_errors=True)
    returncode = run('FC={0} cmake -Ddebug=on -H. -Bbuild'.format(h5pfc))
    assert returncode == 0
    returncode = run('make -s -C build')
    assert returncode == 0
    shutil.move('build/bin/openmc', 'openmc-phdf5-debug')


def test_mpi_hdf5_optimize():
    shutil.rmtree('build', ignore_errors=True)
    returncode = run('FC={0} cmake -Doptimize=on -H. -Bbuild'.format(h5pfc))
    assert returncode == 0
    returncode = run('make -s -C build')
    assert returncode == 0
    shutil.move('build/bin/openmc', 'openmc-phdf5-optimize')


def test_mpi_omp_hdf5():
    shutil.rmtree('build', ignore_errors=True)
    returncode = run('FC={0} cmake -Dopenmp=on -H. -Bbuild'.format(h5pfc))
    assert returncode == 0
    returncode = run('make -s -C build')
    assert returncode == 0
    shutil.move('build/bin/openmc', 'openmc-phdf5-omp')


def test_mpi_omp_hdf5_debug():
    shutil.rmtree('build', ignore_errors=True)
    returncode = run('FC={0} cmake -Dopenmp=on -Ddebug=on '.format(h5pfc) +
                     '-H. -Bbuild')
    assert returncode == 0
    returncode = run('make -s -C build')
    assert returncode == 0
    shutil.move('build/bin/openmc', 'openmc-phdf5-omp-debug')


def test_mpi_omp_hdf5_optimize():
    shutil.rmtree('build', ignore_errors=True)
    returncode = run('FC={0} cmake -Dopenmp=on -Doptimize=on '.format(h5pfc) +
                     '-H. -Bbuild')
    assert returncode == 0
    returncode = run('make -s -C build')
    assert returncode == 0
    shutil.move('build/bin/openmc', 'openmc-phdf5-omp-optimize')


def test_mpi_hdf5_petsc():
    shutil.rmtree('build', ignore_errors=True)
    returncode = run('FC={0} cmake -Dpetsc=on -H. -Bbuild'.format(h5pfc))
    assert returncode == 0
    returncode = run('make -s -C build')
    assert returncode == 0
    shutil.move('build/bin/openmc', 'openmc-phdf5-petsc')


def test_mpi_hdf5_petsc_debug():
    shutil.rmtree('build', ignore_errors=True)
    returncode = run('FC={0} cmake -Dpetsc=on -Ddebug=on '.format(h5pfc) +
                     '-H. -Bbuild')
    assert returncode == 0
    returncode = run('make -s -C build')
    assert returncode == 0
    shutil.move('build/bin/openmc', 'openmc-phdf5-petsc-debug')


def test_mpi_hdf5_petsc_optimize():
    shutil.rmtree('build', ignore_errors=True)
    returncode = run('FC={0} cmake -Dpetsc=on -Doptimize=on '.format(h5pfc) +
                     '-H. -Bbuild')
    assert returncode == 0
    returncode = run('make -s -C build')
    assert returncode == 0
    shutil.move('build/bin/openmc', 'openmc-phdf5-petsc-optimize')


def test_mpi_omp_hdf5_petsc():
    shutil.rmtree('build', ignore_errors=True)
    returncode = run('FC={0} cmake -Dpetsc=on -Dopenmp=on '.format(h5pfc) +
                     '-H. -Bbuild')
    assert returncode == 0
    returncode = run('make -s -C build')
    assert returncode == 0
    shutil.move('build/bin/openmc', 'openmc-omp-phdf5-petsc')


def test_mpi_omp_hdf5_petsc_debug():
    shutil.rmtree('build', ignore_errors=True)
    returncode = run('FC={0} cmake -Dpetsc=on -Dopenmp=on '.format(h5pfc) +
                     '-Ddebug=on -H. -Bbuild')
    assert returncode == 0
    returncode = run('make -s -C build')
    assert returncode == 0
    shutil.move('build/bin/openmc', 'openmc-omp-phdf5-petsc-debug')


def test_mpi_omp_hdf5_petsc_optimize():
    shutil.rmtree('build', ignore_errors=True)
    returncode = run('FC={0} cmake -Dpetsc=on -Dopenmp=on '.format(h5pfc) +
                     '-Doptimize=on -H. -Bbuild')
    assert returncode == 0
    returncode = run('make -s -C build')
    assert returncode == 0
    shutil.move('build/bin/openmc', 'openmc-omp-phdf5-petsc-optimize')


def run(commands):
    proc = Popen(commands, shell=True, stderr=STDOUT, stdout=PIPE)
    print(proc.communicate()[0])
    returncode = proc.returncode
    return returncode


def teardown(commands):
    shutil.rmtree('build', ignore_errors=True)
    shutil.copy('openmc-normal', 'openmc')
