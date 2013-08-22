#!/bin/sh env python

# load in packages
import os
import sys
import nose
import glob
from subprocess import call

OK = '\033[92m'
FAIL = '\033[91m'
ENDC = '\033[0m'
BOLD = '\033[1m'

pwd = os.getcwd()

# clean up all previous executables
openmc_exe = glob.glob(pwd + '/../src/openmc*')
for exe in openmc_exe:
    os.remove(exe)

# run compile test
result_compile = nose.run(argv=['run_tests.py','-v','test_compile'])
if not result_compile:
    sys.stdout.write('\nDid not pass compile tests.')

# run gfortran tests
sys.stdout.write('\n--------------\n')
sys.stdout.write('gfortran tests')
sys.stdout.write('\n--------------\n')
os.chdir(pwd)
os.rename(pwd + '/../src/openmc-gfortran', pwd + '/../src/openmc')
result_gfortran = nose.run(argv=['run_tests.py','-v','--exclude','test_compile'])
if not result_gfortran:
    sys.stdout.write('\nDid not pass gfortran tests\n')

# run gfortran-dbg tests
sys.stdout.write('\n--------------\n')
sys.stdout.write('gfortran-dbg tests')
sys.stdout.write('\n--------------\n')
os.chdir(pwd)
os.rename(pwd + '/../src/openmc-gfortran-dbg', pwd + '/../src/openmc')
result_gfortran_dbg = nose.run(argv=['run_tests.py','-v','--exclude','test_compile'])
if not result_gfortran_dbg:
    sys.stdout.write('\nDid not pass gfortran-dbg tests\n')

# run gfortran-opt tests
sys.stdout.write('\n--------------\n')
sys.stdout.write('gfortran-opt tests')
sys.stdout.write('\n--------------\n')
os.chdir(pwd)
os.rename(pwd + '/../src/openmc-gfortran-opt', pwd + '/../src/openmc')
result_gfortran_opt = nose.run(argv=['run_tests.py','-v','--exclude','test_compile'])
if not result_gfortran_opt:
    sys.stdout.write('\nDid not pass gfortran-opt tests\n')

# run gfortran-hdf5 tests
sys.stdout.write('\n--------------\n')
sys.stdout.write('gfortran-hdf5 tests')
sys.stdout.write('\n--------------\n')
os.chdir(pwd)
os.rename(pwd + '/../src/openmc-gfortran-hdf5', pwd + '/../src/openmc')
result_gfortran_hdf5 = nose.run(argv=['run_tests.py','-v','--exclude','test_compile'])
if not result_gfortran_hdf5:
    sys.stdout.write('\nDid not pass gfortran-hdf5 tests\n')

# run gfortran-mpi tests
sys.stdout.write('\n--------------\n')
sys.stdout.write('gfortran-mpi tests')
sys.stdout.write('\n--------------\n')
os.chdir(pwd)
os.rename(pwd + '/../src/openmc-gfortran-mpi', pwd + '/../src/openmc')
result_gfortran_mpi = nose.run(argv=['run_tests.py','-v','--exclude','test_compile',
                        '--mpi-np','3','--mpi-exec',
                        '/opt/mpich/3.0.4-gnu/bin/mpiexec'])
if not result_gfortran_mpi:
    sys.stdout.write('\nDid not pass gfortran-mpi tests\n')

# run gfortran-phdf5 tests
sys.stdout.write('\n--------------\n')
sys.stdout.write('gfortran-phdf5 tests')
sys.stdout.write('\n--------------\n')
os.chdir(pwd)
os.rename(pwd + '/../src/openmc-gfortran-phdf5', pwd + '/../src/openmc')
result_gfortran_phdf5 = nose.run(argv=['run_tests.py','-v','--exclude','test_compile',
                        '--mpi-np','3','--mpi-exec',
                        '/opt/mpich/3.0.4-gnu/bin/mpiexec'])
if not result_gfortran_phdf5:
    sys.stdout.write('\nDid not pass gfortran-phdf5 tests\n')

# run gfortran-petsc tests
sys.stdout.write('\n--------------\n')
sys.stdout.write('gfortran-petsc tests')
sys.stdout.write('\n--------------\n')
os.chdir(pwd)
os.rename(pwd + '/../src/openmc-gfortran-petsc', pwd + '/../src/openmc')
result_gfortran_petsc = nose.run(argv=['run_tests.py','-v','--exclude','test_compile',
                        '--mpi-np','3','--mpi-exec',
                        '/opt/mpich/3.0.4-gnu/bin/mpiexec'])
if not result_gfortran_petsc:
    sys.stdout.write('\nDid not pass gfortran-petsc tests\n')

# run gfortran-phdf5-petsc tests
sys.stdout.write('\n--------------\n')
sys.stdout.write('gfortran-phdf5-petsc tests')
sys.stdout.write('\n--------------\n')
os.chdir(pwd)
os.rename(pwd + '/../src/openmc-gfortran-phdf5-petsc', pwd + '/../src/openmc')
result_gfortran_phdf5_petsc = nose.run(argv=['run_tests.py','-v','--exclude','test_compile',
                        '--mpi-np','3','--mpi-exec',
                        '/opt/mpich/3.0.4-gnu/bin/mpiexec'])
if not result_gfortran_phdf5_petsc:
    sys.stdout.write('\nDid not pass gfortran-phdf5-petsc tests\n')

# run gfortran-phdf5-petsc-opt tests
sys.stdout.write('\n--------------\n')
sys.stdout.write('gfortran-phdf5-petsc-opt tests')
sys.stdout.write('\n--------------\n')
os.chdir(pwd)
os.rename(pwd + '/../src/openmc-gfortran-phdf5-petsc-opt', pwd + '/../src/openmc')
result_gfortran_phdf5_petsc_opt = nose.run(argv=['run_tests.py','-v','--exclude','test_compile',
                        '--mpi-np','3','--mpi-exec',
                        '/opt/mpich/3.0.4-gnu/bin/mpiexec'])
if not result_gfortran_phdf5_petsc_opt:
    sys.stdout.write('\nDid not pass gfortran-phdf5-petsc-opt tests\n')

sys.stdout.write('\n======================================================\n')
sys.stdout.write('\nSummary of Compilation Option Testing:\n')

sys.stdout.write('\ngfortran tests...........................')
if result_gfortran:
    sys.stdout.write(BOLD + OK + '[OK]' + ENDC)
else:
    sys.stdout.write(BOLD + FAIL + '[FAILED]' + ENDC)

sys.stdout.write('\ngfortran-DEBUG tests.....................')
if result_gfortran_dbg:
    sys.stdout.write(BOLD + OK + '[OK]' + ENDC)
else:
    sys.stdout.write(BOLD + FAIL + '[FAILED]' + ENDC)

sys.stdout.write('\ngfortran-OPTIMIZE tests..................')
if result_gfortran_opt:
    sys.stdout.write(BOLD + OK + '[OK]' + ENDC)
else:
    sys.stdout.write(BOLD + FAIL + '[FAILED]' + ENDC)

sys.stdout.write('\ngfortran-HDF5 tests......................')
if result_gfortran_hdf5:
    sys.stdout.write(BOLD + OK + '[OK]' + ENDC)
else:
    sys.stdout.write(BOLD + FAIL + '[FAILED]' + ENDC)

sys.stdout.write('\ngfortran-MPI tests.......................')
if result_gfortran_mpi:
    sys.stdout.write(BOLD + OK + '[OK]' + ENDC)
else:
    sys.stdout.write(BOLD + FAIL + '[FAILED]' + ENDC)

sys.stdout.write('\ngfortran-MPI-HDF5 tests..................')
if result_gfortran_phdf5:
    sys.stdout.write(BOLD + OK + '[OK]' + ENDC)
else:
    sys.stdout.write(BOLD + FAIL + '[FAILED]' + ENDC)

sys.stdout.write('\ngfortran-MPI-PETSC tests.................')
if result_gfortran_petsc:
    sys.stdout.write(BOLD + OK + '[OK]' + ENDC)
else:
    sys.stdout.write(BOLD + FAIL + '[FAILED]' + ENDC)

sys.stdout.write('\ngfortran-MPI-HDF5-PETSC tests............')
if result_gfortran_phdf5_petsc:
    sys.stdout.write(BOLD + OK + '[OK]' + ENDC)
else:
    sys.stdout.write(BOLD + FAIL + '[FAILED]' + ENDC)

sys.stdout.write('\ngfortran-MPI-PETSC-HDF5-OPTIMIZE tests...')
if result_gfortran_phdf5_petsc_opt:
    sys.stdout.write(BOLD + OK + '[OK]' + ENDC)
else:
    sys.stdout.write(BOLD + FAIL + '[FAILED]' + ENDC + '\n')
