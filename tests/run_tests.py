#!/usr/bin/env python

# load in packages
import os
import sys
import nose
import glob
from subprocess import call
from nose_mpi import NoseMPI

# set mpiexec path
mpiexec = '/opt/mpich/3.0.4-gnu/bin/mpiexec'

# get current working directory
pwd = os.getcwd()
sys.path.append(pwd)

# setup color printing codes
OK = '\033[92m'
FAIL = '\033[91m'
ENDC = '\033[0m'
BOLD = '\033[1m'

# initialize run logicals
run_compile = True
run_gfortran = True
run_gfortran_dbg = True
run_gfortran_opt = True
run_gfortran_hdf5 = True
run_gfortran_mpi = True
run_gfortran_phdf5 = True
run_gfortran_petsc = True
run_gfortran_phdf5_petsc = True
run_gfortran_phdf5_petsc_opt = True

# check for extra command line options
opts = ['all']  # default is run everything
if len(sys.argv) > 1:
    opts = sys.argv[1:]

# if not running everything check for compile tests
# if compile is only present do not run other test cases
if opts[0] != 'all':
    try:
        idx = opts.index('compile')
    except:
        run_compile = False

# check for other test cases to be run
if opts[0] != 'all':
    run_gfortran = False
    run_gfortran_dbg = False
    run_gfortran_opt = False
    run_gfortran_hdf5 = False
    run_gfortran_mpi = False
    run_gfortran_phdf5 = False
    run_gfortran_petsc = False
    run_gfortran_phdf5_petsc = False
    run_gfortran_phdf5_petsc_opt = False
    for item in opts:
        if item == 'all_tests':
            run_gfortran = True
            run_gfortran_dbg = True
            run_gfortran_opt = True
            run_gfortran_hdf5 = True
            run_gfortran_mpi = True
            run_gfortran_phdf5 = True
            run_gfortran_petsc = True
            run_gfortran_phdf5_petsc = True
            run_gfortran_phdf5_petsc_opt = True
            break
        elif item == 'gfortran':
            run_gfortran = True
        elif item == 'gfortran-dbg':
            run_gfortran_dbg = True
        elif item == 'gfortran-opt':
            run_gfortran_opt = True
        elif item == 'gfortran-hdf5':
            run_gfortran_hdf5 = True
        elif item == 'gfortran-mpi':
            run_gfortran_mpi = True
        elif item == 'gfortran-phdf5':
            run_gfortran_phdf5 = True
        elif item == 'gfortran-petsc':
            run_gfortrani_petsc = True
        elif item == 'gfortran-phdf5-petsc':
            run_gfortran_phdf5_petsc = True
        elif item == 'gfortran-phdf5-petsc-opt':
            run_gfortran_phdf5_petsc_opt = True

# begin running tests
if run_compile:
    # clean up all previous executables
    openmc_exe = glob.glob(pwd + '/../src/openmc*')
    for exe in openmc_exe:
        os.remove(exe)

    # run compile test
    result_compile = nose.run(argv=['run_tests.py', '-v', 'test_compile'])
    if not result_compile:
        sys.stdout.write('\nDid not pass compile tests.')

if run_gfortran:
    # run gfortran tests
    sys.stdout.write('\n--------------\n')
    sys.stdout.write('gfortran tests')
    sys.stdout.write('\n--------------\n')
    os.chdir(pwd)
    os.rename(pwd + '/../src/openmc-gfortran', pwd + '/../src/openmc')
    result_gfortran = nose.run(argv=[
        'run_tests.py', '-v', '--exclude', 'test_compile',
        '--exclude', 'cmfd'], addplugins=[NoseMPI()])
    os.rename(pwd + '/../src/openmc', pwd + '/../src/openmc-gfortran')
    if not result_gfortran:
        sys.stdout.write('\nDid not pass gfortran tests\n')

if run_gfortran_dbg:
    # run gfortran-dbg tests
    sys.stdout.write('\n--------------\n')
    sys.stdout.write('gfortran-dbg tests')
    sys.stdout.write('\n--------------\n')
    os.chdir(pwd)
    os.rename(pwd + '/../src/openmc-gfortran-dbg', pwd + '/../src/openmc')
    result_gfortran_dbg = nose.run(argv=[
        'run_tests.py', '-v', '--exclude', 'test_compile',
        '--exclude', 'cmfd'], addplugins=[NoseMPI()])
    os.rename(pwd + '/../src/openmc', pwd + '/../src/openmc-gfortran-dbg')
    if not result_gfortran_dbg:
        sys.stdout.write('\nDid not pass gfortran-dbg tests\n')

if run_gfortran_opt:
    # run gfortran-opt tests
    sys.stdout.write('\n--------------\n')
    sys.stdout.write('gfortran-opt tests')
    sys.stdout.write('\n--------------\n')
    os.chdir(pwd)
    os.rename(pwd + '/../src/openmc-gfortran-opt', pwd + '/../src/openmc')
    result_gfortran_opt = nose.run(argv=[
        'run_tests.py', '-v', '--exclude', 'test_compile',
        '--exclude', 'cmfd'], addplugins=[NoseMPI()])
    os.rename(pwd + '/../src/openmc', pwd + '/../src/openmc-gfortran-opt')
    if not result_gfortran_opt:
        sys.stdout.write('\nDid not pass gfortran-opt tests\n')

if run_gfortran_hdf5:
    # run gfortran-hdf5 tests
    sys.stdout.write('\n--------------\n')
    sys.stdout.write('gfortran-hdf5 tests')
    sys.stdout.write('\n--------------\n')
    os.chdir(pwd)
    os.rename(pwd + '/../src/openmc-gfortran-hdf5', pwd + '/../src/openmc')
    result_gfortran_hdf5 = nose.run(argv=[
        'run_tests.py', '-v', '--exclude', 'test_compile',
        '--exclude', 'cmfd'], addplugins=[NoseMPI()])
    os.rename(pwd + '/../src/openmc', pwd + '/../src/openmc-gfortran-hdf5')
    if not result_gfortran_hdf5:
        sys.stdout.write('\nDid not pass gfortran-hdf5 tests\n')

if run_gfortran_mpi:
    # run gfortran-mpi tests
    sys.stdout.write('\n--------------\n')
    sys.stdout.write('gfortran-mpi tests')
    sys.stdout.write('\n--------------\n')
    os.chdir(pwd)
    os.rename(pwd + '/../src/openmc-gfortran-mpi', pwd + '/../src/openmc')
    result_gfortran_mpi = nose.run(argv=[
        'run_tests.py', '-v', '--exclude', 'test_compile',
        '--exclude', 'cmfd', '--mpi-np', '3', '--mpi-exec', mpiexec],
        addplugins=[NoseMPI()])
    os.rename(pwd + '/../src/openmc', pwd + '/../src/openmc-gfortran-mpi')
    if not result_gfortran_mpi:
        sys.stdout.write('\nDid not pass gfortran-mpi tests\n')

if run_gfortran_phdf5:
    # run gfortran-phdf5 tests
    # set mpiexec path
    sys.stdout.write('\n--------------\n')
    sys.stdout.write('gfortran-phdf5 tests')
    sys.stdout.write('\n--------------\n')
    os.chdir(pwd)
    os.rename(pwd + '/../src/openmc-gfortran-phdf5', pwd + '/../src/openmc')
    result_gfortran_phdf5 = nose.run(argv=[
        'run_tests.py', '-v', '--exclude', 'test_compile', '--exclude',
        'cmfd', '--mpi-np', '3', '--mpi-exec', mpiexec],
        addplugins=[NoseMPI()])
    os.rename(pwd + '/../src/openmc', pwd + '/../src/openmc-gfortran-phdf5')
    if not result_gfortran_phdf5:
        sys.stdout.write('\nDid not pass gfortran-phdf5 tests\n')

if run_gfortran_petsc:
    # run gfortran-petsc tests
    sys.stdout.write('\n--------------\n')
    sys.stdout.write('gfortran-petsc tests')
    sys.stdout.write('\n--------------\n')
    os.chdir(pwd)
    os.rename(pwd + '/../src/openmc-gfortran-petsc', pwd + '/../src/openmc')
    result_gfortran_petsc = nose.run(argv=[
        'run_tests.py', '-v', '--exclude', 'test_compile',
        '--mpi-np', '3', '--mpi-exec', mpiexec], addplugins=[NoseMPI()])
    os.rename(pwd + '/../src/openmc', pwd + '/../src/openmc-gfortran-petsc')
    if not result_gfortran_petsc:
        sys.stdout.write('\nDid not pass gfortran-petsc tests\n')

if run_gfortran_phdf5_petsc:
    # run gfortran-phdf5-petsc tests
    sys.stdout.write('\n--------------\n')
    sys.stdout.write('gfortran-phdf5-petsc tests')
    sys.stdout.write('\n--------------\n')
    os.chdir(pwd)
    os.rename(pwd + '/../src/openmc-gfortran-phdf5-petsc',
              pwd + '/../src/openmc')
    result_gfortran_phdf5_petsc = nose.run(argv=[
        'run_tests.py', '-v', '--exclude', 'test_compile',
        '--mpi-np', '3', '--mpi-exec', mpiexec], addplugins=[NoseMPI()])
    os.rename(pwd + '/../src/openmc',
              pwd + '/../src/openmc-gfortran-phdf5-petsc')
    if not result_gfortran_phdf5_petsc:
        sys.stdout.write('\nDid not pass gfortran-phdf5-petsc tests\n')

if run_gfortran_phdf5_petsc_opt:
    # run gfortran-phdf5-petsc-opt tests
    sys.stdout.write('\n--------------\n')
    sys.stdout.write('gfortran-phdf5-petsc-opt tests')
    sys.stdout.write('\n--------------\n')
    os.chdir(pwd)
    os.rename(pwd + '/../src/openmc-gfortran-phdf5-petsc-opt',
              pwd + '/../src/openmc')
    result_gfortran_phdf5_petsc_opt = nose.run(argv=[
        'run_tests.py', '-v', '--exclude', 'test_compile',
        '--mpi-np', '3', '--mpi-exec', mpiexec], addplugins=[NoseMPI()])
    os.rename(pwd + '/../src/openmc',
              pwd + '/../src/openmc-gfortran-phdf5-petsc-opt')
    if not result_gfortran_phdf5_petsc_opt:
        sys.stdout.write('\nDid not pass gfortran-phdf5-petsc-opt tests\n')

# print out summary of results
sys.stdout.write('\n======================================================\n')
sys.stdout.write('\nSummary of Compilation Option Testing:\n')

if run_compile:
    sys.stdout.write('\ncompilation tests........................')
    if result_compile:
        sys.stdout.write(BOLD + OK + '[OK]' + ENDC)
    else:
        sys.stdout.write(BOLD + FAIL + '[FAILED]' + ENDC)

if run_gfortran:
    sys.stdout.write('\ngfortran tests...........................')
    if result_gfortran:
        sys.stdout.write(BOLD + OK + '[OK]' + ENDC)
    else:
        sys.stdout.write(BOLD + FAIL + '[FAILED]' + ENDC)

if run_gfortran_dbg:
    sys.stdout.write('\ngfortran-DEBUG tests.....................')
    if result_gfortran_dbg:
        sys.stdout.write(BOLD + OK + '[OK]' + ENDC)
    else:
        sys.stdout.write(BOLD + FAIL + '[FAILED]' + ENDC)

if run_gfortran_opt:
    sys.stdout.write('\ngfortran-OPTIMIZE tests..................')
    if result_gfortran_opt:
        sys.stdout.write(BOLD + OK + '[OK]' + ENDC)
    else:
        sys.stdout.write(BOLD + FAIL + '[FAILED]' + ENDC)

if run_gfortran_hdf5:
    sys.stdout.write('\ngfortran-HDF5 tests......................')
    if result_gfortran_hdf5:
        sys.stdout.write(BOLD + OK + '[OK]' + ENDC)
    else:
        sys.stdout.write(BOLD + FAIL + '[FAILED]' + ENDC)

if run_gfortran_mpi:
    sys.stdout.write('\ngfortran-MPI tests.......................')
    if result_gfortran_mpi:
        sys.stdout.write(BOLD + OK + '[OK]' + ENDC)
    else:
        sys.stdout.write(BOLD + FAIL + '[FAILED]' + ENDC)

if run_gfortran_phdf5:
    sys.stdout.write('\ngfortran-MPI-HDF5 tests..................')
    if result_gfortran_phdf5:
        sys.stdout.write(BOLD + OK + '[OK]' + ENDC)
    else:
        sys.stdout.write(BOLD + FAIL + '[FAILED]' + ENDC)

if run_gfortran_petsc:
    sys.stdout.write('\ngfortran-MPI-PETSC tests.................')
    if result_gfortran_petsc:
        sys.stdout.write(BOLD + OK + '[OK]' + ENDC)
    else:
        sys.stdout.write(BOLD + FAIL + '[FAILED]' + ENDC)

if run_gfortran_phdf5_petsc:
    sys.stdout.write('\ngfortran-MPI-HDF5-PETSC tests............')
    if result_gfortran_phdf5_petsc:
        sys.stdout.write(BOLD + OK + '[OK]' + ENDC)
    else:
        sys.stdout.write(BOLD + FAIL + '[FAILED]' + ENDC)

if run_gfortran_phdf5_petsc_opt:
    sys.stdout.write('\ngfortran-MPI-PETSC-HDF5-OPTIMIZE tests...')
    if result_gfortran_phdf5_petsc_opt:
        sys.stdout.write(BOLD + OK + '[OK]' + ENDC + '\n')
    else:
        sys.stdout.write(BOLD + FAIL + '[FAILED]' + ENDC)

sys.stdout.write('\n')
