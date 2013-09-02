#!/usr/bin/env python

from __future__ import print_function

import os
import sys
import nose
import glob
from subprocess import call

from nose_mpi import NoseMPI


def run_compile():
    # clean up all previous executables
    openmc_exe = glob.glob(pwd + '/../src/openmc*')
    for exe in openmc_exe:
        os.remove(exe)

    # run compile test
    result = nose.run(argv=['run_tests.py', 'test_compile'] + flags)
    if not result:
        print('Did not pass compile tests.')
    results.append(('compile', result))


def run_suite(name=None, mpi=False, cmfd=True):
    print('-'*(len(name) + 6))
    print(name + ' tests')
    print('-'*(len(name) + 6))

    # Set arguments list. Note that the first argument is a dummy argument (the
    # script name). It's not actually recursively calling run_tests.py
    argv = ['run_tests.py', '--exclude', 'test_compile'] + flags

    # Add MPI plugin if set
    if mpi:
        plugins = [NoseMPI()]
        argv += ['--mpi-np', '3', '--mpi-exec', mpiexec]
    else:
        plugins = None

    # Exclude CMFD
    if not cmfd:
        argv += ['--exclude', 'cmfd']

    try:
        os.chdir(pwd)
        os.rename(pwd + '/../src/openmc-' + name, pwd + '/../src/openmc')
        result = nose.run(argv=argv, addplugins=plugins)
    finally:
        os.rename(pwd + '/../src/openmc', pwd + '/../src/openmc-' + name)
        if not result:
            print('Did not pass ' + name + ' tests')
        results.append((name, result))

# set mpiexec path
mpiexec = '/opt/mpich/3.0.4-gnu/bin/mpiexec'

# get current working directory
pwd = os.getcwd()
sys.path.append(pwd)

# Set list of tests, either default or from command line
if len(sys.argv) > 1:
    flags = [i for i in sys.argv[1:] if i.startswith('-')]
    tests = [i for i in sys.argv[1:] if not i.startswith('-')]
else:
    flags = []
    tests = ['compile', 'gfortran', 'gfortran-dbg', 'gfortran-opt',
             'gfortran-hdf5', 'gfortran-mpi', 'gfortran-phdf5',
             'gfortran-petsc', 'gfortran-phdf5-petsc',
             'gfortran-phdf5-petsc-opt']

# Run tests
results = []
for name in tests:
    if name == 'compile':
        run_compile()
    elif name in ['gfortran', 'gfortran-dbg', 'gfortran-opt', 'gfortran-hdf5']:
        run_suite(name=name, cmfd=False)
    elif name in ['gfortran-mpi', 'gfortran-phdf5']:
        run_suite(name=name, mpi=True, cmfd=False)
    elif name in ['gfortran-petsc', 'gfortran-phdf5-petsc',
                  'gfortran-phdf5-petsc-opt']:
        run_suite(name=name, mpi=True)

# print out summary of results
print('\n' + '='*54)
print('Summary of Compilation Option Testing:\n')

OK = '\033[92m'
FAIL = '\033[91m'
ENDC = '\033[0m'
BOLD = '\033[1m'
for name, result in results:
    print(name + '.'*(50 - len(name)), end='')
    if result:
        print(BOLD + OK + '[OK]' + ENDC)
    else:
        print(BOLD + FAIL + '[FAILED]' + ENDC)
