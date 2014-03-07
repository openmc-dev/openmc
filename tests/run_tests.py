#!/usr/bin/env python

from __future__ import print_function

import os
import sys
import nose
import glob
from subprocess import call

from nose_mpi import NoseMPI


def run_compile():
    print('-'*17)
    print('Compilation tests')
    print('-'*17)

    # clean up all previous executables
    openmc_exe = glob.glob(pwd + '/../src/openmc*')
    for exe in openmc_exe:
        os.remove(exe)

    # run compile test
    result = nose.run(argv=['nosetests', 'test_compile'] + flags)
    if not result:
        print('Did not pass compile tests.')
    results.append(('compile', result))


def run_suite(name=None, mpi=False):
    print('-'*(len(name) + 6))
    print(name + ' tests')
    print('-'*(len(name) + 6))

    # Set arguments list. Note that the first argument is a dummy argument (the
    # script name). It's not actually recursively calling run_tests.py
    argv = ['nosetests', '--exclude', 'test_compile'] + flags

    # Add MPI plugin if set
    if mpi:
        plugins = [NoseMPI()]
        argv += ['--mpi-np', '3', '--mpi-exec', mpiexec]
    else:
        plugins = None

    try:
        os.chdir(pwd)
        os.rename(pwd + '/../src/openmc-' + name, pwd + '/../src/openmc')
        result = nose.run(argv=argv, addplugins=plugins)
        os.rename(pwd + '/../src/openmc', pwd + '/../src/openmc-' + name)
    except OSError:
        result = False
        print('No OpenMC executable found for ' + name + ' tests')
    if not result:
        print('Did not pass ' + name + ' tests')
    results.append((name, result))

# set mpiexec path
if 'COMPILER' in os.environ:
    compiler = os.environ['COMPILER']
else:
    compiler = 'gnu'
mpiexec = '/opt/mpich/3.0.4-{0}/bin/mpiexec'.format(compiler)

# get current working directory
pwd = os.getcwd()
sys.path.append(pwd)

# Set list of tests, either default or from command line
flags = []
tests = ['compile', 'normal', 'debug', 'optimize',
         'omp', 'omp-debug', 'omp-optimize',
         'hdf5', 'hdf5-debug', 'hdf5-optimize',
         'omp-hdf5', 'omp-hdf5-debug', 'omp-hdf5-optimize',
         'mpi', 'mpi-debug', 'mpi-optimize',
         'mpi-omp', 'mpi-omp-debug', 'mpi-omp-optimize',
         'phdf5', 'phdf5-debug', 'phdf5-optimize',
         'phdf5-omp', 'phdf5-omp-debug', 'phdf5-omp-optimize',
         'petsc', 'petsc-debug', 'petsc-optimize',
         'phdf5-petsc', 'phdf5-petsc-debug', 'phdf5-petsc-optimize',
         'omp-phdf5-petsc', 'omp-phdf5-petsc-debug',
         'omp-phdf5-petsc-optimize']
if len(sys.argv) > 1:
    flags = [i for i in sys.argv[1:] if i.startswith('-')]
    tests_ = [i for i in sys.argv[1:] if not i.startswith('-')]

    # Check for special subsets of tests
    tests__ = []
    for i in tests_:

        # All tests will run all the tests except for compile unless
        # it is also specified on the command line. Note that specifying
        # compile all-tests is the same as not specifying any args
        if i == 'all-tests':
            tests__ = tests
            try:
                idx = tests_.index('compile')  # check for compile test
            except ValueError:
                del tests__[0]
            finally:
                break  # don't need to check for anything else

        # This checks for any subsets of tests. The string after
        # all-XXXX will be used to search through all tests.
        # Specifying XXXX=normal will run tests that don't contain
        # debug or optimize substring.
        if i.startswith('all-'):
            suffix = i.split('all-')[1]
            if suffix == 'normal':
                for j in tests:
                    if j.rfind('debug') == -1 and \
                       j.rfind('optimize') == -1:
                        tests__.append(j)
            else:
                for j in tests:
                    if j.rfind(suffix) != -1:
                        if suffix == 'omp' and j == 'compile':
                            continue
                        if j == 'compile':
                            continue 
                        tests__.append(j)
        else:
            tests__.append(i)  # append specific test (e.g., mpi-debug)
    tests = tests__ if tests__ else tests

# Run tests
results = []
for name in tests:
    if name == 'compile':
        run_compile()
    elif name in ['normal', 'debug', 'optimize',
                  'hdf5', 'hdf5-debug', 'hdf5-optimize',
                  'omp', 'omp-debug', 'omp-optimize',
                  'omp-hdf5', 'omp-hdf5-debug', 'omp-hdf5-optimize']:
        run_suite(name=name)
    elif name in ['mpi', 'mpi-debug', 'mpi-optimize',
                  'mpi-omp', 'mpi-omp-debug', 'mpi-omp-optimize',
                  'phdf5', 'phdf5-debug', 'phdf5-optimize',
                  'phdf5-omp', 'phdf5-omp-debug', 'phdf5-omp-optimize',
                  'petsc', 'petsc-debug', 'petsc-optimize',
                  'phdf5-petsc', 'phdf5-petsc-debug', 'phdf5-petsc-optimize',
                  'omp-phdf5-petsc', 'omp-phdf5-petsc-debug',
                  'omp-phdf5-petsc-optimize']:
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
