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
    result = nose.run(argv=['run_tests.py', 'test_compile'] + flags)
    if not result:
        print('Did not pass compile tests.')
    results.append(('compile', result))


def run_suite(name=None, mpi=False):
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
if os.environ.has_key('COMPILER'):
    compiler = os.environ['COMPILER']
else:
    compiler = 'gnu'
mpiexec = '/opt/mpich/3.0.4-{0}/bin/mpiexec'.format(compiler)

# get current working directory
pwd = os.getcwd()
sys.path.append(pwd)

# Set list of tests, either default or from command line
flags = []
tests = ['compile', 'normal', 'debug', 'optimize', 'hdf5', 'mpi', 'phdf5',
         'petsc', 'phdf5-petsc', 'phdf5-petsc-optimize']
if len(sys.argv) > 1:
    flags = [i for i in sys.argv[1:] if i.startswith('-')]
    tests_ = [i for i in sys.argv[1:] if not i.startswith('-')]
    tests = tests_ if tests_ else tests

# Run tests
results = []
for name in tests:
    if name == 'compile':
        run_compile()
    elif name in ['normal', 'debug', 'optimize', 'hdf5']:
        run_suite(name=name)
    elif name in ['mpi', 'phdf5', 'petsc', 'phdf5-petsc',
                  'phdf5-petsc-optimize']:
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
