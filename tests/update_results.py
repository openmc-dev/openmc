#!/usr/bin/env python

from __future__ import print_function

import os
import re
from subprocess import Popen, call, STDOUT, PIPE
from glob import glob
from optparse import OptionParser
import warnings

parser = OptionParser()
parser.add_option('--mpi_exec', dest='mpi_exec',
                  help="Path to mpiexec for tests that require it.")
parser.add_option('-R', '--tests-regex', dest='regex_tests',
                  help="Run tests matching regular expression. \
                  Test names are the directories present in tests folder.\
                  This uses standard regex syntax to select tests.")
(opts, args) = parser.parse_args()
cwd = os.getcwd()

# Terminal color configurations
OKGREEN = '\033[92m'
FAIL = '\033[91m'
ENDC = '\033[0m'
BOLD = '\033[1m'

# Get a list of all test folders
folders = glob('test_*')

# Check to see if a subset of tests is specified on command line
if opts.regex_tests is not None:
    folders = [item for item in folders if re.search(opts.regex_tests, item)]

# Loop around directories
for adir in sorted(folders):

    # Go into that directory
    os.chdir(adir)
    pwd = os.path.abspath(os.path.dirname('settings.xml'))
    os.putenv('PWD', pwd)

    # Print status to screen
    print(adir, end="")
    sz = len(adir)
    for i in range(35 - sz):
        print('.', end="")

    # Find the test executable
    test_exec = glob('test_*.py')
    assert len(test_exec) == 1, 'There must be only one test executable per ' \
         'test directory'

    # Update the test results
    proc = Popen(['python', test_exec[0], '--update'])
    # Handle DD test separately since it requires running OpenMC twice
    if adir == 'test_domain_decomp':
        if opts.mpi_exec is None:
            warnings.warn('Need to specify an mpi_exec executable for DD test')
            print(BOLD + FAIL + "[FAILED]" + ENDC)
            continue
        else:
            mpi_exec = os.path.abspath(opts.mpi_exec)
        if os.path.exists('results_error.dat'):
            os.remove('results_error.dat')
        proc = Popen(['python', 'test_domain_decomp.py', '--exe', openmc_exe,
                      '--mpi_exec', mpi_exec],
                     stderr=STDOUT, stdout=PIPE)
        returncode = proc.wait()
        if returncode != 0:
            print(BOLD + FAIL + "[FAILED]" + ENDC)
            continue
        os.chdir('1_domain')
        if os.path.exists('results_error.dat'):
            os.rename('results_error.dat', 'results_true.dat')
        os.chdir('..')
        os.chdir('4_domains')
        if os.path.exists('results_error.dat'):
            os.rename('results_error.dat', 'results_true.dat')
        os.chdir('..')
        print(BOLD + OKGREEN + "[OK]" + ENDC)
        os.chdir('..')
        continue

    # Handle distribcell test separately since it requires running OpenMC 3x
    if adir == 'test_filter_distribcell':
        if os.path.exists('results_error.dat'):
            os.remove('results_error.dat')
        proc = Popen(['python', 'test_filter_distribcell.py',
                      '--exe', openmc_exe], stderr=STDOUT, stdout=PIPE)
        returncode = proc.wait()
        if os.path.exists('results_error.dat'):
            print('renamed results_error!!!')
            os.rename('results_error.dat', 'results_true.dat')
        print(BOLD + OKGREEN + "[OK]" + ENDC)
        os.chdir('..')
        continue

    # Run openmc
    proc = Popen([openmc_exe], stderr=STDOUT, stdout=PIPE)
    returncode = proc.wait()
    if returncode == 0:
        print(BOLD + OKGREEN + "[OK]" + ENDC)
    else:
        print(BOLD + FAIL + "[FAILED]" + ENDC)

    # Go back a directory
    os.chdir('..')
