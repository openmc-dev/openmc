#!/usr/bin/env python

from __future__ import print_function

import os
import re
from subprocess import Popen, call, STDOUT, PIPE
from glob import glob
from optparse import OptionParser

parser = OptionParser()
parser.add_option('--exe', dest='exe',
                  help="Path to openmc executable with basic \
                        configuration options (no HDF5, no MPI, etc.)")
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

# Check for valid executable
if opts.exe is None:
    raise Exception('Need to specify an OpenMC executable')
else:
    openmc_exe = os.path.abspath(opts.exe)

# Get a list of all test folders
folders = glob('test_*')

# Check to see if a subset of tests is specified on command line
if opts.regex_tests is not None:
    folders = [item for item in folders if re.search(opts.regex_tests, item)]

# Loop around directories
for adir in sorted(folders):

    # Skip test compile or plot
    if adir.find('plot') != -1:
        continue

    # Go into that directory
    os.chdir(adir)
    pwd = os.path.abspath(os.path.dirname('settings.xml'))
    os.putenv('PWD', pwd)

    # Print status to screen
    print(adir, end="")
    sz = len(adir)
    for i in range(35 - sz):
        print('.', end="")

    # Run openmc
    proc = Popen([openmc_exe], stderr=STDOUT, stdout=PIPE)
    returncode = proc.wait()
    if returncode == 0:
        print(BOLD + OKGREEN + "[OK]" + ENDC)
    else:
        print(BOLD + FAIL + "[FAILED]" + ENDC)

    # Process results
    if returncode == 0:
        call(['python', 'results.py'])
        os.rename('results_test.dat', 'results_true.dat')
        for path in glob("*.binary") + glob("*.out"):
            os.remove(path)

    # Go back a directory
    os.chdir('..')
