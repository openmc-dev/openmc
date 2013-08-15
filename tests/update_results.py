#!/usr/bin/env python

import os
import sys
from subprocess import Popen, call, STDOUT, PIPE
from glob import glob


class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'

    def disable(self):
        self.HEADER = ''
        self.OKBLUE = ''
        self.OKGREEN = ''
        self.WARNING = ''
        self.FAIL = ''
        self.ENDC = ''

# Check for arguments
if len(sys.argv) > 1:
    folders = []
    for i in range(len(sys.argv)):
        if i == 0:
            continue
        folders.append((sys.argv[i], ' ', ' '))
else:
    folders = os.walk('.')

# Loop around directories
for root, dirs, files in folders:

    # Delete first two characters
    if root[0:2] == './':
        root = root[2:]

    # Don't continue if not prefixed by test
    if root[0:4] != 'test':
        continue

    # Skip test compile or plot
    if root == 'test_compile' or root.find('plot') != -1:
        continue

    # Go into that directory
    os.chdir(root)
    pwd = os.path.abspath(os.path.dirname('settings.xml'))
    os.putenv('PWD', pwd)

    # Print status to screen
    sys.stdout.write(root+'... ')

    # Run openmc
    proc = Popen(['../../src/openmc'], stderr=STDOUT, stdout=PIPE)
    returncode = proc.wait()
    if returncode == 0:
        sys.stdout.write(bcolors.OKGREEN + "OK" + bcolors.ENDC + "\n")
    else:
        sys.stdout.write(bcolors.FAIL + "ERROR" + bcolors.ENDC + "\n")

    # Process results
    call(['python', 'results.py'])
    os.rename('results_test.dat', 'results_true.dat')
    for path in glob("*.binary") + glob("*.out"):
        os.remove(path)

    # Go back a directory
    os.chdir('..')
