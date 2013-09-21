#!/usr/bin/env python

import os
import sys
from subprocess import Popen, call, STDOUT, PIPE
from glob import glob

OKGREEN = '\033[92m'
FAIL = '\033[91m'
ENDC = '\033[0m'
BOLD = '\033[1m'

# Check for arguments
if len(sys.argv) > 1:
    folders = []
    for i in range(len(sys.argv)):
        if i == 0:
            continue
        folders.append(sys.argv[i])
else:
    folders = glob('test_*')

# Loop around directories
for adir in sorted(folders):

    # Skip test compile or plot
    if adir == 'test_compile' or adir.find('plot') != -1:
        continue

    # Go into that directory
    os.chdir(adir)
    pwd = os.path.abspath(os.path.dirname('settings.xml'))
    os.putenv('PWD', pwd)

    # Print status to screen
    sys.stdout.write(adir)
    sz = len(adir)
    for i in range(35 - sz):
        sys.stdout.write('.')

    # Run openmc
    proc = Popen(['../../src/openmc'], stderr=STDOUT, stdout=PIPE)
    returncode = proc.wait()
    if returncode == 0:
        sys.stdout.write(BOLD + OKGREEN + "[OK]" + ENDC + "\n")
    else:
        sys.stdout.write(BOLD + FAIL + "[FAILED]" + ENDC + "\n")

    # Process results
    if returncode == 0:
        call(['python', 'results.py'])
        os.rename('results_test.dat', 'results_true.dat')
        for path in glob("*.binary") + glob("*.out"):
            os.remove(path)

    # Go back a directory
    os.chdir('..')
