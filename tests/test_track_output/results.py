#!/usr/bin/env python

import os
import sys
import glob
import shutil
from subprocess import call

# If vtk python module is not available, we can't run track.py so skip this
# test
cwd = os.getcwd()
try:
    import vtk
except ImportError:
    print('----------------Skipping test-------------')
    shutil.copy('results_true.dat', 'results_test.dat')
    exit()

# Run track processing script
call(['../../src/utils/track.py', '-o', 'poly'] +
     glob.glob(''.join((cwd, '/track*'))))
poly = ''.join((cwd, '/poly.pvtp'))
assert os.path.isfile(poly), 'poly.pvtp file not found.'
shutil.copy('poly.pvtp', 'results_test.dat')
