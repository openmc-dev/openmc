#!/usr/bin/env python

import os
import sys
from subprocess import Popen, STDOUT, PIPE

# Loop around directories
for root, dirs, files in os.walk('.'):

  # Don't continue if not prefixed by test
  if root[2:6] != 'test':
    continue
 
  # Go into that directory
  os.chdir(root)
  pwd = os.path.abspath(os.path.dirname('settings.xml'))
  os.putenv('PWD', pwd)

  # Run openmc
  proc = Popen(['../../src/openmc'])
  returncode = proc.wait()

  # Process results
  os.system('python results.py') 
  os.system('mv results_test.dat results_true.dat')
  os.system('rm *.binary')

  # Go back a directory
  os.chdir('..')
