#!/usr/bin/env python

import os
import sys
from subprocess import Popen, STDOUT, PIPE

# Check for arguments
if len(sys.argv) > 1:
    folders = []
    for i in range(len(sys.argv)):
        if i == 0: continue
        folders.append((sys.argv[i],' ',' '))
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
