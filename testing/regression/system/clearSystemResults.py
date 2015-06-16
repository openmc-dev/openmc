#!/usr/bin/env python

import os
from glob import glob

n = len(os.listdir('.'))
ndirs = 0
for i in range(0,n):
  if os.path.isdir(os.listdir('.')[i]):
    ndirs += 1

for i in range(1,ndirs+1):
  os.chdir('./case-'+str(i))
  files = glob('statepoint*')
  for file in files:
    os.remove('./'+file)
  files = glob('error.log')
  for file in files:
    os.remove('./'+file)
  os.chdir('..')

files = glob('system-suite.pbs.*')
for file in files:
  os.remove('./'+file)
files = glob('testSystemResults.dat')
for file in files:
  os.remove('./'+file)
