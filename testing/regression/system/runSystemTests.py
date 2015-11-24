#!/usr/bin/env python

import os
from subprocess import call

n = len(os.listdir('.'))
ndirs = 0
for i in range(0,n):
  if os.path.isdir(os.listdir('.')[i]):
    ndirs += 1

for i in range(1,ndirs+1):
  os.chdir('./case-'+str(i))
  call(['echo', 'Running system test case-'+str(i)+'...'])
  call(['/home/walshjon/dev/shutmc/src/build/bin/shutmc', '.'])
  os.chdir('..')
