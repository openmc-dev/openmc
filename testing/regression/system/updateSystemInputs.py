#!/usr/bin/env python

import os
import fileinput
import sys

case = sys.argv[1]
finput = sys.argv[2]
old = sys.argv[3]
new = sys.argv[4]

def replace(file, sold, snew):
  for line in fileinput.input(file, inplace=True):
    print line.replace(sold, snew).rstrip()

os.chdir('./case-'+case)
replace(finput, old, new) 
