#!/usr/bin/env python

"""
This script reads a cross_sections.out file, adds up the memory usage for each
nuclide and S(a,b) table, and displays the total memory usage.

"""

from __future__ import print_function

import sys
import os

if len(sys.argv) > 1:
    # Get path to cross_sections.out file from command line argument
    filename = sys.argv[-1]
else:
    # Set default path for cross_sections.out
    filename = 'cross_sections.out'
    if not os.path.exists(filename):
        raise OSError('Could not find cross_sections.out file!')

# Open file handle for cross_sections.out file
f = open(filename, 'r')

# Initialize memory size arrays
memory_xs = []
memory_angle = []
memory_energy = []
memory_urr = []
memory_total = []
memory_sab = []

while True:
    # Read next line in file
    line = f.readline()

    # Check for EOF
    if line == '':
        break

    # Look for block listing memory usage for a nuclide
    words = line.split()
    if len(words) == 2 and words[0] == 'Memory':
        memory_xs.append(int(f.readline().split()[-2]))
        memory_angle.append(int(f.readline().split()[-2]))
        memory_energy.append(int(f.readline().split()[-2]))
        memory_urr.append(int(f.readline().split()[-2]))
        memory_total.append(int(f.readline().split()[-2]))

    # Look for memory usage for S(a,b) table
    if len(words) == 5 and words[1] == 'Used':
        memory_sab.append(int(words[-2]))

# Write out summary memory usage
print('Memory Requirements')
print('  Reaction Cross Sections        = ' + str(sum(memory_xs)))
print('  Secondary Angle Distributions  = ' + str(sum(memory_angle)))
print('  Secondary Energy Distributions = ' + str(sum(memory_energy)))
print('  Probability Tables             = ' + str(sum(memory_urr)))
print('  S(a,b) Tables                  = ' + str(sum(memory_sab)))
print('  Total                          = ' + str(sum(memory_total)))
