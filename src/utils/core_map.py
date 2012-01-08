#!/usr/bin/env python

import sys

import numpy as np
import matplotlib.pyplot as pyplot

# Get filename and file object
filename = sys.argv[1]
axial_level = int(sys.argv[2])
if len(sys.argv) > 3:
    stage = int(sys.argv[3])
else:
    stage = 1
fh = open(filename,'r')

# Read size of mesh
words = fh.readline().split()
nx, ny, nz, nstage = [int(item) for item in words]

# Read values
value_dict = {}
for line in fh:
    words = line.split()
    i, j, k, m = [int(item) for item in words[:4]]
    value = float(words[-1])
    value_dict[i,j,k,m] = value

# Set up matrix
matrix = np.array([[value_dict[i+1,j+1,axial_level,stage] for i in range(nx)]
                   for j in range(ny)])

# Make figure
pyplot.pcolor(matrix)
pyplot.colorbar()
pyplot.xlim([0,nx])
pyplot.ylim([0,ny])
pyplot.xticks([])
pyplot.yticks([])
pyplot.show()
