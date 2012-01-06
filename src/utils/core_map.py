#!/usr/bin/env python

import sys

import numpy as np
import matplotlib.pyplot as pyplot

# Get filename and file object
filename = sys.argv[1]
axial_level = int(sys.argv[2])
fh = open(filename,'r')

# Read size of mesh
words = fh.readline().split()
nx, ny, nz = [int(item) for item in words]

# Read values
value_dict = {}
for line in fh:
    words = line.split()
    i, j, k = [int(item) for item in words[:3]]
    value = float(words[3])
    value_dict[i,j,k] = value

# Set up matrix
matrix = np.array([[value_dict[i+1,j+1,axial_level] for i in range(nx)]
                   for j in range(ny)])

# Make figure
pyplot.pcolor(matrix)
pyplot.colorbar()
pyplot.xlim([0,nx])
pyplot.ylim([0,ny])
pyplot.xticks([])
pyplot.yticks([])
pyplot.show()
