#!/usr/bin/env python

from sys import argv
from struct import unpack
from math import sqrt

import numpy as np
import scipy.stats
import matplotlib.pyplot as plt

from statepoint import StatePoint

# Get filename
filename = argv[1]

# Determine axial level
axial_level = int(argv[2]) if len(argv) > 2 else 1
score = int(argv[3]) if len(argv) > 3 else 1

# Create StatePoint object
sp = StatePoint(filename)

# Calculate t-value for 95% two-sided CI
n = sp.current_batch - sp.n_inactive
t_value = scipy.stats.t.ppf(0.975, n - 1)

# Loop over all tallies
print("Reading data...")
for t in sp.tallies:
    n_bins = t.n_score_bins * t.n_filter_bins
    i_mesh = [f.type for f in t.filters].index('mesh')

    # Get Mesh object
    m = sp.meshes[t.filters[i_mesh].bins[0] - 1]
    nx, ny, nz = m.dimension
    ns = t.n_score_bins * t.n_filter_bins / (nx*ny*nz)

    assert n_bins == nx*ny*nz*ns

    # Create lists for tallies
    mean = np.zeros((nx,ny))
    error = np.zeros((nx,ny))
    criteria = np.zeros((nx,ny))

    # Determine starting position for data
    start = sp._f.tell() + (axial_level-1)*ns*16 + (score - 1)*16

    for x in range(nx):
        for y in range(ny):
            # Seek to position of data
            sp._f.seek(start + x*ny*nz*ns*16 + y*nz*ns*16)

            # Read sum and sum-squared
            s, s2 = unpack('=2d', sp._f.read(16))
            s /= n
            mean[x,y] = s
            if s != 0.0:
                error[x,y] = t_value*sqrt((s2/n - s*s)/(n-1))/s
                criteria[x,y] = 1.0 if error[x,y] < 0.05 else 0.0

    # Make figure
    print("Making colorplot...")
    plt.imshow(mean, interpolation="nearest")
    plt.colorbar()
    plt.xlim((0,nx))
    plt.ylim((0,ny))
    plt.xticks(np.linspace(0,nx,5),
               np.linspace(m.lower_left[0],m.upper_right[0],5))
    plt.yticks(np.linspace(0,ny,5),
               np.linspace(m.lower_left[1],m.upper_right[1],5))
    plt.show()
