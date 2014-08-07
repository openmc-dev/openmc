#!/usr/bin/env python

from __future__ import print_function, division

from sys import argv
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

# Read number of realizations for global tallies
sp.n_realizations = sp._get_int()[0]

# Read global tallies
n_global_tallies = sp._get_int()[0]
sp.global_tallies = np.array(sp._get_double(2*n_global_tallies))
sp.global_tallies.shape = (n_global_tallies, 2)

# Flag indicating if tallies are present
tallies_present = sp._get_int()[0]

# Check if tallies are present
if not tallies_present:
    raise Exception("No tally data in state point!")

# Loop over all tallies
print("Reading data...")
for t in sp.tallies:
    # Calculate t-value for 95% two-sided CI
    n = t.n_realizations
    t_value = scipy.stats.t.ppf(0.975, n - 1)

    n_bins = t.total_score_bins * t.total_filter_bins

    # Check for mesh
    if 'mesh' not in t.filters:
        continue

    # Get Mesh object and determine size
    m = sp.meshes[t.filters['mesh'].bins[0] - 1]
    if len(m.dimension) == 2:
        nx, ny = m.dimension
        nz = 1
    else:
        nx, ny, nz = m.dimension

    # Calculate number of score bins
    ns = t.total_score_bins * t.total_filter_bins // (nx*ny*nz)
    assert n_bins == nx*ny*nz*ns

    # Create lists for tallies
    mean = np.zeros((nx, ny))
    error = np.zeros((nx, ny))
    criteria = np.zeros((nx, ny))

    # Determine starting position for data
    start = sp._f.tell() + (axial_level-1)*ns*16 + (score - 1)*16

    for x in range(nx):
        for y in range(ny):
            # Seek to position of data
            sp._f.seek(start + x*ny*nz*ns*16 + y*nz*ns*16)

            # Read sum and sum-squared
            s, s2 = sp._get_double(2)
            s /= n
            mean[x, y] = s
            if s != 0.0:
                error[x, y] = t_value*sqrt((s2/n - s*s)/(n-1))/s
                criteria[x, y] = 1.0 if error[x, y] < 0.05 else 0.0

    # Make figure
    print("Making colorplot...")
    plt.imshow(mean, interpolation="nearest")
    plt.colorbar()
    plt.xlim((0, nx))
    plt.ylim((0, ny))
    plt.xticks(np.linspace(0, nx, 5),
               np.linspace(m.lower_left[0], m.upper_right[0], 5))
    plt.yticks(np.linspace(0, ny, 5),
               np.linspace(m.lower_left[1], m.upper_right[1], 5))
    plt.show()
