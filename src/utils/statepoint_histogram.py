#!/usr/bin/env python

from __future__ import print_function
from sys import argv
from math import sqrt

import numpy as np
import scipy.stats
import matplotlib.pyplot as plt

from openmc.statepoint import StatePoint

# Get filename
filename = argv[1]

# Create StatePoint object
sp = StatePoint(filename)
sp.read_results()
sp.compute_ci()

# Check if tallies are present
if not sp.tallies_present:
    raise Exception("No tally data in state point!")

# Loop over all tallies
for i, t in sp.tallies.items():
    # Determine relative error and fraction of bins with less than 1% half-width
    # of CI
    n_bins = t.mean.size
    relative_error = t.std_dev[t.mean > 0.] / t.mean[t.mean > 0.]
    fraction = float(sum(relative_error < 0.01))/n_bins

    # Display results
    print("Tally " + str(i))
    print("  Fraction under 1% = {0}".format(fraction))
    print("  Min relative error = {0}".format(min(relative_error)))
    print("  Max relative error = {0}".format(max(relative_error)))
    print("  Non-scoring bins = {0}".format(
          1.0 - float(relative_error.size)/n_bins))

    # Plot histogram
    plt.hist(relative_error, 100)
    plt.show()
