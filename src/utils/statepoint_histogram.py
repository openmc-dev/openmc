#!/usr/bin/env python

from __future__ import print_function
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
for i, t in enumerate(sp.tallies):
    # Calculate t-value for 95% two-sided CI
    n = t.n_realizations
    t_value = scipy.stats.t.ppf(0.975, n - 1)

    # Create lists for tallies
    mean = []
    uncertainties = []
    nonzero = []

    n_bins = t.total_score_bins * t.total_filter_bins

    # Loop over filter/score bins
    for j in range(n_bins):
        # Read sum and sum-squared
        s, s2 = sp._get_double(2)
        s /= n
        if s != 0.0:
            relative_error = t_value*sqrt((s2/n - s*s)/(n-1))/s
            uncertainties.append(relative_error)

    # Display number of tallies with less than 1% CIs
    fraction = (sum([1.0 if re < 0.01 else 0.0 for re in uncertainties])
                / n_bins)
    print("Tally " + str(i+1))
    print("  Fraction under 1% = {0}".format(fraction))
    print("  Min relative error = {0}".format(min(uncertainties)))
    print("  Max relative error = {0}".format(max(uncertainties)))
    print("  Non-scoring bins = {0}".format(
          1.0 - float(len(uncertainties))/n_bins))

    # Plot histogram
    plt.hist(uncertainties, 100)
    plt.show()
