#!/usr/bin/env python

from sys import argv
from struct import unpack
from math import sqrt

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
for i, t in enumerate(sp.tallies):
    # Create lists for tallies
    mean = []
    uncertainties = []
    nonzero = []

    n_bins = t.n_score_bins * t.n_filter_bins

    # Loop over filter/score bins
    for j in range(n_bins):
        # Read sum and sum-squared
        s, s2 = unpack('=2d', sp._f.read(16))
        s /= n
        if s != 0.0:
            relative_error = t_value*sqrt((s2/n - s*s)/(n-1))/s
            uncertainties.append(relative_error)

    # Display number of tallies with less than 1% CIs 
    fraction = (sum([1.0 if re < 0.01 else 0.0 for re in uncertainties]) 
                /n_bins)
    print("Tally " + str(i+1))
    print("  Fraction under 1% = {0}".format(fraction))
    print("  Min relative error = {0}".format(min(uncertainties)))
    print("  Max relative error = {0}".format(max(uncertainties)))
    print("  Non-scoring bins = {0}".format(1.0 - float(len(uncertainties))/n_bins))

    # Plot histogram
    plt.hist(uncertainties, 100)
    plt.show()
