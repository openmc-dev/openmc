#!/usr/bin/env python

from math import sqrt
from glob import glob

import scipy.stats
import matplotlib.pyplot as plt

from statepoint import StatePoint

# Find all statepoints in this directory.
files = glob('./statepoint.*.binary')
# Arrange the file list in increasing batch order
files.sort()

# Initialize arrays as needed
mean = [None for x in range(len(files))]
uncert = [None for x in range(len(files))]
#tally_score = []
#tally_filter = []
active_batches = [None for x in range(len(files))]

for i_batch in xrange(len(files)):
    
    # Get filename    
    batch_filename = files[i_batch]
    
    # Create StatePoint object
    sp = StatePoint(batch_filename)
    
    # Check if tallies are present
    if not sp._get_int()[0]:
        raise Exception("No tally data in state point!")
    
    # Increase the dimensionality of our main variables
    mean[i_batch] = [None for x in range(len(sp.tallies))]
    uncert[i_batch] = [None for x in range(len(sp.tallies))]
    
    # Calculate t-value for 95% two-sided CI
    n = sp.current_batch - sp.n_inactive
    t_value = scipy.stats.t.ppf(0.975, n - 1)
    
    # Store the batch count    
    active_batches[i_batch] = n
    
    # Loop over all tallies
    for i_tally, t in enumerate(sp.tallies):
        # Resize the 2nd dimension      
        mean[i_batch][i_tally] = [None for x in range(t.n_score_bins)]
        uncert[i_batch][i_tally] = [None for x in range(t.n_score_bins)]
        
        for i_score in range(t.n_score_bins):
            # Resize the 3rd dimension            
            mean[i_batch][i_tally][i_score] = \
                [None for x in range(t.n_filter_bins)]
            uncert[i_batch][i_tally][i_score] = \
                [None for x in range(t.n_filter_bins)]
            for i_filter in range(t.n_filter_bins):
                s, s2 = sp._get_double(2)
                s /= n
                mean[i_batch][i_tally][i_score][i_filter] = s
                if s != 0.0:
                    relative_error = t_value*sqrt((s2/n - s*s)/(n-1))/s
                else:
                    relative_error = 0.0
                uncert[i_batch][i_tally][i_score][i_filter] = relative_error
       
        # Still need to get some information just to be able to build labels
            

# Reorder the data lists in to a list order more conducive for plotting:
# The indexing should be: [tally][score][filter][batch]
meanPlot = [None for x in range(len(mean[0]))] # Set to the number of tallies
uncertPlot = [None for x in range(len(mean[0]))] # Set to the number of tallies

# Get and set the correct sizes for the rest of the dimensions
for i_tally in range(len(meanPlot)):
    # Set 2nd (score) dimension
    meanPlot[i_tally] = [None for x in range(len(mean[0][i_tally]))]
    uncertPlot[i_tally] = [None for x in range(len(mean[0][i_tally]))]
    
    for i_score in range(len(meanPlot[i_tally])):
        # Set 3rd (filter) dimension
        meanPlot[i_tally][i_score] = \
            [None for x in range(len(mean[0][i_tally][i_score]))]
        uncertPlot[i_tally][i_score] = \
            [None for x in range(len(mean[0][i_tally][i_score]))]
        
        for i_filter in range(len(meanPlot[i_tally][i_score])):
            # Set 4th (batch) dimension
            meanPlot[i_tally][i_score][i_filter] = \
                [None for x in range(len(mean))]
            uncertPlot[i_tally][i_score][i_filter] = \
                [None for x in range(len(mean))]
        
# Now rearrange the data as suitable
for i_batch in range(len(mean)):
    for i_tally in range(len(mean[i_batch])):
        for i_score in range(len(mean[i_batch][i_tally])):
            for i_filter in range(len(mean[i_batch][i_tally][i_score])):
                meanPlot[i_tally][i_score][i_filter][i_batch] = \
                    mean[i_batch][i_tally][i_score][i_filter]
                uncertPlot[i_tally][i_score][i_filter][i_batch] = \
                    uncert[i_batch][i_tally][i_score][i_filter]

# Set plotting constants
xLabel = "Active Batches"

# Begin plotting
for i_tally in range(len(meanPlot)):
    # Set tally string (placeholder until I put tally labels in statePoint)
    tallyStr = "Tally " + str(i_tally)
    
    for i_score in range(len(meanPlot[i_tally])):
        # Set score string
        scoreStr = "" #Not in place until I get that from statePoint above
                
        for i_filter in range(len(meanPlot[i_tally][i_score])):
            # Set filter string
            filterStr = "" #Not in place until I get that from statePoint above
            
            # set Title 
            title = tallyStr + " for " + filterStr
            # set yLabel
            yLabel = scoreStr
            
            # Plot mean with error bars
            plt.errorbar(active_batches, \
                meanPlot[i_tally][i_score][i_filter][:], \
                uncertPlot[i_tally][i_score][i_filter][:])
            plt.xlabel(xLabel)
            plt.ylabel(yLabel)
            plt.title(title)
            plt.show()
            
            # Plot mean with error bars
            plt.plot(active_batches, \
                uncertPlot[i_tally][i_score][i_filter][:])
            plt.xlabel(xLabel)
            plt.ylabel(yLabel)
            plt.title("Error of " + title)
            plt.show()