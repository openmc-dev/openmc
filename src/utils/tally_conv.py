#!/usr/bin/env python

from math import sqrt, pow
from glob import glob

import scipy.stats
import matplotlib.pyplot as plt

from statepoint import StatePoint

# USER OPTIONS

# Set filetype (the file extension desired, without the period.)
# Options are backend dependent, but most backends support png, pdf, ps, eps 
# and svg.  Write "none" if no figures are desired.
fileType = "none"

# Set if cross-sections of reaction rates are desired
printxsFlag = True

# Set if it is desired to show the images
showImg = True

# END USER OPTIONS

# Find all statepoints in this directory.
files = glob('./statepoint.*.binary')
# Arrange the file list in increasing batch order
files.sort()

# Initialize arrays as needed
mean = [None for x in range(len(files))]
uncert = [None for x in range(len(files))]
scoreType = [None for x in range(len(files))]
filterType = [None for x in range(len(files))]
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
    scoreType[i_batch] = [None for x in range(len(sp.tallies))]    
    filterType[i_batch] = [None for x in range(len(sp.tallies))]    
    
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
        scoreType[i_batch][i_tally] = [None for x in range(t.n_score_bins)]
        filterType[i_batch][i_tally] = [None for x in range(t.n_score_bins)]
        
        for i_score in range(t.n_score_bins):
            # Resize the 3rd dimension            
            mean[i_batch][i_tally][i_score] = \
                [None for x in range(t.n_filter_bins)]
            uncert[i_batch][i_tally][i_score] = \
                [None for x in range(t.n_filter_bins)]
            scoreType[i_batch][i_tally][i_score] = t.scores[i_score]  
            filterType[i_batch][i_tally][i_score] = \
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
            #scoreType[i_batch][i_tally][i_score][i_filter] = t.scores[i_score]            
            # Not quite right:
            filterType[i_batch][i_tally][i_score][i_filter] = t.filters[i_filter]


# Reorder the data lists in to a list order more conducive for plotting:
# The indexing should be: [tally][score][filter][batch]
meanPlot = [None for x in range(len(mean[0]))] # Set to the number of tallies
uncertPlot = [None for x in range(len(mean[0]))] # Set to the number of tallies
absUncertPlot = [None for x in range(len(mean[0]))] # Set to the number of tallies
fluxLoc = [None for x in range(len(mean[0]))] # Set to the number of tallies
printxs = [False for x in range(len(mean[0]))] # Set to the number of tallies

# Get and set the correct sizes for the rest of the dimensions
for i_tally in range(len(meanPlot)):
    # Set 2nd (score) dimension
    meanPlot[i_tally] = [None for x in range(len(mean[0][i_tally]))]
    uncertPlot[i_tally] = [None for x in range(len(mean[0][i_tally]))]
    absUncertPlot[i_tally] = [None for x in range(len(mean[0][i_tally]))]
    
    # Initialize flux location so it will be -1 if not found
    fluxLoc[i_tally] = -1
    
    for i_score in range(len(meanPlot[i_tally])):
        # Set 3rd (filter) dimension
        meanPlot[i_tally][i_score] = \
            [None for x in range(len(mean[0][i_tally][i_score]))]
        uncertPlot[i_tally][i_score] = \
            [None for x in range(len(mean[0][i_tally][i_score]))]
        absUncertPlot[i_tally][i_score] = \
            [None for x in range(len(mean[0][i_tally][i_score]))]
        
        for i_filter in range(len(meanPlot[i_tally][i_score])):
            # Set 4th (batch) dimension
            meanPlot[i_tally][i_score][i_filter] = \
                [None for x in range(len(mean))]
            uncertPlot[i_tally][i_score][i_filter] = \
                [None for x in range(len(mean))]
            absUncertPlot[i_tally][i_score][i_filter] = \
                [None for x in range(len(mean))]
        
        # Set flux location if found
        if scoreType[0][i_tally][i_score] == 'flux':
            fluxLoc[i_tally] = i_score

# Set printxs array according to the printXSflag input
if printxsFlag:
    for i_tally in range(len(fluxLoc)):
        if fluxLoc[i_tally] != -1:
            printxs[i_tally] = True
        
# Now rearrange the data as suitable, and perform xs conversion if necessary
for i_batch in range(len(mean)):
    for i_tally in range(len(mean[i_batch])):
        for i_score in range(len(mean[i_batch][i_tally])):
            for i_filter in range(len(mean[i_batch][i_tally][i_score])):
                if (printxs[i_tally] and \
                    ((scoreType[0][i_tally][i_score] != 'flux') and \
                    (scoreType[0][i_tally][i_score] != 'current'))):
                    
                    # Perform rate to xs conversion
                    # mean is mean/fluxmean
                    meanPlot[i_tally][i_score][i_filter][i_batch] = \
                        mean[i_batch][i_tally][i_score][i_filter] / \
                        mean[i_batch][i_tally][fluxLoc[i_tally]][i_filter]
                    
                    # Update the relative uncertainty via error propagation
                    uncertPlot[i_tally][i_score][i_filter][i_batch] = \
                        sqrt(pow(uncert[i_batch][i_tally][i_score][i_filter],2.0) \
                        + pow(uncert[i_batch][i_tally][fluxLoc[i_tally]][i_filter],2.0))
                else: 
                    
                    # Do not perform rate to xs conversion
                    meanPlot[i_tally][i_score][i_filter][i_batch] = \
                        mean[i_batch][i_tally][i_score][i_filter]
                    uncertPlot[i_tally][i_score][i_filter][i_batch] = \
                        uncert[i_batch][i_tally][i_score][i_filter]  
                
                # Both have the same absolute uncertainty calculation
                absUncertPlot[i_tally][i_score][i_filter][i_batch] = \
                    uncert[i_batch][i_tally][i_score][i_filter] * \
                    mean[i_batch][i_tally][i_score][i_filter]    

# Set plotting constants
xLabel = "Active Batches"
xLabel = xLabel.title() # not necessary for now, but is left in to handle if 
# the previous line changes

# Begin plotting
for i_tally in range(len(meanPlot)):
    # Set tally string (placeholder until I put tally labels in statePoint)
    tallyStr = "Tally " + str(i_tally + 1)
    
    for i_score in range(len(meanPlot[i_tally])):
        # Set score string
        scoreStr = scoreType[i_batch][i_tally][i_score]
        scoreStr = scoreStr.title()
        if (printxs[i_tally] and ((scoreStr != 'Flux') and (scoreStr != 'Current'))):
            scoreStr = scoreStr + "-XS"
                
        for i_filter in range(len(meanPlot[i_tally][i_score])):
            # Set filter string
            filterStr = "Filter " + str(i_filter + 1)
            
            # set Title 
            title = "Convergence of " + scoreStr + " in " + tallyStr + " for " + filterStr
            
            # set yLabel
            yLabel = scoreStr
            yLabel = yLabel.title()

            # Set saving filename
            fileName = "tally_" + str(i_tally + 1) + "_" + scoreStr + \
                "_filter_" + str(i_filter+1) + "." + fileType
            REfileName = "tally_" + str(i_tally + 1) + "_" + scoreStr + \
                "RE_filter_" + str(i_filter+1) + "." + fileType
            
            # Plot mean with absolute error bars
            plt.errorbar(active_batches, \
                meanPlot[i_tally][i_score][i_filter][:], \
                absUncertPlot[i_tally][i_score][i_filter][:])
            plt.xlabel(xLabel)
            plt.ylabel(yLabel)
            plt.title(title)
            if (fileType != 'none'):
                plt.savefig(fileName)
            if showImg:            
                plt.show()
            plt.clf()
            
            # Plot relative uncertainty
            plt.plot(active_batches, \
                uncertPlot[i_tally][i_score][i_filter][:])
            plt.xlabel(xLabel)
            plt.ylabel("Relative Error of " + yLabel)
            plt.title("Relative Error of " + title)
            if (fileType != 'none'):
                plt.savefig(REfileName)
            if showImg:            
                plt.show()
            plt.clf()