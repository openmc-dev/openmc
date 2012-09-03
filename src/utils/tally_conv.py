#!/usr/bin/env python

# TODO: 
# - Add savetoCSV functionality, or whatever text output would be best for mathgl, gnuplot or veusz
# - Add smart labels for filters??
# - Add dots to the data points to make them stand out more
# - Smooth lines?
# - Convert my single plt to an array of plots, or perhaps to using the Figure?

from math import sqrt, pow
from glob import glob

import scipy.stats
import matplotlib.pyplot as plt

from statepoint import StatePoint

#from read_inputXML import talliesXML

# USER OPTIONS

# Set filetype (the file extension desired, without the period.)
# Options are backend dependent, but most backends support png, pdf, ps, eps 
# and svg.  Write "none" if no figures are desired.
fileType = "none"

# Set if cross-sections of reaction rates are desired
printxsFlag = True

# Set if it is desired to show the images
showImg = True

# Save to CSV?
savetoCSV = True

# END USER OPTIONS

## Find if tallies.xml exists.
#if glob('./tallies.xml') != None:
#    # It exists
#    tallyData = talliesXML('tallies.xml')
#else: 
#    # It does not exist.
#    tallyData = None

# Find all statepoints in this directory.
files = glob('./statepoint.*.binary')
# Arrange the file list in increasing batch order
files.sort()

# Initialize arrays as needed
mean = [None for x in range(len(files))]
uncert = [None for x in range(len(files))]
scoreType = [None for x in range(len(files))]
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
    
    # Calculate t-value for 95% two-sided CI
    n = sp.current_batch - sp.n_inactive
    t_value = scipy.stats.t.ppf(0.975, n - 1)
    
    # Store the batch count    
    active_batches[i_batch] = n
    
    # Loop over all tallies
    for i_tally, t in enumerate(sp.tallies):
        # Resize the 2nd dimension      
        mean[i_batch][i_tally] = [None for x in range(t.n_filter_bins)]
        uncert[i_batch][i_tally] = [None for x in range(t.n_filter_bins)]
        scoreType[i_batch][i_tally] = [None for x in range(t.n_filter_bins)]
        
        for i_filter in range(t.n_filter_bins):
            # Resize the 3rd dimension            
            mean[i_batch][i_tally][i_filter] = \
                [None for x in range(t.n_score_bins)]
            uncert[i_batch][i_tally][i_filter] = \
                [None for x in range(t.n_score_bins)]
            scoreType[i_batch][i_tally][i_filter] = \
                [None for x in range(t.n_score_bins)]
            
            for i_score in range(t.n_score_bins):
                scoreType[i_batch][i_tally][i_filter][i_score] = t.scores[i_score]                 
                s, s2 = sp._get_double(2)
                s /= n
                mean[i_batch][i_tally][i_filter][i_score] = s
                if s != 0.0:
                    relative_error = t_value*sqrt((s2/n - s*s)/(n-1))/s
                else:
                    relative_error = 0.0
                uncert[i_batch][i_tally][i_filter][i_score] = relative_error

# Reorder the data lists in to a list order more conducive for plotting:
# The indexing should be: [tally][filter][score][batch]
meanPlot = [None for x in range(len(mean[0]))] # Set to the number of tallies
uncertPlot = [None for x in range(len(mean[0]))] # Set to the number of tallies
absUncertPlot = [None for x in range(len(mean[0]))] # Set to the number of tallies
filterLabel = [None for x in range(len(mean[0]))] # Set to the number of tallies
fluxLoc = [None for x in range(len(mean[0]))] # Set to the number of tallies
printxs = [False for x in range(len(mean[0]))] # Set to the number of tallies

# Get and set the correct sizes for the rest of the dimensions
for i_tally in range(len(meanPlot)):
    # Set 2nd (score) dimension
    meanPlot[i_tally] = [None for x in range(len(mean[0][i_tally]))]
    uncertPlot[i_tally] = [None for x in range(len(mean[0][i_tally]))]
    absUncertPlot[i_tally] = [None for x in range(len(mean[0][i_tally]))]
    filterLabel[i_tally] = [None for x in range(len(mean[0][i_tally]))]
    
    # Initialize flux location so it will be -1 if not found
    fluxLoc[i_tally] = -1
    
    for i_filter in range(len(meanPlot[i_tally])):
        # Set 3rd (filter) dimension
        meanPlot[i_tally][i_filter] = \
            [None for x in range(len(mean[0][i_tally][i_filter]))]
        uncertPlot[i_tally][i_filter] = \
            [None for x in range(len(mean[0][i_tally][i_filter]))]
        absUncertPlot[i_tally][i_filter] = \
            [None for x in range(len(mean[0][i_tally][i_filter]))]
        filterLabel[i_tally][i_filter] = \
            [None for x in range(len(mean[0][i_tally][i_filter]))]
        
        for i_score in range(len(meanPlot[i_tally][i_filter])):
            # Set 4th (batch) dimension
            meanPlot[i_tally][i_filter][i_score] = \
                [None for x in range(len(mean))]
            uncertPlot[i_tally][i_filter][i_score] = \
                [None for x in range(len(mean))]
            absUncertPlot[i_tally][i_filter][i_score] = \
                [None for x in range(len(mean))]
                
            # Get filterLabel (this should be moved to its own function)
            #??? How to do?
            
        # Set flux location if found
        # all batches and all tallies will have the same score ordering, hence
        # the 0's in the 1st and 3rd dimensions.
        if scoreType[0][i_tally][0][i_score] == 'flux': 
            fluxLoc[i_tally] = i_score

# Set printxs array according to the printXSflag input
if printxsFlag:
    for i_tally in range(len(fluxLoc)):
        if fluxLoc[i_tally] != -1:
            printxs[i_tally] = True
        
# Now rearrange the data as suitable, and perform xs conversion if necessary
for i_batch in range(len(mean)):
    for i_tally in range(len(mean[i_batch])):
        for i_filter in range(len(mean[i_batch][i_tally])):
            for i_score in range(len(mean[i_batch][i_tally][i_filter])):
                if (printxs[i_tally] and \
                    ((scoreType[0][i_tally][i_filter][i_score] != 'flux') and \
                    (scoreType[0][i_tally][i_filter][i_score] != 'current'))):
                    
                    # Perform rate to xs conversion
                    # mean is mean/fluxmean
                    meanPlot[i_tally][i_filter][i_score][i_batch] = \
                        mean[i_batch][i_tally][i_filter][i_score] / \
                        mean[i_batch][i_tally][i_filter][fluxLoc[i_tally]]
                    
                    # Update the relative uncertainty via error propagation
                    uncertPlot[i_tally][i_filter][i_score][i_batch] = \
                        sqrt(pow(uncert[i_batch][i_tally][i_filter][i_score],2.0) \
                        + pow(uncert[i_batch][i_tally][i_filter][fluxLoc[i_tally]],2.0))
                else: 
                    
                    # Do not perform rate to xs conversion
                    meanPlot[i_tally][i_filter][i_score][i_batch] = \
                        mean[i_batch][i_tally][i_filter][i_score]
                    uncertPlot[i_tally][i_filter][i_score][i_batch] = \
                        uncert[i_batch][i_tally][i_filter][i_score]
                
                # Both have the same absolute uncertainty calculation
                absUncertPlot[i_tally][i_filter][i_score][i_batch] = \
                    uncert[i_batch][i_tally][i_filter][i_score] * \
                    mean[i_batch][i_tally][i_filter][i_score]    

# Set plotting constants
xLabel = "Batches"
xLabel = xLabel.title() # not necessary for now, but is left in to handle if 
# the previous line changes

# Begin plotting
for i_tally in range(len(meanPlot)):
    # Set tally string (placeholder until I put tally labels in statePoint)
    tallyStr = "Tally " + str(i_tally + 1)
    
    for i_filter in range(len(meanPlot[i_tally])):

        # Set filter string
        filterStr = "Filter " + str(i_filter + 1)
                
        for i_score in range(len(meanPlot[i_tally][i_filter])):

            # Set score string
            scoreStr = scoreType[i_batch][i_tally][i_filter][i_score]
            scoreStr = scoreStr.title()
            if (printxs[i_tally] and ((scoreStr != 'Flux') and (scoreStr != 'Current'))):
                scoreStr = scoreStr + "-XS"
            
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
                meanPlot[i_tally][i_filter][i_score][:], \
                absUncertPlot[i_tally][i_filter][i_score][:])
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
                uncertPlot[i_tally][i_filter][i_score][:])
            plt.xlabel(xLabel)
            plt.ylabel("Relative Error of " + yLabel)
            plt.title("Relative Error of " + title)
            if (fileType != 'none'):
                plt.savefig(REfileName)
            if showImg:            
                plt.show()
            plt.clf()
            