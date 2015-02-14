#!/usr/bin/env python

# This program takes OpenMC statepoint binary files and creates a variety of
# outputs from them which should provide the user with an idea of the
# convergence behavior of all the tallies and filters defined by the user in
# tallies.xml.  The program can directly plot the value and errors of each
# tally, filter, score combination; it can save these plots to a file; and 
# it can also save the data used in these plots to a CSV file for importing in
# to other plotting packages such as Excel, gnuplot, MathGL, or Veusz.

# To use the program, run this program from the working directory of the openMC
# problem to analyze.  

# The USER OPTIONS block below provides four options for the user to set:
# fileType, printxs, showImg, and savetoCSV.  See the options block for more 
# information.

from __future__ import print_function
from math import sqrt, pow
from glob import glob

import numpy as np
import scipy.stats
import matplotlib.pyplot as plt

from statepoint import StatePoint

##################################### USER OPTIONS

# Set filetype (the file extension desired, without the period.)
# Options are backend dependent, but most backends support png, pdf, ps, eps 
# and svg.  Write "none" if no saved files are desired.
fileType = "none"

# Set if cross-sections or reaction rates are desired printxs = True means X/S
printxs = False

# Set if the figures should be displayed to screen or not (True means show)
showImg = False

# Save to CSV for use in more advanced plotting programs like GNUPlot, MathGL
savetoCSV = True

##################################### END USER OPTIONS

## Find if tallies.xml exists.
#if glob('./tallies.xml') != None:
#    # It exists
#    tallyData = talliesXML('tallies.xml')
#else: 
#    # It does not exist.
#    tallyData = None

# Find all statepoints in this directory.
files = glob('./statepoint.*.binary')
fileNums = []
begin = 13
# Arrange the file list in increasing batch order
for i in range(len(files)):
    end = files[i].find(".binary")
    fileNums.append(int(files[i][begin:end]))
fileNums.sort()
# Re-make filenames
files = []
for i in range(len(fileNums)):
    files.append("./statepoint." + str(fileNums[i]) + ".binary")

# Initialize arrays as needed
mean = [None for x in range(len(files))]
uncert = [None for x in range(len(files))]
scoreType = [None for x in range(len(files))]
active_batches = [None for x in range(len(files))]

for i_batch in range(len(files)):
    
    # Get filename    
    batch_filename = files[i_batch]
    
    # Create StatePoint object
    sp = StatePoint(batch_filename)

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
    
    # Increase the dimensionality of our main variables
    mean[i_batch] = [None for x in range(len(sp.tallies))]
    uncert[i_batch] = [None for x in range(len(sp.tallies))]
    scoreType[i_batch] = [None for x in range(len(sp.tallies))]   
    
    # Loop over all tallies
    for i_tally, t in enumerate(sp.tallies):
        # Calculate t-value for 95% two-sided CI
        n = t.n_realizations
        t_value = scipy.stats.t.ppf(0.975, n - 1)
    
        # Store the batch count    
        active_batches[i_batch] = n
    
        # Resize the 2nd dimension      
        mean[i_batch][i_tally] = [None for x in range(t.total_filter_bins)]
        uncert[i_batch][i_tally] = [None for x in range(t.total_filter_bins)]
        scoreType[i_batch][i_tally] = [None for x in range(t.total_filter_bins)]
        
        for i_filter in range(t.total_filter_bins):
            # Resize the 3rd dimension            
            mean[i_batch][i_tally][i_filter] = [None for x in range(t.n_nuclides)]
            uncert[i_batch][i_tally][i_filter] = [None for x in range(t.n_nuclides)]
            scoreType[i_batch][i_tally][i_filter] = [None for x in range(t.n_nuclides)]
            print(t.total_filter_bins,t.n_nuclides)            
            for i_nuclide in range(t.n_nuclides):
                mean[i_batch][i_tally][i_filter][i_nuclide] = \
                    [None for x in range(t.n_scores)]
                uncert[i_batch][i_tally][i_filter][i_nuclide] = \
                    [None for x in range(t.n_scores)]
                scoreType[i_batch][i_tally][i_filter][i_nuclide] = \
                    [None for x in range(t.n_scores)]
                i_score = 0
                while i_score < len(t.scores):
                    if ((t.scores[i_score] == 'scatter-pn') or \
                        (t.scores[i_score] == 'ndpp-scatter-pn')):
                        order = t.scatt_order[i_score]
                        for i_s in xrange(order + 1):  
                            scoreType[i_batch][i_tally][i_filter][i_nuclide][i_score] = \
                                t.scores[i_score][:-1]+str(i_s)
                            s, s2 = sp._get_double(2)
                            s /= n
                            mean[i_batch][i_tally][i_filter][i_nuclide][i_score] = s
                            if s != 0.0:
                                relative_error = t_value*sqrt((s2/n - s*s)/(n-1))/s
                            else:
                                relative_error = 0.0
                            uncert[i_batch][i_tally][i_filter][i_nuclide][i_score] = relative_error
                            i_score = i_score + 1
                    elif  (t.scores[i_score] == 'scatter-n'):
                        order = t.scatt_order[i_score]
                        scoreType[i_batch][i_tally][i_filter][i_nuclide][i_score] = \
                            t.scores[i_score][:-1]+str(order)
                        s, s2 = sp._get_double(2)
                        s /= n
                        mean[i_batch][i_tally][i_filter][i_nuclide][i_score] = s
                        if s != 0.0:
                            relative_error = t_value*sqrt((s2/n - s*s)/(n-1))/s
                        else:
                            relative_error = 0.0
                        uncert[i_batch][i_tally][i_filter][i_nuclide][i_score] = relative_error
                        i_score = i_score + 1
                    else:
                        scoreType[i_batch][i_tally][i_filter][i_nuclide][i_score] = \
                            t.scores[i_score]
                        s, s2 = sp._get_double(2)
                        s /= n
                        mean[i_batch][i_tally][i_filter][i_nuclide][i_score] = s
                        if s != 0.0:
                            relative_error = t_value*sqrt((s2/n - s*s)/(n-1))/s
                        else:
                            relative_error = 0.0
                        uncert[i_batch][i_tally][i_filter][i_nuclide][i_score] = relative_error
                        i_score = i_score + 1

# Reorder the data lists in to a list order more conducive for plotting:
# The indexing should be: [tally][filter][score][batch]
meanPlot = [None for x in range(len(mean[0]))] # Set to the number of tallies
uncertPlot = [None for x in range(len(mean[0]))] # Set to the number of tallies
absUncertPlot = [None for x in range(len(mean[0]))] # Set to number of tallies
filterLabel = [None for x in range(len(mean[0]))] #Set to the number of tallies
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
            
        for i_nuclide in range(len(meanPlot[i_tally][i_filter])):
            # Set 4th (nuclide)) dimension
            meanPlot[i_tally][i_filter][i_nuclide] = \
                [None for x in range(len(mean[0][i_tally][i_filter][i_nuclide]))]
            uncertPlot[i_tally][i_filter][i_nuclide] = \
                [None for x in range(len(mean[0][i_tally][i_filter][i_nuclide]))]
            absUncertPlot[i_tally][i_filter][i_nuclide] = \
                [None for x in range(len(mean[0][i_tally][i_filter][i_nuclide]))]
        
            for i_score in range(len(meanPlot[i_tally][i_filter][i_nuclide])):
                # Set 5th (batch) dimension
                meanPlot[i_tally][i_filter][i_nuclide][i_score] = \
                    [None for x in range(len(mean))]
                uncertPlot[i_tally][i_filter][i_nuclide][i_score] = \
                    [None for x in range(len(mean))]
                absUncertPlot[i_tally][i_filter][i_nuclide][i_score] = \
                    [None for x in range(len(mean))]
                    
                # Get filterLabel (this should be moved to its own function)
                #??? How to do?
            
        # Set flux location if found
        # all batches and all tallies will have the same score ordering, hence
        # the 0's in the 1st, 3rd, and 4th dimensions.
        if scoreType[0][i_tally][0][0][i_score] == 'flux': 
            fluxLoc[i_tally] = i_score

# Set printxs array according to the printxs input
if printxs:
    for i_tally in range(len(fluxLoc)):
        if fluxLoc[i_tally] != -1:
            printxs[i_tally] = True
        
# Now rearrange the data as suitable, and perform xs conversion if necessary
for i_batch in range(len(mean)):
    for i_tally in range(len(mean[i_batch])):
        for i_filter in range(len(mean[i_batch][i_tally])):
            for i_nuclide in range(len(mean[i_batch][i_tally][i_filter])):
                for i_score in range(len(mean[i_batch][i_tally][i_filter][i_nuclide])):
                    if (printxs[i_tally] and \
                        ((scoreType[0][i_tally][i_filter][i_nuclide][i_score] != 'flux') and \
                        (scoreType[0][i_tally][i_filter][i_nuclide][i_score] != 'current'))):
                        
                        # Perform rate to xs conversion
                        # mean is mean/fluxmean
                        meanPlot[i_tally][i_filter][i_nuclide][i_score][i_batch] = \
                            mean[i_batch][i_tally][i_filter][i_nuclide][i_score] / \
                            mean[i_batch][i_tally][i_filter][i_nuclide][fluxLoc[i_tally]]
                        
                        # Update the relative uncertainty via error propagation
                        uncertPlot[i_tally][i_filter][i_nuclide][i_score][i_batch] = \
                            sqrt(pow(uncert[i_batch][i_tally][i_filter][i_nuclide][i_score],2) \
                            + pow(uncert[i_batch][i_tally][i_filter][i_nuclide][fluxLoc[i_tally]],2))
                    else: 
                        
                        # Do not perform rate to xs conversion
                        meanPlot[i_tally][i_filter][i_nuclide][i_score][i_batch] = \
                            mean[i_batch][i_tally][i_filter][i_nuclide][i_score]
                        uncertPlot[i_tally][i_filter][i_nuclide][i_score][i_batch] = \
                            uncert[i_batch][i_tally][i_filter][i_nuclide][i_score]

                    # Both have the same absolute uncertainty calculation
                    absUncertPlot[i_tally][i_filter][i_nuclide][i_score][i_batch] = \
                        uncert[i_batch][i_tally][i_filter][i_nuclide][i_score] * \
                        mean[i_batch][i_tally][i_filter][i_nuclide][i_score]    

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
        
        for i_nuclide in range(len(meanPlot[i_tally][i_filter])):
            
            nuclideStr = "Nuclide " + str(i_nuclide + 1)
        
            for i_score in range(len(meanPlot[i_tally][i_filter][i_nuclide])):
    
                # Set score string
                scoreStr = scoreType[i_batch][i_tally][i_filter][i_nuclide][i_score]
                scoreStr = scoreStr.title()
                if (printxs[i_tally] and ((scoreStr != 'Flux') and \
                    (scoreStr != 'Current'))):
                    scoreStr = scoreStr + "-XS"
                
                # set Title 
                title = "Convergence of " + scoreStr + " in " + tallyStr + " for "\
                    + filterStr + " and " + nuclideStr
                
                # set yLabel
                yLabel = scoreStr
                yLabel = yLabel.title()
    
                # Set saving filename
                fileName = "tally_" + str(i_tally + 1) + "_" + scoreStr + \
                    "_filter_" + str(i_filter+1) + "_nuclide_" + str(i_nuclide+1) \
                    + "." + fileType
                REfileName = "tally_" + str(i_tally + 1) + "_" + scoreStr + \
                    "_RE_filter_" + str(i_filter+1) + "_nuclide_" + str(i_nuclide+1) \
                    + "." + fileType
                
                # Plot mean with absolute error bars
                plt.errorbar(active_batches, \
                    meanPlot[i_tally][i_filter][i_nuclide][i_score][:], \
                    absUncertPlot[i_tally][i_filter][i_nuclide][i_score][:],fmt='o-',aa=True)
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
                    uncertPlot[i_tally][i_filter][i_nuclide][i_score][:],'o-',aa=True)
                plt.xlabel(xLabel)
                plt.ylabel("Relative Error of " + yLabel)
                plt.title("Relative Error of " + title)
                if (fileType != 'none'):
                    plt.savefig(REfileName)
                if showImg:            
                    plt.show()
                plt.clf()
            
if savetoCSV:
    # This block loops through each tally, and for each tally:
        # Creates a new file
        # Writes the scores and filters for that tally in csv format.
        # The columns will be: batches,then for each filter: all the scores
        # The rows, of course, are the data points per batch.
    
    for i_tally in range(len(meanPlot)):
        # Set tally string (placeholder until I put tally labels in statePoint)
        tallyStr = "Tally " + str(i_tally + 1)
        CSV_filename = "./tally" + str(i_tally+1)+".csv"
        # Open the file
        f = open(CSV_filename, 'w')  
        
        # Write the header line        
        
        lineText = "Batches" 
        
        for i_filter in range(len(meanPlot[i_tally])):
        
            # Set filter string
            filterStr = "Filter " + str(i_filter + 1)
            
            for i_nuclide in range(len(meanPlot[i_tally][i_filter])):
            
                nuclideStr = "Nuclide " + str(i_nuclide + 1)
                    
                for i_score in range(len(meanPlot[i_tally][i_filter][i_nuclide])):
                    
                    # Set the title
                    scoreStr = scoreType[i_batch][i_tally][i_filter][i_nuclide][i_score]
                    scoreStr = scoreStr.title()
                    if (printxs[i_tally] and ((scoreStr != 'Flux') and \
                        (scoreStr != 'Current'))):
                        scoreStr = scoreStr + "-XS"
                    
                    # set header 
                    headerText = scoreStr + " for " + filterStr + " for " + nuclideStr
                    
                    lineText = lineText + "," + headerText + \
                        ",Abs Unc of " + headerText + \
                        ",Rel Unc of " + headerText
                    
        f.write(lineText + "\n")   

        # Write the data lines, each row is a different batch                     
        
        for i_batch in range(len(meanPlot[i_tally][0][0][0])):
        
            lineText = repr(active_batches[i_batch])        
        
            for i_filter in range(len(meanPlot[i_tally])):
                        
                for i_nuclide in range(len(meanPlot[i_tally][i_filter])):
                        
                    for i_score in range(len(meanPlot[i_tally][i_filter][i_nuclide])):
                        
                        fieldText = \
                            repr(meanPlot[i_tally][i_filter][i_nuclide][i_score][i_batch]) + \
                            "," + \
                            repr(absUncertPlot[i_tally][i_filter][i_nuclide][i_score][i_batch]) +\
                            "," + \
                            repr(uncertPlot[i_tally][i_filter][i_nuclide][i_score][i_batch])
                    
                        lineText = lineText + "," + fieldText
            
            f.write(lineText + "\n")
    
   
