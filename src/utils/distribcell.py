#!/usr/bin/env python2

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

usehdf5 = False
##################################### END USER OPTIONS

## Find if tallies.xml exists.
#if glob('./tallies.xml') != None:
#    # It exists
#    tallyData = talliesXML('tallies.xml')
#else: 
#    # It does not exist.
#    tallyData = None

# Find all statepoints in this directory.
if usehdf5:
  ext = '.h5'
else:
  ext = '.binary'
files = glob('./statepoint.*'+ext)
fileNums = []
begin = 13
# Arrange the file list in increasing batch order
for i in range(len(files)):
    end = files[i].find(ext)
    fileNums.append(int(files[i][begin:end]))
fileNums.sort()
# Re-make filenames
files = []
for i in range(len(fileNums)):
    files.append("./statepoint." + str(fileNums[i]) + ext)

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

    #sp.geom._print_all()

    sp.read_results()
    
    # Read the same value for our cell and distribcell tally
    # This distribcell tally only has one entry
    print sp.get_value(0,[('cell',0)],0)[0]
    print sp.get_value(1,[('distribcell',0)],0)[0]
    print ''

    # Check that the sum of the distribcell tally 
    # is equal to that of the cell tally single entry
    print sp.get_value(2,[('cell',0)],0)[0]
    sum = 0
    for i in range(241):
      #print i
      sum += sp.get_value(3,[('distribcell',i)],0)[0]
    print sum
    print ''


    # What if we don't know the index in the array?
    # Well, we know how to get there. Explained:
    # Universe 0 -> Cell 1 -> Lattice 200 (x,y,z = 4,6,1) -> Universe 6 -> Cell 60
    path = [0,1,(200,4,6,1),6,60]

    # So we have the path to our cell now.
    # Lets get that specific value
    print sp.get_value(3,[('distribcell',path)],0)

