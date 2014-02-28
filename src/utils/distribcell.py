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

usehdf5 = True
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


path1 = [0,1,(5,1,1,1),2,201]
path2 = [0,1,(5,2,1,1),2,201]
path3 = [0,1,(5,1,2,1),2,201]
path4 = [0,1,(5,2,2,1),2,201]

print "Found the following files:",files

for i_batch in range(len(files)):
    
    # Get filename    
    batch_filename = files[i_batch]
    
    # Create StatePoint object
    sp = StatePoint(batch_filename)
    # sp.geom._print_all()
#    sp.geom._get_offset(path1,1)
#    sp.geom._get_offset(path2,1)
#    sp.geom._get_offset(path3,1)
#    sp.geom._get_offset(path4,1)

    sp.read_results()

#    print sp.get_value(0,[('distribcell',path1)],0)
#    print sp.get_value(0,[('distribcell',path2)],0)
#    print sp.get_value(0,[('distribcell',path3)],0)
#    print sp.get_value(0,[('distribcell',path4)],0)
    
    print sp.get_value(0,[('cell',0)],0)
    print sp.get_value(1,[('distribcell',0)],0)
    print sp.get_value(2,[('cell',0)],0)
    sum = [0,0]
    for i in range(241):
      #print i
      sum[0] += sp.get_value(3,[('distribcell',i)],0)[0]
      sum[1] += sp.get_value(3,[('distribcell',i)],0)[1]
    print sum
    print sp.get_value(4,[('cell',0)],0)
    sum = [0,0]
    for i in range(63624):
      #print i
      sum[0] += sp.get_value(5,[('distribcell',i)],0)[0]
      sum[1] += sp.get_value(5,[('distribcell',i)],0)[1]
    print sum
    print sp.get_value(6,[('cell',0)],0)
    sum = [0,0]
    for i in range(241):
      #print i
      sum[0] += sp.get_value(7,[('distribcell',i)],0)[0]
      sum[1] += sp.get_value(7,[('distribcell',i)],0)[1]
    print sum
    print sp.get_value(8,[('cell',0)],0)
    sum = [0,0]
    for i in range(63624):
      #print i
      sum[0] += sp.get_value(9,[('distribcell',i)],0)[0]
      sum[1] += sp.get_value(9,[('distribcell',i)],0)[1]
    print sum
