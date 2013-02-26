#!/usr/bin/python
# Filename: eigenfunction_rms.py

# import packages
import statepoint
import numpy as np
import os

def main():

  # read in statepoint header data
  sp = statepoint.StatePoint('statepoint.ref.binary')

  # read in results
  sp.read_results()

  # extract results
  results = sp.extract_results(2,'nu-fission')

  # extract means and copy
  mean = results['mean'].copy()

  # reshape and integrate over energy
  mean = mean.reshape(results['bin_max'],order='F')
  mean = np.sum(mean,0)
  mean = np.sum(mean,0)
  mean = mean/mean.sum()*(mean > 1.e-8).sum()

  # write gnuplot file
  write_gnuplot('testsrc_pin','Pin Mesh',mean,np.size(mean,0))

  # extract results
  results = sp.extract_results(5,'nu-fission')

  # extract means and copy
  mean = results['mean'].copy()

  # reshape and integrate over energy
  mean = mean.reshape(results['bin_max'],order='F')
  mean = np.sum(mean,0)
  mean = np.sum(mean,0)
  mean = mean/mean.sum()*(mean > 1.e-8).sum()

  # write gnuplot file
  write_gnuplot('testsrc_quart','Quarter Assembly Mesh',mean,np.size(mean,0))

  # extract results
  results = sp.extract_results(8,'nu-fission')

  # extract means and copy
  mean = results['mean'].copy()

  # reshape and integrate over energy
  mean = mean.reshape(results['bin_max'],order='F')
  mean = np.sum(mean,0)
  mean = np.sum(mean,0)
  mean = mean/mean.sum()*(mean > 1.e-8).sum()

  # write gnuplot file
  write_gnuplot('testsrc_assy','Assembly Mesh',mean,np.size(mean,0))

def write_gnuplot(path,name,src,size):

  # Header String for GNUPLOT
  headerstr = """#!/usr/bin/env gnuplot

set terminal pdf enhanced
set output '{output}'
set palette defined (0 '#000090', 1 '#000fff', 2 '#0090ff', 3 '#0fffee', 4 '#90ff70', 5 '#ffee00', 6 '#ff7000', 7 '#ee0000', 8 '#7f0000')
set view map
set size ratio -1
set lmargin at screen 0.10
set rmargin at screen 0.90
set bmargin at screen 0.15
set tmargin at screen 0.90
unset xtics
unset ytics
set title  '{title}'""".format(output=path+'.pdf',title=name)

  # Write out the plot string
  pltstr = "splot '-' matrix with image "

  # Write out the data string
  i = 0
  datastr = ''
  while i < size:
    j = 0
    while j < size:
      datastr = datastr + '{0} '.format(src[i,j][0])
      j += 1
    datastr = datastr + '\n'
    i += 1

  # replace all nan with zero
  datastr = datastr.replace('nan','0.0')

  # Concatenate all
  outstr = headerstr + '\n' + pltstr + '\n' + datastr

  # Write File
  with open(path+".plot",'w') as f:
    f.write(outstr)

  # Run GNUPLOT
  os.system('gnuplot ' + path+".plot")

if __name__ == "__main__":
  main()
