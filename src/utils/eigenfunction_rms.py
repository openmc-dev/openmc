#!/usr/bin/python
# Filename: eigenfunction_rms.py

# import packages
import statepoint
import numpy as np
import os
import sys

def main(tally_id, score_id, batch_start, batch_end, name):

  # read in statepoint header data
  sp = statepoint.StatePoint('statepoint.ref.binary')

  # read in results
  sp.read_results()

  # extract reference mean
  mean_ref = extract_mean(sp, tally_id, score_id)

  # write gnuplot file
  write_src_gnuplot('testsrc_pin','Pin mesh',mean_ref,np.size(mean_ref,0))

  # preallocate arrays
  hists = np.zeros(batch_end - batch_start + 1)
  norms = np.zeros(batch_end - batch_start + 1)

  i = batch_start
  while i <= batch_end:

    # process statepoint
    sp = statepoint.StatePoint('statepoint.'+str(i)+'.binary')
    sp.read_results()

    # extract mean
    mean = extract_mean(sp, tally_id, score_id)

    # calculate L2 norm
    norm = np.linalg.norm(mean - mean_ref)

    # get history information
    n_inactive = sp.n_inactive
    current_batch = sp.current_batch
    n_particles = sp.n_particles
    gen_per_batch = sp.gen_per_batch
    n_histories = (current_batch - n_inactive)*n_particles*gen_per_batch

    # batch in vectors
    hists[i - batch_start] = n_histories
    norms[i - batch_start] = norm

    # print
    print 'Batch: '+str(i)+' Histories: '+str(n_histories)+' Norm: '+str(norm)

    i += 1

  # write out gnuplot file
  write_norm_gnuplot(name,hists,norms,np.size(hists))

def extract_mean(sp, tally_id,score_id):

  # extract results
  results = sp.extract_results(tally_id,score_id)

  # extract means and copy
  mean = results['mean'].copy()

  # reshape and integrate over energy
  mean = mean.reshape(results['bin_max'],order='F')
  mean = np.sum(mean,0)
  mean = np.sum(mean,0)
  mean = mean/mean.sum()*(mean > 1.e-8).sum()

  return mean

def write_norm_gnuplot(path,xdat,ydat,size):

  # Header String for GNUPLOT
  headerstr = """#!/usr/bin/env gnuplot

set terminal pdf enhanced
set output '{output}'
set ylabel 'L-2 norm'
set xlabel 'Histories'
set log x
set log y
""".format(output=path+'.pdf')

  # Write out the plot string
  pltstr = "plot '-' using 1:2 with lines"

  # Write out the data string
  i = 0
  datastr = ''
  while i < size:
    datastr = datastr + '{0} {1}\n'.format(xdat[i],ydat[i])
    i += 1

  # Concatenate all
  outstr = headerstr + '\n' + pltstr + '\n' + datastr

  # Write File
  with open(path+".plot",'w') as f:
    f.write(outstr)

  # Run GNUPLOT
  os.system('gnuplot ' + path+".plot")

def write_src_gnuplot(path,name,src,size):

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
set title '{title}'""".format(output=path+'.pdf',title=name)

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
  tally_id = int(sys.argv[1])
  score_id = sys.argv[2]
  batch_start = int(sys.argv[3])
  batch_end = int(sys.argv[4])
  name = sys.argv[5]
  main(tally_id, score_id, batch_start, batch_end, name)
