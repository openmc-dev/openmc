#!/usr/bin/env python

import os
from subprocess import Popen, STDOUT, PIPE, call
import glob
from optparse import OptionParser

threshold = 0.2
"""Max relative difference in timings that triggers test failure"""

zero_cut = 10
"""If the reference time is zero, this test time in seconds will trigger failure"""

parser = OptionParser()
parser.add_option('--mpi_exec', dest='mpi_exec', default='')
parser.add_option('--mpi_np', dest='mpi_np', default='3')
parser.add_option('--exe', dest='exe')
(opts, args) = parser.parse_args()
cwd = os.getcwd()
logfile = os.path.join(cwd,'stdout.log')

def test_run():
    if opts.mpi_exec != '':
        proc = Popen([opts.mpi_exec, '-np', opts.mpi_np, opts.exe, cwd],
               stderr=STDOUT, stdout=PIPE)
    else:
        proc = Popen([opts.exe, cwd], stderr=STDOUT, stdout=PIPE)

    out,err = proc.communicate()
    print(out)
    with open(logfile,'w') as fh:
      fh.write(out)

    returncode = proc.returncode
    assert returncode == 0, 'OpenMC did not exit successfully.'

def test_results():
    call(['python', 'results.py', logfile])

    with open('results_true.dat') as fh: reftimes = fh.readlines()
    with open('results_test.dat') as fh: testtimes = fh.readlines()

    compare = True
    for ref,test in zip(reftimes,testtimes):
      r = float(ref.split()[1])
      t = float(test.split()[1])
      
      if r > 0.0:
        reldiff = abs(t - r)/r
        if reldiff > threshold:
          compare = False
      elif r < -1e-9 or t < -1e-9:
        raise Exception('Negative time!')
      else:
        if t > zero_cut:
          compare = False

    if not compare:
      os.rename('results_test.dat', 'results_error.dat')
    assert compare, 'Results do not agree.'

def teardown():
    output = glob.glob(os.path.join(cwd,'statepoint.*.*'))
    output.append(logfile)
    output.append(os.path.join(cwd,'results_test.dat'))
    for f in output:
        if os.path.exists(f):
            os.remove(f)

if __name__ == '__main__':

    # test for openmc executable
    if opts.exe is None:
        raise Exception('Must specify OpenMC executable from command line with --exe.')

    # run tests
    test_run()
    test_results()
    teardown()
