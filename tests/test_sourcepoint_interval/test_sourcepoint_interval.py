#!/usr/bin/env python

import os
from subprocess import Popen, STDOUT, PIPE, call
import filecmp
import glob
from optparse import OptionParser

parser = OptionParser()
parser.add_option('--mpi_exec', dest='mpi_exec', default='')
parser.add_option('--mpi_np', dest='mpi_np', default='3')
parser.add_option('--exe', dest='exe')
(opts, args) = parser.parse_args()
cwd = os.getcwd()

def test_run():
    if opts.mpi_exec != '':
        proc = Popen([opts.mpi_exec, '-np', opts.mpi_np, opts.exe, cwd],
               stderr=STDOUT, stdout=PIPE)
    else:
        proc = Popen([opts.exe, cwd], stderr=STDOUT, stdout=PIPE)
    print(proc.communicate()[0])
    returncode = proc.returncode
    assert returncode == 0, 'OpenMC did not exit successfully.'

def test_statepoint_exists():
    statepoint = glob.glob(os.path.join(cwd, 'statepoint.*'))
    assert len(statepoint) == 5, '5 statepoint files must exist.'
    assert statepoint[0].endswith('binary') or statepoint[0].endswith('h5'),\
       'Statepoint file detected that is not binary or hdf5.'

def test_results():
    statepoint = glob.glob(os.path.join(cwd, 'statepoint.8.*'))
    call(['python', 'results.py', statepoint[0]])
    compare = filecmp.cmp('results_test.dat', 'results_true.dat')
    if not compare:
      os.rename('results_test.dat', 'results_error.dat')
    assert compare, 'Results do not agree.'

def teardown():
    output = glob.glob(os.path.join(cwd, 'statepoint.*'))
    output.append(os.path.join(cwd, 'results_test.dat'))
    for f in output:
        if os.path.exists(f):
            os.remove(f)

if __name__ == '__main__':

    # test for openmc executable
    if opts.exe is None:
        raise Exception('Must specify OpenMC executable from command line with --exe.')

    # run tests
    try:
        test_run()
        test_statepoint_exists()
        test_results()
    finally:
        teardown()
