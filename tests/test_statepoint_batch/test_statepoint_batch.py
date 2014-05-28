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

def test_statepoints_exist():
    statepoint = glob.glob(os.path.join(cwd, 'statepoint.3.*'))
    assert len(statepoint) == 1, 'Either multiple or no statepoint.3 files exist.'
    assert statepoint[0].endswith('binary') or statepoint[0].endswith('h5'),\
        'Statepoint.3 file is not a binary or hdf5 file.'
    statepoint = glob.glob(os.path.join(cwd, 'statepoint.6.*'))
    assert len(statepoint) == 1, 'Either multiple or no statepoint.6 files exist.'
    assert statepoint[0].endswith('binary') or statepoint[0].endswith('h5'),\
        'Statepoint.6 file is not a binary or hdf5 file.'
    statepoint = glob.glob(os.path.join(cwd, 'statepoint.9.*'))
    assert len(statepoint) == 1, 'Either multiple or no statepoint.9 files exist.'
    assert statepoint[0].endswith('binary') or statepoint[0].endswith('h5'),\
        'Statepoint.9 file is not a binary or hdf5 file.'

def test_results():
    statepoint = glob.glob(os.path.join(cwd, 'statepoint.9.*'))
    call(['python', 'results.py', statepoint[0]])
    compare = filecmp.cmp('results_test.dat', 'results_true.dat')
    if not compare:
      os.rename('results_test.dat', 'results_error.dat')
    assert compare, 'Results do not agree.'

def teardown():
    output = glob.glob(os.path.join(cwd, 'statepoint*'))
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
        test_statepoints_exist()
        test_results()
    finally:
        teardown()
