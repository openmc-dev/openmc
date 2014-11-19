#!/usr/bin/env python

import os
from subprocess import Popen, STDOUT, PIPE, call
import filecmp
import glob
from optparse import OptionParser

parser = OptionParser()
parser.add_option('--mpi_exec', dest='mpi_exec', default='')
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

def test_output_exists():
    assert os.path.exists('material.m1.binary') or os.path.exists('material.m1.h5'), 'Materials file does not exist.'
    assert os.path.exists('material.m2.binary') or os.path.exists('material.m2.h5'), 'Materials file does not exist.'
    assert os.path.exists('material.m3.binary') or os.path.exists('material.m3.h5'), 'Materials file does not exist.'

def test_results():
    call(['python', 'results.py'])
    compare = filecmp.cmp('results_test.dat', 'results_true.dat')
    if not compare:
      os.rename('results_test.dat', 'results_error.dat')
    assert compare, 'Results do not agree.'

def teardown():
    output = glob.glob(os.path.join(cwd, 'statepoint.20.*'))
    output += glob.glob(os.path.join(cwd, 'material.m*'))
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
        test_output_exists()
        test_results()
    finally:
        teardown()
