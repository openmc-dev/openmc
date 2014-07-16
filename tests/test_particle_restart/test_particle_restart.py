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

def test_created_restart():
    particle = glob.glob(os.path.join(cwd, 'particle_12_192.*'))
    assert len(particle) == 1, 'Either multiple or no particle restart files exist.'
    assert particle[0].endswith('binary') or \
           particle[0].endswith('h5'), 'Particle restart file not a binary or hdf5 file.'

def test_results():
    particle = glob.glob(os.path.join(cwd, 'particle_12_192.*'))
    call(['python', 'results.py', particle[0]])
    compare = filecmp.cmp('results_test.dat', 'results_true.dat')
    if not compare:
      os.rename('results_test.dat', 'results_error.dat')
    assert compare, 'Results do not agree.'

def test_run_restart():
    particle = glob.glob(os.path.join(cwd, 'particle_12_192.*'))
    proc = Popen([opts.exe, '-r', particle[0], cwd], stderr=STDOUT, stdout=PIPE)
    print(proc.communicate()[0])
    returncode = proc.returncode 
    assert returncode == 0, 'Particle restart not successful.'

def teardown():
    output = glob.glob(os.path.join(cwd, 'statepoint.*')) + \
             glob.glob(os.path.join(cwd, 'particle_*')) + \
             [os.path.join(cwd, 'results_test.dat')]
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
        test_created_restart()
        test_results()
        test_run_restart()
    finally:
        teardown()
