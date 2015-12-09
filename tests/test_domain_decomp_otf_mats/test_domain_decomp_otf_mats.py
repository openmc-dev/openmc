#!/usr/bin/env python

import os
from subprocess import Popen, STDOUT, PIPE, call
import filecmp
import glob
from optparse import OptionParser
import shutil

parser = OptionParser()
parser.add_option('--mpi_exec', dest='mpi_exec', default='')
parser.add_option('--mpi_np', dest='mpi_np', default='3')
parser.add_option('--exe', dest='exe')
(opts, args) = parser.parse_args()
cwd = os.getcwd()

def test_run():
    shutil.copy('materials.h5.keep','materials.h5')
    if opts.mpi_exec != '':
        proc = Popen([opts.mpi_exec, '-np', '4', opts.exe, cwd],
               stderr=STDOUT, stdout=PIPE)
        print(proc.communicate()[0])
        returncode = proc.returncode
        assert returncode == 0, 'OpenMC did not exit successfully.'
    else:
        # This test is only relevant if MPI is enabled
        return

def test_output_exists():
    if opts.mpi_exec == '':
        # This test is only relevant if MPI is enabled
        return
    assert os.path.exists('materials-out.h5'), 'Materials file does not exist.'

def test_results():
    if opts.mpi_exec == '':
        # This test is only relevant if MPI is enabled
        return
    call(['python', 'results.py'])
    compare = filecmp.cmp('results_test.dat', 'results_true.dat')
    if not compare:
      os.rename('results_test.dat', 'results_error.dat')
    assert compare, 'Results do not agree.'

def teardown():
    if opts.mpi_exec == '':
        # This test is only relevant if MPI is enabled
        return
    output = glob.glob(os.path.join(cwd, 'statepoint.20.*'))
    output += glob.glob(os.path.join(cwd, 'materials-out.h5'))
    output.append(os.path.join(cwd, 'results_test.dat'))
    for f in output:
        if os.path.exists(f):
            os.remove(f)

if __name__ == '__main__':

    # test for openmc executable
    if opts.exe is None:
        raise Exception('Must specify OpenMC executable from command line with --exe.')

    # test for mpi executable
    if opts.mpi_exec is '':
        raise Exception('Must specify MPI executable from command line with --mpi_exec.')

    # run tests
    try:
        test_run()
        test_output_exists()
        test_results()
    finally:
        teardown()
