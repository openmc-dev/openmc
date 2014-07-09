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

def test_created_outputs():
    outputs = [glob.glob(''.join((cwd, '/track_1_1_1.*')))]
    outputs.append(glob.glob(''.join((cwd, '/track_1_1_2.*'))))
    for files in outputs:
        assert len(files) == 1, 'Multiple or no track files detected.'
        assert files[0].endswith('binary') or files[0].endswith('h5'),\
        'Track files not a binary or hdf5 file'

def test_outputs():
    call(['python', 'results.py'])
    compare = filecmp.cmp('results_test.dat', 'results_true.dat')
    if not compare:
        os.rename('results_test.dat', 'results_error.dat')
    assert compare, 'Results to not agree'

def teardown():
    temp_files = glob.glob('statepoint*')
    temp_files += glob.glob('track*')
    temp_files += glob.glob('*.vtp')
    temp_files += glob.glob('*.pvtp')
    temp_files = temp_files + ['results_test.dat']
    for f in temp_files:
        if os.path.exists(f):
            os.remove(f)

if __name__ == '__main__':

    # test for openmc executable
    if opts.exe is None:
        raise Exception('Must specify OpenMC executable from command line with --exe.')

    # run tests
    try:
        test_run()
        test_created_outputs()
        test_outputs()
    finally:
        teardown()
