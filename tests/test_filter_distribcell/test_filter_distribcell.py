#!/usr/bin/env python

import os
import sys
from subprocess import Popen, STDOUT, PIPE, call
import filecmp
import glob
from optparse import OptionParser
from shutil import copyfile

parser = OptionParser()
parser.add_option('--mpi_exec', dest='mpi_exec', default='')
parser.add_option('--mpi_np', dest='mpi_np', default='3')
parser.add_option('--exe', dest='exe')
(opts, args) = parser.parse_args()

def test_run():

    if opts.mpi_exec != '':
        opts.mpi_exe = os.path.abspath(opts.mpi_exec)

    opts.exe = os.path.abspath(opts.exe)

    # Call 1
    os.chdir('case-1')
    cwd = os.getcwd()
    if opts.mpi_exec != '':
        proc = Popen([opts.mpi_exec, '-np', opts.mpi_np, opts.exe, cwd],
               stderr=STDOUT, stdout=PIPE)
    else:
        proc = Popen([opts.exe, cwd], stderr=STDOUT, stdout=PIPE)
    print(proc.communicate()[0])
    returncode = proc.returncode
    assert returncode == 0

    # Call 2
    os.chdir('../case-2')
    cwd = os.getcwd()
    if opts.mpi_exec != '':
        proc = Popen([opts.mpi_exec, '-np', opts.mpi_np, opts.exe, cwd],
               stderr=STDOUT, stdout=PIPE)
    else:
        proc = Popen([opts.exe, cwd], stderr=STDOUT, stdout=PIPE)
    print(proc.communicate()[0])
    returncode = proc.returncode
    assert returncode == 0

    # Call 3
    os.chdir('../case-3')
    cwd = os.getcwd()
    if opts.mpi_exec != '':
        proc = Popen([opts.mpi_exec, '-np', opts.mpi_np, opts.exe, cwd],
               stderr=STDOUT, stdout=PIPE)
    else:
        proc = Popen([opts.exe, cwd], stderr=STDOUT, stdout=PIPE)
    print(proc.communicate()[0])
    returncode = proc.returncode
    assert returncode == 0

    os.chdir('..')

def test_created_statepoint():
    cwd = os.getcwd()
    statepoint1 = glob.glob(cwd + '/case-1/statepoint.1.*')
    statepoint2 = glob.glob(cwd + '/case-2/statepoint.1.*')
    statepoint3 = glob.glob(cwd + '/case-3/statepoint.3.*')
    assert len(statepoint1) == 1
    assert len(statepoint2) == 1
    assert len(statepoint3) == 1
    string1 = statepoint1.pop()
    string2 = statepoint2.pop()
    string3 = statepoint3.pop()
    assert string1.endswith('binary') or string1.endswith('h5')
    assert string2.endswith('binary') or string2.endswith('h5')
    assert string3.endswith('binary') or string3.endswith('h5')

def test_results():
    cwd = os.getcwd()
    statepoint = list()
    statepoint.append(glob.glob(cwd + '/case-1/statepoint.1.*'))
    statepoint.append(glob.glob(cwd + '/case-2/statepoint.1.*'))
    statepoint.append(glob.glob(cwd + '/case-3/statepoint.3.*'))
    call([sys.executable, 'results.py', statepoint.pop()[0], statepoint.pop()[0], statepoint.pop()[0]])
    compare = filecmp.cmp('results_test.dat', 'results_true.dat')
    if not compare:
      os.rename('results_test.dat', 'results_error.dat')
    assert compare

def teardown():
    cwd = os.getcwd()
    output = glob.glob(cwd + '/statepoint.*')
    output.append(cwd + '/results_test.dat')
    for f in output:
        if os.path.exists(str(f)):
            os.remove(str(f))

if __name__ == '__main__':

    # test for openmc executable
    if opts.exe is None:
        raise Exception('Must specify OpenMC executable from command line with --exe.')

    # run tests
    try:
        test_run()
        test_created_statepoint()
        test_results()
    finally:
        teardown()
