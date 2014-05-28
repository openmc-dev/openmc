#!/usr/bin/env python

import os
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
cwd = os.getcwd()

def test_run():
    openmc_path = cwd + '/../../src/openmc'
    # Call 1
    copyfile("./materials1.xml","./materials.xml")
    copyfile("./geometry1.xml","./geometry.xml")
    copyfile("./settings1.xml","./settings.xml")
    copyfile("./tallies1.xml","./tallies.xml")
    if opts.mpi_exec != '':
        proc = Popen([opts.mpi_exec, '-np', opts.mpi_np, opts.exe, cwd],
               stderr=STDOUT, stdout=PIPE)
    else:
        proc = Popen([opts.exe, cwd], stderr=STDOUT, stdout=PIPE)
    print(proc.communicate()[0])
    returncode = proc.returncode
    assert returncode == 0
    if os.path.isfile("./statepoint.1.binary"):
        os.rename("./statepoint.1.binary","./statepoint.10.binary")
    else:
        os.rename("./statepoint.1.h5","./statepoint.10.h5")
    # Call 2
    copyfile("./materials2.xml","./materials.xml")
    copyfile("./geometry2.xml","./geometry.xml")
    copyfile("./settings2.xml","./settings.xml")
    copyfile("./tallies2.xml","./tallies.xml")
    if opts.mpi_exec != '':
        proc = Popen([opts.mpi_exec, '-np', opts.mpi_np, opts.exe, cwd],
               stderr=STDOUT, stdout=PIPE)
    else:
        proc = Popen([opts.exe, cwd], stderr=STDOUT, stdout=PIPE)
    print(proc.communicate()[0])
    returncode = proc.returncode
    assert returncode == 0
    if os.path.isfile("./statepoint.1.binary"):
        os.rename("./statepoint.1.binary","./statepoint.2.binary")
        os.rename("./statepoint.10.binary","./statepoint.1.binary")
    else:
        os.rename("./statepoint.1.h5","./statepoint.2.h5")
        os.rename("./statepoint.10.h5","./statepoint.1.h5")

    # Call 3
    copyfile("./materials3.xml","./materials.xml")
    copyfile("./geometry3.xml","./geometry.xml")
    copyfile("./settings3.xml","./settings.xml")
    copyfile("./tallies3.xml","./tallies.xml")
    if opts.mpi_exec != '':
        proc = Popen([opts.mpi_exec, '-np', opts.mpi_np, opts.exe, cwd],
               stderr=STDOUT, stdout=PIPE)
    else:
        proc = Popen([opts.exe, cwd], stderr=STDOUT, stdout=PIPE)
    print(proc.communicate()[0])
    returncode = proc.returncode
    assert returncode == 0

def test_created_statepoint():
    statepoint1 = glob.glob(cwd + '/statepoint.1.*')
    statepoint2 = glob.glob(cwd + '/statepoint.2.*')
    statepoint3 = glob.glob(cwd + '/statepoint.3.*')
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
    statepoint = list()
    statepoint.append(glob.glob(cwd + '/statepoint.1.*'))
    statepoint.append(glob.glob(cwd + '/statepoint.2.*'))
    statepoint.append(glob.glob(cwd + '/statepoint.3.*'))
    call(['python', 'results.py', statepoint.pop()[0], statepoint.pop()[0], statepoint.pop()[0]])
    compare = filecmp.cmp('results_test.dat', 'results_true.dat')
    if not compare:
      os.rename('results_test.dat', 'results_error.dat')
    assert compare

def teardown():
    #output = glob.glob(cwd + '/statepoint.*')
    output.append(cwd + '/geometry.xml')
    output.append(cwd + '/materials.xml')
    output.append(cwd + '/settings.xml')
    output.append(cwd + '/tallies.xml')
    output.append(cwd + '/results_test.dat')
    print output
    for f in output:
        print f
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
