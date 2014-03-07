#!/usr/bin/env python

import os
from subprocess import Popen, STDOUT, PIPE, call
import filecmp
from nose_mpi import NoseMPI
import glob
from shutil import copyfile

pwd = os.path.dirname(__file__)

def setup(): 
    os.putenv('PWD', pwd)
    os.chdir(pwd)

def test_run():
    openmc_path = pwd + '/../../src/openmc'
    # Call 1
    copyfile("./materials1.xml","./materials.xml")
    copyfile("./geometry1.xml","./geometry.xml")
    copyfile("./settings1.xml","./settings.xml")
    copyfile("./tallies1.xml","./tallies.xml")
    if int(NoseMPI.mpi_np) > 0:
        proc = Popen([NoseMPI.mpi_exec, '-np', NoseMPI.mpi_np, openmc_path],
               stderr=STDOUT, stdout=PIPE)
    else:
        proc = Popen([openmc_path], stderr=STDOUT, stdout=PIPE)
    print(proc.communicate()[0])
    returncode = proc.returncode
    assert returncode == 0
    # Call 2
    copyfile("./materials2.xml","./materials.xml")
    copyfile("./geometry2.xml","./geometry.xml")
    copyfile("./settings2.xml","./settings.xml")
    copyfile("./tallies2.xml","./tallies.xml")
    if int(NoseMPI.mpi_np) > 0:
        proc = Popen([NoseMPI.mpi_exec, '-np', NoseMPI.mpi_np, openmc_path],
               stderr=STDOUT, stdout=PIPE)
    else:
        proc = Popen([openmc_path], stderr=STDOUT, stdout=PIPE)
    print(proc.communicate()[0])
    returncode = proc.returncode
    assert returncode == 0
    # Call 3
    copyfile("./materials3.xml","./materials.xml")
    copyfile("./geometry3.xml","./geometry.xml")
    copyfile("./settings3.xml","./settings.xml")
    copyfile("./tallies3.xml","./tallies.xml")
    if int(NoseMPI.mpi_np) > 0:
        proc = Popen([NoseMPI.mpi_exec, '-np', NoseMPI.mpi_np, openmc_path],
               stderr=STDOUT, stdout=PIPE)
    else:
        proc = Popen([openmc_path], stderr=STDOUT, stdout=PIPE)
    print(proc.communicate()[0])
    returncode = proc.returncode
    assert returncode == 0

def test_created_statepoint():
    statepoint = glob.glob(pwd + '/statepoint.1.*') + 
                 glob.glob(pwd + '/statepoint.2.*') + 
                 glob.glob(pwd + '/statepoint.3.*') 
    assert len(statepoint) == 3
    assert statepoint[0].endswith('binary') or statepoint[0].endswith('h5')
    assert statepoint[1].endswith('binary') or statepoint[1].endswith('h5')
    assert statepoint[2].endswith('binary') or statepoint[2].endswith('h5')

def test_results():
    statepoint = glob.glob(pwd + '/statepoint.1.*') + 
                 glob.glob(pwd + '/statepoint.2.*') + 
                 glob.glob(pwd + '/statepoint.3.*') 
    call(['python', 'results.py', statepoint)
    compare = filecmp.cmp('results_test.dat', 'results_true.dat')
    if not compare:
      os.rename('results_test.dat', 'results_error.dat')
    assert compare

def teardown():
    statepoint = glob.glob(pwd + '/statepoint.1.*') + 
                 glob.glob(pwd + '/statepoint.2.*') + 
                 glob.glob(pwd + '/statepoint.3.*') 
    output.append(pwd + '/results_test.dat')
    for f in output:
        if os.path.exists(f):
            os.remove(f)
