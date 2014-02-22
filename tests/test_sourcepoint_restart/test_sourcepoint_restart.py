#!/usr/bin/env python

import os
from subprocess import Popen, STDOUT, PIPE, call
import filecmp
from nose_mpi import NoseMPI
import glob

pwd = os.path.dirname(__file__)

def setup(): 
    os.putenv('PWD', pwd)
    os.chdir(pwd)

def test_run():
    openmc_path = pwd + '/../../src/openmc'
    if int(NoseMPI.mpi_np) > 0:
        proc = Popen([NoseMPI.mpi_exec, '-np', NoseMPI.mpi_np, openmc_path],
               stderr=STDOUT, stdout=PIPE)
    else:
        proc = Popen([openmc_path], stderr=STDOUT, stdout=PIPE)
    print(proc.communicate()[0])
    returncode = proc.returncode
    assert returncode == 0

def test_created_statepoint():
    statepoint = glob.glob(pwd + '/statepoint.*')
    assert len(statepoint) == 2
    assert statepoint[0].endswith('binary') or statepoint[0].endswith('h5')
    sourcepoint = glob.glob(pwd + '/source.7.*')
    assert len(sourcepoint) == 1
    assert sourcepoint[0].endswith('binary') or sourcepoint[0].endswith('h5')

def test_results():
    statepoint = glob.glob(pwd + '/statepoint.10.*')
    call(['python', 'results.py', statepoint[0]])
    compare = filecmp.cmp('results_test.dat', 'results_true.dat')
    if not compare:
      os.rename('results_test.dat', 'results_error.dat')
    assert compare
    os.remove(statepoint[0])

def test_restart_form1():
    statepoint = glob.glob(pwd + '/statepoint.7.*')
    sourcepoint = glob.glob(pwd + '/source.7.*')
    openmc_path = pwd + '/../../src/openmc'
    if int(NoseMPI.mpi_np) > 0:
        proc = Popen([NoseMPI.mpi_exec, '-np', NoseMPI.mpi_np, openmc_path, 
                      '-r', statepoint[0], sourcepoint[0]],
               stderr=STDOUT, stdout=PIPE)
    else:
        proc = Popen([openmc_path, '-r', statepoint[0], sourcepoint[0]], stderr=STDOUT, stdout=PIPE)
    print(proc.communicate()[0])
    returncode = proc.returncode
    assert returncode == 0

def test_created_statepoint_form1():
    statepoint = glob.glob(pwd + '/statepoint.10.*')
    assert len(statepoint) == 1
    assert statepoint[0].endswith('binary') or statepoint[0].endswith('h5')

def test_results_form1():
    statepoint = glob.glob(pwd + '/statepoint.10.*')
    call(['python', 'results.py', statepoint[0]])
    compare = filecmp.cmp('results_test.dat', 'results_true.dat')
    if not compare:
      os.rename('results_test.dat', 'results_error.dat')
    assert compare
    os.remove(statepoint[0])

def test_restart_form2():
    statepoint = glob.glob(pwd + '/statepoint.7.*')
    sourcepoint = glob.glob(pwd + '/source.7.*')
    openmc_path = pwd + '/../../src/openmc'
    if int(NoseMPI.mpi_np) > 0:
        proc = Popen([NoseMPI.mpi_exec, '-np', NoseMPI.mpi_np, openmc_path, 
                      '--restart', statepoint[0], sourcepoint[0]],
               stderr=STDOUT, stdout=PIPE)
    else:
        proc = Popen([openmc_path, '--restart', statepoint[0], sourcepoint[0]],
                     stderr=STDOUT, stdout=PIPE)
    print(proc.communicate()[0])
    returncode = proc.returncode
    assert returncode == 0

def test_created_statepoint_form2():
    statepoint = glob.glob(pwd + '/statepoint.10.*')
    assert len(statepoint) == 1
    assert statepoint[0].endswith('binary') or statepoint[0].endswith('h5')

def test_results_form2():
    statepoint = glob.glob(pwd + '/statepoint.10.*')
    call(['python', 'results.py', statepoint[0]])
    compare = filecmp.cmp('results_test.dat', 'results_true.dat')
    if not compare:
      os.rename('results_test.dat', 'results_error.dat')
    assert compare
    os.remove(statepoint[0])

def test_restart_serial():
    statepoint = glob.glob(pwd + '/statepoint.7.*')
    sourcepoint = glob.glob(pwd + '/source.7.*')
    openmc_path = pwd + '/../../src/openmc'
    proc = Popen([openmc_path, '--restart', statepoint[0], sourcepoint[0]], 
                 stderr=STDOUT, stdout=PIPE)
    print(proc.communicate()[0])
    returncode = proc.returncode
    assert returncode == 0

def test_created_statepoint_serial():
    statepoint = glob.glob(pwd + '/statepoint.10.*')
    assert len(statepoint) == 1
    assert statepoint[0].endswith('binary') or statepoint[0].endswith('h5')

def test_results_serial():
    statepoint = glob.glob(pwd + '/statepoint.10.*')
    call(['python', 'results.py', statepoint[0]])
    compare = filecmp.cmp('results_test.dat', 'results_true.dat')
    if not compare:
      os.rename('results_test.dat', 'results_error.dat')
    assert compare

def teardown():
    output = glob.glob(pwd + '/statepoint.*')
    output += glob.glob(pwd + '/source.*')
    output.append(pwd + '/results_test.dat')
    for f in output:
        if os.path.exists(f):
            os.remove(f)
