#!/usr/bin/env python

import os
from subprocess import Popen, STDOUT, PIPE
import filecmp
from nose_mpi import NoseMPI

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
    returncode = proc.wait()
    print(proc.communicate()[0])
    assert returncode == 0

def test_created_statepoint():
    assert os.path.exists(pwd + '/statepoint.10.binary')

def test_results():
    os.system('python results.py')
    compare = filecmp.cmp('results_test.dat', 'results_true.dat')
    if not compare:
      os.rename('results_test.dat', 'results_error.dat')
    assert compare

def teardown():
    output = [pwd + '/statepoint.10.binary', pwd + '/results_test.dat']
    for f in output:
        if os.path.exists(f):
            os.remove(f)
