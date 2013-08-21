#!/usr/bin/env python

import os
from subprocess import Popen, STDOUT, PIPE
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
    returncode = proc.wait()
    print(proc.communicate()[0])
    assert returncode == 0

def test_created_statepoint():
    statepoint = glob.glob('statepoint.10.*')
    assert len(statepoint) == 1
    assert statepoint[0].endswith('binary') or statepoint[0].endswith('h5')

def test_results():
    statepoint = glob.glob('statepoint.10.*')
    os.system('python results.py {0}'.format(statepoint[0]))
    compare = filecmp.cmp('results_test.dat', 'results_true.dat')
    if not compare:
      os.rename('results_test.dat', 'results_error.dat')
    assert compare

def teardown():
    statepoint = glob.glob('statepoint.10.*')
    output = [pwd + '/{0}'.format(statepoint[0]), pwd + '/results_test.dat']
    for f in output:
        if os.path.exists(f):
            os.remove(f)
