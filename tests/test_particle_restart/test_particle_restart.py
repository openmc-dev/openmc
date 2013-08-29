#!/usr/bin/env python

import os
from subprocess import Popen, STDOUT, PIPE
import filecmp
import glob

pwd = os.path.dirname(__file__)

def setup(): 
    os.putenv('PWD', pwd)
    os.chdir(pwd)

def test_run():
    proc = Popen([pwd + '/../../src/openmc'], stderr=PIPE, stdout=PIPE)
    stdout, stderr = proc.communicate()
    assert stderr != ''

def test_created_restart():
    assert os.path.exists(pwd + '/particle_10_394.binary')

def test_results():
    os.system('python results.py')
    compare = filecmp.cmp('results_test.dat', 'results_true.dat')
    if not compare:
      os.rename('results_test.dat', 'results_error.dat')
    assert compare

def test_run_restart():
    proc = Popen([pwd + '/../../src/openmc -r particle_10_394.binary'], 
           stderr=PIPE, stdout=PIPE, shell=True)
    stdout, stderr = proc.communicate()
    assert stderr == ''

def teardown():
    output = glob.glob(pwd + '/particle_*.binary') + [pwd + '/results_test.dat']
    for f in output:
        if os.path.exists(f):
            os.remove(f)
