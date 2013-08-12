#!/usr/bin/env python

import os
from subprocess import Popen, STDOUT, PIPE
import filecmp

pwd = os.path.dirname(__file__)

def setup(): 
    os.putenv('PWD', pwd)
    os.chdir(pwd)

def test_run():
    proc = Popen([pwd + '/../../src/openmc'], stderr=PIPE, stdout=PIPE)
    stdout, stderr = proc.communicate()
    assert stderr != ''

def test_created_restart():
    assert os.path.exists(pwd + '/particle_10_638.binary')

def test_results():
    os.system('python results.py')
    compare = filecmp.cmp('results_test.dat', 'results_true.dat')
    if not compare:
      os.system('cp results_test.dat results_error.dat')
    assert compare

def test_run_restart():
    proc = Popen([pwd + '/../../src/openmc -s particle_10_638.binary'], 
           stderr=PIPE, stdout=PIPE, shell=True)
    stdout, stderr = proc.communicate()
    assert stderr != ''

def teardown():
    output = [pwd + '/particle_10_638.binary', pwd + '/results_test.dat']
    for f in output:
        if os.path.exists(f):
            os.remove(f)
