#!/usr/bin/env python

import os
import glob
from subprocess import Popen, STDOUT, PIPE
import filecmp

pwd = os.path.dirname(__file__)

def setup(): 
    os.putenv('PWD', pwd)
    os.chdir(pwd)

def test_run():
    proc = Popen([pwd + '/../../src/openmc'], stderr=STDOUT, stdout=PIPE)
    returncode = proc.wait()
    print(proc.communicate()[0])
    assert returncode == 0

def test_statepoints_exist():
    assert os.path.exists(pwd + '/statepoint.2.binary')
    assert os.path.exists(pwd + '/statepoint.4.binary')
    assert os.path.exists(pwd + '/statepoint.6.binary')
    assert os.path.exists(pwd + '/statepoint.8.binary')
    assert os.path.exists(pwd + '/statepoint.10.binary')

def test_results():
    os.system('python results.py')
    compare = filecmp.cmp('results_test.dat', 'results_true.dat')
    if not compare:
      os.system('cp results_test.dat results_error.dat')
    assert compare

def teardown():
    output = glob.glob(pwd + '/statepoint.*.binary')
    output.append(pwd + '/results_test.dat')
    for f in output:
        if os.path.exists(f):
            os.remove(f)
