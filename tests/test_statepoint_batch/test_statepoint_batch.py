#!/usr/bin/env python

import os
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
    assert os.path.exists(pwd + '/statepoint.3.binary')
    assert os.path.exists(pwd + '/statepoint.6.binary')
    assert os.path.exists(pwd + '/statepoint.9.binary')

def test_results():
    os.system('python results.py')
    compare = filecmp.cmp('results_test.dat', 'results_true.dat')
    if not compare:
      os.rename('results_test.dat', 'results_error.dat')
    assert compare

def teardown():
    output = [pwd + i for i in ['/statepoint.3.binary',
                                '/statepoint.6.binary',
                                '/statepoint.9.binary',
                                '/results_test.dat']]
    for f in output:
        if os.path.exists(f):
            os.remove(f)
