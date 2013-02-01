#!/usr/bin/env python

import os
from subprocess import Popen, STDOUT, PIPE

pwd = os.path.dirname(__file__)

def setup(): 
    os.putenv('PWD', pwd)
    os.chdir(pwd)

def test_run():
    proc = Popen([pwd + '/../../src/openmc'], stderr=STDOUT, stdout=PIPE)
    returncode = proc.wait()
    print(proc.communicate()[0])
    assert returncode == 0

def test_statepoint_exists():
    assert os.path.exists(pwd + '/statepoint.10.binary')

def test_output_exists():
    assert os.path.exists(pwd + '/tallies.out')

def teardown():
    output = [pwd + '/statepoint.10.binary', pwd + '/tallies.out']
    for f in output:
        if os.path.exists(f):
            os.remove(f)
