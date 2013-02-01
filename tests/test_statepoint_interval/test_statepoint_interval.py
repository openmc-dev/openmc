#!/usr/bin/env python

import os
import glob
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

def test_statepoints_exist():
    assert os.path.exists(pwd + '/statepoint.2.binary')
    assert os.path.exists(pwd + '/statepoint.4.binary')
    assert os.path.exists(pwd + '/statepoint.6.binary')
    assert os.path.exists(pwd + '/statepoint.8.binary')
    assert os.path.exists(pwd + '/statepoint.10.binary')

def teardown():
    output = glob.glob(pwd + '/statepoint.*.binary')
    for f in output:
        if os.path.exists(f):
            os.remove(f)
