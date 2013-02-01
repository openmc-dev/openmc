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
    assert os.path.exists(pwd + '/statepoint.7.binary')

def test_restart_form1():
    proc = Popen([pwd + '/../../src/openmc', '-r',
                  pwd + '/statepoint.7.binary'],
                 stderr=STDOUT, stdout=PIPE)
    returncode = proc.wait()
    print(proc.communicate()[0])
    assert returncode == 0

def test_restart_form2():
    proc = Popen([pwd + '/../../src/openmc', '--restart',
                  pwd + '/statepoint.7.binary'],
                 stderr=STDOUT, stdout=PIPE)
    returncode = proc.wait()
    print(proc.communicate()[0])
    assert returncode == 0

def teardown():
    output = [pwd + '/statepoint.7.binary']
    for f in output:
        if os.path.exists(f):
            os.remove(f)
