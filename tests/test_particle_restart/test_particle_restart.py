#!/usr/bin/env python

import os
from subprocess import Popen, STDOUT, PIPE

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

def test_run_restart():
    proc = Popen([pwd + '/../../src/openmc -s particle_10_638.binary'], 
           stderr=PIPE, stdout=PIPE, shell=True)
    stdout, stderr = proc.communicate()
    assert stderr != ''

def teardown():
    output = [pwd + '/particle_10_638.binary']
    for f in output:
        if os.path.exists(f):
            os.remove(f)
