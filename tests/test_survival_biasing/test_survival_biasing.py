#!/usr/bin/env python

import os
from subprocess import Popen, STDOUT, PIPE

pwd = os.path.dirname(__file__)

def setup(): 
    # Change directory
    os.putenv('PWD', pwd)
    os.chdir(pwd)

def test_run():
    # Run the simulation
    proc = Popen([pwd + '/../../src/openmc'], stderr=STDOUT, stdout=PIPE)
    returncode = proc.wait()

    # Display stdout
    print(proc.communicate()[0])

    # Make sure simulation ran to completion
    assert returncode == 0

def test_created_statepoint():
    # Make sure that statepoint file was created
    assert os.path.exists(pwd + '/statepoint.10.binary')

def teardown():
    # Remove output files
    output = [pwd + '/statepoint.10.binary']
    for f in output:
        if os.path.exists(f):
            os.remove(f)
