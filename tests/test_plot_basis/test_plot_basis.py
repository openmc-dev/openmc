#!/usr/bin/env python

import os
from subprocess import Popen, STDOUT, PIPE

pwd = os.path.dirname(__file__)

def setup(): 
    os.putenv('PWD', pwd)
    os.chdir(pwd)

def test_run():
    proc = Popen([pwd + '/../../src/openmc', '-p'], stderr=STDOUT, stdout=PIPE)
    returncode = proc.wait()
    print(proc.communicate()[0])
    assert returncode == 0

def test_plots_exists():
    assert os.path.exists(pwd + '/1_plot.ppm')
    assert os.path.exists(pwd + '/2_plot.ppm')
    assert os.path.exists(pwd + '/3_plot.ppm')

def teardown():
    output = [pwd + '/1_plot.ppm', pwd + '/2_plot.ppm', pwd + '/3_plot.ppm']
    for f in output:
        if os.path.exists(f):
            os.remove(f)
