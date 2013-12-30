#!/usr/bin/env python

import os
from subprocess import Popen, STDOUT, PIPE
from nose_mpi import NoseMPI

pwd = os.path.dirname(__file__)

def setup(): 
    os.putenv('PWD', pwd)
    os.chdir(pwd)

def test_run():
    openmc_path = pwd + '/../../src/openmc'
    if int(NoseMPI.mpi_np) > 0:
        proc = Popen([NoseMPI.mpi_exec, '-np', NoseMPI.mpi_np, openmc_path, '-p'],
               stderr=STDOUT, stdout=PIPE)
    else:
        proc = Popen([openmc_path, '-p'], stderr=STDOUT, stdout=PIPE)
    print(proc.communicate()[0])
    returncode = proc.returncode
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
