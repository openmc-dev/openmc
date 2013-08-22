#!/usr/bin/env python

import os
from subprocess import Popen, STDOUT, PIPE, call
import filecmp
from nose_mpi import NoseMPI
import glob

pwd = os.path.dirname(__file__)

def setup(): 
    os.putenv('PWD', pwd)
    os.chdir(pwd)

def test_run():
    openmc_path = pwd + '/../../src/openmc'
    if int(NoseMPI.mpi_np) > 0:
        proc = Popen([NoseMPI.mpi_exec, '-np', NoseMPI.mpi_np, openmc_path],
               stderr=PIPE, stdout=PIPE)
    else:
        proc = Popen([openmc_path], stderr=PIPE, stdout=PIPE)
    stdout, stderr = proc.communicate()
    assert stderr != '' 

def test_created_restart():
    particle = glob.glob(pwd + '/particle_10_903.*')
    assert len(particle) == 1
    assert particle[0].endswith('binary') or \
           particle[0].endswith('h5')
    particle = glob.glob(pwd + '/particle_8_777.*')
    assert len(particle) == 1
    assert particle[0].endswith('binary') or \
           particle[0].endswith('h5')

def test_results():
    particle = glob.glob(pwd + '/particle_10_903.*')
    call(['python', 'results.py', particle[0]])
    compare = filecmp.cmp('results_test.dat', 'results_true.dat')
    if not compare:
      os.rename('results_test.dat', 'results_error.dat')
    assert compare

def test_run_restart():
    particle = glob.glob(pwd + '/particle_10_903.*')
    proc = Popen([pwd + '/../../src/openmc -r ' + particle[0]],
           stderr=PIPE, stdout=PIPE, shell=True)
    stdout, stderr = proc.communicate()
    assert stderr == ''

def teardown():
    output = glob.glob(pwd + '/statepoint.*') + \
             glob.glob(pwd + '/particle_*') + \
             [pwd + '/results_test.dat']
    for f in output:
        if os.path.exists(f):
            os.remove(f)
