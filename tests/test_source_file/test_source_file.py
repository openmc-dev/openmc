#!/usr/bin/env python

import os
from subprocess import Popen, STDOUT, PIPE, call
import filecmp
from nose_mpi import NoseMPI
import glob

pwd = os.path.dirname(__file__)

settings1="""<?xml version="1.0"?>
<settings>

  <state_point batches="10" />
  <source_point separate="true" />

  <eigenvalue>
    <batches>10</batches>
    <inactive>5</inactive>
    <particles>1000</particles>
  </eigenvalue>

  <source>
    <space type="box">
      <parameters>-4 -4 -4  4  4  4</parameters>
    </space>
  </source>

</settings>
"""

settings2 = """<?xml version="1.0"?>
<settings>

  <eigenvalue>
    <batches>10</batches>
    <inactive>5</inactive>
    <particles>1000</particles>
  </eigenvalue>

  <source>
    <file> source.10.{0} </file>
  </source>

</settings>
"""

def setup(): 
    os.putenv('PWD', pwd)
    os.chdir(pwd)

def test_run1():
    openmc_path = pwd + '/../../src/openmc'
    if int(NoseMPI.mpi_np) > 0:
        proc = Popen([NoseMPI.mpi_exec, '-np', NoseMPI.mpi_np, openmc_path],
               stderr=STDOUT, stdout=PIPE)
    else:
        proc = Popen([openmc_path], stderr=STDOUT, stdout=PIPE)
    print(proc.communicate()[0])
    returncode = proc.returncode
    assert returncode == 0

def test_statepoint_exists():
    statepoint = glob.glob(pwd + '/statepoint.10.*')
    assert len(statepoint) == 1
    assert statepoint[0].endswith('binary') or statepoint[0].endswith('h5')
    source = glob.glob(pwd + '/source.10.*')
    assert len(statepoint) == 1
    assert source[0].endswith('binary') or source[0].endswith('h5')

def test_run2():
    openmc_path = pwd + '/../../src/openmc'
    source = glob.glob(pwd + '/source.10.*')
    with open('settings.xml','w') as fh:
        fh.write(settings2.format(source[0].split('.')[2]))
    if int(NoseMPI.mpi_np) > 0:
        proc = Popen([NoseMPI.mpi_exec, '-np', NoseMPI.mpi_np, openmc_path],
               stderr=STDOUT, stdout=PIPE)
    else:
        proc = Popen([openmc_path], stderr=STDOUT, stdout=PIPE)
    print(proc.communicate()[0])
    returncode = proc.returncode
    assert returncode == 0

def test_results():
    statepoint = glob.glob(pwd + '/statepoint.10.*')
    call(['python', 'results.py', statepoint[0]])
    compare = filecmp.cmp('results_test.dat', 'results_true.dat')
    if not compare:
      os.rename('results_test.dat', 'results_error.dat')
    assert compare

def teardown():
    with open('settings.xml','w') as fh:
        fh.write(settings1)
    output = glob.glob(pwd + '/statepoint.10.*')
    output += glob.glob(pwd + '/source.10.*')
    output.append(pwd + '/results_test.dat')
    for f in output:
        if os.path.exists(f):
            os.remove(f)
