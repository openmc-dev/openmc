#!/usr/bin/env python

import os
from subprocess import Popen, STDOUT, PIPE, call
import filecmp
import glob

from nose.plugins.skip import SkipTest

from nose_mpi import NoseMPI

pwd = os.path.dirname(__file__)

skipAll = False


def setup():
    os.putenv('PWD', pwd)
    os.chdir(pwd)
    global skipAll
    skipAll = False


def test_run():
    openmc_path = pwd + '/../../src/openmc'
    if int(NoseMPI.mpi_np) > 0:
        proc = Popen([NoseMPI.mpi_exec, '-np', NoseMPI.mpi_np, openmc_path],
                     stderr=STDOUT, stdout=PIPE)
    else:
        proc = Popen([openmc_path], stderr=STDOUT, stdout=PIPE)
    output = proc.communicate()[0]
    print(output)
    returncode = proc.returncode
    if 'CMFD is not available' in output:
        global skipAll
        skipAll = True
        raise SkipTest
    assert returncode == 0


def test_created_statepoint():
    if skipAll:
        raise SkipTest
    statepoint = glob.glob(pwd + '/statepoint.20.*')
    assert len(statepoint) == 1
    assert statepoint[0].endswith('binary') or statepoint[0].endswith('h5')


def test_output_exists():
    if skipAll:
        raise SkipTest
    assert os.path.exists(pwd + '/tallies.out')


def test_results():
    if skipAll:
        raise SkipTest
    statepoint = glob.glob(pwd + '/statepoint.20.*')
    call(['python', 'results.py', statepoint[0]])
    compare = filecmp.cmp('results_test.dat', 'results_true.dat')
    if not compare:
        os.rename('results_test.dat', 'results_error.dat')
    assert compare


def teardown():
    output = glob.glob(pwd + '/statepoint.20.*')
    output.append(pwd + '/tallies.out')
    output.append(pwd + '/results_test.dat')
    for f in output:
        if os.path.exists(f):
            os.remove(f)
