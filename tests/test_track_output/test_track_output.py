#!/usr/bin/env python

import os
from subprocess import Popen, STDOUT, PIPE, call
import filecmp
import glob

from nose.plugins.skip import SkipTest

from nose_mpi import NoseMPI

pwd = os.path.dirname(__file__)


def setup():
    os.putenv('PWD', pwd)
    os.chdir(pwd)


def test_run():
    openmc_path = pwd + '/../../src/openmc'
    if int(NoseMPI.mpi_np) > 0:
        proc = Popen([NoseMPI.mpi_exec, '-np', NoseMPI.mpi_np, openmc_path],
                     stderr=STDOUT, stdout=PIPE)
    else:
        proc = Popen([openmc_path], stderr=STDOUT, stdout=PIPE)
    returncode = proc.wait()
    print(proc.communicate()[0])
    assert returncode == 0


def test_created_outputs():
    outputs = [glob.glob(''.join((pwd, '/track_1_1_1.*')))]
    outputs.append(glob.glob(''.join((pwd, '/track_1_1_2.*'))))
    for files in outputs:
        assert len(files) == 1
        assert files[0].endswith('binary') or files[0].endswith('h5')


def test_outputs():
    # If vtk python module is not available, we can't run track.py so skip this
    # test
    try:
        import vtk
    except ImportError:
        raise SkipTest

    call(['../../src/utils/track.py', '-o', 'poly'] +
         glob.glob(''.join((pwd, '/track*'))))
    poly = ''.join((pwd, '/poly.pvtp'))
    assert os.path.isfile(poly)
    metric = ''.join((pwd, '/true_poly.pvtp'))
    compare = filecmp.cmp(poly, metric)
    if not compare:
        os.rename('poly.pvtp', 'error_poly.pvtp')
    assert compare


def teardown():
    temp_files = glob.glob(''.join((pwd, '/statepoint*')))
    temp_files = temp_files + glob.glob(''.join((pwd, '/track*')))
    temp_files = temp_files + glob.glob(''.join((pwd, '/poly*')))
    for f in temp_files:
        if os.path.exists(f):
            os.remove(f)
