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

##def test_outputs():
##    outputs = [glob.glob(''.join((pwd, '/track_1_1_1.*')))[0]]
##    outputs.append(glob.glob(''.join((pwd, '/track_1_1_2.*')))[0])
##    if outputs[0].endswith('binary'):
##        metrics = [''.join((pwd, '/track_1_1_1_true.binary')),
##                   ''.join((pwd, '/track_1_1_2_true.binary'))]
##    else:
##        metrics = [''.join((pwd, '/track_1_1_1_true.h5')),
##                   ''.join((pwd, '/track_1_1_2_true.h5'))]
##    for i in range(2):
##        output = outputs[i]
##        metric = metrics[i]
##        compare = filecmp.cmp(output, metric)
##        if not compare:
##            if outputs[0].endswith('binary'):
##                extension = 'binary'
##            else:
##                extension = 'h5'
##            os.rename(output, ''.join((pwd, '/track_1_1_{}_error.{}'.format(
##                i+1, extension))))
##        assert compare

def test_outputs():
    call(['../../src/utils/track.py', '-o', 'poly'] +
         glob.glob(''.join((pwd, '/track*'))))
    poly = ''.join((pwd, '/poly.pvtp'))
    assert os.path.isfile(poly)
    metric = ''.join((pwd, '/true_poly.pvtp'))
    compare = filecmp.cmp(poly, metric)
    if not compare:
        os.rename('poly.pvtp', 'poly_error.pvtp')
    assert compare

##def test_results():
##    statepoint = glob.glob(pwd + '/statepoint.10.*')
##    call(['python', 'results.py', statepoint[0]])
##    compare = filecmp.cmp('results_test.dat', 'results_true.dat')
##    if not compare:
##      os.rename('results_test.dat', 'results_error.dat')
##    assert compare

def teardown():
    temp_files = glob.glob(''.join((pwd, '/statepoint*')))
    temp_files = temp_files + glob.glob(''.join((pwd, '/track*')))
    temp_files = temp_files + glob.glob(''.join((pwd, '/poly*')))
    for f in temp_files:
        if os.path.exists(f):
            os.remove(f)
