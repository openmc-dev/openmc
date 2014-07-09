#!/usr/bin/env python

import os
from subprocess import Popen, STDOUT, PIPE, call
import filecmp
import glob
from optparse import OptionParser

parser = OptionParser()
parser.add_option('--mpi_exec', dest='mpi_exec', default='')
parser.add_option('--mpi_np', dest='mpi_np', default='3')
parser.add_option('--exe', dest='exe')
(opts, args) = parser.parse_args()
cwd = os.getcwd()

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

def test_run1():
    if opts.mpi_exec != '':
        proc = Popen([opts.mpi_exec, '-np', opts.mpi_np, opts.exe, cwd],
               stderr=STDOUT, stdout=PIPE)
    else:
        proc = Popen([opts.exe, cwd], stderr=STDOUT, stdout=PIPE)
    print(proc.communicate()[0])
    returncode = proc.returncode
    assert returncode == 0, 'OpenMC did not exit successfully.'

def test_statepoint_exists():
    statepoint = glob.glob(os.path.join(cwd, 'statepoint.10.*'))
    assert len(statepoint) == 1, 'Either multiple or no statepoint files exist.'
    assert statepoint[0].endswith('binary') or statepoint[0].endswith('h5'),\
        'Statepoint file is not a binary or hdf5 file.'
    source = glob.glob(os.path.join(cwd, 'source.10.*'))
    assert len(statepoint) == 1, 'Either multple or no source files exist.'
    assert source[0].endswith('binary') or source[0].endswith('h5'),\
        'Source file is not a binary or hdf5 file.'

def test_run2():
    openmc_path = os.path.join(cwd, '../../src/openmc')
    source = glob.glob(os.path.join(cwd, 'source.10.*'))
    with open('settings.xml','w') as fh:
        fh.write(settings2.format(source[0].split('.')[-1]))
    if opts.mpi_exec != '':
        proc = Popen([opts.mpi_exec, '-np', opts.mpi_np, opts.exe, cwd],
               stderr=STDOUT, stdout=PIPE)
    else:
        proc = Popen([opts.exe, cwd], stderr=STDOUT, stdout=PIPE)
    print(proc.communicate()[0])
    returncode = proc.returncode
    assert returncode == 0, 'OpenMC did not exit successfully.'

def test_results():
    statepoint = glob.glob(os.path.join(cwd, 'statepoint.10.*'))
    call(['python', 'results.py', statepoint[0]])
    compare = filecmp.cmp('results_test.dat', 'results_true.dat')
    if not compare:
      os.rename('results_test.dat', 'results_error.dat')
    assert compare, 'Results do not agree.'

def teardown():
    with open('settings.xml','w') as fh:
        fh.write(settings1)
    output = glob.glob(os.path.join(cwd, 'statepoint.10.*'))
    output += glob.glob(os.path.join(cwd, 'source.10.*'))
    output.append(os.path.join(cwd, 'results_test.dat'))
    for f in output:
        if os.path.exists(f):
            os.remove(f)

if __name__ == '__main__':

    # test for openmc executable
    if opts.exe is None:
        raise Exception('Must specify OpenMC executable from command line with --exe.')

    # run tests
    try:
        test_run1()
        test_statepoint_exists()
        test_run2()
        test_results()
    finally:
        teardown()
