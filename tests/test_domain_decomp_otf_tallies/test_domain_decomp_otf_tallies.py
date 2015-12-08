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

def test_run():
    if opts.mpi_exec != '':
        # Run the 1-domain case
        os.chdir('1_domain')
        rundir = os.path.join(cwd, '1_domain')
        proc = Popen([opts.mpi_exec, '-np', '3', opts.exe, rundir],
               stderr=STDOUT, stdout=PIPE)
        print(proc.communicate()[0])
        returncode = proc.returncode
        assert returncode == 0, 'OpenMC did not exit successfully.'
        os.chdir('..')
        # Run the 4-domain case
        os.chdir('4_domains')
        rundir = os.path.join(cwd, '4_domains')
        proc = Popen([opts.mpi_exec, '-np', '5', opts.exe, rundir],
               stderr=STDOUT, stdout=PIPE)
        print(proc.communicate()[0])
        returncode = proc.returncode
        assert returncode == 0, 'OpenMC did not exit successfully.'
        os.chdir('..')
    else:
        # This test is only relevant if MPI is enabled
        return

def test_created_statepoint():
    if opts.mpi_exec == '':
        # This test is only relevant if MPI is enabled
        return
    rundir = os.path.join(cwd, '1_domain')
    statepoint = glob.glob(os.path.join(rundir, 'statepoint.20.*'))
    assert len(statepoint) == 1, 'Either multiple or no statepoint files exist.'
    assert statepoint[0].endswith('binary') or statepoint[0].endswith('h5'),\
        'Statepoint file is not a binary or hdf5 file.'
    rundir = os.path.join(cwd, '4_domains')
    statepoint = glob.glob(os.path.join(rundir, 'statepoint.20.*'))
    assert len(statepoint) == 4, 'Wrong number of statepoint files: %s' % \
        len(statepoint)
    assert statepoint[0].endswith('binary') or statepoint[0].endswith('h5'),\
        'Statepoint file is not a binary or hdf5 file.'

def test_output_exists():
    if opts.mpi_exec == '':
        # This test is only relevant if MPI is enabled
        return
    rundir = os.path.join(cwd, '1_domain')
    assert os.path.exists(os.path.join(rundir, 'tallies.domain_1.out')), 'Tally output file does not exist.'
    rundir = os.path.join(cwd, '4_domains')
    assert os.path.exists(os.path.join(rundir, 'tallies.domain_1.out')), 'Tally output file does not exist.'
    assert os.path.exists(os.path.join(rundir, 'tallies.domain_2.out')), 'Tally output file does not exist.'
    assert os.path.exists(os.path.join(rundir, 'tallies.domain_3.out')), 'Tally output file does not exist.'
    assert os.path.exists(os.path.join(rundir, 'tallies.domain_4.out')), 'Tally output file does not exist.'

def test_results():
    if opts.mpi_exec == '':
        # This test is only relevant if MPI is enabled
        return
    os.chdir('1_domain')
    rundir = os.path.join(cwd, '1_domain')
    statepoint = glob.glob(os.path.join(rundir, 'statepoint.20.*'))
    call(['python', 'results.py', statepoint[0]])
    compare = filecmp.cmp('results_test.dat', 'results_true.dat')
    if not compare:
      os.rename('results_test.dat', 'results_error.dat')
    assert compare, 'Results do not agree.'
    os.chdir('..')
    os.chdir('4_domains')
    rundir = os.path.join(cwd, '4_domains')
    statepoint = glob.glob(os.path.join(rundir, 'statepoint.20.*'))
    call(['python', 'results.py', statepoint[0]])
    compare = filecmp.cmp('results_test.dat', 'results_true.dat')
    if not compare:
      os.rename('results_test.dat', 'results_error.dat')
    assert compare, 'Results do not agree.'
    os.chdir('..')
    result1 = os.path.join('1_domain','results_test.dat')
    result4 = os.path.join('4_domains','results_test.dat')
    compare = filecmp.cmp(result1, result4)
    assert compare, '4-domain case does not match 1-domain case: ' + \
                    'rng reproducibility is broken'

def teardown():
    if opts.mpi_exec == '':
        # This test is only relevant if MPI is enabled
        return
    rundir = os.path.join(cwd, '1_domain')
    output = glob.glob(os.path.join(rundir, 'statepoint.20.*'))
    output += glob.glob(os.path.join(rundir, 'tallies.domain*'))
    output.append(os.path.join(rundir, 'results_test.dat'))
    rundir = os.path.join(cwd, '4_domains')
    output += glob.glob(os.path.join(rundir, 'statepoint.20.*'))
    output += glob.glob(os.path.join(rundir, 'tallies.domain*'))
    output.append(os.path.join(rundir, 'results_test.dat'))
    for f in output:
        if os.path.exists(f):
            os.remove(f)

if __name__ == '__main__':

    # test for openmc executable
    if opts.exe is None:
        raise Exception('Must specify OpenMC executable from command line with --exe.')

    # test for mpi executable
    if opts.mpi_exec is '':
        raise Exception('Must specify MPI executable from command line with --mpi_exec.')

    # run tests
    try:
        test_run()
        test_created_statepoint()
        test_output_exists()
        test_results()
    finally:
        teardown()
