#!/usr/bin/env python

import os
from subprocess import Popen, STDOUT, PIPE
from optparse import OptionParser

parser = OptionParser()
parser.add_option('--mpi_exec', dest='mpi_exec', default='')
parser.add_option('--mpi_np', dest='mpi_np', default='3')
parser.add_option('--exe', dest='exe')
(opts, args) = parser.parse_args()
cwd = os.getcwd()

def test_run():
    if opts.mpi_exec != '':
        proc = Popen([opts.mpi_exec, '-np', opts.mpi_np, opts.exe, '-p', cwd],
               stderr=STDOUT, stdout=PIPE)
    else:
        proc = Popen([opts.exe, '-p', cwd], stderr=STDOUT, stdout=PIPE)
    print(proc.communicate()[0])
    returncode = proc.returncode
    assert returncode == 0, 'OpenMC did not exit successfully.'

def test_plot_exists():
    assert os.path.exists(os.path.join(cwd ,'1_plot.ppm')), 'Plot ppm file does not exist.'

def teardown():
    output = [os.path.join(cwd ,'1_plot.ppm')]
    for f in output:
        if os.path.exists(f):
            os.remove(f)

if __name__ == '__main__':

    # test for openmc executable
    if opts.exe is None:
        raise Exception('Must specify OpenMC executable from command line with --exe.')

    # run tests
    try:
        test_run()
        test_plot_exists()
    finally:
        teardown()
