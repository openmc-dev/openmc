from __future__ import print_function
import subprocess
from numbers import Integral
import os
import sys

from openmc.checkvalue import check_type

if sys.version_info[0] >= 3:
    basestring = str


class Executor(object):
    """Control execution of OpenMC

    Attributes
    ----------
    working_directory : str
        Path to working directory to run in

    """

    def __init__(self):
        self._working_directory = '.'

    def _run_openmc(self, command, output):
        # Launch a subprocess to run OpenMC
        p = subprocess.Popen(command, shell=True,
                             cwd=self._working_directory,
                             stdout=subprocess.PIPE)

        # Capture and re-print OpenMC output in real-time
        while (True and output):
            line = p.stdout.readline()
            print(line, end='')

            # If OpenMC is finished, break loop
            if line == '' and p.poll() != None:
                break

    @property
    def working_directory(self):
        return self._working_directory

    @working_directory.setter
    def working_directory(self, working_directory):
        check_type("Executor's working directory", working_directory,
                   basestring)
        if not os.path.isdir(working_directory):
            msg = 'Unable to set Executor\'s working directory to {0} ' \
                  'which does not exist'.format(working_directory)
            raise ValueError(msg)

        self._working_directory = working_directory

    def plot_geometry(self, output=True):
        """Run OpenMC in plotting mode"""

        self._run_openmc('openmc -p', output)

    def run_simulation(self, particles=None, threads=None,
                       geometry_debug=False, restart_file=None,
                       tracks=False, mpi_procs=1, output=True):
        """Run an OpenMC simulation.

        Parameters
        ----------
        particles : int
            Number of particles to simulate per generation
        threads : int
            Number of OpenMP threads
        geometry_debug : bool
            Turn on geometry debugging during simulation
        restart_file : str
            Path to restart file to use
        tracks : bool
            Write tracks for all particles
        mpi_procs : int
            Number of MPI processes
        output : bool
            Capture OpenMC output from standard out

        """

        post_args = ' '
        pre_args = ''

        if isinstance(particles, Integral) and particles > 0:
            post_args += '-n {0} '.format(particles)

        if isinstance(threads, Integral) and threads > 0:
            post_args += '-s {0} '.format(threads)

        if geometry_debug:
            post_args += '-g '

        if isinstance(restart_file, basestring):
            post_args += '-r {0} '.format(restart_file)

        if tracks:
            post_args += '-t'

        if isinstance(mpi_procs, Integral) and mpi_procs > 1:
            pre_args += 'mpirun -n {0} '.format(mpi_procs)

        command = pre_args + 'openmc ' + post_args

        self._run_openmc(command, output)
