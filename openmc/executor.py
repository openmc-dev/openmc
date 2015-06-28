import subprocess
import os

from openmc.checkvalue import *


class Executor(object):


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
            print(line),

            # If OpenMC is finished, break loop
            if line == '' and p.poll() != None:
                break

        # Return the returncode (integer, zero if no problems encountered)
        return p.returncode


    @property
    def working_directory(self):
        return self._working_directory


    @working_directory.setter
    def working_directory(self, working_directory):

        if not is_string(working_directory):
            msg = 'Unable to set Executor\'s working directory to {0} ' \
                  'since it is not a string'.format(working_directory)
            raise ValueError(msg)

        elif not os.path.isdir(working_directory):
            msg = 'Unable to set Executor\'s working directory to {0} ' \
                  'which does not exist'.format(working_directory)
            raise ValueError(msg)

        self._working_directory = working_directory


    def plot_geometry(self, output=True):
        return self._run_openmc('openmc -p', output)


    def run_simulation(self, particles=None, threads=None,
                       geometry_debug=False, restart_file=None,
                       tracks=False, mpi_procs=1, output=True):

        post_args = ' '
        pre_args = ''

        if is_integer(particles) and particles > 0:
            post_args += '-n {0} '.format(particles)

        if is_integer(threads) and threads > 0:
            post_args += '-s {0} '.format(threads)

        if geometry_debug:
            post_args += '-g '

        if is_string(restart_file):
            post_args += '-r {0} '.format(restart_file)

        if tracks:
            post_args += '-t'

        if is_integer(mpi_procs) and mpi_procs > 1:
            pre_args += 'mpirun -n {0} '.format(mpi_procs)

        command = pre_args + 'openmc ' + post_args

        return self._run_openmc(command, output)
