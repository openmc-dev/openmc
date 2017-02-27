from __future__ import print_function
import subprocess
from numbers import Integral
import sys

from six import string_types

_summary_indicator = "TIMING STATISTICS"


def _run(command, output, cwd):
    # Launch a subprocess
    p = subprocess.Popen(command, shell=True, cwd=cwd, stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT, universal_newlines=True)

    storage_flag = False

    # Capture and re-print OpenMC output in real-time
    while True:
        # If OpenMC is finished, break loop
        line = p.stdout.readline()
        if not line and p.poll() is not None:
            break

        if output == 'full':
            # If user requested output, print to screen
            print(line, end='')
        elif output == 'summary' and _summary_indicator in line:
            # If they requested a summary, look for the start of the summary
            storage_flag = True

        if storage_flag:
            # If a summary is requested, and we have reached the summary,
            # then print it
            print(line, end='')

    # Return the returncode (integer, zero if no problems encountered)
    return p.returncode


def plot_geometry(output=True, openmc_exec='openmc', cwd='.'):
    """Run OpenMC in plotting mode

    Parameters
    ----------
    output : bool
        Capture OpenMC output from standard out
    openmc_exec : str
        Path to OpenMC executable
    cwd : str, optional
        Path to working directory to run in. Defaults to the current working directory.

    """
    return _run(openmc_exec + ' -p', output, cwd)


def run_volume_calculation(output=True, openmc_exec='openmc', cwd='.'):
    """Run stochastic volume calculations in OpenMC

    Parameters
    ----------
    output : bool
        Capture OpenMC output from standard out
    openmc_exec : str
        Path to OpenMC executable
    cwd : str, optional
        Path to working directory to run in. Defaults to the current working directory.

    """
    return _run(openmc_exec + ' --volume', output, cwd)



def run(particles=None, threads=None, geometry_debug=False,
        restart_file=None, tracks=False, output='full', cwd='.',
        openmc_exec='openmc', mpi_args=None):
    """Run an OpenMC simulation.

    Parameters
    ----------
    particles : int, optional
        Number of particles to simulate per generation.
    threads : int, optional
        Number of OpenMP threads. If OpenMC is compiled with OpenMP threading
        enabled, the default is implementation-dependent but is usually equal to
        the number of hardware threads available (or a value set by the
        :envvar:`OMP_NUM_THREADS` environment variable).
    geometry_debug : bool, optional
        Turn on geometry debugging during simulation. Defaults to False.
    restart_file : str, optional
        Path to restart file to use
    tracks : bool, optional
        Write tracks for all particles. Defaults to False.
    output : {"full", "summary", "none", False}, optional
        Degree of OpenMC output captured from standard out. "full" prints all
        output; "summary" prints only the results summary, and "none" or False
        does not show the output. Defaults to "full".
    cwd : str, optional
        Path to working directory to run in. Defaults to the current working
        directory.
    openmc_exec : str, optional
        Path to OpenMC executable. Defaults to 'openmc'.
    mpi_args : list of str, optional
        MPI execute command and any additional MPI arguments to pass,
        e.g. ['mpiexec', '-n', '8'].

    """

    post_args = ' '
    pre_args = ''

    if isinstance(particles, Integral) and particles > 0:
        post_args += '-n {0} '.format(particles)

    if isinstance(threads, Integral) and threads > 0:
        post_args += '-s {0} '.format(threads)

    if geometry_debug:
        post_args += '-g '

    if isinstance(restart_file, string_types):
        post_args += '-r {0} '.format(restart_file)

    if tracks:
        post_args += '-t'

    if mpi_args is not None:
        pre_args = ' '.join(mpi_args) + ' '

    command = pre_args + openmc_exec + ' ' + post_args

    return _run(command, output, cwd)
