from __future__ import print_function
import subprocess
from numbers import Integral
import sys

if sys.version_info[0] >= 3:
    basestring = str


def _run(command, output, cwd):
    # Launch a subprocess
    p = subprocess.Popen(command, shell=True, cwd=cwd, stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT, universal_newlines=True)

    # Capture and re-print OpenMC output in real-time
    while True:
        # If OpenMC is finished, break loop
        line = p.stdout.readline()
        if not line and p.poll() != None:
            break

        # If user requested output, print to screen
        if output:
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


def run(particles=None, threads=None, geometry_debug=False,
        restart_file=None, tracks=False, mpi_procs=1, output=True,
        openmc_exec='openmc', mpi_exec='mpiexec', cwd='.'):
    """Run an OpenMC simulation.

    Parameters
    ----------
    particles : int, optional
        Number of particles to simulate per generation.
    threads : int, optional
        Number of OpenMP threads. If OpenMC is compiled with OpenMP threading
        enabled, the default is implementation-dependent but is usually equal to
        the number of hardware threads available (or a value set by the
        OMP_NUM_THREADS environment variable).
    geometry_debug : bool, optional
        Turn on geometry debugging during simulation. Defaults to False.
    restart_file : str, optional
        Path to restart file to use
    tracks : bool, optional
        Write tracks for all particles. Defaults to False.
    mpi_procs : int, optional
        Number of MPI processes.
    output : bool, optional
        Capture OpenMC output from standard out. Defaults to True.
    openmc_exec : str, optional
        Path to OpenMC executable. Defaults to 'openmc'.
    mpi_exec : str, optional
        MPI execute command. Defaults to 'mpiexec'.
    cwd : str, optional
        Path to working directory to run in. Defaults to the current working directory.

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
        pre_args += '{} -n {} '.format(mpi_exec, mpi_procs)

    command = pre_args + openmc_exec + ' ' + post_args

    return _run(command, output, cwd)
