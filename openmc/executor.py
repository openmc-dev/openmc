from collections.abc import Iterable
from numbers import Integral
import os
import subprocess

import openmc
from .plots import _get_plot_image


def _process_CLI_arguments(volume=False, geometry_debug=False, particles=None,
                           plot=False, restart_file=None, threads=None,
                           tracks=False, event_based=None,
                           openmc_exec='openmc', mpi_args=None, path_input=None):
    """Converts user-readable flags in to command-line arguments to be run with
    the OpenMC executable via subprocess.

    Parameters
    ----------
    volume : bool, optional
        Run in stochastic volume calculation mode. Defaults to False.
    geometry_debug : bool, optional
        Turn on geometry debugging during simulation. Defaults to False.
    particles : int, optional
        Number of particles to simulate per generation.
    plot : bool, optional
        Run in plotting mode. Defaults to False.
    restart_file : str or PathLike
        Path to restart file to use
    threads : int, optional
        Number of OpenMP threads. If OpenMC is compiled with OpenMP threading
        enabled, the default is implementation-dependent but is usually equal
        to the number of hardware threads available (or a value set by the
        :envvar:`OMP_NUM_THREADS` environment variable).
    tracks : bool, optional
        Enables the writing of particles tracks. The number of particle
        tracks written to tracks.h5 is limited to 1000 unless
        Settings.max_tracks is set. Defaults to False.
    event_based : None or bool, optional
        Turns on event-based parallelism if True. If None, the value in
        the Settings will be used.
    openmc_exec : str, optional
        Path to OpenMC executable. Defaults to 'openmc'.
    mpi_args : list of str, optional
        MPI execute command and any additional MPI arguments to pass,
        e.g., ['mpiexec', '-n', '8'].
    path_input : str or PathLike
        Path to a single XML file or a directory containing XML files for the
        OpenMC executable to read.

    .. versionadded:: 0.13.0

    Returns
    -------
    args : Iterable of str
        The runtime flags converted to CLI arguments of the OpenMC executable

    """

    args = [openmc_exec]

    if volume:
        args.append('--volume')

    if isinstance(particles, Integral) and particles > 0:
        args += ['-n', str(particles)]

    if isinstance(threads, Integral) and threads > 0:
        args += ['-s', str(threads)]

    if geometry_debug:
        args.append('-g')

    if event_based is not None:
        if event_based:
            args.append('-e')

    if isinstance(restart_file, (str, os.PathLike)):
        args += ['-r', str(restart_file)]

    if tracks:
        args.append('-t')

    if plot:
        args.append('-p')

    if mpi_args is not None:
        args = mpi_args + args

    if path_input is not None:
        args += [path_input]

    return args


def _run(args, output, cwd):
    # Optionally print debug information about the invocation when
    # OPENMC_DEBUG is set in the environment. Helps diagnosing cases
    # where OpenMC cannot find XML input files because of a cwd/path mixup.
    if os.environ.get('OPENMC_DEBUG'):
        try:
            print('\n[OPENMC_DEBUG] launching OpenMC')
            print(f'[OPENMC_DEBUG] args: {args}')
            print(f'[OPENMC_DEBUG] cwd: {cwd}')
            print('[OPENMC_DEBUG] files in cwd:')
            for p in sorted(os.listdir(cwd)):
                print('  ', p)
            print('\n')
        except Exception:            
            pass

    # Launch a subprocess
    p = subprocess.Popen(args, cwd=cwd, stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT, universal_newlines=True)

    # Capture and re-print OpenMC output in real-time
    lines = []
    while True:
        # If OpenMC is finished, break loop
        line = p.stdout.readline()
        if not line and p.poll() is not None:
            break

        lines.append(line)
        if output:
            # If user requested output, print to screen
            print(line, end='')

    # Raise an exception if return status is non-zero
    if p.returncode != 0:
        # Get error message from output and simplify whitespace
        output = ''.join(lines)
        if 'ERROR: ' in output:
            _, _, error_msg = output.partition('ERROR: ')
        elif 'what()' in output:
            _, _, error_msg = output.partition('what(): ')
        else:
            error_msg = 'OpenMC aborted unexpectedly.'
        error_msg = ' '.join(error_msg.split())

        raise RuntimeError(error_msg)


def plot_geometry(output=True, openmc_exec='openmc', cwd='.', path_input=None):
    """Run OpenMC in plotting mode

    Parameters
    ----------
    output : bool, optional
        Capture OpenMC output from standard out
    openmc_exec : str, optional
        Path to OpenMC executable
    cwd : str, optional
        Path to working directory to run in
    path_input : str
        Path to a single XML file or a directory containing XML files for the
        OpenMC executable to read.

        .. versionadded:: 0.13.3

    Raises
    ------
    RuntimeError
        If the `openmc` executable returns a non-zero status

    """
    args = [openmc_exec, '-p']
    if path_input is not None:
        args += [path_input]
    _run(args, output, cwd)


def plot_inline(plots, openmc_exec='openmc', cwd='.', path_input=None):
    """Display plots inline in a Jupyter notebook.

    .. versionchanged:: 0.13.0
        The *convert_exec* argument was removed since OpenMC now produces
        .png images directly.


    Parameters
    ----------
    plots : Iterable of openmc.Plot
        Plots to display
    openmc_exec : str
        Path to OpenMC executable
    cwd : str, optional
        Path to working directory to run in
    path_input : str
        Path to a single XML file or a directory containing XML files for the
        OpenMC executable to read.

        .. versionadded:: 0.13.3

    Raises
    ------
    RuntimeError
        If the `openmc` executable returns a non-zero status

    """
    from IPython.display import display

    if not isinstance(plots, Iterable):
        plots = [plots]

    # Create plots.xml
    openmc.Plots(plots).export_to_xml(cwd)

    # Run OpenMC in geometry plotting mode
    plot_geometry(False, openmc_exec, cwd, path_input)

    if plots is not None:
        images = [_get_plot_image(p, cwd) for p in plots]
        display(*images)


def calculate_volumes(threads=None, output=True, cwd='.',
                      openmc_exec='openmc', mpi_args=None,
                      path_input=None):
    """Run stochastic volume calculations in OpenMC.

    This function runs OpenMC in stochastic volume calculation mode. To specify
    the parameters of a volume calculation, one must first create a
    :class:`openmc.VolumeCalculation` instance and assign it to
    :attr:`openmc.Settings.volume_calculations`. For example:

    >>> vol = openmc.VolumeCalculation(domains=[cell1, cell2], samples=100000)
    >>> settings = openmc.Settings()
    >>> settings.volume_calculations = [vol]
    >>> settings.export_to_xml()
    >>> openmc.calculate_volumes()

    Parameters
    ----------
    threads : int, optional
        Number of OpenMP threads. If OpenMC is compiled with OpenMP threading
        enabled, the default is implementation-dependent but is usually equal
        to the number of hardware threads available (or a value set by the
        :envvar:`OMP_NUM_THREADS` environment variable).
    output : bool, optional
        Capture OpenMC output from standard out
    openmc_exec : str, optional
        Path to OpenMC executable. Defaults to 'openmc'.
    mpi_args : list of str, optional
        MPI execute command and any additional MPI arguments to pass,
        e.g., ['mpiexec', '-n', '8'].
    cwd : str, optional
        Path to working directory to run in. Defaults to the current working
        directory.
    path_input : str or PathLike
        Path to a single XML file or a directory containing XML files for the
        OpenMC executable to read.


    Raises
    ------
    RuntimeError
        If the `openmc` executable returns a non-zero status

    See Also
    --------
    openmc.VolumeCalculation

    """

    args = _process_CLI_arguments(volume=True, threads=threads,
                                  openmc_exec=openmc_exec, mpi_args=mpi_args,
                                  path_input=path_input)

    _run(args, output, cwd)


def run(particles=None, threads=None, geometry_debug=False,
        restart_file=None, tracks=False, output=True, cwd='.',
        openmc_exec='openmc', mpi_args=None, event_based=False,
        path_input=None):
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
    restart_file : str or PathLike
        Path to restart file to use
    tracks : bool, optional
        Enables the writing of particles tracks. The number of particle tracks
        written to tracks.h5 is limited to 1000 unless Settings.max_tracks is
        set. Defaults to False.
    output : bool
        Capture OpenMC output from standard out
    cwd : str, optional
        Path to working directory to run in. Defaults to the current working
        directory.
    openmc_exec : str, optional
        Path to OpenMC executable. Defaults to 'openmc'.
    mpi_args : list of str, optional
        MPI execute command and any additional MPI arguments to pass, e.g.,
        ['mpiexec', '-n', '8'].
    event_based : bool, optional
        Turns on event-based parallelism, instead of default history-based

        .. versionadded:: 0.12

    path_input : str or PathLike
        Path to a single XML file or a directory containing XML files for the
        OpenMC executable to read.

        .. versionadded:: 0.13.3

    Raises
    ------
    RuntimeError
        If the `openmc` executable returns a non-zero status

    """

    args = _process_CLI_arguments(
        volume=False, geometry_debug=geometry_debug, particles=particles,
        restart_file=restart_file, threads=threads, tracks=tracks,
        event_based=event_based, openmc_exec=openmc_exec, mpi_args=mpi_args,
        path_input=path_input)

    _run(args, output, cwd)
