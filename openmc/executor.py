from collections.abc import Iterable
from numbers import Integral
import subprocess

import openmc


def _run(args, output, cwd):
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
        raise subprocess.CalledProcessError(p.returncode, ' '.join(args),
                                            ''.join(lines))


def plot_geometry(output=True, openmc_exec='openmc', cwd='.'):
    """Run OpenMC in plotting mode

    Parameters
    ----------
    output : bool, optional
        Capture OpenMC output from standard out
    openmc_exec : str, optional
        Path to OpenMC executable
    cwd : str, optional
        Path to working directory to run in

    Raises
    ------
    subprocess.CalledProcessError
        If the `openmc` executable returns a non-zero status

    """
    _run([openmc_exec, '-p'], output, cwd)


def plot_inline(plots, openmc_exec='openmc', cwd='.', convert_exec='convert'):
    """Display plots inline in a Jupyter notebook.

    This function requires that you have a program installed to convert PPM
    files to PNG files. Typically, that would be `ImageMagick
    <https://www.imagemagick.org>`_ which includes a `convert` command.

    Parameters
    ----------
    plots : Iterable of openmc.Plot
        Plots to display
    openmc_exec : str
        Path to OpenMC executable
    cwd : str, optional
        Path to working directory to run in
    convert_exec : str, optional
        Command that can convert PPM files into PNG files

    Raises
    ------
    subprocess.CalledProcessError
        If the `openmc` executable returns a non-zero status

    """
    from IPython.display import Image, display

    if not isinstance(plots, Iterable):
        plots = [plots]

    # Create plots.xml
    openmc.Plots(plots).export_to_xml()

    # Run OpenMC in geometry plotting mode
    plot_geometry(False, openmc_exec, cwd)

    images = []
    if plots is not None:
        for p in plots:
            if p.filename is not None:
                ppm_file = '{}.ppm'.format(p.filename)
            else:
                ppm_file = 'plot_{}.ppm'.format(p.id)
            png_file = ppm_file.replace('.ppm', '.png')
            subprocess.check_call([convert_exec, ppm_file, png_file])
            images.append(Image(png_file))
        display(*images)


def calculate_volumes(threads=None, output=True, cwd='.',
                      openmc_exec='openmc', mpi_args=None):
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
        enabled, the default is implementation-dependent but is usually equal to
        the number of hardware threads available (or a value set by the
        :envvar:`OMP_NUM_THREADS` environment variable).
    output : bool, optional
        Capture OpenMC output from standard out
    openmc_exec : str, optional
        Path to OpenMC executable. Defaults to 'openmc'.
    mpi_args : list of str, optional
        MPI execute command and any additional MPI arguments to pass,
        e.g. ['mpiexec', '-n', '8'].
    cwd : str, optional
        Path to working directory to run in. Defaults to the current working
        directory.

    Raises
    ------
    subprocess.CalledProcessError
        If the `openmc` executable returns a non-zero status

    See Also
    --------
    openmc.VolumeCalculation

    """
    args = [openmc_exec, '--volume']

    if isinstance(threads, Integral) and threads > 0:
        args += ['-s', str(threads)]
    if mpi_args is not None:
        args = mpi_args + args

    _run(args, output, cwd)


def run(particles=None, threads=None, geometry_debug=False,
        restart_file=None, tracks=False, output=True, cwd='.',
        openmc_exec='openmc', mpi_args=None, event_based=False):
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
    output : bool
        Capture OpenMC output from standard out
    cwd : str, optional
        Path to working directory to run in. Defaults to the current working
        directory.
    openmc_exec : str, optional
        Path to OpenMC executable. Defaults to 'openmc'.
    mpi_args : list of str, optional
        MPI execute command and any additional MPI arguments to pass,
        e.g. ['mpiexec', '-n', '8'].
    event_based : bool, optional
        Turns on event-based parallelism, instead of default history-based

    Raises
    ------
    subprocess.CalledProcessError
        If the `openmc` executable returns a non-zero status

    """
    args = [openmc_exec]

    if isinstance(particles, Integral) and particles > 0:
        args += ['-n', str(particles)]

    if isinstance(threads, Integral) and threads > 0:
        args += ['-s', str(threads)]

    if geometry_debug:
        args.append('-g')

    if event_based:
        args.append('-e')

    if isinstance(restart_file, str):
        args += ['-r', restart_file]

    if tracks:
        args.append('-t')

    if mpi_args is not None:
        args = mpi_args + args

    _run(args, output, cwd)
