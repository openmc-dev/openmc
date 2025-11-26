from contextlib import contextmanager
from ctypes import (c_bool, c_int, c_int32, c_int64, c_double, c_char_p,
                    c_char, POINTER, Structure, c_void_p, create_string_buffer,
                    c_uint64, c_size_t)
import sys
import os
from pathlib import Path
from random import getrandbits
from tempfile import TemporaryDirectory
import traceback as tb

import numpy as np
from numpy.ctypeslib import as_array

from . import _dll
from .error import _error_handler
from ..mpi import comm
from openmc.checkvalue import PathLike
import openmc.lib
import openmc


class _SourceSite(Structure):
    _fields_ = [('r', c_double*3),
                ('u', c_double*3),
                ('E', c_double),
                ('E_parent', c_double),
                ('time', c_double),
                ('wgt', c_double),
                ('delayed_group', c_int),
                ('fission_nuclide', c_int),
                ('surf_id', c_int),
                ('particle', c_int),
                ('parent_nuclide', c_int),
                ('parent_id', c_int64),
                ('progeny_id', c_int64)]


# Define input type for numpy arrays that will be passed into C++ functions
# Must be an int or double array, with single dimension that is contiguous
_array_1d_int = np.ctypeslib.ndpointer(dtype=np.int32, ndim=1,
                                       flags='CONTIGUOUS')
_array_1d_dble = np.ctypeslib.ndpointer(dtype=np.double, ndim=1,
                                        flags='CONTIGUOUS')

_dll.openmc_calculate_volumes.restype = c_int
_dll.openmc_calculate_volumes.errcheck = _error_handler
_dll.openmc_cmfd_reweight.argtypes = c_bool, _array_1d_dble
_dll.openmc_cmfd_reweight.restype = None
_dll.openmc_finalize.restype = c_int
_dll.openmc_finalize.errcheck = _error_handler
_dll.openmc_find_cell.argtypes = [POINTER(c_double*3), POINTER(c_int32),
                                  POINTER(c_int32)]
_dll.openmc_find_cell.restype = c_int
_dll.openmc_find_cell.errcheck = _error_handler
_dll.openmc_hard_reset.restype = c_int
_dll.openmc_hard_reset.errcheck = _error_handler
_dll.openmc_init.argtypes = [c_int, POINTER(POINTER(c_char)), c_void_p]
_dll.openmc_init.restype = c_int
_dll.openmc_init.errcheck = _error_handler
_dll.openmc_get_keff.argtypes = [POINTER(c_double*2)]
_dll.openmc_get_keff.restype = c_int
_dll.openmc_get_keff.errcheck = _error_handler
_dll.openmc_initialize_mesh_egrid.argtypes = [
    c_int, _array_1d_int, c_double
]
_dll.openmc_initialize_mesh_egrid.restype = None
_init_linsolver_argtypes = [_array_1d_int, c_int, _array_1d_int, c_int, c_int,
                            c_double, _array_1d_int, c_bool]
_dll.openmc_initialize_linsolver.argtypes = _init_linsolver_argtypes
_dll.openmc_initialize_linsolver.restype = None
_dll.openmc_is_statepoint_batch.restype = c_bool
_dll.openmc_master.restype = c_bool
_dll.openmc_next_batch.argtypes = [POINTER(c_int)]
_dll.openmc_next_batch.restype = c_int
_dll.openmc_next_batch.errcheck = _error_handler
_dll.openmc_plot_geometry.restype = c_int
_dll.openmc_plot_geometry.errcheck = _error_handler
_dll.openmc_properties_export.argtypes = [c_char_p]
_dll.openmc_properties_export.restype = c_int
_dll.openmc_properties_export.errcheck = _error_handler
_dll.openmc_properties_import.argtypes = [c_char_p]
_dll.openmc_properties_import.restype = c_int
_dll.openmc_properties_import.errcheck = _error_handler
_dll.openmc_run.restype = c_int
_dll.openmc_run.errcheck = _error_handler
_dll.openmc_reset.restype = c_int
_dll.openmc_reset.errcheck = _error_handler
_dll.openmc_reset_timers.restype = c_int
_dll.openmc_reset_timers.errcheck = _error_handler
_run_linsolver_argtypes = [_array_1d_dble, _array_1d_dble, _array_1d_dble,
                           c_double]
_dll.openmc_run_linsolver.argtypes = _run_linsolver_argtypes
_dll.openmc_run_linsolver.restype = c_int
_dll.openmc_source_bank.argtypes = [POINTER(POINTER(_SourceSite)), POINTER(c_int64)]
_dll.openmc_source_bank.restype = c_int
_dll.openmc_source_bank.errcheck = _error_handler
_dll.openmc_simulation_init.restype = c_int
_dll.openmc_simulation_init.errcheck = _error_handler
_dll.openmc_simulation_finalize.restype = c_int
_dll.openmc_simulation_finalize.errcheck = _error_handler
_dll.openmc_statepoint_write.argtypes = [c_char_p, POINTER(c_bool)]
_dll.openmc_statepoint_write.restype = c_int
_dll.openmc_statepoint_write.errcheck = _error_handler
_dll.openmc_statepoint_load.argtypes = [c_char_p]
_dll.openmc_statepoint_load.restype = c_int
_dll.openmc_statepoint_load.errcheck = _error_handler
_dll.openmc_statepoint_write.restype = c_int
_dll.openmc_statepoint_write.errcheck = _error_handler
_dll.openmc_global_bounding_box.argtypes = [POINTER(c_double),
                                            POINTER(c_double)]
_dll.openmc_global_bounding_box.restype = c_int
_dll.openmc_global_bounding_box.errcheck = _error_handler
_dll.openmc_sample_external_source.argtypes = [c_size_t, POINTER(c_uint64), POINTER(_SourceSite)]
_dll.openmc_sample_external_source.restype = c_int
_dll.openmc_sample_external_source.errcheck = _error_handler

def global_bounding_box():
    """Calculate a global bounding box for the model"""
    inf = sys.float_info.max
    llc = np.zeros(3)
    urc = np.zeros(3)
    _dll.openmc_global_bounding_box(llc.ctypes.data_as(POINTER(c_double)),
                                    urc.ctypes.data_as(POINTER(c_double)))
    llc[llc == inf] = np.inf
    urc[urc == inf] = np.inf
    llc[llc == -inf] = -np.inf
    urc[urc == -inf] = -np.inf

    return llc, urc


def calculate_volumes(output=True):
    """Run stochastic volume calculation

    .. versionchanged:: 0.13.0
        The *output* argument was added.

    Parameters
    ----------
    output : bool, optional
        Whether or not to show output. Defaults to showing output

    """

    with quiet_dll(output):
        _dll.openmc_calculate_volumes()


def current_batch():
    """Return the current batch of the simulation.

    Returns
    -------
    int
        Current batch of the simulation

    """
    return c_int.in_dll(_dll, 'current_batch').value


def export_properties(filename=None, output=True):
    """Export physical properties.

    .. versionadded:: 0.13.0

    Parameters
    ----------
    filename : str or None
        Filename to export properties to (defaults to "properties.h5")
    output : bool, optional
        Whether or not to show output. Defaults to showing output

    See Also
    --------
    openmc.lib.import_properties

    """
    if filename is not None:
        filename = c_char_p(filename.encode())

    with quiet_dll(output):
        _dll.openmc_properties_export(filename)


def export_weight_windows(filename="weight_windows.h5", output=True):
    """Export weight windows.

    .. versionadded:: 0.14.0

    Parameters
    ----------
    filename : PathLike or None
        Filename to export weight windows to
    output : bool, optional
        Whether or not to show output.

    See Also
    --------
    openmc.lib.import_weight_windows

    """
    if filename is not None:
        filename = c_char_p(str(filename).encode())

    with quiet_dll(output):
        _dll.openmc_weight_windows_export(filename)


def import_weight_windows(filename='weight_windows.h5', output=True):
    """Import weight windows.

    .. versionadded:: 0.14.0

    Parameters
    ----------
    filename : PathLike or None
        Filename to import weight windows from
    output : bool, optional
        Whether or not to show output.

    See Also
    --------
    openmc.lib.export_weight_windows

    """
    if filename is not None:
        filename = c_char_p(str(filename).encode())

    with quiet_dll(output):
        _dll.openmc_weight_windows_import(filename)


def finalize():
    """Finalize simulation and free memory"""
    _dll.openmc_finalize()
    openmc.lib.is_initialized = False


def find_cell(xyz):
    """Find the cell at a given point

    Parameters
    ----------
    xyz : iterable of float
        Cartesian coordinates of position

    Returns
    -------
    openmc.lib.Cell
        Cell containing the point
    int
        If the cell at the given point is repeated in the geometry, this
        indicates which instance it is, i.e., 0 would be the first instance.

    """
    index = c_int32()
    instance = c_int32()
    _dll.openmc_find_cell((c_double*3)(*xyz), index, instance)
    return openmc.lib.Cell(index=index.value), instance.value


def find_material(xyz):
    """Find the material at a given point

    Parameters
    ----------
    xyz : iterable of float
        Cartesian coordinates of position

    Returns
    -------
    openmc.lib.Material or None
        Material containing the point, or None is no material is found

    """
    index = c_int32()
    instance = c_int32()
    _dll.openmc_find_cell((c_double*3)(*xyz), index, instance)

    mats = openmc.lib.Cell(index=index.value).fill
    if isinstance(mats, (openmc.lib.Material, type(None))):
        return mats
    else:
        return mats[instance.value]


def hard_reset():
    """Reset tallies, timers, and pseudo-random number generator state."""
    _dll.openmc_hard_reset()


def import_properties(filename):
    """Import physical properties.

    .. versionadded:: 0.13.0

    Parameters
    ----------
    filename : str
        Filename to import properties from

    See Also
    --------
    openmc.lib.export_properties

    """
    _dll.openmc_properties_import(filename.encode())


def init(args=None, intracomm=None, output=True):
    """Initialize OpenMC

    .. versionchanged:: 0.13.0
        The *output* argument was added.

    Parameters
    ----------
    args : list of str, optional
        Command-line arguments
    intracomm : mpi4py.MPI.Intracomm or None, optional
        MPI intracommunicator
    output : bool, optional
        Whether or not to show output. Defaults to showing output

    """
    if args is not None:
        args = ['openmc'] + list(args)
    else:
        args = ['openmc']

    argc = len(args)
    # Create the argv array. Note that it is actually expected to be of
    # length argc + 1 with the final item being a null pointer.
    argv = (POINTER(c_char) * (argc + 1))()
    for i, arg in enumerate(args):
        argv[i] = create_string_buffer(arg.encode())

    if intracomm is not None:
        # If an mpi4py communicator was passed, convert it to void* to be passed
        # to openmc_init
        try:
            from mpi4py import MPI
        except ImportError:
            intracomm = None
        else:
            address = MPI._addressof(intracomm)
            intracomm = c_void_p(address)

    with quiet_dll(output):
        _dll.openmc_init(argc, argv, intracomm)
    openmc.lib.is_initialized = True


def is_statepoint_batch():
    """Return whether statepoint will be written in current batch or not.

    Returns
    -------
    bool
        Whether is statepoint batch or not

    """
    return _dll.openmc_is_statepoint_batch()


def iter_batches():
    """Iterator over batches.

    This function returns a generator-iterator that allows Python code to be run
    between batches in an OpenMC simulation. It should be used in conjunction
    with :func:`openmc.lib.simulation_init` and
    :func:`openmc.lib.simulation_finalize`. For example:

    .. code-block:: Python

        with openmc.lib.run_in_memory():
            openmc.lib.simulation_init()
            for _ in openmc.lib.iter_batches():
                # Look at convergence of tallies, for example
                ...
            openmc.lib.simulation_finalize()

    See Also
    --------
    openmc.lib.next_batch

    """
    while True:
        # Run next batch
        status = next_batch()

        # Provide opportunity for user to perform action between batches
        yield

        # End the iteration
        if status != 0:
            break


def keff():
    """Return the calculated k-eigenvalue and its standard deviation.

    Returns
    -------
    tuple
        Mean k-eigenvalue and standard deviation of the mean

    """
    k = (c_double*2)()
    _dll.openmc_get_keff(k)
    return tuple(k)


def master():
    """Return whether processor is master processor or not.

    Returns
    -------
    bool
        Whether is master processor or not

    """
    return _dll.openmc_master()


def next_batch():
    """Run next batch.

    Returns
    -------
    int
        Status after running a batch (0=normal, 1=reached maximum number of
        batches, 2=tally triggers reached)

    """
    status = c_int()
    _dll.openmc_next_batch(status)
    return status.value


def plot_geometry(output=True):
    """Plot geometry

    .. versionchanged:: 0.13.0
        The *output* argument was added.

    Parameters
    ----------
    output : bool, optional
        Whether or not to show output. Defaults to showing output
    """

    with quiet_dll(output):
        _dll.openmc_plot_geometry()


def reset():
    """Reset tally results"""
    _dll.openmc_reset()


def reset_timers():
    """Reset timers."""
    _dll.openmc_reset_timers()


def run(output=True):
    """Run simulation

    .. versionchanged:: 0.13.0
        The *output* argument was added.

    Parameters
    ----------
    output : bool, optional
        Whether or not to show output. Defaults to showing output
    """

    with quiet_dll(output):
        _dll.openmc_run()


def sample_external_source(
        n_samples: int = 1000,
        prn_seed: int | None = None
) -> openmc.ParticleList:
    """Sample external source and return source particles.

    .. versionadded:: 0.13.1

    Parameters
    ----------
    n_samples : int
        Number of samples
    prn_seed : int
        Pseudorandom number generator (PRNG) seed; if None, one will be
        generated randomly.

    Returns
    -------
    openmc.ParticleList
        List of sampled source particles

    """
    if n_samples <= 0:
        raise ValueError("Number of samples must be positive")
    if prn_seed is None:
        prn_seed = getrandbits(63)

    # Call into C API to sample source
    sites_array = (_SourceSite * n_samples)()
    _dll.openmc_sample_external_source(c_size_t(n_samples), c_uint64(prn_seed), sites_array)

    # Convert to list of SourceParticle and return
    return openmc.ParticleList([openmc.SourceParticle(
            r=site.r, u=site.u, E=site.E, time=site.time, wgt=site.wgt,
            delayed_group=site.delayed_group, surf_id=site.surf_id,
            particle=openmc.ParticleType(site.particle)
        )
        for site in sites_array
    ])


def simulation_init():
    """Initialize simulation"""
    _dll.openmc_simulation_init()


def simulation_finalize():
    """Finalize simulation"""
    _dll.openmc_simulation_finalize()


def source_bank():
    """Return source bank as NumPy array

    Returns
    -------
    numpy.ndarray
        Source sites

    """
    # Get pointer to source bank
    ptr = POINTER(_SourceSite)()
    n = c_int64()
    _dll.openmc_source_bank(ptr, n)

    try:
        # Convert to numpy array with appropriate datatype
        bank_dtype = np.dtype(_SourceSite)
        return as_array(ptr, (n.value,)).view(bank_dtype)

    except ValueError as err:
        # If a known numpy error was raised (github.com/numpy/numpy/issues
        # /14214), re-raise with a more helpful error message.
        if len(err.args) == 0:
            raise err
        if err.args[0].startswith('invalid shape in fixed-type tuple'):
            raise ValueError('The source bank is too large to access via '
                'openmc.lib with this version of numpy.  Use a different '
                'version of numpy or reduce the bank size (fewer particles '
                'per MPI process) so that it is smaller than 2 GB.') from err
        else:
            raise err


def statepoint_write(filename=None, write_source=True):
    """Write a statepoint file.

    Parameters
    ----------
    filename : str or None
        Path to the statepoint to write. If None is passed, a default name that
        contains the current batch will be written.
    write_source : bool
        Whether or not to include the source bank in the statepoint.

    """
    if filename is not None:
        filename = c_char_p(filename.encode())
    _dll.openmc_statepoint_write(filename, c_bool(write_source))


def statepoint_load(filename: PathLike):
    """Load a statepoint file.

    Parameters
    ----------
    filename : path-like
        Path to the statepoint to load.

    """
    filename = c_char_p(str(filename).encode())
    _dll.openmc_statepoint_load(filename)


@contextmanager
def run_in_memory(**kwargs):
    """Provides context manager for calling OpenMC shared library functions.

    This function is intended to be used in a 'with' statement and ensures that
    OpenMC is properly initialized/finalized. At the completion of the 'with'
    block, all memory that was allocated during the block is freed. For
    example::

        with openmc.lib.run_in_memory():
            for i in range(n_iters):
                openmc.lib.reset()
                do_stuff()
                openmc.lib.run()

    Parameters
    ----------
    **kwargs
        All keyword arguments are passed to :func:`init`.

    """
    init(**kwargs)
    try:
        yield
    finally:
        finalize()


class TemporarySession:
    """Context manager for running via openmc.lib in a temporary directory.

    This class is useful for accessing functionality from openmc.lib without
    polluting your current working directory with OpenMC files. It is used
    internally as a persistent session to avoid loading cross sections multiple
    times.

    Parameters
    ----------
    model : openmc.Model, optional
        OpenMC model to use for the session. If None, a minimal working model is
        created.
    cwd : PathLike, optional
        Working directory in which to run OpenMC. If None, a temporary directory
        is created and deleted automatically.
    **init_kwargs
        Keyword arguments to pass to :func:`openmc.lib.init`.

    Attributes
    ----------
    model : openmc.Model
        The OpenMC model used for the session.
    comm : mpi4py.MPI.Intracomm
        The MPI intracommunicator used for the session.

    """
    def __init__(self, model=None, cwd=None, **init_kwargs):
        self.init_kwargs = dict(init_kwargs)
        self.cwd = cwd
        if model is None:
            surf = openmc.Sphere(boundary_type="vacuum")
            cell = openmc.Cell(region=-surf)
            model = openmc.Model()
            model.geometry = openmc.Geometry([cell])
            model.settings = openmc.Settings(
                particles=1, batches=1, output={'summary': False})
        self.model = model

        # Determine MPI intercommunicator
        self.init_kwargs.setdefault('intracomm', comm)
        self.comm = self.init_kwargs['intracomm']

    def __enter__(self):
        """Initialize the OpenMC library in a temporary directory."""
        # If already initialized, the context manager is a no-op
        self.already_initialized = openmc.lib.is_initialized
        if self.already_initialized:
            return self

        # Store original working directory
        self.orig_dir = Path.cwd()

        if self.cwd is None:
            # Set up temporary directory on rank 0
            if self.comm.rank == 0:
                self._tmp_dir = TemporaryDirectory()
                self.cwd = self._tmp_dir.name

            # Broadcast the path so that all ranks use the same directory
            self.cwd = self.comm.bcast(self.cwd)

        # Create and change to specified directory
        self.cwd = Path(self.cwd)
        self.cwd.mkdir(parents=True, exist_ok=True)
        os.chdir(self.cwd)

        # Export model on first rank and initialize OpenMC
        if self.comm.rank == 0:
            self.model.export_to_model_xml()
        self.comm.barrier()
        openmc.lib.init(**self.init_kwargs)

        return self

    def __exit__(self, exc_type, exc_value, traceback):
        """Finalize the OpenMC library and clean up temporary directory."""
        if self.already_initialized:
            return

        # If an exception occurred, abort all ranks immediately
        if exc_type is not None:
            # Print exception info on the rank that failed
            tb.print_exception(exc_type, exc_value, traceback)
            sys.stdout.flush()

            # Abort all MPI processes
            self.comm.Abort(1)

        try:
            finalize()
        finally:
            os.chdir(self.orig_dir)

            # Make sure all ranks have finalized before deleting temporary dir
            self.comm.barrier()
            if hasattr(self, '_tmp_dir'):
                self._tmp_dir.cleanup()


class _DLLGlobal:
    """Data descriptor that exposes global variables from libopenmc."""
    def __init__(self, ctype, name):
        self.ctype = ctype
        self.name = name

    def __get__(self, instance, owner):
        return self.ctype.in_dll(_dll, self.name).value

    def __set__(self, instance, value):
        self.ctype.in_dll(_dll, self.name).value = value


class _FortranObject:
    def __repr__(self):
        return f"<{type(self).__name__}(index={self._index})>"


class _FortranObjectWithID(_FortranObject):
    def __init__(self, uid=None, new=True, index=None):
        # Creating the object has already been handled by __new__. In the
        # initializer, all we do is make sure that the object returned has an ID
        # assigned. If the array index of the object is out of bounds, an
        # OutOfBoundsError will be raised here by virtue of referencing self.id
        self.id

    def __repr__(self):
        return f"<{type(self).__name__}(id={self.id})>"


@contextmanager
def quiet_dll(output=True):
    """This context manager allows us to suppress standard output from DLLs

    Parameters
    ----------
    output : bool
        Denotes whether the output should be displayed (True) or not (False)

    .. versionadded:: 0.13.0

    """

    # This contextmanager is modified from that provided here:
    # https://stackoverflow.com/a/14797594

    if output:
        yield
    else:
        sys.stdout.flush()
        # Save the initial file descriptor states
        initial_stdout = sys.stdout
        initial_stdout_fno = os.dup(sys.stdout.fileno())
        # Get a garbage descriptor so we can throw away output
        devnull = os.open(os.devnull, os.O_WRONLY)

        # Get the current stdout stream and make a duplicate of it
        new_stdout = os.dup(1)
        # Copy the garbage output to the stdout stream
        os.dup2(devnull, 1)
        os.close(devnull)
        # Now point stdout to the re-defined stdout
        sys.stdout = os.fdopen(new_stdout, 'w')

        try:
            yield
        finally:
            # Now we just clean up after ourselves and reset the streams
            sys.stdout = initial_stdout
            sys.stdout.flush()
            os.dup2(initial_stdout_fno, 1)
