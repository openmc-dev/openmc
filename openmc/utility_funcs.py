from contextlib import contextmanager
import os
from pathlib import Path
from tempfile import TemporaryDirectory

import h5py

import openmc
from .checkvalue import PathLike


@contextmanager
def change_directory(working_dir: PathLike | None = None, *, tmpdir: bool = False):
    """Context manager for executing in a provided working directory

    Parameters
    ----------
    working_dir : path-like
        Directory to switch to.
    tmpdir : bool
        Whether to use a temporary directory instead of a specific working directory

    """
    orig_dir = Path.cwd()

    # Set up temporary directory if requested
    if tmpdir:
        tmp = TemporaryDirectory()
        working_dir = tmp.name
    elif working_dir is None:
        raise ValueError('Must pass working_dir argument or specify tmpdir=True.')

    working_dir = Path(working_dir)
    working_dir.mkdir(parents=True, exist_ok=True)
    os.chdir(working_dir)
    try:
        yield
    finally:
        os.chdir(orig_dir)
        if tmpdir:
            tmp.cleanup()


def input_path(filename: PathLike) -> Path:
    """Return a path object for an input file based on global configuration

    Parameters
    ----------
    filename : PathLike
        Path to input file

    Returns
    -------
    pathlib.Path
        Path object

    """
    if openmc.config['resolve_paths']:
        return Path(filename).resolve()
    else:
        return Path(filename)


@contextmanager
def h5py_file_or_group(group_or_filename: PathLike | h5py.Group, *args, **kwargs):
    """Context manager for opening an HDF5 file or using an existing group

    Parameters
    ----------
    group_or_filename : path-like or h5py.Group
        Path to HDF5 file, or group from an existing HDF5 file

    """
    if isinstance(group_or_filename, h5py.Group):
        yield group_or_filename
    else:
        with h5py.File(group_or_filename, *args, **kwargs) as f:
            yield f
