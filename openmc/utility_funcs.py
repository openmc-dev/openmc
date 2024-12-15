from contextlib import contextmanager
import os
from pathlib import Path
from tempfile import TemporaryDirectory

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
        raise ValueError("Must pass working_dir argument or specify tmpdir=True.")

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
    if openmc.config["resolve_paths"]:
        return Path(filename).resolve()
    else:
        return Path(filename)
