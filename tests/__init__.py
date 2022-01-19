from contextlib import contextmanager
import os
from shutil import copy
import tempfile


@contextmanager
def cdtemp(files=None):
    """Context manager to change to/return from a tmpdir.

    Parameters
    ----------
    files : Iterable of str or Path-like
        Set of files to copy into the temporary directory
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        cwd = os.getcwd()
        if files:
            for file in files:
                copy(file, tmpdir, follow_symlinks=True)
        try:
            os.chdir(tmpdir)
            yield
        finally:
            os.chdir(cwd)
