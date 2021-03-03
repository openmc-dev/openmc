from contextlib import contextmanager
import os
import tempfile


@contextmanager
def cdtemp():
    """Context manager to change to/return from a tmpdir."""
    with tempfile.TemporaryDirectory() as tmpdir:
        cwd = os.getcwd()
        try:
            os.chdir(tmpdir)
            yield
        finally:
            os.chdir(cwd)
