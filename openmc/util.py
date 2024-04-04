from contextlib import contextmanager
import os
from pathlib import Path


@contextmanager
def _change_directory(working_dir):
    """A context manager for executing in a provided working directory"""
    start_dir = Path.cwd()
    Path.mkdir(working_dir, parents=True, exist_ok=True)
    os.chdir(working_dir)
    try:
        yield
    finally:
        os.chdir(start_dir)