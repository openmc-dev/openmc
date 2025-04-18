import os
import sys
import pytest
import platform
from unittest import mock
from pathlib import Path
from openmc import paths


def test_get_paths_non_recursive():
    """Test get_paths with non-recursive file search."""
    result = paths.get_paths("include", "*", recursive=False)
    assert isinstance(result, list)


def test_get_paths_recursive():
    """Test get_paths with recursive file search."""
    result = paths.get_paths("include", "*", recursive=True)
    assert isinstance(result, list)


def test_get_include_path():
    """Test retrieval of include files and their path."""
    include, include_path = paths.get_include_path()
    assert isinstance(include, list)
    assert isinstance(include_path, list)


def test_get_core_libraries():
    """Test retrieval of core library files and their path."""
    lib, lib_path = paths.get_core_libraries()
    assert isinstance(lib, list)
    assert isinstance(lib_path, list)


def test_get_extra_libraries_existing_dir(tmp_path):
    """Test get_extra_libraries when expected directories exist."""
    # Choose directory based on platform
    if sys.platform == "darwin":
        dylibs_dir = tmp_path / ".dylibs"
        dylibs_dir.mkdir()
        dummy_file = dylibs_dir / "libdummy.dylib"
    else:
        libs_dir = tmp_path.parent / "openmc.libs"
        libs_dir.mkdir()
        dummy_file = libs_dir / "libdummy.so"

    dummy_file.write_text("mock")

    with mock.patch.object(paths, "__path__", [str(tmp_path)]):
        extra_libs, extra_lib_path = paths.get_extra_libraries()
        assert isinstance(extra_libs, list)
        assert isinstance(extra_lib_path, str)
        assert dummy_file in map(Path, extra_libs)
        assert extra_lib_path == str(libs_dir if sys.platform != "darwin" else dylibs_dir)


def test_get_extra_libraries_missing_dir():
    """Test get_extra_libraries when expected directories are missing."""
    with mock.patch.object(paths, "__path__", ["/nonexistent/fakepath"]):
        extra_libs, extra_lib_path = paths.get_extra_libraries()
        assert extra_libs == []
        assert extra_lib_path == []