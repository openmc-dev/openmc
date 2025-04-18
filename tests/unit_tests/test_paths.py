import sys
import pytest
import importlib
from unittest import mock
from pathlib import Path
from openmc import paths


def test_openmc_core_base_path_nameerror(monkeypatch):
    """Test fallback logic when __path__ is not defined."""
    monkeypatch.setitem(sys.modules["openmc.paths"].__dict__, '__path__', None)

    with mock.patch("os.path.exists", return_value=False), \
         mock.patch("sysconfig.get_path", return_value="/mock/path"):
        with pytest.raises(ImportError, match="OpenMC is not installed"):
            importlib.reload(paths)

def test_openmc_core_base_path_importerror(monkeypatch):
    """Test ImportError raised when OpenMC is not installed and no core path is found."""
    monkeypatch.setitem(sys.modules['openmc.paths'].__dict__, '__path__', None)
    with mock.patch("os.path.exists", return_value=False), \
         mock.patch("sysconfig.get_path", return_value="/mock/path"):
        with pytest.raises(ImportError, match="OpenMC is not installed. Please run 'pip install openmc'."):
            importlib.reload(paths)

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

@pytest.mark.parametrize("platform_value, expected_dir_name, expected_ext", [
    ("darwin", ".dylibs", "dylib"),
    ("linux", "openmc.libs", "so"),
])
def test_get_extra_libraries_cross_platform(tmp_path, platform_value, expected_dir_name, expected_ext):
    """Simulate different platforms to test get_extra_libraries logic completely."""
    lib_dir = tmp_path / expected_dir_name
    lib_dir.mkdir()
    dummy_file = lib_dir / f"libdummy.{expected_ext}"
    dummy_file.write_text("mock")

    with mock.patch.object(paths, "__path__", [str(tmp_path)]):
        with mock.patch("sys.platform", new=platform_value):
            extra_libs, extra_lib_path = paths.get_extra_libraries()

            assert isinstance(extra_libs, list)
            if extra_lib_path:
                assert isinstance(extra_lib_path, str)
                assert dummy_file in map(Path, extra_libs)
                assert extra_lib_path == str(lib_dir)
            else:
                assert extra_libs == []

def test_get_extra_libraries_missing_dir():
    """Test get_extra_libraries when expected directories are missing."""
    with mock.patch.object(paths, "__path__", ["/nonexistent/fakepath"]):
        extra_libs, extra_lib_path = paths.get_extra_libraries()
        assert extra_libs == []
        assert extra_lib_path == []

def test_get_extra_libraries_no_match(tmp_path):
    """Ensure get_extra_libraries returns empty if no known lib dirs exist."""
    isolated_path = tmp_path / "nothing_here"
    isolated_path.mkdir()

    with mock.patch.object(paths, "__path__", [str(isolated_path)]):
        extra_libs, extra_lib_path = paths.get_extra_libraries()
        assert extra_libs == []
        assert extra_lib_path == []
