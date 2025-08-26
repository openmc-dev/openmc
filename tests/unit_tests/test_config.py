from collections.abc import Mapping
import os
from pathlib import Path

import openmc
from openmc.config import _default_config
from openmc.data import decay
import pytest


@pytest.fixture(autouse=True, scope='function')
def reset_config_and_env():
    """A fixture to ensure each test has a clean config, env, and CWD."""
    original_env = dict(os.environ)
    original_cwd = os.getcwd()

    # Reset environment variables that affect config
    for key in ['OPENMC_CROSS_SECTIONS', 'OPENMC_MG_CROSS_SECTIONS', 'OPENMC_CHAIN_FILE']:
        if key in os.environ:
            del os.environ[key]

    # Re-initialize the global config object
    openmc.config = _default_config()

    try:
        yield
    finally:
        # Restore environment and CWD
        os.environ.clear()
        os.environ.update(original_env)
        os.chdir(original_cwd)

        # Restore config one last time for safety between modules
        openmc.config = _default_config()


def test_config_basics():
    assert isinstance(openmc.config, Mapping)
    with pytest.warns(UserWarning):
        openmc.config['cross_sections'] = '/path/to/cross_sections.xml'
    del openmc.config['cross_sections']
    assert 'cross_sections' not in openmc.config
    assert 'OPENMC_CROSS_SECTIONS' not in os.environ
    with pytest.raises(KeyError, match="Unrecognized config key: nuke"):
        openmc.config['nuke'] = '/like/to/eat/bacon'
    with pytest.raises(TypeError):
        openmc.config['resolve_paths'] = 'not a bool'


def test_config_path_resolution(tmp_path):
    """Test path resolution logic."""
    os.chdir(tmp_path)
    relative_path = Path("some/file.xml")
    absolute_path = relative_path.resolve()

    # Test with resolve_paths = True (default)
    with pytest.warns(UserWarning):
        openmc.config['cross_sections'] = relative_path
    assert openmc.config['cross_sections'] == absolute_path
    assert openmc.config['cross_sections'].is_absolute()

    # Test with resolve_paths = False
    with openmc.config.patch('resolve_paths', False):
        with pytest.warns(UserWarning):
            openmc.config['chain_file'] = relative_path
        assert openmc.config['chain_file'] == relative_path
        assert not openmc.config['chain_file'].is_absolute()

    assert openmc.config['resolve_paths'] is True


def test_config_patch(tmp_path):
    file_a = tmp_path / "a.xml"; file_a.touch()
    file_b = tmp_path / "b.xml"; file_b.touch()
    openmc.config['cross_sections'] = file_a
    with openmc.config.patch('cross_sections', file_b):
        assert openmc.config['cross_sections'] == file_b.resolve()
    assert openmc.config['cross_sections'] == file_a.resolve()

def test_config_set_envvar(tmp_path):
    """Test that setting config also sets environment variables correctly."""
    os.chdir(tmp_path)
    relative_path = Path("relative.xml")
    with pytest.warns(UserWarning):
        openmc.config['cross_sections'] = relative_path
    expected_path = str(relative_path.resolve())
    assert os.environ['OPENMC_CROSS_SECTIONS'] == expected_path


def test_config_warning_nonexistent_path(tmp_path):
    """Test that a warning is issued for a path that does not exist."""
    bad_path = tmp_path / "a/path/that/does/not/exist.xml"
    with pytest.warns(UserWarning, match=f"Path '{bad_path}' does not exist."):
        openmc.config['chain_file'] = bad_path


def test_config_chain_side_effect(tmp_path):
    """Test that modifying chain_file clears decay data caches."""
    chain_file = tmp_path / "chain.xml"; chain_file.touch()
    decay._DECAY_ENERGY['U235'] = (1.0, 2.0)
    decay._DECAY_PHOTON_ENERGY['PU239'] = {}
    openmc.config['chain_file'] = chain_file
    assert not decay._DECAY_ENERGY and not decay._DECAY_PHOTON_ENERGY
