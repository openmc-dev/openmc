from collections.abc import Mapping
import os

import openmc
import pytest


@pytest.fixture(autouse=True, scope='module')
def reset_config():
    config = dict(openmc.config)
    try:
        yield
    finally:
        openmc.config.clear()
        openmc.config.update(config)


def test_config_basics():
    assert isinstance(openmc.config, Mapping)
    for key, value in openmc.config.items():
        assert isinstance(key, str)
        if key == 'resolve_paths':
            assert isinstance(value, bool)
        else:
            assert isinstance(value, os.PathLike)

    # Set and delete
    openmc.config['cross_sections'] = '/path/to/cross_sections.xml'
    del openmc.config['cross_sections']
    assert 'cross_sections' not in openmc.config
    assert 'OPENMC_CROSS_SECTIONS' not in os.environ

    # Can't use any key
    with pytest.raises(KeyError):
        openmc.config['🐖'] = '/like/to/eat/bacon'


def test_config_patch():
    openmc.config['cross_sections'] = '/path/to/cross_sections.xml'
    with openmc.config.patch('cross_sections', '/path/to/other.xml'):
        assert str(openmc.config['cross_sections']) == '/path/to/other.xml'
    assert str(openmc.config['cross_sections']) == '/path/to/cross_sections.xml'


def test_config_set_envvar():
    openmc.config['cross_sections'] = '/path/to/cross_sections.xml'
    assert os.environ['OPENMC_CROSS_SECTIONS'] == '/path/to/cross_sections.xml'

    openmc.config['mg_cross_sections'] = '/path/to/mg_cross_sections.h5'
    assert os.environ['OPENMC_MG_CROSS_SECTIONS'] == '/path/to/mg_cross_sections.h5'

    openmc.config['chain_file'] = '/path/to/chain_file.xml'
    assert os.environ['OPENMC_CHAIN_FILE'] == '/path/to/chain_file.xml'
