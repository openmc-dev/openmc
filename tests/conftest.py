import os
import pytest
import openmc

from tests.regression_tests import config as regression_config


def pytest_addoption(parser):
    parser.addoption('--exe')
    parser.addoption('--mpi', action='store_true')
    parser.addoption('--mpiexec')
    parser.addoption('--mpi-np')
    parser.addoption('--update', action='store_true')
    parser.addoption('--build-inputs', action='store_true')
    parser.addoption('--event', action='store_true')


def pytest_configure(config):
    opts = ['exe', 'mpi', 'mpiexec', 'mpi_np', 'update', 'build_inputs', 'event']
    for opt in opts:
        if config.getoption(opt) is not None:
            regression_config[opt] = config.getoption(opt)


@pytest.fixture
def run_in_tmpdir(tmpdir):
    orig = tmpdir.chdir()
    try:
        yield
    finally:
        orig.chdir()

@pytest.fixture(scope="module")
def endf_data():
    return os.environ['OPENMC_ENDF_DATA']

@pytest.fixture(scope='session', autouse=True)
def resolve_paths():
    with openmc.config.patch('resolve_paths', False):
        yield


@pytest.fixture(scope='session', autouse=True)
def disable_depletion_multiprocessing_under_mpi():
    """Fork-based depletion multiprocessing may deadlock if MPI is active."""
    if not regression_config['mpi']:
        yield
        return

    from openmc.deplete import pool

    original_setting = pool.USE_MULTIPROCESSING
    pool.USE_MULTIPROCESSING = False
    try:
        yield
    finally:
        pool.USE_MULTIPROCESSING = original_setting
