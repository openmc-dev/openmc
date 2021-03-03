import pytest

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
