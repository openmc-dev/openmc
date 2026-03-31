import os
import hashlib
import pytest
import openmc
import openmc.lib

from tests.regression_tests import config as regression_config

# MD5 hash of the official NNDC HDF5 cross_sections.xml file.
# Generated via: md5sum /path/to/nndc_hdf5/cross_sections.xml
_NNDC_XS_MD5 = "2d00773012eda670bc9f95d96a31c989"

# Collected during pytest_configure, displayed at start and end of session
_environment_warnings = []


def _check_build_environment():
    """Check STRICT_FP and cross section data, collecting any warnings."""
    if not openmc.lib._strict_fp_enabled():
        _environment_warnings.append(
            "OpenMC was NOT built with -DOPENMC_ENABLE_STRICT_FP=on. "
            "Regression test results may not match reference values due to "
            "compiler floating-point optimizations. Rebuild with "
            "-DOPENMC_ENABLE_STRICT_FP=on for reproducible results."
        )

    xs_path = os.environ.get("OPENMC_CROSS_SECTIONS")
    if not xs_path:
        _environment_warnings.append(
            "OPENMC_CROSS_SECTIONS environment variable is not set. "
            "Regression tests require the NNDC HDF5 cross section data."
        )
    elif not os.path.isfile(xs_path):
        _environment_warnings.append(
            f"OPENMC_CROSS_SECTIONS ({xs_path}) is not a valid file path. "
            "Regression tests require the NNDC HDF5 cross section data."
        )
    else:
        with open(xs_path, "rb") as f:
            md5 = hashlib.md5(f.read()).hexdigest()
        if md5 != _NNDC_XS_MD5:
            _environment_warnings.append(
                f"OPENMC_CROSS_SECTIONS ({xs_path}) does not match the "
                "official NNDC HDF5 dataset. Regression tests expect the "
                "NNDC data; results may differ with other cross section "
                "libraries."
            )


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

    _check_build_environment()


def _print_environment_warnings(terminalreporter):
    """Print environment warnings as a visible section."""
    if _environment_warnings:
        terminalreporter.section("OpenMC Environment Warnings")
        for msg in _environment_warnings:
            terminalreporter.line(f"WARNING: {msg}", yellow=True)
        terminalreporter.line("")


def pytest_sessionstart(session):
    """Print environment warnings at the start of the test session."""
    _print_environment_warnings(session.config.pluginmanager.get_plugin(
        "terminalreporter"))


def pytest_terminal_summary(terminalreporter, exitstatus, config):
    """Reprint environment warnings at the end so they aren't missed."""
    _print_environment_warnings(terminalreporter)


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
