import os

import numpy as np
import pytest
import openmc.data


pytestmark = pytest.mark.skipif(
    'OPENMC_MULTIPOLE_LIBRARY' not in os.environ,
    reason='OPENMC_MULTIPOLE_LIBRARY environment variable must be set')


@pytest.fixture(scope='module')
def u235():
    directory = os.environ['OPENMC_MULTIPOLE_LIBRARY']
    filename = os.path.join(directory, '092235.h5')
    return openmc.data.WindowedMultipole.from_hdf5(filename)


@pytest.fixture(scope='module')
def u234():
    directory = os.environ['OPENMC_MULTIPOLE_LIBRARY']
    filename = os.path.join(directory, '092234.h5')
    return openmc.data.WindowedMultipole.from_hdf5(filename)


@pytest.fixture(scope='module')
def fe56():
    directory = os.environ['OPENMC_MULTIPOLE_LIBRARY']
    filename = os.path.join(directory, '026056.h5')
    return openmc.data.WindowedMultipole.from_hdf5(filename)


def test_evaluate_rm(u235):
    """Make sure a Reich-Moore multipole object can be called."""
    energies = [1e-3, 1.0, 10.0, 50.]
    total, absorption, fission = u235(energies, 0.0)
    assert total[1] == pytest.approx(90.64895383)
    total, absorption, fission = u235(energies, 300.0)
    assert total[1] == pytest.approx(91.12534964)


def test_evaluate_mlbw(u234):
    """Make sure a Multi-Level Breit-Wigner multipole object can be called."""
    energies = [1e-3, 1.0, 10.0, 50.]
    total, absorption, fission = u234(energies, 0.0)
    assert total[3] == pytest.approx(15.02827953)
    total, absorption, fission = u234(energies, 300.0)
    assert total[3] == pytest.approx(15.08269143)


def test_high_l(fe56):
    """Test a nuclide (Fe56) with a high l-value (4)."""
    energies = [1e-3, 1.0, 10.0, 1e3, 1e5]
    total, absorption, fission = fe56(energies, 0.0)
    assert total[0] == pytest.approx(25.072619556789267)
    total, absorption, fission = fe56(energies, 300.0)
    assert total[0] == pytest.approx(27.85535792368082)


def test_export_to_hdf5(tmpdir, u235):
    filename = str(tmpdir.join('092235.h5'))
    u235.export_to_hdf5(filename)
    assert os.path.exists(filename)
