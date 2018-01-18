#!/usr/bin/env python

import os

import numpy as np
import pytest
import openmc.data


@pytest.fixture(scope='module')
def u235():
    directory = os.environ['OPENMC_MULTIPOLE_LIBRARY']
    filename = os.path.join(directory, '092235.h5')
    return openmc.data.WindowedMultipole.from_hdf5(filename)


@pytest.fixture(scope='module')
def fe56():
    directory = os.environ['OPENMC_MULTIPOLE_LIBRARY']
    filename = os.path.join(directory, '026056.h5')
    return openmc.data.WindowedMultipole.from_hdf5(filename)


def test_evaluate(u235):
    """Make sure multipole object can be called."""
    energies = [1e-3, 1.0, 10.0, 50.]
    total, absorption, fission = u235(energies, 0.0)
    assert total[1] == pytest.approx(90.64895383)
    total, absorption, fission = u235(energies, 300.0)
    assert total[1] == pytest.approx(91.12534964)


def test_high_l(fe56):
    """Test a nuclide (Fe56) with a high l-value (4)."""
    energies = [1e-3, 1.0, 10.0, 1e3, 1e5]
    total, absorption, fission = fe56(energies, 0.0)
    assert total[0] == pytest.approx(25.072619556789267)
    total, absorption, fission = fe56(energies, 300.0)
    assert total[0] == pytest.approx(27.85535792368082)
