"""Test of the Kalbach-Mann slope calculation when data are
retrieved from ENDF files."""

import os
from pathlib import Path
import pytest

import numpy as np

from openmc.data import IncidentNeutron
from openmc.data.kalbach_mann import _separation_energy, _AtomicRepresentation
from openmc.data import kalbach_slope
from openmc.data import KalbachMann

from . import needs_njoy


@pytest.fixture(scope='module')
def neutron():
    """Neutron AtomicRepresentation."""
    return _AtomicRepresentation(z=0, a=1)


@pytest.fixture(scope='module')
def triton():
    """Triton AtomicRepresentation."""
    return _AtomicRepresentation(z=1, a=3)


@pytest.fixture(scope='module')
def b10():
    """B10 AtomicRepresentation."""
    return _AtomicRepresentation(z=5, a=10)


@pytest.fixture(scope='module')
def c12():
    """C12 AtomicRepresentation."""
    return _AtomicRepresentation(z=6, a=12)


@pytest.fixture(scope='module')
def c13():
    """C13 AtomicRepresentation."""
    return _AtomicRepresentation(z=6, a=13)


@pytest.fixture(scope='module')
def na23():
    """Na23 AtomicRepresentation."""
    return _AtomicRepresentation(z=11, a=23)


def test_atomic_representation(neutron, triton, b10, c12, c13, na23):
    """Test the _AtomicRepresentation class."""
    # Test instantiation from_za
    assert b10 == _AtomicRepresentation.from_za(5010)

    # Test addition
    assert c13 + b10 == na23

    # Test substraction
    assert c13 - c12 == neutron
    assert c13 - b10 == triton

    # Test properties when no information for Kalbach-Mann are given
    assert c13.a == 13
    assert c13.z == 6
    assert c13.n == 7
    assert c13.za == 6013

    # Test properties when information for Kalbach-Mann are given
    assert triton.a == 3
    assert triton.z == 1
    assert triton.n == 2
    assert triton.za == 1003

    # Test instantiation errors
    with pytest.raises(ValueError):
        _AtomicRepresentation(z=5, a=1)
    with pytest.raises(ValueError):
        _AtomicRepresentation(z=-1, a=1)
    with pytest.raises(ValueError):
        _AtomicRepresentation(z=5, a=0)
    with pytest.raises(ValueError):
        _AtomicRepresentation(z=5, a=-2)
    with pytest.raises(ValueError):
        neutron - triton


def test_separation_energy(triton, b10, c13):
    """Comparison to hand-calculations on a simple example."""
    assert _separation_energy(
        compound=c13,
        nucleus=b10,
        particle=triton
    ) == pytest.approx(18.6880713)


def test_kalbach_slope():
    """Comparison to hand-calculations for n + c12 -> c13 -> triton + b10."""
    energy_projectile = 10.2  # [eV]
    energy_emitted = 5.4  # [eV]

    # Check that NotImplementedError is raised if the projectile is not
    # a neutron
    with pytest.raises(NotImplementedError):
        kalbach_slope(
            energy_projectile=energy_projectile,
            energy_emitted=energy_emitted,
            za_projectile=1000,
            za_emitted=1,
            za_target=6012
        )

    assert kalbach_slope(
        energy_projectile=energy_projectile,
        energy_emitted=energy_emitted,
        za_projectile=1,
        za_emitted=1003,
        za_target=6012
    ) == pytest.approx(0.8409921475)


@pytest.mark.parametrize(
    "hdf5_filename, endf_filename", [
        ('O16.h5', 'n-008_O_016.endf'),
        ('Ca46.h5', 'n-020_Ca_046.endf'),
        ('Hg204.h5', 'n-080_Hg_204.endf')
    ]
)
def test_comparison_slope_hdf5(hdf5_filename, endf_filename, endf_data):
    """Test the calculation of the Kalbach-Mann slope done by OpenMC
    by comparing it to HDF5 data. The test is based on the first product
    of MT=5 (neutron). The isotopes tested have been selected because the
    corresponding products in ENDF/B-VII.1 are described using MF=6, LAW=1,
    LANG=2 (i.e., Kalbach-Mann systematics) and the slope is not given
    explicitly.

    If an error occurs during the "validity check", this means that
    the nuclear data evaluation has evolved and the distribution might
    no longer be described using Kalbach-Mann systematics. Another
    isotope needs to be identified and tested.

    Warning: This test is valid as long as ENDF files are not directly
    used to generate the HDF5 files used in the tests.

    """
    # HDF5 data
    hdf5_directory = Path(openmc.config.get('cross_sections')).parent
    hdf5_data = IncidentNeutron.from_hdf5(hdf5_directory / hdf5_filename)
    hdf5_product = hdf5_data[5].products[0]
    hdf5_distribution = hdf5_product.distribution[0]

    # ENDF data
    endf_directory = Path(endf_data)
    endf_path = endf_directory / 'neutrons' / endf_filename
    endf_data = IncidentNeutron.from_endf(endf_path)
    endf_product = endf_data[5].products[0]
    endf_distribution = endf_product.distribution[0]

    # Validity check
    assert isinstance(endf_distribution, KalbachMann)
    assert isinstance(hdf5_distribution, KalbachMann)
    assert endf_product.particle == hdf5_product.particle
    assert len(endf_distribution.slope) == len(hdf5_distribution.slope)

    # Results check
    for i, hdf5_slope in enumerate(hdf5_distribution.slope):
        assert endf_distribution._calculated_slope[i]

        np.testing.assert_array_almost_equal(
            endf_distribution.slope[i].y,
            hdf5_slope.y,
            decimal=5
        )
