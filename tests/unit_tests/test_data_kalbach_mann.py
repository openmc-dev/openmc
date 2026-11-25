"""Test of the Kalbach-Mann slope calculation when data are
retrieved from ENDF files."""

import os
from pathlib import Path
import pytest

import numpy as np

import openmc
from openmc.data import IncidentNeutron
from openmc.data.data import EV_PER_MEV, NEUTRON_MASS_EV
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

    # Only neutron (ZA=1) and photon (ZA=0) projectiles are covered by these
    # systematics. Anything else must say so rather than returning a number or
    # failing as a KeyError from inside za_to_M.
    for za_projectile in (1000, 1001, 1003, 2003, 2004):
        with pytest.raises(NotImplementedError):
            kalbach_slope(
                energy_projectile=energy_projectile,
                energy_emitted=energy_emitted,
                za_projectile=za_projectile,
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


@pytest.mark.parametrize('e_gamma,e_b_cm', [
    (15.0e6, 3.0e6), (40.0e6, 10.0e6), (100.0e6, 90.0e6)])
def test_kalbach_slope_photon(e_gamma, e_b_cm):
    """Eq. 6.5 of the ENDF-6 Formats Manual (BNL-224854-2023, section 6.2):

        a_gamma = a_n(E_gamma, E_b_cm) * sqrt(E_gamma/(2 m_n))
                  * min(4, max(1, 9.3/sqrt(E_b_cm)))

    with E_b_cm and m_n in MeV, and a_n evaluated by plugging E_gamma into the
    incident-neutron slot. Both of those are easy to get subtly wrong -- the
    emission channel energy epsilon_b and the true photon compound system are
    the tempting substitutions, and both are incorrect here -- so the formula
    is pinned against an independent evaluation.
    """
    slope_n = kalbach_slope(e_gamma, e_b_cm, 1, 1, 82208)
    expected = (slope_n*np.sqrt(e_gamma/(2.0*NEUTRON_MASS_EV))
                * min(4.0, max(1.0, 9.3/np.sqrt(e_b_cm/EV_PER_MEV))))

    got = kalbach_slope(e_gamma, e_b_cm, 0, 1, 82208)
    assert got == pytest.approx(expected, rel=1e-12)
    # A photon carries less momentum than a nucleon of the same energy
    assert 0.0 < got < slope_n


def test_kalbach_slope_photon_clip_saturates():
    """The clipping factor saturates at 4 below 5.41 MeV and at 1 above
    86.5 MeV, so the ratio to the unclipped scaling is flat outside that
    window."""
    def ratio(e_b_cm):
        return (kalbach_slope(20.0e6, e_b_cm, 0, 1, 82208)
                / kalbach_slope(20.0e6, e_b_cm, 1, 1, 82208))

    assert ratio(1.0e6) == pytest.approx(ratio(3.0e6), rel=1e-12)
    assert ratio(9.0e7) == pytest.approx(ratio(1.0e8), rel=1e-12)


def test_kalbach_slope_photon_zero_outgoing_energy():
    """A zero outgoing energy is a normal first grid point of an ENDF
    LAW=1/LANG=2 table, and must not raise a divide-by-zero warning."""
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter('error')
        assert kalbach_slope(15.0e6, 0.0, 0, 1, 82208) >= 0.0


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
