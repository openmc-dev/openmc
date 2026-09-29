"""Tests for the level inelastic (ENDF LAW=3) energy distribution.

The HDF5 attribute set of this distribution changed to carry a projectile, so
the tests below pin both directions of that contract: a file written now must
still be readable by an OpenMC built before the change, and a file written
before it must still be readable now.
"""

import h5py
import numpy as np
import pytest

import openmc.data
from openmc.data.data import NEUTRON_MASS_EV
from openmc.data.energy_distribution import LevelInelastic


@pytest.mark.parametrize('particle', ['neutron', 'photon'])
def test_roundtrip(run_in_tmpdir, particle):
    dist = LevelInelastic(-7.368e6, 206.19, particle)
    with h5py.File('level.h5', 'w') as f:
        dist.to_hdf5(f.create_group('d'))
    with h5py.File('level.h5', 'r') as f:
        back = LevelInelastic.from_hdf5(f['d'])

    assert back.q_value == pytest.approx(dist.q_value)
    assert back.mass == pytest.approx(dist.mass)
    assert back.particle == particle
    assert back.threshold == pytest.approx(dist.threshold)


def test_neutron_file_still_readable_by_older_openmc(run_in_tmpdir):
    """A neutron distribution must keep writing the pre-3.1 attributes.

    An OpenMC built before the q_value / mass / particle form existed reads
    `threshold` and `mass_ratio` and nothing else. It does not check the HDF5
    status of that read, so if they are absent it does not fail -- it samples
    from an uninitialised member. Dropping them would silently corrupt every
    neutron calculation run against a regenerated library.
    """
    A, q_value = 206.19, -7.368e6
    dist = LevelInelastic(q_value, A, 'neutron')
    with h5py.File('level.h5', 'w') as f:
        dist.to_hdf5(f.create_group('d'))

    with h5py.File('level.h5', 'r') as f:
        attrs = f['d'].attrs
        assert 'threshold' in attrs
        assert 'mass_ratio' in attrs
        # The two forms describe the same law: E' = mass_ratio*(E - threshold)
        assert attrs['threshold'] == pytest.approx((A + 1.0)/A*abs(q_value))
        assert attrs['mass_ratio'] == pytest.approx((A/(A + 1.0))**2)


def test_photon_file_has_no_legacy_attributes(run_in_tmpdir):
    """There is no pre-3.1 form for a photon projectile, and an older OpenMC
    has no photonuclear support, so it must never see one."""
    with h5py.File('level.h5', 'w') as f:
        LevelInelastic(-7.368e6, 206.19, 'photon').to_hdf5(f.create_group('d'))
    with h5py.File('level.h5', 'r') as f:
        assert 'threshold' not in f['d'].attrs
        assert 'mass_ratio' not in f['d'].attrs


def test_legacy_file_still_readable(run_in_tmpdir):
    """A file written before this change carries only the neutron-only form."""
    A = 206.19
    with h5py.File('old.h5', 'w') as f:
        g = f.create_group('d')
        g.attrs['type'] = np.bytes_('level')
        g.attrs['threshold'] = (A + 1.0)/A*7.368e6
        g.attrs['mass_ratio'] = (A/(A + 1.0))**2
    with h5py.File('old.h5', 'r') as f:
        back = LevelInelastic.from_hdf5(f['d'])

    assert back.mass == pytest.approx(A)
    assert back.q_value == pytest.approx(-7.368e6)
    assert back.particle == 'neutron'


def test_photon_threshold():
    """The reaction opens where the outgoing energy first becomes positive,
    i.e. at the smaller root of E^2/(2 A m_n) - E + |Q| = 0."""
    A, q_value = 206.19, -7.368e6
    dist = LevelInelastic(q_value, A, 'photon')

    b = NEUTRON_MASS_EV*A
    expected = np.sqrt(b)*(np.sqrt(b) - np.sqrt(b - 2.0*abs(q_value)))
    assert dist.threshold == pytest.approx(expected)

    # A photon threshold sits just above |Q|; the neutron one is larger by
    # the recoil factor (A+1)/A
    assert abs(q_value) < dist.threshold < 1.001*abs(q_value)
    assert LevelInelastic(q_value, A, 'neutron').threshold > dist.threshold


def test_hdf5_version_is_bumped():
    """The level attribute set changed, so the format version must say so."""
    assert openmc.data.HDF5_VERSION >= (3, 1)


class _FakeACE:
    """Minimal stand-in for an ACE table holding a single LAW=3 entry.

    `LevelInelastic.from_ace` reads the table type, the two XSS words of the
    distribution and the atomic weight ratio of the header, and nothing else.
    """

    def __init__(self, data_type, threshold, mass_ratio, awr):
        self.name = f'82000.30{data_type}'
        self.data_type = openmc.data.ace.TableType(data_type)
        self.atomic_weight_ratio = awr
        self.xss = np.array([threshold/openmc.data.EV_PER_MEV, mass_ratio])


def _neutron_words(A, q_value):
    """The LAW=3 words as written with neutron kinematics."""
    return (A + 1.0)/A*abs(q_value), (A/(A + 1.0))**2


def _photon_words(A, q_value):
    """The LAW=3 words as written with photonuclear kinematics."""
    return abs(q_value), (A - 1.0)/A


def test_from_ace_neutron():
    A, q_value = 206.19, -7.368e6
    threshold, mass_ratio = _neutron_words(A, q_value)
    dist = LevelInelastic.from_ace(_FakeACE('c', threshold, mass_ratio, A), 0)

    assert dist.particle == 'neutron'
    assert dist.mass == pytest.approx(A)
    assert dist.q_value == pytest.approx(q_value)


def test_from_ace_photonuclear():
    A, q_value = 206.19, -7.368e6
    threshold, mass_ratio = _photon_words(A, q_value)
    dist = LevelInelastic.from_ace(_FakeACE('u', threshold, mass_ratio, A), 0)

    assert dist.particle == 'photon'
    assert dist.mass == pytest.approx(A)
    assert dist.q_value == pytest.approx(q_value)


def test_from_ace_photonuclear_with_neutron_kinematics():
    """Older NJOY wrote photonuclear LAW=3 parameters with neutron kinematics.

    Those tables are what the photonuclear libraries in circulation are made
    of, so reading one must warn and recover |Q| with the convention it was
    written in -- not refuse the table and take library generation down with
    it.
    """
    A, q_value = 206.19, -7.368e6
    threshold, mass_ratio = _neutron_words(A, q_value)
    ace = _FakeACE('u', threshold, mass_ratio, A)

    with pytest.warns(UserWarning, match='neutron kinematics'):
        dist = LevelInelastic.from_ace(ace, 0)

    # The projectile comes from the table type, so the distribution is still
    # photonuclear; only the reading of the two words changes.
    assert dist.particle == 'photon'
    assert dist.mass == pytest.approx(A)
    assert dist.q_value == pytest.approx(q_value)

    # Reading it as though it carried the modern words would put the mass at
    # roughly A/2 and |Q| out by the recoil factor (A + 1)/A.
    assert dist.q_value != pytest.approx(-threshold)


def test_from_ace_inconsistent_mass_raises():
    """A mass ratio matching neither convention is a broken or misparsed
    table, and the caller hears about it."""
    ace = _FakeACE('u', 7.368e6, 0.5, 206.19)
    with pytest.raises(ValueError, match='atomic weight ratio'):
        LevelInelastic.from_ace(ace, 0)

    ace = _FakeACE('c', 7.368e6, 0.5, 206.19)
    with pytest.raises(ValueError, match='atomic weight ratio'):
        LevelInelastic.from_ace(ace, 0)


def test_from_ace_nonphysical_mass_ratio_raises():
    """A mass ratio outside the range either convention can produce cannot be
    inverted at all, and must not raise something opaque like a domain or
    zero-division error."""
    ace = _FakeACE('u', 7.368e6, -1.0, 206.19)
    with pytest.raises(ValueError, match='nothing physical'):
        LevelInelastic.from_ace(ace, 0)


def test_from_ace_unsupported_table_type():
    ace = _FakeACE('h', *_neutron_words(206.19, -7.368e6), 206.19)
    with pytest.raises(NotImplementedError):
        LevelInelastic.from_ace(ace, 0)
