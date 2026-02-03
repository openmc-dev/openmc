"""Unit tests for ParticleType class."""

import pytest
from openmc import ParticleType


# Tests for creating ParticleType instances

def test_create_from_int():
    """Test creation from PDG number."""
    p = ParticleType(2112)
    assert p.pdg_number == 2112
    assert int(p) == 2112


def test_create_from_string_name():
    """Test creation from particle name."""
    p = ParticleType('neutron')
    assert p.pdg_number == 2112

    p = ParticleType('photon')
    assert p.pdg_number == 22


def test_create_from_string_aliases():
    """Test creation from particle aliases."""
    assert ParticleType('n').pdg_number == 2112
    assert ParticleType('gamma').pdg_number == 22
    assert ParticleType('p').pdg_number == 2212
    assert ParticleType('proton').pdg_number == 2212
    assert ParticleType('d').pdg_number == 1000010020
    assert ParticleType('t').pdg_number == 1000010030


def test_create_from_string_nuclide():
    """Test creation from GNDS nuclide name."""
    p = ParticleType('He4')
    assert p.pdg_number == 1000020040

    p = ParticleType('U235')
    assert p.pdg_number == 1000922350

    p = ParticleType('Am242_m1')
    assert p.pdg_number == 1000952421


def test_create_from_string_pdg_prefix():
    """Test creation with pdg: prefix."""
    p = ParticleType('pdg:2112')
    assert p.pdg_number == 2112

    p = ParticleType('PDG:22')
    assert p.pdg_number == 22


def test_create_from_particle_type():
    """Test creation from existing ParticleType."""
    p1 = ParticleType(2112)
    p2 = ParticleType(p1)
    assert p1 == p2
    assert p1 is not p2  # Different instances


def test_legacy_particle_indices():
    """Test backward compatibility with legacy indices 0-3."""
    assert ParticleType(0) == ParticleType.NEUTRON
    assert ParticleType(1) == ParticleType.PHOTON
    assert ParticleType(2) == ParticleType.ELECTRON
    assert ParticleType(3) == ParticleType.POSITRON


def test_create_invalid_type():
    """Test creation with invalid type raises TypeError."""
    with pytest.raises(TypeError):
        ParticleType([2112])

    with pytest.raises(TypeError):
        ParticleType({'pdg': 2112})


def test_create_invalid_string():
    """Test creation with invalid string raises ValueError."""
    with pytest.raises(ValueError):
        ParticleType('')

    with pytest.raises(ValueError):
        ParticleType('pdg:invalid')


def test_create_case_insensitive():
    """Test that string parsing is case insensitive for aliases."""
    assert ParticleType('NEUTRON').pdg_number == 2112
    assert ParticleType('Neutron').pdg_number == 2112
    assert ParticleType('PHOTON').pdg_number == 22


# Tests for equality and comparison

def test_equality_same_pdg():
    """Test that instances with same PDG are equal."""
    p1 = ParticleType(2112)
    p2 = ParticleType(2112)
    assert p1 == p2


def test_equality_int():
    """Test equality comparison with int."""
    p = ParticleType(2112)
    assert p == 2112
    assert 2112 == p


def test_equality_string():
    """Test equality comparison with string."""
    p = ParticleType(2112)
    assert p == 'neutron'
    assert p == 'pdg:2112'
    assert p == 'n'


def test_inequality():
    """Test inequality comparisons."""
    p1 = ParticleType(2112)
    p2 = ParticleType(22)
    assert p1 != p2
    assert p1 != 22
    assert p1 != 'photon'


def test_equality_with_constants():
    """Test equality with class constants."""
    p = ParticleType(2112)
    assert p == ParticleType.NEUTRON
    assert ParticleType.NEUTRON == p


def test_equality_invalid_string():
    """Test equality with invalid string returns False."""
    p = ParticleType(2112)
    assert not (p == 'invalid_particle')
    assert p != 'invalid_particle'


# Tests for hashing behavior

def test_hash_consistency():
    """Test that equal instances hash to the same value."""
    p1 = ParticleType(2112)
    p2 = ParticleType(2112)
    assert hash(p1) == hash(p2)


def test_set_deduplication():
    """Test that equal instances deduplicate in sets."""
    p1 = ParticleType(2112)
    p2 = ParticleType(2112)
    s = {p1, p2}
    assert len(s) == 1


def test_dict_key():
    """Test use as dictionary key."""
    p1 = ParticleType(2112)
    p2 = ParticleType(2112)
    d = {p1: 'neutron'}
    assert d[p2] == 'neutron'


def test_hash_different_particles():
    """Test that different particles have different hashes (usually)."""
    p1 = ParticleType(2112)
    p2 = ParticleType(22)
    # Different PDG numbers should (almost always) have different hashes
    assert hash(p1) != hash(p2)


# Tests for properties and computed attributes

def test_pdg_number_property():
    """Test pdg_number property."""
    p = ParticleType(2112)
    assert p.pdg_number == 2112
    assert isinstance(p.pdg_number, int)


def test_int_conversion():
    """Test __int__ conversion."""
    p = ParticleType(1000020040)
    assert int(p) == 1000020040


def test_zam_elementary():
    """Test zam property for elementary particles."""
    assert ParticleType.NEUTRON.zam is None
    assert ParticleType.PHOTON.zam is None
    assert ParticleType.ELECTRON.zam is None
    assert ParticleType.POSITRON.zam is None


def test_zam_nucleus():
    """Test zam property for nuclear particles."""
    he4 = ParticleType('He4')
    assert he4.zam == (2, 4, 0)

    u235 = ParticleType('U235')
    Z, A, m = u235.zam
    assert Z == 92
    assert A == 235
    assert m == 0


def test_zam_metastable():
    """Test zam property for metastable nuclei."""
    # Am242m has m=1
    am242m = ParticleType('Am242_m1')
    Z, A, m = am242m.zam
    assert Z == 95
    assert A == 242
    assert m == 1


def test_is_nucleus_false():
    """Test is_nucleus for elementary particles."""
    assert not ParticleType.NEUTRON.is_nucleus
    assert not ParticleType.PHOTON.is_nucleus
    assert not ParticleType.ELECTRON.is_nucleus
    assert not ParticleType.POSITRON.is_nucleus


def test_is_nucleus_true():
    """Test is_nucleus for nuclear particles."""
    assert ParticleType.ALPHA.is_nucleus
    assert ParticleType('He4').is_nucleus
    assert ParticleType('U235').is_nucleus
    assert ParticleType.DEUTERON.is_nucleus
    assert ParticleType.TRITON.is_nucleus


# Tests for __str__ and __repr__

def test_str_elementary():
    """Test string representation of elementary particles."""
    assert str(ParticleType.NEUTRON) == 'neutron'
    assert str(ParticleType.PHOTON) == 'photon'
    assert str(ParticleType.ELECTRON) == 'electron'
    assert str(ParticleType.POSITRON) == 'positron'
    assert str(ParticleType.PROTON) == 'H1'


def test_str_nucleus():
    """Test string representation of nuclei."""
    assert str(ParticleType.ALPHA) == 'He4'
    assert str(ParticleType('U235')) == 'U235'
    assert str(ParticleType.DEUTERON) == 'H2'
    assert str(ParticleType.TRITON) == 'H3'


def test_str_arbitrary_pdg():
    """Test string representation of arbitrary PDG number."""
    # PDG number that doesn't match any known particle
    p = ParticleType(12345)
    assert str(p) == 'pdg:12345'


def test_repr():
    """Test repr includes PDG number."""
    p = ParticleType.NEUTRON
    repr_str = repr(p)
    assert 'ParticleType' in repr_str
    assert 'PDG=2112' in repr_str
    assert 'neutron' in repr_str


def test_repr_nucleus():
    """Test repr for nuclear particles."""
    p = ParticleType.ALPHA
    repr_str = repr(p)
    assert 'ParticleType' in repr_str
    assert 'PDG=1000020040' in repr_str
    assert 'He4' in repr_str


def test_create_from_numpy_int():
    """Test creation from numpy integer types."""
    import numpy as np
    p = ParticleType(np.int32(2112))
    assert p.pdg_number == 2112

    p = ParticleType(np.int64(22))
    assert p.pdg_number == 22


def test_equality_numpy_int():
    """Test equality comparison with numpy integer types."""
    import numpy as np
    p = ParticleType(2112)
    assert p == np.int32(2112)
    assert p == np.int64(2112)
