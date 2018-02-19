""" Tests for the AtomNumber class """

import numpy as np
from openmc.deplete import atom_number


def test_indexing():
    """Tests the __getitem__ and __setitem__ routines simultaneously."""

    local_mats = ["10000", "10001"]
    nuclides = ["U238", "U235", "U234"]
    volume = {"10000" : 0.38, "10001" : 0.21}

    number = atom_number.AtomNumber(local_mats, nuclides, volume, 2)

    number["10000", "U238"] = 1.0
    number["10001", "U238"] = 2.0
    number["10000", "U235"] = 3.0
    number["10001", "U235"] = 4.0

    # String indexing
    assert number["10000", "U238"] == 1.0
    assert number["10001", "U238"] == 2.0
    assert number["10000", "U235"] == 3.0
    assert number["10001", "U235"] == 4.0

    # Int indexing
    assert number[0, 0] == 1.0
    assert number[1, 0] == 2.0
    assert number[0, 1] == 3.0
    assert number[1, 1] == 4.0

    number[0, 0] = 5.0

    assert number[0, 0] == 5.0
    assert number["10000", "U238"] == 5.0


def test_properties():
    """Test properties. """
    local_mats = ["10000", "10001"]
    nuclides = ["U238", "U235", "Gd157"]
    volume = {"10000" : 0.38, "10001" : 0.21}

    number = atom_number.AtomNumber(local_mats, nuclides, volume, 2)

    assert list(number.materials) == ["10000", "10001"]
    assert number.n_nuc == 3
    assert list(number.nuclides) == ["U238", "U235", "Gd157"]
    assert number.burnable_nuclides == ["U238", "U235"]


def test_density_indexing():
    """Tests the get and set_atom_density routines simultaneously."""

    local_mats = ["10000", "10001", "10002"]
    nuclides = ["U238", "U235", "U234"]
    volume = {"10000" : 0.38, "10001" : 0.21}

    number = atom_number.AtomNumber(local_mats, nuclides, volume, 2)

    number.set_atom_density("10000", "U238", 1.0)
    number.set_atom_density("10001", "U238", 2.0)
    number.set_atom_density("10002", "U238", 3.0)
    number.set_atom_density("10000", "U235", 4.0)
    number.set_atom_density("10001", "U235", 5.0)
    number.set_atom_density("10002", "U235", 6.0)
    number.set_atom_density("10000", "U234", 7.0)
    number.set_atom_density("10001", "U234", 8.0)
    number.set_atom_density("10002", "U234", 9.0)

    # String indexing
    assert number.get_atom_density("10000", "U238") == 1.0
    assert number.get_atom_density("10001", "U238") == 2.0
    assert number.get_atom_density("10002", "U238") == 3.0
    assert number.get_atom_density("10000", "U235") == 4.0
    assert number.get_atom_density("10001", "U235") == 5.0
    assert number.get_atom_density("10002", "U235") == 6.0
    assert number.get_atom_density("10000", "U234") == 7.0
    assert number.get_atom_density("10001", "U234") == 8.0
    assert number.get_atom_density("10002", "U234") == 9.0

    # Int indexing
    assert number.get_atom_density(0, 0) == 1.0
    assert number.get_atom_density(1, 0) == 2.0
    assert number.get_atom_density(2, 0) == 3.0
    assert number.get_atom_density(0, 1) == 4.0
    assert number.get_atom_density(1, 1) == 5.0
    assert number.get_atom_density(2, 1) == 6.0
    assert number.get_atom_density(0, 2) == 7.0
    assert number.get_atom_density(1, 2) == 8.0
    assert number.get_atom_density(2, 2) == 9.0


    number.set_atom_density(0, 0, 5.0)
    assert number.get_atom_density(0, 0) == 5.0

    # Verify volume is used correctly
    assert number[0, 0] == 5.0 * 0.38
    assert number[1, 0] == 2.0 * 0.21
    assert number[2, 0] == 3.0 * 1.0
    assert number[0, 1] == 4.0 * 0.38
    assert number[1, 1] == 5.0 * 0.21
    assert number[2, 1] == 6.0 * 1.0
    assert number[0, 2] == 7.0 * 0.38
    assert number[1, 2] == 8.0 * 0.21
    assert number[2, 2] == 9.0 * 1.0


def test_get_mat_slice():
    """Tests getting slices."""

    local_mats = ["10000", "10001", "10002"]
    nuclides = ["U238", "U235", "U234"]
    volume = {"10000" : 0.38, "10001" : 0.21}

    number = atom_number.AtomNumber(local_mats, nuclides, volume, 2)

    number.number = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]])

    sl = number.get_mat_slice(0)

    np.testing.assert_array_equal(sl, np.array([1.0, 2.0]))

    sl = number.get_mat_slice("10000")

    np.testing.assert_array_equal(sl, np.array([1.0, 2.0]))


def test_set_mat_slice():
    """Tests getting slices."""

    local_mats = ["10000", "10001", "10002"]
    nuclides = ["U238", "U235", "U234"]
    volume = {"10000" : 0.38, "10001" : 0.21}

    number = atom_number.AtomNumber(local_mats, nuclides, volume, 2)

    number.set_mat_slice(0, [1.0, 2.0])

    assert number[0, 0] == 1.0
    assert number[0, 1] == 2.0

    number.set_mat_slice("10000", [3.0, 4.0])

    assert number[0, 0] == 3.0
    assert number[0, 1] == 4.0
