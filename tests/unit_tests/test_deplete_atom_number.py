""" Tests for atom_number.py. """

import unittest

import numpy as np
from openmc.deplete import atom_number


class TestAtomNumber(unittest.TestCase):
    """Tests for the AtomNumber class."""

    def test_indexing(self):
        """Tests the __getitem__ and __setitem__ routines simultaneously."""

        mat_to_ind = {"10000" : 0, "10001" : 1, "10002" : 2}
        nuc_to_ind = {"U238" : 0, "U235" : 1, "U234" : 2}
        volume = {"10000" : 0.38, "10001" : 0.21}

        number = atom_number.AtomNumber(mat_to_ind, nuc_to_ind, volume, 2, 2)

        number["10000", "U238"] = 1.0
        number["10001", "U238"] = 2.0
        number["10000", "U235"] = 3.0
        number["10001", "U235"] = 4.0

        # String indexing
        self.assertEqual(number["10000", "U238"], 1.0)
        self.assertEqual(number["10001", "U238"], 2.0)
        self.assertEqual(number["10000", "U235"], 3.0)
        self.assertEqual(number["10001", "U235"], 4.0)

        # Int indexing
        self.assertEqual(number[0, 0], 1.0)
        self.assertEqual(number[1, 0], 2.0)
        self.assertEqual(number[0, 1], 3.0)
        self.assertEqual(number[1, 1], 4.0)

        number[0, 0] = 5.0

        self.assertEqual(number[0, 0], 5.0)
        self.assertEqual(number["10000", "U238"], 5.0)

    def test_n_mat(self):
        """Test number of materials property. """
        mat_to_ind = {"10000" : 0, "10001" : 1}
        nuc_to_ind = {"U238" : 0, "U235" : 1, "Gd157" : 2}
        volume = {"10000" : 0.38, "10001" : 0.21}

        number = atom_number.AtomNumber(mat_to_ind, nuc_to_ind, volume, 2, 2)

        self.assertEqual(number.n_mat, 2)

    def test_n_nuc(self):
        """Test number of nuclides property."""
        mat_to_ind = {"10000" : 0, "10001" : 1}
        nuc_to_ind = {"U238" : 0, "U235" : 1, "Gd157" : 2}
        volume = {"10000" : 0.38, "10001" : 0.21}

        number = atom_number.AtomNumber(mat_to_ind, nuc_to_ind, volume, 2, 2)

        self.assertEqual(number.n_nuc, 3)

    def test_burn_nuc_list(self):
        """Test the list of burned nuclides property"""
        mat_to_ind = {"10000" : 0, "10001" : 1}
        nuc_to_ind = {"U238" : 0, "U235" : 1, "Gd157" : 2}
        volume = {"10000" : 0.38, "10001" : 0.21}

        number = atom_number.AtomNumber(mat_to_ind, nuc_to_ind, volume, 2, 2)

        self.assertEqual(number.burn_nuc_list, ["U238", "U235"])

    def test_burn_mat_list(self):
        """Test the list of burned nuclides property"""
        mat_to_ind = {"10000" : 0, "10001" : 1, "10002" : 2}
        nuc_to_ind = {"U238" : 0, "U235" : 1, "Gd157" : 2}
        volume = {"10000" : 0.38, "10001" : 0.21}

        number = atom_number.AtomNumber(mat_to_ind, nuc_to_ind, volume, 2, 2)

        self.assertEqual(number.burn_mat_list, ["10000", "10001"])

    def test_density_indexing(self):
        """Tests the get and set_atom_density routines simultaneously."""

        mat_to_ind = {"10000" : 0, "10001" : 1, "10002" : 2}
        nuc_to_ind = {"U238" : 0, "U235" : 1, "U234" : 2}
        volume = {"10000" : 0.38, "10001" : 0.21}

        number = atom_number.AtomNumber(mat_to_ind, nuc_to_ind, volume, 2, 2)

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
        self.assertEqual(number.get_atom_density("10000", "U238"), 1.0)
        self.assertEqual(number.get_atom_density("10001", "U238"), 2.0)
        self.assertEqual(number.get_atom_density("10002", "U238"), 3.0)
        self.assertEqual(number.get_atom_density("10000", "U235"), 4.0)
        self.assertEqual(number.get_atom_density("10001", "U235"), 5.0)
        self.assertEqual(number.get_atom_density("10002", "U235"), 6.0)
        self.assertEqual(number.get_atom_density("10000", "U234"), 7.0)
        self.assertEqual(number.get_atom_density("10001", "U234"), 8.0)
        self.assertEqual(number.get_atom_density("10002", "U234"), 9.0)

        # Int indexing
        self.assertEqual(number.get_atom_density(0, 0), 1.0)
        self.assertEqual(number.get_atom_density(1, 0), 2.0)
        self.assertEqual(number.get_atom_density(2, 0), 3.0)
        self.assertEqual(number.get_atom_density(0, 1), 4.0)
        self.assertEqual(number.get_atom_density(1, 1), 5.0)
        self.assertEqual(number.get_atom_density(2, 1), 6.0)
        self.assertEqual(number.get_atom_density(0, 2), 7.0)
        self.assertEqual(number.get_atom_density(1, 2), 8.0)
        self.assertEqual(number.get_atom_density(2, 2), 9.0)


        number.set_atom_density(0, 0, 5.0)

        self.assertEqual(number.get_atom_density(0, 0), 5.0)

        # Verify volume is used correctly
        self.assertEqual(number[0, 0], 5.0 * 0.38)
        self.assertEqual(number[1, 0], 2.0 * 0.21)
        self.assertEqual(number[2, 0], 3.0 * 1.0)
        self.assertEqual(number[0, 1], 4.0 * 0.38)
        self.assertEqual(number[1, 1], 5.0 * 0.21)
        self.assertEqual(number[2, 1], 6.0 * 1.0)
        self.assertEqual(number[0, 2], 7.0 * 0.38)
        self.assertEqual(number[1, 2], 8.0 * 0.21)
        self.assertEqual(number[2, 2], 9.0 * 1.0)

    def test_get_mat_slice(self):
        """Tests getting slices."""

        mat_to_ind = {"10000" : 0, "10001" : 1, "10002" : 2}
        nuc_to_ind = {"U238" : 0, "U235" : 1, "U234" : 2}
        volume = {"10000" : 0.38, "10001" : 0.21}

        number = atom_number.AtomNumber(mat_to_ind, nuc_to_ind, volume, 2, 2)

        number.number = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]])

        sl = number.get_mat_slice(0)

        np.testing.assert_array_equal(sl, np.array([1.0, 2.0]))

        sl = number.get_mat_slice("10000")

        np.testing.assert_array_equal(sl, np.array([1.0, 2.0]))

    def test_set_mat_slice(self):
        """Tests getting slices."""

        mat_to_ind = {"10000" : 0, "10001" : 1, "10002" : 2}
        nuc_to_ind = {"U238" : 0, "U235" : 1, "U234" : 2}
        volume = {"10000" : 0.38, "10001" : 0.21}

        number = atom_number.AtomNumber(mat_to_ind, nuc_to_ind, volume, 2, 2)

        number.set_mat_slice(0, [1.0, 2.0])

        self.assertEqual(number[0, 0], 1.0)
        self.assertEqual(number[0, 1], 2.0)

        number.set_mat_slice("10000", [3.0, 4.0])

        self.assertEqual(number[0, 0], 3.0)
        self.assertEqual(number[0, 1], 4.0)


if __name__ == '__main__':
    unittest.main()
