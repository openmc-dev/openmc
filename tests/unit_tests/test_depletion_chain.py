""" Tests for depletion_chain.py"""

from collections import OrderedDict
import os
import unittest
from pathlib import Path

import numpy as np
from openmc.deplete import comm, depletion_chain, reaction_rates, nuclide


_test_filename = str(Path(__file__).parents[2] / 'chains' / 'chain_test.xml')


class TestDepletionChain(unittest.TestCase):
    """ Tests for DepletionChain class."""

    def test__init__(self):
        """ Test depletion chain initialization."""
        dep = depletion_chain.DepletionChain()

        self.assertIsInstance(dep.nuclides, list)
        self.assertIsInstance(dep.nuclide_dict, OrderedDict)
        self.assertIsInstance(dep.react_to_ind, OrderedDict)

    def test_n_nuclides(self):
        """ Test depletion chain n_nuclides parameter. """
        dep = depletion_chain.DepletionChain()

        dep.nuclides = ["NucA", "NucB", "NucC"]

        self.assertEqual(dep.n_nuclides, 3)

    def test_from_endf(self):
        """Test depletion chain building from ENDF. Empty at the moment until we figure
        out a good way to unit-test this."""
        pass

    def test_xml_read(self):
        """ Read chain_test.xml and ensure all values are correct. """
        # Unfortunately, this routine touches a lot of the code, but most of
        # the components external to depletion_chain.py are simple storage
        # types.

        dep = depletion_chain.DepletionChain.xml_read(_test_filename)

        # Basic checks
        self.assertEqual(dep.n_nuclides, 3)

        # A tests
        nuc = dep.nuclides[dep.nuclide_dict["A"]]

        self.assertEqual(nuc.name, "A")
        self.assertEqual(nuc.half_life, 2.36520E+04)
        self.assertEqual(nuc.n_decay_modes, 2)
        modes = nuc.decay_modes
        self.assertEqual([m.target for m in modes], ["B", "C"])
        self.assertEqual([m.type for m in modes], ["beta1", "beta2"])
        self.assertEqual([m.branching_ratio for m in modes], [0.6, 0.4])
        self.assertEqual(nuc.n_reaction_paths, 1)
        self.assertEqual([r.target for r in nuc.reactions], ["C"])
        self.assertEqual([r.type for r in nuc.reactions], ["(n,gamma)"])
        self.assertEqual([r.branching_ratio for r in nuc.reactions], [1.0])

        # B tests
        nuc = dep.nuclides[dep.nuclide_dict["B"]]

        self.assertEqual(nuc.name, "B")
        self.assertEqual(nuc.half_life, 3.29040E+04)
        self.assertEqual(nuc.n_decay_modes, 1)
        modes = nuc.decay_modes
        self.assertEqual([m.target for m in modes], ["A"])
        self.assertEqual([m.type for m in modes], ["beta"])
        self.assertEqual([m.branching_ratio for m in modes], [1.0])
        self.assertEqual(nuc.n_reaction_paths, 1)
        self.assertEqual([r.target for r in nuc.reactions], ["C"])
        self.assertEqual([r.type for r in nuc.reactions], ["(n,gamma)"])
        self.assertEqual([r.branching_ratio for r in nuc.reactions], [1.0])

        # C tests
        nuc = dep.nuclides[dep.nuclide_dict["C"]]

        self.assertEqual(nuc.name, "C")
        self.assertEqual(nuc.n_decay_modes, 0)
        self.assertEqual(nuc.n_reaction_paths, 3)
        self.assertEqual([r.target for r in nuc.reactions], [None, "A", "B"])
        self.assertEqual([r.type for r in nuc.reactions], ["fission", "(n,gamma)", "(n,gamma)"])
        self.assertEqual([r.branching_ratio for r in nuc.reactions], [1.0, 0.7, 0.3])

        # Yield tests
        self.assertEqual(nuc.yield_energies, [0.0253])
        self.assertEqual(list(nuc.yield_data.keys()), [0.0253])
        self.assertEqual(nuc.yield_data[0.0253],
                         [("A", 0.0292737), ("B", 0.002566345)])

    def test_xml_write(self):
        """Test writing a depletion chain to XML."""

        # Prevent different MPI ranks from conflicting
        filename = 'test%u.xml' % comm.rank

        A = nuclide.Nuclide()
        A.name = "A"
        A.half_life = 2.36520e4
        A.decay_modes = [
            nuclide.DecayTuple("beta1", "B", 0.6),
            nuclide.DecayTuple("beta2", "C", 0.4)
        ]
        A.reactions = [nuclide.ReactionTuple("(n,gamma)", "C", 0.0, 1.0)]

        B = nuclide.Nuclide()
        B.name = "B"
        B.half_life = 3.29040e4
        B.decay_modes = [nuclide.DecayTuple("beta", "A", 1.0)]
        B.reactions = [nuclide.ReactionTuple("(n,gamma)", "C", 0.0, 1.0)]

        C = nuclide.Nuclide()
        C.name = "C"
        C.reactions = [
            nuclide.ReactionTuple("fission", None, 2.0e8, 1.0),
            nuclide.ReactionTuple("(n,gamma)", "A", 0.0, 0.7),
            nuclide.ReactionTuple("(n,gamma)", "B", 0.0, 0.3)
        ]
        C.yield_energies = [0.0253]
        C.yield_data = {0.0253: [("A", 0.0292737), ("B", 0.002566345)]}

        chain = depletion_chain.DepletionChain()
        chain.nuclides = [A, B, C]
        chain.xml_write(filename)

        original = open(_test_filename, 'r').read()
        chain_xml = open(filename, 'r').read()
        self.assertEqual(original, chain_xml)

        os.remove(filename)

    def test_form_matrix(self):
        """ Using chain_test, and a dummy reaction rate, compute the matrix. """
        # Relies on test_xml_read passing.

        dep = depletion_chain.DepletionChain.xml_read(_test_filename)

        cell_ind = {"10000": 0, "10001": 1}
        nuc_ind = {"A": 0, "B": 1, "C": 2}
        react_ind = dep.react_to_ind

        react = reaction_rates.ReactionRates(cell_ind, nuc_ind, react_ind)

        dep.nuc_to_react_ind = nuc_ind

        react["10000", "C", "fission"] = 1.0
        react["10000", "A", "(n,gamma)"] = 2.0
        react["10000", "B", "(n,gamma)"] = 3.0
        react["10000", "C", "(n,gamma)"] = 4.0

        mat = dep.form_matrix(react[0, :, :])
        # Loss A, decay, (n, gamma)
        mat00 = -np.log(2) / 2.36520E+04 - 2
        # A -> B, decay, 0.6 branching ratio
        mat10 = np.log(2) / 2.36520E+04 * 0.6
        # A -> C, decay, 0.4 branching ratio + (n,gamma)
        mat20 = np.log(2) / 2.36520E+04 * 0.4 + 2

        # B -> A, decay, 1.0 branching ratio
        mat01 = np.log(2)/3.29040E+04
        # Loss B, decay, (n, gamma)
        mat11 = -np.log(2)/3.29040E+04 - 3
        # B -> C, (n, gamma)
        mat21 = 3

        # C -> A fission, (n, gamma)
        mat02 = 0.0292737 * 1.0 + 4.0 * 0.7
        # C -> B fission, (n, gamma)
        mat12 = 0.002566345 * 1.0 + 4.0 * 0.3
        # Loss C, fission, (n, gamma)
        mat22 = -1.0 - 4.0

        self.assertEqual(mat[0, 0], mat00)
        self.assertEqual(mat[1, 0], mat10)
        self.assertEqual(mat[2, 0], mat20)
        self.assertEqual(mat[0, 1], mat01)
        self.assertEqual(mat[1, 1], mat11)
        self.assertEqual(mat[2, 1], mat21)
        self.assertEqual(mat[0, 2], mat02)
        self.assertEqual(mat[1, 2], mat12)
        self.assertEqual(mat[2, 2], mat22)

    def test_nuc_by_ind(self):
        """ Test nuc_by_ind converter function. """
        dep = depletion_chain.DepletionChain()

        dep.nuclides = ["NucA", "NucB", "NucC"]
        dep.nuclide_dict = {"NucA" : 0, "NucB" : 1, "NucC" : 2}

        self.assertEqual("NucA", dep.nuc_by_ind("NucA"))
        self.assertEqual("NucB", dep.nuc_by_ind("NucB"))
        self.assertEqual("NucC", dep.nuc_by_ind("NucC"))

if __name__ == '__main__':
    unittest.main()
