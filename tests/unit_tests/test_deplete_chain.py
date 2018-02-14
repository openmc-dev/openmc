"""Tests for openmc.deplete.Chain class."""

from collections.abc import Mapping
import os
from pathlib import Path

import numpy as np
from openmc.deplete import comm, Chain, reaction_rates, nuclide


_test_filename = str(Path(__file__).parents[2] / 'chains' / 'chain_test.xml')


def test_init():
    """Test depletion chain initialization."""
    dep = Chain()

    assert isinstance(dep.nuclides, list)
    assert isinstance(dep.nuclide_dict, Mapping)
    assert isinstance(dep.react_to_ind, Mapping)


def test_n_nuclides():
    """Test depletion chain n_nuclides parameter."""
    dep = Chain()
    dep.nuclides = ["NucA", "NucB", "NucC"]

    assert dep.n_nuclides == 3


def test_from_endf():
    """Test depletion chain building from ENDF. Empty at the moment until we figure
    out a good way to unit-test this."""
    pass


def test_from_xml():
    """Read chain_test.xml and ensure all values are correct."""
    # Unfortunately, this routine touches a lot of the code, but most of
    # the components external to depletion_chain.py are simple storage
    # types.

    dep = Chain.from_xml(_test_filename)

    # Basic checks
    assert dep.n_nuclides == 3

    # A tests
    nuc = dep.nuclides[dep.nuclide_dict["A"]]

    assert nuc.name == "A"
    assert nuc.half_life == 2.36520E+04
    assert nuc.n_decay_modes == 2
    modes = nuc.decay_modes
    assert [m.target for m in modes] == ["B", "C"]
    assert [m.type for m in modes] == ["beta1", "beta2"]
    assert [m.branching_ratio for m in modes] == [0.6, 0.4]
    assert nuc.n_reaction_paths == 1
    assert [r.target for r in nuc.reactions] == ["C"]
    assert [r.type for r in nuc.reactions] == ["(n,gamma)"]
    assert [r.branching_ratio for r in nuc.reactions] == [1.0]

    # B tests
    nuc = dep.nuclides[dep.nuclide_dict["B"]]

    assert nuc.name == "B"
    assert nuc.half_life == 3.29040E+04
    assert nuc.n_decay_modes == 1
    modes = nuc.decay_modes
    assert [m.target for m in modes] == ["A"]
    assert [m.type for m in modes] == ["beta"]
    assert [m.branching_ratio for m in modes] == [1.0]
    assert nuc.n_reaction_paths == 1
    assert [r.target for r in nuc.reactions] == ["C"]
    assert [r.type for r in nuc.reactions] == ["(n,gamma)"]
    assert [r.branching_ratio for r in nuc.reactions] == [1.0]

    # C tests
    nuc = dep.nuclides[dep.nuclide_dict["C"]]

    assert nuc.name == "C"
    assert nuc.n_decay_modes == 0
    assert nuc.n_reaction_paths == 3
    assert [r.target for r in nuc.reactions] == [None, "A", "B"]
    assert [r.type for r in nuc.reactions] == ["fission", "(n,gamma)", "(n,gamma)"]
    assert [r.branching_ratio for r in nuc.reactions] == [1.0, 0.7, 0.3]

    # Yield tests
    assert nuc.yield_energies == [0.0253]
    assert list(nuc.yield_data) == [0.0253]
    assert nuc.yield_data[0.0253] == [("A", 0.0292737), ("B", 0.002566345)]


def test_export_to_xml(run_in_tmpdir):
    """Test writing a depletion chain to XML."""

    # Prevent different MPI ranks from conflicting
    filename = 'test{}.xml'.format(comm.rank)

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

    chain = Chain()
    chain.nuclides = [A, B, C]
    chain.export_to_xml(filename)

    original = open(_test_filename, 'r').read()
    chain_xml = open(filename, 'r').read()
    assert original == chain_xml


def test_form_matrix():
    """ Using chain_test, and a dummy reaction rate, compute the matrix. """
    # Relies on test_from_xml passing.

    dep = Chain.from_xml(_test_filename)

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

    assert mat[0, 0] == mat00
    assert mat[1, 0] == mat10
    assert mat[2, 0] == mat20
    assert mat[0, 1] == mat01
    assert mat[1, 1] == mat11
    assert mat[2, 1] == mat21
    assert mat[0, 2] == mat02
    assert mat[1, 2] == mat12
    assert mat[2, 2] == mat22


def test_nuc_by_ind():
    """ Test nuc_by_ind converter function. """
    dep = Chain()
    dep.nuclides = ["NucA", "NucB", "NucC"]
    dep.nuclide_dict = {"NucA" : 0, "NucB" : 1, "NucC" : 2}

    assert "NucA" == dep.nuc_by_ind("NucA")
    assert "NucB" == dep.nuc_by_ind("NucB")
    assert "NucC" == dep.nuc_by_ind("NucC")
