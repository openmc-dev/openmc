"""Tests for openmc.deplete.Chain class."""

from collections.abc import Mapping
import os
from pathlib import Path

import numpy as np
from openmc.data import zam, ATOMIC_SYMBOL
from openmc.deplete import comm, Chain, reaction_rates, nuclide
import pytest

from tests import cdtemp

_ENDF_DATA = Path(os.environ['OPENMC_ENDF_DATA'])

_TEST_CHAIN = """\
<depletion_chain>
  <nuclide name="A" half_life="23652.0" decay_modes="2" decay_energy="0.0" reactions="1">
    <decay type="beta1" target="B" branching_ratio="0.6"/>
    <decay type="beta2" target="C" branching_ratio="0.4"/>
    <reaction type="(n,gamma)" Q="0.0" target="C"/>
  </nuclide>
  <nuclide name="B" half_life="32904.0" decay_modes="1" decay_energy="0.0" reactions="1">
    <decay type="beta" target="A" branching_ratio="1.0"/>
    <reaction type="(n,gamma)" Q="0.0" target="C"/>
  </nuclide>
  <nuclide name="C" reactions="3">
    <reaction type="fission" Q="200000000.0"/>
    <reaction type="(n,gamma)" Q="0.0" target="A" branching_ratio="0.7"/>
    <reaction type="(n,gamma)" Q="0.0" target="B" branching_ratio="0.3"/>
    <neutron_fission_yields>
      <energies>0.0253</energies>
      <fission_yields energy="0.0253">
        <products>A B</products>
        <data>0.0292737 0.002566345</data>
      </fission_yields>
    </neutron_fission_yields>
  </nuclide>
</depletion_chain>
"""


@pytest.fixture(scope='module')
def simple_chain():
    with cdtemp():
        with open('chain_test.xml', 'w') as fh:
            fh.write(_TEST_CHAIN)
        yield Chain.from_xml('chain_test.xml')


def test_init():
    """Test depletion chain initialization."""
    chain = Chain()

    assert isinstance(chain.nuclides, list)
    assert isinstance(chain.nuclide_dict, Mapping)


def test_len():
    """Test depletion chain length."""
    chain = Chain()
    chain.nuclides = ["NucA", "NucB", "NucC"]

    assert len(chain) == 3


def test_from_endf():
    """Test depletion chain building from ENDF files"""
    decay_data = (_ENDF_DATA / 'decay').glob('*.endf')
    fpy_data = (_ENDF_DATA / 'nfy').glob('*.endf')
    neutron_data = (_ENDF_DATA / 'neutrons').glob('*.endf')
    chain = Chain.from_endf(decay_data, fpy_data, neutron_data)

    assert len(chain) == len(chain.nuclides) == len(chain.nuclide_dict) == 3820
    for nuc in chain.nuclides:
        assert nuc == chain[nuc.name]


def test_from_xml(simple_chain):
    """Read chain_test.xml and ensure all values are correct."""
    # Unfortunately, this routine touches a lot of the code, but most of
    # the components external to depletion_chain.py are simple storage
    # types.

    chain = simple_chain

    # Basic checks
    assert len(chain) == 3

    # A tests
    nuc = chain["A"]

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
    nuc = chain["B"]

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
    nuc = chain["C"]

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

    chain_xml = open(filename, 'r').read()
    assert _TEST_CHAIN == chain_xml


def test_form_matrix(simple_chain):
    """ Using chain_test, and a dummy reaction rate, compute the matrix. """
    # Relies on test_from_xml passing.

    chain = simple_chain

    mats = ["10000", "10001"]
    nuclides = ["A", "B", "C"]

    react = reaction_rates.ReactionRates(mats, nuclides, chain.reactions)

    react.set("10000", "C", "fission", 1.0)
    react.set("10000", "A", "(n,gamma)", 2.0)
    react.set("10000", "B", "(n,gamma)", 3.0)
    react.set("10000", "C", "(n,gamma)", 4.0)

    mat = chain.form_matrix(react[0, :, :])
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


def test_getitem():
    """ Test nuc_by_ind converter function. """
    chain = Chain()
    chain.nuclides = ["NucA", "NucB", "NucC"]
    chain.nuclide_dict = {nuc: chain.nuclides.index(nuc)
                        for nuc in chain.nuclides}

    assert "NucA" == chain["NucA"]
    assert "NucB" == chain["NucB"]
    assert "NucC" == chain["NucC"]
