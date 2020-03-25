"""Tests for openmc.deplete.Chain class."""

from collections.abc import Mapping
import os
from pathlib import Path
from itertools import product

import numpy as np
from openmc.deplete import comm, Chain, reaction_rates, nuclide, cram
import pytest

from tests import cdtemp

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
    endf_data = Path(os.environ['OPENMC_ENDF_DATA'])
    decay_data = (endf_data / 'decay').glob('*.endf')
    fpy_data = (endf_data / 'nfy').glob('*.endf')
    neutron_data = (endf_data / 'neutrons').glob('*.endf')
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
    assert nuc.yield_energies == (0.0253,)
    assert list(nuc.yield_data) == [0.0253]
    assert nuc.yield_data[0.0253].products == ("A", "B")
    assert (nuc.yield_data[0.0253].yields == [0.0292737, 0.002566345]).all()


def test_export_to_xml(run_in_tmpdir):
    """Test writing a depletion chain to XML."""

    # Prevent different MPI ranks from conflicting
    filename = 'test{}.xml'.format(comm.rank)

    A = nuclide.Nuclide("A")
    A.half_life = 2.36520e4
    A.decay_modes = [
        nuclide.DecayTuple("beta1", "B", 0.6),
        nuclide.DecayTuple("beta2", "C", 0.4)
    ]
    A.reactions = [nuclide.ReactionTuple("(n,gamma)", "C", 0.0, 1.0)]

    B = nuclide.Nuclide("B")
    B.half_life = 3.29040e4
    B.decay_modes = [nuclide.DecayTuple("beta", "A", 1.0)]
    B.reactions = [nuclide.ReactionTuple("(n,gamma)", "C", 0.0, 1.0)]

    C = nuclide.Nuclide("C")
    C.reactions = [
        nuclide.ReactionTuple("fission", None, 2.0e8, 1.0),
        nuclide.ReactionTuple("(n,gamma)", "A", 0.0, 0.7),
        nuclide.ReactionTuple("(n,gamma)", "B", 0.0, 0.3)
    ]
    C.yield_data = nuclide.FissionYieldDistribution({
        0.0253: {"A": 0.0292737, "B": 0.002566345}})

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

    # Pass equivalent fission yields directly
    # Ensure identical matrix is formed
    f_yields = {"C": {"A": 0.0292737, "B": 0.002566345}}
    new_mat = chain.form_matrix(react[0], f_yields)
    for r, c in product(range(3), range(3)):
        assert new_mat[r, c] == mat[r, c]


def test_getitem():
    """ Test nuc_by_ind converter function. """
    chain = Chain()
    chain.nuclides = ["NucA", "NucB", "NucC"]
    chain.nuclide_dict = {nuc: chain.nuclides.index(nuc)
                        for nuc in chain.nuclides}

    assert "NucA" == chain["NucA"]
    assert "NucB" == chain["NucB"]
    assert "NucC" == chain["NucC"]


def test_set_fiss_q():
    """Make sure new fission q values can be set on the chain"""
    new_q = {"U235": 2.0E8, "U238": 2.0E8, "U234": 5.0E7}
    chain_file = Path(__file__).parents[1] / "chain_simple.xml"
    mod_chain = Chain.from_xml(chain_file, new_q)
    for name, q in new_q.items():
        chain_nuc = mod_chain[name]
        for rx in chain_nuc.reactions:
            if rx.type == 'fission':
                assert rx.Q == q


def test_get_set_chain_br(simple_chain):
    """Test minor modifications to capture branch ratios"""
    expected = {"C": {"A": 0.7, "B": 0.3}}
    assert simple_chain.get_branch_ratios() == expected

    # safely modify
    new_chain = Chain.from_xml("chain_test.xml")
    new_br = {"C": {"A": 0.5, "B": 0.5}, "A": {"C": 0.99, "B": 0.01}}
    new_chain.set_branch_ratios(new_br)
    assert new_chain.get_branch_ratios() == new_br

    # write, re-read
    new_chain.export_to_xml("chain_mod.xml")
    assert Chain.from_xml("chain_mod.xml").get_branch_ratios() == new_br

    # Test non-strict [warn, not error] setting
    bad_br = {"B": {"X": 0.6, "A": 0.4}, "X": {"A": 0.5, "C": 0.5}}
    bad_br.update(new_br)
    new_chain.set_branch_ratios(bad_br, strict=False)
    assert new_chain.get_branch_ratios() == new_br

    # Ensure capture reactions are removed
    rem_br = {"A": {"C": 1.0}}
    new_chain.set_branch_ratios(rem_br)
    # A is not in returned dict because there is no branch
    assert "A" not in new_chain.get_branch_ratios()


def test_capture_branch_infer_ground():
    """Ensure the ground state is infered if not given"""
    # Make up a metastable capture transition:
    infer_br = {"Xe135": {"Xe136_m1": 0.5}}
    set_br = {"Xe135": {"Xe136": 0.5, "Xe136_m1": 0.5}}

    chain_file = Path(__file__).parents[1] / "chain_simple.xml"
    chain = Chain.from_xml(chain_file)

    # Create nuclide to be added into the chain
    xe136m = nuclide.Nuclide("Xe136_m1")

    chain.nuclides.append(xe136m)
    chain.nuclide_dict[xe136m.name] = len(chain.nuclides) - 1

    chain.set_branch_ratios(infer_br, "(n,gamma)")

    assert chain.get_branch_ratios("(n,gamma)") == set_br


def test_capture_branch_no_rxn():
    """Ensure capture reactions that don't exist aren't created"""
    u4br = {"U234": {"U235": 0.5, "U235_m1": 0.5}}

    chain_file = Path(__file__).parents[1] / "chain_simple.xml"
    chain = Chain.from_xml(chain_file)

    u5m = nuclide.Nuclide("U235_m1")

    chain.nuclides.append(u5m)
    chain.nuclide_dict[u5m.name] = len(chain.nuclides) - 1

    with pytest.raises(AttributeError, match="U234"):
        chain.set_branch_ratios(u4br)


def test_capture_branch_failures(simple_chain):
    """Test failure modes for setting capture branch ratios"""

    # Parent isotope not present
    br = {"X": {"A": 0.6, "B": 0.7}}
    with pytest.raises(KeyError, match="X"):
        simple_chain.set_branch_ratios(br)

    # Product isotope not present
    br = {"C": {"X": 0.4, "A": 0.2, "B": 0.4}}
    with pytest.raises(KeyError, match="X"):
        simple_chain.set_branch_ratios(br)

    # Sum of ratios > 1.0
    br = {"C": {"A": 1.0, "B": 1.0}}
    with pytest.raises(ValueError, match=r"Sum of \(n,gamma\).*for C"):
        simple_chain.set_branch_ratios(br, "(n,gamma)")


def test_set_alpha_branches():
    """Test setting of alpha reaction branching ratios"""
    # Build a mock chain
    chain = Chain()

    parent = nuclide.Nuclide()
    parent.name = "A"

    he4 = nuclide.Nuclide()
    he4.name = "He4"

    ground_tgt = nuclide.Nuclide()
    ground_tgt.name = "B"

    meta_tgt = nuclide.Nuclide()
    meta_tgt.name = "B_m1"

    for ix, nuc in enumerate((parent, ground_tgt, meta_tgt, he4)):
        chain.nuclides.append(nuc)
        chain.nuclide_dict[nuc.name] = ix

    # add reactions to parent
    parent.reactions.append(nuclide.ReactionTuple(
        "(n,a)", ground_tgt.name, 1.0, 0.6))
    parent.reactions.append(nuclide.ReactionTuple(
        "(n,a)", meta_tgt.name, 1.0, 0.4))
    parent.reactions.append(nuclide.ReactionTuple(
        "(n,a)", he4.name, 1.0, 1.0))

    expected_ref = {"A": {"B": 0.6, "B_m1": 0.4}}

    assert chain.get_branch_ratios("(n,a)") == expected_ref

    # alter and check again

    altered = {"A": {"B": 0.5, "B_m1": 0.5}}

    chain.set_branch_ratios(altered, "(n,a)")
    assert chain.get_branch_ratios("(n,a)") == altered

    # make sure that alpha particle still produced
    for r in parent.reactions:
        if r.target == he4.name:
            break
    else:
        raise ValueError("Helium has been removed and should not have been")


def test_simple_fission_yields(simple_chain):
    """Check the default fission yields that can be used to form the matrix
    """
    fission_yields = simple_chain.get_default_fission_yields()
    assert fission_yields == {"C": {"A": 0.0292737, "B": 0.002566345}}


def test_fission_yield_attribute(simple_chain):
    """Test the fission_yields property"""
    thermal_yields = simple_chain.get_default_fission_yields()
    # generate default with property
    assert simple_chain.fission_yields[0] == thermal_yields
    empty_chain = Chain()
    empty_chain.fission_yields = thermal_yields
    assert empty_chain.fission_yields[0] == thermal_yields
    empty_chain.fission_yields = [thermal_yields] * 2
    assert empty_chain.fission_yields[0] == thermal_yields
    assert empty_chain.fission_yields[1] == thermal_yields

    # test failure with deplete function
    # number fission yields != number of materials
    dummy_conc = [[1, 2]] * (len(empty_chain.fission_yields) + 1)
    with pytest.raises(
            ValueError, match="fission yield.*not equal.*compositions"):
        cram.deplete(empty_chain, dummy_conc, None, 0.5)


def test_validate(simple_chain):
    """Test the validate method"""

    # current chain is invalid
    # fission yields do not sum to 2.0
    with pytest.raises(ValueError, match="Nuclide C.*fission yields"):
        simple_chain.validate(strict=True, tolerance=0.0)

    with pytest.warns(UserWarning) as record:
        assert not simple_chain.validate(strict=False, quiet=False, tolerance=0.0)
        assert not simple_chain.validate(strict=False, quiet=True, tolerance=0.0)
    assert len(record) == 1
    assert "Nuclide C" in record[0].message.args[0]

    # Fix fission yields but keep to restore later
    old_yields = simple_chain["C"].yield_data
    simple_chain["C"].yield_data = {0.0253: {"A": 1.4, "B": 0.6}}

    assert simple_chain.validate(strict=True, tolerance=0.0)
    with pytest.warns(None) as record:
        assert simple_chain.validate(strict=False, quiet=False, tolerance=0.0)
    assert len(record) == 0

    # Mess up "earlier" nuclide's reactions
    decay_mode = simple_chain["A"].decay_modes.pop()

    with pytest.raises(ValueError, match="Nuclide A.*decay mode"):
        simple_chain.validate(strict=True, tolerance=0.0)

    # restore old fission yields
    simple_chain["C"].yield_data = old_yields

    with pytest.warns(UserWarning) as record:
        assert not simple_chain.validate(strict=False, quiet=False, tolerance=0.0)
    assert len(record) == 2
    assert "Nuclide A" in record[0].message.args[0]
    assert "Nuclide C" in record[1].message.args[0]

    # restore decay modes
    simple_chain["A"].decay_modes.append(decay_mode)


def test_validate_inputs():
    c = Chain()

    with pytest.raises(TypeError, match="tolerance"):
        c.validate(tolerance=None)

    with pytest.raises(ValueError, match="tolerance"):
        c.validate(tolerance=-1)


@pytest.fixture
def gnd_simple_chain():
    chainfile = Path(__file__).parents[1] / "chain_simple.xml"
    return Chain.from_xml(chainfile)


def test_reduce(gnd_simple_chain):
    ref_U5 = gnd_simple_chain["U235"]
    ref_iodine = gnd_simple_chain["I135"]
    ref_U5_yields = ref_U5.yield_data

    no_depth = gnd_simple_chain.reduce(["U235", "I135"], 0)
    # We should get a chain just containing U235 and I135
    assert len(no_depth) == 2
    assert set(no_depth.reactions) == set(gnd_simple_chain.reactions)

    u5_round0 = no_depth["U235"]
    assert u5_round0.n_decay_modes == ref_U5.n_decay_modes
    for newmode, refmode in zip(u5_round0.decay_modes, ref_U5.decay_modes):
        assert newmode.target is None
        assert newmode.type == refmode.type
        assert newmode.branching_ratio == refmode.branching_ratio

    assert u5_round0.n_reaction_paths == ref_U5.n_reaction_paths
    for newrxn, refrxn in zip(u5_round0.reactions, ref_U5.reactions):
        assert newrxn.target is None
        assert newrxn.type == refrxn.type
        assert newrxn.Q == refrxn.Q
        assert newrxn.branching_ratio == refrxn.branching_ratio

    assert u5_round0.yield_data is not None
    assert u5_round0.yield_data.products == ("I135",)
    assert u5_round0.yield_data.yield_matrix == (
            ref_U5_yields.yield_matrix[:, ref_U5_yields.products.index("I135")]
    )

    bareI5 = no_depth["I135"]
    assert bareI5.n_decay_modes == ref_iodine.n_decay_modes
    for newmode, refmode in zip(bareI5.decay_modes, ref_iodine.decay_modes):
        assert newmode.target is None
        assert newmode.type == refmode.type
        assert newmode.branching_ratio == refmode.branching_ratio

    assert bareI5.n_reaction_paths == ref_iodine.n_reaction_paths
    for newrxn, refrxn in zip(bareI5.reactions, ref_iodine.reactions):
        assert newrxn.target is None
        assert newrxn.type == refrxn.type
        assert newrxn.Q == refrxn.Q
        assert newrxn.branching_ratio == refrxn.branching_ratio

    follow_u5 = gnd_simple_chain.reduce(["U235"], 1)
    u5_round1 = follow_u5["U235"]
    assert u5_round1.decay_modes == ref_U5.decay_modes
    assert u5_round1.reactions == ref_U5.reactions
    assert u5_round1.yield_data is not None
    assert (
        u5_round1.yield_data.yield_matrix == ref_U5_yields.yield_matrix
    ).all()

    # Per the chain_simple.xml
    # I135 -> Xe135 -> Cs135
    # I135 -> Xe136
    # No limit on depth
    iodine_chain = gnd_simple_chain.reduce(["I135"])
    truncated_iodine = gnd_simple_chain.reduce(["I135"], 1)
    assert len(iodine_chain) == 4
    assert len(truncated_iodine) == 3
    assert set(iodine_chain.nuclide_dict) == {
        "I135", "Xe135", "Xe136", "Cs135"}
    assert set(truncated_iodine.nuclide_dict) == {"I135", "Xe135", "Xe136"}
    assert iodine_chain.reactions == ["(n,gamma)"]
    assert iodine_chain["I135"].decay_modes == ref_iodine.decay_modes
    assert iodine_chain["I135"].reactions == ref_iodine.reactions
    for mode in truncated_iodine["Xe135"].decay_modes:
        assert mode.target is None

    # Test that no FissionYieldDistribution is made if there are no
    # fission products
    u5_noyields = gnd_simple_chain.reduce(["U235"], 0)["U235"]
    assert u5_noyields.yield_data is None

    # Check early termination if the eventual full chain
    # is specified by using the iodine isotopes
    new_iodine = gnd_simple_chain.reduce(set(iodine_chain.nuclide_dict))
    assert set(iodine_chain.nuclide_dict) == set(new_iodine.nuclide_dict)

    # Failure if some requested isotopes not in chain

    with pytest.raises(IndexError, match=".*not found.*Xx999"):
        gnd_simple_chain.reduce(["U235", "Xx999"])
