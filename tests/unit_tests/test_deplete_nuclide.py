"""Tests for the openmc.deplete.Nuclide class."""

import xml.etree.ElementTree as ET

import numpy
import pytest
from openmc.deplete import nuclide


def test_n_decay_modes():
    """ Test the decay mode count parameter. """

    nuc = nuclide.Nuclide()

    nuc.decay_modes = [
        nuclide.DecayTuple("beta1", "a", 0.5),
        nuclide.DecayTuple("beta2", "b", 0.3),
        nuclide.DecayTuple("beta3", "c", 0.2)
    ]

    assert nuc.n_decay_modes == 3


def test_n_reaction_paths():
    """ Test the reaction path count parameter. """

    nuc = nuclide.Nuclide()

    nuc.reactions = [
        nuclide.ReactionTuple("(n,2n)", "a", 0.0, 1.0),
        nuclide.ReactionTuple("(n,3n)", "b", 0.0, 1.0),
        nuclide.ReactionTuple("(n,4n)", "c", 0.0, 1.0)
    ]

    assert nuc.n_reaction_paths == 3


def test_from_xml():
    """Test reading nuclide data from an XML element."""

    data = """
<nuclide name="U235" reactions="2">
  <decay type="sf" target="U235" branching_ratio="7.2e-11"/>
  <decay type="alpha" target="Th231" branching_ratio="0.999999999928"/>
  <reaction type="(n,2n)" target="U234" Q="-5297781.0"/>
  <reaction type="(n,3n)" target="U233" Q="-12142300.0"/>
  <reaction type="(n,4n)" target="U232" Q="-17885600.0"/>
  <reaction type="(n,gamma)" target="U236" Q="6545200.0"/>
  <reaction type="fission" Q="193405400.0"/>
  <neutron_fission_yields>
    <energies>0.0253</energies>
    <fission_yields energy="0.0253">
      <products>Te134 Zr100 Xe138</products>
      <data>0.062155 0.0497641 0.0481413</data>
    </fission_yields>
  </neutron_fission_yields>
</nuclide>
    """

    element = ET.fromstring(data)
    u235 = nuclide.Nuclide.from_xml(element)

    assert u235.decay_modes == [
        nuclide.DecayTuple('sf', 'U235', 7.2e-11),
        nuclide.DecayTuple('alpha', 'Th231', 1 - 7.2e-11)
    ]
    assert u235.reactions == [
        nuclide.ReactionTuple('(n,2n)', 'U234', -5297781.0, 1.0),
        nuclide.ReactionTuple('(n,3n)', 'U233', -12142300.0, 1.0),
        nuclide.ReactionTuple('(n,4n)', 'U232', -17885600.0, 1.0),
        nuclide.ReactionTuple('(n,gamma)', 'U236', 6545200.0, 1.0),
        nuclide.ReactionTuple('fission', None, 193405400.0, 1.0),
    ]
    expected_yield_data = nuclide.FissionYieldDistribution({
        0.0253: {"Xe138": 0.0481413, "Zr100": 0.0497641, "Te134": 0.062155}})
    assert u235.yield_data == expected_yield_data
    # test accessing the yield energies through the FissionYieldDistribution
    assert u235.yield_energies == (0.0253,)
    assert u235.yield_energies is u235.yield_data.energies
    with pytest.raises(AttributeError):  # not settable
        u235.yield_energies = [0.0253, 5e5]


def test_to_xml_element():
    """Test writing nuclide data to an XML element."""

    C = nuclide.Nuclide("C")
    C.half_life = 0.123
    C.decay_modes = [
        nuclide.DecayTuple('beta-', 'B', 0.99),
        nuclide.DecayTuple('alpha', 'D', 0.01)
    ]
    C.reactions = [
        nuclide.ReactionTuple('fission', None, 2.0e8, 1.0),
        nuclide.ReactionTuple('(n,gamma)', 'A', 0.0, 1.0)
    ]
    C.yield_data = nuclide.FissionYieldDistribution(
        {0.0253: {"A": 0.0292737, "B": 0.002566345}})
    element = C.to_xml_element()

    assert element.get("half_life") == "0.123"

    decay_elems = element.findall("decay")
    assert len(decay_elems) == 2
    assert decay_elems[0].get("type") == "beta-"
    assert decay_elems[0].get("target") == "B"
    assert decay_elems[0].get("branching_ratio") == "0.99"
    assert decay_elems[1].get("type") == "alpha"
    assert decay_elems[1].get("target") == "D"
    assert decay_elems[1].get("branching_ratio") == "0.01"

    rx_elems = element.findall("reaction")
    assert len(rx_elems) == 2
    assert rx_elems[0].get("type") == "fission"
    assert float(rx_elems[0].get("Q")) == 2.0e8
    assert rx_elems[1].get("type") == "(n,gamma)"
    assert rx_elems[1].get("target") == "A"
    assert float(rx_elems[1].get("Q")) == 0.0

    assert element.find('neutron_fission_yields') is not None


def test_fission_yield_distribution():
    """Test an energy-dependent yield distribution"""
    yield_dict = {
        0.0253: {"Xe135": 7.85e-4, "Gd155": 4.08e-12, "Sm149": 1.71e-12},
        1.40e7: {"Xe135": 4.54e-3, "Gd155": 5.83e-8, "Sm149": 2.69e-8},
        5.00e5: {"Xe135": 1.12e-3, "Gd155": 1.32e-12},  # drop Sm149
    }
    yield_dist = nuclide.FissionYieldDistribution(yield_dict)
    assert len(yield_dist) == len(yield_dict)
    assert yield_dist.energies == tuple(sorted(yield_dict.keys()))
    for exp_ene, exp_dist in yield_dict.items():
        act_dist = yield_dict[exp_ene]
        for exp_prod, exp_yield in exp_dist.items():
            assert act_dist[exp_prod] == exp_yield
    exp_yield = numpy.array([
        [4.08e-12, 1.71e-12, 7.85e-4],
        [1.32e-12, 0.0, 1.12e-3],
        [5.83e-8, 2.69e-8, 4.54e-3]])
    assert numpy.array_equal(yield_dist.yield_matrix, exp_yield)

    # Test the operations / special methods for fission yield
    orig_yields = yield_dist[0.0253]
    assert len(orig_yields) == len(yield_dict[0.0253])
    for key, value in yield_dict[0.0253].items():
        assert key in orig_yields
        assert orig_yields[key] == value
    # __getitem__ return yields as a view into yield matrix
    assert orig_yields.yields.base is yield_dist.yield_matrix

    # Scale and increment fission yields
    mod_yields = orig_yields * 2
    assert numpy.array_equal(orig_yields.yields * 2, mod_yields.yields)
    mod_yields += orig_yields
    assert numpy.array_equal(orig_yields.yields * 3, mod_yields.yields)

    mod_yields = 2.0 * orig_yields
    assert numpy.array_equal(orig_yields.yields * 2, mod_yields.yields)

    mod_yields = numpy.float64(2.0) * orig_yields
    assert numpy.array_equal(orig_yields.yields * 2, mod_yields.yields)

    # Failure modes for adding, multiplying yields
    similar = numpy.empty_like(orig_yields.yields)
    with pytest.raises(TypeError):
        orig_yields + similar
    with pytest.raises(TypeError):
        similar + orig_yields
    with pytest.raises(TypeError):
        orig_yields += similar
    with pytest.raises(TypeError):
        orig_yields * similar
    with pytest.raises(TypeError):
        similar * orig_yields
    with pytest.raises(TypeError):
        orig_yields *= similar

    # Test restriction of fission products
    strict_restrict = yield_dist.restrict_products(["Xe135", "Sm149"])
    with_extras = yield_dist.restrict_products(
        ["Xe135", "Sm149", "H1", "U235"])

    assert strict_restrict.products == ("Sm149", "Xe135")
    assert strict_restrict.energies == yield_dist.energies
    assert with_extras.products == ("Sm149", "Xe135")
    assert with_extras.energies == yield_dist.energies
    for ene, new_yields in strict_restrict.items():
        for product in strict_restrict.products:
            assert new_yields[product] == yield_dist[ene][product]
            assert with_extras[ene][product] == yield_dist[ene][product]

    assert yield_dist.restrict_products(["U235"]) is None

def test_validate():

    nuc = nuclide.Nuclide()
    nuc.name = "Test"

    # decay modes: type, target, branching_ratio

    nuc.decay_modes = [
        nuclide.DecayTuple("type 0", "0", 0.5),
        nuclide.DecayTuple("type 1", "1", 0.5),
    ]

    # reactions: type, target, Q, branching_ratio
    nuc.reactions = [
        nuclide.ReactionTuple("0", "0", 1000, 0.3),
        nuclide.ReactionTuple("0", "1", 1000, 0.3),
        nuclide.ReactionTuple("1", "2", 1000, 1.0),
        nuclide.ReactionTuple("0", "3", 1000, 0.4),
    ]

    # fission yields

    nuc.yield_data = {
        0.0253: {"0": 1.5, "1": 0.5},
        1e6: {"0": 1.5, "1": 0.5},
    }

    # nuclide is good and should have no warnings raise
    with pytest.warns(None) as record:
        assert nuc.validate(strict=True, quiet=False, tolerance=0.0)
    assert len(record) == 0

    # invalidate decay modes
    decay = nuc.decay_modes.pop()
    with pytest.raises(ValueError, match="decay mode"):
        nuc.validate(strict=True, quiet=False, tolerance=0.0)

    with pytest.warns(UserWarning) as record:
        assert not nuc.validate(strict=False, quiet=False, tolerance=0.0)
        assert not nuc.validate(strict=False, quiet=True, tolerance=0.0)
    assert len(record) == 1
    assert "decay mode" in record[0].message.args[0]

    # restore decay modes, invalidate reactions
    nuc.decay_modes.append(decay)
    reaction = nuc.reactions.pop()

    with pytest.raises(ValueError, match="0 reaction"):
        nuc.validate(strict=True, quiet=False, tolerance=0.0)

    with pytest.warns(UserWarning) as record:
        assert not nuc.validate(strict=False, quiet=False, tolerance=0.0)
        assert not nuc.validate(strict=False, quiet=True, tolerance=0.0)
    assert len(record) == 1
    assert "0 reaction" in record[0].message.args[0]

    # restore reactions, invalidate fission yields
    nuc.reactions.append(reaction)
    nuc.yield_data[1e6].yields *= 2

    with pytest.raises(ValueError, match=r"fission yields.*1\.0*e"):
        nuc.validate(strict=True, quiet=False, tolerance=0.0)

    with pytest.warns(UserWarning) as record:
        assert not nuc.validate(strict=False, quiet=False, tolerance=0.0)
        assert not nuc.validate(strict=False, quiet=True, tolerance=0.0)
    assert len(record) == 1
    assert "1.0" in record[0].message.args[0]

    # invalidate everything, check that error is raised at decay modes

    decay = nuc.decay_modes.pop()
    reaction = nuc.reactions.pop()

    with pytest.raises(ValueError, match="decay mode"):
        nuc.validate(strict=True, quiet=False, tolerance=0.0)

    # check for warnings
    # should be one warning for decay modes, reactions, fission yields

    with pytest.warns(UserWarning) as record:
        assert not nuc.validate(strict=False, quiet=False, tolerance=0.0)
        assert not nuc.validate(strict=False, quiet=True, tolerance=0.0)
    assert len(record) == 3
    assert "decay mode" in record[0].message.args[0]
    assert "0 reaction" in record[1].message.args[0]
    assert "1.0" in record[2].message.args[0]
