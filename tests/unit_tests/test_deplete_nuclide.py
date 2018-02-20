"""Tests for the openmc.deplete.Nuclide class."""

import xml.etree.ElementTree as ET

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
    assert u235.yield_energies == [0.0253]
    assert u235.yield_data == {
        0.0253: [('Te134', 0.062155), ('Zr100', 0.0497641),
                 ('Xe138', 0.0481413)]
    }


def test_to_xml_element():
    """Test writing nuclide data to an XML element."""

    C = nuclide.Nuclide()
    C.name = "C"
    C.half_life = 0.123
    C.decay_modes = [
        nuclide.DecayTuple('beta-', 'B', 0.99),
        nuclide.DecayTuple('alpha', 'D', 0.01)
    ]
    C.reactions = [
        nuclide.ReactionTuple('fission', None, 2.0e8, 1.0),
        nuclide.ReactionTuple('(n,gamma)', 'A', 0.0, 1.0)
    ]
    C.yield_energies = [0.0253]
    C.yield_data = {0.0253: [("A", 0.0292737), ("B", 0.002566345)]}
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
