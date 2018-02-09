""" Tests for nuclide.py. """

import unittest
import xml.etree.ElementTree as ET

from opendeplete import nuclide


class TestNuclide(unittest.TestCase):
    """ Tests for the nuclide class. """

    def test_n_decay_modes(self):
        """ Test the decay mode count parameter. """

        nuc = nuclide.Nuclide()

        nuc.decay_modes = [
            nuclide.DecayTuple("beta1", "a", 0.5),
            nuclide.DecayTuple("beta2", "b", 0.3),
            nuclide.DecayTuple("beta3", "c", 0.2)
        ]

        self.assertEqual(nuc.n_decay_modes, 3)

    def test_n_reaction_paths(self):
        """ Test the reaction path count parameter. """

        nuc = nuclide.Nuclide()

        nuc.reactions = [
            nuclide.ReactionTuple("(n,2n)", "a", 0.0, 1.0),
            nuclide.ReactionTuple("(n,3n)", "b", 0.0, 1.0),
            nuclide.ReactionTuple("(n,4n)", "c", 0.0, 1.0)
        ]

        self.assertEqual(nuc.n_reaction_paths, 3)

    def test_xml_read(self):
        """Test reading nuclide data from an XML element."""

        data = """
<nuclide_table name="U235" reactions="2">
  <decay_type type="sf" target="U235" branching_ratio="7.2e-11"/>
  <decay_type type="alpha" target="Th231" branching_ratio="0.999999999928"/>
  <reaction_type type="(n,2n)" target="U234" Q="-5297781.0"/>
  <reaction_type type="(n,3n)" target="U233" Q="-12142300.0"/>
  <reaction_type type="(n,4n)" target="U232" Q="-17885600.0"/>
  <reaction_type type="(n,gamma)" target="U236" Q="6545200.0"/>
  <reaction_type type="fission" Q="193405400.0"/>
  <neutron_fission_yields>
    <energies>0.0253</energies>
    <fission_yields energy="0.0253">
      <products>Te134 Zr100 Xe138</products>
      <data>0.062155 0.0497641 0.0481413</data>
    </fission_yields>
  </neutron_fission_yields>
</nuclide_table>
        """

        element = ET.fromstring(data)
        u235 = nuclide.Nuclide.xml_read(element)

        self.assertEqual(u235.decay_modes, [
            nuclide.DecayTuple('sf', 'U235', 7.2e-11),
            nuclide.DecayTuple('alpha', 'Th231', 1 - 7.2e-11)
        ])
        self.assertEqual(u235.reactions, [
            nuclide.ReactionTuple('(n,2n)', 'U234', -5297781.0, 1.0),
            nuclide.ReactionTuple('(n,3n)', 'U233', -12142300.0, 1.0),
            nuclide.ReactionTuple('(n,4n)', 'U232', -17885600.0, 1.0),
            nuclide.ReactionTuple('(n,gamma)', 'U236', 6545200.0, 1.0),
            nuclide.ReactionTuple('fission', None, 193405400.0, 1.0),
        ])
        self.assertEqual(u235.yield_energies, [0.0253])
        self.assertEqual(u235.yield_data, {
            0.0253: [('Te134', 0.062155), ('Zr100', 0.0497641),
                     ('Xe138', 0.0481413)]
        })

    def test_xml_write(self):
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
        element = C.xml_write()

        self.assertEqual(element.get("half_life"), "0.123")

        decay_elems = element.findall("decay_type")
        self.assertEqual(len(decay_elems), 2)
        self.assertEqual(decay_elems[0].get("type"), "beta-")
        self.assertEqual(decay_elems[0].get("target"), "B")
        self.assertEqual(decay_elems[0].get("branching_ratio"), "0.99")
        self.assertEqual(decay_elems[1].get("type"), "alpha")
        self.assertEqual(decay_elems[1].get("target"), "D")
        self.assertEqual(decay_elems[1].get("branching_ratio"), "0.01")

        rx_elems = element.findall("reaction_type")
        self.assertEqual(len(rx_elems), 2)
        self.assertEqual(rx_elems[0].get("type"), "fission")
        self.assertEqual(float(rx_elems[0].get("Q")), 2.0e8)
        self.assertEqual(rx_elems[1].get("type"), "(n,gamma)")
        self.assertEqual(rx_elems[1].get("target"), "A")
        self.assertEqual(float(rx_elems[1].get("Q")), 0.0)

        self.assertIsNotNone(element.find('neutron_fission_yields'))


if __name__ == '__main__':
    unittest.main()
