"""Nuclide module.

Contains the per-nuclide components of a depletion chain.
"""

from collections import namedtuple
try:
    import lxml.etree as ET
except ImportError:
    import xml.etree.ElementTree as ET

DecayTuple = namedtuple('DecayTuple', 'type target branching_ratio')
ReactionTuple = namedtuple('ReactionTuple', 'type target Q branching_ratio')


class Nuclide(object):
    """The Nuclide class.

    Contains everything in a depletion chain relating to a single nuclide.

    Attributes
    ----------
    name : str
        Name of nuclide.
    half_life : float
        Half life of nuclide in s^-1.
    decay_energy : float
        Energy deposited from decay in eV.
    n_decay_modes : int
        Number of decay pathways.
    decay_modes : list of DecayTuple
        Decay mode information. Each element of the list is a named tuple with
        attributes 'type', 'target', and 'branching_ratio'.
    n_reaction_paths : int
        Number of possible reaction pathways.
    reactions : list of ReactionTuple
        Reaction information. Each element of the list is a named tuple with
        attribute 'type', 'target', 'Q', and 'branching_ratio'.
    yield_data : dict of float to list
        Maps tabulated energy to list of (product, yield) for all
        neutron-induced fission products.
    yield_energies : list of float
        Energies at which fission product yiels exist

    """

    def __init__(self):
        # Information about the nuclide
        self.name = None
        self.half_life = None
        self.decay_energy = 0.0

        # Decay paths
        self.decay_modes = []

        # Reaction paths
        self.reactions = []

        # Neutron fission yields, if present
        self.yield_data = {}
        self.yield_energies = []

    @property
    def n_decay_modes(self):
        """Number of decay modes."""
        return len(self.decay_modes)

    @property
    def n_reaction_paths(self):
        """Number of possible reaction pathways."""
        return len(self.reactions)

    @classmethod
    def xml_read(cls, element):
        """Read nuclide from an XML element.

        Parameters
        ----------
        element : xml.etree.ElementTree.Element
            XML element to write nuclide data to

        Returns
        -------
        nuc : Nuclide
            Instance of a nuclide

        """
        nuc = cls()
        nuc.name = element.get('name')

        # Check for half-life
        if 'half_life' in element.attrib:
            nuc.half_life = float(element.get('half_life'))
            nuc.decay_energy = float(element.get('decay_energy', '0'))

        # Check for decay paths
        for decay_elem in element.iter('decay_type'):
            d_type = decay_elem.get('type')
            target = decay_elem.get('target')
            branching_ratio = float(decay_elem.get('branching_ratio'))
            nuc.decay_modes.append(DecayTuple(d_type, target, branching_ratio))

        # Check for reaction paths
        for reaction_elem in element.iter('reaction_type'):
            r_type = reaction_elem.get('type')
            Q = float(reaction_elem.get('Q', '0'))
            branching_ratio = float(reaction_elem.get('branching_ratio', '1'))

            # If the type is not fission, get target and Q value, otherwise
            # just set null values
            if r_type != 'fission':
                target = reaction_elem.get('target')
            else:
                target = None

            # Append reaction
            nuc.reactions.append(ReactionTuple(
                r_type, target, Q, branching_ratio))

        fpy_elem = element.find('neutron_fission_yields')
        if fpy_elem is not None:
            for yields_elem in fpy_elem.iter('fission_yields'):
                E = float(yields_elem.get('energy'))
                products = yields_elem.find('products').text.split()
                yields = [float(y) for y in
                          yields_elem.find('data').text.split()]
                nuc.yield_data[E] = list(zip(products, yields))
            nuc.yield_energies = list(sorted(nuc.yield_data.keys()))

        return nuc

    def xml_write(self):
        """Write nuclide to XML element.

        Returns
        -------
        elem : xml.etree.ElementTree.Element
            XML element to write nuclide data to

        """
        elem = ET.Element('nuclide_table')
        elem.set('name', self.name)

        if self.half_life is not None:
            elem.set('half_life', str(self.half_life))
            elem.set('decay_modes', str(len(self.decay_modes)))
            elem.set('decay_energy', str(self.decay_energy))
            for mode, daughter, br in self.decay_modes:
                mode_elem = ET.SubElement(elem, 'decay_type')
                mode_elem.set('type', mode)
                mode_elem.set('target', daughter)
                mode_elem.set('branching_ratio', str(br))

        elem.set('reactions', str(len(self.reactions)))
        for rx, daughter, Q, br in self.reactions:
            rx_elem = ET.SubElement(elem, 'reaction_type')
            rx_elem.set('type', rx)
            rx_elem.set('Q', str(Q))
            if rx != 'fission':
                rx_elem.set('target', daughter)
            if br != 1.0:
                rx_elem.set('branching_ratio', str(br))

        if self.yield_data:
            fpy_elem = ET.SubElement(elem, 'neutron_fission_yields')
            energy_elem = ET.SubElement(fpy_elem, 'energies')
            energy_elem.text = ' '.join(str(E) for E in self.yield_energies)

            for E in self.yield_energies:
                yields_elem = ET.SubElement(fpy_elem, 'fission_yields')
                yields_elem.set('energy', str(E))

                products_elem = ET.SubElement(yields_elem, 'products')
                products_elem.text = ' '.join(x[0] for x in self.yield_data[E])
                data_elem = ET.SubElement(yields_elem, 'data')
                data_elem.text = ' '.join(str(x[1]) for x in self.yield_data[E])

        return elem
