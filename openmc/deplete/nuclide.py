"""Nuclide module.

Contains the per-nuclide components of a depletion chain.
"""

from collections import namedtuple, defaultdict
from warnings import warn
try:
    import lxml.etree as ET
except ImportError:
    import xml.etree.ElementTree as ET


DecayTuple = namedtuple('DecayTuple', 'type target branching_ratio')
DecayTuple.__doc__ = """\
Decay mode information

Parameters
----------
type : str
    Type of the decay mode, e.g., 'beta-'
target : str
    Nuclide resulting from decay
branching_ratio : float
    Branching ratio of the decay mode

"""
try:
    DecayTuple.type.__doc__ = None
    DecayTuple.target.__doc__ = None
    DecayTuple.branching_ratio.__doc__ = None
except AttributeError:
    # Can't set __doc__ on properties on Python 3.4
    pass


ReactionTuple = namedtuple('ReactionTuple', 'type target Q branching_ratio')
ReactionTuple.__doc__ = """\
Transmutation reaction information

Parameters
----------
type : str
    Type of the reaction, e.g., 'fission'
target : str
    nuclide resulting from reaction
Q : float
    Q value of the reaction in [eV]
branching_ratio : float
    Branching ratio of the reaction

"""
try:
    ReactionTuple.type.__doc__ = None
    ReactionTuple.target.__doc__ = None
    ReactionTuple.Q.__doc__ = None
    ReactionTuple.branching_ratio.__doc__ = None
except AttributeError:
    pass


class Nuclide(object):
    """Decay modes, reactions, and fission yields for a single nuclide.

    Attributes
    ----------
    name : str
        Name of nuclide.
    half_life : float
        Half life of nuclide in [s].
    decay_energy : float
        Energy deposited from decay in [eV].
    n_decay_modes : int
        Number of decay pathways.
    decay_modes : list of openmc.deplete.DecayTuple
        Decay mode information. Each element of the list is a named tuple with
        attributes 'type', 'target', and 'branching_ratio'.
    n_reaction_paths : int
        Number of possible reaction pathways.
    reactions : list of openmc.deplete.ReactionTuple
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
        return len(self.decay_modes)

    @property
    def n_reaction_paths(self):
        return len(self.reactions)

    @classmethod
    def from_xml(cls, element, fission_q=None):
        """Read nuclide from an XML element.

        Parameters
        ----------
        element : xml.etree.ElementTree.Element
            XML element to write nuclide data to
        fission_q : None or float
            User-supplied fission Q value [eV].
            Will be read from the element if not given

        Returns
        -------
        nuc : openmc.deplete.Nuclide
            Instance of a nuclide

        """
        nuc = cls()
        nuc.name = element.get('name')

        # Check for half-life
        if 'half_life' in element.attrib:
            nuc.half_life = float(element.get('half_life'))
            nuc.decay_energy = float(element.get('decay_energy', '0'))

        # Check for decay paths
        for decay_elem in element.iter('decay'):
            d_type = decay_elem.get('type')
            target = decay_elem.get('target')
            branching_ratio = float(decay_elem.get('branching_ratio'))
            nuc.decay_modes.append(DecayTuple(d_type, target, branching_ratio))

        # Check for reaction paths
        for reaction_elem in element.iter('reaction'):
            r_type = reaction_elem.get('type')
            Q = float(reaction_elem.get('Q', '0'))
            branching_ratio = float(reaction_elem.get('branching_ratio', '1'))

            # If the type is not fission, get target and Q value, otherwise
            # just set null values
            if r_type != 'fission':
                target = reaction_elem.get('target')
            else:
                target = None
                if fission_q is not None:
                    Q = fission_q

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

    def to_xml_element(self):
        """Write nuclide to XML element.

        Returns
        -------
        elem : xml.etree.ElementTree.Element
            XML element to write nuclide data to

        """
        elem = ET.Element('nuclide')
        elem.set('name', self.name)

        if self.half_life is not None:
            elem.set('half_life', str(self.half_life))
            elem.set('decay_modes', str(len(self.decay_modes)))
            elem.set('decay_energy', str(self.decay_energy))
            for mode, daughter, br in self.decay_modes:
                mode_elem = ET.SubElement(elem, 'decay')
                mode_elem.set('type', mode)
                mode_elem.set('target', daughter)
                mode_elem.set('branching_ratio', str(br))

        elem.set('reactions', str(len(self.reactions)))
        for rx, daughter, Q, br in self.reactions:
            rx_elem = ET.SubElement(elem, 'reaction')
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

    def validate(self, strict=True, quiet=False, tolerance=1e-4):
        """Search for possible inconsistencies

        The following checks are performed:

            1) for all non-fission reactions and decay modes,
               does the sum of branching ratios equal about one?
            2) for fission reactions, does the sum of fission yield
               fractions equal about two?

        Parameters
        ----------
        strict : bool, optional
            Raise exceptions at the first inconsistency if true.
            Otherwise mark a warning
        quiet : bool, optional
            Flag to suppress warnings and return immediately at
            the first inconsistency. Used only if
            ``strict`` does not evaluate to ``True``.
        tolerance : float, optional
            Absolute tolerance for comparisons. Used to compare computed
            value ``x`` to intended value ``y`` as::

                valid = (y - tolerance <= x <= y + tolerance)

        Returns
        -------
        valid : bool
            True if no inconsistencies were found

        Raises
        ------
        ValueError
            If ``strict`` evaluates to ``True`` and an inconistency was
            found

        See Also
        --------
        openmc.deplete.Chain.validate
        """

        msg_func = ("Nuclide {name} has {prop} that sum to {actual} "
                    "instead of {expected} +/- {tol:7.4e}").format
        valid = True

        # check decay modes
        if self.decay_modes:
            sum_br = sum(m.branching_ratio for m in self.decay_modes)
            stat = 1.0 - tolerance <= sum_br <= 1.0 + tolerance
            if not stat:
                msg = msg_func(
                    name=self.name, actual=sum_br, expected=1.0, tol=tolerance,
                    prop="decay mode branch ratios")
                if strict:
                    raise ValueError(msg)
                elif quiet:
                    return False
                warn(msg)
                valid = False

        if self.reactions:
            type_map = defaultdict(set)
            for reaction in self.reactions:
                type_map[reaction.type].add(reaction)
            for rxn_type, reactions in type_map.items():
                sum_rxn = sum(rx.branching_ratio for rx in reactions)
                stat = 1.0 - tolerance <= sum_rxn <= 1.0 + tolerance
                if stat:
                    continue
                msg = msg_func(
                    name=self.name, actual=sum_br, expected=1.0, tol=tolerance,
                    prop="{} reaction branch ratios".format(rxn_type))
                if strict:
                    raise ValueError(msg)
                elif quiet:
                    return False
                warn(msg)
                valid = False

        if self.yield_data:
            for energy, yield_list in self.yield_data.items():
                sum_yield = sum(y[1] for y in yield_list)
                stat = 2.0 - tolerance <= sum_yield <= 2.0 + tolerance
                if stat:
                    continue
                msg = msg_func(
                    name=self.name, actual=sum_yield,
                    expected=2.0, tol=tolerance,
                    prop="fission yields (E = {:7.4e} eV)".format(energy))
                if strict:
                    raise ValueError(msg)
                elif quiet:
                    return False
                warn(msg)
                valid = False

        return valid
