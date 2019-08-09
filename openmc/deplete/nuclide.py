"""Nuclide module.

Contains the per-nuclide components of a depletion chain.
"""

from numbers import Real
from collections import namedtuple
from collections.abc import Iterable, Mapping

try:
    import lxml.etree as ET
except ImportError:
    import xml.etree.ElementTree as ET

from numpy import asarray, empty

from openmc.checkvalue import check_type


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
        Energies at which fission product yields exist

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
            nuc.yield_data = FissionYieldDistribution.from_xml_element(fpy_elem)
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
            self.yield_data.to_xml_element(fpy_elem)

        return elem


class FissionYieldDistribution(Mapping):
    """Class for storing energy-dependent fission yields for a single nuclide

    Can be used as a dictionary mapping energies and products to fission
    yields::

        >>> fydist = FissionYieldDistribution.from_dict({
        ...     {0.0253: {"Xe135": 0.021}})
        >>> fydist[0.0253]["Xe135"]
        0.021

    Parameters
    ----------
    ordered_energies : iterable of float
        Energies for which fission yield data exist
    orderded_products : iterable of str
        Fission products produced by this parent at all energies
    group_fission_yields : numpy.ndarray or iterable of iterable of float
        Array of shape ``(n_energy, n_products)`` where
        ``group_fission_yields[g][j]`` is the yield of
        ``ordered_products[j]`` due to a fission in energy region ``g``.

    Attributes
    ----------
    energies : tuple
        Energies for which fission yields exist. Converted for
        indexing
    products : tuple
        Fission products produced at all energies. Converted
        for indexing
    yield_matrix : numpy.ndarray
        Array ``(n_energy, n_products)`` where
        ``yield_matrix[g, j]`` is the fission yield of product
        ``j`` for energy group ``g``.

    See Also
    --------
    :meth:`from_xml_element`, :meth:`from_dict`
    """

    def __init__(self, ordered_energies, ordered_products, group_fission_yields):
        check_type("energies", ordered_energies, Iterable, Real)
        self.energies = tuple(ordered_energies)
        check_type("products", ordered_products, Iterable, str)
        self.products = tuple(ordered_products)
        yield_matrix = asarray(group_fission_yields, dtype=float)
        if yield_matrix.shape != (len(self.energies), len(self.products)):
            raise ValueError(
                "Shape of yield matrix inconsistent. "
                "Should be ({}, {}), is {}".format(
                    len(ordered_energies), len(ordered_products),
                    yield_matrix.shape))
        self.yield_matrix = yield_matrix

    def __len__(self):
        return len(self.energies)

    def __getitem__(self, energy):
        if energy not in self.energies:
            raise KeyError(energy)
        return _FissionYield(
            self.products, self.yield_matrix[self.energies.index(energy)])

    def __iter__(self):
        return iter(self.energies)

    @classmethod
    def from_xml_element(cls, element):
        """Construct a distribution from a depletion chain xml file

        Parameters
        ----------
        element : xml.etree.ElementTree.Element
            XML element to pull fission yield data from

        Returns
        -------
        FissionYieldDistribution
        """
        yields = {}
        for elem_index, yield_elem in enumerate(element.iter("fission_yields")):
            energy = float(yield_elem.get("energy"))
            products = yield_elem.find("products").text.split()
            yield_mapobj = map(float, yield_elem.find("data").text.split())
            # Get a map of products to their corresponding yield
            yields[energy] = dict(zip(products, yield_mapobj))

        return cls.from_dict(yields)

    @classmethod
    def from_dict(cls, fission_yields):
        """Construct a distribution from a dictionary of yields

        Parameters
        -----------
        fission_yields : dict
            Dictionary ``{energy: {product: yield}}``

        Returns
        -------
        FissionYieldDistribution
        """
        # mapping {energy: {product: value}}
        energies = tuple(sorted(fission_yields))

        # Get a consistent set of products to produce a matrix of yields
        shared_prod = set()
        for prod_set in map(set, fission_yields.values()):
            shared_prod |= prod_set
        ordered_prod = tuple(sorted(shared_prod))

        yield_matrix = empty((len(energies), len(shared_prod)))

        for g_index, energy in enumerate(energies):
            prod_map = fission_yields[energy]
            for prod_ix, product in enumerate(ordered_prod):
                yield_val = prod_map.get(product)
                if yield_val is None:
                    yield_matrix[g_index, prod_ix] = 0.0
                else:
                    yield_matrix[g_index, prod_ix] = yield_val

        return cls(energies, ordered_prod, yield_matrix)

    def to_xml_element(self, root):
        """Write fission yield data to an xml element

        Parameters
        ----------
        root : xml.etree.ElementTree.Element
            Element to write distribution data to
        """
        for energy, yield_obj in self.items():
            yield_element = ET.SubElement(root, "fission_yields")
            yield_element.set("energy", str(energy))
            product_elem = ET.SubElement(yield_element, "products")
            product_elem.text = " ".join(map(str, yield_obj.products))
            data_elem = ET.SubElement(yield_element, "data")
            data_elem.text = " ".join(map(str, yield_obj.yields))


class _FissionYield(Mapping):
    """Mapping for fission yields of a parent at a specific energy

    Separated to support nested dictionary-like behavior for
    :class:`FissionYieldDistribution`, and allowing math operations
    on a single vector of yields

    Parameters
    ----------
    products : tuple of str
        Products for this specific distribution
    yields : numpy.ndarray
        View into associated :attr:`FissionYieldDistribution.yield_matrix`
    """

    def __init__(self, products, yields):
        self.products = products
        self.yields = yields

    def __getitem__(self, product):
        if product not in self.products:
            raise KeyError(product)
        return self.yields[self.products.index(product)]

    def __len__(self):
        return len(self.products)

    def __iter__(self):
        return iter(self.products)

    def __iadd__(self, other):
        """Increment value from other fission yield"""
        self.yields += other.yields
        return self

    def __mul__(self, value):
        return _FissionYield(self.products, self.yields * value)

    def copy(self):
        """Return an identical yield object, with unique yields"""
        return _FissionYield(self.products, self.yields.copy())
