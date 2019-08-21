"""Nuclide module.

Contains the per-nuclide components of a depletion chain.
"""

from collections import namedtuple
from collections.abc import Mapping

try:
    import lxml.etree as ET
except ImportError:
    import xml.etree.ElementTree as ET

from numpy import empty

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

    Parameters
    ----------
    name : str, optional
        GND name of this nuclide, e.g. ``"He4"``, ``"Am242_m1"``

    Attributes
    ----------
    name : str or None
        Name of nuclide.
    half_life : float or None
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
    yield_data : FissionYieldDistribution or None
        Fission product yields at tabulated energies for this nuclide. Can be
        treated as a nested dictionary ``{energy: {product: yield}}``
    yield_energies : tuple of float or None
        Energies at which fission product yields exist
    """

    def __init__(self, name=None):
        # Information about the nuclide
        self.name = name
        self.half_life = None
        self.decay_energy = 0.0

        # Decay paths
        self.decay_modes = []

        # Reaction paths
        self.reactions = []

        # Neutron fission yields, if present
        self._yield_data = None

    @property
    def n_decay_modes(self):
        return len(self.decay_modes)

    @property
    def n_reaction_paths(self):
        return len(self.reactions)

    @property
    def yield_data(self):
        if self._yield_data is None:
            return {}
        return self._yield_data

    @yield_data.setter
    def yield_data(self, fission_yields):
        check_type("fission_yields", fission_yields, Mapping)
        if isinstance(fission_yields, FissionYieldDistribution):
            self._yield_data = fission_yields
        self._yield_data = FissionYieldDistribution(fission_yields)

    @property
    def yield_energies(self):
        if self._yield_data is None:
            return None
        return self.yield_data.energies

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
    """Energy-dependent fission product yields for a single nuclide

    Can be used as a dictionary mapping energies and products to fission
    yields::

        >>> fydist = FissionYieldDistribution{
        ...     {0.0253: {"Xe135": 0.021}})
        >>> fydist[0.0253]["Xe135"]
        0.021

    Parameters
    ----------
    energies : iterable of float
        Energies for which fission yield data exist. Must be ordered
        by increasing energy
    products : iterable of str
        Fission products produced by this parent at all energies.
        Must be ordered alphabetically
    yield_matrix : numpy.ndarray or iterable of iterable of float
        Array of shape ``(n_energy, n_products)`` where
        ``yield_matrix[g][j]`` is the yield of
        ``ordered_products[j]`` due to a fission in energy region ``g``.

    Attributes
    ----------
    energies : tuple
        Energies for which fission yields exist. Sorted by
        increasing energy
    products : tuple
        Fission products produced at all energies. Sorted by name.
    yield_matrix : numpy.ndarray
        Array ``(n_energy, n_products)`` where
        ``yield_matrix[g, j]`` is the fission yield of product
        ``j`` for energy group ``g``.

    See Also
    --------
    * :meth:`from_xml_element` - Construction methods
    * :class:`FissionYield` - Class used for storing yields at a given energy
    """

    def __init__(self, fission_yields):
        # mapping {energy: {product: value}}
        energies = sorted(fission_yields)

        # Get a consistent set of products to produce a matrix of yields
        shared_prod = set.union(*(set(x) for x in fission_yields.values()))
        ordered_prod = sorted(shared_prod)

        yield_matrix = empty((len(energies), len(shared_prod)))

        for g_index, energy in enumerate(energies):
            prod_map = fission_yields[energy]
            for prod_ix, product in enumerate(ordered_prod):
                yield_val = prod_map.get(product)
                yield_matrix[g_index, prod_ix] = (
                    0.0 if yield_val is None else yield_val)
        self.energies = tuple(energies)
        self.products = tuple(ordered_prod)
        self.yield_matrix = yield_matrix

    def __len__(self):
        return len(self.energies)

    def __getitem__(self, energy):
        if energy not in self.energies:
            raise KeyError(energy)
        return FissionYield(
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
        all_yields = {}
        for elem_index, yield_elem in enumerate(element.iter("fission_yields")):
            energy = float(yield_elem.get("energy"))
            products = yield_elem.find("products").text.split()
            yields = map(float, yield_elem.find("data").text.split())
            # Get a map of products to their corresponding yield
            all_yields[energy] = dict(zip(products, yields))

        return cls(all_yields)

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


class FissionYield(Mapping):
    """Mapping for fission yields of a parent at a specific energy

    Separated to support nested dictionary-like behavior for
    :class:`FissionYieldDistribution`, and allowing math operations
    on a single vector of yields. Can in turn be used like a
    dictionary to fetch fission yields.
    Supports multiplication of a scalar to scale the fission
    yields and addition of another set of yields.

    Does not support resizing / inserting new products that do
    not exist.

    Parameters
    ----------
    products : tuple of str
        Products for this specific distribution
    yields : numpy.ndarray
        Fission product yields for each product in ``products``

    Attributes
    ----------
    products : tuple of str
        Products for this specific distribution
    yields : numpy.ndarray
        Fission product yields for each product in ``products``

    Examples
    --------
    >>> import numpy
    >>> fy_vector = FissionYield(
    ...     ("Xe135", "I129", "Sm149"),
    ...     numpy.array((0.002, 0.001, 0.0003)))
    >>> fy_vector["Xe135"]
    0.002
    >>> new = fy_vector.copy()
    >>> fy_vector *= 2
    >>> fy_vector["Xe135"]
    0.004
    >>> new["Xe135"]
    0.002
    >>> (new + fy_vector)["Sm149"]
    0.0009
    >>> dict(new)
    {"Xe135": 0.002, "I129": 0.001, "Sm149": 0.0003}
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

    def items(self):
        return zip(self.products, self.yields)

    def __add__(self, other):
        new = self.copy()
        new += other
        return new

    def __iadd__(self, other):
        """Increment value from other fission yield"""
        self.yields += other.yields
        return self

    def __imul__(self, scalar):
        self.yields *= scalar
        return self

    def __mul__(self, scalar):
        new = self.copy()
        new *= scalar
        return new

    def __rmul__(self, scalar):
        new = self.copy()
        new *= scalar
        return new

    def __repr__(self):
        return "<{}: {}>".format(self.__class__.__name__, repr(dict(self)))

    def copy(self):
        """Return an identical yield object, with unique yields"""
        return FissionYield(self.products, self.yields.copy())
