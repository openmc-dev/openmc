"""Nuclide module.xml.etree.Ele

Contains the per-nuclide components of a depletion chain.
"""

import bisect
from collections.abc import Mapping
from collections import namedtuple, defaultdict
from warnings import warn
from numbers import Real

import lxml.etree as ET
import numpy as np

from openmc.checkvalue import check_type
from openmc.stats import Univariate

__all__ = [
    "DecayTuple", "ReactionTuple", "Nuclide", "FissionYield",
    "FissionYieldDistribution"]


DecayTuple = namedtuple('DecayTuple', 'type target branching_ratio')
DecayTuple.__doc__ = """\
Decay mode information

Parameters
----------
type : str
    Type of the decay mode, e.g., 'beta-'
target : str or None
    Nuclide resulting from decay. A value of ``None`` implies the
    target does not exist in the currently configured depletion
    chain
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
target : str or None
    Nuclide resulting from reaction. A value of ``None``
    implies either no single target, e.g. from fission,
    or that the target nuclide is not considered
    in the current depletion chain
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


class Nuclide:
    """Decay modes, reactions, and fission yields for a single nuclide.

    Parameters
    ----------
    name : str, optional
        GNDS name of this nuclide, e.g. ``"He4"``, ``"Am242_m1"``

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
    sources : dict
        Dictionary mapping particle type as string to energy distribution of
        decay source represented as :class:`openmc.stats.Univariate`
    yield_data : FissionYieldDistribution or None
        Fission product yields at tabulated energies for this nuclide. Can be
        treated as a nested dictionary ``{energy: {product: yield}}``
    spont_yield_data : FissionYieldDistribution or None
        Spontaneous fission product yields. To keep the same structure as for 
        yield_data this is treated as a nested dictionary ``{energy: {product: yield}}``
        with one energy value = 0.0.
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

        # Decay sources
        self.sources = {}

        # Neutron fission yields, if present
        self._yield_data = None

        # Spontanteous fission yields, if present
        self._spont_yield_data = None

    def __repr__(self):
        n_modes, n_rx = self.n_decay_modes, self.n_reaction_paths
        return f"<Nuclide: {self.name} ({n_modes} modes, {n_rx} reactions)>"

    @property
    def n_decay_modes(self):
        return len(self.decay_modes)

    @property
    def n_reaction_paths(self):
        return len(self.reactions)

    @property
    def yield_data(self):
        if self._yield_data is None:
            return None
        return self._yield_data

    @property
    def spont_yield_data(self):
        if self._spont_yield_data is None:
            return None
        return self._spont_yield_data

    @yield_data.setter
    def yield_data(self, fission_yields):
        if fission_yields is None:
            self._yield_data = None
        else:
            check_type("fission_yields", fission_yields, Mapping)
            if isinstance(fission_yields, FissionYieldDistribution):
                self._yield_data = fission_yields
            else:
                self._yield_data = FissionYieldDistribution(fission_yields)

    @spont_yield_data.setter
    def spont_yield_data(self, fission_yields):
        if fission_yields is None:
            self._spont_yield_data = None
        else:
            check_type("fission_yields", fission_yields, Mapping)
            if isinstance(fission_yields, FissionYieldDistribution):
                self._spont_yield_data = fission_yields
            else:
                self._spont_yield_data = FissionYieldDistribution(fission_yields)

    @property
    def yield_energies(self):
        if self._yield_data is None:
            return None
        return self.yield_data.energies

    def add_decay_mode(self, type, target, branching_ratio):
        """Add decay mode to the nuclide

        Parameters
        ----------
        type : str
            Type of the decay mode, e.g., 'beta-'
        target : str or None
            Nuclide resulting from decay. A value of ``None`` implies the
            target does not exist in the currently configured depletion
            chain
        branching_ratio : float
            Branching ratio of the decay mode

        """
        self.decay_modes.append(
            DecayTuple(type, target, branching_ratio)
        )

    def add_reaction(self, type, target, Q, branching_ratio):
        """Add transmutation reaction to the nuclide

        Parameters
        ----------
        type : str
            Type of the reaction, e.g., 'fission'
        target : str or None
            Nuclide resulting from reaction. A value of ``None``
            implies either no single target, e.g. from fission,
            or that the target nuclide is not considered
            in the current depletion chain
        Q : float
            Q value of the reaction in [eV]
        branching_ratio : float
            Branching ratio of the reaction

        """
        self.reactions.append(
            ReactionTuple(type, target, Q, branching_ratio)
        )

    @classmethod
    def from_xml(cls, element, root=None, fission_q=None):
        """Read nuclide from an XML element.

        Parameters
        ----------
        element : lxml.etree._Element
            XML element to read nuclide data from
        root : lxml.etree._Element, optional
            Root XML element for chain file (only used when fission product
            yields are borrowed from another parent)
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
            if target is not None and target.lower() == "nothing":
                target = None
            branching_ratio = float(decay_elem.get('branching_ratio'))
            nuc.decay_modes.append(DecayTuple(d_type, target, branching_ratio))

        # Check for sources
        for src_elem in element.iter('source'):
            particle = src_elem.get('particle')
            distribution = Univariate.from_xml_element(src_elem)
            nuc.sources[particle] = distribution

        # Check for reaction paths
        for reaction_elem in element.iter('reaction'):
            r_type = reaction_elem.get('type')
            Q = float(reaction_elem.get('Q', '0'))
            branching_ratio = float(reaction_elem.get('branching_ratio', '1'))

            # If the type is not fission, get target and Q value, otherwise
            # just set null values
            if r_type != 'fission':
                target = reaction_elem.get('target')
                if target is not None and target.lower() == "nothing":
                    target = None
            else:
                target = None
                if fission_q is not None:
                    Q = fission_q

            # Append reaction
            nuc.reactions.append(ReactionTuple(
                r_type, target, Q, branching_ratio))

        fpy_elem = element.find('neutron_fission_yields')
        if fpy_elem is not None:
            # Check for use of FPY from other nuclide
            parent = fpy_elem.get('parent')
            if parent is not None:
                assert root is not None
                fpy_elem = root.find(
                    f'.//nuclide[@name="{parent}"]/neutron_fission_yields'
                )
                if fpy_elem is None:
                    raise ValueError(
                        "Fission product yields for {0} borrow from {1}, but {1} is"
                        " not present in the chain file or has no yields.".format(
                            nuc.name, parent
                        ))
                nuc._fpy = parent

            nuc.yield_data = FissionYieldDistribution.from_xml_element(fpy_elem)


        sfy_elem = element.find('spont_fission_yields')
        if sfy_elem is not None:
            # Check for use of FPY from other nuclide
            parent = sfy_elem.get('parent')
            if parent is not None:
                assert root is not None
                sfy_elem = root.find(
                    f'.//nuclide[@name="{parent}"]/spont_fission_yields'
                )
                if sfy_elem is None:
                    raise ValueError(
                        "Spontanteous fission product yields for {0} borrow from {1}, but {1} is"
                        " not present in the chain file or has no yields.".format(
                            nuc.name, parent
                        ))
                nuc._sfy = parent

            nuc.spont_yield_data = FissionYieldDistribution.from_xml_element(sfy_elem)


        return nuc

    def to_xml_element(self):
        """Write nuclide to XML element.

        Returns
        -------
        elem : lxml.etree._Element
            XML element to write nuclide data to

        """
        elem = ET.Element('nuclide')
        elem.set('name', self.name)

        if self.half_life is not None:
            elem.set('half_life', str(self.half_life))
            elem.set('decay_modes', str(len(self.decay_modes)))
            elem.set('decay_energy', str(self.decay_energy))
            for mode_type, daughter, br in self.decay_modes:
                mode_elem = ET.SubElement(elem, 'decay')
                mode_elem.set('type', mode_type)
                if daughter:
                    mode_elem.set('target', daughter)
                mode_elem.set('branching_ratio', str(br))

        # Write decay sources
        if self.sources:
            for particle, source in self.sources.items():
                src_elem = source.to_xml_element('source')
                src_elem.set('particle', particle)
                elem.append(src_elem)

        elem.set('reactions', str(len(self.reactions)))
        for rx, daughter, Q, br in self.reactions:
            rx_elem = ET.SubElement(elem, 'reaction')
            rx_elem.set('type', rx)
            rx_elem.set('Q', str(Q))
            if daughter is not None:
                rx_elem.set('target', daughter)
            if br != 1.0:
                rx_elem.set('branching_ratio', str(br))

        if self.yield_data:
            fpy_elem = ET.SubElement(elem, 'neutron_fission_yields')

            if hasattr(self, '_fpy'):
                # Check for link to other nuclide data
                fpy_elem.set('parent', self._fpy)
            else:
                energy_elem = ET.SubElement(fpy_elem, 'energies')
                energy_elem.text = ' '.join(str(E) for E in self.yield_energies)
                self.yield_data.to_xml_element(fpy_elem)

        if self.spont_yield_data:
            sfy_elem = ET.SubElement(elem, 'spont_fission_yields')

            if hasattr(self, '_sfy'):
                # Check for link to other nuclide data
                sfy_elem.set('parent', self._sfy)
            else:
                self.spont_yield_data.to_xml_element(sfy_elem)

        return elem

    def validate(self, strict=True, quiet=False, tolerance=1e-4):
        """Search for possible inconsistencies

        The following checks are performed:

            1) for all non-fission reactions and decay modes,
               does the sum of branching ratios equal about one?
            2) for fission reactions, does the sum of fission yield
               fractions equal about two?
            3) For spontanous fission, does the sum of fission yield
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
                    prop=f"{rxn_type} reaction branch ratios")
                if strict:
                    raise ValueError(msg)
                elif quiet:
                    return False
                warn(msg)
                valid = False

        if self.yield_data:
            for energy, fission_yield in self.yield_data.items():
                sum_yield = fission_yield.yields.sum()
                stat = 2.0 - tolerance <= sum_yield <= 2.0 + tolerance
                if stat:
                    continue
                msg = msg_func(
                    name=self.name, actual=sum_yield,
                    expected=2.0, tol=tolerance,
                    prop=f"fission yields (E = {energy:7.4e} eV)")
                if strict:
                    raise ValueError(msg)
                elif quiet:
                    return False
                warn(msg)
                valid = False

        if self.spont_yield_data:
            for _, fission_yield in self.spont_yield_data.items():
                sum_yield = fission_yield.yields.sum()
                stat = 2.0 - tolerance <= sum_yield <= 2.0 + tolerance
                if stat:
                    continue
                msg = msg_func(
                    name=self.name, actual=sum_yield,
                    expected=2.0, tol=tolerance,
                    prop=f"spontaneous fission yields")
                if strict:
                    raise ValueError(msg)
                elif quiet:
                    return False
                warn(msg)
                valid = False


        return valid


class FissionYieldDistribution(Mapping):
    """Energy-dependent fission product yields for a single nuclide

    Can be used as a dictionary mapping energies and products to fission
    yields::

        >>> fydist = FissionYieldDistribution(
        ...     {0.0253: {"Xe135": 0.021}})
        >>> fydist[0.0253]["Xe135"]
        0.021

    Parameters
    ----------
    fission_yields : dict
        Dictionary of energies and fission product yields for that energy.
        Expected to be of the form ``{float: {str: float}}``. The first
        float is the energy, typically in eV, that represents this
        distribution. The underlying dictionary maps fission products
        to their respective yields.

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

        yield_matrix = np.empty((len(energies), len(shared_prod)))

        for g_index, energy in enumerate(energies):
            prod_map = fission_yields[energy]
            for prod_ix, product in enumerate(ordered_prod):
                yield_matrix[g_index, prod_ix] = prod_map.get(product, 0.0)
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

    def __repr__(self):
        return "<{} with {} products at {} energies>".format(
            self.__class__.__name__, self.yield_matrix.shape[1],
            len(self.energies))

    @classmethod
    def from_xml_element(cls, element):
        """Construct a distribution from a depletion chain xml file

        Parameters
        ----------
        element : lxml.etree._Element
            XML element to pull fission yield data from

        Returns
        -------
        FissionYieldDistribution
        """
        all_yields = {}
        for yield_elem in element.iter("fission_yields"):
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
        root : lxml.etree._Element
            Element to write distribution data to
        """
        for energy, yield_obj in self.items():
            yield_element = ET.SubElement(root, "fission_yields")
            yield_element.set("energy", str(energy))
            product_elem = ET.SubElement(yield_element, "products")
            product_elem.text = " ".join(map(str, yield_obj.products))
            data_elem = ET.SubElement(yield_element, "data")
            data_elem.text = " ".join(map(str, yield_obj.yields))

    def restrict_products(self, possible_products):
        """Return a new distribution with select products

        .. versionadded:: 0.12

        Parameters
        ----------
        possible_products : iterable of str
            Candidate pool of fission products. Existing products
            not contained here will not exist in the new instance

        Returns
        -------
        FissionYieldDistribution or None
            A value of None indicates no values in
            ``possible_products`` exist in :attr:`products`

        """

        overlap = set(self.products).intersection(possible_products)
        if not overlap:
            return None

        products = sorted(overlap)
        indices = np.searchsorted(self.products, products)

        # coerce back to dictionary to pass back to __init__
        new_yields = {}
        for ene, yields in zip(self.energies, self.yield_matrix.copy()):
            new_yields[ene] = dict(zip(products, yields[indices]))

        return type(self)(new_yields)


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
        Sorted products for this specific distribution
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
    ...     ("I129", "Sm149", "Xe135"),
    ...     numpy.array((0.001, 0.0003, 0.002)))
    >>> fy_vector["Xe135"]
    0.002
    >>> new = FissionYield(fy_vector.products, fy_vector.yields.copy())
    >>> fy_vector *= 2
    >>> fy_vector["Xe135"]
    0.004
    >>> new["Xe135"]
    0.002
    >>> (new + fy_vector)["Sm149"]
    0.0009
    >>> dict(new) == {"Xe135": 0.002, "I129": 0.001, "Sm149": 0.0003}
    True
    """

    def __init__(self, products, yields):
        self.products = products
        self.yields = yields

    def __contains__(self, product):
        ix = bisect.bisect_left(self.products, product)
        return ix != len(self.products) and self.products[ix] == product

    def __getitem__(self, product):
        ix = bisect.bisect_left(self.products, product)
        if ix == len(self.products) or self.products[ix] != product:
            raise KeyError(product)
        return self.yields[ix]

    def __len__(self):
        return len(self.products)

    def __iter__(self):
        return iter(self.products)

    def items(self):
        """Return pairs of product, yield"""
        return zip(self.products, self.yields)

    def __add__(self, other):
        """Add one set of fission yields to this set, return new yields"""
        if not isinstance(other, FissionYield):
            return NotImplemented
        new = FissionYield(self.products, self.yields.copy())
        new += other
        return new

    def __iadd__(self, other):
        """Increment value from other fission yield"""
        if not isinstance(other, FissionYield):
            return NotImplemented
        self.yields += other.yields
        return self

    def __radd__(self, other):
        return self + other

    def __imul__(self, scalar):
        """Scale these fission yields by a real scalar"""
        if not isinstance(scalar, Real):
            return NotImplemented
        self.yields *= scalar
        return self

    def __mul__(self, scalar):
        """Return a new set of yields scaled by a real scalar"""
        if not isinstance(scalar, Real):
            return NotImplemented
        new = FissionYield(self.products, self.yields.copy())
        new *= scalar
        return new

    def __rmul__(self, scalar):
        return self * scalar

    def __repr__(self):
        return f"<{self.__class__.__name__} containing {len(self)} products and yields>"

    def __deepcopy__(self, memo):
        result = FissionYield(self.products, self.yields.copy())
        memo[id(self)] = result
        return result

    # Avoid greedy numpy operations like np.float64 * fission_yield
    # converting this to an array on the fly. Force __rmul__ and
    # __radd__. See issue #1492
    __array_ufunc__ = None
