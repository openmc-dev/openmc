from __future__ import division

from collections import Iterable, OrderedDict
from numbers import Integral
import os
import sys
import copy
import abc

import numpy as np

import openmc
import openmc.checkvalue as cv
from openmc.mgxs import EnergyGroups


if sys.version_info[0] >= 3:
    basestring = str


# Supported cross section types
MGXS_TYPES = ['total',
              'transport',
              'absorption',
              'capture',
              'fission',
              'nu-fission',
              'kappa-fission',
              'scatter',
              'nu-scatter',
              'scatter matrix',
              'nu-scatter matrix',
              'chi']


# Supported domain types
# TODO: Implement Mesh domains
DOMAIN_TYPES = ['cell',
                'distribcell',
                'universe',
                'material']

# Supported domain classes
# TODO: Implement Mesh domains
_DOMAINS = [openmc.Cell,
           openmc.Universe,
           openmc.Material]


class MGXS(object):
    """An abstract multi-group cross section for some energy group structure
    within some spatial domain.

    This class can be used for both OpenMC input generation and tally data
    post-processing to compute spatially-homogenized and energy-integrated
    multi-group cross sections for deterministic neutronics calculations.

    NOTE: Users should instantiate the subclasses of this abstract class.

    Parameters
    ----------
    domain : Material or Cell or Universe
        The domain for spatial homogenization
    domain_type : {'material', 'cell', 'distribcell', 'universe'}
        The domain type for spatial homogenization
    energy_groups : EnergyGroups
        The energy group structure for energy condensation
    by_nuclide : bool
        If true, computes cross sections for each nuclide in domain
    nuclides : Iterable of basestring
        The user-specified nuclides to compute cross sections. If by_nuclide
        is True but nuclides are not specified by the user, all nuclides in the
        spatial domain will be used.
    name : str, optional
        Name of the multi-group cross section. Used as a label to identify
        tallies in OpenMC 'tallies.xml' file.

    Attributes
    ----------
    name : str, optional
        Name of the multi-group cross section
    rxn_type : str
        Reaction type (e.g., 'total', 'nu-fission', etc.)
    by_nuclide : bool
        If true, computes cross sections for each nuclide in domain
    domain : Material or Cell or Universe
        Domain for spatial homogenization
    domain_type : {'material', 'cell', 'distribcell', 'universe'}
        Domain type for spatial homogenization
    energy_groups : EnergyGroups
        Energy group structure for energy condensation
    tally_trigger : Trigger
        An (optional) tally precision trigger given to each tally used to
        compute the cross section
    tallies : OrderedDict
        OpenMC tallies needed to compute the multi-group cross section
    rxn_rate_tally : Tally
        Derived tally for the reaction rate tally used in the numerator to
        compute the multi-group cross section. This attribute is None
        unless the multi-group cross section has been computed.
    xs_tally : Tally
        Derived tally for the multi-group cross section. This attribute
        is None unless the multi-group cross section has been computed.
    num_subdomains : Integral
        The number of subdomains is unity for 'material', 'cell' and 'universe'
        domain types. When the  This is equal to the number of cell instances
        for 'distribcell' domain types (it is equal to unity prior to loading
        tally data from a statepoint file).
    num_nuclides : Integral
        The number of nuclides for which the multi-group cross section is
        being tracked. This is unity if the by_nuclide attribute is False.
    nuclides : list of str or 'sum'
        A list of nuclide string names (e.g., 'U-238', 'O-16') when by_nuclide
        is True and 'sum' when by_nuclide is False.
    sparse : bool
        Whether or not the MGXS' tallies use SciPy's LIL sparse matrix format
        for compressed data storage
    derived : bool
        Whether or not the MGXS is merged from one or more other MGXS

    """

    # This is an abstract class which cannot be instantiated
    __metaclass__ = abc.ABCMeta

    def __init__(self, domain=None, domain_type=None,
                 energy_groups=None, by_nuclide=False, name=''):

        self._name = ''
        self._rxn_type = None
        self._by_nuclide = None
        self._nuclides = None
        self._domain = None
        self._domain_type = None
        self._energy_groups = None
        self._tally_trigger = None
        self._tallies = None
        self._rxn_rate_tally = None
        self._xs_tally = None
        self._sparse = False
        self._derived = False

        self.name = name
        self.by_nuclide = by_nuclide

        if domain_type is not None:
            self.domain_type = domain_type
        if domain is not None:
            self.domain = domain
        if energy_groups is not None:
            self.energy_groups = energy_groups

    def __deepcopy__(self, memo):
        existing = memo.get(id(self))

        # If this is the first time we have tried to copy this object, copy it
        if existing is None:
            clone = type(self).__new__(type(self))
            clone._name = self.name
            clone._rxn_type = self.rxn_type
            clone._by_nuclide = self.by_nuclide
            clone._nuclides = copy.deepcopy(self._nuclides)
            clone._domain = self.domain
            clone._domain_type = self.domain_type
            clone._energy_groups = copy.deepcopy(self.energy_groups, memo)
            clone._tally_trigger = copy.deepcopy(self.tally_trigger, memo)
            clone._rxn_rate_tally = copy.deepcopy(self._rxn_rate_tally, memo)
            clone._xs_tally = copy.deepcopy(self._xs_tally, memo)
            clone._sparse = self.sparse
            clone._derived = self.derived

            clone._tallies = OrderedDict()
            for tally_type, tally in self.tallies.items():
                clone.tallies[tally_type] = copy.deepcopy(tally, memo)

            memo[id(self)] = clone

            return clone

        # If this object has been copied before, return the first copy made
        else:
            return existing

    @property
    def name(self):
        return self._name

    @property
    def rxn_type(self):
        return self._rxn_type

    @property
    def by_nuclide(self):
        return self._by_nuclide

    @property
    def domain(self):
        return self._domain

    @property
    def domain_type(self):
        return self._domain_type

    @property
    def energy_groups(self):
        return self._energy_groups

    @property
    def tally_trigger(self):
        return self._tally_trigger

    @property
    def num_groups(self):
        return self.energy_groups.num_groups

    @property
    def tallies(self):
        return self._tallies

    @property
    def rxn_rate_tally(self):
        return self._rxn_rate_tally

    @property
    def xs_tally(self):
        """Computes multi-group cross section using OpenMC tally arithmetic."""

        if self._xs_tally is None:
            if self.tallies is None:
                msg = 'Unable to get xs_tally since tallies have ' \
                      'not been loaded from a statepoint'
                raise ValueError(msg)

            self._xs_tally = self.rxn_rate_tally / self.tallies['flux']
            self._compute_xs()

        return self._xs_tally

    @property
    def sparse(self):
        return self._sparse

    @property
    def num_subdomains(self):
        domain_filter = self.xs_tally.find_filter(self.domain_type)
        return domain_filter.num_bins

    @property
    def num_nuclides(self):
        if self.by_nuclide:
            return len(self.get_all_nuclides())
        else:
            return 1

    @property
    def nuclides(self):
        if self.by_nuclide:
            return self.get_all_nuclides()
        else:
            return 'sum'

    @property
    def derived(self):
        return self._derived

    @name.setter
    def name(self, name):
        cv.check_type('name', name, basestring)
        self._name = name

    @by_nuclide.setter
    def by_nuclide(self, by_nuclide):
        cv.check_type('by_nuclide', by_nuclide, bool)
        self._by_nuclide = by_nuclide

    @nuclides.setter
    def nuclides(self, nuclides):
        cv.check_iterable_type('nuclides', nuclides, basestring)
        self._nuclides = nuclides

    @domain.setter
    def domain(self, domain):
        cv.check_type('domain', domain, tuple(_DOMAINS))
        self._domain = domain

    @domain_type.setter
    def domain_type(self, domain_type):
        cv.check_value('domain type', domain_type, tuple(DOMAIN_TYPES))
        self._domain_type = domain_type

    @energy_groups.setter
    def energy_groups(self, energy_groups):
        cv.check_type('energy groups', energy_groups, openmc.mgxs.EnergyGroups)
        self._energy_groups = energy_groups

    @tally_trigger.setter
    def tally_trigger(self, tally_trigger):
        cv.check_type('tally trigger', tally_trigger, openmc.Trigger)
        self._tally_trigger = tally_trigger

    @sparse.setter
    def sparse(self, sparse):
        """Convert tally data from NumPy arrays to SciPy list of lists (LIL)
        sparse matrices, and vice versa.

        This property may be used to reduce the amount of data in memory during
        tally data processing. The tally data will be stored as SciPy LIL
        matrices internally within the Tally object. All tally data access
        properties and methods will return data as a dense NumPy array.

        """

        cv.check_type('sparse', sparse, bool)

        # Sparsify or densify the derived MGXS tallies and the base tallies
        if self._xs_tally:
            self.xs_tally.sparse = sparse
        if self._rxn_rate_tally:
            self.rxn_rate_tally.sparse = sparse

        for tally_name in self.tallies:
                self.tallies[tally_name].sparse = sparse

        self._sparse = sparse

    @staticmethod
    def get_mgxs(mgxs_type, domain=None, domain_type=None,
                 energy_groups=None, by_nuclide=False, name=''):
        """Return a MGXS subclass object for some energy group structure within
        some spatial domain for some reaction type.

        This is a factory method which can be used to quickly create MGXS
        subclass objects for various reaction types.

        Parameters
        ----------
        mgxs_type : {'total', 'transport', 'absorption', 'capture', 'fission', 'nu-fission', 'kappa-fission', 'scatter', 'nu-scatter', 'scatter matrix', 'nu-scatter matrix', 'chi'}
            The type of multi-group cross section object to return
        domain : Material or Cell or Universe
            The domain for spatial homogenization
        domain_type : {'material', 'cell', 'distribcell', 'universe'}
            The domain type for spatial homogenization
        energy_groups : EnergyGroups
            The energy group structure for energy condensation
        by_nuclide : bool
            If true, computes cross sections for each nuclide in domain.
            Defaults to False
        name : str, optional
            Name of the multi-group cross section. Used as a label to identify
            tallies in OpenMC 'tallies.xml' file. Defaults to the empty string.

        Returns
        -------
        MGXS
            A subclass of the abstract MGXS class for the multi-group cross
            section type requested by the user

        """

        cv.check_value('mgxs_type', mgxs_type, MGXS_TYPES)

        if mgxs_type == 'total':
            mgxs = TotalXS(domain, domain_type, energy_groups)
        elif mgxs_type == 'transport':
            mgxs = TransportXS(domain, domain_type, energy_groups)
        elif mgxs_type == 'absorption':
            mgxs = AbsorptionXS(domain, domain_type, energy_groups)
        elif mgxs_type == 'capture':
            mgxs = CaptureXS(domain, domain_type, energy_groups)
        elif mgxs_type == 'fission':
            mgxs = FissionXS(domain, domain_type, energy_groups)
        elif mgxs_type == 'nu-fission':
            mgxs = NuFissionXS(domain, domain_type, energy_groups)
        elif mgxs_type == 'kappa-fission':
            mgxs = KappaFissionXS(domain, domain_type, energy_groups)
        elif mgxs_type == 'scatter':
            mgxs = ScatterXS(domain, domain_type, energy_groups)
        elif mgxs_type == 'nu-scatter':
            mgxs = NuScatterXS(domain, domain_type, energy_groups)
        elif mgxs_type == 'scatter matrix':
            mgxs = ScatterMatrixXS(domain, domain_type, energy_groups)
        elif mgxs_type == 'nu-scatter matrix':
            mgxs = NuScatterMatrixXS(domain, domain_type, energy_groups)
        elif mgxs_type == 'chi':
            mgxs = Chi(domain, domain_type, energy_groups)

        mgxs.by_nuclide = by_nuclide
        mgxs.name = name
        return mgxs

    def get_all_nuclides(self):
        """Get all nuclides in the cross section's spatial domain.

        Returns
        -------
        list of str
            A list of the string names for each nuclide in the spatial domain
            (e.g., ['U-235', 'U-238', 'O-16'])

        Raises
        ------
        ValueError
            When this method is called before the spatial domain has been set.

        """

        if self.domain is None:
            raise ValueError('Unable to get all nuclides without a domain')

        # If the user defined nuclides, return them
        if self._nuclides:
            return self._nuclides

        # Otherwise, return all nuclides in the spatial domain
        else:
            nuclides = self.domain.get_all_nuclides()
            return nuclides.keys()

    def get_nuclide_density(self, nuclide):
        """Get the atomic number density in units of atoms/b-cm for a nuclide
        in the cross section's spatial domain.

        Parameters
        ----------
        nuclide : str
            A nuclide name string (e.g., 'U-235')

        Returns
        -------
        Real
            The atomic number density (atom/b-cm) for the nuclide of interest

        Raises
        -------
        ValueError
            When the density is requested for a nuclide which is not found in
            the spatial domain.

        """

        cv.check_type('nuclide', nuclide, basestring)

        # Get list of all nuclides in the spatial domain
        nuclides = self.domain.get_all_nuclides()

        if nuclide not in nuclides:
            msg = 'Unable to get density for nuclide "{0}" which is not in ' \
                  '{1} "{2}"'.format(nuclide, self.domain_type, self.domain.id)
            ValueError(msg)

        density = nuclides[nuclide][1]
        return density

    def get_nuclide_densities(self, nuclides='all'):
        """Get an array of atomic number densities in units of atom/b-cm for all
        nuclides in the cross section's spatial domain.

        Parameters
        ----------
        nuclides : Iterable of str or 'all' or 'sum'
            A list of nuclide name strings (e.g., ['U-235', 'U-238']). The
            special string 'all' will return the atom densities for all nuclides
            in the spatial domain. The special string 'sum' will return the atom
            density summed across all nuclides in the spatial domain. Defaults
            to 'all'.

        Returns
        -------
        ndarray of Real
            An array of the atomic number densities (atom/b-cm) for each of the
            nuclides in the spatial domain

        Raises
        ------
        ValueError
            When this method is called before the spatial domain has been set.

        """

        if self.domain is None:
            raise ValueError('Unable to get nuclide densities without a domain')

        # Sum the atomic number densities for all nuclides
        if nuclides == 'sum':
            nuclides = self.get_all_nuclides()
            densities = np.zeros(1, dtype=np.float)
            for nuclide in nuclides:
                densities[0] += self.get_nuclide_density(nuclide)

        # Tabulate the atomic number densities for all nuclides
        elif nuclides == 'all':
            nuclides = self.get_all_nuclides()
            densities = np.zeros(self.num_nuclides, dtype=np.float)
            for i, nuclide in enumerate(nuclides):
                densities[i] += self.get_nuclide_density(nuclide)

        # Tabulate the atomic number densities for each specified nuclide
        else:
            densities = np.zeros(len(nuclides), dtype=np.float)
            for i, nuclide in enumerate(nuclides):
                densities[i] = self.get_nuclide_density(nuclide)

        return densities

    def _create_tallies(self, scores, all_filters, keys, estimator):
        """Instantiates tallies needed to compute the multi-group cross section.

        This is a helper method for MGXS subclasses to create tallies
        for input file generation. The tallies are stored in the tallies dict.
        This method is called by each subclass' tallies property getter
        which define the parameters given to this parent class method.

        Parameters
        ----------
        scores : Iterable of str
            Scores for each tally
        all_filters : Iterable of tuple of Filter
            Tuples of non-spatial domain filters for each tally
        keys : Iterable of str
            Key string used to store each tally in the tallies dictionary
        estimator : {'analog' or 'tracklength'}
            Type of estimator to use for each tally

        """

        cv.check_iterable_type('scores', scores, basestring)
        cv.check_length('scores', scores, len(keys))
        cv.check_iterable_type('filters', all_filters, openmc.Filter, 1, 2)
        cv.check_type('keys', keys, Iterable, basestring)
        cv.check_value('estimator', estimator, ['analog', 'tracklength'])

        self._tallies = OrderedDict()

        # Create a domain Filter object
        domain_filter = openmc.Filter(self.domain_type, self.domain.id)

        # Create each Tally needed to compute the multi group cross section
        for score, key, filters in zip(scores, keys, all_filters):
            self.tallies[key] = openmc.Tally(name=self.name)
            self.tallies[key].add_score(score)
            self.tallies[key].estimator = estimator
            self.tallies[key].add_filter(domain_filter)

            # If a tally trigger was specified, add it to each tally
            if self.tally_trigger:
                trigger_clone = copy.deepcopy(self.tally_trigger)
                trigger_clone.add_score(score)
                self.tallies[key].add_trigger(trigger_clone)

            # Add all non-domain specific Filters (e.g., 'energy') to the Tally
            for add_filter in filters:
                self.tallies[key].add_filter(add_filter)

            # If this is a by-nuclide cross-section, add all nuclides to Tally
            if self.by_nuclide and score != 'flux':
                all_nuclides = self.domain.get_all_nuclides()
                for nuclide in all_nuclides:
                    self.tallies[key].add_nuclide(nuclide)
            else:
                self.tallies[key].add_nuclide('total')

    def _compute_xs(self):
        """Performs generic cleanup after a subclass' uses tally arithmetic to
        compute a multi-group cross section as a derived tally.

        This method replaces CrossNuclides generated by tally arithmetic with
        the original Nuclide objects in the xs_tally instance attribute. The
        simple Nuclides allow for cleaner output through Pandas DataFrames as
        well as simpler data access through the get_xs(...) class method.

        In addition, this routine resets NaNs in the multi group cross section
        array to 0.0. This may be needed occur if no events were scored in
        certain tally bins, which will lead to a divide-by-zero situation.

        """

        # If computing xs for each nuclide, replace CrossNuclides with originals
        if self.by_nuclide:
            self.xs_tally._nuclides = []
            nuclides = self.get_all_nuclides()
            for nuclide in nuclides:
                self.xs_tally.add_nuclide(openmc.Nuclide(nuclide))

        # Remove NaNs which may have resulted from divide-by-zero operations
        self.xs_tally._mean = np.nan_to_num(self.xs_tally.mean)
        self.xs_tally._std_dev = np.nan_to_num(self.xs_tally.std_dev)
        self.xs_tally.sparse = self.sparse

    def load_from_statepoint(self, statepoint):
        """Extracts tallies in an OpenMC StatePoint with the data needed to
        compute multi-group cross sections.

        This method is needed to compute cross section data from tallies
        in an OpenMC StatePoint object.

        NOTE: The statepoint must first be linked with an OpenMC Summary object.

        Parameters
        ----------
        statepoint : openmc.StatePoint
            An OpenMC StatePoint object with tally data

        Raises
        ------
        ValueError
            When this method is called with a statepoint that has not been
            linked with a summary object.

        """

        cv.check_type('statepoint', statepoint, openmc.statepoint.StatePoint)

        if statepoint.summary is None:
            msg = 'Unable to load data from a statepoint which has not been ' \
                  'linked with a summary file'
            raise ValueError(msg)

        # Override the domain object that loaded from an OpenMC summary file
        # NOTE: This is necessary for micro cross-sections which require
        # the isotopic number densities as computed by OpenMC
        if self.domain_type == 'cell' or self.domain_type == 'distribcell':
            self.domain = statepoint.summary.get_cell_by_id(self.domain.id)
        elif self.domain_type == 'universe':
            self.domain = statepoint.summary.get_universe_by_id(self.domain.id)
        elif self.domain_type == 'material':
            self.domain = statepoint.summary.get_material_by_id(self.domain.id)
        else:
            msg = 'Unable to load data from a statepoint for domain type {0} ' \
                  'which is not yet supported'.format(self.domain_type)
            raise ValueError(msg)

        # Use tally "slicing" to ensure that tallies correspond to our domain
        # NOTE: This is important if tally merging was used
        if self.domain_type != 'distribcell':
            filters = [self.domain_type]
            filter_bins = [(self.domain.id,)]
        # Distribcell filters only accept single cell - neglect it when slicing
        else:
            filters = []
            filter_bins = []

        # Clear any tallies previously loaded from a statepoint
        self._tallies = None
        self._xs_tally = None
        self._rxn_rate_tally = None

        # Find, slice and store Tallies from StatePoint
        # The tally slicing is needed if tally merging was used
        for tally_type, tally in self.tallies.items():
            sp_tally = statepoint.get_tally(tally.scores, tally.filters,
                                            tally.nuclides,
                                            estimator=tally.estimator)
            sp_tally = sp_tally.get_slice(tally.scores, filters,
                                          filter_bins, tally.nuclides)
            sp_tally.sparse = self.sparse
            self.tallies[tally_type] = sp_tally

    def get_xs(self, groups='all', subdomains='all', nuclides='all',
               xs_type='macro', order_groups='increasing', value='mean'):
        """Returns an array of multi-group cross sections.

        This method constructs a 2D NumPy array for the requested multi-group
        cross section data data for one or more energy groups and subdomains.

        Parameters
        ----------
        groups : Iterable of Integral or 'all'
            Energy groups of interest. Defaults to 'all'.
        subdomains : Iterable of Integral or 'all'
            Subdomain IDs of interest. Defaults to 'all'.
        nuclides : Iterable of str or 'all' or 'sum'
            A list of nuclide name strings (e.g., ['U-235', 'U-238']). The
            special string 'all' will return the cross sections for all nuclides
            in the spatial domain. The special string 'sum' will return the
            cross section summed over all nuclides. Defaults to 'all'.
        xs_type: {'macro', 'micro'}
            Return the macro or micro cross section in units of cm^-1 or barns.
            Defaults to 'macro'.
        order_groups: {'increasing', 'decreasing'}
            Return the cross section indexed according to increasing or
            decreasing energy groups (decreasing or increasing energies).
            Defaults to 'increasing'.
        value : str
            A string for the type of value to return - 'mean', 'std_dev' or
            'rel_err' are accepted. Defaults to 'mean'.

        Returns
        -------
        ndarray
            A NumPy array of the multi-group cross section indexed in the order
            each group, subdomain and nuclide is listed in the parameters.

        Raises
        ------
        ValueError
            When this method is called before the multi-group cross section is
            computed from tally data.

        """

        cv.check_value('value', value, ['mean', 'std_dev', 'rel_err'])
        cv.check_value('xs_type', xs_type, ['macro', 'micro'])

        filters = []
        filter_bins = []

        # Construct a collection of the domain filter bins
        if not isinstance(subdomains, basestring):
            cv.check_iterable_type('subdomains', subdomains, Integral, max_depth=2)
            for subdomain in subdomains:
                filters.append(self.domain_type)
                filter_bins.append((subdomain,))

        # Construct list of energy group bounds tuples for all requested groups
        if not isinstance(groups, basestring):
            cv.check_iterable_type('groups', groups, Integral)
            for group in groups:
                filters.append('energy')
                filter_bins.append((self.energy_groups.get_group_bounds(group),))

        # Construct a collection of the nuclides to retrieve from the xs tally
        if self.by_nuclide:
            if nuclides == 'all' or nuclides == 'sum' or nuclides == ['sum']:
                query_nuclides = self.get_all_nuclides()
            else:
                query_nuclides = nuclides
        else:
            query_nuclides = ['total']

        # If user requested the sum for all nuclides, use tally summation
        if nuclides == 'sum' or nuclides == ['sum']:
            xs_tally = self.xs_tally.summation(nuclides=query_nuclides)
            xs = xs_tally.get_values(filters=filters,
                                     filter_bins=filter_bins, value=value)
        else:
            xs = self.xs_tally.get_values(filters=filters, filter_bins=filter_bins,
                                          nuclides=query_nuclides, value=value)

        # Divide by atom number densities for microscopic cross sections
        if xs_type == 'micro':
            if self.by_nuclide:
                densities = self.get_nuclide_densities(nuclides)
            else:
                densities = self.get_nuclide_densities('sum')
            if value == 'mean' or value == 'std_dev':
                xs /= densities[np.newaxis, :, np.newaxis]

        # Reverse data if user requested increasing energy groups since
        # tally data is stored in order of increasing energies
        if order_groups == 'increasing':
            if groups == 'all':
                num_groups = self.num_groups
            else:
                num_groups = len(groups)

            # Reshape tally data array with separate axes for domain and energy
            num_subdomains = int(xs.shape[0] / num_groups)
            new_shape = (num_subdomains, num_groups) + xs.shape[1:]
            xs = np.reshape(xs, new_shape)

            # Reverse energies to align with increasing energy groups
            xs = xs[:, ::-1, :]

            # Eliminate trivial dimensions
            xs = np.squeeze(xs)
            xs = np.atleast_1d(xs)

        return xs

    def get_condensed_xs(self, coarse_groups):
        """Construct an energy-condensed version of this cross section.

        Parameters
        ----------
        coarse_groups : openmc.mgxs.EnergyGroups
            The coarse energy group structure of interest

        Returns
        -------
        MGXS
            A new MGXS condensed to the group structure of interest

        """

        cv.check_type('coarse_groups', coarse_groups, EnergyGroups)
        cv.check_less_than('coarse groups', coarse_groups.num_groups,
                           self.num_groups, equality=True)
        cv.check_value('upper coarse energy', coarse_groups.group_edges[-1],
                       [self.energy_groups.group_edges[-1]])
        cv.check_value('lower coarse energy', coarse_groups.group_edges[0],
                       [self.energy_groups.group_edges[0]])

        # Clone this MGXS to initialize the condensed version
        condensed_xs = copy.deepcopy(self)
        condensed_xs._rxn_rate_tally = None
        condensed_xs._xs_tally = None
        condensed_xs._sparse = False
        condensed_xs._energy_groups = coarse_groups

        # Build energy indices to sum across
        energy_indices = []
        for group in range(coarse_groups.num_groups, 0, -1):
            low, high = coarse_groups.get_group_bounds(group)
            low_index = np.where(self.energy_groups.group_edges == low)[0][0]
            energy_indices.append(low_index)

        fine_edges = self.energy_groups.group_edges

        # Condense each of the tallies to the coarse group structure
        for tally_type, tally in condensed_xs.tallies.items():

            # Make condensed tally derived and null out sum, sum_sq
            tally._derived = True
            tally._sum = None
            tally._sum_sq = None

            # Get tally data arrays reshaped with one dimension per filter
            mean = tally.get_reshaped_data(value='mean')
            std_dev = tally.get_reshaped_data(value='std_dev')

            # Sum across all applicable fine energy group filters
            for i, tally_filter in enumerate(tally.filters):
                if 'energy' not in tally_filter.type:
                    continue
                elif len(tally_filter.bins) != len(fine_edges):
                    continue
                elif not np.allclose(tally_filter.bins, fine_edges):
                    continue
                else:
                    tally_filter.bins = coarse_groups.group_edges
                    mean = np.add.reduceat(mean, energy_indices, axis=i)
                    std_dev = np.add.reduceat(std_dev**2, energy_indices, axis=i)
                    std_dev = np.sqrt(std_dev)

            # Reshape condensed data arrays with one dimension for all filters
            mean = np.reshape(mean, tally.shape)
            std_dev = np.reshape(std_dev, tally.shape)

            # Override tally's data with the new condensed data
            tally._mean = mean
            tally._std_dev = std_dev

        # Compute the energy condensed multi-group cross section
        condensed_xs.sparse = self.sparse
        return condensed_xs

    def get_subdomain_avg_xs(self, subdomains='all'):
        """Construct a subdomain-averaged version of this cross section.

        This method is useful for averaging cross sections across distribcell
        instances. The method performs spatial homogenization to compute the
        scalar flux-weighted average cross section across the subdomains.

        Parameters
        ----------
        subdomains : Iterable of Integral or 'all'
            The subdomain IDs to average across. Defaults to 'all'.

        Returns
        -------
        MGXS
            A new MGXS averaged across the subdomains of interest

        Raises
        ------
        ValueError
            When this method is called before the multi-group cross section is
            computed from tally data.

        """

        # Construct a collection of the subdomain filter bins to average across
        if not isinstance(subdomains, basestring):
            cv.check_iterable_type('subdomains', subdomains, Integral)
        elif self.domain_type == 'distribcell':
            subdomains = np.arange(self.num_subdomains)
        else:
            subdomains = None

        # Clone this MGXS to initialize the subdomain-averaged version
        avg_xs = copy.deepcopy(self)

        if self.derived:
            avg_xs._rxn_rate_tally = avg_xs.rxn_rate_tally.average(
                filter_type=self.domain_type, filter_bins=subdomains)
        else:
            avg_xs._rxn_rate_tally = None
            avg_xs._xs_tally = None

            # Average each of the tallies across subdomains
            for tally_type, tally in avg_xs.tallies.items():
                tally_avg = tally.average(filter_type=self.domain_type,
                                          filter_bins=subdomains)
                avg_xs.tallies[tally_type] = tally_avg

        avg_xs._domain_type = 'avg({0})'.format(self.domain_type)
        avg_xs.sparse = self.sparse
        return avg_xs

    def get_slice(self, nuclides=[], groups=[]):
        """Build a sliced MGXS for the specified nuclides and energy groups.

        This method constructs a new MGXS to encapsulate a subset of the data
        represented by this MGXS. The subset of data to include in the tally
        slice is determined by the nuclides and energy groups specified in
        the input parameters.

        Parameters
        ----------
        nuclides : list of str
            A list of nuclide name strings
            (e.g., ['U-235', 'U-238']; default is [])
        groups : list of Integral
            A list of energy group indices starting at 1 for the high energies
            (e.g., [1, 2, 3]; default is [])

        Returns
        -------
        MGXS
            A new tally which encapsulates the subset of data requested for the
            nuclide(s) and/or energy group(s) requested in the parameters.

        """

        cv.check_iterable_type('nuclides', nuclides, basestring)
        cv.check_iterable_type('energy_groups', groups, Integral)

        # Build lists of filters and filter bins to slice
        if len(groups) == 0:
            filters = []
            filter_bins = []
        else:
            filter_bins = []
            for group in groups:
                group_bounds = self.energy_groups.get_group_bounds(group)
                filter_bins.append(group_bounds)
            filter_bins = [tuple(filter_bins)]
            filters = ['energy']

        # Clone this MGXS to initialize the sliced version
        slice_xs = copy.deepcopy(self)
        slice_xs._rxn_rate_tally = None
        slice_xs._xs_tally = None

        # Slice each of the tallies across nuclides and energy groups
        for tally_type, tally in slice_xs.tallies.items():
            slice_nuclides = [nuc for nuc in nuclides if nuc in tally.nuclides]
            if len(groups) != 0 and tally.contains_filter('energy'):
                tally_slice = tally.get_slice(filters=filters,
                    filter_bins=filter_bins, nuclides=slice_nuclides)
            else:
                tally_slice = tally.get_slice(nuclides=slice_nuclides)
            slice_xs.tallies[tally_type] = tally_slice

        # Assign sliced energy group structure to sliced MGXS
        if groups:
            new_group_edges = []
            for group in groups:
                group_edges = self.energy_groups.get_group_bounds(group)
                new_group_edges.extend(group_edges)
            new_group_edges = np.unique(new_group_edges)
            slice_xs.energy_groups.group_edges = sorted(new_group_edges)

        # Assign sliced nuclides to sliced MGXS
        if nuclides:
            slice_xs.nuclides = nuclides

        slice_xs.sparse = self.sparse
        return slice_xs

    def can_merge(self, other):
        """Determine if another MGXS can be merged with this one

        If results have been loaded from a statepoint, then MGXS are only
        mergeable along one and only one of enegy groups or nuclides.

        Parameters
        ----------
        other : MGXS
            MGXS to check for merging

        """

        if not isinstance(other, type(self)):
            return False

        # Compare reaction type, energy groups, nuclides, domain type
        if self.rxn_type != other.rxn_type:
            return False
        elif not self.energy_groups.can_merge(other.energy_groups):
            return False
        elif self.by_nuclide != other.by_nuclide:
            return False
        elif self.domain_type != other.domain_type:
            return False
        elif 'distribcell' not in self.domain_type and self.domain != other.domain:
            return False
        elif not self.xs_tally.can_merge(other.xs_tally):
            return False
        elif not self.rxn_rate_tally.can_merge(other.rxn_rate_tally):
            return False

        # If all conditionals pass then MGXS are mergeable
        return True

    def merge(self, other):
        """Merge another MGXS with this one

        MGXS are only mergeable if their energy groups and nuclides are either
        identical or mutually exclusive. If results have been loaded from a
        statepoint, then MGXS are only mergeable along one and only one of
        energy groups or nuclides.

        Parameters
        ----------
        other : MGXS
            MGXS to merge with this one

        Returns
        -------
        merged_mgxs : MGXS
            Merged MGXS

        """

        if not self.can_merge(other):
            raise ValueError('Unable to merge MGXS')

        # Create deep copy of tally to return as merged tally
        merged_mgxs = copy.deepcopy(self)
        merged_mgxs._derived = True

        # Merge energy groups
        if self.energy_groups != other.energy_groups:
            merged_groups = self.energy_groups.merge(other.energy_groups)
            merged_mgxs.energy_groups = merged_groups

        # Merge nuclides
        if self.nuclides != other.nuclides:

            # The nuclides must be mutually exclusive
            for nuclide in self.nuclides:
                if nuclide in other.nuclides:
                    msg = 'Unable to merge MGXS with shared nuclides'
                    raise ValueError(msg)

            # Concatenate lists of nuclides for the merged MGXS
            merged_mgxs.nuclides = self.nuclides + other.nuclides

        # Null base tallies but merge reaction rate and cross section tallies
        merged_mgxs._tallies = OrderedDict()
        merged_mgxs._rxn_rate_tally = self.rxn_rate_tally.merge(other.rxn_rate_tally)
        merged_mgxs._xs_tally = self.xs_tally.merge(other.xs_tally)

        return merged_mgxs

    def print_xs(self, subdomains='all', nuclides='all', xs_type='macro'):
        """Print a string representation for the multi-group cross section.

        Parameters
        ----------
        subdomains : Iterable of Integral or 'all'
            The subdomain IDs of the cross sections to include in the report.
            Defaults to 'all'.
        nuclides : Iterable of str or 'all' or 'sum'
            The nuclides of the cross-sections to include in the report. This
            may be a list of nuclide name strings (e.g., ['U-235', 'U-238']).
            The special string 'all' will report the cross sections for all
            nuclides in the spatial domain. The special string 'sum' will report
            the cross sections summed over all nuclides. Defaults to 'all'.
        xs_type: {'macro', 'micro'}
            Return the macro or micro cross section in units of cm^-1 or barns.
            Defaults to 'macro'.

        """

        # Construct a collection of the subdomains to report
        if not isinstance(subdomains, basestring):
            cv.check_iterable_type('subdomains', subdomains, Integral)
        elif self.domain_type == 'distribcell':
            subdomains = np.arange(self.num_subdomains, dtype=np.int)
        else:
            subdomains = [self.domain.id]

        # Construct a collection of the nuclides to report
        if self.by_nuclide:
            if nuclides == 'all':
                nuclides = self.get_all_nuclides()
            elif nuclides == 'sum':
                nuclides = ['sum']
            else:
                cv.check_iterable_type('nuclides', nuclides, basestring)
        else:
            nuclides = ['sum']

        cv.check_value('xs_type', xs_type, ['macro', 'micro'])

        # Build header for string with type and domain info
        string = 'Multi-Group XS\n'
        string += '{0: <16}=\t{1}\n'.format('\tReaction Type', self.rxn_type)
        string += '{0: <16}=\t{1}\n'.format('\tDomain Type', self.domain_type)
        string += '{0: <16}=\t{1}\n'.format('\tDomain ID', self.domain.id)

        # If cross section data has not been computed, only print string header
        if self.tallies is None:
            print(string)
            return

        # Loop over all subdomains
        for subdomain in subdomains:

            if self.domain_type == 'distribcell':
                string += '{0: <16}=\t{1}\n'.format('\tSubdomain', subdomain)

            # Loop over all Nuclides
            for nuclide in nuclides:

                # Build header for nuclide type
                if nuclide != 'sum':
                    string += '{0: <16}=\t{1}\n'.format('\tNuclide', nuclide)

                # Build header for cross section type
                if xs_type == 'macro':
                    string += '{0: <16}\n'.format('\tCross Sections [cm^-1]:')
                else:
                    string += '{0: <16}\n'.format('\tCross Sections [barns]:')

                template = '{0: <12}Group {1} [{2: <10} - {3: <10}MeV]:\t'

                # Loop over energy groups ranges
                for group in range(1, self.num_groups+1):
                    bounds = self.energy_groups.get_group_bounds(group)
                    string += template.format('', group, bounds[0], bounds[1])
                    average = self.get_xs([group], [subdomain], [nuclide],
                                          xs_type=xs_type, value='mean')
                    rel_err = self.get_xs([group], [subdomain], [nuclide],
                                          xs_type=xs_type, value='rel_err')
                    average = average.flatten()[0]
                    rel_err = rel_err.flatten()[0] * 100.
                    string += '{:.2e} +/- {:1.2e}%'.format(average, rel_err)
                    string += '\n'
                string += '\n'
            string += '\n'

        print(string)

    def build_hdf5_store(self, filename='mgxs.h5', directory='mgxs',
                         subdomains='all', nuclides='all',
                         xs_type='macro', append=True):
        """Export the multi-group cross section data to an HDF5 binary file.

        This method constructs an HDF5 file which stores the multi-group
        cross section data. The data is stored in a hierarchy of HDF5 groups
        from the domain type, domain id, subdomain id (for distribcell domains),
        nuclides and cross section type. Two datasets for the mean and standard
        deviation are stored for each subdomain entry in the HDF5 file.

        NOTE: This requires the h5py Python package.

        Parameters
        ----------
        filename : str
            Filename for the HDF5 file. Defaults to 'mgxs.h5'.
        directory : str
            Directory for the HDF5 file. Defaults to 'mgxs'.
        subdomains : Iterable of Integral or 'all'
            The subdomain IDs of the cross sections to include in the report.
            Defaults to 'all'.
        nuclides : Iterable of str or 'all' or 'sum'
            The nuclides of the cross-sections to include in the report. This
            may be a list of nuclide name strings (e.g., ['U-235', 'U-238']).
            The special string 'all' will report the cross sections for all
            nuclides in the spatial domain. The special string 'sum' will report
            the cross sections summed over all nuclides. Defaults to 'all'.
        xs_type: {'macro', 'micro'}
            Store the macro or micro cross section in units of cm^-1 or barns.
            Defaults to 'macro'.
        append : bool
            If true, appends to an existing HDF5 file with the same filename
            directory (if one exists). Defaults to True.

        Raises
        ------
        ValueError
            When this method is called before the multi-group cross section is
            computed from tally data.
        ImportError
            When h5py is not installed.

        """

        import h5py

        # Make directory if it does not exist
        if not os.path.exists(directory):
            os.makedirs(directory)

        filename = os.path.join(directory, filename)
        filename = filename.replace(' ', '-')

        if append and os.path.isfile(filename):
            xs_results = h5py.File(filename, 'a')
        else:
            xs_results = h5py.File(filename, 'w')

        # Construct a collection of the subdomains to report
        if not isinstance(subdomains, basestring):
            cv.check_iterable_type('subdomains', subdomains, Integral)
        elif self.domain_type == 'distribcell':
            subdomains = np.arange(self.num_subdomains, dtype=np.int)
        elif self.domain_type == 'avg(distribcell)':
            domain_filter = self.xs_tally.find_filter('avg(distribcell)')
            subdomains = domain_filter.bins
        else:
            subdomains = [self.domain.id]

        # Construct a collection of the nuclides to report
        if self.by_nuclide:
            if nuclides == 'all':
                nuclides = self.get_all_nuclides()
                densities = np.zeros(len(nuclides), dtype=np.float)
            elif nuclides == 'sum':
                nuclides = ['sum']
            else:
                cv.check_iterable_type('nuclides', nuclides, basestring)
        else:
            nuclides = ['sum']

        cv.check_value('xs_type', xs_type, ['macro', 'micro'])

        # Create an HDF5 group within the file for the domain
        domain_type_group = xs_results.require_group(self.domain_type)
        domain_group = domain_type_group.require_group(str(self.domain.id))

        # Determine number of digits to pad subdomain group keys
        num_digits = len(str(self.num_subdomains))

        # Create a separate HDF5 group for each subdomain
        for i, subdomain in enumerate(subdomains):

            # Create an HDF5 group for the subdomain
            if self.domain_type == 'distribcell':
                group_name = ''.zfill(num_digits)
                subdomain_group = domain_group.require_group(group_name)
            else:
                subdomain_group = domain_group

            # Create a separate HDF5 group for the rxn type
            rxn_group = subdomain_group.require_group(self.rxn_type)

            # Create a separate HDF5 group for each nuclide
            for j, nuclide in enumerate(nuclides):

                if nuclide != 'sum':
                    density = densities[j]
                    nuclide_group = rxn_group.require_group(nuclide)
                    nuclide_group.require_dataset('density', dtype=np.float64,
                                                  data=[density], shape=(1,))
                else:
                    nuclide_group = rxn_group

                # Extract the cross section for this subdomain and nuclide
                average = self.get_xs(subdomains=[subdomain], nuclides=[nuclide],
                                      xs_type=xs_type, value='mean')
                std_dev = self.get_xs(subdomains=[subdomain], nuclides=[nuclide],
                                      xs_type=xs_type, value='std_dev')
                average = average.squeeze()
                std_dev = std_dev.squeeze()

                # Add MGXS results data to the HDF5 group
                nuclide_group.require_dataset('average', dtype=np.float64,
                                              shape=average.shape, data=average)
                nuclide_group.require_dataset('std. dev.', dtype=np.float64,
                                              shape=std_dev.shape, data=std_dev)

        # Close the results HDF5 file
        xs_results.close()

    def export_xs_data(self, filename='mgxs', directory='mgxs',
                       format='csv', groups='all', xs_type='macro'):
        """Export the multi-group cross section data to a file.

        This method leverages the functionality in the Pandas library to export
        the multi-group cross section data in a variety of output file formats
        for storage and/or post-processing.

        Parameters
        ----------
        filename : str
            Filename for the exported file. Defaults to 'mgxs'.
        directory : str
            Directory for the exported file. Defaults to 'mgxs'.
        format : {'csv', 'excel', 'pickle', 'latex'}
            The format for the exported data file. Defaults to 'csv'.
        groups : Iterable of Integral or 'all'
            Energy groups of interest. Defaults to 'all'.
        xs_type: {'macro', 'micro'}
            Store the macro or micro cross section in units of cm^-1 or barns.
            Defaults to 'macro'.

        """

        cv.check_type('filename', filename, basestring)
        cv.check_type('directory', directory, basestring)
        cv.check_value('format', format, ['csv', 'excel', 'pickle', 'latex'])
        cv.check_value('xs_type', xs_type, ['macro', 'micro'])

        # Make directory if it does not exist
        if not os.path.exists(directory):
            os.makedirs(directory)

        filename = os.path.join(directory, filename)
        filename = filename.replace(' ', '-')

        # Get a Pandas DataFrame for the data
        df = self.get_pandas_dataframe(groups=groups, xs_type=xs_type)

        # Capitalize column label strings
        df.columns = df.columns.astype(str)
        df.columns = map(str.title, df.columns)

        # Export the data using Pandas IO API
        if format == 'csv':
            df.to_csv(filename + '.csv', index=False)
        elif format == 'excel':
            df.to_excel(filename + '.xls', index=False)
        elif format == 'pickle':
            df.to_pickle(filename + '.pkl')
        elif format == 'latex':
            if self.domain_type == 'distribcell':
                msg = 'Unable to export distribcell multi-group cross section' \
                      'data to a LaTeX table'
                raise NotImplementedError(msg)

            df.to_latex(filename + '.tex', bold_rows=True,
                        longtable=True, index=False)

            # Surround LaTeX table with code needed to run pdflatex
            with open(filename + '.tex','r') as original:
                data = original.read()
            with open(filename + '.tex','w') as modified:
                modified.write(
                    '\\documentclass[preview, 12pt, border=1mm]{standalone}\n')
                modified.write('\\usepackage{caption}\n')
                modified.write('\\usepackage{longtable}\n')
                modified.write('\\usepackage{booktabs}\n')
                modified.write('\\begin{document}\n\n')
                modified.write(data)
                modified.write('\n\\end{document}')

    def get_pandas_dataframe(self, groups='all', nuclides='all',
                             xs_type='macro', summary=None):
        """Build a Pandas DataFrame for the MGXS data.

        This method leverages the Tally.get_pandas_dataframe(...) method, but
        renames the columns with terminology appropriate for cross section data.

        Parameters
        ----------
        groups : Iterable of Integral or 'all'
            Energy groups of interest. Defaults to 'all'.
        nuclides : Iterable of str or 'all' or 'sum'
            The nuclides of the cross-sections to include in the dataframe. This
            may be a list of nuclide name strings (e.g., ['U-235', 'U-238']).
            The special string 'all' will include the cross sections for all
            nuclides in the spatial domain. The special string 'sum' will
            include the cross sections summed over all nuclides. Defaults
            to 'all'.
        xs_type: {'macro', 'micro'}
            Return macro or micro cross section in units of cm^-1 or barns.
            Defaults to 'macro'.
        summary : None or Summary
            An optional Summary object to be used to construct columns for
            distribcell tally filters (default is None). The geometric
            information in the Summary object is embedded into a multi-index
            column with a geometric "path" to each distribcell intance.
            NOTE: This option requires the OpenCG Python package.

        Returns
        -------
        pandas.DataFrame
            A Pandas DataFrame for the cross section data.

        Raises
        ------
        ValueError
            When this method is called before the multi-group cross section is
            computed from tally data.

        """

        if not isinstance(groups, basestring):
            cv.check_iterable_type('groups', groups, Integral)
        if nuclides != 'all' and nuclides != 'sum':
            cv.check_iterable_type('nuclides', nuclides, basestring)
        cv.check_value('xs_type', xs_type, ['macro', 'micro'])

        # Get a Pandas DataFrame from the derived xs tally
        if self.by_nuclide and nuclides == 'sum':

            # Use tally summation to sum across all nuclides
            query_nuclides = self.get_all_nuclides()
            xs_tally = self.xs_tally.summation(nuclides=query_nuclides)
            df = xs_tally.get_pandas_dataframe(summary=summary)

            # Remove nuclide column since it is homogeneous and redundant
            df.drop('nuclide', axis=1, inplace=True)

        # If the user requested a specific set of nuclides
        elif self.by_nuclide and nuclides != 'all':
            xs_tally = self.xs_tally.get_slice(nuclides=nuclides)
            df = xs_tally.get_pandas_dataframe(summary=summary)

        # If the user requested all nuclides, keep nuclide column in dataframe
        else:
            df = self.xs_tally.get_pandas_dataframe(summary=summary)

        # Remove the score column since it is homogeneous and redundant
        if summary and 'distribcell' in self.domain_type:
            df = df.drop('score', level=0, axis=1)
        else:
            df = df.drop('score', axis=1)

        # Override energy groups bounds with indices
        all_groups = np.arange(self.num_groups, 0, -1, dtype=np.int)
        all_groups = np.repeat(all_groups, self.num_nuclides)
        if 'energy low [MeV]' in df and 'energyout low [MeV]' in df:
            df.rename(columns={'energy low [MeV]': 'group in'},
                      inplace=True)
            in_groups = np.tile(all_groups, self.num_subdomains)
            in_groups = np.repeat(in_groups, df.shape[0] / in_groups.size)
            df['group in'] = in_groups
            del df['energy high [MeV]']

            df.rename(columns={'energyout low [MeV]': 'group out'},
                      inplace=True)
            out_groups = np.tile(all_groups, df.shape[0] / all_groups.size)
            df['group out'] = out_groups
            del df['energyout high [MeV]']
            columns = ['group in', 'group out']

        elif 'energyout low [MeV]' in df:
            df.rename(columns={'energyout low [MeV]': 'group out'},
                      inplace=True)
            in_groups = np.tile(all_groups, self.num_subdomains)
            df['group out'] = in_groups
            del df['energyout high [MeV]']
            columns = ['group out']

        elif 'energy low [MeV]' in df:
            df.rename(columns={'energy low [MeV]': 'group in'}, inplace=True)
            in_groups = np.tile(all_groups, self.num_subdomains)
            df['group in'] = in_groups
            del df['energy high [MeV]']
            columns = ['group in']

        # Select out those groups the user requested
        if not isinstance(groups, basestring):
            if 'group in' in df:
                df = df[df['group in'].isin(groups)]
            if 'group out' in df:
                df = df[df['group out'].isin(groups)]

        # If user requested micro cross sections, divide out the atom densities
        if xs_type == 'micro':
            if self.by_nuclide:
                densities = self.get_nuclide_densities(nuclides)
            else:
                densities = self.get_nuclide_densities('sum')
            tile_factor = df.shape[0] / len(densities)
            df['mean'] /= np.tile(densities, tile_factor)
            df['std. dev.'] /= np.tile(densities, tile_factor)

        # Sort the dataframe by domain type id (e.g., distribcell id) and
        # energy groups such that data is from fast to thermal
        df.sort_values(by=[self.domain_type] + columns, inplace=True)
        return df


class TotalXS(MGXS):
    """A total multi-group cross section."""

    def __init__(self, domain=None, domain_type=None,
                 groups=None, by_nuclide=False, name=''):
        super(TotalXS, self).__init__(domain, domain_type,
                                      groups, by_nuclide, name)
        self._rxn_type = 'total'

    @property
    def tallies(self):
        """Construct the OpenMC tallies needed to compute this cross section.

        This method constructs two tracklength tallies to compute the 'flux'
        and 'total' reaction rates in the spatial domain and energy groups
        of interest.

        """

        # Instantiate tallies if they do not exist
        if self._tallies is None:

            # Create a list of scores for each Tally to be created
            scores = ['flux', 'total']
            estimator = 'tracklength'
            keys = scores

            # Create the non-domain specific Filters for the Tallies
            group_edges = self.energy_groups.group_edges
            energy_filter = openmc.Filter('energy', group_edges)
            filters = [[energy_filter], [energy_filter]]

            # Initialize the Tallies
            self._create_tallies(scores, filters, keys, estimator)

        return self._tallies

    @property
    def rxn_rate_tally(self):
        if self._rxn_rate_tally is None :
            self._rxn_rate_tally = self.tallies['total']
            self._rxn_rate_tally.sparse = self.sparse
        return self._rxn_rate_tally


class TransportXS(MGXS):
    """A transport-corrected total multi-group cross section."""

    def __init__(self, domain=None, domain_type=None,
                 groups=None, by_nuclide=False, name=''):
        super(TransportXS, self).__init__(domain, domain_type,
                                          groups, by_nuclide, name)
        self._rxn_type = 'transport'

    @property
    def tallies(self):
        """Construct the OpenMC tallies needed to compute this cross section.

        This method constructs three analog tallies to compute the 'flux',
        'total' and 'scatter-P1' reaction rates in the spatial domain and
        energy groups of interest.

        """

        # Instantiate tallies if they do not exist
        if self._tallies is None:

            # Create a list of scores for each Tally to be created
            scores = ['flux', 'total', 'scatter-P1']
            estimator = 'analog'
            keys = scores

            # Create the non-domain specific Filters for the Tallies
            group_edges = self.energy_groups.group_edges
            energy_filter = openmc.Filter('energy', group_edges)
            energyout_filter = openmc.Filter('energyout', group_edges)
            filters = [[energy_filter], [energy_filter], [energyout_filter]]

            # Initialize the Tallies
            self._create_tallies(scores, filters, keys, estimator)

        return self._tallies

    @property
    def rxn_rate_tally(self):
        if self._rxn_rate_tally is None:
            scatter_p1 = copy.deepcopy(self.tallies['scatter-P1'])

            # Use tally slicing to remove scatter-P0 data from scatter-P1 tally
            self.tallies['scatter-P1'] = \
                scatter_p1.get_slice(scores=['scatter-P1'])

            self.tallies['scatter-P1'].filters[-1].type = 'energy'
            self._rxn_rate_tally = \
                self.tallies['total'] - self.tallies['scatter-P1']
            self._rxn_rate_tally.sparse = self.sparse

        return self._rxn_rate_tally


class AbsorptionXS(MGXS):
    """An absorption multi-group cross section."""

    def __init__(self, domain=None, domain_type=None,
                 groups=None, by_nuclide=False, name=''):
        super(AbsorptionXS, self).__init__(domain, domain_type,
                                           groups, by_nuclide, name)
        self._rxn_type = 'absorption'

    @property
    def tallies(self):
        """Construct the OpenMC tallies needed to compute this cross section.

        This method constructs two tracklength tallies to compute the 'flux'
        and 'absorption' reaction rates in the spatial domain and energy
        groups of interest.

        """

        # Instantiate tallies if they do not exist
        if self._tallies is None:

            # Create a list of scores for each Tally to be created
            scores = ['flux', 'absorption']
            estimator = 'tracklength'
            keys = scores

            # Create the non-domain specific Filters for the Tallies
            group_edges = self.energy_groups.group_edges
            energy_filter = openmc.Filter('energy', group_edges)
            filters = [[energy_filter], [energy_filter]]

            # Initialize the Tallies
            self._create_tallies(scores, filters, keys, estimator)

        return self._tallies

    @property
    def rxn_rate_tally(self):
        if self._rxn_rate_tally is None:
            self._rxn_rate_tally = self.tallies['absorption']
            self._rxn_rate_tally.sparse = self.sparse
        return self._rxn_rate_tally


class CaptureXS(MGXS):
    """A capture multi-group cross section.

    The neutron capture reaction rate is defined as the difference between
    OpenMC's 'absorption' and 'fission' reaction rate score types. This includes
    not only radiative capture, but all forms of neutron disappearance aside
    from fission (e.g., MT > 100).

    """

    def __init__(self, domain=None, domain_type=None,
                 groups=None, by_nuclide=False, name=''):
        super(CaptureXS, self).__init__(domain, domain_type,
                                        groups, by_nuclide, name)
        self._rxn_type = 'capture'

    @property
    def tallies(self):
        """Construct the OpenMC tallies needed to compute this cross section.

        This method constructs two tracklength tallies to compute the 'flux'
        and 'capture' reaction rates in the spatial domain and energy
        groups of interest.

        """

        # Instantiate tallies if they do not exist
        if self._tallies is None:

            # Create a list of scores for each Tally to be created
            scores = ['flux', 'absorption', 'fission']
            estimator = 'tracklength'
            keys = scores

            # Create the non-domain specific Filters for the Tallies
            group_edges = self.energy_groups.group_edges
            energy_filter = openmc.Filter('energy', group_edges)
            filters = [[energy_filter], [energy_filter], [energy_filter]]

            # Initialize the Tallies
            self._create_tallies(scores, filters, keys, estimator)

        return self._tallies

    @property
    def rxn_rate_tally(self):
        if self._rxn_rate_tally is None:
            self._rxn_rate_tally = \
                self.tallies['absorption'] - self.tallies['fission']
            self._rxn_rate_tally.sparse = self.sparse
        return self._rxn_rate_tally

class FissionXSBase(MGXS):
    """A fission production multi-group cross section base class
       for NuFission and KappaFission
    """

    # This is an abstract class which cannot be instantiated
    __metaclass__ = abc.ABCMeta

    def __init__(self, rxn_type, domain=None, domain_type=None,
                 groups=None, by_nuclide=False, name=''):
        super(FissionXSBase, self).__init__(domain, domain_type,
                                            groups, by_nuclide, name)
        self._rxn_type = rxn_type

    @property
    def tallies(self):
        """Construct the OpenMC tallies needed to compute this cross section.

        This method constructs two tracklength tallies to compute the 'flux'
        and 'rxn_type' reaction rates in the spatial domain and energy
        groups of interest.

        """

        # Instantiate tallies if they do not exist
        if self._tallies is None:

            # Create a list of scores for each Tally to be created
            scores = ['flux', self._rxn_type]
            estimator = 'tracklength'
            keys = scores

            # Create the non-domain specific Filters for the Tallies
            group_edges = self.energy_groups.group_edges
            energy_filter = openmc.Filter('energy', group_edges)
            filters = [[energy_filter], [energy_filter]]

            # Initialize the Tallies
            self._create_tallies(scores, filters, keys, estimator)

        return self._tallies

    @property
    def rxn_rate_tally(self):
        if self._rxn_rate_tally is None:
            self._rxn_rate_tally = self.tallies[self._rxn_type]
            self._rxn_rate_tally.sparse = self.sparse
        return self._rxn_rate_tally


class FissionXS(FissionXSBase):
    """A fission multi-group cross section."""

    def __init__(self, domain=None, domain_type=None,
                 groups=None, by_nuclide=False, name=''):
        super(FissionXS, self).__init__('fission', domain, domain_type,
                                        groups, by_nuclide, name)


class NuFissionXS(FissionXSBase):
    """A fission production multi-group cross section."""

    def __init__(self, domain=None, domain_type=None,
                 groups=None, by_nuclide=False, name=''):
        super(NuFissionXS, self).__init__('nu-fission', domain, domain_type,
                                          groups, by_nuclide, name)

class KappaFissionXS(FissionXSBase):
    """A recoverable fission energy production rate multi-group cross section."""

    def __init__(self, domain=None, domain_type=None,
                 groups=None, by_nuclide=False, name=''):
        super(KappaFissionXS, self).__init__('kappa-fission', domain, domain_type,
                                             groups, by_nuclide, name)

class ScatterXS(MGXS):
    """A scatter multi-group cross section."""

    def __init__(self, domain=None, domain_type=None,
                 groups=None, by_nuclide=False, name=''):
        super(ScatterXS, self).__init__(domain, domain_type,
                                        groups, by_nuclide, name)
        self._rxn_type = 'scatter'

    @property
    def tallies(self):
        """Construct the OpenMC tallies needed to compute this cross section.

        This method constructs two tracklength tallies to compute the 'flux'
        and 'scatter' reaction rates in the spatial domain and energy
        groups of interest.

        """

        # Instantiate tallies if they do not exist
        if self._tallies is None:

            # Create a list of scores for each Tally to be created
            scores = ['flux', 'scatter']
            estimator = 'tracklength'
            keys = scores

            # Create the non-domain specific Filters for the Tallies
            group_edges = self.energy_groups.group_edges
            energy_filter = openmc.Filter('energy', group_edges)
            filters = [[energy_filter], [energy_filter]]

            # Intialize the Tallies
            self._create_tallies(scores, filters, keys, estimator)

        return self._tallies

    @property
    def rxn_rate_tally(self):
        if self._rxn_rate_tally is None:
            self._rxn_rate_tally = self.tallies['scatter']
            self._rxn_rate_tally.sparse = self.sparse
        return self._rxn_rate_tally


class NuScatterXS(MGXS):
    """A nu-scatter multi-group cross section."""

    def __init__(self, domain=None, domain_type=None,
                 groups=None, by_nuclide=False, name=''):
        super(NuScatterXS, self).__init__(domain, domain_type,
                                          groups, by_nuclide, name)
        self._rxn_type = 'nu-scatter'

    @property
    def tallies(self):
        """Construct the OpenMC tallies needed to compute this cross section.

        This method constructs two analog tallies to compute the 'flux'
        and 'nu-scatter' reaction rates in the spatial domain and energy
        groups of interest.

        """

        # Instantiate tallies if they do not exist
        if self._tallies is None:

            # Create a list of scores for each Tally to be created
            scores = ['flux', 'nu-scatter']
            estimator = 'analog'
            keys = scores

            # Create the non-domain specific Filters for the Tallies
            group_edges = self.energy_groups.group_edges
            energy_filter = openmc.Filter('energy', group_edges)
            filters = [[energy_filter], [energy_filter]]

            # Initialize the Tallies
            self._create_tallies(scores, filters, keys, estimator)

        return self._tallies

    @property
    def rxn_rate_tally(self):
        if self._rxn_rate_tally is None:
            self._rxn_rate_tally = self.tallies['nu-scatter']
            self._rxn_rate_tally.sparse = self.sparse
        return self._rxn_rate_tally


class ScatterMatrixXS(MGXS):
    """A scattering matrix multi-group cross section.

    Attributes
    ----------
    correction : 'P0' or None
        Apply the P0 correction to scattering matrices if set to 'P0'

    """

    def __init__(self, domain=None, domain_type=None,
                 groups=None, by_nuclide=False, name=''):
        super(ScatterMatrixXS, self).__init__(domain, domain_type,
                                              groups, by_nuclide, name)
        self._rxn_type = 'scatter matrix'
        self._correction = 'P0'

    def __deepcopy__(self, memo):
        clone = super(ScatterMatrixXS, self).__deepcopy__(memo)
        clone._correction = self.correction
        return clone

    @property
    def correction(self):
        return self._correction

    @property
    def tallies(self):
        """Construct the OpenMC tallies needed to compute this cross section.

        This method constructs three analog tallies to compute the 'flux',
        'scatter' and 'scatter-P1' reaction rates in the spatial domain and
        energy groups of interest.

        """

        # Instantiate tallies if they do not exist
        if self._tallies is None:

            group_edges = self.energy_groups.group_edges
            energy = openmc.Filter('energy', group_edges)
            energyout = openmc.Filter('energyout', group_edges)

            # Create a list of scores for each Tally to be created
            if self.correction == 'P0':
                scores = ['flux', 'scatter', 'scatter-P1']
                filters = [[energy], [energy, energyout], [energyout]]
            else:
                scores = ['flux', 'scatter']
                filters = [[energy], [energy, energyout]]

            estimator = 'analog'
            keys = scores

            # Initialize the Tallies
            self._create_tallies(scores, filters, keys, estimator)

        return self._tallies

    @property
    def rxn_rate_tally(self):

        if self._rxn_rate_tally is None:
            # If using P0 correction subtract scatter-P1 from the diagonal
            if self.correction == 'P0':
                scatter_p1 = self.tallies['scatter-P1']
                scatter_p1 = scatter_p1.get_slice(scores=['scatter-P1'])
                energy_filter = self.tallies['scatter'].find_filter('energy')
                energy_filter = copy.deepcopy(energy_filter)
                scatter_p1 = scatter_p1.diagonalize_filter(energy_filter)
                self._rxn_rate_tally = self.tallies['scatter'] - scatter_p1
            else:
                self._rxn_rate_tally = self.tallies['scatter']

            self._rxn_rate_tally.sparse = self.sparse

        return self._rxn_rate_tally

    @correction.setter
    def correction(self, correction):
        cv.check_value('correction', correction, ('P0', None))
        self._correction = correction

    def get_slice(self, nuclides=[], in_groups=[], out_groups=[]):
        """Build a sliced ScatterMatrix for the specified nuclides and
        energy groups.

        This method constructs a new MGXS to encapsulate a subset of the data
        represented by this MGXS. The subset of data to include in the tally
        slice is determined by the nuclides and energy groups specified in
        the input parameters.

        Parameters
        ----------
        nuclides : list of str
            A list of nuclide name strings
            (e.g., ['U-235', 'U-238']; default is [])
        in_groups : list of Integral
            A list of incoming energy group indices starting at 1 for the high
            energies (e.g., [1, 2, 3]; default is [])
        out_groups : list of Integral
            A list of outgoing energy group indices starting at 1 for the high
            energies (e.g., [1, 2, 3]; default is [])

        Returns
        -------
        MGXS
            A new tally which encapsulates the subset of data requested for the
            nuclide(s) and/or energy group(s) requested in the parameters.

        """

        # Call super class method and null out derived tallies
        slice_xs = super(ScatterMatrixXS, self).get_slice(nuclides, in_groups)
        slice_xs._rxn_rate_tally = None
        slice_xs._xs_tally = None

        # Slice outgoing energy groups if needed
        if len(out_groups) != 0:
            filter_bins = []
            for group in out_groups:
                group_bounds = self.energy_groups.get_group_bounds(group)
                filter_bins.append(group_bounds)
            filter_bins = [tuple(filter_bins)]

            # Slice each of the tallies across energyout groups
            for tally_type, tally in slice_xs.tallies.items():
                if tally.contains_filter('energyout'):
                    tally_slice = tally.get_slice(filters=['energyout'],
                                                  filter_bins=filter_bins)
                    slice_xs.tallies[tally_type] = tally_slice

        slice_xs.sparse = self.sparse
        return slice_xs

    def get_xs(self, in_groups='all', out_groups='all',
               subdomains='all', nuclides='all', xs_type='macro',
               order_groups='increasing', value='mean'):
        """Returns an array of multi-group cross sections.

        This method constructs a 2D NumPy array for the requested scattering
        matrix data data for one or more energy groups and subdomains.

        Parameters
        ----------
        in_groups : Iterable of Integral or 'all'
            Incoming energy groups of interest. Defaults to 'all'.
        out_groups : Iterable of Integral or 'all'
            Outgoing energy groups of interest. Defaults to 'all'.
        subdomains : Iterable of Integral or 'all'
            Subdomain IDs of interest. Defaults to 'all'.
        nuclides : Iterable of str or 'all' or 'sum'
            A list of nuclide name strings (e.g., ['U-235', 'U-238']). The
            special string 'all' will return the cross sections for all nuclides
            in the spatial domain. The special string 'sum' will return the
            cross section summed over all nuclides. Defaults to 'all'.
        xs_type: {'macro', 'micro'}
            Return the macro or micro cross section in units of cm^-1 or barns.
            Defaults to 'macro'.
        order_groups: {'increasing', 'decreasing'}
            Return the cross section indexed according to increasing or
            decreasing energy groups (decreasing or increasing energies).
            Defaults to 'increasing'.
        value : str
            A string for the type of value to return - 'mean', 'std_dev', or
            'rel_err' are accepted. Defaults to the empty string.

        Returns
        -------
        ndarray
            A NumPy array of the multi-group cross section indexed in the order
            each group and subdomain is listed in the parameters.

        Raises
        ------
        ValueError
            When this method is called before the multi-group cross section is
            computed from tally data.

        """

        cv.check_value('value', value, ['mean', 'std_dev', 'rel_err'])
        cv.check_value('xs_type', xs_type, ['macro', 'micro'])

        filters = []
        filter_bins = []

        # Construct a collection of the domain filter bins
        if not isinstance(subdomains, basestring):
            cv.check_iterable_type('subdomains', subdomains, Integral, max_depth=2)
            for subdomain in subdomains:
                filters.append(self.domain_type)
                filter_bins.append((subdomain,))

        # Construct list of energy group bounds tuples for all requested groups
        if not isinstance(in_groups, basestring):
            cv.check_iterable_type('groups', in_groups, Integral)
            for group in in_groups:
                filters.append('energy')
                filter_bins.append((self.energy_groups.get_group_bounds(group),))

        # Construct list of energy group bounds tuples for all requested groups
        if not isinstance(out_groups, basestring):
            cv.check_iterable_type('groups', out_groups, Integral)
            for group in out_groups:
                filters.append('energyout')
                filter_bins.append((self.energy_groups.get_group_bounds(group),))

        # Construct a collection of the nuclides to retrieve from the xs tally
        if self.by_nuclide:
            if nuclides == 'all' or nuclides == 'sum' or nuclides == ['sum']:
                query_nuclides = self.get_all_nuclides()
            else:
                query_nuclides = nuclides
        else:
            query_nuclides = ['total']

        # Use tally summation if user requested the sum for all nuclides
        if nuclides == 'sum' or nuclides == ['sum']:
            xs_tally = self.xs_tally.summation(nuclides=query_nuclides)
            xs = xs_tally.get_values(filters=filters,
                                     filter_bins=filter_bins, value=value)
        else:
            xs = self.xs_tally.get_values(filters=filters,
                                          filter_bins=filter_bins,
                                          nuclides=query_nuclides, value=value)

        xs = np.nan_to_num(xs)

        # Divide by atom number densities for microscopic cross sections
        if xs_type == 'micro':
            if self.by_nuclide:
                densities = self.get_nuclide_densities(nuclides)
            else:
                densities = self.get_nuclide_densities('sum')
            if value == 'mean' or value == 'std_dev':
                xs /= densities[np.newaxis, :, np.newaxis]

        # Reverse data if user requested increasing energy groups since
        # tally data is stored in order of increasing energies
        if order_groups == 'increasing':
            if in_groups == 'all':
                num_in_groups = self.num_groups
            else:
                num_in_groups = len(in_groups)
            if out_groups == 'all':
                num_out_groups = self.num_groups
            else:
                num_out_groups = len(out_groups)

            # Reshape tally data array with separate axes for domain and energy
            num_subdomains = int(xs.shape[0] / (num_in_groups * num_out_groups))
            new_shape = (num_subdomains, num_in_groups, num_out_groups)
            new_shape += xs.shape[1:]
            xs = np.reshape(xs, new_shape)

            # Reverse energies to align with increasing energy groups
            xs = xs[:, ::-1, ::-1, :]

            # Eliminate trivial dimensions
            xs = np.squeeze(xs)
            xs = np.atleast_2d(xs)

        return xs

    def print_xs(self, subdomains='all', nuclides='all', xs_type='macro'):
        """Prints a string representation for the multi-group cross section.

        Parameters
        ----------
        subdomains : Iterable of Integral or 'all'
            The subdomain IDs of the cross sections to include in the report.
            Defaults to 'all'.
        nuclides : Iterable of str or 'all' or 'sum'
            The nuclides of the cross-sections to include in the report. This
            may be a list of nuclide name strings (e.g., ['U-235', 'U-238']).
            The special string 'all' will report the cross sections for all
            nuclides in the spatial domain. The special string 'sum' will report
            the cross sections summed over all nuclides. Defaults to 'all'.
        xs_type: {'macro', 'micro'}
            Return the macro or micro cross section in units of cm^-1 or barns.
            Defaults to 'macro'.

        """

        # Construct a collection of the subdomains to report
        if not isinstance(subdomains, basestring):
            cv.check_iterable_type('subdomains', subdomains, Integral)
        elif self.domain_type == 'distribcell':
            subdomains = np.arange(self.num_subdomains, dtype=np.int)
        else:
            subdomains = [self.domain.id]

        # Construct a collection of the nuclides to report
        if self.by_nuclide:
            if nuclides == 'all':
                nuclides = self.get_all_nuclides()
            if nuclides == 'sum':
                nuclides = ['sum']
            else:
                cv.check_iterable_type('nuclides', nuclides, basestring)
        else:
            nuclides = ['sum']

        cv.check_value('xs_type', xs_type, ['macro', 'micro'])

        # Build header for string with type and domain info
        string = 'Multi-Group XS\n'
        string += '{0: <16}=\t{1}\n'.format('\tReaction Type', self.rxn_type)
        string += '{0: <16}=\t{1}\n'.format('\tDomain Type', self.domain_type)
        string += '{0: <16}=\t{1}\n'.format('\tDomain ID', self.domain.id)

        # If cross section data has not been computed, only print string header
        if self.tallies is None:
            print(string)
            return

        string += '{0: <16}\n'.format('\tEnergy Groups:')
        template = '{0: <12}Group {1} [{2: <10} - {3: <10}MeV]\n'

        # Loop over energy groups ranges
        for group in range(1, self.num_groups+1):
            bounds = self.energy_groups.get_group_bounds(group)
            string += template.format('', group, bounds[0], bounds[1])

        # Loop over all subdomains
        for subdomain in subdomains:

            if self.domain_type == 'distribcell':
                string += \
                    '{0: <16}=\t{1}\n'.format('\tSubdomain', subdomain)

            # Loop over all Nuclides
            for nuclide in nuclides:

                # Build header for nuclide type
                if xs_type != 'sum':
                    string += '{0: <16}=\t{1}\n'.format('\tNuclide', nuclide)

                # Build header for cross section type
                if xs_type == 'macro':
                    string += '{0: <16}\n'.format('\tCross Sections [cm^-1]:')
                else:
                    string += '{0: <16}\n'.format('\tCross Sections [barns]:')

                template = '{0: <12}Group {1} -> Group {2}:\t\t'

                # Loop over incoming/outgoing energy groups ranges
                for in_group in range(1, self.num_groups+1):
                    for out_group in range(1, self.num_groups+1):
                        string += template.format('', in_group, out_group)
                        average = \
                            self.get_xs([in_group], [out_group],
                                        [subdomain], [nuclide],
                                        xs_type=xs_type, value='mean')
                        rel_err = \
                            self.get_xs([in_group], [out_group],
                                        [subdomain], [nuclide],
                                        xs_type=xs_type, value='rel_err')
                        average = average.flatten()[0]
                        rel_err = rel_err.flatten()[0] * 100.
                        string += '{:1.2e} +/- {:1.2e}%'.format(average, rel_err)
                        string += '\n'
                    string += '\n'
                string += '\n'
            string += '\n'

        print(string)


class NuScatterMatrixXS(ScatterMatrixXS):
    """A scattering production matrix multi-group cross section."""

    def __init__(self, domain=None, domain_type=None,
                 groups=None, by_nuclide=False, name=''):
        super(NuScatterMatrixXS, self).__init__(domain, domain_type,
                                                groups, by_nuclide, name)
        self._rxn_type = 'nu-scatter matrix'

    @property
    def tallies(self):
        """Construct the OpenMC tallies needed to compute this cross section.

        This method constructs three analog tallies to compute the 'flux',
        'nu-scatter' and 'scatter-P1' reaction rates in the spatial domain and
        energy groups of interest.

        """

        # Instantiate tallies if they do not exist
        if self._tallies is None:

            # Create the non-domain specific Filters for the Tallies
            group_edges = self.energy_groups.group_edges
            energy = openmc.Filter('energy', group_edges)
            energyout = openmc.Filter('energyout', group_edges)

            # Create a list of scores for each Tally to be created
            if self.correction == 'P0':
                scores = ['flux', 'nu-scatter', 'scatter-P1']
                estimator = 'analog'
                keys = ['flux', 'scatter', 'scatter-P1']
                filters = [[energy], [energy, energyout], [energyout]]
            else:
                scores = ['flux', 'nu-scatter']
                estimator = 'analog'
                keys = ['flux', 'scatter']
                filters = [[energy], [energy, energyout]]

            # Intialize the Tallies
            self._create_tallies(scores, filters, keys, estimator)

        return self._tallies

class Chi(MGXS):
    """The fission spectrum."""

    def __init__(self, domain=None, domain_type=None,
                 groups=None, by_nuclide=False, name=''):
        super(Chi, self).__init__(domain, domain_type, groups, by_nuclide, name)
        self._rxn_type = 'chi'

    @property
    def tallies(self):
        """Construct the OpenMC tallies needed to compute this cross section.

        This method constructs two analog tallies to compute 'nu-fission'
        reaction rates with 'energy' and 'energyout' filters in the spatial
        domain and energy groups of interest.

        """

        # Instantiate tallies if they do not exist
        if self._tallies is None:

            # Create a list of scores for each Tally to be created
            scores = ['nu-fission', 'nu-fission']
            estimator = 'analog'
            keys = ['nu-fission-in', 'nu-fission-out']

            # Create the non-domain specific Filters for the Tallies
            group_edges = self.energy_groups.group_edges
            energyout = openmc.Filter('energyout', group_edges)
            energyin = openmc.Filter('energy', [group_edges[0], group_edges[-1]])
            filters = [[energyin], [energyout]]

            # Intialize the Tallies
            self._create_tallies(scores, filters, keys, estimator)

        return self._tallies

    @property
    def rxn_rate_tally(self):
        if self._rxn_rate_tally is None:
            self._rxn_rate_tally = self.tallies['nu-fission-out']
            self._rxn_rate_tally.sparse = self.sparse
        return self._rxn_rate_tally

    @property
    def xs_tally(self):
        """Computes chi fission spectrum using OpenMC tally arithmetic."""

        if self._xs_tally is None:
            nu_fission_in = self.tallies['nu-fission-in']

            # Remove coarse energy filter to keep it out of tally arithmetic
            energy_filter = nu_fission_in.find_filter('energy')
            nu_fission_in.remove_filter(energy_filter)

            # Compute chi
            self._xs_tally = self.rxn_rate_tally / nu_fission_in
            super(Chi, self)._compute_xs()

            # Add the coarse energy filter back to the nu-fission tally
            nu_fission_in.add_filter(energy_filter)

        return self._xs_tally

    def get_slice(self, nuclides=[], groups=[]):
        """Build a sliced Chi for the specified nuclides and energy groups.

        This method constructs a new MGXS to encapsulate a subset of the data
        represented by this MGXS. The subset of data to include in the tally
        slice is determined by the nuclides and energy groups specified in
        the input parameters.

        Parameters
        ----------
        nuclides : list of str
            A list of nuclide name strings
            (e.g., ['U-235', 'U-238']; default is [])
        groups : list of Integral
            A list of energy group indices starting at 1 for the high energies
            (e.g., [1, 2, 3]; default is [])

        Returns
        -------
        MGXS
            A new tally which encapsulates the subset of data requested for the
            nuclide(s) and/or energy group(s) requested in the parameters.

        """

        # Temporarily remove energy filter from nu-fission-in since its
        # group structure will work in super MGXS.get_slice(...) method
        nu_fission_in = self.tallies['nu-fission-in']
        energy_filter = nu_fission_in.find_filter('energy')
        nu_fission_in.remove_filter(energy_filter)

        # Call super class method and null out derived tallies
        slice_xs = super(Chi, self).get_slice(nuclides, groups)
        slice_xs._rxn_rate_tally = None
        slice_xs._xs_tally = None

        # Slice energy groups if needed
        if len(groups) != 0:
            filter_bins = []
            for group in groups:
                group_bounds = self.energy_groups.get_group_bounds(group)
                filter_bins.append(group_bounds)
            filter_bins = [tuple(filter_bins)]

            # Slice nu-fission-out tally along energyout filter
            nu_fission_out = slice_xs.tallies['nu-fission-out']
            tally_slice = nu_fission_out.get_slice(filters=['energyout'],
                                                   filter_bins=filter_bins)
            slice_xs._tallies['nu-fission-out'] = tally_slice

        # Add energy filter back to nu-fission-in tallies
        self.tallies['nu-fission-in'].add_filter(energy_filter)
        slice_xs._tallies['nu-fission-in'].add_filter(energy_filter)

        slice_xs.sparse = self.sparse
        return slice_xs

    def merge(self, other):
        """Merge another Chi with this one

        If results have been loaded from a statepoint, then Chi are only
        mergeable along one and only one of energy groups or nuclides.

        Parameters
        ----------
        other : MGXS
            MGXS to merge with this one

        Returns
        -------
        merged_mgxs : MGXS
            Merged MGXS
        """

        if not self.can_merge(other):
            raise ValueError('Unable to merge Chi')

        # Create deep copy of tally to return as merged tally
        merged_mgxs = copy.deepcopy(self)
        merged_mgxs._derived = True
        merged_mgxs._rxn_rate_tally = None
        merged_mgxs._xs_tally = None

        # Merge energy groups
        if self.energy_groups != other.energy_groups:
            merged_groups = self.energy_groups.merge(other.energy_groups)
            merged_mgxs.energy_groups = merged_groups

        # Merge nuclides
        if self.nuclides != other.nuclides:

            # The nuclides must be mutually exclusive
            for nuclide in self.nuclides:
                if nuclide in other.nuclides:
                    msg = 'Unable to merge Chi with shared nuclides'
                    raise ValueError(msg)

            # Concatenate lists of nuclides for the merged MGXS
            merged_mgxs.nuclides = self.nuclides + other.nuclides

        # Merge tallies
        for tally_key in self.tallies:
            merged_tally = self.tallies[tally_key].merge(other.tallies[tally_key])
            merged_mgxs.tallies[tally_key] = merged_tally

        return merged_mgxs

    def get_xs(self, groups='all', subdomains='all', nuclides='all',
               xs_type='macro', order_groups='increasing', value='mean'):
        """Returns an array of the fission spectrum.

        This method constructs a 2D NumPy array for the requested multi-group
        cross section data data for one or more energy groups and subdomains.

        Parameters
        ----------
        groups : Iterable of Integral or 'all'
            Energy groups of interest. Defaults to 'all'.
        subdomains : Iterable of Integral or 'all'
            Subdomain IDs of interest. Defaults to 'all'.
        nuclides : Iterable of str or 'all' or 'sum'
            A list of nuclide name strings (e.g., ['U-235', 'U-238']). The
            special string 'all' will return the cross sections for all nuclides
            in the spatial domain. The special string 'sum' will return the
            cross section summed over all nuclides. Defaults to 'all'.
        xs_type: {'macro', 'micro'}
            This parameter is not relevant for chi but is included here to
            mirror the parent MGXS.get_xs(...) class method
        order_groups: {'increasing', 'decreasing'}
            Return the cross section indexed according to increasing or
            decreasing energy groups (decreasing or increasing energies).
            Defaults to 'increasing'.
        value : str
            A string for the type of value to return - 'mean', 'std_dev', or
            'rel_err' are accepted. Defaults to 'mean'.

        Returns
        -------
        ndarray
            A NumPy array of the multi-group cross section indexed in the order
            each group, subdomain and nuclide is listed in the parameters.

        Raises
        ------
        ValueError
            When this method is called before the multi-group cross section is
            computed from tally data.

        """

        cv.check_value('value', value, ['mean', 'std_dev', 'rel_err'])
        cv.check_value('xs_type', xs_type, ['macro', 'micro'])

        filters = []
        filter_bins = []

        # Construct a collection of the domain filter bins
        if not isinstance(subdomains, basestring):
            cv.check_iterable_type('subdomains', subdomains, Integral, max_depth=2)
            for subdomain in subdomains:
                filters.append(self.domain_type)
                filter_bins.append((subdomain,))

        # Construct list of energy group bounds tuples for all requested groups
        if not isinstance(groups, basestring):
            cv.check_iterable_type('groups', groups, Integral)
            for group in groups:
                filters.append('energyout')
                filter_bins.append((self.energy_groups.get_group_bounds(group),))

        # If chi was computed for each nuclide in the domain
        if self.by_nuclide:

            # Get the sum as the fission source weighted average chi for all
            # nuclides in the domain
            if nuclides == 'sum' or nuclides == ['sum']:

                # Retrieve the fission production tallies
                nu_fission_in = self.tallies['nu-fission-in']
                nu_fission_out = self.tallies['nu-fission-out']

                # Sum out all nuclides
                nuclides = self.get_all_nuclides()
                nu_fission_in = nu_fission_in.summation(nuclides=nuclides)
                nu_fission_out = nu_fission_out.summation(nuclides=nuclides)

                # Remove coarse energy filter to keep it out of tally arithmetic
                energy_filter = nu_fission_in.find_filter('energy')
                nu_fission_in.remove_filter(energy_filter)

                # Compute chi and store it as the xs_tally attribute so we can
                # use the generic get_xs(...) method
                xs_tally = nu_fission_out / nu_fission_in

                # Add the coarse energy filter back to the nu-fission tally
                nu_fission_in.add_filter(energy_filter)

                xs = xs_tally.get_values(filters=filters,
                                         filter_bins=filter_bins, value=value)

            # Get chi for all nuclides in the domain
            elif nuclides == 'all':
                nuclides = self.get_all_nuclides()
                xs = self.xs_tally.get_values(filters=filters,
                                              filter_bins=filter_bins,
                                              nuclides=nuclides, value=value)

            # Get chi for user-specified nuclides in the domain
            else:
                cv.check_iterable_type('nuclides', nuclides, basestring)
                xs = self.xs_tally.get_values(filters=filters,
                                              filter_bins=filter_bins,
                                              nuclides=nuclides, value=value)

        # If chi was computed as an average of nuclides in the domain
        else:
            xs = self.xs_tally.get_values(filters=filters,
                                          filter_bins=filter_bins, value=value)

        # Reverse data if user requested increasing energy groups since
        # tally data is stored in order of increasing energies
        if order_groups == 'increasing':

            # Reshape tally data array with separate axes for domain and energy
            if groups == 'all':
                num_groups = self.num_groups
            else:
                num_groups = len(groups)
            num_subdomains = int(xs.shape[0] / num_groups)
            new_shape = (num_subdomains, num_groups) + xs.shape[1:]
            xs = np.reshape(xs, new_shape)

            # Reverse energies to align with increasing energy groups
            xs = xs[:, ::-1, :]

            # Eliminate trivial dimensions
            xs = np.squeeze(xs)
            xs = np.atleast_1d(xs)

        xs = np.nan_to_num(xs)
        return xs

    def get_pandas_dataframe(self, groups='all', nuclides='all',
                             xs_type='macro', summary=None):
        """Build a Pandas DataFrame for the MGXS data.

        This method leverages the Tally.get_pandas_dataframe(...) method, but
        renames the columns with terminology appropriate for cross section data.

        Parameters
        ----------
        groups : Iterable of Integral or 'all'
            Energy groups of interest. Defaults to 'all'.
        nuclides : Iterable of str or 'all' or 'sum'
            The nuclides of the cross-sections to include in the dataframe. This
            may be a list of nuclide name strings (e.g., ['U-235', 'U-238']).
            The special string 'all' will include the cross sections for all
            nuclides in the spatial domain. The special string 'sum' will
            include the cross sections summed over all nuclides. Defaults to
            'all'.
        xs_type: {'macro', 'micro'}
            Return macro or micro cross section in units of cm^-1 or barns.
            Defaults to 'macro'.
        summary : None or Summary
            An optional Summary object to be used to construct columns for
            distribcell tally filters (default is None). The geometric
            information in the Summary object is embedded into a multi-index
            column with a geometric "path" to each distribcell intance.
            NOTE: This option requires the OpenCG Python package.

        Returns
        -------
        pandas.DataFrame
            A Pandas DataFrame for the cross section data.

        Raises
        ------
        ValueError
            When this method is called before the multi-group cross section is
            computed from tally data.

        """

        # Build the dataframe using the parent class method
        df = super(Chi, self).get_pandas_dataframe(groups, nuclides,
                                                   xs_type, summary)

        # If user requested micro cross sections, multiply by the atom
        # densities to cancel out division made by the parent class method
        if xs_type == 'micro':
            if self.by_nuclide:
                densities = self.get_nuclide_densities(nuclides)
            else:
                densities = self.get_nuclide_densities('sum')
            tile_factor = df.shape[0] / len(densities)
            df['mean'] *= np.tile(densities, tile_factor)
            df['std. dev.'] *= np.tile(densities, tile_factor)

        return df
