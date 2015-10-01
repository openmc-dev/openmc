from collections import Iterable
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


# Supported domain types
DOMAIN_TYPES = ['cell',
                'distribcell',
                'universe',
                'material']

# Supported domain objects
DOMAINS = [openmc.Cell,
           openmc.Universe,
           openmc.Material]


class MultiGroupXS(object):
    """A multi-group cross section for some energy group structure within
    some spatial domain.

    This class can be used for both OpenMC input generation and tally data
    post-processing to compute spatially-homogenized and energy-integrated
    multi-group cross sections for deterministic neutronics calculations.

    Parameters
    ----------
    domain : Material or Cell or Universe
        The domain for spatial homogenization
    domain_type : {'material', 'cell', 'distribcell', 'universe'}
        The domain type for spatial homogenization
    energy_groups : EnergyGroups
        The energy group structure for energy condensation
    by_nuclide : bool
        If true, computes multi-group cross sections for each nuclide in domain
    name : str, optional
        Name of the multi-group cross section. Used as a label to identify
        tallies in OpenMC tallies.xml file.

    Attributes
    ----------
    name : str, optional
        Name of the multi-group cross section
    rxn_type : str
        Reaction type (e.g., 'total', 'nu-fission', etc.)
    by_nuclide : bool
        If true, computes multi-group cross sections for each nuclide in domain
    domain : Material or Cell or Universe
        Domain for spatial homogenization
    domain_type : {'material', 'cell', 'distribcell', 'universe'}
        Domain type for spatial homogenization
    energy_groups : EnergyGroups
        Energy group structure for energy condensation
    num_groups : Integral
        Number of energy groups
    tallies : dict
        OpenMC tallies needed to compute the multi-group cross section
    xs_tally : Tally
        Derived tally for the multi-group cross section. This attribute
        is None unless the multi-group cross section has been computed.

    """

    # This is an abstract class which cannot be instantiated
    __metaclass__ = abc.ABCMeta

    def __init__(self, domain=None, domain_type=None,
                 energy_groups=None, by_nuclide=False, name=''):

        self._name = ''
        self._rxn_type = None
        self._by_nuclide = None
        self._domain = None
        self._domain_type = None
        self._energy_groups = None
        self._num_groups = None
        self._tallies = dict()
        self._xs_tally = None

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
            clone._domain = self.domain
            clone._domain_type = self.domain_type
            clone._energy_groups = copy.deepcopy(self.energy_groups, memo)
            clone._num_groups = self.num_groups
            clone._xs_tally = copy.deepcopy(self.xs_tally, memo)

            clone._tallies = dict()
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
    def num_groups(self):
        return self._num_groups

    @property
    def tallies(self):
        return self._tallies

    @property
    def xs_tally(self):
        return self._xs_tally

    @property
    def num_subdomains(self):
        tally = self.tallies.values()[0]
        domain_filter = tally.find_filter(self.domain_type)
        return domain_filter.num_bins

    @name.setter
    def name(self, name):
        cv.check_type('name', name, basestring)
        self._name = name

    @by_nuclide.setter
    def by_nuclide(self, by_nuclide):
        cv.check_type('by_nuclide', by_nuclide, bool)
        self._by_nuclide = by_nuclide

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

    @domain.setter
    def domain(self, domain):
        cv.check_type('domain', domain, tuple(DOMAINS))
        self._domain = domain

    @domain_type.setter
    def domain_type(self, domain_type):
        cv.check_value('domain type', domain_type, tuple(DOMAIN_TYPES))
        self._domain_type = domain_type

    @energy_groups.setter
    def energy_groups(self, energy_groups):
        cv.check_type('energy groups', energy_groups, openmc.mgxs.EnergyGroups)
        self._energy_groups = energy_groups
        self._num_groups = energy_groups.num_groups

    def get_all_nuclides(self):
        """Get all nuclides in the cross section's spatial domain.

        Returns
        -------
        nuclides : list of str
            A list of the string names for each nuclide in the problem domain
            (e.g., ['U-235', 'U-238', 'O-16'])

        Raises
        ------
        ValueError
            When this method is called before the spatial domain has been set.

        """

        if self.domain is None:
            raise ValueError('Unable to get all nuclides without a domain')

        nuclides = self.domain.get_all_nuclides()
        return nuclides.keys()

    def get_nuclide_density(self, nuclide):
        """Get the atomic number density in units of atoms/b-cm for a nuclide
        in the cross section's spatial domain.

        Paramters
        ---------
        nuclide : str
            A nuclide name string (e.g., 'U-235')

        Returns
        -------
        density : Real
            The atomic number density (atom/b-cm) for the nuclide of interest

        Raises
        ------
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

        Paramters
        ---------
        nuclides : Iterable of str or 'all' or 'sum'
            A list of nuclide name strings (e.g., ['U-235', 'U-238']). The
            special string 'all' will return the atom densities for all nuclides
            in the spatial domain. The special string 'sum' will return the atom
            density summed across all nuclides in the spatial domain.

        Returns
        -------
        densities : ndarray of float
            An array of the atomic number densities (atom/b-cm) for each of the
            nuclides in the problem domain

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

        # Sum the atomic number densities for all nuclides
        elif nuclides == 'all':
            nuclides = self.get_all_nuclides()
            densities = np.zeros(self.num_nuclides, dtype=np.float)
            for i, nuclide in enumerate(nuclides):
                densities[i] += self.get_nuclide_density(nuclide)

        # Store each nuclide's atomic number density in an array
        else:
            densities = np.zeros(len(nuclides), dtype=np.float)
            for i, nuclide in enumerate(nuclides):
                densities[i] = self.get_nuclide_density(nuclide)

        return densities

    @abc.abstractmethod
    def create_tallies(self, scores, all_filters, keys, estimator):
        """Instantiates tallies needed to compute the multi-group cross section.

        This is a helper method for MultiGroupXS subclasses to create tallies
        for input file generation. The tallies are stored in the tallies dict.

        Parameters
        ----------
        scores : Iterable of str
            Scores for each tally

        filters : Iterable of tuple of Filter
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

        # Create a domain Filter object
        domain_filter = openmc.Filter(self.domain_type, self.domain.id)
        domain_filter.num_bins = 1

        for score, key, filters in zip(scores, keys, all_filters):
            self.tallies[key] = openmc.Tally(name=self.name)
            self.tallies[key].add_score(score)
            self.tallies[key].estimator = estimator
            self.tallies[key].add_filter(domain_filter)

            # Add all non-domain specific Filters (e.g., 'energy') to the Tally
            for filter in filters:
                self.tallies[key].add_filter(filter)

            # If this is a by nuclide cross-section, add all nuclides to Tally
            if self.by_nuclide and score != 'flux':
                all_nuclides = self.domain.get_all_nuclides()
                for nuclide in all_nuclides:
                    self.tallies[key].add_nuclide(nuclide)
            else:
                self.tallies[key].add_nuclide('total')

    @abc.abstractmethod
    def compute_xs(self):
        """Performs generic cleanup after a subclass' uses tally arithmetic to
        compute a multi-group cross section as a derived tally."""

        # If computing xs for each nuclide, replace CrossNuclides with originals
        if self.by_nuclide:
            self.xs_tally._nuclides = []
            nuclides = self.domain.get_all_nuclides()
            for nuclide in nuclides:
                self.xs_tally.add_nuclide(openmc.Nuclide(nuclide))

        # Remove NaNs which may have resulted from divide-by-zero operations
        self._xs_tally._mean = np.nan_to_num(self.xs_tally.mean)
        self._xs_tally._std_dev = np.nan_to_num(self.xs_tally.std_dev)

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

        """

        cv.check_type('statepoint', statepoint, openmc.statepoint.StatePoint)

        if not statepoint.with_summary:
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
            msg = 'Unable to load data from a statepoint for domain type {} ' \
                  'which is not yet supported'.format(self.domain_type)
            raise ValueError(msg)

        # Create Tallies to search for in StatePoint
        self.create_tallies()

        if self.domain_type == 'distribcell':
            filters = []
            filter_bins = []
        else:
            filters = [self.domain_type]
            filter_bins = [(self.domain.id,)]

        # Find, slice and store Tallies from StatePoint
        # The tally slicing is needed if tally merging was used
        for tally_type, tally in self.tallies.items():
            sp_tally = statepoint.get_tally(tally.scores, tally.filters,
                                            tally.nuclides,
                                            estimator=tally.estimator)
            sp_tally = sp_tally.get_slice(tally.scores, filters,
                                          filter_bins, tally.nuclides)
            self.tallies[tally_type] = sp_tally

    def get_xs(self, groups='all', subdomains='all', nuclides='all',
               xs_type='macro', order_groups='increasing', value='mean'):
        """Returns an array of multi-group cross sections.

        This method constructs a 2D NumPy array for the requested multi-group
        cross section data data for one or more energy groups and subdomains.

        Parameters
        ----------
        groups : Iterable of Integral or 'all'
            Energy groups of interest

        subdomains : Iterable of Integral or 'all'
            Subdomain IDs of interest

        nuclides : Iterable of str or 'all' or 'sum'
            A list of nuclide name strings (e.g., ['U-235', 'U-238']). The
            special string 'all' (default) will return the cross sections for
            all nuclides in the spatial domain. The special string 'sum' will
            return the cross section summed over all nuclides.

        xs_type: {'macro' or 'micro'}
            Return the macro or micro cross section in units of cm^-1 or barns

        order_groups: {'increasing', 'decreasing'}
            Return the cross section indexed according to increasing (default)
            or decreasing energy groups (decreasing or increasing energies)

        value : str
            A string for the type of value to return - 'mean' (default),
            'std_dev' or 'rel_err' are accepted

        Returns
        -------
        xs : ndarray
            A NumPy array of the multi-group cross section indexed in the order
            each group, subdomain and nuclide is listed in the parameters.

        Raises
        ------
        ValueError
            When this method is called before the multi-group cross section is
            computed from tally data.

        """

        if self.xs_tally is None:
            msg = 'Unable to get cross section since it has not been computed'
            raise ValueError(msg)

        cv.check_value('value', value, ['mean', 'std_dev', 'rel_err'])
        cv.check_value('xs_type', xs_type, ['macro', 'micro'])

        filters = []
        filter_bins = []

        # Construct a collection of the domain filter bins
        if subdomains != 'all':
            cv.check_iterable_type('subdomains', subdomains, Integral)
            for subdomain in subdomains:
                filters.append(self.domain_type)
                filter_bins.append((subdomain,))

        # Construct list of energy group bounds tuples for all requested groups
        if groups != 'all':
            cv.check_iterable_type('groups', groups, Integral)
            for group in groups:
                filters.append('energy')
                filter_bins.append((self.energy_groups.get_group_bounds(group),))

        # Construct a collection of the nuclides to retrieve from the xs tally
        # NOTE: We must not override the "nuclides" parameter since it is used
        # to retrieve atomic number densities for micro xs
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
            # Reshape tally data array with separate axes for domain and energy
            if groups == 'all':
                num_groups = self.num_groups
            else:
                num_groups = len(groups)
            num_subdomains = xs.shape[0] / num_groups
            new_shape = (num_subdomains, num_groups) + xs.shape[1:]
            xs = np.reshape(xs, new_shape)

            # Reverse energies to align with increasing energy groups
            xs = xs[:, ::-1, :]

            # Reshape array to original axes (filters, nuclides, scores)
            new_shape = (num_subdomains * num_groups,) + xs.shape[2:]
            xs = np.reshape(xs, new_shape)

        return xs

    def get_condensed_xs(self, coarse_groups):
        """Construct an energy-condensed version of this cross section.

        Parameters
        ----------
        coarse_groups : openmc.mgxs.EnergyGroups
            The coarse energy group structure of interest

        Returns
        -------
        MultiGroupXS
            A new MultiGroupXS condensed to the group structure of interest

        """

        if self.xs_tally is None:
            msg = 'Unable to get a condensed coarse group cross section ' \
                  'since the fine group cross section has not been computed'
            raise ValueError(msg)

        cv.check_type('coarse_groups', coarse_groups, EnergyGroups)
        cv.check_less_than('coarse groups', coarse_groups.num_groups,
                           self.num_groups, equality=True)
        cv.check_value('upper coarse energy', coarse_groups.group_edges[-1],
                       [self.energy_groups.group_edges[-1]])
        cv.check_value('lower coarse energy', coarse_groups.group_edges[0],
                       [self.energy_groups.group_edges[0]])

        # Clone this MultiGroupXS to initialize the condensed version
        condensed_xs = copy.deepcopy(self)
        condensed_xs.energy_groups = coarse_groups

        # Build indices to sum up over
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
            for i, filter in enumerate(tally.filters):
                if 'energy' in filter.type and all(filter.bins == fine_edges):
                    filter.bins = coarse_groups.group_edges
                    mean = np.add.reduceat(mean, energy_indices, axis=i)
                    std_dev = np.add.reduceat(std_dev**2, energy_indices, axis=i)
                    std_dev = np.sqrt(std_dev)

            # Reshape condensed data arrays with one dimension for all filters
            new_shape = \
                (tally.num_filter_bins, tally.num_nuclides, tally.num_score_bins,)
            mean = np.reshape(mean, new_shape)
            std_dev = np.reshape(std_dev, new_shape)

            # Override tally's data with the new condensed data
            tally._mean = mean
            tally._std_dev = std_dev

        # Compute the energy condensed multi-group cross section
        condensed_xs.compute_xs()
        return condensed_xs

    def get_subdomain_avg_xs(self, subdomains='all'):
        """Construct a subdomain-averaged version of this cross section.

        This is primarily useful for averaging across distribcell instances.
        This routine performs spatial homogenization to compute the scalar
        flux-weighted average cross section across the subdomains.

        Parameters
        ----------
        subdomains : Iterable of Integral or 'all'
            The subdomain IDs to average across

        Returns
        -------
        MultiGroupXS
            A new MultiGroupXS averaged across the subdomains of interest

        Raises
        ------
        ValueError
            When this method is called before the multi-group cross section is
            computed from tally data.

        """

        if self.xs_tally is None:
            msg = 'Unable to get subdomain-averaged cross section since the ' \
                  'subdomain-distributed cross section has not been computed'
            raise ValueError(msg)

        # Construct a collection of the subdomain filter bins to average across
        if subdomains != 'all':
            cv.check_iterable_type('subdomains', subdomains, Integral)
        elif self.domain_type == 'distribcell':
            subdomains = np.arange(self.num_subdomains)
        else:
            subdomains = [0]

        # Clone this MultiGroupXS to initialize the subdomain-averaged version
        avg_xs = copy.deepcopy(self)
        avg_xs.domain_type = 'cell'

        # Average each of the tallies across subdomains
        for tally_type, tally in avg_xs.tallies.items():

            # Make condensed tally derived and null out sum, sum_sq
            tally._derived = True
            tally._sum = None
            tally._sum_sq = None

            # Get tally data arrays reshaped with one dimension per filter
            mean = tally.get_reshaped_data(value='mean')
            std_dev = tally.get_reshaped_data(value='std_dev')

            # Get the mean of the mean, std. dev. across requested subdomains
            mean = np.mean(mean[subdomains, ...], axis=0)
            std_dev = np.mean(std_dev[subdomains, ...]**2, axis=0)
            std_dev = np.sqrt(std_dev)

            # If domain is distribcell, make subdomain-averaged a 'cell' domain
            domain_filter = tally.find_filter(self._domain_type)
            if domain_filter.type == 'distribcell':
                domain_filter.type = 'cell'
                domain_filter.num_bins = 1

            # Reshape averaged data arrays with one dimension for all filters
            new_shape = \
                (tally.num_filter_bins, tally.num_nuclides, tally.num_score_bins,)
            mean = np.reshape(mean, new_shape)
            std_dev = np.reshape(std_dev, new_shape)

            # Override tally's data with the new condensed data
            tally._mean = mean
            tally._std_dev = std_dev

        # Compute the subdomain-averaged multi-group cross section
        avg_xs.compute_xs()

        return avg_xs

    def print_xs(self, subdomains='all', nuclides='all', xs_type='macro'):
        """Prints a string representation for the multi-group cross section.

        Parameters
        ----------
        subdomains : Iterable of Integral or 'all'
            The subdomain IDs of the cross sections to include in the report

        nuclides : Iterable of str or 'all' or 'sum'
            The nuclides of the cross-sections to include in the report. This
            may be a list of nuclide name strings (e.g., ['U-235', 'U-238']).
            The special string 'all' (default) will report the cross sections
            for all nuclides in the spatial domain. The special string 'sum'
            will report the cross sections summed over all nuclides.

        xs_type: {'macro' or 'micro'}
            Return the macro or micro cross section in units of cm^-1 or barns

        """

        # Construct a collection of the subdomains to report
        if subdomains != 'all':
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
        if self.xs_tally is None:
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
                                          xs_type=xs_type, value='rel_err')*100
                    average = np.nan_to_num(average.flatten())[0]
                    rel_err = np.nan_to_num(rel_err.flatten())[0]
                    string += '{:.2e} +/- {:1.2e}%'.format(average, rel_err)
                    string += '\n'
                string += '\n'
            string += '\n'

        print(string)

    def build_hdf5_store(self, filename='mgxs', directory='mgxs',
                         xs_type='macro', append=True):
        """Export the multi-group cross section data into an HDF5 binary file.

        This routine constructs an HDF5 file which stores the multi-group
        cross section data. The data is be stored in a hierarchy of HDF5 groups
        from the domain type, domain id, subdomain id (for distribcell domains),
        and cross section type. Two datasets for the mean and standard deviation
        are stored for each subddomain entry in the HDF5 file.

        NOTE: This requires the h5py Python package.

        Parameters
        ----------
        filename : str
            Filename for the HDF5 file (default is 'mgxs')

        directory : str
            Directory for the HDF5 file (default is 'mgxs')

        xs_type: {'macro' or 'micro'}
            Store the macro or micro cross section in units of cm^-1 or barns

        append : boolean
            If true, appends to an existing HDF5 file with the same filename
            directory (if one exists)

        Raises
        ------
        ValueError
            When this method is called before the multi-group cross section is
            computed from tally data.
        ImportError
            When h5py is not installed.

        """

        if self.xs_tally is None:
            msg = 'Unable to get build HDF5 store since the ' \
                  'cross section has not been computed'
            raise ValueError(msg)

        # Attempt to import h5py
        try:
            import h5py
        except ImportError:
            msg = 'The h5py Python package must be installed on your system'
            raise ImportError(msg)

        # Make directory if it does not exist
        if not os.path.exists(directory):
            os.makedirs(directory)

        filename = directory + '/' + filename + '.h5'
        filename = filename.replace(' ', '-')

        if append and os.path.isfile(filename):
            xs_results = h5py.File(filename, 'a')
        else:
            xs_results = h5py.File(filename, 'w')

        cv.check_value('xs_type', xs_type, ['macro', 'micro'])

        if self.by_nuclide:
            nuclides = self.domain.get_all_nuclides()
            densities = np.zeros(len(nuclides), dtype=np.float)
            for i, nuclide in enumerate(nuclides):
                densities[i] = nuclides[nuclide][1]
        else:
            nuclides = ['sum']

        # Create an HDF5 group within the file for the domain
        domain_type_group = xs_results.require_group(self.domain_type)
        group_name = '{0} {1}'.format(self.domain_type, self.domain.id)
        domain_group = domain_type_group.require_group(group_name)

        if self.domain_type == 'distribcell':
            subdomains = np.arange(self.num_subdomains, dtype=np.int)
        else:
            subdomains = [self.domain.id]

        # Determine number of digits to pad subdomain group keys
        num_digits = len(str(self.num_subdomains))

        # Create a separate HDF5 group for each subdomain
        for i, subdomain in enumerate(subdomains):

            # Create an HDF5 group for the subdomain
            if self.domain_type == 'distribcell':
                group_name = str(subdomain).zfill(num_digits)
                subdomain_group = domain_group.require_group(group_name)
            else:
                subdomain_group = domain_group

            # Create a separate HDF5 group for the rxn type
            rxn_group = subdomain_group.require_group(self.rxn_type)

            # Create a separate HDF5 group for each nuclide
            for j, nuclide in enumerate(nuclides):

                if nuclide != 'sum':
                    nuclide_group = rxn_group.require_group(nuclide)
                    nuclide_group.require_dataset('density', dtype=np.float64,
                                                  data=[densities[j]], shape=(1,))
                else:
                    nuclide_group = rxn_group

                # Extract the cross section for this subdomain and nuclide
                average = self.get_xs(subdomains=[subdomain], nuclides=[nuclide],
                                      xs_type=xs_type, value='mean')
                std_dev = self.get_xs(subdomains=[subdomain], nuclides=[nuclide],
                                      xs_type=xs_type, value='std_dev')
                average = average.squeeze()
                std_dev = std_dev.squeeze()

                # Add MultiGroupXS results data to the HDF5 group
                nuclide_group.require_dataset('average', dtype=np.float64,
                                              shape=average.shape, data=average)
                nuclide_group.require_dataset('std. dev.', dtype=np.float64,
                                              shape=std_dev.shape, data=std_dev)

        # Close the MultiGroup results HDF5 file
        xs_results.close()

    def export_xs_data(self, filename='mgxs', directory='mgxs',
                       format='csv', groups='all', xs_type='macro'):
        """Export the multi-group cross section data to a file.

        This routine leverages the functionality in the Pandas library to
        export the multi-group cross section data in a variety of output
        file formats for storage and/or post-processing.

        Parameters
        ----------
        filename : str
            Filename for the exported file (default is 'mgxs')

        directory : str
            Directory for the exported file (default is 'mgxs')

        format : {'csv', 'excel', 'pickle', 'latex'}
            The format for the exported data file

        groups : Iterable of Integral or 'all'
            Energy groups of interest

        xs_type: {'macro' or 'micro'}
            Store the macro or micro cross section in units of cm^-1 or barns

        """

        cv.check_type('filename', filename, basestring)
        cv.check_type('directory', directory, basestring)
        cv.check_value('format', format, ['csv', 'excel', 'pickle', 'latex'])
        cv.check_value('xs_type', xs_type, ['macro', 'micro'])

        # Make directory if it does not exist
        if not os.path.exists(directory):
            os.makedirs(directory)

        filename = directory + '/' + filename
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
        """Build a Pandas DataFrame for the MultiGroupXS data.

        This routine leverages the Tally.get_pandas_dataframe(...) routine, but
        renames the columns with terminology appropriate for cross section data.

        Parameters
        ----------
        groups : Iterable of Integral or 'all'
            Energy groups of interest

        nuclides : Iterable of str or 'all' or 'sum'
            The nuclides of the cross-sections to include in the dataframe. This
            may be a list of nuclide name strings (e.g., ['U-235', 'U-238']).
            The special string 'all' (default) will include the cross sections
            for all nuclides in the spatial domain. The special string 'sum'
            will include the cross sections summed over all nuclides.

        xs_type: {'macro' or 'micro'}
            Return macro or micro cross section in units of cm^-1 or barns

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

        if self.xs_tally is None:
            msg = 'Unable to get Pandas DataFrame since the ' \
                  'cross section has not been computed'
            raise ValueError(msg)

        if groups != 'all':
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
        if summary and self.domain_type == 'distribcell':
            df = df.drop('score', level=0, axis=1)
        else:
            df = df.drop('score', axis=1)

        # Rename energy(out) columns
        columns = []
        if 'energy [MeV]' in df:
            df.rename(columns={'energy [MeV]': 'group in'}, inplace=True)
            columns.append('group in')
        if 'energyout [MeV]' in df:
            df.rename(columns={'energyout [MeV]': 'group out'}, inplace=True)
            columns.append('group out')

        # Loop over all energy groups and override the bounds with indices
        template = '({0:.1e} - {1:.1e})'
        bins = self.energy_groups.group_edges
        for column in columns:
            for i in range(self.num_groups):
                group = template.format(bins[i], bins[i+1])
                row_indices = df[column] == group
                df.loc[row_indices, column] = self.num_groups - i

        # Select out those groups the user requested
        if groups != 'all':
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
        df.sort([self.domain_type] + columns, inplace=True)

        return df


class TotalXS(MultiGroupXS):

    def __init__(self, domain=None, domain_type=None,
                 groups=None, by_nuclide=False, name=''):
        super(TotalXS, self).__init__(domain, domain_type, groups, by_nuclide, name)
        self._rxn_type = 'total'

    def create_tallies(self):
        """Construct the OpenMC tallies needed to compute this cross section."""

        # Create a list of scores for each Tally to be created
        scores = ['flux', 'total']
        estimator = 'tracklength'
        keys = scores

        # Create the non-domain specific Filters for the Tallies
        group_edges = self.energy_groups.group_edges
        energy_filter = openmc.Filter('energy', group_edges)
        filters = [[energy_filter], [energy_filter]]

        # Initialize the Tallies
        super(TotalXS, self).create_tallies(scores, filters, keys, estimator)

    def compute_xs(self):
        """Computes the multi-group total cross sections using OpenMC
        tally arithmetic."""

        self._xs_tally = self.tallies['total'] / self.tallies['flux']
        super(TotalXS, self).compute_xs()


class TransportXS(MultiGroupXS):

    def __init__(self, domain=None, domain_type=None,
                 groups=None, by_nuclide=False, name=''):
        super(TransportXS, self).__init__(domain, domain_type, groups, by_nuclide, name)
        self._rxn_type = 'transport'

    def create_tallies(self):
        """Construct the OpenMC tallies needed to compute this cross section."""

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
        super(TransportXS, self).create_tallies(scores, filters, keys, estimator)

    def load_from_statepoint(self, statepoint):
        super(TransportXS, self).load_from_statepoint(statepoint)
        scatter_p1 = self.tallies['scatter-P1']
        self.tallies['scatter-P1'] = scatter_p1.get_slice(scores=['scatter-P1'])
        self.tallies['scatter-P1'].filters[-1].type = 'energy'

    def compute_xs(self):
        """Computes the multi-group transport cross sections using OpenMC
        tally arithmetic."""

        self._xs_tally = self.tallies['total'] - self.tallies['scatter-P1']
        self._xs_tally /= self.tallies['flux']
        super(TransportXS, self).compute_xs()


class AbsorptionXS(MultiGroupXS):

    def __init__(self, domain=None, domain_type=None,
                 groups=None, by_nuclide=False, name=''):
        super(AbsorptionXS, self).__init__(domain, domain_type, groups, by_nuclide, name)
        self._rxn_type = 'absorption'

    def create_tallies(self):
        """Construct the OpenMC tallies needed to compute this cross section."""

        # Create a list of scores for each Tally to be created
        scores = ['flux', 'absorption']
        estimator = 'tracklength'
        keys = scores

        # Create the non-domain specific Filters for the Tallies
        group_edges = self.energy_groups.group_edges
        energy_filter = openmc.Filter('energy', group_edges)
        filters = [[energy_filter], [energy_filter]]

        # Initialize the Tallies
        super(AbsorptionXS, self).create_tallies(scores, filters, keys, estimator)

    def compute_xs(self):
        """Computes the multi-group absorption cross sections using OpenMC
        tally arithmetic."""

        self._xs_tally = self.tallies['absorption'] / self.tallies['flux']
        super(AbsorptionXS, self).compute_xs()


class CaptureXS(MultiGroupXS):

    def __init__(self, domain=None, domain_type=None,
                 groups=None, by_nuclide=False, name=''):
        super(CaptureXS, self).__init__(domain, domain_type, groups, by_nuclide, name)
        self._rxn_type = 'capture'

    def create_tallies(self):
        """Construct the OpenMC tallies needed to compute this cross section."""

        # Create a list of scores for each Tally to be created
        scores = ['flux', 'absorption', 'fission']
        estimator = 'tracklength'
        keys = scores

        # Create the non-domain specific Filters for the Tallies
        group_edges = self.energy_groups.group_edges
        energy_filter = openmc.Filter('energy', group_edges)
        filters = [[energy_filter], [energy_filter], [energy_filter]]

        # Initialize the Tallies
        super(CaptureXS, self).create_tallies(scores, filters, keys, estimator)

    def compute_xs(self):
        """Computes the multi-group capture cross sections using OpenMC
        tally arithmetic."""

        self._xs_tally = self.tallies['absorption'] - self.tallies['fission']
        self._xs_tally /= self.tallies['flux']
        super(CaptureXS, self).compute_xs()


class FissionXS(MultiGroupXS):

    def __init__(self, domain=None, domain_type=None,
                 groups=None, by_nuclide=False, name=''):
        super(FissionXS, self).__init__(domain, domain_type, groups, by_nuclide, name)
        self._rxn_type = 'fission'

    def create_tallies(self):
        """Construct the OpenMC tallies needed to compute this cross section."""

        # Create a list of scores for each Tally to be created
        scores = ['flux', 'fission']
        estimator = 'tracklength'
        keys = scores

        # Create the non-domain specific Filters for the Tallies
        group_edges = self.energy_groups.group_edges
        energy_filter = openmc.Filter('energy', group_edges)
        filters = [[energy_filter], [energy_filter]]

        # Initialize the Tallies
        super(FissionXS, self).create_tallies(scores, filters, keys, estimator)

    def compute_xs(self):
        """Computes the multi-group fission cross sections using OpenMC
        tally arithmetic."""

        self._xs_tally = self.tallies['fission'] / self.tallies['flux']
        super(FissionXS, self).compute_xs()


class NuFissionXS(MultiGroupXS):

    def __init__(self, domain=None, domain_type=None,
                 groups=None, by_nuclide=False, name=''):
        super(NuFissionXS, self).__init__(domain, domain_type, groups, by_nuclide, name)
        self._rxn_type = 'nu-fission'

    def create_tallies(self):
        """Construct the OpenMC tallies needed to compute this cross section."""

        # Create a list of scores for each Tally to be created
        scores = ['flux', 'nu-fission']
        estimator = 'tracklength'
        keys = scores

        # Create the non-domain specific Filters for the Tallies
        group_edges = self.energy_groups.group_edges
        energy_filter = openmc.Filter('energy', group_edges)
        filters = [[energy_filter], [energy_filter]]

        # Initialize the Tallies
        super(NuFissionXS, self).create_tallies(scores, filters, keys, estimator)

    def compute_xs(self):
        """Computes the multi-group nu-fission cross sections using OpenMC
        tally arithmetic."""

        self._xs_tally = self.tallies['nu-fission'] / self.tallies['flux']
        super(NuFissionXS, self).compute_xs()


class ScatterXS(MultiGroupXS):

    def __init__(self, domain=None, domain_type=None,
                 groups=None, by_nuclide=False, name=''):
        super(ScatterXS, self).__init__(domain, domain_type, groups, by_nuclide, name)
        self._rxn_type = 'scatter'

    def create_tallies(self):
        """Construct the OpenMC tallies needed to compute this cross section."""

        # Create a list of scores for each Tally to be created
        scores = ['flux', 'scatter']
        estimator = 'tracklength'
        keys = scores

        # Create the non-domain specific Filters for the Tallies
        group_edges = self.energy_groups.group_edges
        energy_filter = openmc.Filter('energy', group_edges)
        filters = [[energy_filter], [energy_filter]]

        # Intialize the Tallies
        super(ScatterXS, self).create_tallies(scores, filters, keys, estimator)

    def compute_xs(self):
        """Computes the scattering multi-group cross sections using
        OpenMC tally arithmetic."""

        self._xs_tally = self.tallies['scatter'] / self.tallies['flux']
        super(ScatterXS, self).compute_xs()


class NuScatterXS(MultiGroupXS):

    def __init__(self, domain=None, domain_type=None,
                 groups=None, by_nuclide=False, name=''):
        super(NuScatterXS, self).__init__(domain, domain_type, groups, by_nuclide, name)
        self._rxn_type = 'nu-scatter'

    def create_tallies(self):
        """Construct the OpenMC tallies needed to compute this cross section."""

        # Create a list of scores for each Tally to be created
        scores = ['flux', 'nu-scatter']
        estimator = 'analog'
        keys = scores

        # Create the non-domain specific Filters for the Tallies
        group_edges = self.energy_groups.group_edges
        energy_filter = openmc.Filter('energy', group_edges)
        filters = [[energy_filter], [energy_filter]]

        # Initialize the Tallies
        super(NuScatterXS, self).create_tallies(scores, filters, keys, estimator)

    def compute_xs(self):
        """Computes the nu-scattering multi-group cross section using OpenMC
        tally arithmetic."""

        self._xs_tally = self.tallies['nu-scatter'] / self.tallies['flux']
        super(NuScatterXS, self).compute_xs()


class ScatterMatrixXS(MultiGroupXS):

    def __init__(self, domain=None, domain_type=None,
                 groups=None, by_nuclide=False, name=''):
        super(ScatterMatrixXS, self).__init__(domain, domain_type, groups, by_nuclide, name)
        self._rxn_type = 'scatter matrix'

    def create_tallies(self):
        """Construct the OpenMC tallies needed to compute this cross section."""

        group_edges = self.energy_groups.group_edges
        energy = openmc.Filter('energy', group_edges)
        energyout = openmc.Filter('energyout', group_edges)

        # Create a list of scores for each Tally to be created
        scores = ['flux', 'scatter', 'scatter-P1']
        filters = [[energy], [energy, energyout], [energyout]]

        estimator = 'analog'
        keys = scores

        # Initialize the Tallies
        super(ScatterMatrixXS, self).create_tallies(scores, filters, keys, estimator)

    def compute_xs(self, correction='P0'):
        """Computes the multi-group scattering matrix using OpenMC
        tally arithmetic.

        Parameters
        ----------
        correction : {'P0' or None}
            If 'P0', applies the P0 transport correction to the diagonal of the
            scattering matrix.

        """

        # If using P0 correction subtract scatter-P1 from the diagonal
        if correction == 'P0':
            scatter_p1 = self.tallies['scatter-P1']
            scatter_p1 = scatter_p1.get_slice(scores=['scatter-P1'])
            energy_filter = openmc.Filter(type='energy')
            energy_filter.bins = self.energy_groups.group_edges
            scatter_p1 = scatter_p1.diagonalize_filter(energy_filter)
            rxn_tally = self.tallies['scatter'] - scatter_p1
        else:
            rxn_tally = self.tallies['scatter']

        self._xs_tally = rxn_tally / self.tallies['flux']
        super(ScatterMatrixXS, self).compute_xs()

    def get_xs(self, in_groups='all', out_groups='all', subdomains='all',
               nuclides='all', order_groups='increasing',
               xs_type='macro', value='mean'):
        """Returns an array of multi-group cross sections.

        This method constructs a 2D NumPy array for the requested scattering
        matrix data data for one or more energy groups and subdomains.

        Parameters
        ----------
        in_groups : Iterable of Integral or 'all'
            Incoming energy groups of interest

        out_groups : Iterable of Integral or 'all'
            Outgoing energy groups of interest

        subdomains : Iterable of Integral or 'all'
            Subdomain IDs of interest

        nuclides : Iterable of str or 'all' or 'sum'
            A list of nuclide name strings (e.g., ['U-235', 'U-238']). The
            special string 'all' (default) will return the cross sections for
            all nuclides in the spatial domain. The special string 'sum' will
            return the cross section summed over all nuclides.

        xs_type: {'macro' or 'micro'}
            Return the macro or micro cross section in units of cm^-1 or barns

        xs_type: {'macro' or 'micro'}
            Return the macro or micro cross section in units of cm^-1 or barns

        value : str
            A string for the type of value to return - 'mean' (default),
            'std_dev' or 'rel_err' are accepted

        Returns
        -------
        xs : ndarray
            A NumPy array of the multi-group cross section indexed in the order
            each group and subdomain is listed in the parameters.

        Raises
        ------
        ValueError
            When this method is called before the multi-group cross section is
            computed from tally data.

        """

        if self.xs_tally is None:
            msg = 'Unable to get cross section since it has not been computed'
            raise ValueError(msg)

        cv.check_value('value', value, ['mean', 'std_dev', 'rel_err'])
        cv.check_value('xs_type', xs_type, ['macro', 'micro'])

        filters = []
        filter_bins = []

        # Construct a collection of the domain filter bins
        if subdomains != 'all':
            cv.check_iterable_type('subdomains', subdomains, Integral)
            for subdomain in subdomains:
                filters.append(self.domain_type)
                filter_bins.append((subdomain,))

        # Construct list of energy group bounds tuples for all requested groups
        if in_groups != 'all':
            cv.check_iterable_type('groups', in_groups, Integral)
            for group in in_groups:
                filters.append('energy')
                filter_bins.append((self.energy_groups.get_group_bounds(group),))

        # Construct list of energy group bounds tuples for all requested groups
        if out_groups != 'all':
            cv.check_iterable_type('groups', out_groups, Integral)
            for group in out_groups:
                filters.append('energyout')
                filter_bins.append((self.energy_groups.get_group_bounds(group),))

        # Construct a collection of the nuclides to retrieve from the xs tally
        # NOTE: We must not override the "nuclides" parameter since it is used
        # to retrieve atomic number densities for micro xs
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
            xs = self.xs_tally.get_values(filters=filters, filter_bins=filter_bins,
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
            # Reshape tally data array with separate axes for domain and energy
            if in_groups == 'all':
                num_in_groups = self.num_groups
            else:
                num_in_groups = len(in_groups)
            if out_groups == 'all':
                num_out_groups = self.num_groups
            else:
                num_out_groups = len(out_groups)
            num_subdomains = xs.shape[0] / (num_in_groups * num_out_groups)
            new_shape = (num_subdomains, num_in_groups, num_out_groups)
            new_shape += xs.shape[1:]
            xs = np.reshape(xs, new_shape)

            # Reverse energies to align with increasing energy groups
            xs = xs[:, ::-1, ::-1, :]

            # Reshape array to original axes (filters, nuclides, scores)
            new_shape = (num_subdomains * num_in_groups * num_out_groups,)
            new_shape += xs.shape[3:]
            xs = np.reshape(xs, new_shape)

        return xs

    def print_xs(self, subdomains='all', nuclides='all', xs_type='macro'):
        """Prints a string representation for the multi-group cross section.

        Parameters
        ----------
        subdomains : Iterable of Integral or 'all'
            The subdomain IDs of the cross sections to include in the report

        nuclides : Iterable of str or 'all' or 'sum'
            The nuclides of the cross-sections to include in the report. This
            may be a list of nuclide name strings (e.g., ['U-235', 'U-238']).
            The special string 'all' (default) will report the cross sections
            for all nuclides in the spatial domain. The special string 'sum'
            will report the cross sections summed over all nuclides.

        xs_type: {'macro' or 'micro'}
            Return the macro or micro cross section in units of cm^-1 or barns

        """

        # Construct a collection of the subdomains to report
        if subdomains != 'all':
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
        if self.xs_tally is None:
            print(string)
            return

        string += '{0: <16}\n'.format('\tEnergy Groups:')
        template = '{0: <12}Group {1} [{2: <10} - {3: <10}MeV]\n'

        # Loop over energy groups ranges
        for group in range(1, self.num_groups+1):
            bounds = self.energy_groups.get_group_bounds(group)
            string += template.format('', group, bounds[0], bounds[1])

        if subdomains == 'all':
            if self.domain_type == 'distribcell':
                subdomains = np.arange(self.num_subdomains, dtype=np.int)
            else:
                subdomains = [self.domain.id]

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
                            self.get_xs([in_group], [out_group], [subdomain],
                                        [nuclide], xs_type=xs_type, value='mean')
                        rel_err = \
                            self.get_xs([in_group], [out_group], [subdomain],
                                        [nuclide], xs_type=xs_type, value='rel_err') * 100
                        average = np.nan_to_num(average.flatten())[0]
                        rel_err = np.nan_to_num(rel_err.flatten())[0]
                        string += '{:1.2e} +/- {:1.2e}%'.format(average, rel_err)
                        string += '\n'
                    string += '\n'
                string += '\n'
            string += '\n'

        print(string)


class NuScatterMatrixXS(ScatterMatrixXS):

    def __init__(self, domain=None, domain_type=None,
                 groups=None, by_nuclide=False, name=''):
        super(NuScatterMatrixXS, self).__init__(domain, domain_type, groups, by_nuclide, name)
        self._rxn_type = 'nu-scatter matrix'

    def create_tallies(self):
        """Construct the OpenMC tallies needed to compute this cross section."""

        # Create a list of scores for each Tally to be created
        scores = ['flux', 'scatter', 'scatter-P1']
        estimator = 'analog'
        keys = scores

        # Create the non-domain specific Filters for the Tallies
        group_edges = self.energy_groups.group_edges
        energy = openmc.Filter('energy', group_edges)
        energyout = openmc.Filter('energyout', group_edges)
        filters = [[energy], [energy, energyout], [energyout]]

        # Intialize the Tallies
        super(ScatterMatrixXS, self).create_tallies(scores, filters, keys, estimator)

class Chi(MultiGroupXS):

    def __init__(self, domain=None, domain_type=None,
                 groups=None, by_nuclide=False, name=''):
        super(Chi, self).__init__(domain, domain_type, groups, by_nuclide, name)
        self._rxn_type = 'chi'

    def create_tallies(self):
        """Construct the OpenMC tallies needed to compute this cross section."""

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
        super(Chi, self).create_tallies(scores, filters, keys, estimator)

    def compute_xs(self):
        """Computes chi fission spectrum using OpenMC tally arithmetic."""

        # Retrieve the fission production tallies
        nu_fission_in = self.tallies['nu-fission-in']
        nu_fission_out = self.tallies['nu-fission-out']

        # Remove the coarse energy filter to keep it out of tally arithmetic
        nu_fission_in.remove_filter(nu_fission_in.filters[-1])

        # Compute chi
        self._xs_tally = nu_fission_out / nu_fission_in
        super(Chi, self).compute_xs()

    def get_xs(self, groups='all', subdomains='all', nuclides='all',
               order_groups='increasing', xs_type='macro', value='mean'):
        """Returns an array of multi-group cross sections.

        This method constructs a 2D NumPy array for the requested multi-group
        cross section data data for one or more energy groups and subdomains.

        Parameters
        ----------
        groups : Iterable of Integral or 'all'
            Energy groups of interest

        subdomains : Iterable of Integral or 'all'
            Subdomain IDs of interest

        nuclides : Iterable of str or 'all' or 'sum'
            A list of nuclide name strings (e.g., ['U-235', 'U-238']). The
            special string 'all' (default) will return the cross sections for
            all nuclides in the spatial domain. The special string 'sum' will
            return the cross section summed over all nuclides.

        xs_type: {'macro' or 'micro'}
            Return the macro or micro cross section in units of cm^-1 or barns

        xs_type: {'macro' or 'micro'}
            This parameter is not relevant for chi but is included here to
            mirror the parent MultiGroupXS.get_xs(...) class method

        value : str
            A string for the type of value to return - 'mean' (default),
            'std_dev' or 'rel_err' are accepted

        Returns
        -------
        xs : ndarray
            A NumPy array of the multi-group cross section indexed in the order
            each group, subdomain and nuclide is listed in the parameters.

        Raises
        ------
        ValueError
            When this method is called before the multi-group cross section is
            computed from tally data.

        """

        if self.xs_tally is None:
            msg = 'Unable to get cross section since it has not been computed'
            raise ValueError(msg)

        cv.check_value('value', value, ['mean', 'std_dev', 'rel_err'])
        cv.check_value('xs_type', xs_type, ['macro', 'micro'])

        filters = []
        filter_bins = []

        # Construct a collection of the domain filter bins
        if subdomains != 'all':
            cv.check_iterable_type('subdomains', subdomains, Integral)
            for subdomain in subdomains:
                filters.append(self.domain_type)
                filter_bins.append((subdomain,))

        # Construct list of energy group bounds tuples for all requested groups
        if groups != 'all':
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

                # Compute chi and store it as the xs_tally attribute so we can use
                # the generic get_xs routine
                xs_tally = nu_fission_out / nu_fission_in
                xs = xs_tally.get_values(filters=filters,
                                         filter_bins=filter_bins, value=value)

            # Get chi for all nuclides in the domain
            elif nuclides == 'all':
                nuclides = self.get_all_nuclides()
                xs = self.xs_tally.get_values(filters=filters, filter_bins=filter_bins,
                                              nuclides=nuclides, value=value)

            # Get chi for user-specified nuclides in the domain
            else:
                cv.check_iterable_type('nuclides', nuclides, basestring)
                xs = self.xs_tally.get_values(filters=filters, filter_bins=filter_bins,
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
            num_subdomains = xs.shape[0] / num_groups
            new_shape = (num_subdomains, num_groups) + xs.shape[1:]
            xs = np.reshape(xs, new_shape)

            # Reverse energies to align with increasing energy groups
            xs = xs[:, ::-1, :]

            # Reshape array to original axes (filters, nuclides, scores)
            new_shape = (num_subdomains * num_groups,) + new_shape[2:]
            xs = np.reshape(xs, new_shape)

        xs = np.nan_to_num(xs)
        return xs

    def get_pandas_dataframe(self, groups='all', nuclides='all',
                             xs_type='macro', summary=None):
        """Build a Pandas DataFrame for the MultiGroupXS data.

        This routine leverages the Tally.get_pandas_dataframe(...) routine, but
        renames the columns with terminology appropriate for cross section data.

        Parameters
        ----------
        groups : Iterable of Integral or 'all'
            Energy groups of interest

        nuclides : Iterable of str or 'all' or 'sum'
            The nuclides of the cross-sections to include in the dataframe. This
            may be a list of nuclide name strings (e.g., ['U-235', 'U-238']).
            The special string 'all' (default) will include the cross sections
            for all nuclides in the spatial domain. The special string 'sum'
            will include the cross sections summed over all nuclides.

        xs_type: {'macro' or 'micro'}
            Return macro or micro cross section in units of cm^-1 or barns

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

        # Build the dataframe using the parent class routine
        df = super(Chi, self).get_pandas_dataframe(groups, nuclides,
                                                   xs_type, summary)

        # If user requested micro cross sections, multiply by the atom
        # densities to cancel out division made by the parent class routine
        if xs_type == 'micro':
            if self.by_nuclide:
                densities = self.get_nuclide_densities(nuclides)
            else:
                densities = self.get_nuclide_densities('sum')
            tile_factor = df.shape[0] / len(densities)
            df['mean'] *= np.tile(densities, tile_factor)
            df['std. dev.'] *= np.tile(densities, tile_factor)

        return df