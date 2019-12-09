from collections import OrderedDict
from collections.abc import Iterable
import itertools
from numbers import Integral
import warnings
import os
import sys
import copy
from abc import ABCMeta

import numpy as np

import openmc
from openmc.mgxs import MGXS
from openmc.mgxs.mgxs import _DOMAIN_TO_FILTER
import openmc.checkvalue as cv


# Supported cross section types
MDGXS_TYPES = ['delayed-nu-fission',
               'chi-delayed',
               'beta',
               'decay-rate',
               'delayed-nu-fission matrix']

# Maximum number of delayed groups, from src/constants.F90
MAX_DELAYED_GROUPS = 8


class MDGXS(MGXS):
    """An abstract multi-delayed-group cross section for some energy and delayed
    group structures within some spatial domain.

    This class can be used for both OpenMC input generation and tally data
    post-processing to compute spatially-homogenized and energy-integrated
    multi-group and multi-delayed-group cross sections for downstream
    neutronics calculations.

    NOTE: Users should instantiate the subclasses of this abstract class.

    Parameters
    ----------
    domain : openmc.Material or openmc.Cell or openmc.Universe or openmc.RegularMesh
        The domain for spatial homogenization
    domain_type : {'material', 'cell', 'distribcell', 'universe', 'mesh'}
        The domain type for spatial homogenization
    energy_groups : openmc.mgxs.EnergyGroups
        The energy group structure for energy condensation
    by_nuclide : bool
        If true, computes cross sections for each nuclide in domain
    name : str, optional
        Name of the multi-group cross section. Used as a label to identify
        tallies in OpenMC 'tallies.xml' file.
    delayed_groups : list of int
        Delayed groups to filter out the xs
    num_polar : Integral, optional
        Number of equi-width polar angle bins for angle discretization;
        defaults to one bin
    num_azimuthal : Integral, optional
        Number of equi-width azimuthal angle bins for angle discretization;
        defaults to one bin

    Attributes
    ----------
    name : str, optional
        Name of the multi-group cross section
    rxn_type : str
        Reaction type (e.g., 'chi-delayed', 'beta', etc.)
    by_nuclide : bool
        If true, computes cross sections for each nuclide in domain
    domain : openmc.Material or openmc.Cell or openmc.Universe or openmc.RegularMesh
        Domain for spatial homogenization
    domain_type : {'material', 'cell', 'distribcell', 'universe', 'mesh'}
        Domain type for spatial homogenization
    energy_groups : openmc.mgxs.EnergyGroups
        Energy group structure for energy condensation
    delayed_groups : list of int
        Delayed groups to filter out the xs
    num_polar : Integral
        Number of equi-width polar angle bins for angle discretization
    num_azimuthal : Integral
        Number of equi-width azimuthal angle bins for angle discretization
    tally_trigger : openmc.Trigger
        An (optional) tally precision trigger given to each tally used to
        compute the cross section
    scores : list of str
        The scores in each tally used to compute the multi-group cross section
    filters : list of openmc.Filter
        The filters in each tally used to compute the multi-group cross section
    tally_keys : list of str
        The keys into the tallies dictionary for each tally used to compute
        the multi-group cross section
    estimator : {'tracklength', 'analog'}
        The tally estimator used to compute the multi-group cross section
    tallies : collections.OrderedDict
        OpenMC tallies needed to compute the multi-group cross section
    rxn_rate_tally : openmc.Tally
        Derived tally for the reaction rate tally used in the numerator to
        compute the multi-group cross section. This attribute is None
        unless the multi-group cross section has been computed.
    xs_tally : openmc.Tally
        Derived tally for the multi-group cross section. This attribute
        is None unless the multi-group cross section has been computed.
    num_subdomains : int
        The number of subdomains is unity for 'material', 'cell', and 'universe'
        domain types. This is equal to the number of cell instances for
        'distribcell' domain types (it is equal to unity prior to loading
        tally data from a statepoint file) and the number of mesh cells for
        'mesh' domain types.
    num_nuclides : int
        The number of nuclides for which the multi-group cross section is
        being tracked. This is unity if the by_nuclide attribute is False.
    nuclides : Iterable of str or 'sum'
        The optional user-specified nuclides for which to compute cross
        sections (e.g., 'U-238', 'O-16'). If by_nuclide is True but nuclides
        are not specified by the user, all nuclides in the spatial domain
        are included. This attribute is 'sum' if by_nuclide is false.
    sparse : bool
        Whether or not the MGXS' tallies use SciPy's LIL sparse matrix format
        for compressed data storage
    loaded_sp : bool
        Whether or not a statepoint file has been loaded with tally data
    derived : bool
        Whether or not the MGXS is merged from one or more other MGXS
    hdf5_key : str
        The key used to index multi-group cross sections in an HDF5 data store

    """

    def __init__(self, domain=None, domain_type=None, energy_groups=None,
                 delayed_groups=None, by_nuclide=False, name='',
                 num_polar=1, num_azimuthal=1):
        super().__init__(domain, domain_type, energy_groups, by_nuclide, name,
                         num_polar, num_azimuthal)

        self._delayed_groups = None

        if delayed_groups is not None:
            self.delayed_groups = delayed_groups

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
            clone._delayed_groups = copy.deepcopy(self.delayed_groups, memo)
            clone._num_polar = self.num_polar
            clone._num_azimuthal = self.num_azimuthal
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
    def _dont_squeeze(self):
        """Create a tuple of axes which should not be removed during the get_xs
        process
        """
        if self.num_polar > 1 or self.num_azimuthal > 1:
            return (0, 1, 3, 4)
        else:
            return (1, 2)

    @property
    def delayed_groups(self):
        return self._delayed_groups

    @property
    def num_delayed_groups(self):
        if self.delayed_groups is None:
            return 1
        else:
            return len(self.delayed_groups)

    @delayed_groups.setter
    def delayed_groups(self, delayed_groups):

        if delayed_groups is not None:

            cv.check_type('delayed groups', delayed_groups, list, int)
            cv.check_greater_than('num delayed groups', len(delayed_groups), 0)

            # Check that the groups are within [1, MAX_DELAYED_GROUPS]
            for group in delayed_groups:
                cv.check_greater_than('delayed group', group, 0)
                cv.check_less_than('delayed group', group, MAX_DELAYED_GROUPS,
                                   equality=True)

        self._delayed_groups = delayed_groups

    @property
    def filters(self):

        # Create the non-domain specific Filters for the Tallies
        group_edges = self.energy_groups.group_edges
        energy_filter = openmc.EnergyFilter(group_edges)

        if self.delayed_groups != None:
            delayed_filter = openmc.DelayedGroupFilter(self.delayed_groups)
            filters = [[energy_filter], [delayed_filter, energy_filter]]
        else:
            filters = [[energy_filter], [energy_filter]]
        return self._add_angle_filters(filters)

    @staticmethod
    def get_mgxs(mdgxs_type, domain=None, domain_type=None, energy_groups=None,
                 delayed_groups=None, by_nuclide=False, name='',
                 num_polar=1, num_azimuthal=1):
        """Return a MDGXS subclass object for some energy group structure within
        some spatial domain for some reaction type.

        This is a factory method which can be used to quickly create MDGXS
        subclass objects for various reaction types.

        Parameters
        ----------
        mdgxs_type : {'delayed-nu-fission', 'chi-delayed', 'beta', 'decay-rate', 'delayed-nu-fission matrix'}
            The type of multi-delayed-group cross section object to return
        domain : openmc.Material or openmc.Cell or openmc.Universe or
            openmc.RegularMesh
            The domain for spatial homogenization
        domain_type : {'material', 'cell', 'distribcell', 'universe', 'mesh'}
            The domain type for spatial homogenization
        energy_groups : openmc.mgxs.EnergyGroups
            The energy group structure for energy condensation
        by_nuclide : bool
            If true, computes cross sections for each nuclide in domain.
            Defaults to False
        name : str, optional
            Name of the multi-group cross section. Used as a label to identify
            tallies in OpenMC 'tallies.xml' file. Defaults to the empty string.
        delayed_groups : list of int
            Delayed groups to filter out the xs
        num_polar : Integral, optional
            Number of equi-width polar angle bins for angle discretization;
            defaults to one bin
        num_azimuthal : Integral, optional
            Number of equi-width azimuthal angle bins for angle discretization;
            defaults to one bin
        Returns
        -------
        openmc.mgxs.MDGXS
            A subclass of the abstract MDGXS class for the multi-delayed-group
            cross section type requested by the user

        """

        cv.check_value('mdgxs_type', mdgxs_type, MDGXS_TYPES)

        if mdgxs_type == 'delayed-nu-fission':
            mdgxs = DelayedNuFissionXS(domain, domain_type, energy_groups,
                                       delayed_groups)
        elif mdgxs_type == 'chi-delayed':
            mdgxs = ChiDelayed(domain, domain_type, energy_groups,
                               delayed_groups)
        elif mdgxs_type == 'beta':
            mdgxs = Beta(domain, domain_type, energy_groups, delayed_groups)
        elif mdgxs_type == 'decay-rate':
            mdgxs = DecayRate(domain, domain_type, energy_groups, delayed_groups)
        elif mdgxs_type == 'delayed-nu-fission matrix':
            mdgxs = DelayedNuFissionMatrixXS(domain, domain_type, energy_groups,
                                             delayed_groups)

        mdgxs.by_nuclide = by_nuclide
        mdgxs.name = name
        mdgxs.num_polar = num_polar
        mdgxs.num_azimuthal = num_azimuthal
        return mdgxs

    def get_xs(self, groups='all', subdomains='all', nuclides='all',
               xs_type='macro', order_groups='increasing',
               value='mean', delayed_groups='all', squeeze=True, **kwargs):
        """Returns an array of multi-delayed-group cross sections.

        This method constructs a 4D NumPy array for the requested
        multi-delayed-group cross section data for one or more
        subdomains (1st dimension), delayed groups (2nd demension),
        energy groups (3rd dimension), and nuclides (4th dimension).

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
        value : {'mean', 'std_dev', 'rel_err'}
            A string for the type of value to return. Defaults to 'mean'.
        delayed_groups : list of int or 'all'
            Delayed groups of interest. Defaults to 'all'.
        squeeze : bool
            A boolean representing whether to eliminate the extra dimensions
            of the multi-dimensional array to be returned. Defaults to True.

        Returns
        -------
        numpy.ndarray
            A NumPy array of the multi-group cross section indexed in the order
            each group, subdomain and nuclide is listed in the parameters.

        Raises
        ------
        ValueError
            When this method is called before the multi-delayed-group cross
            section is computed from tally data.

        """

        cv.check_value('value', value, ['mean', 'std_dev', 'rel_err'])
        cv.check_value('xs_type', xs_type, ['macro', 'micro'])

        # FIXME: Unable to get microscopic xs for mesh domain because the mesh
        # cells do not know the nuclide densities in each mesh cell.
        if self.domain_type == 'mesh' and xs_type == 'micro':
            msg = 'Unable to get micro xs for mesh domain since the mesh ' \
                  'cells do not know the nuclide densities in each mesh cell.'
            raise ValueError(msg)

        filters = []
        filter_bins = []

        # Construct a collection of the domain filter bins
        if not isinstance(subdomains, str):
            cv.check_iterable_type('subdomains', subdomains, Integral,
                                   max_depth=3)
            for subdomain in subdomains:
                filters.append(_DOMAIN_TO_FILTER[self.domain_type])
                filter_bins.append((subdomain,))

        # Construct list of energy group bounds tuples for all requested groups
        if not isinstance(groups, str):
            cv.check_iterable_type('groups', groups, Integral)
            for group in groups:
                filters.append(openmc.EnergyFilter)
                filter_bins.append(
                    (self.energy_groups.get_group_bounds(group),))

        # Construct list of delayed group tuples for all requested groups
        if not isinstance(delayed_groups, str):
            cv.check_type('delayed groups', delayed_groups, list, int)
            for delayed_group in delayed_groups:
                filters.append(openmc.DelayedGroupFilter)
                filter_bins.append((delayed_group,))

        # Construct a collection of the nuclides to retrieve from the xs tally
        if self.by_nuclide:
            if nuclides == 'all' or nuclides == 'sum' or nuclides == ['sum']:
                query_nuclides = self.get_nuclides()
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
            xs = self.xs_tally.get_values(filters=filters,
                                          filter_bins=filter_bins,
                                          nuclides=query_nuclides, value=value)

        # Divide by atom number densities for microscopic cross sections
        if xs_type == 'micro':
            if self.by_nuclide:
                densities = self.get_nuclide_densities(nuclides)
            else:
                densities = self.get_nuclide_densities('sum')
            if value == 'mean' or value == 'std_dev':
                xs /= densities[np.newaxis, :, np.newaxis]

        # Eliminate the trivial score dimension
        xs = np.squeeze(xs, axis=len(xs.shape) - 1)
        xs = np.nan_to_num(xs)

        if groups == 'all':
            num_groups = self.num_groups
        else:
            num_groups = len(groups)

        if delayed_groups == 'all':
            num_delayed_groups = self.num_delayed_groups
        else:
            num_delayed_groups = len(delayed_groups)

        # Reshape tally data array with separate axes for domain,
        # energy groups, delayed groups, and nuclides
        # Accommodate the polar and azimuthal bins if needed
        num_subdomains = \
            int(xs.shape[0] / (num_groups * num_delayed_groups *
                               self.num_polar * self.num_azimuthal))
        if self.num_polar > 1 or self.num_azimuthal > 1:
            new_shape = (self.num_polar, self.num_azimuthal, num_subdomains,
                         num_delayed_groups, num_groups)
        else:
            new_shape = (num_subdomains, num_delayed_groups, num_groups)
        new_shape += xs.shape[1:]
        xs = np.reshape(xs, new_shape)

        # Reverse data if user requested increasing energy groups since
        # tally data is stored in order of increasing energies
        if order_groups == 'increasing':
            xs = xs[..., ::-1, :]

        if squeeze:
            # We want to squeeze out everything but the polar, azimuthal,
            # delayed group, and energy group data.
            xs = self._squeeze_xs(xs)

        return xs

    def get_slice(self, nuclides=[], groups=[], delayed_groups=[]):
        """Build a sliced MDGXS for the specified nuclides, energy groups,
        and delayed groups.

        This method constructs a new MDGXS to encapsulate a subset of the data
        represented by this MDGXS. The subset of data to include in the tally
        slice is determined by the nuclides, energy groups, delayed groups
        specified in the input parameters.

        Parameters
        ----------
        nuclides : list of str
            A list of nuclide name strings
            (e.g., ['U-235', 'U-238']; default is [])
        groups : list of int
            A list of energy group indices starting at 1 for the high energies
            (e.g., [1, 2, 3]; default is [])
        delayed_groups : list of int
            A list of delayed group indices
            (e.g., [1, 2, 3]; default is [])

        Returns
        -------
        openmc.mgxs.MDGXS
            A new MDGXS object which encapsulates the subset of data requested
            for the nuclide(s) and/or energy group(s) and/or delayed group(s)
            requested in the parameters.

        """

        cv.check_iterable_type('nuclides', nuclides, str)
        cv.check_iterable_type('energy_groups', groups, Integral)
        cv.check_type('delayed groups', delayed_groups, list, int)

        # Build lists of filters and filter bins to slice
        filters = []
        filter_bins = []

        if len(groups) != 0:
            energy_bins = []
            for group in groups:
                group_bounds = self.energy_groups.get_group_bounds(group)
                energy_bins.append(group_bounds)
            filter_bins.append(tuple(energy_bins))
            filters.append(openmc.EnergyFilter)

        if len(delayed_groups) != 0:
            filter_bins.append(tuple(delayed_groups))
            filters.append(openmc.DelayedGroupFilter)

        # Clone this MGXS to initialize the sliced version
        slice_xs = copy.deepcopy(self)
        slice_xs._rxn_rate_tally = None
        slice_xs._xs_tally = None

        # Slice each of the tallies across nuclides and energy groups
        for tally_type, tally in slice_xs.tallies.items():
            slice_nuclides = [nuc for nuc in nuclides if nuc in tally.nuclides]
            if filters != []:
                tally_slice = tally.get_slice(filters=filters,
                                              filter_bins=filter_bins,
                                              nuclides=slice_nuclides)
            else:
                tally_slice = tally.get_slice(nuclides=slice_nuclides)
            slice_xs.tallies[tally_type] = tally_slice

        # Assign sliced energy group structure to sliced MDGXS
        if groups:
            new_group_edges = []
            for group in groups:
                group_edges = self.energy_groups.get_group_bounds(group)
                new_group_edges.extend(group_edges)
            new_group_edges = np.unique(new_group_edges)
            slice_xs.energy_groups.group_edges = sorted(new_group_edges)

        # Assign sliced delayed group structure to sliced MDGXS
        if delayed_groups:
            slice_xs.delayed_groups = delayed_groups

        # Assign sliced nuclides to sliced MGXS
        if nuclides:
            slice_xs.nuclides = nuclides

        slice_xs.sparse = self.sparse
        return slice_xs

    def merge(self, other):
        """Merge another MGXS with this one

        MGXS are only mergeable if their energy groups and nuclides are either
        identical or mutually exclusive. If results have been loaded from a
        statepoint, then MGXS are only mergeable along one and only one of
        energy groups or nuclides.

        Parameters
        ----------
        other : openmc.mgxs.MDGXS
            MDGXS to merge with this one

        Returns
        -------
        merged_mdgxs : openmc.mgxs.MDGXS
            Merged MDGXS

        """

        merged_mdgxs = super().merge(other)

        # Merge delayed groups
        if self.delayed_groups != other.delayed_groups:
            merged_mdgxs.delayed_groups = list(set(self.delayed_groups +
                                                   other.delayed_groups))

        return merged_mdgxs

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

        if self.delayed_groups is None:
            super().print_xs(subdomains, nuclides, xs_type)
            return

        # Construct a collection of the subdomains to report
        if not isinstance(subdomains, str):
            cv.check_iterable_type('subdomains', subdomains, Integral)
        elif self.domain_type == 'distribcell':
            subdomains = np.arange(self.num_subdomains, dtype=np.int)
        elif self.domain_type == 'mesh':
            xyz = [range(1, x + 1) for x in self.domain.dimension]
            subdomains = list(itertools.product(*xyz))
        else:
            subdomains = [self.domain.id]

        # Construct a collection of the nuclides to report
        if self.by_nuclide:
            if nuclides == 'all':
                nuclides = self.get_nuclides()
            elif nuclides == 'sum':
                nuclides = ['sum']
            else:
                cv.check_iterable_type('nuclides', nuclides, str)
        else:
            nuclides = ['sum']

        cv.check_value('xs_type', xs_type, ['macro', 'micro'])

        # Build header for string with type and domain info
        string = 'Multi-Delayed-Group XS\n'
        string += '{0: <16}=\t{1}\n'.format('\tReaction Type', self.rxn_type)
        string += '{0: <16}=\t{1}\n'.format('\tDomain Type', self.domain_type)
        string += '{0: <16}=\t{1}\n'.format('\tDomain ID', self.domain.id)

        # Generate the header for an individual XS
        xs_header = '\tCross Sections [{0}]:'.format(self.get_units(xs_type))

        # If cross section data has not been computed, only print string header
        if self.tallies is None:
            print(string)
            return

        # Set polar/azimuthal bins
        if self.num_polar > 1 or self.num_azimuthal > 1:
            polar_bins = np.linspace(0., np.pi, num=self.num_polar + 1,
                                     endpoint=True)
            azimuthal_bins = np.linspace(-np.pi, np.pi,
                                         num=self.num_azimuthal + 1,
                                         endpoint=True)

        # Loop over all subdomains
        for subdomain in subdomains:

            if self.domain_type == 'distribcell' or self.domain_type == 'mesh':
                string += '{0: <16}=\t{1}\n'.format('\tSubdomain', subdomain)

            # Loop over all Nuclides
            for nuclide in nuclides:

                # Build header for nuclide type
                if nuclide != 'sum':
                    string += '{0: <16}=\t{1}\n'.format('\tNuclide', nuclide)

                # Add the cross section header
                string += '{0: <16}\n'.format(xs_header)

                for delayed_group in self.delayed_groups:

                    template = '{0: <12}Delayed Group {1}:\t'
                    string += template.format('', delayed_group)
                    string += '\n'

                    template = '{0: <12}Group {1} [{2: <10} - {3: <10}eV]:\t'

                    average_xs = self.get_xs(nuclides=[nuclide],
                                             subdomains=[subdomain],
                                             xs_type=xs_type, value='mean',
                                             delayed_groups=[delayed_group])
                    rel_err_xs = self.get_xs(nuclides=[nuclide],
                                             subdomains=[subdomain],
                                             xs_type=xs_type, value='rel_err',
                                             delayed_groups=[delayed_group])
                    rel_err_xs = rel_err_xs * 100.

                    if self.num_polar > 1 or self.num_azimuthal > 1:
                        # Loop over polar, azimuthal, and energy group ranges
                        for pol in range(len(polar_bins) - 1):
                            pol_low, pol_high = polar_bins[pol: pol + 2]
                            for azi in range(len(azimuthal_bins) - 1):
                                azi_low, azi_high = azimuthal_bins[azi: azi + 2]
                                string += '\t\tPolar Angle: [{0:5f} - {1:5f}]'.format(
                                    pol_low, pol_high) + \
                                    '\tAzimuthal Angle: [{0:5f} - {1:5f}]'.format(
                                    azi_low, azi_high) + '\n'
                                for group in range(1, self.num_groups + 1):
                                    bounds = \
                                        self.energy_groups.get_group_bounds(group)
                                    string += '\t' + template.format('', group,
                                                                     bounds[0],
                                                                     bounds[1])
                                    string += '{0:.2e} +/- {1:.2e}%'.format(
                                        average_xs[pol, azi, group - 1],
                                        rel_err_xs[pol, azi, group - 1])
                                    string += '\n'
                                string += '\n'
                    else:
                        # Loop over energy groups ranges
                        for group in range(1, self.num_groups+1):
                            bounds = self.energy_groups.get_group_bounds(group)
                            string += template.format('', group, bounds[0], bounds[1])
                            string += '{0:.2e} +/- {1:.2e}%'.format(
                                average_xs[group - 1], rel_err_xs[group - 1])
                            string += '\n'
                    string += '\n'
                string += '\n'

        print(string)

    def export_xs_data(self, filename='mgxs', directory='mgxs',
                       format='csv', groups='all', xs_type='macro',
                       delayed_groups='all'):
        """Export the multi-delayed-group cross section data to a file.

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
        delayed_groups : list of int or 'all'
            Delayed groups of interest. Defaults to 'all'.

        """

        cv.check_type('filename', filename, str)
        cv.check_type('directory', directory, str)
        cv.check_value('format', format, ['csv', 'excel', 'pickle', 'latex'])
        cv.check_value('xs_type', xs_type, ['macro', 'micro'])

        # Make directory if it does not exist
        if not os.path.exists(directory):
            os.makedirs(directory)

        filename = os.path.join(directory, filename)
        filename = filename.replace(' ', '-')

        # Get a Pandas DataFrame for the data
        df = self.get_pandas_dataframe(groups=groups, xs_type=xs_type,
                                       delayed_groups=delayed_groups)

        # Export the data using Pandas IO API
        if format == 'csv':
            df.to_csv(filename + '.csv', index=False)
        elif format == 'excel':
            if self.domain_type == 'mesh':
                df.to_excel(filename + '.xls')
            else:
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
                             xs_type='macro', paths=True,
                             delayed_groups='all'):
        """Build a Pandas DataFrame for the MDGXS data.

        This method leverages :meth:`openmc.Tally.get_pandas_dataframe`, but
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
        paths : bool, optional
            Construct columns for distribcell tally filters (default is True).
            The geometric information in the Summary object is embedded into
            a Multi-index column with a geometric "path" to each distribcell
            instance.
        delayed_groups : list of int or 'all'
            Delayed groups of interest. Defaults to 'all'.

        Returns
        -------
        pandas.DataFrame
            A Pandas DataFrame for the cross section data.

        Raises
        ------
        ValueError
            When this method is called before the multi-delayed-group cross
            section is computed from tally data.

        """

        if not isinstance(groups, str):
            cv.check_iterable_type('groups', groups, Integral)
        if nuclides != 'all' and nuclides != 'sum':
            cv.check_iterable_type('nuclides', nuclides, str)
        if not isinstance(delayed_groups, str):
            cv.check_type('delayed groups', delayed_groups, list, int)

        cv.check_value('xs_type', xs_type, ['macro', 'micro'])

        # Get a Pandas DataFrame from the derived xs tally
        if self.by_nuclide and nuclides == 'sum':

            # Use tally summation to sum across all nuclides
            xs_tally = self.xs_tally.summation(nuclides=self.get_nuclides())
            df = xs_tally.get_pandas_dataframe(paths=paths)

            # Remove nuclide column since it is homogeneous and redundant
            if self.domain_type == 'mesh':
                df.drop('sum(nuclide)', axis=1, level=0, inplace=True)
            else:
                df.drop('sum(nuclide)', axis=1, inplace=True)

        # If the user requested a specific set of nuclides
        elif self.by_nuclide and nuclides != 'all':
            xs_tally = self.xs_tally.get_slice(nuclides=nuclides)
            df = xs_tally.get_pandas_dataframe(paths=paths)

        # If the user requested all nuclides, keep nuclide column in dataframe
        else:
            df = self.xs_tally.get_pandas_dataframe(paths=paths)

        # Remove the score column since it is homogeneous and redundant
        if self.domain_type == 'mesh':
            df = df.drop('score', axis=1, level=0)
        else:
            df = df.drop('score', axis=1)

        # Convert azimuthal, polar, energy in and energy out bin values in to
        # bin indices
        columns = self._df_convert_columns_to_bins(df)

        # Select out those groups the user requested
        if not isinstance(groups, str):
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
            densities = np.repeat(densities, len(self.rxn_rate_tally.scores))
            tile_factor = int(df.shape[0] / len(densities))
            df['mean'] /= np.tile(densities, tile_factor)
            df['std. dev.'] /= np.tile(densities, tile_factor)

        # Sort the dataframe by domain type id (e.g., distribcell id) and
        # energy groups such that data is from fast to thermal
        if self.domain_type == 'mesh':
            mesh_str = 'mesh {0}'.format(self.domain.id)
            df.sort_values(by=[(mesh_str, 'x'), (mesh_str, 'y'),
                               (mesh_str, 'z')] + columns, inplace=True)
        else:
            df.sort_values(by=[self.domain_type] + columns, inplace=True)

        return df


class ChiDelayed(MDGXS):
    r"""The delayed fission spectrum.

    This class can be used for both OpenMC input generation and tally data
    post-processing to compute spatially-homogenized and energy-integrated
    multi-group and multi-delayed-group cross sections for multi-group
    neutronics calculations. At a minimum, one needs to set the
    :attr:`ChiDelayed.energy_groups` and :attr:`ChiDelayed.domain` properties.
    Tallies for the flux and appropriate reaction rates over the specified
    domain are generated automatically via the :attr:`ChiDelayed.tallies`
    property, which can then be appended to a :class:`openmc.Tallies` instance.

    For post-processing, the :meth:`MGXS.load_from_statepoint` will pull in the
    necessary data to compute multi-group cross sections from a
    :class:`openmc.StatePoint` instance. The derived multi-group cross
    section can then be obtained from the :attr:`ChiDelayed.xs_tally` property.

    For a spatial domain :math:`V`, energy group :math:`[E_g,E_{g-1}]`, and
    delayed group :math:`d`, the delayed fission spectrum is calculated as:

    .. math::

       \begin{aligned}
       \langle \nu^d \sigma_{f,g' \rightarrow g} \phi \rangle &= \int_{r \in V}
       dr \int_{4\pi} d\Omega' \int_0^\infty dE' \int_{E_g}^{E_{g-1}} dE \;
       \chi(E) \nu^d \sigma_f (r, E') \psi(r, E', \Omega')\\
       \langle \nu^d \sigma_f \phi \rangle &= \int_{r \in V} dr \int_{4\pi}
       d\Omega' \int_0^\infty dE' \int_0^\infty dE \; \chi(E) \nu^d \sigma_f (r,
       E') \psi(r, E', \Omega') \\
       \chi_g^d &= \frac{\langle \nu^d \sigma_{f,g' \rightarrow g} \phi \rangle}
       {\langle \nu^d \sigma_f \phi \rangle}
       \end{aligned}

    Parameters
    ----------
    domain : openmc.Material or openmc.Cell or openmc.Universe or openmc.RegularMesh
        The domain for spatial homogenization
    domain_type : {'material', 'cell', 'distribcell', 'universe', 'mesh'}
        The domain type for spatial homogenization
    groups : openmc.mgxs.EnergyGroups
        The energy group structure for energy condensation
    by_nuclide : bool
        If true, computes cross sections for each nuclide in domain
    name : str, optional
        Name of the multi-group cross section. Used as a label to identify
        tallies in OpenMC 'tallies.xml' file.
    delayed_groups : list of int
        Delayed groups to filter out the xs
    num_polar : Integral, optional
        Number of equi-width polar angle bins for angle discretization;
        defaults to one bin
    num_azimuthal : Integral, optional
        Number of equi-width azimuthal angle bins for angle discretization;
        defaults to one bin

    Attributes
    ----------
    name : str, optional
        Name of the multi-group cross section
    rxn_type : str
        Reaction type (e.g., 'total', 'nu-fission', etc.)
    by_nuclide : bool
        If true, computes cross sections for each nuclide in domain
    domain : openmc.Material or openmc.Cell or openmc.Universe or openmc.RegularMesh
        Domain for spatial homogenization
    domain_type : {'material', 'cell', 'distribcell', 'universe', 'mesh'}
        Domain type for spatial homogenization
    energy_groups : openmc.mgxs.EnergyGroups
        Energy group structure for energy condensation
    delayed_groups : list of int
        Delayed groups to filter out the xs
    num_polar : Integral
        Number of equi-width polar angle bins for angle discretization
    num_azimuthal : Integral
        Number of equi-width azimuthal angle bins for angle discretization
    tally_trigger : openmc.Trigger
        An (optional) tally precision trigger given to each tally used to
        compute the cross section
    scores : list of str
        The scores in each tally used to compute the multi-group cross section
    filters : list of openmc.Filter
        The filters in each tally used to compute the multi-group cross section
    tally_keys : list of str
        The keys into the tallies dictionary for each tally used to compute
        the multi-group cross section
    estimator : 'analog'
        The tally estimator used to compute the multi-group cross section
    tallies : collections.OrderedDict
        OpenMC tallies needed to compute the multi-group cross section. The keys
        are strings listed in the :attr:`ChiDelayed.tally_keys` property and
        values are instances of :class:`openmc.Tally`.
    rxn_rate_tally : openmc.Tally
        Derived tally for the reaction rate tally used in the numerator to
        compute the multi-group cross section. This attribute is None
        unless the multi-group cross section has been computed.
    xs_tally : openmc.Tally
        Derived tally for the multi-group cross section. This attribute
        is None unless the multi-group cross section has been computed.
    num_subdomains : int
        The number of subdomains is unity for 'material', 'cell' and 'universe'
        domain types. This is equal to the number of cell instances
        for 'distribcell' domain types (it is equal to unity prior to loading
        tally data from a statepoint file).
    num_nuclides : int
        The number of nuclides for which the multi-group cross section is
        being tracked. This is unity if the by_nuclide attribute is False.
    nuclides : Iterable of str or 'sum'
        The optional user-specified nuclides for which to compute cross
        sections (e.g., 'U-238', 'O-16'). If by_nuclide is True but nuclides
        are not specified by the user, all nuclides in the spatial domain
        are included. This attribute is 'sum' if by_nuclide is false.
    sparse : bool
        Whether or not the MGXS' tallies use SciPy's LIL sparse matrix format
        for compressed data storage
    loaded_sp : bool
        Whether or not a statepoint file has been loaded with tally data
    derived : bool
        Whether or not the MGXS is merged from one or more other MGXS
    hdf5_key : str
        The key used to index multi-group cross sections in an HDF5 data store

    """

    def __init__(self, domain=None, domain_type=None, energy_groups=None,
                 delayed_groups=None, by_nuclide=False, name='',
                 num_polar=1, num_azimuthal=1):
        super().__init__(domain, domain_type, energy_groups, delayed_groups,
                         by_nuclide, name, num_polar, num_azimuthal)
        self._rxn_type = 'chi-delayed'
        self._estimator = 'analog'

    @property
    def scores(self):
        return ['delayed-nu-fission', 'delayed-nu-fission']

    @property
    def filters(self):
        # Create the non-domain specific Filters for the Tallies
        group_edges = self.energy_groups.group_edges
        energyout = openmc.EnergyoutFilter(group_edges)
        energyin = openmc.EnergyFilter([group_edges[0], group_edges[-1]])
        if self.delayed_groups is not None:
            delayed_filter = openmc.DelayedGroupFilter(self.delayed_groups)
            filters = [[delayed_filter, energyin], [delayed_filter, energyout]]
        else:
            filters = [[energyin], [energyout]]

        return self._add_angle_filters(filters)

    @property
    def tally_keys(self):
        return ['delayed-nu-fission-in', 'delayed-nu-fission-out']

    @property
    def rxn_rate_tally(self):
        if self._rxn_rate_tally is None:
            self._rxn_rate_tally = self.tallies['delayed-nu-fission-out']
            self._rxn_rate_tally.sparse = self.sparse
        return self._rxn_rate_tally

    @property
    def xs_tally(self):

        if self._xs_tally is None:
            delayed_nu_fission_in = self.tallies['delayed-nu-fission-in']

            # Remove coarse energy filter to keep it out of tally arithmetic
            energy_filter = delayed_nu_fission_in.find_filter(
                openmc.EnergyFilter)
            delayed_nu_fission_in.remove_filter(energy_filter)

            # Compute chi
            self._xs_tally = self.rxn_rate_tally / delayed_nu_fission_in
            super()._compute_xs()

            # Add the coarse energy filter back to the nu-fission tally
            delayed_nu_fission_in.filters.append(energy_filter)

        return self._xs_tally

    def get_homogenized_mgxs(self, other_mgxs):
        """Construct a homogenized MGXS with other MGXS objects.

        This method constructs a new MGXS object that is the flux-weighted
        combination of two MGXS objects. It is equivalent to what one would
        obtain if the tally spatial domain were designed to encompass the
        individual domains for both MGXS objects.

        Parameters
        ----------
        other_mgxs : openmc.mgxs.MGXS or Iterable of openmc.mgxs.MGXS
            The MGXS to homogenize with this one.

        Returns
        -------
        openmc.mgxs.MGXS
            A new homogenized MGXS

        Raises
        ------
        ValueError
            If the other_mgxs is of a different type.

        """

        return self._get_homogenized_mgxs(other_mgxs, 'delayed-nu-fission-in')


    def get_slice(self, nuclides=[], groups=[], delayed_groups=[]):
        """Build a sliced ChiDelayed for the specified nuclides and energy
           groups.

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
        delayed_groups : list of int
            A list of delayed group indices
            (e.g., [1, 2, 3]; default is [])

        Returns
        -------
        openmc.mgxs.MDGXS
            A new MDGXS which encapsulates the subset of data requested
            for the nuclide(s) and/or energy group(s) and/or delayed group(s)
            requested in the parameters.

        """

        # Temporarily remove energy filter from delayed-nu-fission-in since its
        # group structure will work in super MGXS.get_slice(...) method
        delayed_nu_fission_in = self.tallies['delayed-nu-fission-in']
        energy_filter = delayed_nu_fission_in.find_filter(openmc.EnergyFilter)
        delayed_nu_fission_in.remove_filter(energy_filter)

        # Call super class method and null out derived tallies
        slice_xs = super().get_slice(nuclides, groups, delayed_groups)
        slice_xs._rxn_rate_tally = None
        slice_xs._xs_tally = None

        # Slice energy groups if needed
        filters = []
        filter_bins = []

        if len(groups) != 0:
            energy_bins = []
            for group in groups:
                group_bounds = self.energy_groups.get_group_bounds(group)
                energy_bins.append(group_bounds)
            filter_bins.append(tuple(energy_bins))
            filters.append(openmc.EnergyoutFilter)

        if len(delayed_groups) != 0:
            filter_bins.append(tuple(delayed_groups))
            filters.append(openmc.DelayedGroupFilter)

        if filters != []:

            # Slice nu-fission-out tally along energyout filter
            delayed_nu_fission_out = slice_xs.tallies['delayed-nu-fission-out']
            tally_slice = delayed_nu_fission_out.get_slice \
                          (filters=filters, filter_bins=filter_bins)
            slice_xs._tallies['delayed-nu-fission-out'] = tally_slice

        # Add energy filter back to nu-fission-in tallies
        self.tallies['delayed-nu-fission-in'].add_filter(energy_filter)
        slice_xs._tallies['delayed-nu-fission-in'].add_filter(energy_filter)

        slice_xs.sparse = self.sparse
        return slice_xs

    def merge(self, other):
        """Merge another ChiDelayed with this one

        If results have been loaded from a statepoint, then ChiDelayed are only
        mergeable along one and only one of energy groups or nuclides.

        Parameters
        ----------
        other : openmc.mdgxs.MGXS
            MGXS to merge with this one

        Returns
        -------
        merged_mdgxs : openmc.mgxs.MDGXS
            Merged MDGXS
        """

        if not self.can_merge(other):
            raise ValueError('Unable to merge ChiDelayed')

        # Create deep copy of tally to return as merged tally
        merged_mdgxs = copy.deepcopy(self)
        merged_mdgxs._derived = True
        merged_mdgxs._rxn_rate_tally = None
        merged_mdgxs._xs_tally = None

        # Merge energy groups
        if self.energy_groups != other.energy_groups:
            merged_groups = self.energy_groups.merge(other.energy_groups)
            merged_mdgxs.energy_groups = merged_groups

        # Merge delayed groups
        if self.delayed_groups != other.delayed_groups:
            merged_mdgxs.delayed_groups = list(set(self.delayed_groups +
                                                   other.delayed_groups))

        # Merge nuclides
        if self.nuclides != other.nuclides:

            # The nuclides must be mutually exclusive
            for nuclide in self.nuclides:
                if nuclide in other.nuclides:
                    msg = 'Unable to merge Chi Delayed with shared nuclides'
                    raise ValueError(msg)

            # Concatenate lists of nuclides for the merged MGXS
            merged_mdgxs.nuclides = self.nuclides + other.nuclides

        # Merge tallies
        for tally_key in self.tallies:
            merged_tally = self.tallies[tally_key].merge\
                           (other.tallies[tally_key])
            merged_mdgxs.tallies[tally_key] = merged_tally

        return merged_mdgxs

    def get_xs(self, groups='all', subdomains='all', nuclides='all',
               xs_type='macro', order_groups='increasing',
               value='mean', delayed_groups='all', squeeze=True, **kwargs):
        """Returns an array of the delayed fission spectrum.

        This method constructs a 4D NumPy array for the requested
        multi-delayed-group cross section data for one or more
        subdomains (1st dimension), delayed groups (2nd demension),
        energy groups (3rd dimension), and nuclides (4th dimension).

        Parameters
        ----------
        groups : Iterable of Integral or 'all'
            Energy groups of interest. Defaults to 'all'.
        delayed_groups : list of int or 'all'
            Delayed groups of interest. Defaults to 'all'.
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
        value : {'mean', 'std_dev', 'rel_err'}
            A string for the type of value to return. Defaults to 'mean'.
        squeeze : bool
            A boolean representing whether to eliminate the extra dimensions
            of the multi-dimensional array to be returned. Defaults to True.

        Returns
        -------
        numpy.ndarray
            A NumPy array of the multi-group and multi-delayed-group cross
            section indexed in the order each group, subdomain and nuclide is
            listed in the parameters.

        Raises
        ------
        ValueError
            When this method is called before the multi-group cross section is
            computed from tally data.

        """

        cv.check_value('value', value, ['mean', 'std_dev', 'rel_err'])
        cv.check_value('xs_type', xs_type, ['macro', 'micro'])

        # FIXME: Unable to get microscopic xs for mesh domain because the mesh
        # cells do not know the nuclide densities in each mesh cell.
        if self.domain_type == 'mesh' and xs_type == 'micro':
            msg = 'Unable to get micro xs for mesh domain since the mesh ' \
                  'cells do not know the nuclide densities in each mesh cell.'
            raise ValueError(msg)

        filters = []
        filter_bins = []

        # Construct a collection of the domain filter bins
        if not isinstance(subdomains, str):
            cv.check_iterable_type('subdomains', subdomains, Integral,
                                   max_depth=3)
            for subdomain in subdomains:
                filters.append(_DOMAIN_TO_FILTER[self.domain_type])
                filter_bins.append((subdomain,))

        # Construct list of energy group bounds tuples for all requested groups
        if not isinstance(groups, str):
            cv.check_iterable_type('groups', groups, Integral)
            for group in groups:
                filters.append(openmc.EnergyoutFilter)
                filter_bins.append(
                    (self.energy_groups.get_group_bounds(group),))

        # Construct list of delayed group tuples for all requested groups
        if not isinstance(delayed_groups, str):
            cv.check_type('delayed groups', delayed_groups, list, int)
            for delayed_group in delayed_groups:
                filters.append(openmc.DelayedGroupFilter)
                filter_bins.append((delayed_group,))

        # If chi delayed was computed for each nuclide in the domain
        if self.by_nuclide:

            # Get the sum as the fission source weighted average chi for all
            # nuclides in the domain
            if nuclides == 'sum' or nuclides == ['sum']:

                # Retrieve the fission production tallies
                delayed_nu_fission_in = self.tallies['delayed-nu-fission-in']
                delayed_nu_fission_out = self.tallies['delayed-nu-fission-out']

                # Sum out all nuclides
                nuclides = self.get_nuclides()
                delayed_nu_fission_in = delayed_nu_fission_in.summation\
                                        (nuclides=nuclides)
                delayed_nu_fission_out = delayed_nu_fission_out.summation\
                                         (nuclides=nuclides)

                # Remove coarse energy filter to keep it out of tally arithmetic
                energy_filter = delayed_nu_fission_in.find_filter(
                    openmc.EnergyFilter)
                delayed_nu_fission_in.remove_filter(energy_filter)

                # Compute chi and store it as the xs_tally attribute so we can
                # use the generic get_xs(...) method
                xs_tally = delayed_nu_fission_out / delayed_nu_fission_in

                # Add the coarse energy filter back to the nu-fission tally
                delayed_nu_fission_in.filters.append(energy_filter)

                xs = xs_tally.get_values(filters=filters,
                                         filter_bins=filter_bins, value=value)

            # Get chi delayed for all nuclides in the domain
            elif nuclides == 'all':
                nuclides = self.get_nuclides()
                xs = self.xs_tally.get_values(filters=filters,
                                              filter_bins=filter_bins,
                                              nuclides=nuclides, value=value)

            # Get chi delayed for user-specified nuclides in the domain
            else:
                cv.check_iterable_type('nuclides', nuclides, str)
                xs = self.xs_tally.get_values(filters=filters,
                                              filter_bins=filter_bins,
                                              nuclides=nuclides, value=value)

        # If chi delayed was computed as an average of nuclides in the domain
        else:
            xs = self.xs_tally.get_values(filters=filters,
                                          filter_bins=filter_bins, value=value)

        # Eliminate the trivial score dimension
        xs = np.squeeze(xs, axis=len(xs.shape) - 1)
        xs = np.nan_to_num(xs)

        # Reshape tally data array with separate axes for domain and energy
        if groups == 'all':
            num_groups = self.num_groups
        else:
            num_groups = len(groups)

        if delayed_groups == 'all':
            num_delayed_groups = self.num_delayed_groups
        else:
            num_delayed_groups = len(delayed_groups)

        # Reshape tally data array with separate axes for domain, energy
        # groups, and accomodate the polar and azimuthal bins if needed
        num_subdomains = int(xs.shape[0] / (num_delayed_groups *
                                            num_groups * self.num_polar *
                                            self.num_azimuthal))
        if self.num_polar > 1 or self.num_azimuthal > 1:
            new_shape = (self.num_polar, self.num_azimuthal, num_subdomains,
                         num_delayed_groups, num_groups)
        else:
            new_shape = (num_subdomains, num_delayed_groups, num_groups)
        new_shape += xs.shape[1:]
        xs = np.reshape(xs, new_shape)

        # Reverse data if user requested increasing energy groups since
        # tally data is stored in order of increasing energies
        if order_groups == 'increasing':
            xs = xs[..., ::-1, :]

        if squeeze:
            # We want to squeeze out everything but the polar, azimuthal,
            # and energy group data.
            xs = self._squeeze_xs(xs)

        return xs


class DelayedNuFissionXS(MDGXS):
    r"""A fission delayed neutron production multi-group cross section.

    This class can be used for both OpenMC input generation and tally data
    post-processing to compute spatially-homogenized and energy-integrated
    multi-group fission neutron production cross sections for multi-group
    neutronics calculations. At a minimum, one needs to set the
    :attr:`DelayedNuFissionXS.energy_groups` and :attr:`DelayedNuFissionXS.domain`
    properties. Tallies for the flux and appropriate reaction rates over the
    specified domain are generated automatically via the
    :attr:`DelayedNuFissionXS.tallies` property, which can then be appended to a
    :class:`openmc.Tallies` instance.

    For post-processing, the :meth:`MGXS.load_from_statepoint` will pull in the
    necessary data to compute multi-group cross sections from a
    :class:`openmc.StatePoint` instance. The derived multi-group cross section
    can then be obtained from the :attr:`DelayedNuFissionXS.xs_tally` property.

    For a spatial domain :math:`V`, energy group :math:`[E_g,E_{g-1}]`, and
    delayed group :math:`d`, the fission delayed neutron production cross
    section is calculated as:

    .. math::

       \frac{\int_{r \in V} dr \int_{4\pi} d\Omega \int_{E_g}^{E_{g-1}} dE \;
       \nu^d \sigma_f (r, E) \psi (r, E, \Omega)}{\int_{r \in V} dr \int_{4\pi}
       d\Omega \int_{E_g}^{E_{g-1}} dE \; \psi (r, E, \Omega)}.


    Parameters
    ----------
    domain : openmc.Material or openmc.Cell or openmc.Universe or openmc.RegularMesh
        The domain for spatial homogenization
    domain_type : {'material', 'cell', 'distribcell', 'universe', 'mesh'}
        The domain type for spatial homogenization
    groups : openmc.mgxs.EnergyGroups
        The energy group structure for energy condensation
    by_nuclide : bool
        If true, computes cross sections for each nuclide in domain
    name : str, optional
        Name of the multi-group cross section. Used as a label to identify
        tallies in OpenMC 'tallies.xml' file.
    delayed_groups : list of int
        Delayed groups to filter out the xs
    num_polar : Integral, optional
        Number of equi-width polar angle bins for angle discretization;
        defaults to one bin
    num_azimuthal : Integral, optional
        Number of equi-width azimuthal angle bins for angle discretization;
        defaults to one bin

    Attributes
    ----------
    name : str, optional
        Name of the multi-group cross section
    rxn_type : str
        Reaction type (e.g., 'total', 'nu-fission', etc.)
    by_nuclide : bool
        If true, computes cross sections for each nuclide in domain
    domain : openmc.Material or openmc.Cell or openmc.Universe or openmc.RegularMesh
        Domain for spatial homogenization
    domain_type : {'material', 'cell', 'distribcell', 'universe', 'mesh'}
        Domain type for spatial homogenization
    energy_groups : openmc.mgxs.EnergyGroups
        Energy group structure for energy condensation
    delayed_groups : list of int
        Delayed groups to filter out the xs
    num_polar : Integral
        Number of equi-width polar angle bins for angle discretization
    num_azimuthal : Integral
        Number of equi-width azimuthal angle bins for angle discretization
    tally_trigger : openmc.Trigger
        An (optional) tally precision trigger given to each tally used to
        compute the cross section
    scores : list of str
        The scores in each tally used to compute the multi-group cross section
    filters : list of openmc.Filter
        The filters in each tally used to compute the multi-group cross section
    tally_keys : list of str
        The keys into the tallies dictionary for each tally used to compute
        the multi-group cross section
    estimator : {'tracklength', 'analog'}
        The tally estimator used to compute the multi-group cross section
    tallies : collections.OrderedDict
        OpenMC tallies needed to compute the multi-group cross section. The keys
        are strings listed in the :attr:`DelayedNuFissionXS.tally_keys` property
        and values are instances of :class:`openmc.Tally`.
    rxn_rate_tally : openmc.Tally
        Derived tally for the reaction rate tally used in the numerator to
        compute the multi-group cross section. This attribute is None
        unless the multi-group cross section has been computed.
    xs_tally : openmc.Tally
        Derived tally for the multi-group cross section. This attribute
        is None unless the multi-group cross section has been computed.
    num_subdomains : int
        The number of subdomains is unity for 'material', 'cell' and 'universe'
        domain types. This is equal to the number of cell instances
        for 'distribcell' domain types (it is equal to unity prior to loading
        tally data from a statepoint file).
    num_nuclides : int
        The number of nuclides for which the multi-group cross section is
        being tracked. This is unity if the by_nuclide attribute is False.
    nuclides : Iterable of str or 'sum'
        The optional user-specified nuclides for which to compute cross
        sections (e.g., 'U-238', 'O-16'). If by_nuclide is True but nuclides
        are not specified by the user, all nuclides in the spatial domain
        are included. This attribute is 'sum' if by_nuclide is false.
    sparse : bool
        Whether or not the MGXS' tallies use SciPy's LIL sparse matrix format
        for compressed data storage
    loaded_sp : bool
        Whether or not a statepoint file has been loaded with tally data
    derived : bool
        Whether or not the MGXS is merged from one or more other MGXS
    hdf5_key : str
        The key used to index multi-group cross sections in an HDF5 data store

    """

    def __init__(self, domain=None, domain_type=None, energy_groups=None,
                 delayed_groups=None, by_nuclide=False, name='',
                 num_polar=1, num_azimuthal=1):
        super().__init__(domain, domain_type, energy_groups, delayed_groups,
                         by_nuclide, name, num_polar, num_azimuthal)
        self._rxn_type = 'delayed-nu-fission'


class Beta(MDGXS):
    r"""The delayed neutron fraction.

    This class can be used for both OpenMC input generation and tally data
    post-processing to compute spatially-homogenized and energy-integrated
    multi-group and multi-delayed group cross sections for multi-group
    neutronics calculations. At a minimum, one needs to set the
    :attr:`Beta.energy_groups` and :attr:`Beta.domain` properties. Tallies for
    the flux and appropriate reaction rates over the specified domain are
    generated automatically via the :attr:`Beta.tallies` property, which can
    then be appended to a :class:`openmc.Tallies` instance.

    For post-processing, the :meth:`MGXS.load_from_statepoint` will pull in the
    necessary data to compute multi-group cross sections from a
    :class:`openmc.StatePoint` instance. The derived multi-group cross section
    can then be obtained from the :attr:`Beta.xs_tally` property.

    For a spatial domain :math:`V`, energy group :math:`[E_g,E_{g-1}]`, and
    delayed group :math:`d`, the delayed neutron fraction is calculated as:

    .. math::

       \begin{aligned}
       \langle \nu^d \sigma_f \phi \rangle &= \int_{r \in V} dr \int_{4\pi}
       d\Omega' \int_0^\infty dE' \int_0^\infty dE \; \chi(E) \nu^d
       \sigma_f (r, E') \psi(r, E', \Omega') \\
       \langle \nu \sigma_f \phi \rangle &= \int_{r \in V} dr \int_{4\pi}
       d\Omega' \int_0^\infty dE' \int_0^\infty dE \; \chi(E) \nu
       \sigma_f (r, E') \psi(r, E', \Omega') \\
       \beta_{d,g} &= \frac{\langle \nu^d \sigma_f \phi \rangle}
       {\langle \nu \sigma_f \phi \rangle}
       \end{aligned}

    NOTE: The Beta MGXS is the delayed neutron fraction computed directly from
    the nuclear data. Often the delayed neutron fraction is
    "importance-weighted" by the adjoint flux and called "beta-effective". It
    is important to make clear that this Beta is not importance-weighted.

    Parameters
    ----------
    domain : openmc.Material or openmc.Cell or openmc.Universe or openmc.RegularMesh
        The domain for spatial homogenization
    domain_type : {'material', 'cell', 'distribcell', 'universe', 'mesh'}
        The domain type for spatial homogenization
    groups : openmc.mgxs.EnergyGroups
        The energy group structure for energy condensation
    by_nuclide : bool
        If true, computes cross sections for each nuclide in domain
    name : str, optional
        Name of the multi-group cross section. Used as a label to identify
        tallies in OpenMC 'tallies.xml' file.
    delayed_groups : list of int
        Delayed groups to filter out the xs
    num_polar : Integral, optional
        Number of equi-width polar angle bins for angle discretization;
        defaults to one bin
    num_azimuthal : Integral, optional
        Number of equi-width azimuthal angle bins for angle discretization;
        defaults to one bin

    Attributes
    ----------
    name : str, optional
        Name of the multi-group cross section
    rxn_type : str
        Reaction type (e.g., 'total', 'nu-fission', etc.)
    by_nuclide : bool
        If true, computes cross sections for each nuclide in domain
    domain : openmc.Material or openmc.Cell or openmc.Universe or openmc.RegularMesh
        Domain for spatial homogenization
    domain_type : {'material', 'cell', 'distribcell', 'universe', 'mesh'}
        Domain type for spatial homogenization
    energy_groups : openmc.mgxs.EnergyGroups
        Energy group structure for energy condensation
    delayed_groups : list of int
        Delayed groups to filter out the xs
    num_polar : Integral
        Number of equi-width polar angle bins for angle discretization
    num_azimuthal : Integral
        Number of equi-width azimuthal angle bins for angle discretization
    tally_trigger : openmc.Trigger
        An (optional) tally precision trigger given to each tally used to
        compute the cross section
    scores : list of str
        The scores in each tally used to compute the multi-group cross section
    filters : list of openmc.Filter
        The filters in each tally used to compute the multi-group cross section
    tally_keys : list of str
        The keys into the tallies dictionary for each tally used to compute
        the multi-group cross section
    estimator : {'tracklength', 'analog'}
        The tally estimator used to compute the multi-group cross section
    tallies : collections.OrderedDict
        OpenMC tallies needed to compute the multi-group cross section. The keys
        are strings listed in the :attr:`Beta.tally_keys` property and
        values are instances of :class:`openmc.Tally`.
    rxn_rate_tally : openmc.Tally
        Derived tally for the reaction rate tally used in the numerator to
        compute the multi-group cross section. This attribute is None
        unless the multi-group cross section has been computed.
    xs_tally : openmc.Tally
        Derived tally for the multi-group cross section. This attribute
        is None unless the multi-group cross section has been computed.
    num_subdomains : int
        The number of subdomains is unity for 'material', 'cell' and 'universe'
        domain types. This is equal to the number of cell instances
        for 'distribcell' domain types (it is equal to unity prior to loading
        tally data from a statepoint file).
    num_nuclides : int
        The number of nuclides for which the multi-group cross section is
        being tracked. This is unity if the by_nuclide attribute is False.
    nuclides : Iterable of str or 'sum'
        The optional user-specified nuclides for which to compute cross
        sections (e.g., 'U-238', 'O-16'). If by_nuclide is True but nuclides
        are not specified by the user, all nuclides in the spatial domain
        are included. This attribute is 'sum' if by_nuclide is false.
    sparse : bool
        Whether or not the MGXS' tallies use SciPy's LIL sparse matrix format
        for compressed data storage
    loaded_sp : bool
        Whether or not a statepoint file has been loaded with tally data
    derived : bool
        Whether or not the MGXS is merged from one or more other MGXS
    hdf5_key : str
        The key used to index multi-group cross sections in an HDF5 data store

    """

    def __init__(self, domain=None, domain_type=None, energy_groups=None,
                 delayed_groups=None, by_nuclide=False, name='',
                 num_polar=1, num_azimuthal=1):
        super().__init__(domain, domain_type, energy_groups, delayed_groups,
                         by_nuclide, name, num_polar, num_azimuthal)
        self._rxn_type = 'beta'

    @property
    def scores(self):
        return ['nu-fission', 'delayed-nu-fission']

    @property
    def tally_keys(self):
        return ['nu-fission', 'delayed-nu-fission']

    @property
    def rxn_rate_tally(self):
        if self._rxn_rate_tally is None:
            self._rxn_rate_tally = self.tallies['delayed-nu-fission']
            self._rxn_rate_tally.sparse = self.sparse
        return self._rxn_rate_tally

    @property
    def xs_tally(self):

        if self._xs_tally is None:
            nu_fission = self.tallies['nu-fission']

            # Compute beta
            self._xs_tally = self.rxn_rate_tally / nu_fission
            super()._compute_xs()

        return self._xs_tally

    def get_homogenized_mgxs(self, other_mgxs):
        """Construct a homogenized MGXS with other MGXS objects.

        This method constructs a new MGXS object that is the flux-weighted
        combination of two MGXS objects. It is equivalent to what one would
        obtain if the tally spatial domain were designed to encompass the
        individual domains for both MGXS objects.

        Parameters
        ----------
        other_mgxs : openmc.mgxs.MGXS or Iterable of openmc.mgxs.MGXS
            The MGXS to homogenize with this one.

        Returns
        -------
        openmc.mgxs.MGXS
            A new homogenized MGXS

        Raises
        ------
        ValueError
            If the other_mgxs is of a different type.

        """

        return self._get_homogenized_mgxs(other_mgxs, 'nu-fission')


class DecayRate(MDGXS):
    r"""The decay rate for delayed neutron precursors.

    This class can be used for both OpenMC input generation and tally data
    post-processing to compute spatially-homogenized and energy-integrated
    multi-group and multi-delayed group cross sections for multi-group
    neutronics calculations. At a minimum, one needs to set the
    :attr:`DecayRate.energy_groups` and :attr:`DecayRate.domain` properties.
    Tallies for the flux and appropriate reaction rates over the specified
    domain are generated automatically via the :attr:`DecayRate.tallies`
    property, which can then be appended to a :class:`openmc.Tallies` instance.

    For post-processing, the :meth:`MGXS.load_from_statepoint` will pull in the
    necessary data to compute multi-group cross sections from a
    :class:`openmc.StatePoint` instance. The derived multi-group cross section
    can then be obtained from the :attr:`DecayRate.xs_tally` property.

    For a spatial domain :math:`V`, energy group :math:`[E_g,E_{g-1}]`, and
    delayed group :math:`d`, the decay rate is calculated as:

    .. math::

       \begin{aligned}
       \langle \lambda_d \nu^d \sigma_f \phi \rangle &= \int_{r \in V} dr
       \int_{4\pi} d\Omega' \int_0^\infty dE' \int_0^\infty dE \; \lambda_d \nu^d
       \sigma_f (r, E') \psi(r, E', \Omega') \\
       \langle \nu^d \sigma_f \phi \rangle &= \int_{r \in V} dr \int_{4\pi}
       d\Omega' \int_0^\infty dE' \int_0^\infty dE \; \chi(E) \nu^d
       \sigma_f (r, E') \psi(r, E', \Omega') \\
       \lambda_d &= \frac{\langle \lambda_d \nu^d \sigma_f \phi \rangle}
       {\langle \nu^d \sigma_f \phi \rangle}
       \end{aligned}

    Parameters
    ----------
    domain : openmc.Material or openmc.Cell or openmc.Universe or openmc.RegularMesh
        The domain for spatial homogenization
    domain_type : {'material', 'cell', 'distribcell', 'universe', 'mesh'}
        The domain type for spatial homogenization
    groups : openmc.mgxs.EnergyGroups
        The energy group structure for energy condensation
    by_nuclide : bool
        If true, computes cross sections for each nuclide in domain
    name : str, optional
        Name of the multi-group cross section. Used as a label to identify
        tallies in OpenMC 'tallies.xml' file.
    delayed_groups : list of int
        Delayed groups to filter out the xs
    num_polar : Integral, optional
        Number of equi-width polar angle bins for angle discretization;
        defaults to one bin
    num_azimuthal : Integral, optional
        Number of equi-width azimuthal angle bins for angle discretization;
        defaults to one bin

    Attributes
    ----------
    name : str, optional
        Name of the multi-group cross section
    rxn_type : str
        Reaction type (e.g., 'total', 'nu-fission', etc.)
    by_nuclide : bool
        If true, computes cross sections for each nuclide in domain
    domain : openmc.Material or openmc.Cell or openmc.Universe or openmc.RegularMesh
        Domain for spatial homogenization
    domain_type : {'material', 'cell', 'distribcell', 'universe', 'mesh'}
        Domain type for spatial homogenization
    energy_groups : openmc.mgxs.EnergyGroups
        Energy group structure for energy condensation
    delayed_groups : list of int
        Delayed groups to filter out the xs
    num_polar : Integral
        Number of equi-width polar angle bins for angle discretization
    num_azimuthal : Integral
        Number of equi-width azimuthal angle bins for angle discretization
    tally_trigger : openmc.Trigger
        An (optional) tally precision trigger given to each tally used to
        compute the cross section
    scores : list of str
        The scores in each tally used to compute the multi-group cross section
    filters : list of openmc.Filter
        The filters in each tally used to compute the multi-group cross section
    tally_keys : list of str
        The keys into the tallies dictionary for each tally used to compute
        the multi-group cross section
    estimator : {'tracklength', 'analog'}
        The tally estimator used to compute the multi-group cross section
    tallies : collections.OrderedDict
        OpenMC tallies needed to compute the multi-group cross section. The keys
        are strings listed in the :attr:`DecayRate.tally_keys` property and
        values are instances of :class:`openmc.Tally`.
    rxn_rate_tally : openmc.Tally
        Derived tally for the reaction rate tally used in the numerator to
        compute the multi-group cross section. This attribute is None
        unless the multi-group cross section has been computed.
    xs_tally : openmc.Tally
        Derived tally for the multi-group cross section. This attribute
        is None unless the multi-group cross section has been computed.
    num_subdomains : int
        The number of subdomains is unity for 'material', 'cell' and 'universe'
        domain types. This is equal to the number of cell instances
        for 'distribcell' domain types (it is equal to unity prior to loading
        tally data from a statepoint file).
    num_nuclides : int
        The number of nuclides for which the multi-group cross section is
        being tracked. This is unity if the by_nuclide attribute is False.
    nuclides : Iterable of str or 'sum'
        The optional user-specified nuclides for which to compute cross
        sections (e.g., 'U-238', 'O-16'). If by_nuclide is True but nuclides
        are not specified by the user, all nuclides in the spatial domain
        are included. This attribute is 'sum' if by_nuclide is false.
    sparse : bool
        Whether or not the MGXS' tallies use SciPy's LIL sparse matrix format
        for compressed data storage
    loaded_sp : bool
        Whether or not a statepoint file has been loaded with tally data
    derived : bool
        Whether or not the MGXS is merged from one or more other MGXS
    hdf5_key : str
        The key used to index multi-group cross sections in an HDF5 data store

    """

    def __init__(self, domain=None, domain_type=None, energy_groups=None,
                 delayed_groups=None, by_nuclide=False, name='',
                 num_polar=1, num_azimuthal=1):
        super().__init__(domain, domain_type, energy_groups, delayed_groups,
                         by_nuclide, name, num_polar, num_azimuthal)
        self._rxn_type = 'decay-rate'

    @property
    def scores(self):
        return ['delayed-nu-fission', 'decay-rate']

    @property
    def tally_keys(self):
        return ['delayed-nu-fission', 'decay-rate']

    @property
    def filters(self):

        # Create the non-domain specific Filters for the Tallies
        group_edges = self.energy_groups.group_edges
        energy_filter = openmc.EnergyFilter(group_edges)

        if self.delayed_groups is not None:
            delayed_filter = openmc.DelayedGroupFilter(self.delayed_groups)
            filters = [[delayed_filter, energy_filter], [delayed_filter,
                                                         energy_filter]]
        else:
            filters = [[energy_filter], [energy_filter]]

        return self._add_angle_filters(filters)

    @property
    def xs_tally(self):

        if self._xs_tally is None:
            delayed_nu_fission = self.tallies['delayed-nu-fission']

            # Compute the decay rate
            self._xs_tally = self.rxn_rate_tally / delayed_nu_fission
            super()._compute_xs()

        return self._xs_tally

    def get_homogenized_mgxs(self, other_mgxs):
        """Construct a homogenized MGXS with other MGXS objects.

        This method constructs a new MGXS object that is the flux-weighted
        combination of two MGXS objects. It is equivalent to what one would
        obtain if the tally spatial domain were designed to encompass the
        individual domains for both MGXS objects.

        Parameters
        ----------
        other_mgxs : openmc.mgxs.MGXS or Iterable of openmc.mgxs.MGXS
            The MGXS to homogenize with this one.

        Returns
        -------
        openmc.mgxs.MGXS
            A new homogenized MGXS

        Raises
        ------
        ValueError
            If the other_mgxs is of a different type.

        """

        return self._get_homogenized_mgxs(other_mgxs, 'delayed-nu-fission')


class MatrixMDGXS(MDGXS):
    """An abstract multi-delayed-group cross section for some energy group and
    delayed group structure within some spatial domain. This class is
    specifically intended for cross sections which depend on both the incoming
    and outgoing energy groups and are therefore represented by matrices.
    An example of this is the delayed-nu-fission matrix.

    This class can be used for both OpenMC input generation and tally data
    post-processing to compute spatially-homogenized and energy-integrated
    multi-group and multi-delayed-group cross sections for downstream neutronics
    calculations.

    NOTE: Users should instantiate the subclasses of this abstract class.

    Parameters
    ----------
    domain : openmc.Material or openmc.Cell or openmc.Universe or openmc.RegularMesh
        The domain for spatial homogenization
    domain_type : {'material', 'cell', 'distribcell', 'universe', 'mesh'}
        The domain type for spatial homogenization
    energy_groups : openmc.mgxs.EnergyGroups
        The energy group structure for energy condensation
    by_nuclide : bool
        If true, computes cross sections for each nuclide in domain
    name : str, optional
        Name of the multi-group cross section. Used as a label to identify
        tallies in OpenMC 'tallies.xml' file.
    delayed_groups : list of int
        Delayed groups to filter out the xs
    num_polar : Integral, optional
        Number of equi-width polar angle bins for angle discretization;
        defaults to one bin
    num_azimuthal : Integral, optional
        Number of equi-width azimuthal angle bins for angle discretization;
        defaults to one bin

    Attributes
    ----------
    name : str, optional
        Name of the multi-group cross section
    rxn_type : str
        Reaction type (e.g., 'total', 'nu-fission', etc.)
    by_nuclide : bool
        If true, computes cross sections for each nuclide in domain
    domain : openmc.Material or openmc.Cell or openmc.Universe or openmc.RegularMesh
        Domain for spatial homogenization
    domain_type : {'material', 'cell', 'distribcell', 'universe', 'mesh'}
        Domain type for spatial homogenization
    energy_groups : openmc.mgxs.EnergyGroups
        Energy group structure for energy condensation
    delayed_groups : list of int
        Delayed groups to filter out the xs
    num_polar : Integral
        Number of equi-width polar angle bins for angle discretization
    num_azimuthal : Integral
        Number of equi-width azimuthal angle bins for angle discretization
    tally_trigger : openmc.Trigger
        An (optional) tally precision trigger given to each tally used to
        compute the cross section
    scores : list of str
        The scores in each tally used to compute the multi-group cross section
    filters : list of openmc.Filter
        The filters in each tally used to compute the multi-group cross section
    tally_keys : list of str
        The keys into the tallies dictionary for each tally used to compute
        the multi-group cross section
    estimator : {'tracklength', 'collision', 'analog'}
        The tally estimator used to compute the multi-group cross section
    tallies : collections.OrderedDict
        OpenMC tallies needed to compute the multi-group cross section
    rxn_rate_tally : openmc.Tally
        Derived tally for the reaction rate tally used in the numerator to
        compute the multi-group cross section. This attribute is None
        unless the multi-group cross section has been computed.
    xs_tally : openmc.Tally
        Derived tally for the multi-group cross section. This attribute
        is None unless the multi-group cross section has been computed.
    num_subdomains : int
        The number of subdomains is unity for 'material', 'cell' and 'universe'
        domain types. This is equal to the number of cell instances
        for 'distribcell' domain types (it is equal to unity prior to loading
        tally data from a statepoint file) and the number of mesh cells for
        'mesh' domain types.
    num_nuclides : int
        The number of nuclides for which the multi-group cross section is
        being tracked. This is unity if the by_nuclide attribute is False.
    nuclides : Iterable of str or 'sum'
        The optional user-specified nuclides for which to compute cross
        sections (e.g., 'U238', 'O16'). If by_nuclide is True but nuclides
        are not specified by the user, all nuclides in the spatial domain
        are included. This attribute is 'sum' if by_nuclide is false.
    sparse : bool
        Whether or not the MGXS' tallies use SciPy's LIL sparse matrix format
        for compressed data storage
    loaded_sp : bool
        Whether or not a statepoint file has been loaded with tally data
    derived : bool
        Whether or not the MGXS is merged from one or more other MGXS
    hdf5_key : str
        The key used to index multi-group cross sections in an HDF5 data store

    """

    @property
    def _dont_squeeze(self):
        """Create a tuple of axes which should not be removed during the get_xs
        process
        """
        if self.num_polar > 1 or self.num_azimuthal > 1:
            return (0, 1, 3, 4, 5)
        else:
            return (1, 2, 3)

    @property
    def filters(self):
        # Create the non-domain specific Filters for the Tallies
        group_edges = self.energy_groups.group_edges
        energy = openmc.EnergyFilter(group_edges)
        energyout = openmc.EnergyoutFilter(group_edges)

        if self.delayed_groups is not None:
            delayed = openmc.DelayedGroupFilter(self.delayed_groups)
            filters = [[energy], [delayed, energy, energyout]]
        else:
            filters = [[energy], [energy, energyout]]

        return self._add_angle_filters(filters)

    def get_xs(self, in_groups='all', out_groups='all',
               subdomains='all', nuclides='all',
               xs_type='macro', order_groups='increasing',
               row_column='inout', value='mean', delayed_groups='all',
               squeeze=True, **kwargs):
        """Returns an array of multi-group cross sections.

        This method constructs a 4D NumPy array for the requested
        multi-group cross section data for one or more subdomains
        (1st dimension), delayed groups (2nd dimension), energy groups in
        (3rd dimension), energy groups out (4th dimension), and nuclides
        (5th dimension).

        Parameters
        ----------
        in_groups : Iterable of Integral or 'all'
            Incoming energy groups of interest. Defaults to 'all'.
        out_groups : Iterable of Integral or 'all'
            Outgoing energy groups of interest. Defaults to 'all'.
        subdomains : Iterable of Integral or 'all'
            Subdomain IDs of interest. Defaults to 'all'.
        nuclides : Iterable of str or 'all' or 'sum'
            A list of nuclide name strings (e.g., ['U235', 'U238']). The
            special string 'all' will return the cross sections for all
            nuclides in the spatial domain. The special string 'sum' will
            return the cross section summed over all nuclides. Defaults to
            'all'.
        xs_type: {'macro', 'micro'}
            Return the macro or micro cross section in units of cm^-1 or barns.
            Defaults to 'macro'.
        order_groups: {'increasing', 'decreasing'}
            Return the cross section indexed according to increasing or
            decreasing energy groups (decreasing or increasing energies).
            Defaults to 'increasing'.
        row_column: {'inout', 'outin'}
            Return the cross section indexed first by incoming group and
            second by outgoing group ('inout'), or vice versa ('outin').
            Defaults to 'inout'.
        value : {'mean', 'std_dev', 'rel_err'}
            A string for the type of value to return. Defaults to 'mean'.
        delayed_groups : list of int or 'all'
            Delayed groups of interest. Defaults to 'all'.
        squeeze : bool
            A boolean representing whether to eliminate the extra dimensions
            of the multi-dimensional array to be returned. Defaults to True.

        Returns
        -------
        numpy.ndarray
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

        # FIXME: Unable to get microscopic xs for mesh domain because the mesh
        # cells do not know the nuclide densities in each mesh cell.
        if self.domain_type == 'mesh' and xs_type == 'micro':
            msg = 'Unable to get micro xs for mesh domain since the mesh ' \
                  'cells do not know the nuclide densities in each mesh cell.'
            raise ValueError(msg)

        filters = []
        filter_bins = []

        # Construct a collection of the domain filter bins
        if not isinstance(subdomains, str):
            cv.check_iterable_type('subdomains', subdomains, Integral,
                                   max_depth=3)
            for subdomain in subdomains:
                filters.append(_DOMAIN_TO_FILTER[self.domain_type])
                filter_bins.append((subdomain,))

        # Construct list of energy group bounds tuples for all requested groups
        if not isinstance(in_groups, str):
            cv.check_iterable_type('groups', in_groups, Integral)
            for group in in_groups:
                filters.append(openmc.EnergyFilter)
                filter_bins.append((
                    self.energy_groups.get_group_bounds(group),))

        # Construct list of energy group bounds tuples for all requested groups
        if not isinstance(out_groups, str):
            cv.check_iterable_type('groups', out_groups, Integral)
            for group in out_groups:
                filters.append(openmc.EnergyoutFilter)
                filter_bins.append((
                    self.energy_groups.get_group_bounds(group),))

        # Construct list of delayed group tuples for all requested groups
        if not isinstance(delayed_groups, str):
            cv.check_type('delayed groups', delayed_groups, list, int)
            for delayed_group in delayed_groups:
                filters.append(openmc.DelayedGroupFilter)
                filter_bins.append((delayed_group,))

        # Construct a collection of the nuclides to retrieve from the xs tally
        if self.by_nuclide:
            if nuclides == 'all' or nuclides == 'sum' or nuclides == ['sum']:
                query_nuclides = self.get_nuclides()
            else:
                query_nuclides = nuclides
        else:
            query_nuclides = ['total']

        # Use tally summation if user requested the sum for all nuclides
        if nuclides == 'sum' or nuclides == ['sum']:
            xs_tally = self.xs_tally.summation(nuclides=query_nuclides)
            xs = xs_tally.get_values(filters=filters, filter_bins=filter_bins,
                                     value=value)
        else:
            xs = self.xs_tally.get_values(filters=filters,
                                          filter_bins=filter_bins,
                                          nuclides=query_nuclides, value=value)

        # Divide by atom number densities for microscopic cross sections
        if xs_type == 'micro':
            if self.by_nuclide:
                densities = self.get_nuclide_densities(nuclides)
            else:
                densities = self.get_nuclide_densities('sum')
            if value == 'mean' or value == 'std_dev':
                xs /= densities[np.newaxis, :, np.newaxis]

        # Eliminate the trivial score dimension
        xs = np.squeeze(xs, axis=len(xs.shape) - 1)
        xs = np.nan_to_num(xs)

        if in_groups == 'all':
            num_in_groups = self.num_groups
        else:
            num_in_groups = len(in_groups)

        if out_groups == 'all':
            num_out_groups = self.num_groups
        else:
            num_out_groups = len(out_groups)

        if delayed_groups == 'all':
            num_delayed_groups = self.num_delayed_groups
        else:
            num_delayed_groups = len(delayed_groups)

        # Reshape tally data array with separate axes for domain and energy
        # Accomodate the polar and azimuthal bins if needed
        num_subdomains = int(xs.shape[0] / (num_delayed_groups *
                                            num_in_groups * num_out_groups *
                                            self.num_polar *
                                            self.num_azimuthal))
        if self.num_polar > 1 or self.num_azimuthal > 1:
            new_shape = (self.num_polar, self.num_azimuthal, num_subdomains,
                         num_delayed_groups, num_in_groups, num_out_groups)
            new_shape += xs.shape[1:]
            xs = np.reshape(xs, new_shape)

            # Transpose the matrix if requested by user
            if row_column == 'outin':
                xs = np.swapaxes(xs, 4, 5)
        else:
            new_shape = (num_subdomains, num_delayed_groups, num_in_groups,
                         num_out_groups)
            new_shape += xs.shape[1:]
            xs = np.reshape(xs, new_shape)

            # Transpose the matrix if requested by user
            if row_column == 'outin':
                xs = np.swapaxes(xs, 2, 3)

        # Reverse data if user requested increasing energy groups since
        # tally data is stored in order of increasing energies
        if order_groups == 'increasing':
            xs = xs[..., ::-1, ::-1, :]

        if squeeze:
            # We want to squeeze out everything but the polar, azimuthal,
            # and in/out energy group data.
            xs = self._squeeze_xs(xs)

        return xs

    def get_slice(self, nuclides=[], in_groups=[], out_groups=[],
                  delayed_groups=[]):
        """Build a sliced MatrixMDGXS object for the specified nuclides and
        energy groups.

        This method constructs a new MdGXS to encapsulate a subset of the data
        represented by this MdGXS. The subset of data to include in the tally
        slice is determined by the nuclides, energy groups, and delayed groups
        specified in the input parameters.

        Parameters
        ----------
        nuclides : list of str
            A list of nuclide name strings
            (e.g., ['U235', 'U238']; default is [])
        in_groups : list of int
            A list of incoming energy group indices starting at 1 for the high
            energies (e.g., [1, 2, 3]; default is [])
        out_groups : list of int
            A list of outgoing energy group indices starting at 1 for the high
            energies (e.g., [1, 2, 3]; default is [])
        delayed_groups : list of int
            A list of delayed group indices
            (e.g., [1, 2, 3]; default is [])

        Returns
        -------
        openmc.mgxs.MatrixMDGXS
            A new MatrixMDGXS object which encapsulates the subset of data
            requested for the nuclide(s) and/or energy group(s) requested in
            the parameters.

        """

        # Call super class method and null out derived tallies
        slice_xs = super().get_slice(nuclides, in_groups, delayed_groups)
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
                if tally.contains_filter(openmc.EnergyoutFilter):
                    tally_slice = tally.get_slice(
                        filters=[openmc.EnergyoutFilter],
                        filter_bins=filter_bins)
                    slice_xs.tallies[tally_type] = tally_slice

        slice_xs.sparse = self.sparse
        return slice_xs

    def print_xs(self, subdomains='all', nuclides='all', xs_type='macro'):
        """Prints a string representation for the multi-group cross section.

        Parameters
        ----------
        subdomains : Iterable of Integral or 'all'
            The subdomain IDs of the cross sections to include in the report.
            Defaults to 'all'.
        nuclides : Iterable of str or 'all' or 'sum'
            The nuclides of the cross-sections to include in the report. This
            may be a list of nuclide name strings (e.g., ['U235', 'U238']).
            The special string 'all' will report the cross sections for all
            nuclides in the spatial domain. The special string 'sum' will
            report the cross sections summed over all nuclides. Defaults to
            'all'.
        xs_type: {'macro', 'micro'}
            Return the macro or micro cross section in units of cm^-1 or barns.
            Defaults to 'macro'.

        """

        # Construct a collection of the subdomains to report
        if not isinstance(subdomains, str):
            cv.check_iterable_type('subdomains', subdomains, Integral)
        elif self.domain_type == 'distribcell':
            subdomains = np.arange(self.num_subdomains, dtype=np.int)
        elif self.domain_type == 'mesh':
            xyz = [range(1, x + 1) for x in self.domain.dimension]
            subdomains = list(itertools.product(*xyz))
        else:
            subdomains = [self.domain.id]

        # Construct a collection of the nuclides to report
        if self.by_nuclide:
            if nuclides == 'all':
                nuclides = self.get_nuclides()
            if nuclides == 'sum':
                nuclides = ['sum']
            else:
                cv.check_iterable_type('nuclides', nuclides, str)
        else:
            nuclides = ['sum']

        cv.check_value('xs_type', xs_type, ['macro', 'micro'])

        # Build header for string with type and domain info
        string = 'Multi-Delayed-Group XS\n'
        string += '{0: <16}=\t{1}\n'.format('\tReaction Type', self.rxn_type)
        string += '{0: <16}=\t{1}\n'.format('\tDomain Type', self.domain_type)
        string += '{0: <16}=\t{1}\n'.format('\tDomain ID', self.domain.id)

        # Generate the header for an individual XS
        xs_header = '\tCross Sections [{0}]:'.format(self.get_units(xs_type))

        # If cross section data has not been computed, only print string header
        if self.tallies is None:
            print(string)
            return

        string += '{0: <16}\n'.format('\tEnergy Groups:')
        template = '{0: <12}Group {1} [{2: <10} - {3: <10}eV]\n'

        # Loop over energy groups ranges
        for group in range(1, self.num_groups + 1):
            bounds = self.energy_groups.get_group_bounds(group)
            string += template.format('', group, bounds[0], bounds[1])

        # Set polar and azimuthal bins if necessary
        if self.num_polar > 1 or self.num_azimuthal > 1:
            pol_bins = np.linspace(0., np.pi, num=self.num_polar + 1,
                                   endpoint=True)
            azi_bins = np.linspace(-np.pi, np.pi, num=self.num_azimuthal + 1,
                                   endpoint=True)

        # Loop over all subdomains
        for subdomain in subdomains:

            if self.domain_type == 'distribcell':
                string += '{: <16}=\t{}\n'.format('\tSubdomain', subdomain)

            # Loop over all Nuclides
            for nuclide in nuclides:

                # Build header for nuclide type
                if xs_type != 'sum':
                    string += '{: <16}=\t{}\n'.format('\tNuclide', nuclide)

                # Build header for cross section type
                string += '{: <16}\n'.format(xs_header)

                if self.delayed_groups is not None:

                    for delayed_group in self.delayed_groups:

                        template = '{0: <12}Delayed Group {1}:\t'
                        string += template.format('', delayed_group)
                        string += '\n'

                        template = '{0: <12}Group {1} -> Group {2}:\t\t'

                        average_xs = self.get_xs(nuclides=[nuclide],
                                                 subdomains=[subdomain],
                                                 xs_type=xs_type, value='mean',
                                                 delayed_groups=[delayed_group])
                        rel_err_xs = self.get_xs(nuclides=[nuclide],
                                                 subdomains=[subdomain],
                                                 xs_type=xs_type,
                                                 value='rel_err',
                                                 delayed_groups=[delayed_group])
                        rel_err_xs = rel_err_xs * 100.

                        if self.num_polar > 1 or self.num_azimuthal > 1:
                            # Loop over polar, azi, and in/out group ranges
                            for pol in range(len(pol_bins) - 1):
                                pol_low, pol_high = pol_bins[pol: pol + 2]
                                for azi in range(len(azi_bins) - 1):
                                    azi_low, azi_high = azi_bins[azi: azi + 2]
                                    string += '\t\tPolar Angle: [{0:5f} - {1:5f}]'.format(
                                        pol_low, pol_high) + \
                                        '\tAzimuthal Angle: [{0:5f} - {1:5f}]'.format(
                                        azi_low, azi_high) + '\n'
                                    for in_group in range(1, self.num_groups + 1):
                                        for out_group in range(1, self.num_groups + 1):
                                            string += '\t' + template.format(
                                                '', in_group, out_group)
                                            string += '{0:.2e} +/- {1:.2e}%'.format(
                                                average_xs[pol, azi, in_group - 1,
                                                           out_group - 1],
                                                rel_err_xs[pol, azi, in_group - 1,
                                                           out_group - 1])
                                            string += '\n'
                                        string += '\n'
                                    string += '\n'
                        else:
                            # Loop over incoming/outgoing energy groups ranges
                            for in_group in range(1, self.num_groups + 1):
                                for out_group in range(1, self.num_groups + 1):
                                    string += template.format(
                                        '', in_group, out_group)
                                    string += '{:.2e} +/- {:.2e}%'.format(
                                        average_xs[in_group-1, out_group-1],
                                        rel_err_xs[in_group-1, out_group-1])
                                    string += '\n'
                                string += '\n'
                        string += '\n'
                else:

                    template = '{0: <12}Group {1} -> Group {2}:\t\t'

                    average_xs = self.get_xs(nuclides=[nuclide],
                                             subdomains=[subdomain],
                                             xs_type=xs_type, value='mean')
                    rel_err_xs = self.get_xs(nuclides=[nuclide],
                                             subdomains=[subdomain],
                                             xs_type=xs_type, value='rel_err')
                    rel_err_xs = rel_err_xs * 100.

                    if self.num_polar > 1 or self.num_azimuthal > 1:
                        # Loop over polar, azi, and in/out energy group ranges
                        for pol in range(len(pol_bins) - 1):
                            pol_low, pol_high = pol_bins[pol: pol + 2]
                            for azi in range(len(azi_bins) - 1):
                                azi_low, azi_high = azi_bins[azi: azi + 2]
                                string += '\t\tPolar Angle: [{0:5f} - {1:5f}]'.format(
                                    pol_low, pol_high) + \
                                    '\tAzimuthal Angle: [{0:5f} - {1:5f}]'.format(
                                    azi_low, azi_high) + '\n'
                                for in_group in range(1, self.num_groups + 1):
                                    for out_group in range(1, self.num_groups + 1):
                                        string += '\t' + template.format(
                                            '', in_group, out_group)
                                        string += '{0:.2e} +/- {1:.2e}%'.format(
                                            average_xs[pol, azi, in_group - 1,
                                                       out_group - 1],
                                            rel_err_xs[pol, azi, in_group - 1,
                                                       out_group - 1])
                                        string += '\n'
                                    string += '\n'
                                string += '\n'
                    else:
                        # Loop over incoming/outgoing energy groups ranges
                        for in_group in range(1, self.num_groups + 1):
                            for out_group in range(1, self.num_groups + 1):
                                string += template.format('', in_group,
                                                          out_group)
                                string += '{0:.2e} +/- {1:.2e}%'.format(
                                    average_xs[in_group - 1, out_group - 1],
                                    rel_err_xs[in_group - 1, out_group - 1])
                                string += '\n'
                            string += '\n'
                    string += '\n'
                string += '\n'
            string += '\n'

        print(string)


class DelayedNuFissionMatrixXS(MatrixMDGXS):
    r"""A fission delayed neutron production matrix multi-group cross section.

    This class can be used for both OpenMC input generation and tally data
    post-processing to compute spatially-homogenized and energy-integrated
    multi-group fission neutron production cross sections for multi-group
    neutronics calculations. At a minimum, one needs to set the
    :attr:`DelayedNuFissionMatrixXS.energy_groups` and
    :attr:`DelayedNuFissionMatrixXS.domain` properties. Tallies for the flux and
    appropriate reaction rates over the specified domain are generated
    automatically via the :attr:`DelayedNuFissionMatrixXS.tallies` property,
    which can then be appended to a :class:`openmc.Tallies` instance.

    For post-processing, the :meth:`MGXS.load_from_statepoint` will pull in the
    necessary data to compute multi-group cross sections from a
    :class:`openmc.StatePoint` instance. The derived multi-group cross section
    can then be obtained from the :attr:`DelayedNuFissionMatrixXS.xs_tally`
    property.

    For a spatial domain :math:`V`, energy group :math:`[E_g,E_{g-1}]`, and
    delayed group :math:`d`, the fission delayed neutron production cross
    section is calculated as:

    .. math::

       \begin{aligned}
       \langle \nu\sigma_{f,g'\rightarrow g} \phi \rangle &= \int_{r \in V} dr
       \int_{4\pi} d\Omega' \int_{E_{g'}}^{E_{g'-1}} dE' \int_{E_g}^{E_{g-1}} dE
       \; \chi(E) \nu\sigma_f^d (r, E') \psi(r, E', \Omega')\\
       \langle \phi \rangle &= \int_{r \in V} dr \int_{4\pi} d\Omega
       \int_{E_g}^{E_{g-1}} dE \; \psi (r, E, \Omega) \\
       \nu\sigma_{f,g'\rightarrow g} &= \frac{\langle \nu\sigma_{f,g'\rightarrow
       g}^d \phi \rangle}{\langle \phi \rangle}
       \end{aligned}


    Parameters
    ----------
    domain : openmc.Material or openmc.Cell or openmc.Universe or openmc.RegularMesh
        The domain for spatial homogenization
    domain_type : {'material', 'cell', 'distribcell', 'universe', 'mesh'}
        The domain type for spatial homogenization
    groups : openmc.mgxs.EnergyGroups
        The energy group structure for energy condensation
    by_nuclide : bool
        If true, computes cross sections for each nuclide in domain
    name : str, optional
        Name of the multi-group cross section. Used as a label to identify
        tallies in OpenMC 'tallies.xml' file.
    delayed_groups : list of int
        Delayed groups to filter out the xs
    num_polar : Integral, optional
        Number of equi-width polar angle bins for angle discretization;
        defaults to one bin
    num_azimuthal : Integral, optional
        Number of equi-width azimuthal angle bins for angle discretization;
        defaults to one bin

    Attributes
    ----------
    name : str, optional
        Name of the multi-group cross section
    rxn_type : str
        Reaction type (e.g., 'total', 'nu-fission', etc.)
    by_nuclide : bool
        If true, computes cross sections for each nuclide in domain
    domain : openmc.Material or openmc.Cell or openmc.Universe or openmc.RegularMesh
        Domain for spatial homogenization
    domain_type : {'material', 'cell', 'distribcell', 'universe', 'mesh'}
        Domain type for spatial homogenization
    energy_groups : openmc.mgxs.EnergyGroups
        Energy group structure for energy condensation
    delayed_groups : list of int
        Delayed groups to filter out the xs
    num_polar : Integral
        Number of equi-width polar angle bins for angle discretization
    num_azimuthal : Integral
        Number of equi-width azimuthal angle bins for angle discretization
    tally_trigger : openmc.Trigger
        An (optional) tally precision trigger given to each tally used to
        compute the cross section
    scores : list of str
        The scores in each tally used to compute the multi-group cross section
    filters : list of openmc.Filter
        The filters in each tally used to compute the multi-group cross section
    tally_keys : list of str
        The keys into the tallies dictionary for each tally used to compute
        the multi-group cross section
    estimator : 'analog'
        The tally estimator used to compute the multi-group cross section
    tallies : collections.OrderedDict
        OpenMC tallies needed to compute the multi-group cross section. The keys
        are strings listed in the :attr:`DelayedNuFissionXS.tally_keys` property
        and values are instances of :class:`openmc.Tally`.
    rxn_rate_tally : openmc.Tally
        Derived tally for the reaction rate tally used in the numerator to
        compute the multi-group cross section. This attribute is None
        unless the multi-group cross section has been computed.
    xs_tally : openmc.Tally
        Derived tally for the multi-group cross section. This attribute
        is None unless the multi-group cross section has been computed.
    num_subdomains : int
        The number of subdomains is unity for 'material', 'cell' and 'universe'
        domain types. This is equal to the number of cell instances
        for 'distribcell' domain types (it is equal to unity prior to loading
        tally data from a statepoint file).
    num_nuclides : int
        The number of nuclides for which the multi-group cross section is
        being tracked. This is unity if the by_nuclide attribute is False.
    nuclides : Iterable of str or 'sum'
        The optional user-specified nuclides for which to compute cross
        sections (e.g., 'U-238', 'O-16'). If by_nuclide is True but nuclides
        are not specified by the user, all nuclides in the spatial domain
        are included. This attribute is 'sum' if by_nuclide is false.
    sparse : bool
        Whether or not the MGXS' tallies use SciPy's LIL sparse matrix format
        for compressed data storage
    loaded_sp : bool
        Whether or not a statepoint file has been loaded with tally data
    derived : bool
        Whether or not the MGXS is merged from one or more other MGXS
    hdf5_key : str
        The key used to index multi-group cross sections in an HDF5 data store

    """

    def __init__(self, domain=None, domain_type=None, energy_groups=None,
                 delayed_groups=None, by_nuclide=False, name='',
                 num_polar=1, num_azimuthal=1):
        super().__init__(domain, domain_type, energy_groups, delayed_groups,
                         by_nuclide, name, num_polar, num_azimuthal)
        self._rxn_type = 'delayed-nu-fission'
        self._hdf5_key = 'delayed-nu-fission matrix'
        self._estimator = 'analog'
        self._valid_estimators = ['analog']
