from __future__ import division

from collections import Iterable, OrderedDict
from numbers import Integral
import warnings
import os
import sys
import copy
import abc

import numpy as np

from mgxs import MGXS, MGXS_TYPES, DOMAIN_TYPES, _DOMAINS
from openmc.mgxs import EnergyGroups, DelayedGroups
from openmc import Mesh
import openmc
import openmc.checkvalue as cv

# Supported cross section types
MDGXS_TYPES = ['delayed-nu-fission',
               'chi-delayed',
               'beta']

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
    domain : openmc.Material or openmc.Cell or openmc.Universe or openmc.Mesh
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
    delayed_groups : openmc.mgxs.DelayedGroups
        Delayed groups to filter out the xs

    Attributes
    ----------
    name : str, optional
        Name of the multi-group cross section
    rxn_type : str
        Reaction type (e.g., 'chi-delayed', 'beta', etc.)
    by_nuclide : bool
        If true, computes cross sections for each nuclide in domain
    domain : Material or Cell or Universe or Mesh
        Domain for spatial homogenization
    domain_type : {'material', 'cell', 'distribcell', 'universe', 'mesh'}
        Domain type for spatial homogenization
    energy_groups : openmc.mgxs.EnergyGroups
        Energy group structure for energy condensation
    delayed_groups : openmc.mgxs.DelayedGroups
        Delayed groups to filter out the xs
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
        Whether or not the MDGXS' tallies use SciPy's LIL sparse matrix format
        for compressed data storage
    loaded_sp : bool
        Whether or not a statepoint file has been loaded with tally data
    derived : bool
        Whether or not the MDGXS is merged from one or more other MDGXS
    hdf5_key : str
        The key used to index multi-group cross sections in an HDF5 data store

    """

    # This is an abstract class which cannot be instantiated
    __metaclass__ = abc.ABCMeta

    def __init__(self, domain=None, domain_type=None, energy_groups=None,
                 by_nuclide=False, name='', delayed_groups=None):
        super(MDGXS, self).__init__(domain, domain_type, energy_groups,
                                    by_nuclide, name)
        self._delayed_groups = None

        if delayed_groups is not None:
            self.delayed_groups = delayed_groups

    def __deepcopy__(self, memo):
        super(MDGXS, self).__deepcopy__(memo)
        existing = memo.get(id(self))

        # If this is the first time we have tried to copy this object, copy it
        if existing is None:
            clone._delayed_groups = copy.deepcopy(self.delayed_groups, memo)

            return clone

        # If this object has been copied before, return the first copy made
        else:
            return existing

    @property
    def delayed_groups(self):
        return self._delayed_groups

    @property
    def num_delayed_groups(self):
        return self.delayed_groups.num_groups

    @delayed_groups.setter
    def delayed_groups(self, delayed_groups):
        cv.check_type('delayed groups', delayed_groups,
                      openmc.mgxs.DelayedGroups)
        self._delayed_groups = delayed_groups

    @property
    def filters(self):

        # Create the non-domain specific Filters for the Tallies
        group_edges = self.energy_groups.group_edges
        energy_filter = openmc.Filter('energy', group_edges)

        if self.delayed_groups != None:
            delayed_groups = self.delayed_groups.groups
            delayed_filter = openmc.Filter('delayedgroup', delayed_groups)
            return [[energy_filter], [delayed_filter, energy_filter]]
        else:
            return [[energy_filter], [energy_filter]]

    @staticmethod
    def get_mgxs(mdgxs_type, domain=None, domain_type=None,
                 energy_groups=None, by_nuclide=False, name='',
                 delayed_groups=None):
        """Return a MDGXS subclass object for some energy group structure within
        some spatial domain for some reaction type.

        This is a factory method which can be used to quickly create MDGXS
        subclass objects for various reaction types.

        Parameters
        ----------
        mdgxs_type : {'delayed-nu-fission', 'chi-delayed', 'beta'}
            The type of multi-delayed-group cross section object to return
        domain : openmc.Material or openmc.Cell or openmc.Universe or
            openmc.Mesh
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
        delayed_groups : openmc.mgxs.DelayedGroups
            Delayed groups to filter out the xs

        Returns
        -------
        openmc.mgxs.MDGXS
            A subclass of the abstract MDGXS class for the multi-delayed-group
            cross section type requested by the user

        """

        cv.check_value('mdgxs_type', mdgxs_type, MDGXS_TYPES)

        if mdgxs_type == 'delayed-nu-fission':
            mdgxs = DelayedNuFissionXS(domain, domain_type, energy_groups)
        elif mdgxs_type == 'chi-delayed':
            mdgxs = ChiDelayed(domain, domain_type, energy_groups)
        elif mdgxs_type == 'beta':
            mdgxs = Beta(domain, domain_type, energy_groups)

        mdgxs.by_nuclide = by_nuclide
        mdgxs.name = name
        mdgxs.delayed_groups = delayed_groups
        return mdgxs

    def get_xs(self, groups='all', subdomains='all', nuclides='all',
               xs_type='macro', order_groups='increasing',
               value='mean', delayed_groups='all', **kwargs):
        """Returns an array of multi-delayed-group cross sections.

        This method constructs a 2D NumPy array for the requested
        multi-delayed-group cross section data data for one or more energy
        groups, delayed groups, and subdomains.

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
        delayed_groups : Iterable of Integral or 'all'
            Delayed groups of interest. Defaults to 'all'.

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
        if not isinstance(subdomains, basestring):
            cv.check_iterable_type('subdomains', subdomains, Integral,
                                   max_depth=3)
            for subdomain in subdomains:
                filters.append(self.domain_type)
                filter_bins.append((subdomain,))

        # Construct list of energy group bounds tuples for all requested groups
        if not isinstance(groups, basestring):
            cv.check_iterable_type('groups', groups, Integral)
            for group in groups:
                filters.append('energy')
                filter_bins.append(
                    (self.energy_groups.get_group_bounds(group),))

        # Construct list of delayed group tuples for all requested groups
        if not isinstance(delayed_groups, basestring):
            cv.check_iterable_type('delayed_groups', delayed_groups, Integral)
            for delayed_group in delayed_groups:
                filters.append('delayedgroup')
                filter_bins.append((delayed_group,))

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

        cv.check_iterable_type('nuclides', nuclides, basestring)
        cv.check_iterable_type('energy_groups', groups, Integral)
        cv.check_iterable_type('delayed_groups', delayed_groups, Integral)

        # Build lists of filters and filter bins to slice
        filters = []
        filter_bins = []

        if len(groups) != 0:
            energy_bins = []
            for group in groups:
                group_bounds = self.energy_groups.get_group_bounds(group)
                energy_bins.append(group_bounds)
            filter_bins.append(tuple(energy_bins))
            filters.append('energy')

        if len(delayed_groups) != 0:
            filter_bins.append(tuple(delayed_groups))
            filters.append('delayedgroup')

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
            slice_xs.delayed_groups.groups = delayed_groups

        # Assign sliced nuclides to sliced MGXS
        if nuclides:
            slice_xs.nuclides = nuclides

        slice_xs.sparse = self.sparse
        return slice_xs

    def can_merge(self, other):
        """Determine if another MDGXS can be merged with this one

        If results have been loaded from a statepoint, then MGXS are only
        mergeable along one and only one of enegy groups or nuclides.

        Parameters
        ----------
        other : openmc.mgxs.MGXS
            MGXS to check for merging

        """

        can_merge = super(MDGXS, self).can_merge(other)

        # Compare delayed groups
        if not self.delayed_groups.can_merge(other.delayed_groups):
            can_merge = False

        # If all conditionals pass then MDGXS are mergeable
        return can_merge

    def merge(self, other):
        """Merge another MDGXS with this one

        MDGXS are only mergeable if their energy groups and nuclides are either
        identical or mutually exclusive. If results have been loaded from a
        statepoint, then MDGXS are only mergeable along one and only one of
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

        merged_mdgxs = super(MDGXS, self).merge(other)

        # Merge delayed groups
        if self.delayed_groups != other.delayed_groups:
            merged_delayed_groups = self.delayed_groups.merge(
                other.delayed_groups)
            merged_mdgxs.delayed_groups = merged_delayed_groups

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

        if self.delayed_groups == None:
            super(MDGXS, self).print_xs(subdomains, nuclides, xs_type)
            return

        # Construct a collection of the subdomains to report
        if not isinstance(subdomains, basestring):
            cv.check_iterable_type('subdomains', subdomains, Integral)
        elif self.domain_type == 'distribcell':
            subdomains = np.arange(self.num_subdomains, dtype=np.int)
        elif self.domain_type == 'mesh':
            xyz = [range(1, x+1) for x in self.domain.dimension]
            subdomains = list(itertools.product(*xyz))
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

                for delayed_group in self.delayed_groups.groups:

                    template = '{0: <12}Delayed Group {1}:\t'
                    string += template.format('', delayed_group)
                    string += '\n'

                    template = '{0: <12}Group {1} [{2: <10} - {3: <10}MeV]:\t'

                    # Loop over energy groups ranges
                    for group in range(1, self.num_groups+1):
                        bounds = self.energy_groups.get_group_bounds(group)
                        string += template.format('', group, bounds[0], bounds[1])
                        average = self.get_xs([group], [subdomain], [nuclide],
                                              xs_type=xs_type, value='mean',
                                              delayed_groups=[delayed_group])
                        rel_err = self.get_xs([group], [subdomain], [nuclide],
                                              xs_type=xs_type, value='rel_err',
                                              delayed_groups=[delayed_group])
                        average = average.flatten()[0]
                        rel_err = rel_err.flatten()[0] * 100.
                        string += '{:.2e} +/- {:1.2e}%'.format(average, rel_err)
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
        delayed_groups : Iterable of Integral or 'all'
            Delayed groups of interest. Defaults to 'all'.

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
                             xs_type='macro', distribcell_paths=True,
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
        distribcell_paths : bool, optional
            Construct columns for distribcell tally filters (default is True).
            The geometric information in the Summary object is embedded into
            a Multi-index column with a geometric "path" to each distribcell
            instance.
        delayed_groups : Iterable of Integral or 'all'
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

        df = super(MDGXS, self).get_pandas_dataframe(groups, nuclides, xs_type,
                                                     distribcell_paths)

        if not isinstance(delayed_groups, basestring):
            cv.check_iterable_type('delayed groups', delayed_groups, Integral)

        # Select out those delayed groups the user requested
        if not isinstance(delayed_groups, basestring):
            if 'delayedgroup' in df:
                df = df[df['delayedgroup'].isin(delayed_groups)]

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

    For post-processing, the :meth:`MDGXS.load_from_statepoint` will pull in the
    necessary data to compute multi-group cross sections from a
    :class:`openmc.StatePoint` instance. The derived multi-group cross
    section can then be obtained from the :attr:`ChiDelayed.xs_tally` property.

    For a spatial domain :math:`V`, energy group :math:`[E_g,E_{g-1}]`, and
    delayed group :math:`d`, the delayed fission spectrum is calculated as:

    .. math::

       \langle \nu^d \sigma_{f,g' \rightarrow g} \phi \rangle &= \int_{r \in V}
       dr \int_{4\pi} d\Omega' \int_0^\infty dE' \int_{E_g}^{E_{g-1}} dE \;
       \chi(E) \nu^d \sigma_f (r, E') \psi(r, E', \Omega')\\
       \langle \nu^d \sigma_f \phi \rangle &= \int_{r \in V} dr \int_{4\pi}
       d\Omega' \int_0^\infty dE' \int_0^\infty dE \; \chi(E) \nu^d \sigma_f (r,
       E') \psi(r, E', \Omega') \\
       \chi_g^d &= \frac{\langle \nu^d \sigma_{f,g' \rightarrow g} \phi \rangle}
       {\langle \nu^d \sigma_f \phi \rangle}

    Parameters
    ----------
    domain : openmc.Material or openmc.Cell or openmc.Universe or openmc.Mesh
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
    delayed_groups : openmc.mgxs.DelayedGroups
        Delayed groups to filter out the xs

    Attributes
    ----------
    name : str, optional
        Name of the multi-group cross section
    rxn_type : str
        Reaction type (e.g., 'total', 'nu-fission', etc.)
    by_nuclide : bool
        If true, computes cross sections for each nuclide in domain
    domain : Material or Cell or Universe or Mesh
        Domain for spatial homogenization
    domain_type : {'material', 'cell', 'distribcell', 'universe', 'mesh'}
        Domain type for spatial homogenization
    energy_groups : openmc.mgxs.EnergyGroups
        Energy group structure for energy condensation
    delayed_groups : openmc.mgxs.DelayedGroups
        Delayed groups to filter out the xs
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
        domain types. When the  This is equal to the number of cell instances
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
                 by_nuclide=False, name='', delayed_groups=None):
        super(ChiDelayed, self).__init__(domain, domain_type, energy_groups,
                                         by_nuclide, name, delayed_groups)
        self._rxn_type = 'chi-delayed'

    @property
    def scores(self):
        return ['delayed-nu-fission', 'delayed-nu-fission']

    @property
    def filters(self):
        # Create the non-domain specific Filters for the Tallies
        group_edges = self.energy_groups.group_edges
        energyout = openmc.Filter('energyout', group_edges)
        energyin = openmc.Filter('energy', [group_edges[0], group_edges[-1]])
        if self.delayed_groups != None:
            delayed_groups = self.delayed_groups.groups
            delayed_filter = openmc.Filter('delayedgroup', delayed_groups)
            return [[delayed_filter, energyin], [delayed_filter, energyout]]
        else:
            return [[energyin], [energyout]]

    @property
    def tally_keys(self):
        return ['delayed-nu-fission-in', 'delayed-nu-fission-out']

    @property
    def estimator(self):
        return 'analog'

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
            energy_filter = delayed_nu_fission_in.find_filter('energy')
            delayed_nu_fission_in.remove_filter(energy_filter)

            # Compute chi
            self._xs_tally = self.rxn_rate_tally / delayed_nu_fission_in
            super(ChiDelayed, self)._compute_xs()

            # Add the coarse energy filter back to the nu-fission tally
            delayed_nu_fission_in.filters.append(energy_filter)

        return self._xs_tally

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
        energy_filter = delayed_nu_fission_in.find_filter('energy')
        delayed_nu_fission_in.remove_filter(energy_filter)

        # Call super class method and null out derived tallies
        slice_xs = super(ChiDelayed, self).get_slice(nuclides, groups,
                                                     delayed_groups)
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
            filters.append('energyout')

        if len(delayed_groups) != 0:
            filter_bins.append(tuple(delayed_groups))
            filters.append('delayedgroup')

        if filters != []:

            # Slice nu-fission-out tally along energyout filter
            delayed_nu_fission_out = slice_xs.tallies['delayed-nu-fission-out']
            tally_slice = delayed_nu_fission_out.get_slice\
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
        other : openmc.mdgxs.MDGXS
            MDGXS to merge with this one

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
            merged_delayed_groups = self.delayed_groups.merge\
                                    (other.delayed_groups)
            merged_mdgxs.delayed_groups = merged_delayed_groups

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
               value='mean', delayed_groups='all', **kwargs):
        """Returns an array of the delayed fission spectrum.

        This method constructs a 2D NumPy array for the requested multi-group
        and multi-delayed group cross section data data for one or more energy
        groups and subdomains.

        Parameters
        ----------
        groups : Iterable of Integral or 'all'
            Energy groups of interest. Defaults to 'all'.
        delayed_groups : Iterable of Integral or 'all'
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
            mirror the parent MDGXS.get_xs(...) class method
        order_groups: {'increasing', 'decreasing'}
            Return the cross section indexed according to increasing or
            decreasing energy groups (decreasing or increasing energies).
            Defaults to 'increasing'.
        value : {'mean', 'std_dev', 'rel_err'}
            A string for the type of value to return. Defaults to 'mean'.

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
        if not isinstance(subdomains, basestring):
            cv.check_iterable_type('subdomains', subdomains, Integral,
                                   max_depth=3)
            for subdomain in subdomains:
                filters.append(self.domain_type)
                filter_bins.append((subdomain,))

        # Construct list of energy group bounds tuples for all requested groups
        if not isinstance(groups, basestring):
            cv.check_iterable_type('groups', groups, Integral)
            for group in groups:
                filters.append('energyout')
                filter_bins.append(
                    (self.energy_groups.get_group_bounds(group),))

        # Construct list of delayed group tuples for all requested groups
        if not isinstance(delayed_groups, basestring):
            cv.check_iterable_type('delayed_groups', delayed_groups, Integral)
            for delayed_group in delayed_groups:
                filters.append('delayedgroup')
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
                nuclides = self.get_all_nuclides()
                delayed_nu_fission_in = delayed_nu_fission_in.summation\
                                        (nuclides=nuclides)
                delayed_nu_fission_out = delayed_nu_fission_out.summation\
                                         (nuclides=nuclides)

                # Remove coarse energy filter to keep it out of tally arithmetic
                energy_filter = delayed_nu_fission_in.find_filter('energy')
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
                nuclides = self.get_all_nuclides()
                xs = self.xs_tally.get_values(filters=filters,
                                              filter_bins=filter_bins,
                                              nuclides=nuclides, value=value)

            # Get chi delayed for user-specified nuclides in the domain
            else:
                cv.check_iterable_type('nuclides', nuclides, basestring)
                xs = self.xs_tally.get_values(filters=filters,
                                              filter_bins=filter_bins,
                                              nuclides=nuclides, value=value)

        # If chi delayed was computed as an average of nuclides in the domain
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

    For post-processing, the :meth:`MDGXS.load_from_statepoint` will pull in the
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
    domain : openmc.Material or openmc.Cell or openmc.Universe or openmc.Mesh
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
    delayed_groups : openmc.mgxs.DelayedGroups
        Delayed groups to filter out the xs

    Attributes
    ----------
    name : str, optional
        Name of the multi-group cross section
    rxn_type : str
        Reaction type (e.g., 'total', 'nu-fission', etc.)
    by_nuclide : bool
        If true, computes cross sections for each nuclide in domain
    domain : Material or Cell or Universe or Mesh
        Domain for spatial homogenization
    domain_type : {'material', 'cell', 'distribcell', 'universe', 'mesh'}
        Domain type for spatial homogenization
    energy_groups : openmc.mgxs.EnergyGroups
        Energy group structure for energy condensation
    delayed_groups : openmc.mgxs.DelayedGroups
        Delayed groups to filter out the xs
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
        domain types. When the  This is equal to the number of cell instances
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
                 by_nuclide=False, name='', delayed_groups=None):
        super(DelayedNuFissionXS, self).__init__(domain, domain_type,
                                                 energy_groups, by_nuclide,
                                                 name, delayed_groups)
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

    For post-processing, the :meth:`MDGXS.load_from_statepoint` will pull in the
    necessary data to compute multi-group cross sections from a
    :class:`openmc.StatePoint` instance. The derived multi-group cross section
    can then be obtained from the :attr:`Beta.xs_tally` property.

    For a spatial domain :math:`V`, energy group :math:`[E_g,E_{g-1}]`, and
    delayed group :math:`d`, the delayed neutron fraction is calculated as:

    .. math::

       \langle \nu^d \sigma_f \phi \rangle &= \int_{r \in V} dr \int_{4\pi}
       d\Omega' \int_0^\infty dE' \int_0^\infty dE \; \chi(E) \nu^d
       \sigma_f (r, E') \psi(r, E', \Omega') \\
       \langle \nu \sigma_f \phi \rangle &= \int_{r \in V} dr \int_{4\pi}
       d\Omega' \int_0^\infty dE' \int_0^\infty dE \; \chi(E) \nu
       \sigma_f (r, E') \psi(r, E', \Omega') \\
       \beta_{d,g} &= \frac{\langle \nu^d \sigma_f \phi \rangle}
       {\langle \nu \sigma_f \phi \rangle}

    Parameters
    ----------
    domain : openmc.Material or openmc.Cell or openmc.Universe or openmc.Mesh
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
    delayed_groups : openmc.mgxs.DelayedGroups
        Delayed groups to filter out the xs

    Attributes
    ----------
    name : str, optional
        Name of the multi-group cross section
    rxn_type : str
        Reaction type (e.g., 'total', 'nu-fission', etc.)
    by_nuclide : bool
        If true, computes cross sections for each nuclide in domain
    domain : Material or Cell or Universe or Mesh
        Domain for spatial homogenization
    domain_type : {'material', 'cell', 'distribcell', 'universe', 'mesh'}
        Domain type for spatial homogenization
    energy_groups : openmc.mgxs.EnergyGroups
        Energy group structure for energy condensation
    delayed_groups : openmc.mgxs.DelayedGroups
        Delayed groups to filter out the xs
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
        domain types. When the  This is equal to the number of cell instances
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
                 by_nuclide=False, name='', delayed_groups=None):
        super(Beta, self).__init__(domain, domain_type, energy_groups,
                                   by_nuclide, name, delayed_groups)
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
            super(Beta, self)._compute_xs()

        return self._xs_tally
