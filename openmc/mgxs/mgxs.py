from collections import Iterable
from numbers import Integral
import os
import sys
import copy
import abc
import pickle

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
                'material',
                'mesh']

# Supported domain objects
DOMAINS = [openmc.Cell,
           openmc.Universe,
           openmc.Material,
           openmc.Mesh]


class MultiGroupXS(object):
    """A multi-group cross-section for some energy group structure within
    some spatial domain.

    This class can be used for both OpenMC input generation and tally data
    post-processing to compute spatially-homogenized and energy-integrated
    multi-group cross-sections for deterministic neutronics calculations.

    Parameters
    ----------
    domain : Material or Cell or Universe or Mesh
        The domain for spatial homogenization
    domain_type : {'material', 'cell', 'distribcell', 'universe' or 'mesh'}
        The domain type for spatial homogenization
    energy_groups : EnergyGroups
        The energy group structure for energy condensation
    name : str, optional
        Name of the multi-group cross-section. Used as a label to identify
        tallies in OpenMC tallies.xml file.

    Attributes
    ----------
    name : str, optional
        Name of the multi-group cross-section
    xs_type : str
        Cross-section type (e.g., 'total', 'nu-fission', etc.)
    domain : Material or Cell or Universe or Mesh
        Domain for spatial homogenization
    domain_type : {'material', 'cell', 'distribcell', 'universe' or 'mesh'}
        Domain type for spatial homogenization
    energy_groups : EnergyGroups
        Energy group structure for energy condensation
    num_groups : Integral
        Number of energy groups
    tallies : dict
        OpenMC tallies needed to compute the multi-group cross-section
    xs_tally : Tally
        Derived tally for the multi-group cross-section. This attribute
        is None unless the multi-group cross-section has been computed.

    """

    # This is an abstract class which cannot be instantiated
    __metaclass__ = abc.ABCMeta

    def __init__(self, domain=None, domain_type=None,
                 energy_groups=None, name=''):

        self._name = ''
        self._xs_type = None
        self._domain = None
        self._domain_type = None
        self._energy_groups = None
        self._num_groups = None
        self._tallies = dict()
        self._xs_tally = None

        self.name = name
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
            clone._xs_type = self.xs_type
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
    def xs_type(self):
        return self._xs_type

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

    @abc.abstractmethod
    def create_tallies(self, scores, all_filters, keys, estimator):
        """Instantiates tallies needed to compute the multi-group cross-section.

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
        cv.check_iterable_type('filters', all_filters, openmc.Filter, 1, 2)
        cv.check_type('keys', keys, Iterable, basestring)
        cv.check_length('scores', scores, len(keys))
        cv.check_value('estimator', estimator, ['analog', 'tracklength'])

        # Create a domain Filter object
        domain_filter = openmc.Filter(self.domain_type, self.domain.id)
        domain_filter.num_bins = 1

        for score, key, filters in zip(scores, keys, all_filters):
            self.tallies[key] = openmc.Tally(name=self.name)
            self.tallies[key].add_score(score)
            self.tallies[key].estimator = estimator
            self.tallies[key].add_filter(domain_filter)

            # Add all non-domain specific Filters (i.e., 'energy') to the Tally
            for filter in filters:
                self.tallies[key].add_filter(filter)

    @abc.abstractmethod
    def compute_xs(self):
        """Computes multi-group cross-sections using OpenMC tally arithmetic."""

        return

    def load_from_statepoint(self, statepoint):
        """Extracts tallies in an OpenMC StatePoint with the data needed to
        compute multi-group cross-sections.

        This method is needed to compute cross-section data from tallies
        in an OpenMC StatePoint object.

        Parameters
        ----------
        statepoint : openmc.StatePoint
            An OpenMC StatePoint object with tally data

        """

        cv.check_type('statepoint', statepoint, openmc.statepoint.StatePoint)

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

    def get_xs(self, groups='all', subdomains='all', value='mean'):
        """Returns an array of multi-group cross-sections.

        This method constructs a 2D NumPy array for the requested multi-group
        cross-section data data for one or more energy groups and subdomains.

        Parameters
        ----------
        groups : Iterable of Integral or 'all'
            Energy groups of interest

        subdomains : Iterable of Integral or 'all'
            Subdomain IDs of interest

        value : str
            A string for the type of value to return - 'mean' (default),
            'std_dev' or 'rel_err' are accepted

        Returns
        -------
        xs : ndarray
            A NumPy array of the multi-group cross-section indexed in the order
            each group and subdomain is listed in the parameters.

        Raises
        ------
        ValueError
            When this method is called before the multi-group cross-section is
            computed from tally data.

        """

        if self.xs_tally is None:
            msg = 'Unable to get cross-section since it has not been computed'
            raise ValueError(msg)

        cv.check_value('value', value, ['mean', 'std_dev', 'rel_err'])

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

        # Query the multi-group cross-section tally for the data
        xs = self.xs_tally.get_values(filters=filters,
                                      filter_bins=filter_bins, value=value)
        return xs

    def get_condensed_xs(self, coarse_groups):
        """Construct an energy-condensed version of this cross-section.

        Parameters
        ----------
        coarse_groups : openmc.mgxs.EnergyGroups
            The coarse energy group structure of interest

        Returns
        -------
        MultiGroupXS
            A new MultiGroupXS condensed to the group structure of interest

        """

        raise NotImplementedError('Energy condensation is not yet implemented')

    def get_subdomain_avg_xs(self, subdomains='all'):
        """Construct a subdomain-averaged version of this cross-section.

        This is primarily useful for averaging across distribcell instances.

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
            When this method is called before the multi-group cross-section is
            computed from tally data.

        """

        if self.xs_tally is None:
            msg = 'Unable to get cross-section since it has not been computed'
            raise ValueError(msg)

        # Construct a collection of the subdomain filter bins to average across
        if subdomains == 'all':
            if self.domain_type == 'distribcell':
                subdomains = np.arange(self.num_subdomains)
            else:
                subdomains = [self.domain.id]
        else:
            cv.check_iterable_type('subdomains', subdomains, Integral)

        # Clone this MultiGroupXS to initialize the condensed version
        avg_xs = copy.deepcopy(self)

        # Reset subdomain indices and offsets for distribcell domains
        if self.domain_type == 'distribcell':
            avg_xs._offset = 0

        # Overwrite tallies with new subdomain-averaged versions
        avg_xs._tallies = {}
        for tally_type, tally in self.tallies.items():
            tally_sum = tally.summation(filter_type=self.domain_type,
                                        filter_bins=subdomains)
            tally_sum /= len(subdomains)
            avg_xs.tallies[tally_type] = tally_sum

        # Compute the condensed single group cross-section
        avg_xs.compute_xs()

        return avg_xs

    def print_xs(self, subdomains='all'):
        """Prints a string representation for the multi-group cross-section.

        Parameters
        ----------
        subdomains : Iterable of Integral or 'all'
            The subdomain IDs of the cross-sections to include in the report

        """

        if subdomains != 'all':
            cv.check_iterable_type('subdomains', subdomains, Integral)

        # Build header for string with type and domain info
        string = 'Multi-Group XS\n'
        string += '{0: <16}=\t{1}\n'.format('\tType', self.xs_type)
        string += '{0: <16}=\t{1}\n'.format('\tDomain Type', self.domain_type)
        string += '{0: <16}=\t{1}\n'.format('\tDomain ID', self.domain.id)

        # Append cross-section data if it has been computed
        if self.xs_tally is not None:
            if subdomains == 'all':
                if self.domain_type == 'distribcell':
                    subdomains = np.arange(self.num_subdomains, dtype=np.int)
                else:
                    subdomains = [self.domain.id]

            # Loop over all subdomains
            for subdomain in subdomains:

                if self.domain_type == 'distribcell':
                    string += '{0: <16}=\t{1}\n'.format('\tSubdomain', subdomain)

                string += '{0: <16}\n'.format('\tCross-Sections [cm^-1]:')
                template = '{0: <12}Group {1} [{2: <10} - {3: <10}MeV]:\t'

                # Loop over energy groups ranges
                for group in range(1, self.num_groups+1):
                    bounds = self.energy_groups.get_group_bounds(group)
                    string += template.format('', group, bounds[0], bounds[1])
                    average = self.get_xs([group], [subdomain], 'mean')
                    rel_err = self.get_xs([group], [subdomain], 'rel_err')*100.
                    average = np.nan_to_num(average.flatten())[0]
                    rel_err = np.nan_to_num(rel_err.flatten())[0]
                    string += '{:.2e} +/- {:1.2e}%'.format(average, rel_err)
                    string += '\n'
                string += '\n'

        print(string)

    def pickle(self, filename='mgxs', directory='mgxs'):
        """Store the MultiGroupXS as a pickled binary file.

        Parameters
        ----------
        filename : str
            Filename for the pickled binary file (default is 'mgxs')

        directory : str
            Directory for the pickled binary file (default is 'mgxs')

        """

        cv.check_type('filename', filename, basestring)
        cv.check_type('directory', directory, basestring)

        # Make directory if it does not exist
        if not os.path.exists(directory):
            os.makedirs(directory)

        # Create an empty dictionary to store the data
        xs_results = dict()

        # Store all of this MultiGroupXS' class attributes in the dictionary
        xs_results['name'] = self.name
        xs_results['xs_type'] = self.xs_type
        xs_results['domain_type'] = self.domain_type
        xs_results['domain'] = self.domain
        xs_results['energy_groups'] = self.energy_groups
        xs_results['tallies'] = self.tallies
        xs_results['xs_tally'] = self.xs_tally
        xs_results['offset'] = self.offset
        xs_results['subdomain_indices'] = self.subdomain_indices

        # Pickle the MultiGroupXS results to a binary file
        filename = directory + '/' + filename + '.pkl'
        filename = filename.replace(' ', '-')
        pickle.dump(xs_results, open(filename, 'wb'))

    def unpickle(self, filename='mgxs', directory='mgxs'):
        """Restore the MultiGroupXS from a pickled binary file.

        Parameters
        ----------
        filename : str
            Filename for the pickled binary file (default is 'mgxs')

        directory : str
            Directory for the pickled binary file (default is 'mgxs')

        Raises
        ------
        ValueError
            When the requested filename does not exist.

        """

        cv.check_type('filename', filename, basestring)
        cv.check_type('directory', directory, basestring)

        filename = directory + '/' + filename + '.pkl'
        filename = filename.replace(' ', '-')

        # Check that the file exists
        if not os.path.exists(filename):
            msg = 'Unable to import from filename="{0}"'.format(filename)
            raise ValueError(msg)

        # Load the pickle file into a dictionary
        xs_results = pickle.load(open(filename, 'rb'))

        # Store the MultiGroupXS class attributes
        self.name = xs_results['name']
        self._xs_type = xs_results['xs_type']
        self.domain_type = xs_results['domain_type']
        self.domain = xs_results['domain']
        self.energy_groups = xs_results['energy_groups']
        self._tallies = xs_results['tallies']
        self._xs_tally = xs_results['xs_tally']
        self._offset = xs_results['offset']

    def build_hdf5_store(self, filename='mgxs', directory='mgxs', append=True):
        """Export the multi-group cross-section data into an HDF5 binary file.

        This routine constructs an HDF5 file which stores the multi-group
        cross-section data. The data is be stored in a hierarchy of HDF5 groups
        from the domain type, domain id, subdomain id (for distribcell domains),
        and cross-section type. Two datasets for the mean and standard deviation
        are stored for each subddomain entry in the HDF5 file.

        NOTE: This requires the h5py Python package.

        Parameters
        ----------
        filename : str
            Filename for the HDF5 file (default is 'mgxs')

        directory : str
            Directory for the HDF5 file (default is 'mgxs')

        append : boolean
            If true, appends to an existing HDF5 file with the same filename
            directory (if one exists)

        Raises
        ------
        ValueError
            When this method is called before the multi-group cross-section is
            computed from tally data.
        ImportError
            When h5py is not installed.

        """

        if self.xs_tally is None:
            msg = 'Unable to get build HDF5 store since the ' \
                  'cross-section has not been computed'
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

        # Create a separate HDF5 dataset for each subdomain
        for i, subdomain in enumerate(subdomains):

            # Create an HDF5 group for the subdomain
            if self.domain_type == 'distribcell':
                group_name = str(subdomain).zfill(num_digits)
                subdomain_group = domain_group.require_group(group_name)
            else:
                subdomain_group = domain_group

            # Create a separate HDF5 group for the xs type
            xs_group = subdomain_group.require_group(self.xs_type)

            # Extract the cross-section for this
            average = self.get_xs(subdomains=[subdomain], value='mean')
            std_dev = self.get_xs(subdomains=[subdomain], value='std_dev')
            average = average.squeeze()
            std_dev = std_dev.squeeze()

            # Add MultiGroupXS results data to the HDF5 group
            xs_group.require_dataset('average', dtype=np.float64,
                                     shape=average.shape, data=average)
            xs_group.require_dataset('std. dev.', dtype=np.float64,
                                     shape=std_dev.shape, data=std_dev)

        # Close the MultiGroup results HDF5 file
        xs_results.close()

    def export_xs_data(self, filename='mgxs', directory='mgxs', format='csv'):
        """Export the multi-group cross-section data to a file.

        This routine leverages the functionality in the Pandas library to
        export the multi-group cross-section data in a variety of output
        file formats for storage and/or post-processing.

        Parameters
        ----------
        filename : str
            Filename for the exported file (default is 'mgxs')

        directory : str
            Directory for the exported file (default is 'mgxs')

        format : {'csv', 'excel', 'pickle', 'latex'}
            The format for the exported data file

        groups : {'indices' or 'bounds'}
            When set to 'indices' (default), integer group indices are inserted
            in the energy in/out column(s) of the DataFrame. When it is 'bounds'
            the lower and upper energy bounds are used.

        summary : None or Summary
            An optional Summary object to be used to construct columns for
            distribcell tally filters (default is None). The geometric
            information in the Summary object is embedded into a Multi-index
            column with a geometric "path" to each distribcell intance.
            NOTE: This option requires the OpenCG Python package.

        """

        cv.check_type('filename', filename, basestring)
        cv.check_type('directory', directory, basestring)
        cv.check_value('format', format, ['csv', 'excel', 'pickle', 'latex'])

        # Make directory if it does not exist
        if not os.path.exists(directory):
            os.makedirs(directory)

        filename = directory + '/' + filename
        filename = filename.replace(' ', '-')

        # Get a Pandas DataFrame for the data
        df = self.get_pandas_dataframe()

        # Capitalize column label strings
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
                msg = 'Unable to export distribcell multi-group cross-section' \
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


    def get_pandas_dataframe(self, groups='indices', summary=None):
        """Build a Pandas DataFrame for the MultiGroupXS data.

        This routine leverages the Tally.get_pandas_dataframe(...) routine, but
        renames the columns with terminology appropriate for cross-section data.

        Parameters
        ----------
        groups : {'indices' or 'bounds'}
            When set to 'indices' (default), integer group indices are inserted
            in the energy in/out column(s) of the DataFrame. When it is 'bounds'
            the lower and upper energy bounds are used.

        summary : None or Summary
            An optional Summary object to be used to construct columns for
            distribcell tally filters (default is None). The geometric
            information in the Summary object is embedded into a Multi-index
            column with a geometric "path" to each distribcell intance.
            NOTE: This option requires the OpenCG Python package.

        Returns
        -------
        pandas.DataFrame
            A Pandas DataFrame for the cross-section data.

        Raises
        ------
        ValueError
            When this method is called before the multi-group cross-section is
            computed from tally data.

        """

        if self.xs_tally is None:
            msg = 'Unable to get Pandas DataFrame since the ' \
                  'cross-section has not been computed'
            raise ValueError(msg)

        # Get a Pandas DataFrame from the derived xs tally
        df = self.xs_tally.get_pandas_dataframe(summary=summary)

        # Remove the score column since it is homogeneous and redundant
        if summary and self.domain_type == 'distribcell':
            df = df.drop('score', level=0, axis=1)
        else:
            df = df.drop('score', axis=1)

        # Use group indices in place of energy bounds ("1" for fastest group)
        if groups == 'indices':

            # Rename the column label for energy in the dataframe
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

            # Sort the dataframe by domain type id (e.g., distribcell id) and
            # energy groups such that data is from fast to thermal
            df.sort([self.domain_type] + columns, inplace=True)

        return df


class TotalXS(MultiGroupXS):

    def __init__(self, domain=None, domain_type=None, groups=None, name=''):
        super(TotalXS, self).__init__(domain, domain_type, groups, name)
        self._xs_type = 'total'

    def create_tallies(self):
        """Construct the OpenMC tallies needed to compute this cross-section."""

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
        """Computes the multi-group total cross-sections using OpenMC
        tally arithmetic."""

        self._xs_tally = self.tallies['total'] / self.tallies['flux']
        self._xs_tally._mean = np.nan_to_num(self.xs_tally.mean)
        self._xs_tally._std_dev = np.nan_to_num(self.xs_tally.std_dev)


class TransportXS(MultiGroupXS):

    def __init__(self, domain=None, domain_type=None, groups=None, name=''):
        super(TransportXS, self).__init__(domain, domain_type, groups, name)
        self._xs_type = 'transport'

    def create_tallies(self):
        """Construct the OpenMC tallies needed to compute this cross-section."""

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

    def compute_xs(self):
        """Computes the multi-group transport cross-sections using OpenMC
        tally arithmetic."""

        self._xs_tally = self.tallies['total'] - self.tallies['scatter-P1']
        self._xs_tally /= self.tallies['flux']
        self._xs_tally._mean = np.nan_to_num(self.xs_tally.mean)
        self._xs_tally._std_dev = np.nan_to_num(self.xs_tally.std_dev)


class AbsorptionXS(MultiGroupXS):

    def __init__(self, domain=None, domain_type=None, groups=None, name=''):
        super(AbsorptionXS, self).__init__(domain, domain_type, groups, name)
        self._xs_type = 'absorption'

    def create_tallies(self):
        """Construct the OpenMC tallies needed to compute this cross-section."""

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
        """Computes the multi-group absorption cross-sections using OpenMC
        tally arithmetic."""

        self._xs_tally = self.tallies['absorption'] / self.tallies['flux']
        self._xs_tally._mean = np.nan_to_num(self.xs_tally.mean)
        self._xs_tally._std_dev = np.nan_to_num(self.xs_tally.std_dev)


class CaptureXS(MultiGroupXS):

    def __init__(self, domain=None, domain_type=None, groups=None, name=''):
        super(CaptureXS, self).__init__(domain, domain_type, groups, name)
        self._xs_type = 'capture'

    def create_tallies(self):
        """Construct the OpenMC tallies needed to compute this cross-section."""

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
        """Computes the multi-group capture cross-sections using OpenMC
        tally arithmetic."""

        self._xs_tally = self.tallies['absorption'] - self.tallies['fission']
        self._xs_tally /= self.tallies['flux']
        self._xs_tally._mean = np.nan_to_num(self.xs_tally.mean)
        self._xs_tally._std_dev = np.nan_to_num(self.xs_tally.std_dev)


class FissionXS(MultiGroupXS):

    def __init__(self, domain=None, domain_type=None, groups=None, name=''):
        super(FissionXS, self).__init__(domain, domain_type, groups, name)
        self._xs_type = 'fission'

    def create_tallies(self):
        """Construct the OpenMC tallies needed to compute this cross-section."""

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
        """Computes the multi-group fission cross-sections using OpenMC
        tally arithmetic."""

        self._xs_tally = self.tallies['fission'] / self.tallies['flux']
        self._xs_tally._mean = np.nan_to_num(self.xs_tally.mean)
        self._xs_tally._std_dev = np.nan_to_num(self.xs_tally.std_dev)


class NuFissionXS(MultiGroupXS):

    def __init__(self, domain=None, domain_type=None, groups=None, name=''):
        super(NuFissionXS, self).__init__(domain, domain_type, groups, name)
        self._xs_type = 'nu-fission'

    def create_tallies(self):
        """Construct the OpenMC tallies needed to compute this cross-section."""

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
        """Computes the multi-group nu-fission cross-sections using OpenMC
        tally arithmetic."""

        self._xs_tally = self.tallies['nu-fission'] / self.tallies['flux']
        self._xs_tally._mean = np.nan_to_num(self.xs_tally.mean)
        self._xs_tally._std_dev = np.nan_to_num(self.xs_tally.std_dev)


class ScatterXS(MultiGroupXS):

    def __init__(self, domain=None, domain_type=None, groups=None, name=''):
        super(ScatterXS, self).__init__(domain, domain_type, groups, name)
        self._xs_type = 'scatter'

    def create_tallies(self):
        """Construct the OpenMC tallies needed to compute this cross-section."""

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
        """Computes the scattering multi-group cross-sections using
        OpenMC tally arithmetic."""

        self._xs_tally = self.tallies['scatter'] / self.tallies['flux']
        self._xs_tally._mean = np.nan_to_num(self.xs_tally.mean)
        self._xs_tally._std_dev = np.nan_to_num(self.xs_tally.std_dev)


class NuScatterXS(MultiGroupXS):

    def __init__(self, domain=None, domain_type=None, groups=None, name=''):
        super(NuScatterXS, self).__init__(domain, domain_type, groups, name)
        self._xs_type = 'nu-scatter'

    def create_tallies(self):
        """Construct the OpenMC tallies needed to compute this cross-section."""

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
        """Computes the nu-scattering multi-group cross-section using OpenMC
        tally arithmetic."""

        self._xs_tally = self.tallies['nu-scatter'] / self.tallies['flux']
        self._xs_tally._mean = np.nan_to_num(self.xs_tally.mean)
        self._xs_tally._std_dev = np.nan_to_num(self.xs_tally.std_dev)


class ScatterMatrixXS(MultiGroupXS):

    def __init__(self, domain=None, domain_type=None, groups=None, name=''):
        super(ScatterMatrixXS, self).__init__(domain, domain_type, groups, name)
        self._xs_type = 'scatter matrix'

    def create_tallies(self):
        """Construct the OpenMC tallies needed to compute this cross-section."""

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

    def compute_xs(self, correction='None'):
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
            scatter_p1 = self.tallies['scatter-1']
            scatter_p1 = scatter_p1.get_slice(scores=['scatter-P1'])
            energy_filter = openmc.Filter(type='energy')
            energy_filter.bins = self.energy_groups.group_edges
            scatter_p1 = scatter_p1.diagonalize_filter(energy_filter)
            rxn_tally = self.tallies['scatter'] - scatter_p1
        else:
            rxn_tally = self.tallies['scatter']

        self._xs_tally = rxn_tally / self.tallies['flux']
        self._xs_tally._mean = np.nan_to_num(self.xs_tally.mean)
        self._xs_tally._std_dev = np.nan_to_num(self.xs_tally.std_dev)

    def get_xs(self, in_groups='all', out_groups='all',
               subdomains='all', value='mean'):
        """Returns an array of multi-group cross-sections.

        This method constructs a 2D NumPy array for the requested multi-group
        cross-section data data for one or more energy groups and subdomains.

        Parameters
        ----------
        in_groups : Iterable of Integral or 'all'
            Incoming energy groups of interest

        out_groups : Iterable of Integral or 'all'
            Outgoing energy groups of interest

        subdomains : Iterable of Integral or 'all'
            Subdomain IDs of interest

        value : str
            A string for the type of value to return - 'mean' (default),
            'std_dev' or 'rel_err' are accepted

        Returns
        -------
        xs : ndarray
            A NumPy array of the multi-group cross-section indexed in the order
            each group and subdomain is listed in the parameters.

        Raises
        ------
        ValueError
            When this method is called before the multi-group cross-section is
            computed from tally data.

        """

        if self.xs_tally is None:
            msg = 'Unable to get cross-section since it has not been computed'
            raise ValueError(msg)

        cv.check_value('value', value, ['mean', 'std_dev', 'rel_err'])

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

        # Query the multi-group cross-section tally for the data
        xs = self.xs_tally.get_values(filters=filters,
                                      filter_bins=filter_bins, value=value)
        xs = np.nan_to_num(xs)
        return xs

    def print_xs(self, subdomains='all'):
        """Prints a string representation for the multi-group cross-section.

        Parameters
        ----------
        subdomains : Iterable of Integral or 'all'
            The subdomain IDs of the cross-sections to include in the report

        """

        if subdomains != 'all':
            cv.check_iterable_type('subdomains', subdomains, Integral)

        # Build header for string with type and domain info
        string = 'Multi-Group XS\n'
        string += '{0: <16}=\t{1}\n'.format('\tType', self.xs_type)
        string += '{0: <16}=\t{1}\n'.format('\tDomain Type', self.domain_type)
        string += '{0: <16}=\t{1}\n'.format('\tDomain ID', self.domain.id)

        # Append cross-section data if it has been computed
        if self.xs_tally is not None:
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

                string += '{0: <16}\n'.format('\tCross-Sections [cm^-1]:')
                template = '{0: <12}Group {1} -> Group {2}:\t\t'

                # Loop over incoming/outgoing energy groups ranges
                for in_group in range(1, self.num_groups+1):
                    for out_group in range(1, self.num_groups+1):
                        string += template.format('', in_group, out_group)
                        average = self.get_xs([in_group], [out_group],
                                              [subdomain], 'mean')
                        rel_err = self.get_xs([in_group], [out_group],
                                              [subdomain], 'rel_err') * 100.
                        average = np.nan_to_num(average.flatten())[0]
                        rel_err = np.nan_to_num(rel_err.flatten())[0]
                        string += '{:1.2e} +/- {:1.2e}%'.format(average, rel_err)
                        string += '\n'
                    string += '\n'

        print(string)


class NuScatterMatrixXS(ScatterMatrixXS):

    def __init__(self, domain=None, domain_type=None, groups=None,  name=''):
        super(NuScatterMatrixXS, self).__init__(domain, domain_type, groups, name)
        self._xs_type = 'nu-scatter matrix'

    def create_tallies(self):
        """Construct the OpenMC tallies needed to compute this cross-section."""

        # Create a list of scores for each Tally to be created
        scores = ['flux', 'nu-scatter', 'scatter-P1']
        estimator = 'analog'
        keys = scores

        # Create the non-domain specific Filters for the Tallies
        group_edges = self.energy_groups.group_edges
        energy = openmc.Filter('energy', group_edges)
        energyout = openmc.Filter('energyout', group_edges)
        filters = [[energy], [energy, energyout], [energyout]]

        # Intialize the Tallies
        super(ScatterMatrixXS, self).create_tallies(scores, filters, keys, estimator)

    def compute_xs(self, correction='None'):
        """Computes the multi-group nu-scattering matrix using OpenMC
        tally arithmetic.

        Parameters
        ----------
        correction : {'P0' or None}
            If 'P0', applies the P0 transport correction to the diagonal of the
            scattering matrix

        """

        # If using P0 correction subtract scatter-P1 from the diagonal
        if correction == 'P0':
            scatter_p1 = self.tallies['scatter-1']
            scatter_p1 = scatter_p1.get_slice(scores=['scatter-P1'])
            energy_filter = openmc.Filter(type='energy')
            energy_filter.bins = self.energy_groups.group_edges
            scatter_p1 = scatter_p1.diagonalize_filter(energy_filter)
            rxn_tally = self.tallies['nu-scatter'] - scatter_p1
        else:
            rxn_tally = self.tallies['nu-scatter']

        self._xs_tally = rxn_tally / self.tallies['flux']
        self._xs_tally._mean = np.nan_to_num(self.xs_tally.mean)
        self._xs_tally._std_dev = np.nan_to_num(self.xs_tally.std_dev)

class Chi(MultiGroupXS):

    def __init__(self, domain=None, domain_type=None, groups=None, name=''):
        super(Chi, self).__init__(domain, domain_type, groups, name)
        self._xs_type = 'chi'

    def create_tallies(self):
        """Construct the OpenMC tallies needed to compute this cross-section."""

        # Create a list of scores for each Tally to be created
        scores = ['nu-fission', 'nu-fission']
        estimator = 'analog'
        keys = ['nu-fission-in', 'nu-fission-out']

        # Create the non-domain specific Filters for the Tallies
        group_edges = self.energy_groups.group_edges
        energyout_filter = openmc.Filter('energyout', group_edges)
        filters = [[], [energyout_filter]]

        # Intialize the Tallies
        super(Chi, self).create_tallies(scores, filters, keys, estimator)

    def compute_xs(self):
        """Computes chi fission spectrum using OpenMC tally arithmetic."""

        nu_fission_in = self.tallies['nu-fission-in']
        nu_fission_out = self.tallies['nu-fission-out']
        self._xs_tally = nu_fission_out / nu_fission_in
        self._xs_tally._mean = np.nan_to_num(self.xs_tally.mean)
        self._xs_tally._std_dev = np.nan_to_num(self.xs_tally.std_dev)
