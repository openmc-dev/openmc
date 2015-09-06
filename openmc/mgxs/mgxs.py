from collections import Iterable
from numbers import Integral
import os
import sys
import copy
import abc
import pickle
import subprocess

import numpy as np

import openmc
import openmc.checkvalue as cv
from openmc.mgxs import EnergyGroups


if sys.version_info[0] >= 3:
    basestring = str


# Supported cross-section types
XS_TYPES = ['total',
            'transport',
            'absorption',
            'capture',
            'scatter',
            'nu-scatter',
            'scatter matrix',
            'nu-scatter matrix',
            'fission',
            'nu-fission',
            'chi']

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

# LaTeX Greek symbols for each cross-section type
GREEK = dict()
GREEK['total'] = '$\\Sigma_{t}$'
GREEK['transport'] = '$\\Sigma_{tr}$'
GREEK['absorption'] = '$\\Sigma_{a}$'
GREEK['capture'] = '$\\Sigma_{c}$'
GREEK['scatter'] = '$\\Sigma_{s}$'
GREEK['nu-scatter'] = '$\\nu\\Sigma_{s}$'
GREEK['scatter matrix'] = '$\\Sigma_{s}$'
GREEK['nu-scatter matrix'] = '$\\nu\\Sigma_{s}$'
GREEK['fission'] = '$\\Sigma_{f}$'
GREEK['nu-fission'] = '$\\nu\\Sigma_{f}$'
GREEK['chi'] = '$\\chi$'


class MultiGroupXS(object):
    """A multi-group cross-section for some energy groups structure within
    some spatial domain.

    This class can be used for both OpenMC input generation and tally data
    post-processing to compute spatially-homogenized and energy-integrated
    multi-group cross-sections for deterministic neutronics calculations.

    Parameters
    ----------
    name : str, optional
        Name of the multi-group cross-section. If not specified, the name is
        the empty string.
    domain : Material or Cell or Universe or Mesh
        The domain for spatial homogenization
    domain_type : {'material', 'cell', 'distribcell', 'universe' or 'mesh'}
        The domain type for spatial homogenization
    energy_groups : EnergyGroups
        The energy group structure for energy condensation

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
        Tallies needed to compute the multi-group cross-section
    xs : Tally
        Derived tally for the multi-group cross-section. This attribute
        is None unless the multi-group cross-section has been computed.
    subdomain_offsets : dict
        Integral subdomain IDs (keys) mapped to integral tally data array
        offsets (values). When the domain_type is 'distribcell', each subdomain
        ID corresponds to an instance of the cell domain. For all other domain
        types, there is only one subdomain for the domain and this dictionary
        will trivially map zero to zero.
    offset : Integral
        The filter offset for the domain filter

    """

    # This is an abstract class which cannot be instantiated
    metaclass__ = abc.ABCMeta

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

        # A dictionary used to compute indices into the xs array
        # Keys   - Domain ID (ie, Material ID, Region ID for districell, etc)
        # Values - Offset/stride into xs array
        # NOTE: This is primarily used for distribcell domain types
        self._subdomain_indices = dict()
        self._offset = None

        self.name = name
        if not domain_type is None:
            self.domain_type = domain_type
        if not domain is None:
            self.domain = domain
        if not energy_groups is None:
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
            clone._subdomain_indices = copy.deepcopy(self.subdomain_indices, memo)
            clone._offset = copy.deepcopy(self.offset, memo)

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
    def offset(self):
        return self._offset

    @property
    def subdomain_indices(self):
        return self._subdomain_indices

    @name.setter
    def name(self, name):
        cv.check_type('name', name, basestring)
        self._name = name

    @domain.setter
    def domain(self, domain):
        cv.check_type('domain', domain, DOMAINS)
        self._domain = domain
        if self._domain_type in ['material', 'cell', 'universe', 'mesh']:
            self._subdomain_indices[domain.id] = 0

    @domain_type.setter
    def domain_type(self, domain_type):
        cv.check_type('domain type', domain_type, DOMAIN_TYPES)
        self._domain_type = domain_type

    @energy_groups.setter
    def energy_groups(self, energy_groups):
        cv.check_type('energy groups', energy_groups, openmc.mgxs.EnergyGroups)
        self._energy_groups = energy_groups
        self._num_groups = energy_groups._num_groups

    def _find_domain_offset(self):
        """Finds and stores the offset of the domain tally filter"""
        tally = self.tallies.values()[0]
        filter = tally.find_filter(self.domain_type, [self.domain.id])
        self._offset = filter.offset

    def set_subdomain_index(self, subdomain_id, index):
        """Set the filter bin index for a subdomain of the domain.

        This is primarily useful when the domain type is 'distribcell', in
        which case one may wish to map each subdomain (a cell instance) to its
        filter bin in the derived multi-group cross-section tally data array.

        Parameters
        ----------
        subdomain_id : Integral
            The ID for the subdomain
        index : Integral
            The filter bin index for the subdomain

        """

        cv.check_type('subdomain id', subdomain_id, Integral)
        cv.check_type('subdomain offset', index, Integral)
        cv.check_greater_than('subdomain id', subdomain_id, 0, True)
        cv.check_greater_than('subdomain offset', subdomain_id, 0, True)
        self._subdomain_indices[subdomain_id] = index

    @abc.abstractmethod
    def _create_tallies(self, scores, all_filters, keys, estimator):
        """Instantiates tallies needed to compute the multi-group cross-section

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

        cv.check_value('scores', scores, openmc.SCORE_TYPES)
        cv.check_iterable_type('filters', all_filters, openmc.Filter, 1, 2)
        cv.check_type('keys', keys, Iterable, basestring)
        cv.check_length('scores', scores, len(keys))
        cv.check_value('estimator', estimator, ['analog', 'tracklength'])

        # Create a domain Filter object
        domain_filter = openmc.Filter(self.domain_type, self.domain.id)

        for score, key, filters in zip(scores, keys, all_filters):
            self.tallies[key] = openmc.Tally(name=self.name)
            self.tallies[key].add_score(score)
            self.tallies[key].estimator = estimator
            self.tallies[key].add_filter(domain_filter)

            # Add all non-domain specific Filters (i.e., 'energy') to the Tally
            for filter in filters:
                self.tallies[key].add_filter(filter)

    def get_subdomain_indices(self, subdomains='all'):
        """Get the indices for one or more subdomains.

        This method can be used to extract the indices into the multi-group
        cross-section tally data array for a subdomain (i.e., cell instance).

        See also : get_subdomains

        Parameters
        ----------
        subdomains : Iterable of Integral or 'all'
            Subdomain IDs of interest

        Returns
        ----------
        indices : ndarray
            The subdomain indices indexed in the order of the subdomains

        Raises
        ------
        ValueError
            When one of the subdomains is not a valid subdomain ID.

        """

        if subdomains != 'all':
            cv.check_type('subdomains', subdomains, Iterable, Integral)

        if subdomains == 'all':
            num_subdomains = len(self.subdomain_indices)
            indices = np.arange(num_subdomains)
        else:
            indices = np.zeros(len(subdomains), dtype=np.int64)

            for i, subdomain in enumerate(subdomains):
                if subdomain in self.subdomain_indices:
                    indices[i] = self.subdomain_indices[subdomain]
                else:
                    msg = 'Unable to get index for subdomain "{0}" since it ' \
                          'is not a valid subdomain'.format(subdomain)
                    raise ValueError(msg)

        return indices

    def get_subdomains(self, indices='all'):
        """Get the subdomain IDs for one or more indices.

        This method can be used to extract the subdomains for the multi-group
        cross-section from their indices in the tally data array.

        See also : get_subdomain_indices

        Parameters
        ----------
        indices : Iterable of Integral or 'all'
            Subdomain indices of interest

        Returns
        ----------
        subdomains : ndarray
            Array of subdomain IDs indexed in the order of the indices

        Raises
        ------
        ValueError
            When one of the indices is not a valid subdomain index.

        """

        if indices != 'all':
            cv.check_type('offsets', indices, Iterable, Integral)

        if indices == 'all':
            indices = self.get_subdomain_indices()

        subdomains = np.zeros(len(indices), dtype=np.int64)
        keys = self.subdomain_indices.keys()
        values = self.subdomain_indices.values()

        for i, index in enumerate(indices):
            if index in values:
                subdomains[i] = keys[values.index(index)]
            else:
                msg = 'Unable to get subdomain for index "{0}" since it ' \
                      'is not a valid index'.format(index)
                raise ValueError(msg)

        return subdomains

    def get_xs(self, groups='all', subdomains='all', value='mean'):
        """

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

        filters = []
        filter_bins = []

        # Construct a collection of the domain filter bins
        if subdomains != 'all':
            cv.check_value('subdomains', subdomains, Iterable, Integral)
            filters.append(self.domain_type)
            filter_bins.append(tuple(subdomains))

        if groups != 'all':
            cv.check_value('groups', groups, Iterable, Integral)
            filters.append('energy')
            filter_bins.append(self.energy_groups.get_group_bounds(groups))

        # Query the multi-group cross-section tally for the data
        xs = self.xs_tally.get_values(filters=filters,
                                      filter_bins=filter_bins, value=value)
        return xs

    def get_condensed_xs(self, coarse_groups):
        """This routine takes in a collection of 2-tuples of energy groups"""

        cv.check_value('coarse groups', coarse_groups, EnergyGroups)

        # FIXME: this should use the Tally.slice(...) routine

    def get_subdomain_avg_xs(self, subdomains='all'):
        """Construct a subdomain-averaged version of this cross-section.

        Parameters
        ----------
        subdomains : Iterable of Integral or 'all'
            The subdomain IDs to average across

        Returns
        -------
        MultiGroupXS
            This MultiGroupXS averaged across subdomains of interest

        """

        if subdomains != 'all':
            cv.check_value('subdomains', subdomains, Iterable, Integral)

        avg_xs = copy.deepcopy(self)

        if self.domain_type == 'distribcell':
            avg_xs.domain_type = 'cell'
            avg_xs._subdomain_indices = {}
            avg_xs._offset = 0

            # Spatially average each tally
            for key, old_tally in avg_xs.tallies.items():
                # FIXME: Need to create Tally.mean(...)
                slice_tally = old_tally.slice(filters=[avg_xs.domain_type],
                                              filter_bins=[tuple(subdomains)])
                avg_tally = slice_tally.mean(filters=[avg_xs.domain_type],
                                             filter_bins=[tuple(subdomains)])
                avg_xs.tallies[key] = avg_tally

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
            cv.check_value('subdomains', subdomains, Iterable, Integral)

        string = 'Multi-Group XS\n'
        string += '{0: <16}=\t{1}\n'.format('\tType', self.xs_type)
        string += '{0: <16}=\t{1}\n'.format('\tDomain Type', self.domain_type)
        string += '{0: <16}=\t{1}\n'.format('\tDomain ID', self.domain.id)

        if self.xs_tally is not None:
            if subdomains == 'all':
                subdomains = self.get_subdomain_indices()

            # Loop over all subdomains
            for subdomain in subdomains:

                if self.domain_type == 'distribcell':
                    string += '{0: <16}=\t{1}\n'.format('\tSubDomain', subdomain)

                string += '{0: <16}\n'.format('\tCross-Sections [cm^-1]:')

                # Loop over energy groups ranges
                for group in range(1,self.num_groups+1):
                    bounds = self.energy_groups.get_group_bounds(group)
                    string += '{0: <12}Group {1} [{2: <10} - ' \
                         '{3: <10}MeV]:\t'.format('', group, bounds[0], bounds[1])
                    average = self.get_xs([group], [subdomain], 'mean')
                    rel_err = self.get_xs([group], [subdomain], 'rel_err')*100.
                    string += '{:.2e}+/-{:1.2e}%'.format(average, rel_err)
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
        xs_results['offset'] = self._offset
        xs_results['subdomain_indices'] = self._subdomain_indices

        # Pickle the MultiGroupXS results to a binary file
        filename = directory + '/' + filename + '.pkl'
        filename = filename.replace(' ', '-')
        pickle.dump(xs_results, open(filename, 'wb'))

    def restore_from_file(self, filename='mgxs', directory='mgxs'):
        """Restore the MultiGroupXS from a pickled binary file.

        Parameters
        ----------
        filename : str
            Filename for the pickled binary file (default is 'mgxs')
        directory : str
            Directory for the pickled binary file (default is 'mgxs')
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
        self.xs_type = xs_results['xs_type']
        self.domain_type = xs_results['domain_type']
        self.domain = xs_results['domain']
        self.energy_groups = xs_results['energy_groups']
        self.tallies = xs_results['tallies']
        self.xs_tally = xs_results['xs_tally']
        self._offset = xs_results['offset']
        self._subdomain_indices = xs_results['subdomain_indices']

    def export_xs_data(self, subdomains='all', filename='mgxs',
                       directory='mgxs', format='hdf5', append=True):
        """Export the multi-group cross-secttion data to a file.

        This routine leverages the functionality in the Pandas library to
        export DataFrames to CSV, HDF5, LaTeX and PDF files.

        Parameters
        ----------
        subdomains : Iterable of Integral or 'all'
        filename : str
            Filename for the exported file (default is 'mgxs')
        directory : str
            Directory for the exported file (default is 'mgxs')
        format : {'csv', 'hdf5', 'latex', 'pdf'}
            The format for the exported data file
        append : bool
            If True (default), appends to an existing file if possible
        """

        if subdomains != 'all':
            cv.check_type('submdomains', subdomains, Iterable, Integral)
        cv.check_type('filename', filename, basestring)
        cv.check_type('directory', directory, basestring)
        cv.check_values('format', format, ['hdf5', 'pickle'])
        cv.check_type('append', append, bool)

        # Make directory if it does not exist
        if not os.path.exists(directory):
            os.makedirs(directory)

        # FIXME: Use pandas dataframes!!

class TotalXS(MultiGroupXS):

    def __init__(self, name='', domain=None, domain_type=None, groups=None):
        super(TotalXS, self).__init__(name, domain, domain_type, groups)
        self.xs_type = 'total'

    def create_tallies(self):

        # Create a list of scores for each Tally to be created
        scores = ['flux', 'total']
        estimator = 'tracklength'
        keys = scores

        # Create the non-domain specific Filters for the Tallies
        group_edges = self.energy_groups.group_edges
        energy_filter = openmc.Filter('energy', group_edges)
        filters = [[energy_filter], [energy_filter]]

        # Intialize the Tallies
        super(TotalXS, self)._create_tallies(scores, filters, keys, estimator)

    def compute_xs(self):
        self.xs_tally = self.tallies['total'] / self.tallies['flux']


class TransportXS(MultiGroupXS):

    def __init__(self, name='', domain=None, domain_type=None, groups=None):
        super(TransportXS, self).__init__(name, domain, domain_type, groups)
        self.xs_type = 'transport'

    def create_tallies(self):

        # Create a list of scores for each Tally to be created
        scores = ['flux', 'total', 'scatter-1']
        estimator = 'analog'
        keys = scores

        # Create the non-domain specific Filters for the Tallies
        group_edges = self.energy_groups.group_edges
        energy_filter = openmc.Filter('energy', group_edges)
        energyout_filter = openmc.Filter('energyout', group_edges)
        filters = [[energy_filter], [energy_filter], [energyout_filter]]

        # Initialize the Tallies
        super(TransportXS, self)._create_tallies(scores, filters, keys, estimator)

    def compute_xs(self):
        self.xs_tally = self.tallies['total'] - self.tallies['scatter-1']
        self.xs_tally /= self.tallies['flux']


class AbsorptionXS(MultiGroupXS):

    def __init__(self, name='', domain=None, domain_type=None, groups=None):
        super(AbsorptionXS, self).__init__(name, domain, domain_type, groups)
        self.xs_type = 'absorption'

    def create_tallies(self):

        # Create a list of scores for each Tally to be created
        scores = ['flux', 'absorption']
        estimator = 'tracklength'
        keys = scores

        # Create the non-domain specific Filters for the Tallies
        group_edges = self.energy_groups.group_edges
        energy_filter = openmc.Filter('energy', group_edges)
        filters = [[energy_filter], [energy_filter]]

        # Intialize the Tallies
        super(AbsorptionXS, self)._create_tallies(scores, filters, keys, estimator)

    def compute_xs(self):
        self.xs_tally = self.tallies['absorption'] / self.tallies['flux']


class CaptureXS(MultiGroupXS):

    def __init__(self, name='', domain=None, domain_type=None, groups=None):
        super(CaptureXS, self).__init__(name, domain, domain_type, groups)
        self._xs_type = 'capture'

    def create_tallies(self):

        # Create a list of scores for each Tally to be created
        scores = ['flux', 'absorption', 'fission']
        estimator = 'tracklength'
        keys = scores

        # Create the non-domain specific Filters for the Tallies
        group_edges = self.energy_groups.group_edges
        energy_filter = openmc.Filter('energy', group_edges)
        filters = [[energy_filter], [energy_filter], [energy_filter]]

        # Intialize the Tallies
        super(CaptureXS, self)._create_tallies(scores, filters, keys, estimator)

    def compute_xs(self):
        self.xs_tally = self.tallies['absorption'] - self.tallies['fission']
        self.xs_tally /= self.tallies['flux']


class FissionXS(MultiGroupXS):

    def __init__(self, name='', domain=None, domain_type=None, energy_groups=None):
        super(FissionXS, self).__init__(name, domain, domain_type, energy_groups)
        self._xs_type = 'fission'

    def create_tallies(self):

        # Create a list of scores for each Tally to be created
        scores = ['flux', 'fission']
        estimator = 'tracklength'
        keys = scores

        # Create the non-domain specific Filters for the Tallies
        group_edges = self._energy_groups._group_edges
        energy_filter = openmc.Filter('energy', group_edges)
        filters = [[energy_filter], [energy_filter]]

        # Intialize the Tallies
        super(FissionXS, self)._create_tallies(scores, filters, keys, estimator)

    def compute_xs(self):
        self.xs_tally = self.tallies['fission'] / self.tallies['flux']


class NuFissionXS(MultiGroupXS):

    def __init__(self, name='', domain=None, domain_type=None, groups=None):
        super(NuFissionXS, self).__init__(name, domain, domain_type, groups)
        self._xs_type = 'nu-fission'

    def create_tallies(self):

        # Create a list of scores for each Tally to be created
        scores = ['flux', 'nu-fission']
        estimator = 'tracklength'
        keys = scores

        # Create the non-domain specific Filters for the Tallies
        group_edges = self.energy_groups.group_edges
        energy_filter = openmc.Filter('energy', group_edges)
        filters = [[energy_filter], [energy_filter]]

        # Intialize the Tallies
        super(NuFissionXS, self)._create_tallies(scores, filters, keys, estimator)

    def compute_xs(self):
        self.xs_tally = self.tallies['nu-fission'] / self.tallies['flux']


class ScatterXS(MultiGroupXS):

    def __init__(self, name='', domain=None, domain_type=None, energy_groups=None):
        super(ScatterXS, self).__init__(name, domain, domain_type, energy_groups)
        self._xs_type = 'scatter'

    def create_tallies(self):

        # Create a list of scores for each Tally to be created
        scores = ['flux', 'scatter']
        estimator = 'tracklength'
        keys = scores

        # Create the non-domain specific Filters for the Tallies
        group_edges = self.energy_groups.group_edges
        energy_filter = openmc.Filter('energy', group_edges)
        filters = [[energy_filter], [energy_filter]]

        # Intialize the Tallies
        super(ScatterXS, self)._create_tallies(scores, filters, keys, estimator)

    def compute_xs(self):
        self.xs_tally = self.tallies['scatter'] / self.tallies['flux']


class NuScatterXS(MultiGroupXS):

    def __init__(self, name='', domain=None, domain_type=None, groups=None):
        super(NuScatterXS, self).__init__(name, domain, domain_type, groups)
        self._xs_type = 'nu-scatter'

    def create_tallies(self):

        # Create a list of scores for each Tally to be created
        scores = ['flux', 'nu-scatter']
        estimator = 'analog'
        keys = scores

        # Create the non-domain specific Filters for the Tallies
        group_edges = self.energy_groups.group_edges
        energy_filter = openmc.Filter('energy', group_edges)
        filters = [[energy_filter], [energy_filter]]

        # Intialize the Tallies
        super(NuScatterXS, self)._create_tallies(scores, filters, keys, estimator)

    def compute_xs(self):
        self.xs_tally = self.tallies['nu-scatter'] / self.tallies['flux']


class ScatterMatrixXS(MultiGroupXS):

    def __init__(self, name='', domain=None, domain_type=None, groups=None):
        super(ScatterMatrixXS, self).__init__(name, domain, domain_type, groups)
        self._xs_type = 'scatter matrix'

    def create_tallies(self):

        # Create a list of scores for each Tally to be created
        scores = ['flux', 'scatter', 'scatter-1']
        estimator = 'analog'
        keys = scores

        # Create the non-domain specific Filters for the Tallies
        group_edges = self.energy_groups.group_edges
        energy_filter = openmc.Filter('energy', group_edges)
        energyout_filter = openmc.Filter('energyout', group_edges)
        filters = [[energy_filter], [energy_filter, energyout_filter], [energyout_filter]]

        # Intialize the Tallies
        super(ScatterMatrixXS, self)._create_tallies(scores, filters, keys, estimator)

    def compute_xs(self):
        self.xs_tally = self.tallies['scatter'] - self.tallies['scatter-1']
        self.xs_tally /= self.tallies['flux']

    def get_condensed_xs(self, coarse_groups):
        """This routine takes in a collection of 2-tuples of energy groups"""

        cv.check_value('coarse groups', coarse_groups, EnergyGroups)

        # FIXME: this should use the Tally.slice(...) routine

        # Error checking for the group bounds is done here
        new_groups = self.energy_groups.getCondensedGroups(coarse_groups)
        num_coarse_groups = new_groups._num_groups

    def get_xs(self, in_groups='all', out_groups='all',
              subdomains='all', value='mean'):

        if self.xs_tally is None:
            msg = 'Unable to get cross-section since it has not been computed'
            raise ValueError(msg)

        cv.check_value('value', value, ['mean', 'std. dev.', 'rel. err.'])
        if in_groups != 'all':
            cv.check_value('in groups', in_groups, Iterable, Integral)
        if out_groups != 'all':
            cv.check_value('out groups', out_groups, Iterable, Integral)
        if subdomains != 'all':
            cv.check_value('subdomains', subdomains, Iterable, Integral)

        # FIXME: Make this use Tally.get_values()

    def print_xs(self, subdomains='all'):

        if subdomains != 'all':
            cv.check_value('subdomains', subdomains, Iterable, Integral)

        string = 'Multi-Group XS\n'
        string += '{0: <16}{1}{2}\n'.format('\tType', '=\t', self.xs_type)
        string += '{0: <16}{1}{2}\n'.format('\tDomain Type', '=\t', self.domain_type)
        string += '{0: <16}{1}{2}\n'.format('\tDomain ID', '=\t', self.domain.id)

        string += '{0: <16}\n'.format('\tEnergy Groups:')

        # Loop over energy groups ranges
        for group in range(1,self.num_groups+1):
            bounds = self.energy_groups.get_group_bounds(group)
            string += '{0: <12}Group {1} [{2: <10} - ' \
                      '{3: <10}MeV]\n'.format('', group, bounds[0], bounds[1])

        if subdomains == 'all':
            subdomains = self._subdomain_indices.keys()

        for subdomain in subdomains:

            if self.domain_type == 'distribcell':
                string += '{0: <16}{1}{2}\n'.format('\tSubDomain', '=\t', subdomain)

            string += '{0: <16}\n'.format('\tCross-Sections [cm^-1]:')

            # Loop over energy groups ranges
            for in_group in range(1,self.num_groups+1):
                for out_group in range(1,self.num_groups+1):
                    string += '{0: <12}Group {1} -> Group {2}:\t\t'.format('', in_group, out_group)
                    average = self.get_xs([in_group], [out_group], [subdomain], 'mean')
                    rel_err = self.get_xs([in_group], [out_group], [subdomain], 'rel. err.')
                    string += '{:.2e}+/-{:1.2e}%'.format(average[0,0,0], rel_err[0,0,0])
                    string += '\n'
                string += '\n'
        print(string)


class NuScatterMatrixXS(ScatterMatrixXS):

    def __init__(self, name='', domain=None, domain_type=None, groups=None):
        super(NuScatterMatrixXS, self).__init__(name, domain, domain_type, groups)
        self.xs_type = 'nu-scatter matrix'

    def create_tallies(self):

        # Create a list of scores for each Tally to be created
        scores = ['flux', 'nu-scatter', 'scatter-1']
        estimator = 'analog'
        keys = scores

        # Create the non-domain specific Filters for the Tallies
        group_edges = self.energy_groups.group_edges
        energy_filter = openmc.Filter('energy', group_edges)
        energyout_filter = openmc.Filter('energyout', group_edges)
        filters = [[energy_filter], [energy_filter, energyout_filter], [energyout_filter]]

        # Intialize the Tallies
        super(ScatterMatrixXS, self)._create_tallies(scores, filters, keys, estimator)

    def compute_xs(self):
        self.xs_tally = self.tallies['nu-scatter'] - self.tallies['scatter-1']
        self.xs_tally /= self.tallies['flux']


class Chi(MultiGroupXS):

    def __init__(self, name='', domain=None, domain_type=None, groups=None):
        super(Chi, self).__init__(name, domain, domain_type, groups)
        self._xs_type = 'chi'

    def create_tallies(self):

        # Create a list of scores for each Tally to be created
        scores = ['nu-fission', 'nu-fission']
        estimator = 'analog'
        keys = ['nu-fission-in', 'nu-fission-out']

        # Create the non-domain specific Filters for the Tallies
        group_edges = self._energy_groups._group_edges
        energy_filter = openmc.Filter('energy', group_edges)
        energyout_filter = openmc.Filter('energyout', group_edges)
        filters = [[energy_filter], [energyout_filter]]

        # Intialize the Tallies
        super(Chi, self)._create_tallies(scores, filters, keys, estimator)

    def compute_xs(self):

    # Extract and clean the Tally data
    tally_data, zero_indices = super(Chi, self).getAllTallyData()
    nu_fission_in = tally_data['nu-fission-in']
    nu_fission_out = tally_data['nu-fission-out']

    # Set any zero reaction rates to -1
    nu_fission_in[0, zero_indices['nu-fission-in']] = -1.

    # FIXME - uncertainty propagation
    self._xs_tally = infermc.error_prop.arithmetic.divide_by_scalar(nu_fission_out,
                                   nu_fission_in.sum(2)[0, :, np.newaxis, ...],
                                   corr, False)

    # Compute the total across all groups per subdomain
    norm = self._xs_tally.sum(2)[0, :, np.newaxis, ...]

    # Set any zero norms (in non-fissionable domains) to -1
    norm_indices = norm == 0.
    norm[norm_indices] = -1.

    # Normalize chi to 1.0
    # FIXME - uncertainty propagation
    self._xs_tally = infermc.error_prop.arithmetic.divide_by_scalar(self._xs_tally, norm,
                                                              corr, False)

    # For any region without flux or reaction rate, convert xs to zero
    self._xs_tally[:, norm_indices] = 0.

    # FIXME - uncertainty propagation - this is just a temporary fix
    self._xs_tally[1, ...] = 0.

    # Correct -0.0 to +0.0
    self._xs_tally += 0.