from collections import Iterable
from numbers import Integral, Real
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
GREEK['diffusion'] = '$D$'


def flip_axis(arr, axis=0):
    """Flip contents of `axis` in array 'arr'
    Taken verbatim from:
    https://github.com/nipy/nibabel/blob/master/nibabel/orientations.py
    """
    arr = np.asanyarray(arr)
    arr = arr.swapaxes(0, axis)
    arr = np.flipud(arr)
    return arr.swapaxes(axis, 0)


class MultiGroupXS(object):
    """

    """

    # This is an abstract class which cannot be instantiated
    metaclass__ = abc.ABCMeta

    def __init__(self, name='', domain=None,
                 domain_type=None, energy_groups=None):
        """
        :param name:
        :param domain:
        :param domain_type:
        :param energy_groups:
        :return:
        """

        self._name = ''
        self._xs_type = None
        self._domain = None
        self._domain_type = None
        self._energy_groups = None
        self._num_groups = None
        self._tallies = dict()
        self._xs = None

        # A dictionary used to compute indices into the xs array
        # Keys   - Domain ID (ie, Material ID, Region ID for districell, etc)
        # Values - Offset/stride into xs array
        self._subdomain_offsets = dict()
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

        # If this is the first time we have tried to copy this object, create a copy
        if existing is None:
            clone = type(self).__new__(type(self))
            clone._name = self._name
            clone._xs_type = self._xs_type
            clone._domain = self._domain
            clone._domain_type = self._domain_type
            clone._energy_groups = copy.deepcopy(self._energy_groups, memo)
            clone._num_groups = self._num_groups
            clone._xs = copy.deepcopy(self._xs, memo)
            clone._subdomain_offsets = copy.deepcopy(self._subdomain_offsets, memo)
            clone._offset = copy.deepcopy(self._offset, memo)

            clone._tallies = dict()
            for tally_type, tally in self._tallies.items():
                clone._tallies[tally_type] = copy.deepcopy(tally, memo)

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

    @name.setter
    def name(self, name):
        cv.check_type('MultiGroupXS name', name, basestring)
        self._name = name

    @domain.setter
    def domain(self, domain):
        cv.check_type('MultiGroupXS domain', domain, DOMAINS)
        self._domain = domain
        if self._domain_type in ['material', 'cell', 'universe', 'mesh']:
            self._subdomain_offsets[domain.id] = 0

    @energy_groups.setter
    def energy_groups(self, energy_groups):
        cv.check_type('MultiGroupXS energy groups', energy_groups,
                      openmc.mgxs.EnergyGroups)
        self._energy_groups = energy_groups
        self._num_groups = energy_groups._num_groups

    @domain_type.setter
    def domain_type(self, domain_type):
        cv.check_type('MultiGroupXS domain type', domain_type, DOMAIN_TYPES)
        self._domain_type = domain_type

    def find_domain_offset(self):
        tally = self.tallies[self.tallies.keys()[0]]
        filter = tally.find_filter(self.domain_type, [self.domain.id])
        self._offset = filter.offset

    def set_subdomain_offset(self, domain_id, offset):
        """
        :param domain_id:
        :param offset:
        :return:
        """

        cv.check_type('subdomain id', domain_id, Integral)
        cv.check_type('subdomain offset', offset, Integral)
        self._subdomain_offsets[domain_id] = offset

    @abc.abstractmethod
    def _create_tallies(self, scores, filters, keys, estimator):
        """

        :param scores:
        :param filters:
        :param keys:
        :param estimator:
        :return:
        """

        if self.energy_groups is None:
            raise ValueError('Unable to create Tallies without energy groups')
        elif self.domain is None:
             raise ValueError('Unable to create Tallies without a domain')
        elif self.domain_type is None:
             raise ValueError('Unable to create Tallies without a domain type')

        cv.check_type('scores', scores, Iterable, basestring)
        cv.check_value('scores', scores, openmc.SCORE_TYPES)
        cv.check_type('filters', scores, Iterable, openmc.Filter)
        cv.check_type('keys', keys, Iterable, basestring)
        cv.check_value('# scores', len(scores), len(keys))
        cv.check_type('estimator', estimator, basestring)
        cv.check_value('estimator', estimator, ['analog', 'tracklength'])

        # Create a domain Filter object
        domain_filter = openmc.Filter(self.domain_type, self.domain.id)

        for score, key, filters in zip(scores, keys, filters):
            self.tallies[key] = openmc.Tally(name=self.name)
            self.tallies[key].add_score(score)
            self.tallies[key].estimator = estimator
            self.tallies[key].add_filter(domain_filter)

            # Add all non-domain specific Filters (ie, energy) to the Tally
            for filter in filters:
                self.tallies[key].add_filter(filter)

    def get_subdomain_offsets(self, subdomains='all'):
        """

        :param subdomains:
        :return:
        """

        if subdomains != 'all':
            cv.check_type('subdomains', subdomains, Iterable, Integral)

        if subdomains == 'all':
            offsets = np.arange(self.xs.shape[1])
        else:
            offsets = np.zeros(len(subdomains), dtype=np.int64)

            for i, subdomain in enumerate(subdomains):
                if subdomain in self._subdomain_offsets:
                    offsets[i] = self._subdomain_offsets[subdomain]
                else:
                    msg = 'Unable to get index for subdomain "{0}" since it ' \
                          'is not a subdomain in cross-section'.format(subdomain)
                    raise ValueError(msg)

        return offsets

    def get_subdomains(self, offsets='all'):

        if offsets != 'all':
            cv.check_type('offsets', offsets, Iterable, Integral)

        if offsets == 'all':
            offsets = self.get_subdomain_offsets()

        subdomains = np.zeros(len(offsets), dtype=np.int64)
        keys = self._subdomain_offsets.keys()
        values = self._subdomain_offsets.values()

        for i, offset in enumerate(offsets):
            if offset in values:
                subdomains[i] = keys[values.index(offset)]
            else:
                msg = 'Unable to get subdomain for offset "{0}" since it ' \
                      'is not an offset in the cross-section'.format(offset)
                raise ValueError(msg)

        return subdomains

    def get_xs(self, groups='all', subdomains='all', metric='mean'):

        if self.xs is None:
            msg = 'Unable to get cross-section since it has not been computed'
            raise ValueError(msg)

        cv.check_value('metric', metric, ['mean', 'std. dev.', 'rel. err.'])
        if groups != 'all':
            cv.check_value('groups', groups, Iterable, Integral)
        if subdomains != 'all':
            cv.check_value('subdomains', subdomains, Iterable, Integral)

        # FIXME: Make this use Tally.get_values()

    def get_condensed_xs(self, coarse_groups):
        """This routine takes in a collection of 2-tuples of energy groups"""

        cv.check_value('coarse groups', coarse_groups, EnergyGroups)

        # FIXME: this should use the Tally.slice(...) routine

    def get_domain_vg_xs(self, subdomains='all'):

        if self.domain_type != 'distribcell':
            msg = 'Unable to compute domain averaged "{0}" xs for "{1}"' \
                  '"{2}" since it is not a distribcell'.format(self._xs_type,
                  self._domain_type, self._domain.id)
            raise ValueError(msg)

        if subdomains != 'all':
            cv.check_value('subdomains', subdomains, Iterable, Integral)

        # FIXME: This should use tally arithmetic

    def print_xs(self, subdomains='all'):

        if subdomains != 'all':
            cv.check_value('subdomains', subdomains, Iterable, Integral)

        string = 'Multi-Group XS\n'
        string += '{0: <16}{1}{2}\n'.format('\tType', '=\t', self.xs_type)
        string += '{0: <16}{1}{2}\n'.format('\tDomain Type', '=\t', self.domain_type)
        string += '{0: <16}{1}{2}\n'.format('\tDomain ID', '=\t', self.domain.id)

        if subdomains == 'all':
            subdomains = self._subdomain_offsets.keys()

        # Loop over all subdomains
        for subdomain in subdomains:

            if self.domain_type == 'distribcell':
                string += '{0: <16}{1}{2}\n'.format('\tSubDomain', '=\t', subdomain)

            string += '{0: <16}\n'.format('\tCross-Sections [cm^-1]:')

            # Loop over energy groups ranges
            for group in range(1,self.num_groups+1):
                bounds = self._energy_groups.getGroupBounds(group)
                string += '{0: <12}Group {1} [{2: <10} - ' \
                          '{3: <10}MeV]:\t'.format('', group, bounds[0], bounds[1])
                average = self.get_xs([group], [subdomain], 'mean')
                rel_err = self.get_xs([group], [subdomain], 'rel. err.')
                string += '{:.2e}+/-{:1.2e}%'.format(average[0,0,0], rel_err[0,0,0])
                string += '\n'
            string += '\n'

        print(string)

    def dump_to_file(self, filename='multigroupxs', directory='multigroupxs'):

        cv.check_type('filename', filename, basestring)
        cv.check_type('directory', directory, basestring)

        # Make directory if it does not exist
        if not os.path.exists(directory):
            os.makedirs(directory)

        # Create an empty dictionary to store the data
        xs_results = dict()

        # Store all of this MultiGroupXS' class attributes in the dictionary
        xs_results['name'] = self.name
        xs_results['xs type'] = self.xs_type
        xs_results['domain type'] = self.domain_type
        xs_results['domain'] = self.domain
        xs_results['energy groups'] = self.energy_groups
        xs_results['tallies'] = self.tallies
        xs_results['xs'] = self.xs
        xs_results['offset'] = self._offset
        xs_results['subdomain offsets'] = self._subdomain_offsets

        # Pickle the MultiGroupXS results to a file
        filename = directory + '/' + filename + '.pkl'
        filename = filename.replace(' ', '-')
        pickle.dump(xs_results, open(filename, 'wb'))

    def restore_from_file(self, filename='multigroupxs', directory='multigroupxs'):

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
        self.xs_type = xs_results['xs type']
        self.domain_type = xs_results['domain type']
        self.domain = xs_results['domain']
        self.energy_groups = xs_results['energy groups']
        self.tallies = xs_results['tallies']
        self.xs = xs_results['xs']
        self._offset = xs_results['offset']
        self._subdomain_offsets = xs_results['subdomain offsets']

    def exportResults(self, subdomains='all', filename='multigroupxs',
                      directory='multigroupxs', format='hdf5', append=True):

        if subdomains != 'all':
            cv.check_type('submdomains', subdomains, Iterable, Integral)
        cv.check_type('filename', filename, basestring)
        cv.check_type('directory', directory, basestring)
        cv.check_values('format', format, ['hdf5', 'pickle'])
        cv.check_type('append', append, bool)

        # Make directory if it does not exist
        if not os.path.exists(directory):
            os.makedirs(directory)

        # FIXME: Use tally arithmetic!!!

    def print_pdf(self, subdomains='all', filename='multigroupxs',
                 directory='multigroupxs'):

        if subdomains != 'all':
            cv.check_type('submdomains', subdomains, Iterable, Integral)
        cv.check_type('filename', filename, basestring)
        cv.check_type('directory', directory, basestring)

        # Make directory if it does not exist
        if not os.path.exists(directory):
            os.makedirs(directory)

        filename = filename.replace(' ', '-')

        # Generate LaTeX file
        self.exportResults(subdomains, filename, '.', 'latex', False)

        # Compile LaTeX to PDF
        FNULL = open(os.devnull, 'w')
        subprocess.check_call('pdflatex {0}.tex'.format(filename),
                              shell=True, stdout=FNULL)

        # Move PDF to requested directory and cleanup temporary LaTeX files
        if directory != '.':
            os.system('mv {0}.pdf {1}'.format(filename, directory))
        os.system('rm {0}.tex {0}.aux {0}.log'.format(filename))


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
        self.xs = self.tallies['total'] / self.tallies['flux']


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
        self.xs = self.tallies['total'] - self.tallies['scatter-1']
        self.xs /= self.tallies['flux']


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
        self.xs = self.tallies['absorption'] / self.tallies['flux']


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
        self.xs = self.tallies['absorption'] - self.tallies['fission']
        self.xs /= self.tallies['flux']


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
        self.xs = self.tallies['fission'] / self.tallies['flux']


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
        self.xs = self.tallies['nu-fission'] / self.tallies['flux']


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
        self.xs = self.tallies['scatter'] / self.tallies['flux']


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
        self.xs = self.tallies['nu-scatter'] / self.tallies['flux']


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
        self.xs = self.tallies['scatter'] - self.tallies['scatter-1']
        self.xs /= self.tallies['flux']

    def get_condensed_xs(self, coarse_groups):
        """This routine takes in a collection of 2-tuples of energy groups"""

        cv.check_value('coarse groups', coarse_groups, EnergyGroups)

        # FIXME: this should use the Tally.slice(...) routine

        # Error checking for the group bounds is done here
        new_groups = self.energy_groups.getCondensedGroups(coarse_groups)
        num_coarse_groups = new_groups._num_groups

    def get_xs(self, in_groups='all', out_groups='all',
              subdomains='all', metric='mean'):

        if self.xs is None:
            msg = 'Unable to get cross-section since it has not been computed'
            raise ValueError(msg)

        cv.check_value('metric', metric, ['mean', 'std. dev.', 'rel. err.'])
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
            subdomains = self._subdomain_offsets.keys()

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
        self.xs = self.tallies['nu-scatter'] - self.tallies['scatter-1']
        self.xs /= self.tallies['flux']


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
    self._xs = infermc.error_prop.arithmetic.divide_by_scalar(nu_fission_out,
                                   nu_fission_in.sum(2)[0, :, np.newaxis, ...],
                                   corr, False)

    # Compute the total across all groups per subdomain
    norm = self._xs.sum(2)[0, :, np.newaxis, ...]

    # Set any zero norms (in non-fissionable domains) to -1
    norm_indices = norm == 0.
    norm[norm_indices] = -1.

    # Normalize chi to 1.0
    # FIXME - uncertainty propagation
    self._xs = infermc.error_prop.arithmetic.divide_by_scalar(self._xs, norm,
                                                              corr, False)

    # For any region without flux or reaction rate, convert xs to zero
    self._xs[:, norm_indices] = 0.

    # FIXME - uncertainty propagation - this is just a temporary fix
    self._xs[1, ...] = 0.

    # Correct -0.0 to +0.0
    self._xs += 0.