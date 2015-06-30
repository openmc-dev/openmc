from collections import Iterable
import copy
import os
import pickle
import itertools
from numbers import Integral
from xml.etree import ElementTree as ET
import sys

import numpy as np

from openmc import Mesh, Filter, Trigger, Nuclide
from openmc.summary import Summary
from openmc.checkvalue import check_type, check_value, check_greater_than
from openmc.clean_xml import *


if sys.version_info[0] >= 3:
    basestring = str

# "Static" variable for auto-generated Tally IDs
AUTO_TALLY_ID = 10000


def reset_auto_tally_id():
    global AUTO_TALLY_ID
    AUTO_TALLY_ID = 10000


class Tally(object):
    """A tally defined by a set of scores that are accumulated for a list of
    nuclides given a set of filters.

    Parameters
    ----------
    tally_id : int, optional
        Unique identifier for the tally. If none is specified, an identifier
        will automatically be assigned
    name : str, optional
        Name of the tally. If not specified, the name is the empty string.

    Attributes
    ----------
    id : int
        Unique identifier for the tally
    name : str
        Name of the tally
    filters : list of openmc.filter.Filter
        List of specified filters for the tally
    nuclides : list of openmc.nuclide.Nuclide
        List of nuclides to score results for
    scores : list of str
        List of defined scores, e.g. 'flux', 'fission', etc.
    estimator : {'analog', 'tracklength'}
        Type of estimator for the tally
    triggers : list of openmc.trigger.Trigger
        List of tally triggers
    num_score_bins : int
        Total number of scores, accounting for the fact that a single
        user-specified score, e.g. scatter-P3 or flux-Y2,2, might have multiple
        bins
    num_scores : int
        Total number of user-specified scores
    num_filter_bins : int
        Total number of filter bins accounting for all filters
    num_bins : int
        Total number of bins for the tally
    num_realizations : int
        Total number of realizations
    with_summary : bool
        Whether or not a Summary has been linked
    sum : ndarray
        An array containing the sum of each independent realization for each bin
    sum_sq : ndarray
        An array containing the sum of each independent realization squared for
        each bin
    mean : ndarray
        An array containing the sample mean for each bin
    std_dev : ndarray
        An array containing the sample standard deviation for each bin

    """

    def __init__(self, tally_id=None, name=''):
        # Initialize Tally class attributes
        self.id = tally_id
        self.name = name
        self._filters = []
        self._nuclides = []
        self._scores = []
        self._estimator = None
        self._triggers = []

        self._num_score_bins = 0
        self._num_realizations = 0
        self._with_summary = False

        self._sum = None
        self._sum_sq = None
        self._mean = None
        self._std_dev = None

    def __deepcopy__(self, memo):
        existing = memo.get(id(self))

        # If this is the first time we have tried to copy this object, create a copy
        if existing is None:
            clone = type(self).__new__(type(self))
            clone.id = self.id
            clone.name = self.name
            clone.estimator = self.estimator
            clone.num_score_bins = self.num_score_bins
            clone.num_realizations = self.num_realizations
            clone._sum = copy.deepcopy(self.sum, memo)
            clone._sum_sq = copy.deepcopy(self.sum_sq, memo)
            clone._mean = copy.deepcopy(self.mean, memo)
            clone._std_dev = copy.deepcopy(self.std_dev, memo)

            clone.filters = []
            for filter in self.filters:
                clone.add_filter(copy.deepcopy(filter, memo))

            clone.nuclides = []
            for nuclide in self.nuclides:
                clone.add_nuclide(copy.deepcopy(nuclide, memo))

            clone.scores = []
            for score in self.scores:
                clone.add_score(score)

            clone.triggers = []
            for trigger in self.triggers:
                clone.add_trigger(trigger)

            memo[id(self)] = clone

            return clone

        # If this object has been copied before, return the first copy made
        else:
            return existing

    def __eq__(self, tally2):
        # Check all filters
        if len(self.filters) != len(tally2.filters):
            return False

        for filter in self.filters:
            if filter not in tally2.filters:
                return False

        # Check all nuclides
        if len(self.nuclides) != len(tally2.nuclides):
            return False

        for nuclide in self.nuclides:
            if nuclide not in tally2.nuclides:
                return False

        # Check all scores
        if len(self.scores) != len(tally2.scores):
            return False

        for score in self.scores:
            if score not in tally2.scores:
                return False

        if self.estimator != tally2.estimator:
            return False

        return True

    def __hash__(self):
        hashable = []

        for filter in self.filters:
            hashable.append((filter.type, tuple(filter.bins)))

        for nuclide in self.nuclides:
            hashable.append(nuclide.name)

        for score in self.scores:
            hashable.append(score)

        hashable.append(self.estimator)
        hashable.append(self.name)

        return hash(tuple(hashable))

    def __add__(self, other):
        # FIXME: Error checking: must check that results has been
        # set and that # bins is the same

        new_tally = Tally()
        new_tally._mean = self._mean + other._mean
        new_tally._std_dev = np.sqrt(self.std_dev**2 + other.std_dev**2)

    @property
    def id(self):
        return self._id

    @property
    def name(self):
        return self._name

    @property
    def filters(self):
        return self._filters

    @property
    def nuclides(self):
        return self._nuclides

    @property
    def num_nuclides(self):
        return len(self._nuclides)

    @property
    def scores(self):
        return self._scores

    @property
    def num_scores(self):
        return len(self._scores)

    @property
    def num_score_bins(self):
        return self._num_score_bins

    @property
    def num_filter_bins(self):
        num_bins = 1

        for filter in self._filters:
            num_bins *= filter.num_bins

        return num_bins

    @property
    def num_bins(self):
        num_bins = self.num_filter_bins
        num_bins *= self.num_nuclides
        num_bins *= self.num_score_bins
        return num_bins

    @property
    def estimator(self):
        return self._estimator

    @property
    def triggers(self):
        return self._triggers

    @property
    def num_realizations(self):
        return self._num_realizations

    @property
    def with_summary(self):
        return self._with_summary

    @property
    def sum(self):
        return self._sum

    @property
    def sum_sq(self):
        return self._sum_sq

    @property
    def mean(self):
        return self._mean

    @property
    def std_dev(self):
        return self._std_dev

    @estimator.setter
    def estimator(self, estimator):
        check_value('estimator', estimator, ['analog', 'tracklength'])
        self._estimator = estimator

    def add_trigger(self, trigger):
        """Add a tally trigger to the tally

        Parameters
        ----------
        trigger : openmc.trigger.Trigger
            Trigger to add

        """

        if not isinstance(trigger, Trigger):
            msg = 'Unable to add a tally trigger for Tally ID={0} to ' \
                  'since "{1}" is not a Trigger'.format(self.id, trigger)
            raise ValueError(msg)

        self._triggers.append(trigger)

    @id.setter
    def id(self, tally_id):
        if tally_id is None:
            global AUTO_TALLY_ID
            self._id = AUTO_TALLY_ID
            AUTO_TALLY_ID += 1
        else:
            check_type('tally ID', tally_id, Integral)
            check_greater_than('tally ID', tally_id, 0)
            self._id = tally_id

    @name.setter
    def name(self, name):
        check_type('tally name', name, basestring)
        self._name = name

    def add_filter(self, filter):
        """Add a filter to the tally

        Parameters
        ----------
        filter : openmc.filter.Filter
            Filter to add

        """

        if not isinstance(filter, Filter):
            msg = 'Unable to add Filter "{0}" to Tally ID={1} since it is ' \
                  'not a Filter object'.format(filter, self.id)
            raise ValueError(msg)

        self._filters.append(filter)

    def add_nuclide(self, nuclide):
        """Specify that scores for a particular nuclide should be accumulated

        Parameters
        ----------
        nuclide : openmc.nuclide.Nuclide
            Nuclide to add

        """

        self._nuclides.append(nuclide)

    def add_score(self, score):
        """Specify a quantity to be scored

        Parameters
        ----------
        score : str
            Score to be accumulated, e.g. 'flux'

        """

        if not isinstance(score, basestring):
            msg = 'Unable to add score "{0}" to Tally ID={1} since it is ' \
                  'not a string'.format(score, self.id)
            raise ValueError(msg)

        # If the score is already in the Tally, don't add it again
        if score in self.scores:
            return
        else:
            self._scores.append(score)

    @num_score_bins.setter
    def num_score_bins(self, num_score_bins):
        self._num_score_bins = num_score_bins

    @num_realizations.setter
    def num_realizations(self, num_realizations):
        check_type('number of realizations', num_realizations, Integral)
        check_greater_than('number of realizations', num_realizations, 0, True)
        self._num_realizations = num_realizations

    @with_summary.setter
    def with_summary(self, with_summary):
        check_type('with_summary', with_summary, bool)
        self._with_summary = with_summary

    @sum.setter
    def sum(self, sum):
        check_type('sum', sum, Iterable)
        self._sum = sum

    @sum_sq.setter
    def sum_sq(self, sum_sq):
        check_type('sum_sq', sum_sq, Iterable)
        self._sum_sq = sum_sq

    def remove_score(self, score):
        """Remove a score from the tally

        Parameters
        ----------
        score : str
            Score to remove

        """

        if score not in self.scores:
            msg = 'Unable to remove score "{0}" from Tally ID={1} since the ' \
                  'Tally does not contain this score'.format(score, self.id)
            ValueError(msg)

        self._scores.remove(score)

    def remove_filter(self, filter):
        """Remove a filter from the tally

        Parameters
        ----------
        filter : openmc.filter.Filter
            Filter to remove

        """

        if filter not in self.filters:
            msg = 'Unable to remove filter "{0}" from Tally ID={1} since the ' \
                  'Tally does not contain this filter'.format(filter, self.id)
            ValueError(msg)

        self._filters.remove(filter)

    def remove_nuclide(self, nuclide):
        """Remove a nuclide from the tally

        Parameters
        ----------
        nuclide : openmc.nuclide.Nuclide
            Nuclide to remove

        """

        if nuclide not in self.nuclides:
            msg = 'Unable to remove nuclide "{0}" from Tally ID={1} since the ' \
                  'Tally does not contain this nuclide'.format(nuclide, self.id)
            ValueError(msg)

        self._nuclides.remove(nuclide)

    def compute_std_dev(self, t_value=1.0):
        """Compute the sample mean and standard deviation for each bin in the tally

        Parameters
        ----------
        t_value : float, optional
            Student's t-value applied to the uncertainty. Defaults to 1.0,
            meaning the reported value is the sample standard deviation.

        """

        # Calculate sample mean and standard deviation
        self._mean = self.sum / self.num_realizations
        self._std_dev = np.sqrt((self.sum_sq / self.num_realizations -
                                 self.mean**2) / (self.num_realizations - 1))
        self._std_dev *= t_value

    def __repr__(self):
        string = 'Tally\n'
        string += '{0: <16}{1}{2}\n'.format('\tID', '=\t', self.id)
        string += '{0: <16}{1}{2}\n'.format('\tName', '=\t', self.name)

        string += '{0: <16}\n'.format('\tFilters')

        for filter in self.filters:
            string += '{0: <16}\t\t{1}\t{2}\n'.format('', filter.type,
                                                          filter.bins)

        string += '{0: <16}{1}'.format('\tNuclides', '=\t')

        for nuclide in self.nuclides:
            if isinstance(nuclide, Nuclide):
                string += '{0} '.format(nuclide.name)
            else:
                string += '{0} '.format(nuclide)

        string += '\n'

        string += '{0: <16}{1}{2}\n'.format('\tScores', '=\t', self.scores)
        string += '{0: <16}{1}{2}\n'.format('\tEstimator', '=\t', self.estimator)

        return string

    def can_merge(self, tally):
        """Determine if another tally can be merged with this one

        Parameters
        ----------
        tally : Tally
            Tally to check for merging

        """

        if not isinstance(tally, Tally):
            return False

        # Must have same estimator
        if self.estimator != tally.estimator:
            return False

        # Must have same nuclides
        if len(self.nuclides) != len(tally.nuclides):
            return False

        for nuclide in self.nuclides:
            if nuclide not in tally.nuclides:
                return False

        # Must have same or mergeable filters
        if len(self.filters) != len(tally.filters):
            return False

        # Look to see if all filters are the same, or one or more can be merged
        for filter1 in self.filters:
            mergeable_filter = False

            for filter2 in tally.filters:
                if filter1 == filter2 or filter1.can_merge(filter2):
                    mergeable_filter = True
                    break

            # If no mergeable filter was found, the tallies are not mergeable
            if not mergeable_filter:
                return False

        # Tallies are mergeable if all conditional checks passed
        return True

    def merge(self, tally):
        """Merge another tally with this one

        Parameters
        ----------
        tally : Tally
            Tally to merge with this one

        Returns
        -------
        merged_tally : Tally
            Merged tallies

        """

        if not self.can_merge(tally):
            msg = 'Unable to merge tally ID={0} with {1}'.format(tally.id, self.id)
            raise ValueError(msg)

        # Create deep copy of tally to return as merged tally
        merged_tally = copy.deepcopy(self)

        # Differentiate Tally with a new auto-generated Tally ID
        merged_tally.id = None

        # Merge filters
        for i, filter1 in enumerate(merged_tally.filters):
            for filter2 in tally.filters:
                if filter1 != filter2 and filter1.can_merge(filter2):
                    merged_filter = filter1.merge(filter2)
                    merged_tally.filters[i] = merged_filter
                    break

        # Add scores from second tally to merged tally
        for score in tally.scores:
            merged_tally.add_score(score)

        # Add triggers from second tally to merged tally
        for trigger in tally.triggers:
            merged_tally.add_trigger(trigger)

        return merged_tally

    def get_tally_xml(self):
        """Return XML representation of the tally

        Returns
        -------
        element : xml.etree.ElementTree.Element
            XML element containing tally data

        """

        element = ET.Element("tally")

        # Tally ID
        element.set("id", str(self.id))

        # Optional Tally name
        if self.name != '':
            element.set("name", self.name)

        # Optional Tally filters
        for filter in self.filters:
            subelement = ET.SubElement(element, "filter")
            subelement.set("type", str(filter.type))

            if filter.bins is not None:
                bins = ''
                for bin in filter.bins:
                    bins += '{0} '.format(bin)

                subelement.set("bins", bins.rstrip(' '))

        # Optional Nuclides
        if len(self.nuclides) > 0:
            nuclides = ''
            for nuclide in self.nuclides:
                if isinstance(nuclide, Nuclide):
                    nuclides += '{0} '.format(nuclide.name)
                else:
                    nuclides += '{0} '.format(nuclide)

            subelement = ET.SubElement(element, "nuclides")
            subelement.text = nuclides.rstrip(' ')

        # Scores
        if len(self.scores) == 0:
            msg = 'Unable to get XML for Tally ID={0} since it does not ' \
                  'contain any scores'.format(self.id)
            raise ValueError(msg)

        else:
            scores = ''
            for score in self.scores:
                scores += '{0} '.format(score)

            subelement = ET.SubElement(element,    "scores")
            subelement.text = scores.rstrip(' ')

        # Tally estimator type
        if self.estimator is not None:
            subelement = ET.SubElement(element, "estimator")
            subelement.text = self.estimator

        # Optional Triggers
        for trigger in self.triggers:
            trigger.get_trigger_xml(element)

        return element

    def find_filter(self, filter_type):
        """Return a filter in the tally that matches a specified type

        Parameters
        ----------
        filter_type : str
            Type of the filter, e.g. 'mesh'

        Returns
        -------
        filter : openmc.filter.Filter
            Filter from this tally with matching type, or None if no matching
            Filter is found

        Raises
        ------
        ValueError
            If no matching Filter is found

        """

        filter = None

        # Look through all of this Tally's Filters for the type requested
        for test_filter in self.filters:
            if test_filter.type == filter_type:
                filter = test_filter
                break

        # If we did not find the Filter, throw an Exception
        if filter is None:
            msg = 'Unable to find filter type "{0}" in ' \
                  'Tally ID={1}'.format(filter_type, self.id)
            raise ValueError(msg)

        return filter

    def get_filter_index(self, filter_type, filter_bin):
        """Returns the index in the Tally's results array for a Filter bin

        Parameters
        ----------
        filter_type : str
            The type of Filter (e.g., 'cell', 'energy', etc.)

        filter_bin : int, list
            The bin is an integer ID for 'material', 'surface', 'cell',
            'cellborn', and 'universe' Filters. The bin is an integer for the
            cell instance ID for 'distribcell' Filters. The bin is a 2-tuple of
            floats for 'energy' and 'energyout' filters corresponding to the
            energy boundaries of the bin of interest.  The bin is a (x,y,z)
            3-tuple for 'mesh' filters corresponding to the mesh cell of
            interest.

        Returns
        -------
             The index in the Tally data array for this filter bin.

        """

        # Find the equivalent Filter in this Tally's list of Filters
        filter = self.find_filter(filter_type)

        # Get the index for the requested bin from the Filter and return it
        filter_index = filter.get_bin_index(filter_bin)
        return filter_index

    def get_nuclide_index(self, nuclide):
        """Returns the index in the Tally's results array for a Nuclide bin

        Parameters
        ----------
        nuclide : str
            The name of the Nuclide (e.g., 'H-1', 'U-238')

        Returns
        -------
        nuclide_index : int
            The index in the Tally data array for this nuclide.

        Raises
        ------
        KeyError
            When the argument passed to the 'nuclide' parameter cannot be found
            in the Tally.

        """

        nuclide_index = -1

        # Look for the user-requested nuclide in all of the Tally's Nuclides
        for i, test_nuclide in enumerate(self.nuclides):

            # If the Summary was linked, then values are Nuclide objects
            if isinstance(test_nuclide, Nuclide):
                if test_nuclide._name == nuclide:
                    nuclide_index = i
                    break

            # If the Summary has not been linked, then values are ZAIDs
            else:
                if test_nuclide == nuclide:
                    nuclide_index = i
                    break

        if nuclide_index == -1:
            msg = 'Unable to get the nuclide index for Tally since "{0}" ' \
                  'is not one of the nuclides'.format(nuclide)
            raise KeyError(msg)
        else:
            return nuclide_index

    def get_score_index(self, score):
        """Returns the index in the Tally's results array for a score bin

        Parameters
        ----------
        score : str
            The score string (e.g., 'absorption', 'nu-fission')

        Returns
        -------
        score_index : int
            The index in the Tally data array for this score.

        Raises
        ------
        ValueError
            When the argument passed to the 'score' parameter cannot be found in
            the Tally.

        """

        try:
            score_index = self.scores.index(score)

        except ValueError:
            msg = 'Unable to get the score index for Tally since "{0}" ' \
                  'is not one of the scores'.format(score)
            raise ValueError(msg)

        return score_index

    def get_values(self, scores=[], filters=[], filter_bins=[],
                   nuclides=[], value='mean'):
        """Returns a tally score value given a list of filters to satisfy.

        This routine constructs a 3D NumPy array for the requested Tally data
        indexed by filter bin, nuclide bin, and score index. The routine will
        order the data in the array

        Parameters
        ----------
        scores : list
            A list of one or more score strings
            (e.g., ['absorption', 'nu-fission']; default is [])

        filters : list
            A list of filter type strings
            (e.g., ['mesh', 'energy']; default is [])

        filter_bins : list
            A list of the filter bins corresponding to the filter_types
            parameter (e.g., [1, (0., 0.625e-6)]; default is []). Each bin in
            the list is the integer ID for 'material', 'surface', 'cell',
            'cellborn', and 'universe' Filters. Each bin is an integer for the
            cell instance ID for 'distribcell Filters. Each bin is a 2-tuple of
            floats for 'energy' and 'energyout' filters corresponding to the
            energy boundaries of the bin of interest.  The bin is a (x,y,z)
            3-tuple for 'mesh' filters corresponding to the mesh cell of
            interest. The order of the bins in the list must correspond of the
            filter_types parameter.

        nuclides : list
            A list of nuclide name strings
            (e.g., ['U-235', 'U-238']; default is [])

        value : str
            A string for the type of value to return  - 'mean' (default),
            'std_dev', 'rel_err', 'sum', or 'sum_sq' are accepted

        Returns
        -------
        float or ndarray
            A scalar or NumPy array of the Tally data indexed in the order
            each filter, nuclide and score is listed in the parameters.

        Raises
        ------
        ValueError
            When this routine is called before the Tally is populated with data
            by the StatePoint.read_results() routine. ValueError is also thrown
            if the input parameters do not correspond to the Tally's attributes,
            e.g., if the score(s) do not match those in the Tally.

        """

        # Ensure that StatePoint.read_results() was called first
        if self._sum is None or self._sum_sq is None:
            msg = 'The Tally ID={0} has no data to return. Call the ' \
                  'StatePoint.read_results() routine before using ' \
                  'Tally.get_values(...)'.format(self.id)
            raise ValueError(msg)

        # Compute batch statistics if not yet computed
        self.compute_std_dev()

        ############################      FILTERS      #########################
        # Determine the score indices from any of the requested scores
        if filters:
            # Initialize empty list of indices for each bin in each Filter
            filter_indices = []

            # Loop over all of the Tally's Filters
            for i, filter in enumerate(self.filters):
                # Initialize empty list of indices for this Filter's bins
                filter_indices.append([])

                user_filter = False

                # If a user-requested Filter, get the user-requested bins
                for j, test_filter in enumerate(filters):
                    if filter.type == test_filter:
                        bins = filter_bins[j]
                        user_filter = True
                        break

                # If not a user-requested Filter, get all bins
                if not user_filter:
                    # Create list of 2- or 3-tuples tuples for mesh cell bins
                    if filter.type == 'mesh':
                        dimension = filter.mesh.dimension
                        xyz = map(lambda x: np.arange(1, x+1), dimension)
                        bins = list(itertools.product(*xyz))

                    # Create list of 2-tuples for energy boundary bins
                    elif filter.type in ['energy', 'energyout']:
                        bins = []
                        for i in range(filter.num_bins):
                            bins.append((filter.bins[i], filter.bins[i+1]))

                    # Create list of IDs for bins for all other Filter types
                    else:
                        bins = filter.bins

                # Add indices for each bin in this Filter to the list
                for bin in bins:
                    filter_indices[i].append(
                        self.get_filter_index(filter.type, bin))

            # Apply cross-product sum between all filter bin indices
            filter_indices = list(map(sum, itertools.product(*filter_indices)))

        # If user did not specify any specific Filters, use them all
        else:
            filter_indices = np.arange(self.num_filter_bins)

        ############################      NUCLIDES      ########################
        # Determine the score indices from any of the requested scores
        if nuclides:
            nuclide_indices = np.zeros(len(nuclides), dtype=np.int)
            for i, nuclide in enumerate(nuclides):
                nuclide_indices[i] = self.get_nuclide_index(nuclide)

        # If user did not specify any specific Nuclides, use them all
        else:
            nuclide_indices = np.arange(self.num_nuclides)

        #############################      SCORES      #########################
        # Determine the score indices from any of the requested scores
        if scores:
            score_indices = np.zeros(len(scores), dtype=np.int)
            for i, score in enumerate(scores):
                score_indices[i] = self.get_score_index(score)

        # If user did not specify any specific scores, use them all
        else:
            score_indices = np.arange(self.num_scores)

        # Construct cross-product of all three index types with each other
        indices = np.ix_(filter_indices, nuclide_indices, score_indices)

        # Return the desired result from Tally
        if value == 'mean':
            data = self.mean[indices]
        elif value == 'std_dev':
            data = self.std_dev[indices]
        elif value == 'rel_err':
            data = self.std_dev[indices] / self.mean[indices]
        elif value == 'sum':
            data = self.sum[indices]
        elif value == 'sum_sq':
            data = self.sum_sq[indices]
        else:
            msg = 'Unable to return results from Tally ID={0} since the ' \
                  'the requested value "{1}" is not \'mean\', \'std_dev\', ' \
                  '\rel_err\', \'sum\', or \'sum_sq\''.format(self.id, value)
            raise LookupError(msg)

        return data.squeeze()

    def get_pandas_dataframe(self, filters=True, nuclides=True,
                             scores=True, summary=None):
        """Build a Pandas DataFrame for the Tally data.

        This routine constructs a Pandas DataFrame object for the Tally data
        with columns annotated by filter, nuclide and score bin information.
        This capability has been tested for Pandas >=v0.13.1. However, if p
        possible, it is recommended to use the v0.16 or newer versions of
        Pandas since this this routine uses the Multi-index Pandas feature.

        Parameters
        ----------
        filters : bool
            Include columns with filter bin information (default is True).

        nuclides : bool
            Include columns with nuclide bin information (default is True).

        scores : bool
            Include columns with score bin information (default is True).

        summary : None or Summary
            An optional Summary object to be used to construct columns for
            distribcell tally filters (default is None). The geometric
            information in the Summary object is embedded into a Multi-index
            column with a geometric "path" to each distribcell intance.  NOTE:
            This option requires the OpenCG Python package.

        Returns
        -------
        pandas.DataFrame
            A Pandas DataFrame with each column annotated by filter, nuclide and
            score bin information (if these parameters are True), and the mean
            and standard deviation of the Tally's data.

        Raises
        ------
        KeyError
            When this routine is called before the Tally is populated with data
            by the StatePoint.read_results() routine.

        """

        # Ensure that StatePoint.read_results() was called first
        if self._sum is None or self._sum_sq is None:
            msg = 'The Tally ID={0} has no data to return. Call the ' \
                  'StatePoint.read_results() routine before using ' \
                  'Tally.get_pandas_dataframe(...)'.format(self.id)
            raise KeyError(msg)

        # If using Summary, ensure StatePoint.link_with_summary(...) was called
        if summary and not self.with_summary:
            msg = 'The Tally ID={0} has not been linked with the Summary. ' \
                  'Call the StatePoint.link_with_summary(...) routine ' \
                  'before using Tally.get_pandas_dataframe(...) with ' \
                  'Summary info'.format(self.id)
            raise KeyError(msg)

        # Attempt to import the pandas package
        try:
            import pandas as pd
        except ImportError:
            msg = 'The pandas Python package must be installed on your system'
            raise ImportError(msg)

        # Compute batch statistics if not yet computed
        self.compute_std_dev()

        # Initialize a pandas dataframe for the tally data
        df = pd.DataFrame()

        # Find the total length of the tally data array
        data_size = self.sum.size

        # Build DataFrame columns for filters if user requested them
        if filters:
            for filter in self.filters:
                # mesh filters
                if filter.type == 'mesh':

                    # Initialize dictionary to build Pandas Multi-index column
                    filter_dict = {}

                    # Append Mesh ID as outermost index of mult-index
                    mesh_id = filter.mesh.id
                    mesh_key = 'mesh {0}'.format(mesh_id)

                    # Find mesh dimensions - use 3D indices for simplicity
                    if (len(filter.mesh.dimension) == 3):
                        nx, ny, nz = filter.mesh.dimension
                    else:
                        nx, ny = filter.mesh.dimension
                        nz = 1

                    # Generate multi-index sub-column for x-axis
                    filter_bins = np.arange(1, nx+1)
                    repeat_factor = ny * nz * filter.stride
                    filter_bins = np.repeat(filter_bins, repeat_factor)
                    tile_factor = data_size / len(filter_bins)
                    filter_bins = np.tile(filter_bins, tile_factor)
                    filter_dict[(mesh_key, 'x')] = filter_bins

                    # Generate multi-index sub-column for y-axis
                    filter_bins = np.arange(1, ny+1)
                    repeat_factor = nz * filter.stride
                    filter_bins = np.repeat(filter_bins, repeat_factor)
                    tile_factor = data_size / len(filter_bins)
                    filter_bins = np.tile(filter_bins, tile_factor)
                    filter_dict[(mesh_key, 'y')] = filter_bins

                    # Generate multi-index sub-column for z-axis
                    filter_bins = np.arange(1, nz+1)
                    repeat_factor = filter.stride
                    filter_bins = np.repeat(filter_bins, repeat_factor)
                    tile_factor = data_size / len(filter_bins)
                    filter_bins = np.tile(filter_bins, tile_factor)
                    filter_dict[(mesh_key, 'z')] = filter_bins

                    # Append the multi-index column to the DataFrame
                    df = pd.concat([df, pd.DataFrame(filter_dict)], axis=1)

                # distribcell filters
                elif filter.type == 'distribcell':
                    if isinstance(summary, Summary):
                        # Attempt to import the OpenCG package
                        try:
                            import opencg
                        except ImportError:
                            msg = 'The OpenCG package must be installed ' \
                                  'to use a Summary for distribcell dataframes'
                            raise ImportError(msg)

                        # Create and extract the OpenCG geometry the Summary
                        summary.make_opencg_geometry()
                        opencg_geometry = summary.opencg_geometry
                        openmc_geometry = summary.openmc_geometry

                        # Use OpenCG to compute the number of regions
                        opencg_geometry.initializeCellOffsets()
                        num_regions = opencg_geometry._num_regions

                        # Initialize a dictionary mapping OpenMC distribcell
                        # offsets to OpenCG LocalCoords linked lists
                        offsets_to_coords = {}

                        # Use OpenCG to compute LocalCoords linked list for
                        # each region and store in dictionary
                        for region in range(num_regions):
                            coords = opencg_geometry.findRegion(region)
                            path = opencg.get_path(coords)
                            cell_id = path[-1]

                            # If this region is in Cell corresponding to the
                            # distribcell filter bin, store it in dictionary
                            if cell_id == filter.bins[0]:
                                offset = openmc_geometry.get_offset(path,
                                     filter.offset)
                                offsets_to_coords[offset] = coords

                        # Each distribcell offset is a DataFrame bin
                        # Unravel the paths into DataFrame columns
                        num_offsets = len(offsets_to_coords)

                        # Initialize termination condition for while loop
                        levels_remain = True
                        counter = 0

                        # Iterate over each level in the CSG tree hierarchy
                        while levels_remain:
                            levels_remain = False

                            # Initialize dictionary to build Pandas Multi-index
                            # column for this level in the CSG tree hierarchy
                            level_dict = {}

                            # Initialize prefix Multi-index keys
                            counter += 1
                            level_key = 'level {0}'.format(counter)
                            univ_key = (level_key, 'univ', 'id')
                            cell_key = (level_key, 'cell', 'id')
                            lat_id_key = (level_key, 'lat', 'id')
                            lat_x_key = (level_key, 'lat', 'x')
                            lat_y_key = (level_key, 'lat', 'y')
                            lat_z_key = (level_key, 'lat', 'z')

                            # Allocate NumPy arrays for each CSG level and
                            # each Multi-index column in the DataFrame
                            level_dict[univ_key] = np.empty(num_offsets)
                            level_dict[cell_key] = np.empty(num_offsets)
                            level_dict[lat_id_key] = np.empty(num_offsets)
                            level_dict[lat_x_key] = np.empty(num_offsets)
                            level_dict[lat_y_key] = np.empty(num_offsets)
                            level_dict[lat_z_key] = np.empty(num_offsets)

                            # Initialize Multi-index columns to NaN - this is
                            # necessary since some distribcell instances may
                            # have very different LocalCoords linked lists
                            level_dict[univ_key][:] = np.nan
                            level_dict[cell_key][:] = np.nan
                            level_dict[lat_id_key][:] = np.nan
                            level_dict[lat_x_key][:] = np.nan
                            level_dict[lat_y_key][:] = np.nan
                            level_dict[lat_z_key][:] = np.nan

                            # Iterate over all regions (distribcell instances)
                            for offset in range(num_offsets):
                                coords = offsets_to_coords[offset]

                                # If entire LocalCoords has been unraveled into
                                # Multi-index columns already, continue
                                if coords is None:
                                    continue

                                # Assign entry to Universe Multi-index column
                                if coords._type == 'universe':
                                    univ_id = coords._universe._id
                                    cell_id = coords._cell._id
                                    level_dict[univ_key][offset] = univ_id
                                    level_dict[cell_key][offset] = cell_id

                                # Assign entry to Lattice Multi-index column
                                else:
                                    lat_id = coords._lattice._id
                                    lat_x = coords._lat_x
                                    lat_y = coords._lat_y
                                    lat_z = coords._lat_z
                                    level_dict[lat_id_key][offset] = lat_id
                                    level_dict[lat_x_key][offset] = lat_x
                                    level_dict[lat_y_key][offset] = lat_y
                                    level_dict[lat_z_key][offset] = lat_z

                                # Move to next node in LocalCoords linked list
                                if coords._next is None:
                                    offsets_to_coords[offset] = None
                                else:
                                    offsets_to_coords[offset] = coords._next
                                    levels_remain = True

                            # Tile the Multi-index columns
                            for level_key, level_bins in level_dict.items():
                                level_bins = \
                                     np.repeat(level_bins, filter.stride)
                                tile_factor = data_size / len(level_bins)
                                level_bins = np.tile(level_bins, tile_factor)
                                level_dict[level_key] = level_bins

                            # Append the multi-index column to the DataFrame
                            df = pd.concat([df, pd.DataFrame(level_dict)],
                                           axis=1)

                    # Create DataFrame column for distribcell instances IDs
                    # NOTE: This is performed regardless of whether the user
                    # requests Summary geomeric information
                    filter_bins = np.arange(filter.num_bins)
                    filter_bins = np.repeat(filter_bins, filter.stride)
                    tile_factor = data_size / len(filter_bins)
                    filter_bins = np.tile(filter_bins, tile_factor)
                    df[filter.type] = filter_bins

                # energy, energyout filters
                elif 'energy' in filter.type:
                    bins = filter.bins
                    num_bins = filter.num_bins

                    # Create strings for
                    template = '{0:.1e} - {1:.1e}'
                    filter_bins = []
                    for i in range(num_bins):
                        filter_bins.append(template.format(bins[i], bins[i+1]))

                    # Tile the energy bins into a DataFrame column
                    filter_bins = np.repeat(filter_bins, filter.stride)
                    tile_factor = data_size / len(filter_bins)
                    filter_bins = np.tile(filter_bins, tile_factor)
                    df[filter.type + ' [MeV]'] = filter_bins

                # universe, material, surface, cell, and cellborn filters
                else:
                    filter_bins = np.repeat(filter.bins, filter.stride)
                    tile_factor = data_size / len(filter_bins)
                    filter_bins = np.tile(filter_bins, tile_factor)
                    df[filter.type] = filter_bins

        # Include DataFrame column for nuclides if user requested it
        if nuclides:
            nuclides = []

            for nuclide in self.nuclides:
                # Write Nuclide name if Summary info was linked with StatePoint
                if isinstance(nuclide, Nuclide):
                    nuclides.append(nuclide.name)
                else:
                    nuclides.append(nuclide)

            # Tile the nuclide bins into a DataFrame column
            nuclides = np.repeat(nuclides, len(self.scores))
            tile_factor = data_size / len(nuclides)
            df['nuclide'] = np.tile(nuclides, tile_factor)

        # Include column for scores if user requested it
        if scores:
            tile_factor = data_size / len(self.scores)
            df['score'] = np.tile(self.scores, tile_factor)

        # Append columns with mean, std. dev. for each tally bin
        df['mean'] = self.mean.ravel()
        df['std. dev.'] = self.std_dev.ravel()

        df.index.name = 'bin'
        df = df.dropna(axis=1)
        return df

    def export_results(self, filename='tally-results', directory='.',
                      format='hdf5', append=True):
        """Exports tallly results to an HDF5 or Python pickle binary file.

        Parameters
        ----------
        filename : str
            The name of the file for the results (default is 'tally-results')

        directory : str
            The name of the directory for the results (default is '.')

        format : str
            The format for the exported file - HDF5 ('hdf5', default) and
            Python pickle ('pkl') files are supported

        append : bool
            Whether or not to append the results to the file (default is True)

        Raises
        ------
        KeyError
            When this routine is called before the Tally is populated with data
            by the StatePoint.read_results() routine.

        """

        # Ensure that StatePoint.read_results() was called first
        if self._sum is None or self._sum_sq is None:
            msg = 'The Tally ID={0} has no data to export. Call the ' \
                  'StatePoint.read_results() routine before using ' \
                  'Tally.export_results(...)'.format(self.id)
            raise KeyError(msg)

        if not isinstance(filename, basestring):
            msg = 'Unable to export the results for Tally ID={0} to ' \
                  'filename="{1}" since it is not a ' \
                  'string'.format(self.id, filename)
            raise ValueError(msg)

        elif not isinstance(directory, basestring):
            msg = 'Unable to export the results for Tally ID={0} to ' \
                  'directory="{1}" since it is not a ' \
                  'string'.format(self.id, directory)
            raise ValueError(msg)

        elif format not in ['hdf5', 'pkl', 'csv']:
            msg = 'Unable to export the results for Tally ID={0} to format ' \
                  '"{1}" since it is not supported'.format(self.id, format)
            raise ValueError(msg)

        elif not isinstance(append, bool):
            msg = 'Unable to export the results for Tally ID={0} since the ' \
                  'append parameter is not True/False'.format(self.id, append)
            raise ValueError(msg)

        # Make directory if it does not exist
        if not os.path.exists(directory):
            os.makedirs(directory)

        # HDF5 binary file
        if format == 'hdf5':
            import h5py

            filename = directory + '/' + filename + '.h5'

            if append:
                tally_results = h5py.File(filename, 'a')
            else:
                tally_results = h5py.File(filename, 'w')

            # Create an HDF5 group within the file for this particular Tally
            tally_group = tally_results.create_group('Tally-{0}'.format(self.id))

            # Add basic Tally data to the HDF5 group
            tally_group.create_dataset('id', data=self.id)
            tally_group.create_dataset('name', data=self.name)
            tally_group.create_dataset('estimator', data=self.estimator)
            tally_group.create_dataset('scores', data=np.array(self.scores))

            # Add a string array of the nuclides to the HDF5 group
            nuclides = []

            for nuclide in self.nuclides:
                nuclides.append(nuclide.name)

            tally_group.create_dataset('nuclides', data=np.array(nuclides))

            # Create an HDF5 sub-group for the Filters
            filter_group = tally_group.create_group('filters')

            for filter in self.filters:
                filter_group.create_dataset(filter.type, data=filter.bins)

            # Add all results to the main HDF5 group for the Tally
            tally_group.create_dataset('sum', data=self.sum)
            tally_group.create_dataset('sum_sq', data=self.sum_sq)
            tally_group.create_dataset('mean', data=self.mean)
            tally_group.create_dataset('std_dev', data=self.std_dev)

            # Close the Tally results HDF5 file
            tally_results.close()

        # Python pickle binary file
        elif format == 'pkl':
            # Load the dictionary from the Pickle file
            filename = directory + '/' + filename + '.pkl'

            if os.path.exists(filename) and append:
                tally_results = pickle.load(file(filename, 'rb'))
            else:
                tally_results = {}

            # Create a nested dictionary within the file for this particular Tally
            tally_results['Tally-{0}'.format(self.id)] = {}
            tally_group = tally_results['Tally-{0}'.format(self.id)]

            # Add basic Tally data to the nested dictionary
            tally_group['id'] = self.id
            tally_group['name'] = self.name
            tally_group['estimator'] = self.estimator
            tally_group['scores'] = np.array(self.scores)

            # Add a string array of the nuclides to the HDF5 group
            nuclides = []

            for nuclide in self.nuclides:
                nuclides.append(nuclide.name)

            tally_group['nuclides'] = np.array(nuclides)

            # Create a nested dictionary for the Filters
            tally_group['filters'] = {}
            filter_group = tally_group['filters']

            for filter in self.filters:
                filter_group[filter.type] = filter.bins

            # Add all results to the main sub-dictionary for the Tally
            tally_group['sum'] = self.sum
            tally_group['sum_sq'] = self.sum_sq
            tally_group['mean'] = self.mean
            tally_group['std_dev'] = self.std_dev

            # Pickle the Tally results to a file
            pickle.dump(tally_results, open(filename, 'wb'))


class TalliesFile(object):
    """Tallies file used for an OpenMC simulation. Corresponds directly to the
    tallies.xml input file.

    """

    def __init__(self):
        # Initialize TalliesFile class attributes
        self._tallies = []
        self._meshes = []
        self._tallies_file = ET.Element("tallies")

    def add_tally(self, tally, merge=False):
        """Add a tally to the file

        Parameters
        ----------
        tally : Tally
            Tally to add to file
        merge : bool
            Indicate whether the tally should be merged with an existing tally,
            if possible. Defaults to False.

        """

        if not isinstance(tally, Tally):
            msg = 'Unable to add a non-Tally {0} to the TalliesFile'.format(tally)
            raise ValueError(msg)

        if merge:
            merged = False

            # Look for a tally to merge with this one
            for i, tally2 in enumerate(self._tallies):

                # If a mergeable tally is found
                if tally2.can_merge(tally):
                    # Replace tally 2 with the merged tally
                    merged_tally = tally2.merge(tally)
                    self._tallies[i] = merged_tally
                    merged = True
                    break

            # If not mergeable tally was found, simply add this tally
            if not merged:
                self._tallies.append(tally)

        else:
            self._tallies.append(tally)

    def remove_tally(self, tally):
        """Remove a tally from the file

        Parameters
        ----------
        tally : Tally
            Tally to remove

        """

        self._tallies.remove(tally)

    def merge_tallies(self):
        """Merge any mergeable tallies together. Note that n-way merges are
        possible.

        """

        for i, tally1 in enumerate(self._tallies):
            for j, tally2 in enumerate(self._tallies):
                # Do not merge the same tally with itself
                if i == j:
                    continue

                # If the two tallies are mergeable
                if tally1.can_merge(tally2):
                    # Replace tally 1 with the merged tally
                    merged_tally = tally1.merge(tally2)
                    self._tallies[i] = merged_tally

                    # Remove tally 2 since it is no longer needed
                    self._tallies.pop(j)

                    # Continue iterating from the first loop
                    break

    def add_mesh(self, mesh):
        """Add a mesh to the file

        Parameters
        ----------
        mesh : openmc.mesh.Mesh
            Mesh to add to the file

        """

        if not isinstance(mesh, Mesh):
            msg = 'Unable to add a non-Mesh {0} to the TalliesFile'.format(mesh)
            raise ValueError(msg)

        self._meshes.append(mesh)

    def remove_mesh(self, mesh):
        """Remove a mesh from the file

        Parameters
        ----------
        mesh : openmc.mesh.Mesh
            Mesh to remove from the file

        """

        self._meshes.remove(mesh)

    def _create_tally_subelements(self):
        for tally in self._tallies:
            xml_element = tally.get_tally_xml()
            self._tallies_file.append(xml_element)

    def _create_mesh_subelements(self):
        for mesh in self._meshes:
            if len(mesh._name) > 0:
                self._tallies_file.append(ET.Comment(mesh._name))

            xml_element = mesh.get_mesh_xml()
            self._tallies_file.append(xml_element)

    def export_to_xml(self):
        """Create a tallies.xml file that can be used for a simulation.

        """

        self._create_mesh_subelements()
        self._create_tally_subelements()

        # Clean the indentation in the file to be user-readable
        clean_xml_indentation(self._tallies_file)

        # Write the XML Tree to the tallies.xml file
        tree = ET.ElementTree(self._tallies_file)
        tree.write("tallies.xml", xml_declaration=True,
                             encoding='utf-8', method="xml")
