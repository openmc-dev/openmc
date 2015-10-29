from __future__ import division

from collections import Iterable, defaultdict
import copy
import os
import pickle
import itertools
from numbers import Integral, Real
from xml.etree import ElementTree as ET
import sys

import numpy as np

from openmc import Mesh, Filter, Trigger, Nuclide
from openmc.cross import CrossScore, CrossNuclide, CrossFilter
from openmc.filter import _FILTER_TYPES
import openmc.checkvalue as cv
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
    tally_id : Integral, optional
        Unique identifier for the tally. If none is specified, an identifier
        will automatically be assigned
    name : str, optional
        Name of the tally. If not specified, the name is the empty string.

    Attributes
    ----------
    id : Integral
        Unique identifier for the tally
    name : str
        Name of the tally
    filters : list of openmc.filter.Filter
        List of specified filters for the tally
    nuclides : list of openmc.nuclide.Nuclide
        List of nuclides to score results for
    scores : list of str
        List of defined scores, e.g. 'flux', 'fission', etc.
    estimator : {'analog', 'tracklength', 'collision'}
        Type of estimator for the tally
    triggers : list of openmc.trigger.Trigger
        List of tally triggers
    num_score_bins : Integral
        Total number of scores, accounting for the fact that a single
        user-specified score, e.g. scatter-P3 or flux-Y2,2, might have multiple
        bins
    num_scores : Integral
        Total number of user-specified scores
    num_filter_bins : Integral
        Total number of filter bins accounting for all filters
    num_bins : Integral
        Total number of bins for the tally
    num_realizations : Integral
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
        self._with_batch_statistics = False
        self._derived = False

        self._sp_filename = None
        self._results_read = False

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
            clone._sum = copy.deepcopy(self._sum, memo)
            clone._sum_sq = copy.deepcopy(self._sum_sq, memo)
            clone._mean = copy.deepcopy(self._mean, memo)
            clone._std_dev = copy.deepcopy(self._std_dev, memo)
            clone._with_summary = self.with_summary
            clone._with_batch_statistics = self.with_batch_statistics
            clone._derived = self.derived
            clone._sp_filename = self._sp_filename
            clone._results_read = self._results_read

            clone._filters = []
            for filter in self.filters:
                clone.add_filter(copy.deepcopy(filter, memo))

            clone._nuclides = []
            for nuclide in self.nuclides:
                clone.add_nuclide(copy.deepcopy(nuclide, memo))

            clone._scores = []
            for score in self.scores:
                clone.add_score(score)

            clone._triggers = []
            for trigger in self.triggers:
                clone.add_trigger(trigger)

            memo[id(self)] = clone

            return clone

        # If this object has been copied before, return the first copy made
        else:
            return existing

    def __eq__(self, other):
        if not isinstance(other, Tally):
            return False

        # Check all filters
        if len(self.filters) != len(other.filters):
            return False

        for filter in self.filters:
            if filter not in other.filters:
                return False

        # Check all nuclides
        if len(self.nuclides) != len(other.nuclides):
            return False

        for nuclide in self.nuclides:
            if nuclide not in other.nuclides:
                return False

        # Check all scores
        if len(self.scores) != len(other.scores):
            return False

        for score in self.scores:
            if score not in other.scores:
                return False

        if self.estimator != other.estimator:
            return False

        return True

    def __ne__(self, other):
        return not self == other

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

    def __repr__(self):
        string = 'Tally\n'
        string += '{0: <16}{1}{2}\n'.format('\tID', '=\t', self.id)
        string += '{0: <16}{1}{2}\n'.format('\tName', '=\t', self.name)

        string += '{0: <16}{1}\n'.format('\tFilters', '=\t')

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

        for filter in self.filters:
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
        if not self._sp_filename:
            return None

        if not self._results_read:
            import h5py

            # Open the HDF5 statepoint file
            f = h5py.File(self._sp_filename, 'r')

            # Extract Tally data from the file
            data = f['tallies/tally {0}/results'.format(
                self.id)].value
            sum = data['sum']
            sum_sq = data['sum_sq']

            # Define a routine to convert 0 to 1
            def nonzero(val):
                return 1 if not val else val

            # Reshape the results arrays
            new_shape = (nonzero(self.num_filter_bins),
                         nonzero(self.num_nuclides),
                         nonzero(self.num_score_bins))

            sum = np.reshape(sum, new_shape)
            sum_sq = np.reshape(sum_sq, new_shape)

            # Set the data for this Tally
            self._sum = sum
            self._sum_sq = sum_sq

            # Indicate that Tally results have been read
            self._results_read = True

            # Close the HDF5 statepoint file
            f.close()

        return self._sum

    @property
    def sum_sq(self):
        if not self._sp_filename:
            return None

        if not self._results_read:
            # Force reading of sum and sum_sq
            self.sum

        return self._sum_sq

    @property
    def mean(self):
        if self._mean is None:
            if not self._sp_filename:
                return None

            self._mean = self.sum / self.num_realizations
        return self._mean

    @property
    def std_dev(self):
        if self._std_dev is None:
            if not self._sp_filename:
                return None

            n = self.num_realizations
            nonzero = np.abs(self.mean) > 0
            self._std_dev = np.zeros_like(self.mean)
            self._std_dev[nonzero] = np.sqrt((self.sum_sq[nonzero]/n -
                                              self.mean[nonzero]**2)/(n - 1))
            self.with_batch_statistics = True
        return self._std_dev

    @property
    def with_batch_statistics(self):
        return self._with_batch_statistics

    @property
    def derived(self):
        return self._derived

    @estimator.setter
    def estimator(self, estimator):
        cv.check_value('estimator', estimator,
                    ['analog', 'tracklength', 'collision'])
        self._estimator = estimator

    def add_trigger(self, trigger):
        """Add a tally trigger to the tally

        Parameters
        ----------
        trigger : openmc.trigger.Trigger
            Trigger to add

        """

        if not isinstance(trigger, Trigger):
            msg = 'Unable to add a tally trigger for Tally ID="{0}" to ' \
                  'since "{1}" is not a Trigger'.format(self.id, trigger)
            raise ValueError(msg)

        if trigger not in self.triggers:
            self.triggers.append(trigger)

    @id.setter
    def id(self, tally_id):
        if tally_id is None:
            global AUTO_TALLY_ID
            self._id = AUTO_TALLY_ID
            AUTO_TALLY_ID += 1
        else:
            cv.check_type('tally ID', tally_id, Integral)
            cv.check_greater_than('tally ID', tally_id, 0, equality=True)
            self._id = tally_id

    @name.setter
    def name(self, name):
        if name is not None:
            cv.check_type('tally name', name, basestring)
            self._name = name
        else:
            self._name = None

    def add_filter(self, filter):
        """Add a filter to the tally

        Parameters
        ----------
        filter : openmc.filter.Filter
            Filter to add

        """

        if not isinstance(filter, (Filter, CrossFilter)):
            msg = 'Unable to add Filter "{0}" to Tally ID="{1}" since it is ' \
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

        if not isinstance(score, (basestring, CrossScore)):
            msg = 'Unable to add score "{0}" to Tally ID="{1}" since it is ' \
                  'not a string'.format(score, self.id)
            raise ValueError(msg)

        # If the score is already in the Tally, don't add it again
        if score in self.scores:
            return
        # Normal score strings
        if isinstance(score, basestring):
            self._scores.append(score.strip())
        # CrossScores
        else:
            self._scores.append(score)

    @num_score_bins.setter
    def num_score_bins(self, num_score_bins):
        self._num_score_bins = num_score_bins

    @num_realizations.setter
    def num_realizations(self, num_realizations):
        cv.check_type('number of realizations', num_realizations, Integral)
        cv.check_greater_than('number of realizations', num_realizations, 0, True)
        self._num_realizations = num_realizations

    @with_summary.setter
    def with_summary(self, with_summary):
        cv.check_type('with_summary', with_summary, bool)
        self._with_summary = with_summary

    @with_batch_statistics.setter
    def with_batch_statistics(self, with_batch_statistics):
        cv.check_type('with_batch_statistics', with_batch_statistics, bool)
        self._with_batch_statistics = with_batch_statistics

    @sum.setter
    def sum(self, sum):
        cv.check_type('sum', sum, Iterable)
        self._sum = sum

    @sum_sq.setter
    def sum_sq(self, sum_sq):
        cv.check_type('sum_sq', sum_sq, Iterable)
        self._sum_sq = sum_sq

    def remove_score(self, score):
        """Remove a score from the tally

        Parameters
        ----------
        score : str
            Score to remove

        """

        if score not in self.scores:
            msg = 'Unable to remove score "{0}" from Tally ID="{1}" since ' \
                  'the Tally does not contain this score'.format(score, self.id)
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
            msg = 'Unable to remove filter "{0}" from Tally ID="{1}" since the ' \
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
            msg = 'Unable to remove nuclide "{0}" from Tally ID="{1}" since the ' \
                  'Tally does not contain this nuclide'.format(nuclide, self.id)
            ValueError(msg)

        self._nuclides.remove(nuclide)

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

        # Check if only one tally contains a delayed group filter
        tally1_dg = False
        for filter1 in self.filters:
            if filter1.type == 'delayedgroup':
                tally1_dg = True

        tally2_dg = False
        for filter2 in tally.filters:
            if filter2.type == 'delayedgroup':
                tally2_dg = True

        # Return False if only one tally has a delayed group filter
        if (tally1_dg or tally2_dg) and not (tally1_dg and tally2_dg):
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
            msg = 'Unable to merge tally ID="{0}" with "{1}"'.format(tally.id, self.id)
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
            msg = 'Unable to get XML for Tally ID="{0}" since it does not ' \
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
                  'Tally ID="{1}"'.format(filter_type, self.id)
            raise ValueError(msg)

        return filter

    def get_filter_index(self, filter_type, filter_bin):
        """Returns the index in the Tally's results array for a Filter bin

        Parameters
        ----------
        filter_type : str
            The type of Filter (e.g., 'cell', 'energy', etc.)
        filter_bin : Integral or tuple
            The bin is an integer ID for 'material', 'surface', 'cell',
            'cellborn', and 'universe' Filters. The bin is an integer for the
            cell instance ID for 'distribcell' Filters. The bin is a 2-tuple of
            floats for 'energy' and 'energyout' filters corresponding to the
            energy boundaries of the bin of interest.  The bin is a (x,y,z)
            3-tuple for 'mesh' filters corresponding to the mesh cell of
            interest.

        Returns
        -------
             The index in the Tally data array for this filter bin

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

    def get_filter_indices(self, filters=[], filter_bins=[]):
        """Get indices into the filter axis of this tally's data arrays.

        This is a helper method for the Tally.get_values(...) method to
        extract tally data. This method returns the indices into the filter
        axis of the tally's data array (axis=0) for particular combinations
        of filters and their corresponding bins.

        Parameters
        ----------
        filters : list of str
            A list of filter type strings
            (e.g., ['mesh', 'energy']; default is [])
        filter_bins : list of Iterables
            A list of the filter bins corresponding to the filter_types
            parameter (e.g., [(1,), (0., 0.625e-6)]; default is []). Each bin
            in the list is the integer ID for 'material', 'surface', 'cell',
            'cellborn', and 'universe' Filters. Each bin is an integer for the
            cell instance ID for 'distribcell' Filters. Each bin is a 2-tuple of
            floats for 'energy' and 'energyout' filters corresponding to the
            energy boundaries of the bin of interest.  The bin is a (x,y,z)
            3-tuple for 'mesh' filters corresponding to the mesh cell of
            interest. The order of the bins in the list must correspond to the
            filter_types parameter.

        Returns
        -------
        ndarray
            A NumPy array of the filter indices

        """

        cv.check_iterable_type('filters', filters, basestring)
        cv.check_iterable_type('filter_bins', filter_bins, tuple)

        # Determine the score indices from any of the requested scores
        if filters:
            # Initialize empty list of indices for each bin in each Filter
            filter_indices = []

            # Loop over all of the Tally's Filters
            for i, filter in enumerate(self.filters):
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
                        for k in range(filter.num_bins):
                            bins.append((filter.bins[k], filter.bins[k+1]))

                    # Create list of cell instance IDs for distribcell Filters
                    elif filter.type == 'distribcell':
                        bins = np.arange(filter.num_bins)

                    # Create list of IDs for bins for all other filter types
                    else:
                        bins = filter.bins

                # Initialize a NumPy array for the Filter bin indices
                filter_indices.append(np.zeros(len(bins), dtype=np.int))

                # Add indices for each bin in this Filter to the list
                for j, bin in enumerate(bins):
                    filter_index = self.get_filter_index(filter.type, bin)
                    filter_indices[i][j] = filter_index

                # Account for stride in each of the previous filters
                for indices in filter_indices[:i]:
                    indices *= filter.num_bins

            # Apply outer product sum between all filter bin indices
            filter_indices = list(map(sum, itertools.product(*filter_indices)))

        # If user did not specify any specific Filters, use them all
        else:
            filter_indices = np.arange(self.num_filter_bins)

        return filter_indices

    def get_nuclide_indices(self, nuclides):
        """Get indices into the nuclide axis of this tally's data arrays.

        This is a helper method for the Tally.get_values(...) method to
        extract tally data. This method returns the indices into the nuclide
        axis of the tally's data array (axis=1) for one or more nuclides.

        Parameters
        ----------
        nuclides : list of str
            A list of nuclide name strings
            (e.g., ['U-235', 'U-238']; default is [])

        Returns
        -------
        ndarray
            A NumPy array of the nuclide indices

        """

        cv.check_iterable_type('nuclides', nuclides, basestring)

        # Determine the score indices from any of the requested scores
        if nuclides:
            nuclide_indices = np.zeros(len(nuclides), dtype=np.int)
            for i, nuclide in enumerate(nuclides):
                nuclide_indices[i] = self.get_nuclide_index(nuclide)

        # If user did not specify any specific Nuclides, use them all
        else:
            nuclide_indices = np.arange(self.num_nuclides)

        return nuclide_indices

    def get_score_indices(self, scores):
        """Get indices into the score axis of this tally's data arrays.

        This is a helper method for the Tally.get_values(...) method to
        extract tally data. This method returns the indices into the score
        axis of the tally's data array (axis=2) for one or more scores.

        Parameters
        ----------
        scores : list of str
            A list of one or more score strings
            (e.g., ['absorption', 'nu-fission']; default is [])

        Returns
        -------
        ndarray
            A NumPy array of the score indices

        """

        cv.check_iterable_type('scores', scores, basestring)

        # Determine the score indices from any of the requested scores
        if scores:
            score_indices = np.zeros(len(scores), dtype=np.int)
            for i, score in enumerate(scores):
                score_indices[i] = self.get_score_index(score)

        # If user did not specify any specific scores, use them all
        else:
            score_indices = np.arange(self.num_scores)

        return score_indices

    def get_values(self, scores=[], filters=[], filter_bins=[],
                   nuclides=[], value='mean'):
        """Returns one or more tallied values given a list of scores, filters,
        filter bins and nuclides.

        This method constructs a 3D NumPy array for the requested Tally data
        indexed by filter bin, nuclide bin, and score index. The method will
        order the data in the array as specified in the parameter lists.

        Parameters
        ----------
        scores : list of str
            A list of one or more score strings
            (e.g., ['absorption', 'nu-fission']; default is [])
        filters : list of str
            A list of filter type strings
            (e.g., ['mesh', 'energy']; default is [])
        filter_bins : list of Iterables
            A list of the filter bins corresponding to the filter_types
            parameter (e.g., [(1,), (0., 0.625e-6)]; default is []). Each bin
            in the list is the integer ID for 'material', 'surface', 'cell',
            'cellborn', and 'universe' Filters. Each bin is an integer for the
            cell instance ID for 'distribcell' Filters. Each bin is a 2-tuple of
            floats for 'energy' and 'energyout' filters corresponding to the
            energy boundaries of the bin of interest.  The bin is a (x,y,z)
            3-tuple for 'mesh' filters corresponding to the mesh cell of
            interest. The order of the bins in the list must correspond to the
            filter_types parameter.
        nuclides : list of str
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
            When this method is called before the Tally is populated with data
            by the StatePoint.read_results() method. ValueError is also thrown
            if the input parameters do not correspond to the Tally's attributes,
            e.g., if the score(s) do not match those in the Tally.

        """

        # Ensure that StatePoint.read_results() was called first
        if (value == 'mean' and self.mean is None) or \
           (value == 'std_dev' and self.std_dev is None) or \
           (value == 'rel_err' and self.mean is None) or \
           (value == 'sum' and self.sum is None) or \
           (value == 'sum_sq' and self.sum_sq is None):
            msg = 'The Tally ID="{0}" has no data to return. Call the ' \
                  'StatePoint.read_results() method before using ' \
                  'Tally.get_values(...)'.format(self.id)
            raise ValueError(msg)

        # Get filter, nuclide and score indices
        filter_indices = self.get_filter_indices(filters, filter_bins)
        nuclide_indices = self.get_nuclide_indices(nuclides)
        score_indices = self.get_score_indices(scores)

        # Construct outer product of all three index types with each other
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
            msg = 'Unable to return results from Tally ID="{0}" since the ' \
                  'the requested value "{1}" is not \'mean\', \'std_dev\', ' \
                  '\'rel_err\', \'sum\', or \'sum_sq\''.format(self.id, value)
            raise LookupError(msg)

        return data

    def get_pandas_dataframe(self, filters=True, nuclides=True,
                             scores=True, summary=None):
        """Build a Pandas DataFrame for the Tally data.

        This method constructs a Pandas DataFrame object for the Tally data
        with columns annotated by filter, nuclide and score bin information.

        This capability has been tested for Pandas >=0.13.1. However, it is
        recommended to use v0.16 or newer versions of Pandas since this method
        uses the Multi-index Pandas feature.

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
            column with a geometric "path" to each distribcell intance.
            NOTE: This option requires the OpenCG Python package.

        Returns
        -------
        pandas.DataFrame
            A Pandas DataFrame with each column annotated by filter, nuclide and
            score bin information (if these parameters are True), and the mean
            and standard deviation of the Tally's data.

        Raises
        ------
        KeyError
            When this method is called before the Tally is populated with data
            by the StatePoint.read_results() method.
        ImportError
            When Pandas can not be found on the caller's system

        """

        # Ensure that StatePoint.read_results() was called first
        if self.mean is None or self.std_dev is None:
            msg = 'The Tally ID="{0}" has no data to return. Call the ' \
                  'StatePoint.read_results() method before using ' \
                  'Tally.get_pandas_dataframe(...)'.format(self.id)
            raise KeyError(msg)

        # If using Summary, ensure StatePoint.link_with_summary(...) was called
        if summary and not self.with_summary:
            msg = 'The Tally ID="{0}" has not been linked with the Summary. ' \
                  'Call the StatePoint.link_with_summary(...) method ' \
                  'before using Tally.get_pandas_dataframe(...) with ' \
                  'Summary info'.format(self.id)
            raise KeyError(msg)

        # Attempt to import Pandas
        try:
            import pandas as pd
        except ImportError:
            msg = 'The Pandas Python package must be installed on your system'
            raise ImportError(msg)

        # Initialize a pandas dataframe for the tally data
        df = pd.DataFrame()

        # Find the total length of the tally data array
        data_size = self.mean.size

        # Build DataFrame columns for filters if user requested them
        if filters:

            # Append each Filter's DataFrame to the overall DataFrame
            for filter in self.filters:
                filter_df = filter.get_pandas_dataframe(data_size, summary)
                df = pd.concat([df, filter_df], axis=1)

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

        df = df.dropna(axis=1)

        # Expand the columns into Pandas MultiIndices for readability
        if pd.__version__ >= '0.16':
            columns = copy.deepcopy(df.columns.values)

            # Convert all elements in columns list to tuples
            for i, column in enumerate(columns):
                if not isinstance(column, tuple):
                    columns[i] = (column,)

            # Make each tuple the same length
            max_len_column = len(max(columns, key=len))
            for i, column in enumerate(columns):
                delta_len = max_len_column - len(column)
                if delta_len > 0:
                    new_column = list(column)
                    new_column.extend(['']*delta_len)
                    columns[i] = tuple(new_column)

            # Create and set a MultiIndex for the DataFrame's columns
            df.columns = pd.MultiIndex.from_tuples(columns)

        return df

    def get_reshaped_data(self, value='mean'):
        """Returns an array of tally data with one dimension per filter.

        The tally data in OpenMC is stored as a 3D array with the dimensions
        corresponding to filters, nuclides and scores. As a result, tally data
        can be opaque for a user to directly index (i.e., without use of the
        Tally.get_values(...) method) since one must know how to properly use
        the number of bins and strides for each filter to index into the first
        (filter) dimension.

        This builds and returns a reshaped version of the tally data array with
        unique dimensions corresponding to each tally filter. For example,
        suppose this tally has arrays of data with shape (8,5,5) corresponding
        to two filters (2 and 4 bins, respectively), five nuclides and five
        scores. This method will return a version of the data array with the
        with a new shape of (2,4,5,5) such that the first two dimensions
        correspond directly to the two filters with two and four bins.

        Parameters
        ----------
        value : str
            A string for the type of value to return  - 'mean' (default),
            'std_dev', 'rel_err', 'sum', or 'sum_sq' are accepted

        Returns
        -------
        ndarray
            The tally data array indexed by filters, nuclides and scores.

        """

        # Get the 3D array of data in filters, nuclides and scores
        data = self.get_values(value=value)

        # Build a new array shape with one dimension per filter
        new_shape = ()
        for filter in self.filters:
            new_shape += (filter.num_bins, )
        new_shape += (self.num_nuclides,)
        new_shape += (self.num_score_bins,)

        # Reshape the data with one dimension for each filter
        data = np.reshape(data, new_shape)
        return data

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
            When this method is called before the Tally is populated with data
            by the StatePoint.read_results() method.

        """

        # Ensure that StatePoint.read_results() was called first
        if self._sum is None or self._sum_sq is None and not self.derived:
            msg = 'The Tally ID="{0}" has no data to export. Call the ' \
                  'StatePoint.read_results() method before using ' \
                  'Tally.export_results(...)'.format(self.id)
            raise KeyError(msg)

        if not isinstance(filename, basestring):
            msg = 'Unable to export the results for Tally ID="{0}" to ' \
                  'filename="{1}" since it is not a ' \
                  'string'.format(self.id, filename)
            raise ValueError(msg)

        elif not isinstance(directory, basestring):
            msg = 'Unable to export the results for Tally ID="{0}" to ' \
                  'directory="{1}" since it is not a ' \
                  'string'.format(self.id, directory)
            raise ValueError(msg)

        elif format not in ['hdf5', 'pkl', 'csv']:
            msg = 'Unable to export the results for Tally ID="{0}" to format ' \
                  '"{1}" since it is not supported'.format(self.id, format)
            raise ValueError(msg)

        elif not isinstance(append, bool):
            msg = 'Unable to export the results for Tally ID="{0}" since the ' \
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

    def _outer_product(self, other, binary_op):
        """Combines filters, scores and nuclides with another tally.

        This is a helper method for the tally arithmetic methods. The filters,
        scores and nuclides from both tallies are enumerated into all possible
        combinations and expressed as CrossFilter, CrossScore and
        CrossNuclide objects in the new derived tally.

        Parameters
        ----------
        other : Tally
            The tally on the right hand side of the outer product
        binary_op : {'+', '-', '*', '/', '^'}
            The binary operation in the outer product

        Returns
        -------
        Tally
            A new Tally that is the outer product with this one.

        Raises
        ------
        ValueError
            When this method is called before the other tally is populated
            with data by the StatePoint.read_results() method.

        """

        # Check that results have been read
        if not other.derived and other.sum is None:
            msg = 'Unable to use tally arithmetic with Tally ID="{0}" ' \
                  'since it does not contain any results.'.format(other.id)
            raise ValueError(msg)

        new_tally = Tally()
        new_tally.with_batch_statistics = True
        new_tally._derived = True

        # Construct a combined derived name from the two tally operands
        if self.name != '' and other.name != '':
            new_name = '({0} {1} {2})'.format(self.name, binary_op, other.name)
            new_tally.name = new_name

        # Find any shared filters between the two tallies
        self_filters = set(self.filters)
        other_filters = set(other.filters)
        filter_intersect = self_filters.intersection(other_filters)

        # Align the shared filters to follow in each tally operand
        for i, filter in enumerate(filter_intersect):
            self_index = self.filters.index(filter)
            other_filter = other.filters[self_index]
            if other_filter != filter:
                other = other.swap_filters(filter, other_filter)

        data = self._align_tally_data(other)

        if binary_op == '+':
            new_tally._mean = data['self']['mean'] + data['other']['mean']
            new_tally._std_dev = np.sqrt(data['self']['std. dev.']**2 +
                                         data['other']['std. dev.']**2)
        elif binary_op == '-':
            new_tally._mean = data['self']['mean'] - data['other']['mean']
            new_tally._std_dev = np.sqrt(data['self']['std. dev.']**2 +
                                         data['other']['std. dev.']**2)
        elif binary_op == '*':
            self_rel_err = data['self']['std. dev.'] / data['self']['mean']
            other_rel_err = data['other']['std. dev.'] / data['other']['mean']
            new_tally._mean = data['self']['mean'] * data['other']['mean']
            new_tally._std_dev = np.abs(new_tally.mean) * \
                                 np.sqrt(self_rel_err**2 + other_rel_err**2)
        elif binary_op == '/':
            self_rel_err = data['self']['std. dev.'] / data['self']['mean']
            other_rel_err = data['other']['std. dev.'] / data['other']['mean']
            new_tally._mean = data['self']['mean'] / data['other']['mean']
            new_tally._std_dev = np.abs(new_tally.mean) * \
                                 np.sqrt(self_rel_err**2 + other_rel_err**2)
        elif binary_op == '^':
            mean_ratio = data['other']['mean'] / data['self']['mean']
            first_term = mean_ratio * data['self']['std. dev.']
            second_term = \
                np.log(data['self']['mean']) * data['other']['std. dev.']
            new_tally._mean = data['self']['mean'] ** data['other']['mean']
            new_tally._std_dev = np.abs(new_tally.mean) * \
                                 np.sqrt(first_term**2 + second_term**2)

        if self.estimator == other.estimator:
            new_tally.estimator = self.estimator
        if self.with_summary and other.with_summary:
            new_tally.with_summary = self.with_summary
        if self.num_realizations == other.num_realizations:
            new_tally.num_realizations = self.num_realizations

        # If filters are identical, simply reuse them in derived tally
        if self.filters == other.filters:
            for self_filter in self.filters:
                new_tally.add_filter(self_filter)

        # Generate filter "outer products" for non-identical filters
        else:

            # Find the common longest sequence of shared filters
            match = 0
            for self_filter, other_filter in zip(self.filters, other.filters):
                if self_filter == other_filter:
                    match += 1
                else:
                    break

            match_filters = self.filters[:match]
            cross_filters = [self.filters[match:], other.filters[match:]]

            # Simply reuse shared filters in derived tally
            for filter in match_filters:
                new_tally.add_filter(filter)

            # Use cross filters to combine non-shared filters in derived tally
            if len(self.filters) != match and len(other.filters) == match:
                for filter in cross_filters[0]:
                    new_tally.add_filter(filter)
            elif len(other.filters) == match and len(other.filters) != match:
                for filter in cross_filters[1]:
                    new_tally.add_filter(filter)
            else:
                for self_filter, other_filter in itertools.product(*cross_filters):
                    new_filter = CrossFilter(self_filter, other_filter, binary_op)
                    new_tally.add_filter(new_filter)

        # Generate score "outer products"
        if self.scores == other.scores:
            new_tally.num_score_bins = self.num_score_bins
            for self_score in self.scores:
                new_tally.add_score(self_score)
        else:
            new_tally.num_score_bins = self.num_score_bins * other.num_score_bins
            all_scores = [self.scores, other.scores]
            for self_score, other_score in itertools.product(*all_scores):
                new_score = CrossScore(self_score, other_score, binary_op)
                new_tally.add_score(new_score)

        # Generate nuclide "outer products"
        if self.nuclides == other.nuclides:
            for self_nuclide in self.nuclides:
                new_tally.nuclides.append(self_nuclide)
        else:
            all_nuclides = [self.nuclides, other.nuclides]
            for self_nuclide, other_nuclide in itertools.product(*all_nuclides):
                new_nuclide = CrossNuclide(self_nuclide, other_nuclide, binary_op)
                new_tally.add_nuclide(new_nuclide)

        # Correct each Filter's stride
        stride = new_tally.num_nuclides * new_tally.num_score_bins
        for filter in reversed(new_tally.filters):
            filter.stride = stride
            stride *= filter.num_bins

        return new_tally

    def _align_tally_data(self, other):
        """Aligns data from two tallies for tally arithmetic.

        This is a helper method to construct a dict of dicts of the "aligned"
        data arrays from each tally for tally arithmetic. The method analyzes
        the filters, scores and nuclides in both tally's and determines how to
        appropriately align the data for vectorized arithmetic. For example,
        if the two tallies have different filters, this method will use NumPy
        'tile' and 'repeat' operations to the new data arrays such that all
        possible combinations of the data in each tally's bins will be made
        when the arithmetic operation is applied to the arrays.

        Parameters
        ----------
        other : Tally
            The tally to outer product with this tally

        Returns
        -------
        dict
            A dictionary of dictionaries to "aligned" 'mean' and 'std. dev'
            NumPy arrays for each tally's data.

        """

        self_mean = copy.deepcopy(self.mean)
        self_std_dev = copy.deepcopy(self.std_dev)
        other_mean = copy.deepcopy(other.mean)
        other_std_dev = copy.deepcopy(other.std_dev)

        if self.filters != other.filters:

            # Determine the number of paired combinations of filter bins
            # between the two tallies and repeat arrays along filter axes
            diff1 = list(set(self.filters).difference(set(other.filters)))
            diff2 = list(set(other.filters).difference(set(self.filters)))

            # Determine the factors by which each tally operands' data arrays
            # must be tiled or repeated for the tally outer product
            other_tile_factor = 1
            self_repeat_factor = 1
            for filter in diff1:
                other_tile_factor *= filter.num_bins
            for filter in diff2:
                self_repeat_factor *= filter.num_bins

            # Tile / repeat the tally data for the tally outer product
            self_shape = list(self.mean.shape)
            other_shape = list(other.mean.shape)
            self_shape[0] *= self_repeat_factor
            self_mean = np.repeat(self_mean, self_repeat_factor)
            self_std_dev = np.repeat(self_std_dev, self_repeat_factor)

            if self_repeat_factor == 1:
                other_shape[0] *= other_tile_factor
                other_mean = np.repeat(other_mean, other_tile_factor, axis=0)
                other_std_dev = np.repeat(other_std_dev, other_tile_factor, axis=0)
            else:
                other_mean = np.tile(other_mean, (other_tile_factor, 1, 1))
                other_std_dev = np.tile(other_std_dev, (other_tile_factor, 1, 1))

            # NumPy repeat and tile routines return 1D flattened arrays
            # Reshape arrays as 3D with filters, nuclides and scores axes
            self_mean.shape = tuple(self_shape)
            self_std_dev.shape = tuple(self_shape)
            other_mean.shape = tuple(other_shape)
            other_std_dev.shape = tuple(other_shape)

        if self.nuclides != other.nuclides:

            # Determine the number of paired combinations of nuclides
            # between the two tallies and repeat arrays along nuclide axes
            self_repeat_factor = other.num_nuclides
            other_tile_factor = self.num_nuclides

            # Replicate the data
            self_mean = np.repeat(self_mean, self_repeat_factor, axis=1)
            other_mean = np.tile(other_mean, (1, other_tile_factor, 1))
            self_std_dev = np.repeat(self_std_dev, self_repeat_factor, axis=1)
            other_std_dev = np.tile(other_std_dev, (1, other_tile_factor, 1))

            # NumPy repeat and tile routines return 1D flattened arrays
            # Reshape arrays as 3D with filters, nuclides and scores axes
            self_shape = list(self.mean.shape)
            self_shape[1] *= self_repeat_factor
            self_mean.shape = tuple(self_shape)
            self_std_dev.shape = tuple(self_shape)

        if self.scores != other.scores:

            # Determine the number of paired combinations of score bins
            # between the two tallies and repeat arrays along score axes
            self_repeat_factor = other.num_score_bins
            other_tile_factor = self.num_score_bins

            # Replicate the data
            self_mean = np.repeat(self_mean, self_repeat_factor, axis=2)
            other_mean = np.tile(other_mean, (1, 1, other_tile_factor))
            self_std_dev = np.repeat(self_std_dev, self_repeat_factor, axis=2)
            other_std_dev = np.tile(other_std_dev, (1, 1, other_tile_factor))

            # NumPy repeat and tile routines return 1D flattened arrays
            # Reshape arrays as 3D with filters, nuclides and scores axes
            self_shape = list(self.mean.shape)
            self_shape[2] *= self_repeat_factor
            self_mean.shape = tuple(self_shape)
            self_std_dev.shape = tuple(self_shape)

        data = {}
        data['self'] = {}
        data['other'] = {}
        data['self']['mean'] = self_mean
        data['other']['mean'] = other_mean
        data['self']['std. dev.'] = self_std_dev
        data['other']['std. dev.'] = other_std_dev
        return data

    def swap_filters(self, filter1, filter2):
        """Reverse the ordering of two filters in this tally

        This is a helper method for tally arithmetic which helps align the data
        in two tallies with shared filters. This method copies this tally and
        reverses the order of the two filters.

        Parameters
        ----------
        filter1 : Filter
            The filter to swap with filter2

        filter2 : Filter
            The filter to swap with filter1

        Returns
        -------
        swap_tally
            A copy of this tally with the filters swapped

        Raises
        ------
        ValueError
            If this is a derived tally or this method is called before the tally
            is populated with data by the StatePoint.read_results() method.

        """

        # Check that results have been read
        if not self.derived and self.sum is None:
            msg = 'Unable to use tally arithmetic with Tally ID="{0}" ' \
                  'since it does not contain any results.'.format(self.id)
            raise ValueError(msg)

        cv.check_type('filter1', filter1, Filter)
        cv.check_type('filter2', filter2, Filter)

        if filter1 == filter2:
            msg = 'Unable to swap a filter with itself'
            raise ValueError(msg)
        elif filter1 not in self.filters:
            msg = 'Unable to swap "{0}" filter1 in Tally ID="{1}" since it ' \
                  'does not contain such a filter'.format(filter1.type, self.id)
            raise ValueError(msg)
        elif filter2 not in self.filters:
            msg = 'Unable to swap "{0}" filter2 in Tally ID="{1}" since it ' \
                  'does not contain such a filter'.format(filter2.type, self.id)
            raise ValueError(msg)

        swap_tally = copy.deepcopy(self)

        # Swap the filters in the copied version of this Tally
        filter1_index = swap_tally.filters.index(filter1)
        filter2_index = swap_tally.filters.index(filter2)
        swap_tally.filters[filter1_index] = filter2
        swap_tally.filters[filter2_index] = filter1

        # Update the strides for each of the filters
        stride = swap_tally.num_nuclides * swap_tally.num_score_bins
        for filter in reversed(swap_tally.filters):
            filter.stride = stride
            stride *= filter.num_bins

        # Construct lists of tuples for the bins in each of the two filters
        filters = [filter1.type, filter2.type]
        if filter1.type == 'distribcell':
            filter1_bins = np.arange(filter.num_bins)
        else:
            filter1_bins = [(filter1.get_bin(i)) for i in range(filter1.num_bins)]

        if filter1.type == 'distribcell':
            filter2_bins = np.arange(filter2.num_bins)
        else:
            filter2_bins = [filter2.get_bin(i) for i in range(filter2.num_bins)]

        # Adjust the sum data array to relect the new filter order
        if self.sum is not None:
            for bin1, bin2 in itertools.product(filter1_bins, filter2_bins):
                filter_bins = [(bin1,), (bin2,)]
                data = self.get_values(filters=filters,
                                       filter_bins=filter_bins, value='sum')
                indices = swap_tally.get_filter_indices(filters, filter_bins)
                swap_tally.sum[indices, :, :] = data

        # Adjust the sum_sq data array to relect the new filter order
        if self.sum_sq is not None:
            for bin1, bin2 in itertools.product(filter1_bins, filter2_bins):
                filter_bins = [(bin1,), (bin2,)]
                data = self.get_values(filters=filters,
                                       filter_bins=filter_bins, value='sum_sq')
                indices = swap_tally.get_filter_indices(filters, filter_bins)
                swap_tally.sum_sq[indices, :, :] = data

        # Adjust the mean data array to relect the new filter order
        if self.mean is not None:
            for bin1, bin2 in itertools.product(filter1_bins, filter2_bins):
                filter_bins = [(bin1,), (bin2,)]
                data = self.get_values(filters=filters,
                                       filter_bins=filter_bins, value='mean')
                indices = swap_tally.get_filter_indices(filters, filter_bins)
                swap_tally._mean[indices, :, :] = data

        # Adjust the std_dev data array to relect the new filter order
        if self.std_dev is not None:
            for bin1, bin2 in itertools.product(filter1_bins, filter2_bins):
                filter_bins = [(bin1,), (bin2,)]
                data = self.get_values(filters=filters,
                                       filter_bins=filter_bins, value='std_dev')
                indices = swap_tally.get_filter_indices(filters, filter_bins)
                swap_tally._std_dev[indices, :, :] = data

        return swap_tally

    def __add__(self, other):
        """Adds this tally to another tally or scalar value.

        This method builds a new tally with data that is the sum of this
        tally's data and that from the other tally or scalar value. If the
        filters, scores and nuclides in the two tallies are not the same, then
        they are combined in all possible ways in the new derived tally.

        Uncertainty propagation is used to compute the standard deviation
        for the new tally's data. It is important to note that this makes
        the assumption that the tally data is independently distributed.
        In most use cases, this is *not* true and may lead to under-prediction
        of the uncertainty. The uncertainty propagation model is from the
        following source:

        https://en.wikipedia.org/wiki/Propagation_of_uncertainty

        Parameters
        ----------
        other : Tally or Real
            The tally or scalar value to add to this tally

        Returns
        -------
        Tally
            A new derived tally which is the sum of this tally and the other
            tally or scalar value in the addition.

        Raises
        ------
        ValueError
            When this method is called before the Tally is populated with data
            by the StatePoint.read_results() method.

        """

        # Check that results have been read
        if not self.derived and self.sum is None:
            msg = 'Unable to use tally arithmetic with Tally ID="{0}" ' \
                  'since it does not contain any results.'.format(self.id)
            raise ValueError(msg)

        if isinstance(other, Tally):
            new_tally = self._outer_product(other, binary_op='+')

        elif isinstance(other, Real):
            new_tally = Tally(name='derived')
            new_tally._derived = True
            new_tally.with_batch_statistics = True
            new_tally.name = self.name
            new_tally._mean = self.mean + other
            new_tally._std_dev = self.std_dev
            new_tally.estimator = self.estimator
            new_tally.with_summary = self.with_summary
            new_tally.num_realization = self.num_realizations
            new_tally.num_score_bins = self.num_score_bins

            for filter in self.filters:
                new_tally.add_filter(filter)
            for nuclide in self.nuclides:
                new_tally.add_nuclide(nuclide)
            for score in self.scores:
                new_tally.add_score(score)

        else:
            msg = 'Unable to add "{0}" to Tally ID="{1}"'.format(other, self.id)
            raise ValueError(msg)

        return new_tally

    def __sub__(self, other):
        """Subtracts another tally or scalar value from this tally.

        This method builds a new tally with data that is the difference of
        this tally's data and that from the other tally or scalar value. If the
        filters, scores and nuclides in the two tallies are not the same, then
        they are combined in all possible ways in the new derived tally.

        Uncertainty propagation is used to compute the standard deviation
        for the new tally's data. It is important to note that this makes
        the assumption that the tally data is independently distributed.
        In most use cases, this is *not* true and may lead to under-prediction
        of the uncertainty. The uncertainty propagation model is from the
        following source:

        https://en.wikipedia.org/wiki/Propagation_of_uncertainty

        Parameters
        ----------
        other : Tally or Real
            The tally or scalar value to subtract from this tally

        Returns
        -------
        Tally
            A new derived tally which is the difference of this tally and the
            other tally or scalar value in the subtraction.

        Raises
        ------
        ValueError
            When this method is called before the Tally is populated with data
            by the StatePoint.read_results() method.

        """

        # Check that results have been read
        if not self.derived and self.sum is None:
            msg = 'Unable to use tally arithmetic with Tally ID="{0}" ' \
                  'since it does not contain any results.'.format(self.id)
            raise ValueError(msg)

        if isinstance(other, Tally):
            new_tally = self._outer_product(other, binary_op='-')

        elif isinstance(other, Real):
            new_tally = Tally(name='derived')
            new_tally._derived = True
            new_tally.name = self.name
            new_tally._mean = self.mean - other
            new_tally._std_dev = self.std_dev
            new_tally.estimator = self.estimator
            new_tally.with_summary = self.with_summary
            new_tally.num_realization = self.num_realizations
            new_tally.num_score_bins = self.num_score_bins

            for filter in self.filters:
                new_tally.add_filter(filter)
            for nuclide in self.nuclides:
                new_tally.add_nuclide(nuclide)
            for score in self.scores:
                new_tally.add_score(score)

        else:
            msg = 'Unable to subtract "{0}" from Tally ' \
                  'ID="{1}"'.format(other, self.id)
            raise ValueError(msg)

        return new_tally

    def __mul__(self, other):
        """Multiplies this tally with another tally or scalar value.

        This method builds a new tally with data that is the product of
        this tally's data and that from the other tally or scalar value. If the
        filters, scores and nuclides in the two tallies are not the same, then
        they are combined in all possible ways in the new derived tally.

        Uncertainty propagation is used to compute the standard deviation
        for the new tally's data. It is important to note that this makes
        the assumption that the tally data is independently distributed.
        In most use cases, this is *not* true and may lead to under-prediction
        of the uncertainty. The uncertainty propagation model is from the
        following source:

        https://en.wikipedia.org/wiki/Propagation_of_uncertainty

        Parameters
        ----------
        other : Tally or Real
            The tally or scalar value to multiply with this tally

        Returns
        -------
        Tally
            A new derived tally which is the product of this tally and the
            other tally or scalar value in the multiplication.

        Raises
        ------
        ValueError
            When this method is called before the Tally is populated with data
            by the StatePoint.read_results() method.

        """

        # Check that results have been read
        if not self.derived and self.sum is None:
            msg = 'Unable to use tally arithmetic with Tally ID="{0}" ' \
                  'since it does not contain any results.'.format(self.id)
            raise ValueError(msg)

        if isinstance(other, Tally):
            new_tally = self._outer_product(other, binary_op='*')

        elif isinstance(other, Real):
            new_tally = Tally(name='derived')
            new_tally._derived = True
            new_tally.name = self.name
            new_tally._mean = self.mean * other
            new_tally._std_dev = self.std_dev * np.abs(other)
            new_tally.estimator = self.estimator
            new_tally.with_summary = self.with_summary
            new_tally.num_realization = self.num_realizations
            new_tally.num_score_bins = self.num_score_bins

            for filter in self.filters:
                new_tally.add_filter(filter)
            for nuclide in self.nuclides:
                new_tally.add_nuclide(nuclide)
            for score in self.scores:
                new_tally.add_score(score)

        else:
            msg = 'Unable to multiply Tally ID="{0}" ' \
                  'by "{1}"'.format(self.id, other)
            raise ValueError(msg)

        return new_tally

    def __truediv__(self, other):
        """Divides this tally by another tally or scalar value.

        This method builds a new tally with data that is the dividend of
        this tally's data and that from the other tally or scalar value. If the
        filters, scores and nuclides in the two tallies are not the same, then
        they are combined in all possible ways in the new derived tally.

        Uncertainty propagation is used to compute the standard deviation
        for the new tally's data. It is important to note that this makes
        the assumption that the tally data is independently distributed.
        In most use cases, this is *not* true and may lead to under-prediction
        of the uncertainty. The uncertainty propagation model is from the
        following source:

        https://en.wikipedia.org/wiki/Propagation_of_uncertainty

        Parameters
        ----------
        other : Tally or Real
            The tally or scalar value to divide this tally by

        Returns
        -------
        Tally
            A new derived tally which is the dividend of this tally and the
            other tally or scalar value in the division.

        Raises
        ------
        ValueError
            When this method is called before the Tally is populated with data
            by the StatePoint.read_results() method.

        """

        # Check that results have been read
        if not self.derived and self.sum is None:
            msg = 'Unable to use tally arithmetic with Tally ID="{0}" ' \
                  'since it does not contain any results.'.format(self.id)
            raise ValueError(msg)

        if isinstance(other, Tally):
            new_tally = self._outer_product(other, binary_op='/')

        elif isinstance(other, Real):
            new_tally = Tally(name='derived')
            new_tally._derived = True
            new_tally.name = self.name
            new_tally._mean = self.mean / other
            new_tally._std_dev = self.std_dev * np.abs(1. / other)
            new_tally.estimator = self.estimator
            new_tally.with_summary = self.with_summary
            new_tally.num_realization = self.num_realizations
            new_tally.num_score_bins = self.num_score_bins

            for filter in self.filters:
                new_tally.add_filter(filter)
            for nuclide in self.nuclides:
                new_tally.add_nuclide(nuclide)
            for score in self.scores:
                new_tally.add_score(score)

        else:
            msg = 'Unable to divide Tally ID="{0}" ' \
                  'by "{1}"'.format(self.id, other)
            raise ValueError(msg)

        return new_tally

    def __div__(self, other):
        return self.__truediv__(other)

    def __pow__(self, power):
        """Raises this tally to another tally or scalar value power.

        This method builds a new tally with data that is the power of
        this tally's data to that from the other tally or scalar value. If the
        filters, scores and nuclides in the two tallies are not the same, then
        they are combined in all possible ways in the new derived tally.

        Uncertainty propagation is used to compute the standard deviation
        for the new tally's data. It is important to note that this makes
        the assumption that the tally data is independently distributed.
        In most use cases, this is *not* true and may lead to under-prediction
        of the uncertainty. The uncertainty propagation model is from the
        following source:

        https://en.wikipedia.org/wiki/Propagation_of_uncertainty

        Parameters
        ----------
        power : Tally or Real
            The tally or scalar value exponent

        Returns
        -------
        Tally
            A new derived tally which is this tally raised to the power of the
            other tally or scalar value in the exponentiation.

        Raises
        ------
        ValueError
            When this method is called before the Tally is populated with data
            by the StatePoint.read_results() method.

        """

        # Check that results have been read
        if not self.derived and self.sum is None:
            msg = 'Unable to use tally arithmetic with Tally ID="{0}" ' \
                  'since it does not contain any results.'.format(self.id)
            raise ValueError(msg)

        if isinstance(power, Tally):
            new_tally = self._outer_product(power, binary_op='^')

        elif isinstance(power, Real):
            new_tally = Tally(name='derived')
            new_tally._derived = True
            new_tally.name = self.name
            new_tally._mean = self._mean ** power
            self_rel_err = self.std_dev / self.mean
            new_tally._std_dev = np.abs(new_tally._mean * power * self_rel_err)
            new_tally.estimator = self.estimator
            new_tally.with_summary = self.with_summary
            new_tally.num_realization = self.num_realizations
            new_tally.num_score_bins = self.num_score_bins

            for filter in self.filters:
                new_tally.add_filter(filter)
            for nuclide in self.nuclides:
                new_tally.add_nuclide(nuclide)
            for score in self.scores:
                new_tally.add_score(score)

        else:
            msg = 'Unable to raise Tally ID="{0}" to ' \
                  'power "{1}"'.format(self.id, power)
            raise ValueError(msg)

        return new_tally

    def __radd__(self, other):
        """Right addition with a scalar value.

        This reverses the operands and calls the __add__ method.

        Parameters
        ----------
        other : Integer or Real
            The scalar value to add to this tally

        Returns
        -------
        Tally
            A new derived tally of this tally added with the scalar value.

        """

        return self + other

    def __rsub__(self, other):
        """Right subtraction from a scalar value.

        This reverses the operands and calls the __sub__ method.

        Parameters
        ----------
        other : Integer or Real
            The scalar value to subtract this tally from

        Returns
        -------
        Tally
            A new derived tally of this tally subtracted from the scalar value.

        """

        return -1. * self + other

    def __rmul__(self, other):
        """Right multiplication with a scalar value.

        This reverses the operands and calls the __mul__ method.

        Parameters
        ----------
        other : Integer or Real
            The scalar value to multiply with this tally

        Returns
        -------
        Tally
            A new derived tally of this tally multiplied by the scalar value.

        """

        return self * other

    def __rdiv__(self, other):
        """Right division with a scalar value.

        This reverses the operands and calls the __div__ method.

        Parameters
        ----------
        other : Integer or Real
            The scalar value to divide by this tally

        Returns
        -------
        Tally
            A new derived tally of the scalar value divided by this tally.

        """

        return other * self**-1

    def __pos__(self):
        """The absolute value of this tally.

        Returns
        -------
        Tally
            A new derived tally which is the absolute value of this tally.

        """

        new_tally = copy.deepcopy(self)
        new_tally._mean = np.abs(new_tally.mean)
        return new_tally

    def __neg__(self):
        """The negated value of this tally.

        Returns
        -------
        Tally
            A new derived tally which is the negated value of this tally.

        """

        new_tally = self * -1
        return new_tally

    def get_slice(self, scores=[], filters=[], filter_bins=[], nuclides=[]):
        """Build a sliced tally for the specified filters, scores and nuclides.

        This method constructs a new tally to encapsulate a subset of the data
        represented by this tally. The subset of data to include in the tally
        slice is determined by the scores, filters and nuclides specified in
        the input parameters.

        Parameters
        ----------
        scores : list of str
            A list of one or more score strings
            (e.g., ['absorption', 'nu-fission']; default is [])
        filters : list of str
            A list of filter type strings
            (e.g., ['mesh', 'energy']; default is [])
        filter_bins : list of Iterables
            A list of the filter bins corresponding to the filter_types
            parameter (e.g., [(1,), (0., 0.625e-6)]; default is []). Each bin
            in the list is the integer ID for 'material', 'surface', 'cell',
            'cellborn', and 'universe' Filters. Each bin is an integer for the
            cell instance ID for 'distribcell Filters. Each bin is a 2-tuple of
            floats for 'energy' and 'energyout' filters corresponding to the
            energy boundaries of the bin of interest.  The bin is a (x,y,z)
            3-tuple for 'mesh' filters corresponding to the mesh cell of
            interest. The order of the bins in the list must correspond to the
            filter_types parameter.
        nuclides : list of str
            A list of nuclide name strings
            (e.g., ['U-235', 'U-238']; default is [])

        Returns
        -------
        Tally
            A new tally which encapsulates the subset of data requested in the
            order each filter, nuclide and score is listed in the parameters.

        Raises
        ------
        ValueError
            When this method is called before the Tally is populated with data
            by the StatePoint.read_results() method.

        """

        # Ensure that StatePoint.read_results() was called first
        if not self.derived and self.sum is None:
            msg = 'Unable to use tally arithmetic with Tally ID="{0}" ' \
                  'since it does not contain any results.'.format(self.id)
            raise ValueError(msg)

        new_tally = copy.deepcopy(self)

        if self.sum is not None:
            new_sum = self.get_values(scores, filters, filter_bins,
                                      nuclides, 'sum')
            new_tally.sum = new_sum
        if self.sum_sq is not None:
            new_sum_sq = self.get_values(scores, filters, filter_bins,
                                         nuclides, 'sum_sq')
            new_tally.sum_sq = new_sum_sq
        if self.mean is not None:
            new_mean = self.get_values(scores, filters, filter_bins,
                                       nuclides, 'mean')
            new_tally._mean = new_mean
        if self.std_dev is not None:
            new_std_dev = self.get_values(scores, filters, filter_bins,
                                       nuclides, 'std_dev')
            new_tally._std_dev = new_std_dev

        # SCORES
        if scores:
            score_indices = []

            # Determine the score indices from any of the requested scores
            for score in self.scores:
                if score not in scores:
                    score_index = self.get_score_index(score)
                    score_indices.append(score_index)

            # Loop over indices in reverse to remove excluded scores
            for score_index in reversed(score_indices):
                new_tally.remove_score(self.scores[score_index])
                new_tally.num_score_bins -= 1

        # NUCLIDES
        if nuclides:
            nuclide_indices = []

            # Determine the nuclide indices from any of the requested nuclides
            for nuclide in self.nuclides:
                if nuclide.name not in nuclides:
                    nuclide_index = self.get_nuclide_index(nuclide.name)
                    nuclide_indices.append(nuclide_index)

            # Loop over indices in reverse to remove excluded Nuclides
            for nuclide_index in reversed(nuclide_indices):
                new_tally.remove_nuclide(self.nuclides[nuclide_index])

        # FILTERS
        if filters:

            # Determine the filter indices from any of the requested filters
            for i, filter_type in enumerate(filters):
                filter = new_tally.find_filter(filter_type)

                # Remove and/or reorder filter bins to user specifications
                bin_indices = []

                for filter_bin in filter_bins[i]:
                    bin_index = filter.get_bin_index(filter_bin)
                    if filter_type in ['energy', 'energyout']:
                        bin_indices.extend([bin_index, bin_index+1])
                    elif filter_type == 'distribcell':
                        indices = [(bin,) for bin in range(filter.num_bins)]
                        bin_indices.extend(indices)
                    else:
                        bin_indices.append(bin_index)

                filter.bins = filter.bins[bin_indices]
                filter.num_bins = len(filter_bins[i])

        # Correct each Filter's stride
        stride = new_tally.num_nuclides * new_tally.num_score_bins
        for filter in reversed(new_tally.filters):
            filter.stride = stride
            stride *= filter.num_bins

        return new_tally

    def summation(self, scores=[], filter_type=None,
                  filter_bins=[], nuclides=[]):
        """Vectorized sum of tally data across scores, filter bins and/or
        nuclides using tally addition.

        This method constructs a new tally to encapsulate the sum of the data
        represented by the summation of the data in this tally. The tally data
        sum is determined by the scores, filter bins and nuclides specified
        in the input parameters.

        Parameters
        ----------
        scores : list of str
            A list of one or more score strings to sum across
            (e.g., ['absorption', 'nu-fission']; default is [])
        filter_type : str
            A filter type string (e.g., 'cell', 'energy') corresponding to the
            filter bins to sum across
        filter_bins : Iterable of Integral or tuple
            A list of the filter bins corresponding to the filter_type parameter
            Each bin in the list is the integer ID for 'material', 'surface',
            'cell', 'cellborn', and 'universe' Filters. Each bin is an integer
            for the cell instance ID for 'distribcell Filters. Each bin is a
            2-tuple of floats for 'energy' and 'energyout' filters corresponding
            to the energy boundaries of the bin of interest. Each bin is an
            (x,y,z) 3-tuple for 'mesh' filters corresponding to the mesh cell of
            interest.
        nuclides : list of str
            A list of nuclide name strings to sum across
            (e.g., ['U-235', 'U-238']; default is [])

        Returns
        -------
        Tally
            A new tally which encapsulates the sum of data requested.
        """

        # If user did not specify any scores, do not sum across scores
        if len(scores) == 0:
            scores = [[]]
        # Sum across any scores specified by the user
        else:
            scores = [[score] for score in scores]

        # If user did not specify any nuclides, do not sum across nuclides
        if len(nuclides) == 0:
            nuclides = [[]]
        # Sum across any nuclides specified by the user
        else:
            nuclides = [[nuclide] for nuclide in nuclides]

        # Sum across any filter bins specified by the user
        if filter_type in _FILTER_TYPES:
            filter_bins = [[(filter_bin,)] for filter_bin in filter_bins]
            filters = [[filter_type]]
        # If user did not specify a filter type, do not sum across filter bins
        else:
            filter_bins = [[]]
            filters = [[]]

        # Initialize Tally sum
        tally_sum = 0

        # Iterate over all Tally slice operands in summation
        prod = [scores, filters, filter_bins, nuclides]
        summed_filters = defaultdict(list)
        for scores, filters, filter_bins, nuclides in itertools.product(*prod):
            tally_slice = self.get_slice(scores, filters, filter_bins, nuclides)

            # Remove filters summed across to avoid bulky CrossFilters
            if filter_type:
                filter = tally_slice.find_filter(filter_type)
                tally_slice.remove_filter(filter)
                summed_filters[filter_type].append(filter)

            # Accumulate this Tally slice into the Tally sum
            tally_sum += tally_slice

        # Add back the filter(s) which were summed across to derived tally
        for filter_type in summed_filters:
            filters = summed_filters[filter_type]
            for i in range(1, len(filters)):
                filters[i] = CrossFilter(filters[i-1], filters[i], '+')
            tally_sum.add_filter(filters[-1])

        return tally_sum

    def diagonalize_filter(self, new_filter):
        """Diagonalize the tally data array along a new axis of filter bins.

        This is a helper method for the tally arithmetic methods. This method
        adds the new filter to a derived tally constructed copied from this one.
        The data in the derived tally arrays is "diagonalized" along the bins in
        the new filter. This functionality is used by the openmc.mgxs module; to
        transport-correct scattering matrices by subtracting a 'scatter-P1'
        reaction rate tally with an energy filter from an 'scatter' reaction
        rate tally with both energy and energyout filters.

        Parameters
        ----------
        new_filter : Filter
            The filter along which to diagonalize the data in the new

        Returns
        -------
        Tally
            A new derived Tally with data diagaonalized along the new filter.

        """

        cv.check_type('new_filter', new_filter, Filter)

        if new_filter in self.filters:
            msg = 'Unable to diagonalize Tally ID="{0}" which already ' \
                  'contains a "{1}" filter'.format(self.id, new_filter.type)
            raise ValueError(msg)

        # Add the new filter to a copy of this Tally
        new_tally = copy.deepcopy(self)
        new_tally.add_filter(new_filter)

        # Determine the shape of data in the new diagonalized Tally
        num_filter_bins = new_tally.num_filter_bins
        num_nuclides = new_tally.num_nuclides
        num_score_bins = new_tally.num_score_bins
        new_shape = (num_filter_bins, num_nuclides, num_score_bins)

        # Determine "base" indices along the new "diagonal", and the factor
        # by which the "base" indices should be repeated to account for all
        # other filter bins in the diagonalized tally
        indices = np.arange(0, new_filter.num_bins**2, new_filter.num_bins+1)
        diag_factor = int(self.num_filter_bins / new_filter.num_bins)
        diag_indices = np.zeros(self.num_filter_bins, dtype=np.int)

        # Determine the filter indices along the new "diagonal"
        for i in range(diag_factor):
            start = i * new_filter.num_bins
            end = (i+1) * new_filter.num_bins
            diag_indices[start:end] = indices + (i * new_filter.num_bins**2)

        # Inject this Tally's data along the diagonal of the diagonalized Tally
        if self.sum is not None:
            new_tally._sum = np.zeros(new_shape, dtype=np.float64)
            new_tally._sum[diag_indices, :, :] = self.sum
        if self.sum_sq is not None:
            new_tally._sum_sq = np.zeros(new_shape, dtype=np.float64)
            new_tally._sum_sq[diag_indices, :, :] = self.sum_sq
        if self.mean is not None:
            new_tally._mean = np.zeros(new_shape, dtype=np.float64)
            new_tally._mean[diag_indices, :, :] = self.mean
        if self.std_dev is not None:
            new_tally._std_dev = np.zeros(new_shape, dtype=np.float64)
            new_tally._std_dev[diag_indices, :, :] = self.std_dev

        # Correct each Filter's stride
        stride = new_tally.num_nuclides * new_tally.num_score_bins
        for filter in reversed(new_tally.filters):
            filter.stride = stride
            stride *= filter.num_bins

        return new_tally


class TalliesFile(object):
    """Tallies file used for an OpenMC simulation. Corresponds directly to the
    tallies.xml input file.

    """

    def __init__(self):
        # Initialize TalliesFile class attributes
        self._tallies = []
        self._meshes = []
        self._tallies_file = ET.Element("tallies")

    @property
    def tallies(self):
        return self._tallies

    @property
    def meshes(self):
        return self._meshes

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
            msg = 'Unable to add a non-Tally "{0}" to the TalliesFile'.format(tally)
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
            msg = 'Unable to add a non-Mesh "{0}" to the TalliesFile'.format(mesh)
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

        # Reset xml element tree
        self._tallies_file.clear()

        self._create_mesh_subelements()
        self._create_tally_subelements()

        # Clean the indentation in the file to be user-readable
        clean_xml_indentation(self._tallies_file)

        # Write the XML Tree to the tallies.xml file
        tree = ET.ElementTree(self._tallies_file)
        tree.write("tallies.xml", xml_declaration=True,
                             encoding='utf-8', method="xml")
