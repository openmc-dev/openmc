from collections.abc import Iterable, MutableSequence
import copy
import re
from functools import partial, reduce
from itertools import product
from numbers import Integral, Real
import operator
from pathlib import Path
import warnings
from xml.etree import ElementTree as ET

import numpy as np
import pandas as pd
import scipy.sparse as sps
import h5py

import openmc
import openmc.checkvalue as cv
from openmc._xml import clean_indentation
from .mixin import IDManagerMixin


# The tally arithmetic product types. The tensor product performs the full
# cross product of the data in two tallies with respect to a specified axis
# (filters, nuclides, or scores). The entrywise product performs the arithmetic
# operation entrywise across the entries in two tallies with respect to a
# specified axis.
_PRODUCT_TYPES = ['tensor', 'entrywise']

# The following indicate acceptable types when setting Tally.scores,
# Tally.nuclides, and Tally.filters
_SCORE_CLASSES = (str, openmc.CrossScore, openmc.AggregateScore)
_NUCLIDE_CLASSES = (str, openmc.CrossNuclide, openmc.AggregateNuclide)
_FILTER_CLASSES = (openmc.Filter, openmc.CrossFilter, openmc.AggregateFilter)

# Valid types of estimators
ESTIMATOR_TYPES = ['tracklength', 'collision', 'analog']


class Tally(IDManagerMixin):
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
    filters : list of openmc.Filter
        List of specified filters for the tally
    nuclides : list of openmc.Nuclide
        List of nuclides to score results for
    scores : list of str
        List of defined scores, e.g. 'flux', 'fission', etc.
    estimator : {'analog', 'tracklength', 'collision'}
        Type of estimator for the tally
    triggers : list of openmc.Trigger
        List of tally triggers
    num_scores : int
        Total number of scores
    num_filter_bins : int
        Total number of filter bins accounting for all filters
    num_bins : int
        Total number of bins for the tally
    shape : 3-tuple of int
        The shape of the tally data array ordered as the number of filter bins,
        nuclide bins and score bins
    filter_strides : list of int
        Stride in memory for each filter
    num_realizations : int
        Total number of realizations
    with_summary : bool
        Whether or not a Summary has been linked
    sum : numpy.ndarray
        An array containing the sum of each independent realization for each bin
    sum_sq : numpy.ndarray
        An array containing the sum of each independent realization squared for
        each bin
    mean : numpy.ndarray
        An array containing the sample mean for each bin
    std_dev : numpy.ndarray
        An array containing the sample standard deviation for each bin
    derived : bool
        Whether or not the tally is derived from one or more other tallies
    sparse : bool
        Whether or not the tally uses SciPy's LIL sparse matrix format for
        compressed data storage
    derivative : openmc.TallyDerivative
        A material perturbation derivative to apply to all scores in the tally.

    """

    next_id = 1
    used_ids = set()

    def __init__(self, tally_id=None, name=''):
        # Initialize Tally class attributes
        self.id = tally_id
        self.name = name
        self._filters = cv.CheckedList(_FILTER_CLASSES, 'tally filters')
        self._nuclides = cv.CheckedList(_NUCLIDE_CLASSES, 'tally nuclides')
        self._scores = cv.CheckedList(_SCORE_CLASSES, 'tally scores')
        self._estimator = None
        self._triggers = cv.CheckedList(openmc.Trigger, 'tally triggers')
        self._derivative = None

        self._num_realizations = 0
        self._with_summary = False

        self._sum = None
        self._sum_sq = None
        self._mean = None
        self._std_dev = None
        self._with_batch_statistics = False
        self._derived = False
        self._sparse = False

        self._sp_filename = None
        self._results_read = False

    def __repr__(self):
        parts = ['Tally']
        parts.append('{: <15}=\t{}'.format('ID', self.id))
        parts.append('{: <15}=\t{}'.format('Name', self.name))
        if self.derivative is not None:
            parts.append('{: <15}=\t{}'.format('Derivative ID', self.derivative.id))
        filters = ', '.join(type(f).__name__ for f in self.filters)
        parts.append('{: <15}=\t{}'.format('Filters', filters))
        nuclides = ' '.join(str(nuclide) for nuclide in self.nuclides)
        parts.append('{: <15}=\t{}'.format('Nuclides', nuclides))
        parts.append('{: <15}=\t{}'.format('Scores', self.scores))
        parts.append('{: <15}=\t{}'.format('Estimator', self.estimator))
        return '\n\t'.join(parts)

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
    def num_filters(self):
        return len(self.filters)

    @property
    def num_filter_bins(self):
        return reduce(operator.mul, (f.num_bins for f in self.filters), 1)

    @property
    def num_bins(self):
        return self.num_filter_bins * self.num_nuclides * self.num_scores

    @property
    def shape(self):
        return (self.num_filter_bins, self.num_nuclides, self.num_scores)

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

    def _read_results(self):
        if self._results_read:
            return

        # Open the HDF5 statepoint file
        with h5py.File(self._sp_filename, 'r') as f:
            # Extract Tally data from the file
            data = f['tallies/tally {}/results'.format(self.id)]
            sum_ = data[:, :, 0]
            sum_sq = data[:, :, 1]

            # Reshape the results arrays
            sum_ = np.reshape(sum_, self.shape)
            sum_sq = np.reshape(sum_sq, self.shape)

            # Set the data for this Tally
            self._sum = sum_
            self._sum_sq = sum_sq

            # Convert NumPy arrays to SciPy sparse LIL matrices
            if self.sparse:
                self._sum = sps.lil_matrix(self._sum.flatten(), self._sum.shape)
                self._sum_sq = sps.lil_matrix(self._sum_sq.flatten(), self._sum_sq.shape)

        # Indicate that Tally results have been read
        self._results_read = True

    @property
    def sum(self):
        if not self._sp_filename or self.derived:
            return None

        # Make sure results have been read
        self._read_results()

        if self.sparse:
            return np.reshape(self._sum.toarray(), self.shape)
        else:
            return self._sum

    @property
    def sum_sq(self):
        if not self._sp_filename or self.derived:
            return None

        # Make sure results have been read
        self._read_results()

        if self.sparse:
            return np.reshape(self._sum_sq.toarray(), self.shape)
        else:
            return self._sum_sq

    @property
    def mean(self):
        if self._mean is None:
            if not self._sp_filename:
                return None

            self._mean = self.sum / self.num_realizations

            # Convert NumPy array to SciPy sparse LIL matrix
            if self.sparse:
                self._mean = sps.lil_matrix(self._mean.flatten(),
                                            self._mean.shape)

        if self.sparse:
            return np.reshape(self._mean.toarray(), self.shape)
        else:
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

            # Convert NumPy array to SciPy sparse LIL matrix
            if self.sparse:
                self._std_dev = sps.lil_matrix(self._std_dev.flatten(),
                                               self._std_dev.shape)

            self.with_batch_statistics = True

        if self.sparse:
            return np.reshape(self._std_dev.toarray(), self.shape)
        else:
            return self._std_dev

    @property
    def with_batch_statistics(self):
        return self._with_batch_statistics

    @property
    def derived(self):
        return self._derived

    @property
    def derivative(self):
        return self._derivative

    @property
    def sparse(self):
        return self._sparse

    @estimator.setter
    def estimator(self, estimator):
        cv.check_value('estimator', estimator, ESTIMATOR_TYPES)
        self._estimator = estimator

    @triggers.setter
    def triggers(self, triggers):
        cv.check_type('tally triggers', triggers, MutableSequence)
        self._triggers = cv.CheckedList(openmc.Trigger, 'tally triggers',
                                        triggers)

    @name.setter
    def name(self, name):
        cv.check_type('tally name', name, str, none_ok=True)
        self._name = name

    @derivative.setter
    def derivative(self, deriv):
        cv.check_type('tally derivative', deriv, openmc.TallyDerivative,
                      none_ok=True)
        self._derivative = deriv

    @filters.setter
    def filters(self, filters):
        cv.check_type('tally filters', filters, MutableSequence)

        # If the filter is already in the Tally, raise an error
        visited_filters = set()
        for f in filters:
            if f in visited_filters:
                msg = 'Unable to add a duplicate filter "{}" to Tally ID="{}" ' \
                      'since duplicate filters are not supported in the OpenMC ' \
                      'Python API'.format(f, self.id)
                raise ValueError(msg)
            visited_filters.add(f)

        self._filters = cv.CheckedList(_FILTER_CLASSES, 'tally filters', filters)

    @nuclides.setter
    def nuclides(self, nuclides):
        cv.check_type('tally nuclides', nuclides, MutableSequence)

        # If the nuclide is already in the Tally, raise an error
        visited_nuclides = set()
        for nuc in nuclides:
            if nuc in visited_nuclides:
                msg = 'Unable to add a duplicate nuclide "{}" to Tally ID="{}" ' \
                      'since duplicate nuclides are not supported in the OpenMC ' \
                      'Python API'.format(nuclide, self.id)
                raise ValueError(msg)
            visited_nuclides.add(nuc)

        self._nuclides = cv.CheckedList(_NUCLIDE_CLASSES, 'tally nuclides',
                                        nuclides)

    @scores.setter
    def scores(self, scores):
        cv.check_type('tally scores', scores, MutableSequence)

        visited_scores = set()
        for i, score in enumerate(scores):
            # If the score is already in the Tally, raise an error
            if score in visited_scores:
                msg = 'Unable to add a duplicate score "{}" to Tally ID="{}" ' \
                      'since duplicate scores are not supported in the OpenMC ' \
                      'Python API'.format(score, self.id)
                raise ValueError(msg)
            visited_scores.add(score)

            # If score is a string, strip whitespace
            if isinstance(score, str):
                # Check to see if scores are deprecated before storing
                for deprecated in ['scatter-', 'nu-scatter-', 'scatter-p',
                                   'nu-scatter-p', 'scatter-y', 'nu-scatter-y',
                                    'flux-y', 'total-y']:
                    if score.strip().startswith(deprecated):
                        msg = score.strip() + ' is no longer supported.'
                        raise ValueError(msg)
                scores[i] = score.strip()

        self._scores = cv.CheckedList(_SCORE_CLASSES, 'tally scores', scores)

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

        # Convert NumPy arrays to SciPy sparse LIL matrices
        if sparse and not self.sparse:
            if self._sum is not None:
                self._sum = sps.lil_matrix(self._sum.flatten(), self._sum.shape)
            if self._sum_sq is not None:
                self._sum_sq = sps.lil_matrix(self._sum_sq.flatten(),
                                              self._sum_sq.shape)
            if self._mean is not None:
                self._mean = sps.lil_matrix(self._mean.flatten(),
                                            self._mean.shape)
            if self._std_dev is not None:
                self._std_dev = sps.lil_matrix(self._std_dev.flatten(),
                                               self._std_dev.shape)

            self._sparse = True

        # Convert SciPy sparse LIL matrices to NumPy arrays
        elif not sparse and self.sparse:
            if self._sum is not None:
                self._sum = np.reshape(self._sum.toarray(), self.shape)
            if self._sum_sq is not None:
                self._sum_sq = np.reshape(self._sum_sq.toarray(), self.shape)
            if self._mean is not None:
                self._mean = np.reshape(self._mean.toarray(), self.shape)
            if self._std_dev is not None:
                self._std_dev = np.reshape(self._std_dev.toarray(), self.shape)
            self._sparse = False

    def remove_score(self, score):
        """Remove a score from the tally

        Parameters
        ----------
        score : str
            Score to remove

        """

        if score not in self.scores:
            msg = 'Unable to remove score "{}" from Tally ID="{}" since ' \
                  'the Tally does not contain this score'.format(score, self.id)
            raise ValueError(msg)

        self._scores.remove(score)

    def remove_filter(self, old_filter):
        """Remove a filter from the tally

        Parameters
        ----------
        old_filter : openmc.Filter
            Filter to remove

        """

        if old_filter not in self.filters:
            msg = 'Unable to remove filter "{}" from Tally ID="{}" since the ' \
                  'Tally does not contain this filter'.format(old_filter, self.id)
            raise ValueError(msg)

        self._filters.remove(old_filter)

    def remove_nuclide(self, nuclide):
        """Remove a nuclide from the tally

        Parameters
        ----------
        nuclide : openmc.Nuclide
            Nuclide to remove

        """

        if nuclide not in self.nuclides:
            msg = 'Unable to remove nuclide "{}" from Tally ID="{}" since the ' \
                  'Tally does not contain this nuclide'.format(nuclide, self.id)
            raise ValueError(msg)

        self._nuclides.remove(nuclide)

    def _can_merge_filters(self, other):
        """Determine if another tally's filters can be merged with this one's

        The types of filters between the two tallies must match identically.
        The bins in all of the filters must match identically, or be mergeable
        in only one filter. This is a helper method for the can_merge(...)
        and merge(...) methods.

        Parameters
        ----------
        other : openmc.Tally
            Tally to check for mergeable filters

        """

        # Two tallies must have the same number of filters
        if len(self.filters) != len(other.filters):
            return False

        # Return False if only one tally has a delayed group filter
        tally1_dg = self.contains_filter(openmc.DelayedGroupFilter)
        tally2_dg = other.contains_filter(openmc.DelayedGroupFilter)
        if tally1_dg != tally2_dg:
            return False

        # Look to see if all filters are the same, or one or more can be merged
        for filter1 in self.filters:
            mergeable = False

            for filter2 in other.filters:
                if filter1 == filter2 or filter1.can_merge(filter2):
                    mergeable = True
                    break

            # If no mergeable filter was found, the tallies are not mergeable
            if not mergeable:
                return False

        # Tally filters are mergeable if all conditional checks passed
        return True

    def _can_merge_nuclides(self, other):
        """Determine if another tally's nuclides can be merged with this one's

        The nuclides between the two tallies must be mutually exclusive or
        identically matching. This is a helper method for the can_merge(...)
        and merge(...) methods.

        Parameters
        ----------
        other : openmc.Tally
            Tally to check for mergeable nuclides

        """

        no_nuclides_match = True
        all_nuclides_match = True

        # Search for each of this tally's nuclides in the other tally
        for nuclide in self.nuclides:
            if nuclide not in other.nuclides:
                all_nuclides_match = False
            else:
                no_nuclides_match = False

        # Search for each of the other tally's nuclides in this tally
        for nuclide in other.nuclides:
            if nuclide not in self.nuclides:
                all_nuclides_match = False
            else:
                no_nuclides_match = False

        # Either all nuclides should match, or none should
        return no_nuclides_match or all_nuclides_match

    def _can_merge_scores(self, other):
        """Determine if another tally's scores can be merged with this one's

        The scores between the two tallies must be mutually exclusive or
        identically matching. This is a helper method for the can_merge(...)
        and merge(...) methods.

        Parameters
        ----------
        other : openmc.Tally
            Tally to check for mergeable scores

        """

        no_scores_match = True
        all_scores_match = True

        # Search for each of this tally's scores in the other tally
        for score in self.scores:
            if score in other.scores:
                no_scores_match = False

        # Search for each of the other tally's scores in this tally
        for score in other.scores:
            if score not in self.scores:
                all_scores_match = False
            else:
                no_scores_match = False

            if score == 'current' and score not in self.scores:
                return False

        # Nuclides cannot be specified on 'flux' scores
        if 'flux' in self.scores or 'flux' in other.scores:
            if self.nuclides != other.nuclides:
                return False

        # Either all scores should match, or none should
        return no_scores_match or all_scores_match

    def can_merge(self, other):
        """Determine if another tally can be merged with this one

        If results have been loaded from a statepoint, then tallies are only
        mergeable along one and only one of filter bins, nuclides or scores.

        Parameters
        ----------
        other : openmc.Tally
            Tally to check for merging

        """

        if not isinstance(other, Tally):
            return False

        # Must have same estimator
        if self.estimator != other.estimator:
            return False

        equal_filters = sorted(self.filters) == sorted(other.filters)
        equal_nuclides = sorted(self.nuclides) == sorted(other.nuclides)
        equal_scores = sorted(self.scores) == sorted(other.scores)
        equality = [equal_filters, equal_nuclides, equal_scores]

        # If all filters, nuclides and scores match then tallies are mergeable
        if all(equality):
            return True

        # Variables to indicate filter bins, nuclides, and scores that can be merged
        can_merge_filters = self._can_merge_filters(other)
        can_merge_nuclides = self._can_merge_nuclides(other)
        can_merge_scores = self._can_merge_scores(other)
        mergeability = [can_merge_filters, can_merge_nuclides, can_merge_scores]

        if not all(mergeability):
            return False

        # If the tally results have been read from the statepoint, at least two
        # of filters, nuclides and scores must match
        else:
            return not self._results_read or sum(equality) >= 2

    def merge(self, other):
        """Merge another tally with this one

        If results have been loaded from a statepoint, then tallies are only
        mergeable along one and only one of filter bins, nuclides or scores.

        Parameters
        ----------
        other : openmc.Tally
            Tally to merge with this one

        Returns
        -------
        merged_tally : openmc.Tally
            Merged tallies

        """

        if not self.can_merge(other):
            msg = 'Unable to merge tally ID="{}" with "{}"'.format(
                other.id, self.id)
            raise ValueError(msg)

        # Create deep copy of tally to return as merged tally
        merged_tally = copy.deepcopy(self)

        # Differentiate Tally with a new auto-generated Tally ID
        merged_tally.id = None

        # Create deep copy of other tally to use for array concatenation
        other_copy = copy.deepcopy(other)

        # Identify if filters, nuclides and scores are mergeable and/or equal
        merge_filters = self._can_merge_filters(other)
        merge_nuclides = self._can_merge_nuclides(other)
        merge_scores = self._can_merge_scores(other)
        equal_filters = sorted(self.filters) == sorted(other.filters)
        equal_nuclides = sorted(self.nuclides) == sorted(other.nuclides)
        equal_scores = sorted(self.scores) == sorted(other.scores)

        # If two tallies can be merged along a filter's bins
        if merge_filters and not equal_filters:

            # Search for mergeable filters
            for i, filter1 in enumerate(self.filters):
                for filter2 in other.filters:
                    if filter1 != filter2 and filter1.can_merge(filter2):
                        other_copy._swap_filters(other_copy.filters[i], filter2)
                        merged_tally.filters[i] = filter1.merge(filter2)
                        join_right = filter1 < filter2
                        merge_axis = i
                        break

        # If two tallies can be merged along nuclide bins
        if merge_nuclides and not equal_nuclides:
            merge_axis = self.num_filters
            join_right = True

            # Add unique nuclides from other tally to merged tally
            for nuclide in other.nuclides:
                if nuclide not in merged_tally.nuclides:
                    merged_tally.nuclides.append(nuclide)

        # If two tallies can be merged along score bins
        if merge_scores and not equal_scores:
            merge_axis = self.num_filters + 1
            join_right = True

            # Add unique scores from other tally to merged tally
            for score in other.scores:
                if score not in merged_tally.scores:
                    merged_tally.scores.append(score)

        # Add triggers from other tally to merged tally
        for trigger in other.triggers:
            merged_tally.triggers.append(trigger)

        # If results have not been read, then return tally for input generation
        if self._results_read is None:
            return merged_tally
        # Otherwise, this is a derived tally which needs merged results arrays
        else:
            self._derived = True

        # Concatenate sum arrays if present in both tallies
        if self.sum is not None and other_copy.sum is not None:
            self_sum = self.get_reshaped_data(value='sum')
            other_sum = other_copy.get_reshaped_data(value='sum')

            if join_right:
                merged_sum = np.concatenate((self_sum, other_sum),
                                            axis=merge_axis)
            else:
                merged_sum = np.concatenate((other_sum, self_sum),
                                            axis=merge_axis)

            merged_tally._sum = np.reshape(merged_sum, merged_tally.shape)

        # Concatenate sum_sq arrays if present in both tallies
        if self.sum_sq is not None and other.sum_sq is not None:
            self_sum_sq = self.get_reshaped_data(value='sum_sq')
            other_sum_sq = other_copy.get_reshaped_data(value='sum_sq')

            if join_right:
                merged_sum_sq = np.concatenate((self_sum_sq, other_sum_sq),
                                               axis=merge_axis)
            else:
                merged_sum_sq = np.concatenate((other_sum_sq, self_sum_sq),
                                               axis=merge_axis)

            merged_tally._sum_sq = np.reshape(merged_sum_sq, merged_tally.shape)

        # Concatenate mean arrays if present in both tallies
        if self.mean is not None and other.mean is not None:
            self_mean = self.get_reshaped_data(value='mean')
            other_mean = other_copy.get_reshaped_data(value='mean')

            if join_right:
                merged_mean = np.concatenate((self_mean, other_mean),
                                             axis=merge_axis)
            else:
                merged_mean = np.concatenate((other_mean, self_mean),
                                             axis=merge_axis)

            merged_tally._mean = np.reshape(merged_mean, merged_tally.shape)

        # Concatenate std. dev. arrays if present in both tallies
        if self.std_dev is not None and other.std_dev is not None:
            self_std_dev = self.get_reshaped_data(value='std_dev')
            other_std_dev = other_copy.get_reshaped_data(value='std_dev')

            if join_right:
                merged_std_dev = np.concatenate((self_std_dev, other_std_dev),
                                                axis=merge_axis)
            else:
                merged_std_dev = np.concatenate((other_std_dev, self_std_dev),
                                                axis=merge_axis)

            merged_tally._std_dev = np.reshape(merged_std_dev, merged_tally.shape)

        # Sparsify merged tally if both tallies are sparse
        merged_tally.sparse = self.sparse and other.sparse

        return merged_tally

    def to_xml_element(self):
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
        if len(self.filters) > 0:
            subelement = ET.SubElement(element, "filters")
            subelement.text = ' '.join(str(f.id) for f in self.filters)

        # Optional Nuclides
        if self.nuclides:
            subelement = ET.SubElement(element, "nuclides")
            subelement.text = ' '.join(str(n) for n in self.nuclides)

        # Scores
        if len(self.scores) == 0:
            msg = 'Unable to get XML for Tally ID="{}" since it does not ' \
                  'contain any scores'.format(self.id)
            raise ValueError(msg)

        else:
            scores = ''
            for score in self.scores:
                scores += '{} '.format(score)

            subelement = ET.SubElement(element, "scores")
            subelement.text = scores.rstrip(' ')

        # Tally estimator type
        if self.estimator is not None:
            subelement = ET.SubElement(element, "estimator")
            subelement.text = self.estimator

        # Optional Triggers
        for trigger in self.triggers:
            trigger.get_trigger_xml(element)

        # Optional derivatives
        if self.derivative is not None:
            subelement = ET.SubElement(element, "derivative")
            subelement.text = str(self.derivative.id)

        return element

    def contains_filter(self, filter_type):
        """Looks for a filter in the tally that matches a specified type

        Parameters
        ----------
        filter_type : openmc.FilterMeta
            Type of the filter, e.g. MeshFilter

        Returns
        -------
        filter_found : bool
            True if the tally contains a filter of the requested type;
            otherwise false

        """
        for test_filter in self.filters:
            if type(test_filter) is filter_type:
                return True
        return False

    def find_filter(self, filter_type):
        """Return a filter in the tally that matches a specified type

        Parameters
        ----------
        filter_type : openmc.FilterMeta
            Type of the filter, e.g. MeshFilter

        Returns
        -------
        filter_found : openmc.Filter
            Filter from this tally with matching type, or None if no matching
            Filter is found

        Raises
        ------
        ValueError
            If no matching Filter is found

        """

        # Look through all of this Tally's Filters for the type requested
        for test_filter in self.filters:
            if type(test_filter) is filter_type:
                return test_filter

            # Also check to see if the desired filter is wrapped up in an
            # aggregate
            elif isinstance(test_filter, openmc.AggregateFilter):
                if isinstance(test_filter.aggregate_filter, filter_type):
                    return test_filter

        # If we did not find the Filter, throw an Exception
        msg = 'Unable to find filter type "{}" in Tally ID="{}"'.format(
            filter_type, self.id)
        raise ValueError(msg)

    def get_nuclide_index(self, nuclide):
        """Returns the index in the Tally's results array for a Nuclide bin

        Parameters
        ----------
        nuclide : str
            The name of the Nuclide (e.g., 'H1', 'U238')

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
        # Look for the user-requested nuclide in all of the Tally's Nuclides
        for i, test_nuclide in enumerate(self.nuclides):
            # If the Summary was linked, then values are Nuclide objects
            if isinstance(test_nuclide, openmc.Nuclide):
                if test_nuclide.name == nuclide:
                    return i

            # If the Summary has not been linked, then values are ZAIDs
            else:
                if test_nuclide == nuclide:
                    return i

        msg = ('Unable to get the nuclide index for Tally since "{}" '
               'is not one of the nuclides'.format(nuclide))
        raise KeyError(msg)

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
            msg = 'Unable to get the score index for Tally since "{}" ' \
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
        filters : Iterable of openmc.FilterMeta
            An iterable of filter types
            (e.g., [MeshFilter, EnergyFilter]; default is [])
        filter_bins : Iterable of tuple
            A list of tuples of filter bins corresponding to the filter_types
            parameter (e.g., [(1,), ((0., 0.625e-6),)]; default is []). Each
            tuple contains bins for the corresponding filter type in the filters
            parameter. Each bin is an integer ID for Material-, Surface-,
            Cell-, Cellborn-, and Universe- Filters. Each bin is an integer
            for the cell instance ID for DistribcellFilters. Each bin is a
            2-tuple of floats for Energy- and Energyout- Filters corresponding
            to the energy boundaries of the bin of interest. The bin is an
            (x,y,z) 3-tuple for MeshFilters corresponding to the mesh cell
            of interest. The order of the bins in the list must correspond to
            the filter_types parameter.

        Returns
        -------
        numpy.ndarray
            A NumPy array of the filter indices

        """

        cv.check_type('filters', filters, Iterable, openmc.FilterMeta)
        cv.check_type('filter_bins', filter_bins, Iterable, tuple)

        # If user did not specify any specific Filters, use them all
        if not filters:
            return np.arange(self.num_filter_bins)

        # Initialize empty list of indices for each bin in each Filter
        filter_indices = []

        # Loop over all of the Tally's Filters
        for i, self_filter in enumerate(self.filters):
            # If a user-requested Filter, get the user-requested bins
            for j, test_filter in enumerate(filters):
                if type(self_filter) is test_filter:
                    bins = filter_bins[j]
                    break
            else:
                # If not a user-requested Filter, get all bins
                if isinstance(self_filter, openmc.DistribcellFilter):
                    # Create list of cell instance IDs for distribcell Filters
                    bins = list(range(self_filter.num_bins))

                elif isinstance(self_filter, openmc.EnergyFunctionFilter):
                    # EnergyFunctionFilters don't have bins so just add a None
                    bins = [None]

                else:
                    # Create list of IDs for bins for all other filter types
                    bins = self_filter.bins

            # Add indices for each bin in this Filter to the list
            indices = np.array([self_filter.get_bin_index(b) for b in bins])
            filter_indices.append(indices)

            # Account for stride in each of the previous filters
            for indices in filter_indices[:i]:
                indices *= self_filter.num_bins

        # Apply outer product sum between all filter bin indices
        return list(map(sum, product(*filter_indices)))

    def get_nuclide_indices(self, nuclides):
        """Get indices into the nuclide axis of this tally's data arrays.

        This is a helper method for the Tally.get_values(...) method to
        extract tally data. This method returns the indices into the nuclide
        axis of the tally's data array (axis=1) for one or more nuclides.

        Parameters
        ----------
        nuclides : list of str
            A list of nuclide name strings
            (e.g., ['U235', 'U238']; default is [])

        Returns
        -------
        numpy.ndarray
            A NumPy array of the nuclide indices

        """

        cv.check_iterable_type('nuclides', nuclides, str)

        # If user did not specify any specific Nuclides, use them all
        if not nuclides:
            return np.arange(self.num_nuclides)

        # Determine the score indices from any of the requested scores
        nuclide_indices = np.zeros_like(nuclides, dtype=int)
        for i, nuclide in enumerate(nuclides):
            nuclide_indices[i] = self.get_nuclide_index(nuclide)
        return nuclide_indices

    def get_score_indices(self, scores):
        """Get indices into the score axis of this tally's data arrays.

        This is a helper method for the Tally.get_values(...) method to
        extract tally data. This method returns the indices into the score
        axis of the tally's data array (axis=2) for one or more scores.

        Parameters
        ----------
        scores : list of str or openmc.CrossScore
            A list of one or more score strings
            (e.g., ['absorption', 'nu-fission']; default is [])

        Returns
        -------
        numpy.ndarray
            A NumPy array of the score indices

        """

        for score in scores:
            if not isinstance(score, (str, openmc.CrossScore)):
                msg = 'Unable to get score indices for score "{}" in Tally ' \
                      'ID="{}" since it is not a string or CrossScore'\
                      .format(score, self.id)
                raise ValueError(msg)

        # Determine the score indices from any of the requested scores
        if scores:
            score_indices = np.zeros(len(scores), dtype=int)
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
        filters : Iterable of openmc.FilterMeta
            An iterable of filter types
            (e.g., [MeshFilter, EnergyFilter]; default is [])
        filter_bins : list of Iterables
            A list of tuples of filter bins corresponding to the filter_types
            parameter (e.g., [(1,), ((0., 0.625e-6),)]; default is []). Each
            tuple contains bins for the corresponding filter type in the filters
            parameter. Each bins is the integer ID for 'material', 'surface',
            'cell', 'cellborn', and 'universe' Filters. Each bin is an integer
            for the cell instance ID for 'distribcell' Filters. Each bin is a
            2-tuple of floats for 'energy' and 'energyout' filters corresponding
            to the energy boundaries of the bin of interest. The bin is an
            (x,y,z) 3-tuple for 'mesh' filters corresponding to the mesh cell
            of interest. The order of the bins in the list must correspond to
            the filter_types parameter.
        nuclides : list of str
            A list of nuclide name strings
            (e.g., ['U235', 'U238']; default is [])
        value : str
            A string for the type of value to return  - 'mean' (default),
            'std_dev', 'rel_err', 'sum', or 'sum_sq' are accepted

        Returns
        -------
        float or numpy.ndarray
            A scalar or NumPy array of the Tally data indexed in the order
            each filter, nuclide and score is listed in the parameters.

        Raises
        ------
        ValueError
            When this method is called before the Tally is populated with data,
            or the input parameters do not correspond to the Tally's attributes,
            e.g., if the score(s) do not match those in the Tally.

        """

        # Ensure that the tally has data
        if (value == 'mean' and self.mean is None) or \
           (value == 'std_dev' and self.std_dev is None) or \
           (value == 'rel_err' and self.mean is None) or \
           (value == 'sum' and self.sum is None) or \
           (value == 'sum_sq' and self.sum_sq is None):
            msg = 'The Tally ID="{}" has no data to return'.format(self.id)
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
            msg = 'Unable to return results from Tally ID="{}" since the ' \
                  'the requested value "{}" is not \'mean\', \'std_dev\', ' \
                  '\'rel_err\', \'sum\', or \'sum_sq\''.format(self.id, value)
            raise LookupError(msg)

        return data

    def get_pandas_dataframe(self, filters=True, nuclides=True, scores=True,
                             derivative=True, paths=True, float_format='{:.2e}'):
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
        derivative : bool
            Include columns with differential tally info (default is True).
        paths : bool, optional
            Construct columns for distribcell tally filters (default is True).
            The geometric information in the Summary object is embedded into a
            Multi-index column with a geometric "path" to each distribcell
            instance.
        float_format : str
            All floats in the DataFrame will be formatted using the given
            format string before printing.

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

        """

        # Ensure that the tally has data
        if self.mean is None or self.std_dev is None:
            msg = 'The Tally ID="{}" has no data to return'.format(self.id)
            raise KeyError(msg)

        # Initialize a pandas dataframe for the tally data
        df = pd.DataFrame()

        # Find the total length of the tally data array
        data_size = self.mean.size

        # Build DataFrame columns for filters if user requested them
        if filters:
            # Append each Filter's DataFrame to the overall DataFrame
            for f, stride in zip(self.filters, self.filter_strides):
                filter_df = f.get_pandas_dataframe(
                    data_size, stride, paths=paths)
                df = pd.concat([df, filter_df], axis=1)

        # Include DataFrame column for nuclides if user requested it
        if nuclides:
            nuclides = []
            column_name = 'nuclide'

            for nuclide in self.nuclides:
                if isinstance(nuclide, openmc.Nuclide):
                    nuclides.append(nuclide.name)
                elif isinstance(nuclide, openmc.AggregateNuclide):
                    nuclides.append(nuclide.name)
                    column_name = '{}(nuclide)'.format(nuclide.aggregate_op)
                else:
                    nuclides.append(nuclide)

            # Tile the nuclide bins into a DataFrame column
            nuclides = np.repeat(nuclides, len(self.scores))
            tile_factor = data_size / len(nuclides)
            df[column_name] = np.tile(nuclides, int(tile_factor))

        # Include column for scores if user requested it
        if scores:
            scores = []
            column_name = 'score'

            for score in self.scores:
                if isinstance(score, (str, openmc.CrossScore)):
                    scores.append(str(score))
                elif isinstance(score, openmc.AggregateScore):
                    scores.append(score.name)
                    column_name = '{}(score)'.format(score.aggregate_op)

            tile_factor = data_size / len(self.scores)
            df[column_name] = np.tile(scores, int(tile_factor))

        # Include columns for derivatives if user requested it
        if derivative and (self.derivative is not None):
            df['d_variable'] = self.derivative.variable
            if self.derivative.material is not None:
                df['d_material'] = self.derivative.material
            if self.derivative.nuclide is not None:
                df['d_nuclide'] = self.derivative.nuclide

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

            # Create and set a MultiIndex for the DataFrame's columns, but only
            # if any column actually is multi-level (e.g., a mesh filter)
            if any(len(c) > 1 for c in columns):
                df.columns = pd.MultiIndex.from_tuples(columns)

        # Modify the df.to_string method so that it prints formatted strings.
        # Credit to http://stackoverflow.com/users/3657742/chrisb for this trick
        df.to_string = partial(df.to_string, float_format=float_format.format)

        return df

    def get_reshaped_data(self, value='mean'):
        """Returns an array of tally data with one dimension per filter.

        The tally data in OpenMC is stored as a 3D array with the dimensions
        corresponding to filters, nuclides and scores. As a result, tally data
        can be opaque for a user to directly index (i.e., without use of
        :meth:`openmc.Tally.get_values`) since one must know how to properly use
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
        numpy.ndarray
            The tally data array indexed by filters, nuclides and scores.

        """

        # Get the 3D array of data in filters, nuclides and scores
        data = self.get_values(value=value)

        # Build a new array shape with one dimension per filter
        new_shape = tuple(f.num_bins for f in self.filters)
        new_shape += (self.num_nuclides, self.num_scores)

        # Reshape the data with one dimension for each filter
        data = np.reshape(data, new_shape)
        return data

    def hybrid_product(self, other, binary_op, filter_product=None,
                       nuclide_product=None, score_product=None):
        """Combines filters, scores and nuclides with another tally.

        This is a helper method for the tally arithmetic operator overloaded
        methods. It is called a "hybrid product" because it performs a
        combination of tensor (or Kronecker) and entrywise (or Hadamard)
        products. The filters from both tallies are combined using an entrywise
        (or Hadamard) product on matching filters. By default, if all nuclides
        are identical in the two tallies, the entrywise product is performed
        across nuclides; else the tensor product is performed. By default, if
        all scores are identical in the two tallies, the entrywise product is
        performed across scores; else the tensor product is performed. Users
        can also call the method explicitly and specify the desired product.

        Parameters
        ----------
        other : openmc.Tally
            The tally on the right hand side of the hybrid product
        binary_op : {'+', '-', '*', '/', '^'}
            The binary operation in the hybrid product
        filter_product : {'tensor', 'entrywise' or None}
            The type of product (tensor or entrywise) to be performed between
            filter data. The default is the entrywise product. Currently only
            the entrywise product is supported since a tally cannot contain
            two of the same filter.
        nuclide_product : {'tensor', 'entrywise' or None}
            The type of product (tensor or entrywise) to be performed between
            nuclide data. The default is the entrywise product if all nuclides
            between the two tallies are the same; otherwise the default is
            the tensor product.
        score_product : {'tensor', 'entrywise' or None}
            The type of product (tensor or entrywise) to be performed between
            score data. The default is the entrywise product if all scores
            between the two tallies are the same; otherwise the default is
            the tensor product.

        Returns
        -------
        openmc.Tally
            A new Tally that is the hybrid product with this one.

        Raises
        ------
        ValueError
            When this method is called before the other tally is populated
            with data.

        """

        # Set default value for filter product if it was not set
        if filter_product is None:
            filter_product = 'entrywise'
        elif filter_product == 'tensor':
            msg = 'Unable to perform Tally arithmetic with a tensor product' \
                  'for the filter data as this is not currently supported.'
            raise ValueError(msg)

        # Set default value for nuclide product if it was not set
        if nuclide_product is None:
            if self.nuclides == other.nuclides:
                nuclide_product = 'entrywise'
            else:
                nuclide_product = 'tensor'

        # Set default value for score product if it was not set
        if score_product is None:
            if self.scores == other.scores:
                score_product = 'entrywise'
            else:
                score_product = 'tensor'

        # Check product types
        cv.check_value('filter product', filter_product, _PRODUCT_TYPES)
        cv.check_value('nuclide product', nuclide_product, _PRODUCT_TYPES)
        cv.check_value('score product', score_product, _PRODUCT_TYPES)

        # Check that results have been read
        if not other.derived and other.sum is None:
            msg = 'Unable to use tally arithmetic with Tally ID="{}" ' \
                  'since it does not contain any results.'.format(other.id)
            raise ValueError(msg)

        new_tally = Tally()
        new_tally._derived = True
        new_tally.with_batch_statistics = True
        new_tally._num_realizations = self.num_realizations
        new_tally._estimator = self.estimator
        new_tally._with_summary = self.with_summary
        new_tally._sp_filename = self._sp_filename

        # Construct a combined derived name from the two tally operands
        if self.name != '' and other.name != '':
            new_name = '({} {} {})'.format(self.name, binary_op, other.name)
            new_tally.name = new_name

        # Query the mean and std dev so the tally data is read in from file
        # if it has not already been read in.
        self.mean, self.std_dev, other.mean, other.std_dev

        # Create copies of self and other tallies to rearrange for tally
        # arithmetic
        self_copy = copy.deepcopy(self)
        other_copy = copy.deepcopy(other)

        self_copy.sparse = False
        other_copy.sparse = False

        # Align the tally data based on desired hybrid product
        data = self_copy._align_tally_data(other_copy, filter_product,
                                           nuclide_product, score_product)

        # Perform tally arithmetic operation
        if binary_op == '+':
            new_tally._mean = data['self']['mean'] + data['other']['mean']
            new_tally._std_dev = np.sqrt(data['self']['std. dev.']**2 +
                                         data['other']['std. dev.']**2)
        elif binary_op == '-':
            new_tally._mean = data['self']['mean'] - data['other']['mean']
            new_tally._std_dev = np.sqrt(data['self']['std. dev.']**2 +
                                         data['other']['std. dev.']**2)
        elif binary_op == '*':
            with np.errstate(divide='ignore', invalid='ignore'):
                self_rel_err = data['self']['std. dev.'] / data['self']['mean']
                other_rel_err = data['other']['std. dev.'] / data['other']['mean']
            new_tally._mean = data['self']['mean'] * data['other']['mean']
            new_tally._std_dev = np.abs(new_tally.mean) * \
                                 np.sqrt(self_rel_err**2 + other_rel_err**2)
        elif binary_op == '/':
            with np.errstate(divide='ignore', invalid='ignore'):
                self_rel_err = data['self']['std. dev.'] / data['self']['mean']
                other_rel_err = data['other']['std. dev.'] / data['other']['mean']
                new_tally._mean = data['self']['mean'] / data['other']['mean']
            new_tally._std_dev = np.abs(new_tally.mean) * \
                                 np.sqrt(self_rel_err**2 + other_rel_err**2)
        elif binary_op == '^':
            with np.errstate(divide='ignore', invalid='ignore'):
                mean_ratio = data['other']['mean'] / data['self']['mean']
            first_term = mean_ratio * data['self']['std. dev.']
            second_term = \
                np.log(data['self']['mean']) * data['other']['std. dev.']
            new_tally._mean = data['self']['mean'] ** data['other']['mean']
            new_tally._std_dev = np.abs(new_tally.mean) * \
                                 np.sqrt(first_term**2 + second_term**2)

        # Convert any infs and nans to zero
        new_tally._mean[np.isinf(new_tally._mean)] = 0
        new_tally._mean = np.nan_to_num(new_tally._mean)
        new_tally._std_dev[np.isinf(new_tally._std_dev)] = 0
        new_tally._std_dev = np.nan_to_num(new_tally._std_dev)

        # Set tally attributes
        if self_copy.estimator == other_copy.estimator:
            new_tally.estimator = self_copy.estimator
        if self_copy.with_summary and other_copy.with_summary:
            new_tally.with_summary = self_copy.with_summary
        if self_copy.num_realizations == other_copy.num_realizations:
            new_tally.num_realizations = self_copy.num_realizations

        # Add filters to the new tally
        if filter_product == 'entrywise':
            for self_filter in self_copy.filters:
                new_tally.filters.append(self_filter)
        else:
            all_filters = [self_copy.filters, other_copy.filters]
            for self_filter, other_filter in product(*all_filters):
                new_filter = openmc.CrossFilter(self_filter, other_filter,
                                                binary_op)
                new_tally.filters.append(new_filter)

        # Add nuclides to the new tally
        if nuclide_product == 'entrywise':
            for self_nuclide in self_copy.nuclides:
                new_tally.nuclides.append(self_nuclide)
        else:
            all_nuclides = [self_copy.nuclides, other_copy.nuclides]
            for self_nuclide, other_nuclide in product(*all_nuclides):
                new_nuclide = openmc.CrossNuclide(self_nuclide, other_nuclide,
                                                  binary_op)
                new_tally.nuclides.append(new_nuclide)

        # Define helper function that handles score units appropriately
        # depending on the binary operator
        def cross_score(score1, score2, binary_op):
            if binary_op == '+' or binary_op == '-':
                if score1 == score2:
                    return score1
                else:
                    return openmc.CrossScore(score1, score2, binary_op)
            else:
                return openmc.CrossScore(score1, score2, binary_op)

        # Add scores to the new tally
        if score_product == 'entrywise':
            for self_score in self_copy.scores:
                new_score = cross_score(self_score, self_score, binary_op)
                new_tally.scores.append(new_score)
        else:
            all_scores = [self_copy.scores, other_copy.scores]
            for self_score, other_score in product(*all_scores):
                new_score = cross_score(self_score, other_score, binary_op)
                new_tally.scores.append(new_score)

        return new_tally

    @property
    def filter_strides(self):
        all_strides = []
        stride = self.num_nuclides * self.num_scores
        for self_filter in reversed(self.filters):
            all_strides.append(stride)
            stride *= self_filter.num_bins
        return all_strides[::-1]

    def _align_tally_data(self, other, filter_product, nuclide_product,
                          score_product):
        """Aligns data from two tallies for tally arithmetic.

        This is a helper method to construct a dict of dicts of the "aligned"
        data arrays from each tally for tally arithmetic. The method analyzes
        the filters, scores and nuclides in both tallies and determines how to
        appropriately align the data for vectorized arithmetic. For example,
        if the two tallies have different filters, this method will use NumPy
        'tile' and 'repeat' operations to the new data arrays such that all
        possible combinations of the data in each tally's bins will be made
        when the arithmetic operation is applied to the arrays.

        Parameters
        ----------
        other : openmc.Tally
            The tally to outer product with this tally
        filter_product : {'entrywise'}
            The type of product to be performed between filter data. Currently,
            only the entrywise product is supported for the filter product.
        nuclide_product : {'tensor', 'entrywise'}
            The type of product (tensor or entrywise) to be performed between
            nuclide data.
        score_product : {'tensor', 'entrywise'}
            The type of product (tensor or entrywise) to be performed between
            score data.

        Returns
        -------
        dict
            A dictionary of dictionaries to "aligned" 'mean' and 'std. dev'
            NumPy arrays for each tally's data.

        """

        # Get the set of filters that each tally is missing
        other_missing_filters = set(self.filters) - set(other.filters)
        self_missing_filters = set(other.filters) - set(self.filters)

        # Add filters present in self but not in other to other
        for other_filter in other_missing_filters:
            filter_copy = copy.deepcopy(other_filter)
            other._mean = np.repeat(other.mean, filter_copy.num_bins, axis=0)
            other._std_dev = np.repeat(other.std_dev, filter_copy.num_bins, axis=0)
            other.filters.append(filter_copy)

        # Add filters present in other but not in self to self
        for self_filter in self_missing_filters:
            filter_copy = copy.deepcopy(self_filter)
            self._mean = np.repeat(self.mean, filter_copy.num_bins, axis=0)
            self._std_dev = np.repeat(self.std_dev, filter_copy.num_bins, axis=0)
            self.filters.append(filter_copy)

        # Align other filters with self filters
        for i, self_filter in enumerate(self.filters):
            other_index = other.filters.index(self_filter)

            # If necessary, swap other filter
            if other_index != i:
                other._swap_filters(self_filter, other.filters[i])

        # Repeat and tile the data by nuclide in preparation for performing
        # the tensor product across nuclides.
        if nuclide_product == 'tensor':
            self._mean = np.repeat(self.mean, other.num_nuclides, axis=1)
            self._std_dev = np.repeat(self.std_dev, other.num_nuclides, axis=1)
            other._mean = np.tile(other.mean, (1, self.num_nuclides, 1))
            other._std_dev = np.tile(other.std_dev, (1, self.num_nuclides, 1))

        # Add nuclides to each tally such that each tally contains the complete
        # set of nuclides necessary to perform an entrywise product. New
        # nuclides added to a tally will have all their scores set to zero.
        else:

            # Get the set of nuclides that each tally is missing
            other_missing_nuclides = set(self.nuclides) - set(other.nuclides)
            self_missing_nuclides = set(other.nuclides) - set(self.nuclides)

            # Add nuclides present in self but not in other to other
            for nuclide in other_missing_nuclides:
                other._mean = np.insert(other.mean, other.num_nuclides, 0, axis=1)
                other._std_dev = np.insert(other.std_dev, other.num_nuclides, 0,
                                           axis=1)
                other.nuclides.append(nuclide)

            # Add nuclides present in other but not in self to self
            for nuclide in self_missing_nuclides:
                self._mean = np.insert(self.mean, self.num_nuclides, 0, axis=1)
                self._std_dev = np.insert(self.std_dev, self.num_nuclides, 0,
                                          axis=1)
                self.nuclides.append(nuclide)

            # Align other nuclides with self nuclides
            for i, nuclide in enumerate(self.nuclides):
                other_index = other.get_nuclide_index(nuclide)

                # If necessary, swap other nuclide
                if other_index != i:
                    other._swap_nuclides(nuclide, other.nuclides[i])

        # Repeat and tile the data by score in preparation for performing
        # the tensor product across scores.
        if score_product == 'tensor':
            self._mean = np.repeat(self.mean, other.num_scores, axis=2)
            self._std_dev = np.repeat(self.std_dev, other.num_scores, axis=2)
            other._mean = np.tile(other.mean, (1, 1, self.num_scores))
            other._std_dev = np.tile(other.std_dev, (1, 1, self.num_scores))

        # Add scores to each tally such that each tally contains the complete set
        # of scores necessary to perform an entrywise product. New scores added
        # to a tally will be set to zero.
        else:

            # Get the set of scores that each tally is missing
            other_missing_scores = set(self.scores) - set(other.scores)
            self_missing_scores = set(other.scores) - set(self.scores)

            # Add scores present in self but not in other to other
            for score in other_missing_scores:
                other._mean = np.insert(other.mean, other.num_scores, 0, axis=2)
                other._std_dev = np.insert(other.std_dev, other.num_scores, 0, axis=2)
                other.scores.append(score)

            # Add scores present in other but not in self to self
            for score in self_missing_scores:
                self._mean = np.insert(self.mean, self.num_scores, 0, axis=2)
                self._std_dev = np.insert(self.std_dev, self.num_scores, 0, axis=2)
                self.scores.append(score)

            # Align other scores with self scores
            for i, score in enumerate(self.scores):
                other_index = other.scores.index(score)

                # If necessary, swap other score
                if other_index != i:
                    other._swap_scores(score, other.scores[i])

        data = {}
        data['self'] = {}
        data['other'] = {}
        data['self']['mean'] = self.mean
        data['other']['mean'] = other.mean
        data['self']['std. dev.'] = self.std_dev
        data['other']['std. dev.'] = other.std_dev
        return data

    def _swap_filters(self, filter1, filter2):
        """Reverse the ordering of two filters in this tally

        This is a helper method for tally arithmetic which helps align the data
        in two tallies with shared filters. This method reverses the order of
        the two filters in place.

        Parameters
        ----------
        filter1 : Filter
            The filter to swap with filter2
        filter2 : Filter
            The filter to swap with filter1

        Raises
        ------
        ValueError
            If this is a derived tally or this method is called before the tally
            is populated with data.

        """

        cv.check_type('filter1', filter1, _FILTER_CLASSES)
        cv.check_type('filter2', filter2, _FILTER_CLASSES)

        # Check that the filters exist in the tally and are not the same
        if filter1 == filter2:
            return
        elif filter1 not in self.filters:
            msg = 'Unable to swap "{}" filter1 in Tally ID="{}" since it ' \
                  'does not contain such a filter'.format(filter1.type, self.id)
            raise ValueError(msg)
        elif filter2 not in self.filters:
            msg = 'Unable to swap "{}" filter2 in Tally ID="{}" since it ' \
                  'does not contain such a filter'.format(filter2.type, self.id)
            raise ValueError(msg)

        # Construct lists of tuples for the bins in each of the two filters
        filters = [type(filter1), type(filter2)]
        if isinstance(filter1, openmc.DistribcellFilter):
            filter1_bins = [b for b in range(filter1.num_bins)]
        elif isinstance(filter1, openmc.EnergyFunctionFilter):
            filter1_bins = [None]
        else:
            filter1_bins = filter1.bins

        if isinstance(filter2, openmc.DistribcellFilter):
            filter2_bins = [b for b in range(filter2.num_bins)]
        elif isinstance(filter2, openmc.EnergyFunctionFilter):
            filter2_bins = [None]
        else:
            filter2_bins = filter2.bins

        # Create variables to store views of data in the misaligned structure
        mean = {}
        std_dev = {}

        # Store the data from the misaligned structure
        for i, (bin1, bin2) in enumerate(product(filter1_bins, filter2_bins)):
            filter_bins = [(bin1,), (bin2,)]

            if self.mean is not None:
                mean[i] = self.get_values(
                    filters=filters, filter_bins=filter_bins, value='mean')

            if self.std_dev is not None:
                std_dev[i] = self.get_values(
                    filters=filters, filter_bins=filter_bins, value='std_dev')

        # Swap the filters in the copied version of this Tally
        filter1_index = self.filters.index(filter1)
        filter2_index = self.filters.index(filter2)
        self.filters[filter1_index] = filter2
        self.filters[filter2_index] = filter1

        # Realign the data
        for i, (bin1, bin2) in enumerate(product(filter1_bins, filter2_bins)):
            filter_bins = [(bin1,), (bin2,)]
            indices = self.get_filter_indices(filters, filter_bins)

            if self.mean is not None:
                self.mean[indices, :, :] = mean[i]

            if self.std_dev is not None:
                self.std_dev[indices, :, :] = std_dev[i]

    def _swap_nuclides(self, nuclide1, nuclide2):
        """Reverse the ordering of two nuclides in this tally

        This is a helper method for tally arithmetic which helps align the data
        in two tallies with shared nuclides. This method reverses the order of
        the two nuclides in place.

        Parameters
        ----------
        nuclide1 : Nuclide
            The nuclide to swap with nuclide2

        nuclide2 : Nuclide
            The nuclide to swap with nuclide1

        Raises
        ------
        ValueError
            If this is a derived tally or this method is called before the tally
            is populated with data.

        """

        # Check that results have been read
        if not self.derived and self.sum is None:
            msg = 'Unable to use tally arithmetic with Tally ID="{}" ' \
                  'since it does not contain any results.'.format(self.id)
            raise ValueError(msg)

        cv.check_type('nuclide1', nuclide1, _NUCLIDE_CLASSES)
        cv.check_type('nuclide2', nuclide2, _NUCLIDE_CLASSES)

        # Check that the nuclides exist in the tally and are not the same
        if nuclide1 == nuclide2:
            msg = 'Unable to swap a nuclide with itself'
            raise ValueError(msg)
        elif nuclide1 not in self.nuclides:
            msg = 'Unable to swap nuclide1 "{}" in Tally ID="{}" since it ' \
                  'does not contain such a nuclide'\
                  .format(nuclide1.name, self.id)
            raise ValueError(msg)
        elif nuclide2 not in self.nuclides:
            msg = 'Unable to swap "{}" nuclide2 in Tally ID="{}" since it ' \
                  'does not contain such a nuclide'\
                  .format(nuclide2.name, self.id)
            raise ValueError(msg)

        # Swap the nuclides in the Tally
        nuclide1_index = self.get_nuclide_index(nuclide1)
        nuclide2_index = self.get_nuclide_index(nuclide2)
        self.nuclides[nuclide1_index] = nuclide2
        self.nuclides[nuclide2_index] = nuclide1

        # Adjust the mean data array to relect the new nuclide order
        if self.mean is not None:
            nuclide1_mean = self.mean[:, nuclide1_index, :].copy()
            nuclide2_mean = self.mean[:, nuclide2_index, :].copy()
            self.mean[:, nuclide2_index, :] = nuclide1_mean
            self.mean[:, nuclide1_index, :] = nuclide2_mean

        # Adjust the std_dev data array to relect the new nuclide order
        if self.std_dev is not None:
            nuclide1_std_dev = self.std_dev[:, nuclide1_index, :].copy()
            nuclide2_std_dev = self.std_dev[:, nuclide2_index, :].copy()
            self.std_dev[:, nuclide2_index, :] = nuclide1_std_dev
            self.std_dev[:, nuclide1_index, :] = nuclide2_std_dev

    def _swap_scores(self, score1, score2):
        """Reverse the ordering of two scores in this tally

        This is a helper method for tally arithmetic which helps align the data
        in two tallies with shared scores. This method reverses the order
        of the two scores in place.

        Parameters
        ----------
        score1 : str or CrossScore
            The score to swap with score2

        score2 : str or CrossScore
            The score to swap with score1

        Raises
        ------
        ValueError
            If this is a derived tally or this method is called before the tally
            is populated with data.

        """

        # Check that results have been read
        if not self.derived and self.sum is None:
            msg = 'Unable to use tally arithmetic with Tally ID="{}" ' \
                  'since it does not contain any results.'.format(self.id)
            raise ValueError(msg)

        # Check that the scores are valid
        if not isinstance(score1, (str, openmc.CrossScore)):
            msg = 'Unable to swap score1 "{}" in Tally ID="{}" since it is ' \
                  'not a string or CrossScore'.format(score1, self.id)
            raise ValueError(msg)
        elif not isinstance(score2, (str, openmc.CrossScore)):
            msg = 'Unable to swap score2 "{}" in Tally ID="{}" since it is ' \
                  'not a string or CrossScore'.format(score2, self.id)
            raise ValueError(msg)

        # Check that the scores exist in the tally and are not the same
        if score1 == score2:
            msg = 'Unable to swap a score with itself'
            raise ValueError(msg)
        elif score1 not in self.scores:
            msg = 'Unable to swap score1 "{}" in Tally ID="{}" since it ' \
                  'does not contain such a score'.format(score1, self.id)
            raise ValueError(msg)
        elif score2 not in self.scores:
            msg = 'Unable to swap score2 "{}" in Tally ID="{}" since it ' \
                  'does not contain such a score'.format(score2, self.id)
            raise ValueError(msg)

        # Swap the scores in the Tally
        score1_index = self.get_score_index(score1)
        score2_index = self.get_score_index(score2)
        self.scores[score1_index] = score2
        self.scores[score2_index] = score1

        # Adjust the mean data array to relect the new nuclide order
        if self.mean is not None:
            score1_mean = self.mean[:, :, score1_index].copy()
            score2_mean = self.mean[:, :, score2_index].copy()
            self.mean[:, :, score2_index] = score1_mean
            self.mean[:, :, score1_index] = score2_mean

        # Adjust the std_dev data array to relect the new nuclide order
        if self.std_dev is not None:
            score1_std_dev = self.std_dev[:, :, score1_index].copy()
            score2_std_dev = self.std_dev[:, :, score2_index].copy()
            self.std_dev[:, :, score2_index] = score1_std_dev
            self.std_dev[:, :, score1_index] = score2_std_dev

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
        other : openmc.Tally or float
            The tally or scalar value to add to this tally

        Returns
        -------
        openmc.Tally
            A new derived tally which is the sum of this tally and the other
            tally or scalar value in the addition.

        Raises
        ------
        ValueError
            When this method is called before the Tally is populated with data.

        """

        # Check that results have been read
        if not self.derived and self.sum is None:
            msg = 'Unable to use tally arithmetic with Tally ID="{}" ' \
                  'since it does not contain any results.'.format(self.id)
            raise ValueError(msg)

        if isinstance(other, Tally):
            new_tally = self.hybrid_product(other, binary_op='+')

            # If both tally operands were sparse, sparsify the new tally
            if self.sparse and other.sparse:
                new_tally.sparse = True

        elif isinstance(other, Real):
            new_tally = Tally(name='derived')
            new_tally._derived = True
            new_tally.with_batch_statistics = True
            new_tally.name = self.name
            new_tally._mean = self.mean + other
            new_tally._std_dev = self.std_dev
            new_tally.estimator = self.estimator
            new_tally.with_summary = self.with_summary
            new_tally.num_realizations = self.num_realizations

            new_tally.filters = copy.deepcopy(self.filters)
            new_tally.nuclides = copy.deepcopy(self.nuclides)
            new_tally.scores = copy.deepcopy(self.scores)

            # If this tally operand is sparse, sparsify the new tally
            new_tally.sparse = self.sparse

        else:
            msg = 'Unable to add "{}" to Tally ID="{}"'.format(other, self.id)
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
        other : openmc.Tally or float
            The tally or scalar value to subtract from this tally

        Returns
        -------
        openmc.Tally
            A new derived tally which is the difference of this tally and the
            other tally or scalar value in the subtraction.

        Raises
        ------
        ValueError
            When this method is called before the Tally is populated with data.

        """

        # Check that results have been read
        if not self.derived and self.sum is None:
            msg = 'Unable to use tally arithmetic with Tally ID="{}" ' \
                  'since it does not contain any results.'.format(self.id)
            raise ValueError(msg)

        if isinstance(other, Tally):
            new_tally = self.hybrid_product(other, binary_op='-')

            # If both tally operands were sparse, sparsify the new tally
            if self.sparse and other.sparse:
                new_tally.sparse = True

        elif isinstance(other, Real):
            new_tally = Tally(name='derived')
            new_tally._derived = True
            new_tally.name = self.name
            new_tally._mean = self.mean - other
            new_tally._std_dev = self.std_dev
            new_tally.estimator = self.estimator
            new_tally.with_summary = self.with_summary
            new_tally.num_realizations = self.num_realizations

            new_tally.filters = copy.deepcopy(self.filters)
            new_tally.nuclides = copy.deepcopy(self.nuclides)
            new_tally.scores = copy.deepcopy(self.scores)

            # If this tally operand is sparse, sparsify the new tally
            new_tally.sparse = self.sparse

        else:
            msg = 'Unable to subtract "{}" from Tally ID="{}"'.format(other, self.id)
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
        other : openmc.Tally or float
            The tally or scalar value to multiply with this tally

        Returns
        -------
        openmc.Tally
            A new derived tally which is the product of this tally and the
            other tally or scalar value in the multiplication.

        Raises
        ------
        ValueError
            When this method is called before the Tally is populated with data.

        """

        # Check that results have been read
        if not self.derived and self.sum is None:
            msg = 'Unable to use tally arithmetic with Tally ID="{}" ' \
                  'since it does not contain any results.'.format(self.id)
            raise ValueError(msg)

        if isinstance(other, Tally):
            new_tally = self.hybrid_product(other, binary_op='*')

            # If original tally operands were sparse, sparsify the new tally
            if self.sparse and other.sparse:
                new_tally.sparse = True

        elif isinstance(other, Real):
            new_tally = Tally(name='derived')
            new_tally._derived = True
            new_tally.name = self.name
            new_tally._mean = self.mean * other
            new_tally._std_dev = self.std_dev * np.abs(other)
            new_tally.estimator = self.estimator
            new_tally.with_summary = self.with_summary
            new_tally.num_realizations = self.num_realizations

            new_tally.filters = copy.deepcopy(self.filters)
            new_tally.nuclides = copy.deepcopy(self.nuclides)
            new_tally.scores = copy.deepcopy(self.scores)

            # If this tally operand is sparse, sparsify the new tally
            new_tally.sparse = self.sparse

        else:
            msg = 'Unable to multiply Tally ID="{}" by "{}"'.format(self.id, other)
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
        other : openmc.Tally or float
            The tally or scalar value to divide this tally by

        Returns
        -------
        openmc.Tally
            A new derived tally which is the dividend of this tally and the
            other tally or scalar value in the division.

        Raises
        ------
        ValueError
            When this method is called before the Tally is populated with data.

        """

        # Check that results have been read
        if not self.derived and self.sum is None:
            msg = 'Unable to use tally arithmetic with Tally ID="{}" ' \
                  'since it does not contain any results.'.format(self.id)
            raise ValueError(msg)

        if isinstance(other, Tally):
            new_tally = self.hybrid_product(other, binary_op='/')

            # If original tally operands were sparse, sparsify the new tally
            if self.sparse and other.sparse:
                new_tally.sparse = True

        elif isinstance(other, Real):
            new_tally = Tally(name='derived')
            new_tally._derived = True
            new_tally.name = self.name
            new_tally._mean = self.mean / other
            new_tally._std_dev = self.std_dev * np.abs(1. / other)
            new_tally.estimator = self.estimator
            new_tally.with_summary = self.with_summary
            new_tally.num_realizations = self.num_realizations

            new_tally.filters = copy.deepcopy(self.filters)
            new_tally.nuclides = copy.deepcopy(self.nuclides)
            new_tally.scores = copy.deepcopy(self.scores)

            # If this tally operand is sparse, sparsify the new tally
            new_tally.sparse = self.sparse

        else:
            msg = 'Unable to divide Tally ID="{}" by "{}"'.format(self.id, other)
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
        power : openmc.Tally or float
            The tally or scalar value exponent

        Returns
        -------
        openmc.Tally
            A new derived tally which is this tally raised to the power of the
            other tally or scalar value in the exponentiation.

        Raises
        ------
        ValueError
            When this method is called before the Tally is populated with data.

        """

        # Check that results have been read
        if not self.derived and self.sum is None:
            msg = 'Unable to use tally arithmetic with Tally ID="{}" ' \
                  'since it does not contain any results.'.format(self.id)
            raise ValueError(msg)

        if isinstance(power, Tally):
            new_tally = self.hybrid_product(power, binary_op='^')

            # If original tally operand was sparse, sparsify the new tally
            if self.sparse:
                new_tally.sparse = True

        elif isinstance(power, Real):
            new_tally = Tally(name='derived')
            new_tally._derived = True
            new_tally.name = self.name
            new_tally._mean = self._mean ** power
            self_rel_err = self.std_dev / self.mean
            new_tally._std_dev = np.abs(new_tally._mean * power * self_rel_err)
            new_tally.estimator = self.estimator
            new_tally.with_summary = self.with_summary
            new_tally.num_realizations = self.num_realizations

            new_tally.filters = copy.deepcopy(self.filters)
            new_tally.nuclides = copy.deepcopy(self.nuclides)
            new_tally.scores = copy.deepcopy(self.scores)

            # If original tally was sparse, sparsify the exponentiated tally
            new_tally.sparse = self.sparse

        else:
            msg = 'Unable to raise Tally ID="{}" to power "{}"'.format(self.id, power)
            raise ValueError(msg)

        return new_tally

    def __radd__(self, other):
        """Right addition with a scalar value.

        This reverses the operands and calls the __add__ method.

        Parameters
        ----------
        other : float
            The scalar value to add to this tally

        Returns
        -------
        openmc.Tally
            A new derived tally of this tally added with the scalar value.

        """

        return self + other

    def __rsub__(self, other):
        """Right subtraction from a scalar value.

        This reverses the operands and calls the __sub__ method.

        Parameters
        ----------
        other : float
            The scalar value to subtract this tally from

        Returns
        -------
        openmc.Tally
            A new derived tally of this tally subtracted from the scalar value.

        """

        return -1. * self + other

    def __rmul__(self, other):
        """Right multiplication with a scalar value.

        This reverses the operands and calls the __mul__ method.

        Parameters
        ----------
        other : float
            The scalar value to multiply with this tally

        Returns
        -------
        openmc.Tally
            A new derived tally of this tally multiplied by the scalar value.

        """

        return self * other

    def __rdiv__(self, other):
        """Right division with a scalar value.

        This reverses the operands and calls the __div__ method.

        Parameters
        ----------
        other : float
            The scalar value to divide by this tally

        Returns
        -------
        openmc.Tally
            A new derived tally of the scalar value divided by this tally.

        """

        return other * self**-1

    def __abs__(self):
        """The absolute value of this tally.

        Returns
        -------
        openmc.Tally
            A new derived tally which is the absolute value of this tally.

        """

        new_tally = copy.deepcopy(self)
        new_tally._mean = np.abs(new_tally.mean)
        return new_tally

    def __neg__(self):
        """The negated value of this tally.

        Returns
        -------
        openmc.Tally
            A new derived tally which is the negated value of this tally.

        """

        new_tally = self * -1
        return new_tally

    def get_slice(self, scores=[], filters=[], filter_bins=[], nuclides=[],
                  squeeze=False):
        """Build a sliced tally for the specified filters, scores and nuclides.

        This method constructs a new tally to encapsulate a subset of the data
        represented by this tally. The subset of data to include in the tally
        slice is determined by the scores, filters and nuclides specified in
        the input parameters.

        Parameters
        ----------
        scores : list of str
            A list of one or more score strings (e.g., ['absorption',
            'nu-fission']
        filters : Iterable of openmc.FilterMeta
            An iterable of filter types (e.g., [MeshFilter, EnergyFilter])
        filter_bins : list of Iterables
            A list of iterables of filter bins corresponding to the specified
            filter types (e.g., [(1,), ((0., 0.625e-6),)]). Each iterable
            contains bins to slice for the corresponding filter type in the
            filters parameter. Each bin is the integer ID for 'material',
            'surface', 'cell', 'cellborn', and 'universe' Filters. Each bin is
            an integer for the cell instance ID for 'distribcell' Filters. Each
            bin is a 2-tuple of floats for 'energy' and 'energyout' filters
            corresponding to the energy boundaries of the bin of interest. The
            bin is an (x,y,z) 3-tuple for 'mesh' filters corresponding to the
            mesh cell of interest. The order of the bins in the list must
            correspond to the `filters` argument.
        nuclides : list of str
            A list of nuclide name strings (e.g., ['U235', 'U238'])
        squeeze : bool
            Whether to remove filters with only a single bin in the sliced tally

        Returns
        -------
        openmc.Tally
            A new tally which encapsulates the subset of data requested in the
            order each filter, nuclide and score is listed in the parameters.

        Raises
        ------
        ValueError
            When this method is called before the Tally is populated with data.

        """

        # Ensure that the tally has data
        if not self.derived and self.sum is None:
            msg = 'Unable to use tally arithmetic with Tally ID="{}" ' \
                  'since it does not contain any results.'.format(self.id)
            raise ValueError(msg)

        # Create deep copy of tally to return as sliced tally
        new_tally = copy.deepcopy(self)
        new_tally._derived = True

        # Differentiate Tally with a new auto-generated Tally ID
        new_tally.id = None

        new_tally.sparse = False

        if not self.derived and self.sum is not None:
            new_sum = self.get_values(scores, filters, filter_bins,
                                      nuclides, 'sum')
            new_tally.sum = new_sum
        if not self.derived and self.sum_sq is not None:
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
                f = new_tally.find_filter(filter_type)

                # Remove filters with only a single bin if requested
                if squeeze:
                    if len(filter_bins[i]) == 1:
                        new_tally.filters.remove(f)
                        continue
                    else:
                        raise RuntimeError('Cannot remove sliced filter with '
                                           'more than one bin.')

                # Remove and/or reorder filter bins to user specifications
                bin_indices = [f.get_bin_index(b)
                               for b in filter_bins[i]]
                bin_indices = np.unique(bin_indices)

                # Set bins for sliced filter
                new_filter = copy.copy(f)
                new_filter.bins = [f.bins[i] for i in bin_indices]

                # Set number of bins manually for mesh/distribcell filters
                if filter_type is openmc.DistribcellFilter:
                    new_filter._num_bins = f._num_bins

                # Replace existing filter with new one
                for j, test_filter in enumerate(new_tally.filters):
                    if isinstance(test_filter, filter_type):
                        new_tally.filters[j] = new_filter

        # If original tally was sparse, sparsify the sliced tally
        new_tally.sparse = self.sparse
        return new_tally

    def summation(self, scores=[], filter_type=None,
                  filter_bins=[], nuclides=[], remove_filter=False):
        """Vectorized sum of tally data across scores, filter bins and/or
        nuclides using tally aggregation.

        This method constructs a new tally to encapsulate the sum of the data
        represented by the summation of the data in this tally. The tally data
        sum is determined by the scores, filter bins and nuclides specified
        in the input parameters.

        Parameters
        ----------
        scores : list of str
            A list of one or more score strings to sum across
            (e.g., ['absorption', 'nu-fission']; default is [])
        filter_type : openmc.FilterMeta
            Type of the filter, e.g. MeshFilter
        filter_bins : Iterable of int or tuple
            A list of the filter bins corresponding to the filter_type parameter
            Each bin in the list is the integer ID for 'material', 'surface',
            'cell', 'cellborn', and 'universe' Filters. Each bin is an integer
            for the cell instance ID for 'distribcell' Filters. Each bin is a
            2-tuple of floats for 'energy' and 'energyout' filters corresponding
            to the energy boundaries of the bin of interest. Each bin is an
            (x,y,z) 3-tuple for 'mesh' filters corresponding to the mesh cell of
            interest.
        nuclides : list of str
            A list of nuclide name strings to sum across
            (e.g., ['U235', 'U238']; default is [])
        remove_filter : bool
            If a filter is being summed over, this bool indicates whether to
            remove that filter in the returned tally. Default is False.

        Returns
        -------
        openmc.Tally
            A new tally which encapsulates the sum of data requested.
        """

        # Create new derived Tally for summation
        tally_sum = Tally()
        tally_sum._derived = True
        tally_sum._estimator = self.estimator
        tally_sum._num_realizations = self.num_realizations
        tally_sum._with_batch_statistics = self.with_batch_statistics
        tally_sum._with_summary = self.with_summary
        tally_sum._sp_filename = self._sp_filename
        tally_sum._results_read = self._results_read

        # Get tally data arrays reshaped with one dimension per filter
        mean = self.get_reshaped_data(value='mean')
        std_dev = self.get_reshaped_data(value='std_dev')

        # Sum across any filter bins specified by the user
        if isinstance(filter_type, openmc.FilterMeta):
            find_filter = self.find_filter(filter_type)

            # If user did not specify filter bins, sum across all bins
            if len(filter_bins) == 0:
                bin_indices = np.arange(find_filter.num_bins)

                if isinstance(find_filter, openmc.DistribcellFilter):
                    filter_bins = np.arange(find_filter.num_bins)
                elif isinstance(find_filter, openmc.EnergyFunctionFilter):
                    filter_bins = [None]
                else:
                    filter_bins = find_filter.bins

            # Only sum across bins specified by the user
            else:
                bin_indices = \
                    [find_filter.get_bin_index(bin) for bin in filter_bins]

            # Sum across the bins in the user-specified filter
            for i, self_filter in enumerate(self.filters):
                if type(self_filter) == filter_type:
                    shape = mean.shape
                    mean = np.take(mean, indices=bin_indices, axis=i)
                    std_dev = np.take(std_dev, indices=bin_indices, axis=i)

                    # NumPy take introduces a new dimension in output array
                    # for some special cases that must be removed
                    if len(mean.shape) > len(shape):
                        mean = np.squeeze(mean, axis=i)
                        std_dev = np.squeeze(std_dev, axis=i)

                    mean = np.sum(mean, axis=i, keepdims=True)
                    std_dev = np.sum(std_dev**2, axis=i, keepdims=True)
                    std_dev = np.sqrt(std_dev)

                    # Add AggregateFilter to the tally sum
                    if not remove_filter:
                        filter_sum = openmc.AggregateFilter(self_filter,
                            [tuple(filter_bins)], 'sum')
                        tally_sum.filters.append(filter_sum)

                # Add a copy of each filter not summed across to the tally sum
                else:
                    tally_sum.filters.append(copy.deepcopy(self_filter))

        # Add a copy of this tally's filters to the tally sum
        else:
            tally_sum._filters = copy.deepcopy(self.filters)

        # Sum across any nuclides specified by the user
        if len(nuclides) != 0:
            nuclide_bins = [self.get_nuclide_index(nuclide) for nuclide in nuclides]
            axis_index = self.num_filters
            mean = np.take(mean, indices=nuclide_bins, axis=axis_index)
            std_dev = np.take(std_dev, indices=nuclide_bins, axis=axis_index)
            mean = np.sum(mean, axis=axis_index, keepdims=True)
            std_dev = np.sum(std_dev**2, axis=axis_index, keepdims=True)
            std_dev = np.sqrt(std_dev)

            # Add AggregateNuclide to the tally sum
            nuclide_sum = openmc.AggregateNuclide(nuclides, 'sum')
            tally_sum.nuclides.append(nuclide_sum)

        # Add a copy of this tally's nuclides to the tally sum
        else:
            tally_sum._nuclides = copy.deepcopy(self.nuclides)

        # Sum across any scores specified by the user
        if len(scores) != 0:
            score_bins = [self.get_score_index(score) for score in scores]
            axis_index = self.num_filters + 1
            mean = np.take(mean, indices=score_bins, axis=axis_index)
            std_dev = np.take(std_dev, indices=score_bins, axis=axis_index)
            mean = np.sum(mean, axis=axis_index, keepdims=True)
            std_dev = np.sum(std_dev**2, axis=axis_index, keepdims=True)
            std_dev = np.sqrt(std_dev)

            # Add AggregateScore to the tally sum
            score_sum = openmc.AggregateScore(scores, 'sum')
            tally_sum.scores.append(score_sum)

        # Add a copy of this tally's scores to the tally sum
        else:
            tally_sum._scores = copy.deepcopy(self.scores)

        # Reshape condensed data arrays with one dimension for all filters
        mean = np.reshape(mean, tally_sum.shape)
        std_dev = np.reshape(std_dev, tally_sum.shape)

        # Assign tally sum's data with the new arrays
        tally_sum._mean = mean
        tally_sum._std_dev = std_dev

        # If original tally was sparse, sparsify the tally summation
        tally_sum.sparse = self.sparse
        return tally_sum

    def average(self, scores=[], filter_type=None,
                filter_bins=[], nuclides=[], remove_filter=False):
        """Vectorized average of tally data across scores, filter bins and/or
        nuclides using tally aggregation.

        This method constructs a new tally to encapsulate the average of the
        data represented by the average of the data in this tally. The tally
        data average is determined by the scores, filter bins and nuclides
        specified in the input parameters.

        Parameters
        ----------
        scores : list of str
            A list of one or more score strings to average across
            (e.g., ['absorption', 'nu-fission']; default is [])
        filter_type : openmc.FilterMeta
            Type of the filter, e.g. MeshFilter
        filter_bins : Iterable of int or tuple
            A list of the filter bins corresponding to the filter_type parameter
            Each bin in the list is the integer ID for 'material', 'surface',
            'cell', 'cellborn', and 'universe' Filters. Each bin is an integer
            for the cell instance ID for 'distribcell' Filters. Each bin is a
            2-tuple of floats for 'energy' and 'energyout' filters corresponding
            to the energy boundaries of the bin of interest. Each bin is an
            (x,y,z) 3-tuple for 'mesh' filters corresponding to the mesh cell of
            interest.
        nuclides : list of str
            A list of nuclide name strings to average across
            (e.g., ['U235', 'U238']; default is [])
        remove_filter : bool
            If a filter is being averaged over, this bool indicates whether to
            remove that filter in the returned tally. Default is False.

        Returns
        -------
        openmc.Tally
            A new tally which encapsulates the average of data requested.
        """

        # Create new derived Tally for average
        tally_avg = Tally()
        tally_avg._derived = True
        tally_avg._estimator = self.estimator
        tally_avg._num_realizations = self.num_realizations
        tally_avg._with_batch_statistics = self.with_batch_statistics
        tally_avg._with_summary = self.with_summary
        tally_avg._sp_filename = self._sp_filename
        tally_avg._results_read = self._results_read

        # Get tally data arrays reshaped with one dimension per filter
        mean = self.get_reshaped_data(value='mean')
        std_dev = self.get_reshaped_data(value='std_dev')

        # Average across any filter bins specified by the user
        if isinstance(filter_type, openmc.FilterMeta):
            find_filter = self.find_filter(filter_type)

            # If user did not specify filter bins, average across all bins
            if len(filter_bins) == 0:
                bin_indices = np.arange(find_filter.num_bins)

                if isinstance(find_filter, openmc.DistribcellFilter):
                    filter_bins = np.arange(find_filter.num_bins)
                elif isinstance(find_filter, openmc.EnergyFunctionFilter):
                    filter_bins = [None]
                else:
                    filter_bins = find_filter.bins

            # Only average across bins specified by the user
            else:
                bin_indices = \
                    [find_filter.get_bin_index(bin) for bin in filter_bins]

            # Average across the bins in the user-specified filter
            for i, self_filter in enumerate(self.filters):
                if isinstance(self_filter, filter_type):
                    shape = mean.shape
                    mean = np.take(mean, indices=bin_indices, axis=i)
                    std_dev = np.take(std_dev, indices=bin_indices, axis=i)

                    # NumPy take introduces a new dimension in output array
                    # for some special cases that must be removed
                    if len(mean.shape) > len(shape):
                        mean = np.squeeze(mean, axis=i)
                        std_dev = np.squeeze(std_dev, axis=i)

                    mean = np.nanmean(mean, axis=i, keepdims=True)
                    std_dev = np.nanmean(std_dev**2, axis=i, keepdims=True)
                    std_dev /= len(bin_indices)
                    std_dev = np.sqrt(std_dev)

                    # Add AggregateFilter to the tally avg
                    if not remove_filter:
                        filter_sum = openmc.AggregateFilter(self_filter,
                            [tuple(filter_bins)], 'avg')
                        tally_avg.filters.append(filter_sum)

                # Add a copy of each filter not averaged across to the tally avg
                else:
                    tally_avg.filters.append(copy.deepcopy(self_filter))

        # Add a copy of this tally's filters to the tally avg
        else:
            tally_avg._filters = copy.deepcopy(self.filters)

        # Sum across any nuclides specified by the user
        if len(nuclides) != 0:
            nuclide_bins = [self.get_nuclide_index(nuclide) for nuclide in nuclides]
            axis_index = self.num_filters
            mean = np.take(mean, indices=nuclide_bins, axis=axis_index)
            std_dev = np.take(std_dev, indices=nuclide_bins, axis=axis_index)
            mean = np.nanmean(mean, axis=axis_index, keepdims=True)
            std_dev = np.nanmean(std_dev**2, axis=axis_index, keepdims=True)
            std_dev /= len(nuclide_bins)
            std_dev = np.sqrt(std_dev)

            # Add AggregateNuclide to the tally avg
            nuclide_avg = openmc.AggregateNuclide(nuclides, 'avg')
            tally_avg.nuclides.append(nuclide_avg)

        # Add a copy of this tally's nuclides to the tally avg
        else:
            tally_avg._nuclides = copy.deepcopy(self.nuclides)

        # Sum across any scores specified by the user
        if len(scores) != 0:
            score_bins = [self.get_score_index(score) for score in scores]
            axis_index = self.num_filters + 1
            mean = np.take(mean, indices=score_bins, axis=axis_index)
            std_dev = np.take(std_dev, indices=score_bins, axis=axis_index)
            mean = np.nanmean(mean, axis=axis_index, keepdims=True)
            std_dev = np.nanmean(std_dev**2, axis=axis_index, keepdims=True)
            std_dev /= len(score_bins)
            std_dev = np.sqrt(std_dev)

            # Add AggregateScore to the tally avg
            score_sum = openmc.AggregateScore(scores, 'avg')
            tally_avg.scores.append(score_sum)

        # Add a copy of this tally's scores to the tally avg
        else:
            tally_avg._scores = copy.deepcopy(self.scores)

        # Reshape condensed data arrays with one dimension for all filters
        mean = np.reshape(mean, tally_avg.shape)
        std_dev = np.reshape(std_dev, tally_avg.shape)

        # Assign tally avg's data with the new arrays
        tally_avg._mean = mean
        tally_avg._std_dev = std_dev

        # If original tally was sparse, sparsify the tally average
        tally_avg.sparse = self.sparse
        return tally_avg

    def diagonalize_filter(self, new_filter, filter_position=-1):
        """Diagonalize the tally data array along a new axis of filter bins.

        This is a helper method for the tally arithmetic methods. This method
        adds the new filter to a derived tally constructed copied from this one.
        The data in the derived tally arrays is "diagonalized" along the bins in
        the new filter. This functionality is used by the openmc.mgxs module; to
        transport-correct scattering matrices by subtracting a 'scatter-P1'
        reaction rate tally with an energy filter from a 'scatter' reaction
        rate tally with both energy and energyout filters.

        Parameters
        ----------
        new_filter : Filter
            The filter along which to diagonalize the data in the new
        filter_position : int
            Where to place the new filter in the Tally.filters list. Defaults
            to last position.

        Returns
        -------
        openmc.Tally
            A new derived Tally with data diagaonalized along the new filter.

        """

        cv.check_type('new_filter', new_filter, _FILTER_CLASSES)
        cv.check_type('filter_position', filter_position, Integral)

        if new_filter in self.filters:
            msg = 'Unable to diagonalize Tally ID="{}" which already ' \
                  'contains a "{}" filter'.format(self.id, type(new_filter))
            raise ValueError(msg)

        # Add the new filter to a copy of this Tally
        new_tally = copy.deepcopy(self)
        new_tally.filters.insert(filter_position, new_filter)

        # Determine "base" indices along the new "diagonal", and the factor
        # by which the "base" indices should be repeated to account for all
        # other filter bins in the diagonalized tally
        indices = np.arange(0, new_filter.num_bins**2, new_filter.num_bins+1)
        diag_factor = self.num_filter_bins // new_filter.num_bins
        diag_indices = np.zeros(self.num_filter_bins, dtype=int)

        # Determine the filter indices along the new "diagonal"
        for i in range(diag_factor):
            start = i * new_filter.num_bins
            end = (i+1) * new_filter.num_bins
            diag_indices[start:end] = indices + (i * new_filter.num_bins**2)

        # Inject this Tally's data along the diagonal of the diagonalized Tally
        if not self.derived and self.sum is not None:
            new_tally._sum = np.zeros(new_tally.shape, dtype=np.float64)
            new_tally._sum[diag_indices, :, :] = self.sum
        if not self.derived and self.sum_sq is not None:
            new_tally._sum_sq = np.zeros(new_tally.shape, dtype=np.float64)
            new_tally._sum_sq[diag_indices, :, :] = self.sum_sq
        if self.mean is not None:
            new_tally._mean = np.zeros(new_tally.shape, dtype=np.float64)
            new_tally._mean[diag_indices, :, :] = self.mean
        if self.std_dev is not None:
            new_tally._std_dev = np.zeros(new_tally.shape, dtype=np.float64)
            new_tally._std_dev[diag_indices, :, :] = self.std_dev

        # If original tally was sparse, sparsify the diagonalized tally
        new_tally.sparse = self.sparse
        return new_tally


class Tallies(cv.CheckedList):
    """Collection of Tallies used for an OpenMC simulation.

    This class corresponds directly to the tallies.xml input file. It can be
    thought of as a normal Python list where each member is a :class:`Tally`. It
    behaves like a list as the following example demonstrates:

    >>> t1 = openmc.Tally()
    >>> t2 = openmc.Tally()
    >>> t3 = openmc.Tally()
    >>> tallies = openmc.Tallies([t1])
    >>> tallies.append(t2)
    >>> tallies += [t3]

    Parameters
    ----------
    tallies : Iterable of openmc.Tally
        Tallies to add to the collection

    """

    def __init__(self, tallies=None):
        super().__init__(Tally, 'tallies collection')
        if tallies is not None:
            self += tallies

    def append(self, tally, merge=False):
        """Append tally to collection

        Parameters
        ----------
        tally : openmc.Tally
            Tally to append
        merge : bool
            Indicate whether the tally should be merged with an existing tally,
            if possible. Defaults to False.

        """
        if not isinstance(tally, Tally):
            msg = 'Unable to add a non-Tally "{}" to the ' \
                  'Tallies instance'.format(tally)
            raise TypeError(msg)

        if merge:
            merged = False

            # Look for a tally to merge with this one
            for i, tally2 in enumerate(self):

                # If a mergeable tally is found
                if tally2.can_merge(tally):
                    # Replace tally2 with the merged tally
                    merged_tally = tally2.merge(tally)
                    self[i] = merged_tally
                    merged = True
                    break

            # If no mergeable tally was found, simply add this tally
            if not merged:
                super().append(tally)

        else:
            super().append(tally)

    def insert(self, index, item):
        """Insert tally before index

        Parameters
        ----------
        index : int
            Index in list
        item : openmc.Tally
            Tally to insert

        """
        super().insert(index, item)

    def merge_tallies(self):
        """Merge any mergeable tallies together. Note that n-way merges are
        possible.

        """

        for i, tally1 in enumerate(self):
            for j, tally2 in enumerate(self):
                # Do not merge the same tally with itself
                if i == j:
                    continue

                # If the two tallies are mergeable
                if tally1.can_merge(tally2):
                    # Replace tally 1 with the merged tally
                    merged_tally = tally1.merge(tally2)
                    self[i] = merged_tally

                    # Remove tally 2 since it is no longer needed
                    self.pop(j)

                    # Continue iterating from the first loop
                    break

    def _create_tally_subelements(self, root_element):
        for tally in self:
            root_element.append(tally.to_xml_element())

    def _create_mesh_subelements(self, root_element):
        already_written = set()
        for tally in self:
            for f in tally.filters:
                if isinstance(f, openmc.MeshFilter):
                    if f.mesh.id not in already_written:
                        if len(f.mesh.name) > 0:
                            root_element.append(ET.Comment(f.mesh.name))

                        root_element.append(f.mesh.to_xml_element())
                        already_written.add(f.mesh.id)

    def _create_filter_subelements(self, root_element):
        already_written = dict()
        for tally in self:
            for f in tally.filters:
                if f not in already_written:
                    root_element.append(f.to_xml_element())
                    already_written[f] = f.id
                elif f.id != already_written[f]:
                    # Set the IDs of identical filters with different
                    # user-defined IDs to the same value
                    f.id = already_written[f]

    def _create_derivative_subelements(self, root_element):
        # Get a list of all derivatives referenced in a tally.
        derivs = []
        for tally in self:
            deriv = tally.derivative
            if deriv is not None and deriv not in derivs:
                derivs.append(deriv)

        # Add the derivatives to the XML tree.
        for d in derivs:
            root_element.append(d.to_xml_element())

    def export_to_xml(self, path='tallies.xml'):
        """Create a tallies.xml file that can be used for a simulation.

        Parameters
        ----------
        path : str
            Path to file to write. Defaults to 'tallies.xml'.

        """

        root_element = ET.Element("tallies")
        self._create_mesh_subelements(root_element)
        self._create_filter_subelements(root_element)
        self._create_tally_subelements(root_element)
        self._create_derivative_subelements(root_element)

        # Clean the indentation in the file to be user-readable
        clean_indentation(root_element)

        # Check if path is a directory
        p = Path(path)
        if p.is_dir():
            p /= 'tallies.xml'

        # Write the XML Tree to the tallies.xml file
        tree = ET.ElementTree(root_element)
        tree.write(str(p), xml_declaration=True, encoding='utf-8')
