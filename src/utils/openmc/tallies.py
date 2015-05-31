import copy
import os
import itertools
from xml.etree import ElementTree as ET

import numpy as np

from openmc import Mesh, Filter, Trigger, Nuclide
from openmc.summary import Summary
from openmc.clean_xml import *
from openmc.checkvalue import *


# "Static" variable for auto-generated Tally IDs
AUTO_TALLY_ID = 10000


def reset_auto_tally_id():
    global AUTO_TALLY_ID
    AUTO_TALLY_ID = 10000


class Tally(object):

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


    def _align_tally_data(self, other):

        self_mean = copy.deepcopy(self.mean)
        self_std_dev = copy.deepcopy(self.std_dev)
        other_mean = copy.deepcopy(other.mean)
        other_std_dev = copy.deepcopy(other.std_dev)

        if self.filters != other.filters:

            # Determine the number of paired combinations of filter bins
            # between the two tallies and repeat arrays along filter axes
            self_num_filter_bins = self.mean.shape[0]
            other_num_filter_bins = other.mean.shape[0]
            num_filter_bins = self_num_filter_bins * other_num_filter_bins
            self_new_bins = max(num_filter_bins - self_num_filter_bins, 1)
            other_new_bins = max(num_filter_bins - other_num_filter_bins, 1)

            # Replicate the data
            self_mean = np.repeat(self_mean, self_new_bins, axis=0)
            other_mean = np.tile(other_mean, (other_new_bins, 1, 1))
            self_std_dev = np.repeat(self_std_dev, self_new_bins, axis=0)
            other_std_dev = np.tile(other_std_dev, (other_new_bins, 1, 1))

        if self.nuclides != other.nuclides:

            # Determine the number of paired combinations of nuclides
            # between the two tallies and repeat arrays along nuclide axes
            self_num_nuclide_bins = self.mean.shape[1]
            other_num_nuclide_bins = other.mean.shape[1]
            num_nuclide_bins = self_num_nuclide_bins * other_num_nuclide_bins
            self_new_bins = max(num_nuclide_bins - self_num_nuclide_bins, 1)
            other_new_bins = max(num_nuclide_bins - other_num_nuclide_bins, 1)

            # Replicate the data
            self_mean = np.repeat(self_mean, self_new_bins, axis=1)
            other_mean = np.tile(other_mean, (1, other_new_bins, 1))
            self_std_dev = np.repeat(self_std_dev, self_new_bins, axis=1)
            other_std_dev = np.tile(other_std_dev, (1, other_new_bins, 1))

        if self.scores != other.scores:

            # Determine the number of paired combinations of score bins
            # between the two tallies and repeat arrays along score axes
            self_num_score_bins = self.mean.shape[2]
            other_num_score_bins = other.mean.shape[2]
            num_score_bins = self_num_score_bins * other_num_score_bins
            self_new_bins = max(num_score_bins - self_num_score_bins, 1)
            other_new_bins = max(num_score_bins - other_num_score_bins, 1)

            # Replicate the data
            self_mean = np.repeat(self_mean, self_new_bins, axis=2)
            other_mean = np.tile(other_mean, (1, 1, other_new_bins))
            self_std_dev = np.repeat(self_std_dev, self_new_bins, axis=2)
            other_std_dev = np.tile(other_std_dev, (1, 1, other_new_bins))

        data = {}
        data['self'] = {}
        data['other'] = {}
        data['self']['mean'] = self_mean
        data['other']['mean'] = other_mean
        data['self']['std. dev.'] = self_std_dev
        data['other']['std. dev.'] = other_std_dev
        return data


    def __add__(self, other):

        # Check that results have been read
        if self.mean is None:
            msg = 'Unable to use tally arithmetic with Tally ID={0} ' \
                  'since it does not contain any results.'.format(self.id)
            raise ValueError(msg)

        new_tally = Tally(name='derived')
        new_tally.with_batch_statistics = True

        if isinstance(other, Tally):

            # Check that results have been read
            if other.mean is None:
                msg = 'Unable to use tally arithmetic with Tally ID={0} ' \
                      'since it does not contain any results.'.format(other.id)
                raise ValueError(msg)

            # FIXME: Need new CrossFilter class

            data = self._align_tally_data(other)

            new_tally._mean = data['self']['mean'] + data['other']['mean']
            new_tally._std_dev = np.sqrt(data['self']['std. dev.']**2 + \
                                         data['other']['std. dev.']**2)

            if self.estimator == other.estimator:
                new_tally.estimator = self.estimator

            if self.with_summary and other.with_summary:
                new_tally.with_summary = self.with_summary

            if self.num_realizations == other.num_realizations:
                new_tally.num_realizations = self.num_realizations

            # If the two Tallies have same filters, replicate them in new Tally
            if self.filters == other.filters:
                for filter in self.filters:
                    new_tally.add_filter(filter)

            # Generate filter "cross product"
            else:
                for self_filter in self.filters:
                    for other_filter in other.filters:
                        new_filter = _CrossFilter(self_filter, other_filter, '+')
                        new_tally.add_filter(new_filter)

            # If the two Tallies have same nuclides, replicate them in new Tally
            if self.nuclides == other.nuclides:
                for nuclide in self.nuclides:
                    new_tally.add_nuclide(nuclide)

            # Generate nuclide "cross product"
            else:
                for self_nuclide in self.nuclides:
                    for other_nuclide in other.nuclides:
                        new_nuclide = _CrossNuclide(self_nuclide, other_nuclide, '+')
                        new_tally.add_nuclide(new_nuclide)

            # If the two Tallies have same scores, replicate them in new Tally
            if self.scores == other.scores:
                for score in self.scores:
                    new_tally.add_score(score)

            # Generate score "cross product"
            else:
                for self_score in self.scores:
                    for other_score in other.scores:
                        new_score = _CrossScore(self_score, other_score, '+')
                        new_tally.add_score(new_score)

        elif is_integer(other) or is_float(other):

            new_tally._mean = self._mean + other
            new_tally._std_dev = self._std_dev

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
            msg = 'Unable to add {0} to Tally ID={1}'.format(other, self.id)
            raise ValueError(msg)

        return new_tally


    def __radd__(self, other):
        return self + other


    def __sub__(self, other):

        # Check that results have been read
        if self.mean is None:
            msg = 'Unable to use tally arithmetic with Tally ID={0} ' \
                  'since it does not contain any results.'.format(self.id)
            raise ValueError(msg)

        new_tally = Tally(name='derived')
        new_tally.with_batch_statistics = True

        if isinstance(other, Tally):

            # Check that results have been read
            if other.mean is None:
                msg = 'Unable to use tally arithmetic with Tally ID={0} ' \
                      'since it does not contain any results.'.format(other.id)
                raise ValueError(msg)

            # FIXME: Need new CrossFilter class

            data = self._align_tally_data(other)

            new_tally._mean = data['self']['mean'] - data['other']['mean']
            new_tally._std_dev = np.sqrt(data['self']['std. dev.']**2 + \
                                         data['other']['std. dev.']**2)

            if self.estimator == other.estimator:
                new_tally.estimator = self.estimator

            if self.with_summary and other.with_summary:
                new_tally.with_summary = self.with_summary

            if self.num_realizations == other.num_realizations:
                new_tally.num_realizations = self.num_realizations

            # If the two Tallies have same filters, replicate them in new Tally
            if self.filters == other.filters:
                for filter in self.filters:
                    new_tally.add_filter(filter)

            # Generate filter "cross product"
            else:
                for self_filter in self.filters:
                    for other_filter in other.filters:
                        new_filter = _CrossFilter(self_filter, other_filter, '-')
                        new_tally.add_filter(new_filter)

            # If the two Tallies have same nuclides, replicate them in new Tally
            if self.nuclides == other.nuclides:
                for nuclide in self.nuclides:
                    new_tally.add_nuclide(nuclide)

            # Generate nuclide "cross product"
            else:
                for self_nuclide in self.nuclides:
                    for other_nuclide in other.nuclides:
                        new_nuclide = _CrossNuclide(self_nuclide, other_nuclide, '-')
                        new_tally.add_nuclide(new_nuclide)

            # If the two Tallies have same scores, replicate them in new Tally
            if self.scores == other.scores:
                for score in self.scores:
                    new_tally.add_score(score)

            # Generate score "cross product"
            else:
                for self_score in self.scores:
                    for other_score in other.scores:
                        new_score = _CrossScore(self_score, other_score, '-')
                        new_tally.add_score(new_score)

        elif is_integer(other) or is_float(other):

            new_tally._mean = self._mean - other
            new_tally._std_dev = self._std_dev

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
            msg = 'Unable to subtract {0} from Tally ' \
                  'ID={1}'.format(other, self.id)
            raise ValueError(msg)

        return new_tally


    def __rsub__(self, other):
        return -1. * self + other


    def __mul__(self, other):

        # Check that results have been read
        if self.mean is None:
            msg = 'Unable to use tally arithmetic with Tally ID={0} ' \
                  'since it does not contain any results.'.format(self.id)
            raise ValueError(msg)

        new_tally = Tally(name='derived')
        new_tally.with_batch_statistics = True

        if isinstance(other, Tally):

            # Check that results have been read
            if other.mean is None:
                msg = 'Unable to use tally arithmetic with Tally ID={0} ' \
                      'since it does not contain any results.'.format(other.id)
                raise ValueError(msg)

            data = self._align_tally_data(other)

            self_rel_err = data['self']['std. dev.'] / data['self']['mean']
            other_rel_err = data['other']['std. dev.'] / data['other']['mean']
            new_tally._mean = data['self']['mean'] * data['other']['mean']
            new_tally._std_dev = np.abs(new_tally.mean) * \
                                 np.sqrt(self_rel_err**2 + other_rel_err**2)

            if self.estimator == other.estimator:
                new_tally.estimator = self.estimator

            if self.with_summary and other.with_summary:
                new_tally.with_summary = self.with_summary

            if self.num_realizations == other.num_realizations:
                new_tally.num_realizations = self.num_realizations

            # If the two Tallies have same filters, replicate them in new Tally
            if self.filters == other.filters:
                for filter in self.filters:
                    new_tally.add_filter(filter)

            # Generate filter "cross product"
            else:
                for self_filter in self.filters:
                    for other_filter in other.filters:
                        new_filter = _CrossFilter(self_filter, other_filter, '*')
                        new_tally.add_filter(new_filter)

            # If the two Tallies have same nuclides, replicate them in new Tally
            if self.nuclides == other.nuclides:
                for nuclide in self.nuclides:
                    new_tally.add_nuclide(nuclide)

            # Generate nuclide "cross product"
            else:
                for self_nuclide in self.nuclides:
                    for other_nuclide in other.nuclides:
                        new_nuclide = _CrossNuclide(self_nuclide, other_nuclide, '*')
                        new_tally.add_nuclide(new_nuclide)

            # If the two Tallies have same scores, replicate them in new Tally
            if self.scores == other.scores:
                for score in self.scores:
                    new_tally.add_score(score)

            # Generate score "cross product"
            else:
                for self_score in self.scores:
                    for other_score in other.scores:
                        new_score = _CrossScore(self_score, other_score, '*')
                        new_tally.add_score(new_score)

        elif is_integer(other) or is_float(other):

            new_tally._mean = self._mean * other
            new_tally._std_dev = self._std_dev * np.abs(other)

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
            msg = 'Unable to multiply Tally ID={1} ' \
                  'by {0}'.format(self.id, other)
            raise ValueError(msg)

        return new_tally


    def __rmul__(self, other):
        return self * other


    def __div__(self, other):

        # Check that results have been read
        if self.mean is None:
            msg = 'Unable to use tally arithmetic with Tally ID={0} ' \
                  'since it does not contain any results.'.format(self.id)
            raise ValueError(msg)

        new_tally = Tally(name='derived')
        new_tally.with_batch_statistics = True

        if isinstance(other, Tally):

            # Check that results have been read
            if other.mean is None:
                msg = 'Unable to use tally arithmetic with Tally ID={0} ' \
                      'since it does not contain any results.'.format(other.id)
                raise ValueError(msg)

            # FIXME: Need new CrossFilter class

            data = self._align_tally_data(other)

            self_rel_err = data['self']['std. dev.'] / data['self']['mean']
            other_rel_err = data['other']['std. dev.'] / data['other']['mean']
            new_tally._mean = data['self']['mean'] / data['other']['mean']
            new_tally._std_dev = np.abs(new_tally.mean) * \
                                 np.sqrt(self_rel_err**2 + other_rel_err**2)

            if self.estimator == other.estimator:
                new_tally.estimator = self.estimator

            if self.with_summary and other.with_summary:
                new_tally.with_summary = self.with_summary

            if self.num_realizations == other.num_realizations:
                new_tally.num_realizations = self.num_realizations

            # If the two Tallies have same filters, replicate them in new Tally
            if self.filters == other.filters:
                for filter in self.filters:
                    new_tally.add_filter(filter)

            # Generate filter "cross product"
            else:
                for self_filter in self.filters:
                    for other_filter in other.filters:
                        new_filter = _CrossFilter(self_filter, other_filter, '/')
                        new_tally.add_filter(new_filter)

            # If the two Tallies have same nuclides, replicate them in new Tally
            if self.nuclides == other.nuclides:
                for nuclide in self.nuclides:
                    new_tally.add_nuclide(nuclide)

            # Generate nuclide "cross product"
            else:
                for self_nuclide in self.nuclides:
                    for other_nuclide in other.nuclides:
                        new_nuclide = _CrossNuclide(self_nuclide, other_nuclide, '/')
                        new_tally.add_nuclide(new_nuclide)

            # If the two Tallies have same scores, replicate them in new Tally
            if self.scores == other.scores:
                for score in self.scores:
                    new_tally.add_score(score)

            # Generate score "cross product"
            else:
                for self_score in self.scores:
                    for other_score in other.scores:
                        new_score = _CrossScore(self_score, other_score, '/')
                        new_tally.add_score(new_score)

        elif is_integer(other) or is_float(other):

            new_tally._mean = self._mean / other
            new_tally._std_dev = self._std_dev * np.abs(1. / other)

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
            msg = 'Unable to divide Tally ID={0} ' \
                  'by {1}'.format(self.id, other)
            raise ValueError(msg)

        return new_tally


    def __rdiv__(self, other):
        return self * (1. / other)


    def __pow__(self, power):

        # Check that results have been read
        if self.mean is None:
            msg = 'Unable to use tally arithmetic with Tally ID={0} ' \
                  'since it does not contain any results.'.format(self.id)
            raise ValueError(msg)

        new_tally = Tally(name='derived')
        new_tally.with_batch_statistics = True

        if isinstance(power, Tally):

            # Check that results have been read
            if power.mean is None:
                msg = 'Unable to use tally arithmetic with Tally ID={0} ' \
                      'since it does not contain any results.'.format(power.id)
                raise ValueError(msg)

            # FIXME: Need new CrossFilter class

            data = self._align_tally_data(power)

            mean_ratio = data['other']['mean'] / data['self']['mean']
            first_term = mean_ratio * data['self']['std. dev.']
            second_term = np.log(data['self']['mean']) * data['other']['std. dev.']
            new_tally._mean = data['self']['mean'] ** data['other']['mean']
            new_tally._std_dev = np.abs(new_tally.mean) * \
                                 np.sqrt(first_term**2 + second_term**2)

            if self.estimator == power.estimator:
                new_tally.estimator = self.estimator

            if self.with_summary and power.with_summary:
                new_tally.with_summary = self.with_summary

            if self.num_realizations == power.num_realizations:
                new_tally.num_realizations = self.num_realizations

            # If the two Tallies have same filters, replicate them in new Tally
            if self.filters == power.filters:
                for filter in self.filters:
                    new_tally.add_filter(filter)

            # Generate filter "cross product"
            else:
                for self_filter in self.filters:
                    for power_filter in power.filters:
                        new_filter = _CrossFilter(self_filter, power_filter, '^')
                        new_tally.add_filter(new_filter)

            # If the two Tallies have same nuclides, replicate them in new Tally
            if self.nuclides == power.nuclides:
                for nuclide in self.nuclides:
                    new_tally.add_nuclide(nuclide)

            # Generate nuclide "cross product"
            else:
                for self_nuclide in self.nuclides:
                    for power_nuclide in power.nuclides:
                        new_nuclide = _CrossNuclide(self_nuclide, power_nuclide, '^')
                        new_tally.add_nuclide(new_nuclide)

            # If the two Tallies have same scores, replicate them in new Tally
            if self.scores == power.scores:
                for score in self.scores:
                    new_tally.add_score(score)

            # Generate score "cross product"
            else:
                for self_score in self.scores:
                    for power_score in power.scores:
                        new_score = _CrossScore(self_score, power_score, '^')
                        new_tally.add_score(new_score)

        elif is_integer(power) or is_float(power):

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
            msg = 'Unable to raise Tally ID={0} to ' \
                  'power {1}'.format(self.id, power)
            raise ValueError(msg)

        return new_tally

    '''
    def sum(self, axis=None):

    def slice(self, filters=[], nuclides=[], scores=[])
    '''


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
            if not filter in tally2.filters:
                return False

        # Check all nuclides
        if len(self.nuclides) != len(tally2.nuclides):
            return False

        for nuclide in self.nuclides:
            if not nuclide in tally2.nuclides:
                return False

        # Check all scores
        if len(self.scores) != len(tally2.scores):
            return False

        for score in self.scores:
            if not score in tally2.scores:
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
    def num_filters(self):
        return len(self._filters)


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


    @property
    def with_batch_statistics(self):
        return self._with_batch_statistics


    @estimator.setter
    def estimator(self, estimator):

        if not estimator in ['analog', 'tracklength']:
            msg = 'Unable to set the estimator for Tally ID={0} to {1} since ' \
                  'it is not a valid estimator type'.format(self.id, estimator)
            raise ValueError(msg)

        self._estimator = estimator


    def add_trigger(self, trigger):

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

        # Check that the ID is an integer and wasn't already used
        elif not is_integer(tally_id):
            msg = 'Unable to set a non-integer Tally ID {0}'.format(tally_id)
            raise ValueError(msg)

        elif tally_id < 0:
            msg = 'Unable to set Tally ID to {0} since it must be a ' \
                  'non-negative integer'.format(tally_id)
            raise ValueError(msg)

        else:
            self._id = tally_id


    @name.setter
    def name(self, name):

        if not is_string(name):
            msg = 'Unable to set name for Tally ID={0} with a non-string ' \
                  'value "{1}"'.format(self.id, name)
            raise ValueError(msg)

        else:
            self._name = name


    def add_filter(self, filter):

        global filters

        if not isinstance(filter, (Filter, _CrossFilter)):
            msg = 'Unable to add Filter "{0}" to Tally ID={1} since it is ' \
                  'not a Filter object'.format(filter, self.id)
            raise ValueError(msg)

        self._filters.append(filter)


    def add_nuclide(self, nuclide):
        self._nuclides.append(nuclide)


    def add_score(self, score):

        if not is_string(score) and not isinstance(score, _CrossScore):
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

        if not is_integer(num_realizations):
            msg = 'Unable to set the number of realizations to "{0}" for ' \
                  'Tally ID={1} since it is not an ' \
                  'integer'.format(num_realizations)
            raise ValueError(msg)

        elif num_realizations < 0:
            msg = 'Unable to set the number of realizations to "{0}" for ' \
                  'Tally ID={1} since it is a negative ' \
                  'value'.format(num_realizations)
            raise ValueError(msg)

        self._num_realizations = num_realizations


    @with_summary.setter
    def with_summary(self, with_summary):

        if not isinstance(with_summary, bool):
            msg = 'Unable to set with_summary to a non-boolean ' \
                  'value "{0}"'.format(with_summary)
            raise ValueError(msg)

        self._with_summary = with_summary


    @with_batch_statistics.setter
    def with_batch_statistics(self, with_batch_statistics):

        if not isinstance(with_batch_statistics, bool):
            msg = 'Unable to set with_batch_statistics to a non-boolean ' \
                  'value "{0}"'.format(with_batch_statistics)
            raise ValueError(msg)

        self._with_batch_statistics = with_batch_statistics


    def set_results(self, sum, sum_sq):

        if not isinstance(sum, (tuple, list, np.ndarray)):
            msg = 'Unable to set the sum to "{0}" for Tally ID={1} since ' \
                  'it is not a Python tuple/list or NumPy ' \
                  'array'.format(sum, self.id)
            raise ValueError(msg)

        if not isinstance(sum_sq, (tuple, list, np.ndarray)):
            msg = 'Unable to set the sum to "{0}" for Tally ID={1} since ' \
                  'it is not a Python tuple/list or NumPy ' \
                  'array'.format(sum_sq, self.id)
            raise ValueError(msg)

        self._sum = sum
        self._sum_sq = sum_sq


    def remove_score(self, score):

        if not score in self.scores:
            msg = 'Unable to remove score "{0}" from Tally ID={1} since the ' \
                  'Tally does not contain this score'.format(score, self.id)
            ValueError(msg)

        self._scores.remove(score)


    def remove_filter(self, filter):

        if not filter in self.filters:
            msg = 'Unable to remove filter "{0}" from Tally ID={1} since the ' \
                  'Tally does not contain this filter'.format(filter, self.id)
            ValueError(msg)

        self._filters.remove(filter)


    def remove_nuclide(self, nuclide):

        if not nuclide in self.nuclides:
            msg = 'Unable to remove nuclide "{0}" from Tally ID={1} since the ' \
                  'Tally does not contain this nuclide'.format(nuclide, self.id)
            ValueError(msg)

        self._nuclides.remove(nuclide)


    def compute_std_dev(self, t_value=1.0):

        # Calculate sample mean and standard deviation
        self._mean = self.sum / self.num_realizations
        self._std_dev = np.sqrt((self.sum_sq / self.num_realizations - \
                                 self.mean**2) / (self.num_realizations - 1))
        self._std_dev *= t_value
        self.with_batch_statistics = True


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

        if not isinstance(tally, Tally):
            return False

        # Must have same estimator
        if self.estimator != tally.estimator:
            return False

        # Must have same nuclides
        if len(self.nuclides) != len(tally.nuclides):
            return False

        for nuclide in self.nuclides:
            if not nuclide in tally.nuclides:
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

            # If no mergeable filter was found, the tallies are not mergable
            if not mergeable_filter:
                return False

        # Tallies are mergeable if all conditional checks passed
        return True


    def merge(self, tally):

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

            if not filter.bins is None:

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
        if not self.estimator is None:
            subelement = ET.SubElement(element, "estimator")
            subelement.text = self.estimator

        # Optional Triggers
        for trigger in self.triggers:
            trigger.get_trigger_xml(element)

        return element


    def find_filter(self, filter_type):

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
              'cellborn', and 'universe' Filters. The bin is an integer for
              the cell instance ID for 'distribcell' Filters. The bin is
              a 2-tuple of floats for 'energy' and 'energyout' filters
              corresponding to the energy boundaries of the bin of interest.
              The bin is a (x,y,z) 3-tuple for 'mesh' filters corresponding to
              the mesh cell of interest.

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
             The index in the Tally data array for this nuclide.

        Raises
        ------
             KeyError : An error when the argument passed to the 'nuclide' 
             parameter cannot be found in the Tally.
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
             The index in the Tally data array for this score.

        Raises
        ------
             ValueError: An error when the argument passed to the 'score' 
             parameter cannot be found in the Tally.
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
             parameter (e.g., [1, (0., 0.625e-6)]; default is []). Each bin
             in the list is the integer ID for 'material', 'surface', 'cell',
             'cellborn', and 'universe' Filters. Each bin is an integer for
             the cell instance ID for 'distribcell Filters. Each bin is
             a 2-tuple of floats for 'energy' and 'energyout' filters
             corresponding to the energy boundaries of the bin of interest.
             The bin is a (x,y,z) 3-tuple for 'mesh' filters corresponding
             to the mesh cell of interest. The order of the bins in the list
             must correspond of the filter_types parameter.

        nuclides : list
             A list of nuclide name strings
             (e.g., ['U-235', 'U-238']; default is [])

        value : str
             A string for the type of value to return  - 'mean' (default),
             'std_dev', 'rel_err', 'sum', or 'sum_sq' are accepted

        Returns
        -------
             A scalar or NumPy array of the Tally data indexed in the order
             each filter, nuclide and score is listed in the parameters.

        Raises
        ------
             ValueError : An error when this routine is called before the Tally
             is populated with data by the StatePoint.read_results() routine.
        """

        # Ensure that StatePoint.read_results() was called first
        if (value == 'mean' and self.mean is None) or \
           (value == 'std_dev' and self.std_dev is None) or \
           (value == 'rel_err' and self.mean is None) or \
           (value == 'sum' and self.sum is None) or \
           (value == 'sum_sq' and self.sum_sq is None):
            msg = 'The Tally ID={0} has no data to return. Call the ' \
                  'StatePoint.read_results() routine before using ' \
                  'Tally.get_values(...)'.format(self.id)
            raise ValueError(msg)


        # Compute batch statistics if not yet computed
        if not self.with_batch_statistics:
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
                        xyz = map(lambda x: np.arange(1,x+1), dimension)
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
            filter_indices = map(sum, itertools.product(*filter_indices))

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
             for distribcell tally filters (default is None). The geometric
             information in the Summary object is embedded into a Multi-index
             column with a geometric "path" to each distribcell intance.
             NOTE: This option requires the OpenCG Python package.

        Returns
        -------
             A Pandas DataFrame with each column annotated by filter, nuclide
             and score bin information (if these parameters are True), and the
             mean and standard deviation of the Tally's data.

        Raises
        ------
             KeyError : An error when this routine is called before the Tally
             is populated with data by the StatePoint.read_results() routine.
        """

        # Ensure that StatePoint.read_results() was called first
        if self.mean is None or self.std_dev is None:
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
        if not self.with_batch_statistics:
            self.compute_std_dev()

        # Initialize a pandas dataframe for the tally data
        df = pd.DataFrame()

        # Find the total length of the tally data array
        data_size = self.mean.size

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
                                if coords == None:
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
                                if coords._next == None:
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
             KeyError : An error when this routine is called before the Tally
             is populated with data by the StatePoint.read_results() routine.
        """

        # Ensure that StatePoint.read_results() was called first
        if self._sum is None or self._sum_sq is None:
            msg = 'The Tally ID={0} has no data to export. Call the ' \
                  'StatePoint.read_results() routine before using ' \
                  'Tally.export_results(...)'.format(self.id)
            raise KeyError(msg)

        if not is_string(filename):
            msg = 'Unable to export the results for Tally ID={0} to ' \
                  'filename="{1}" since it is not a ' \
                  'string'.format(self.id, filename)
            raise ValueError(msg)

        elif not is_string(directory):
            msg = 'Unable to export the results for Tally ID={0} to ' \
                  'directory="{1}" since it is not a ' \
                  'string'.format(self.id, directory)
            raise ValueError(msg)

        elif not format in ['hdf5', 'pkl', 'csv']:
            msg = 'Unable to export the results for Tally ID={0} to format ' \
                  '"{1}" since it is not supported'.format(self.id, format)
            raise ValueError(msg)

        elif not isinstance(append, (bool, np.bool)):
            msg = 'Unable to export the results for Tally ID={0} since the ' \
                  'append parameters is not True/False'.format(self.id, append)
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

            import pickle

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

            tally_group['nuclides']= np.array(nuclides)

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

    def __init__(self):

        # Initialize TalliesFile class attributes
        self._tallies = []
        self._meshes = []
        self._tallies_file = ET.Element("tallies")


    def add_tally(self, tally, merge=False):

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
        self._tallies.remove(tally)


    def merge_tallies(self):

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

        if not isinstance(mesh, Mesh):
            msg = 'Unable to add a non-Mesh {0} to the TalliesFile'.format(mesh)
            raise ValueError(msg)

        self._meshes.append(mesh)


    def remove_mesh(self, mesh):
        self._meshes.remove(mesh)


    def create_tally_subelements(self):

        for tally in self._tallies:
            xml_element = tally.get_tally_xml()
            self._tallies_file.append(xml_element)


    def create_mesh_subelements(self):

        for mesh in self._meshes:

            if len(mesh._name) > 0:
                self._tallies_file.append(ET.Comment(mesh._name))

            xml_element = mesh.get_mesh_xml()
            self._tallies_file.append(xml_element)


    def export_to_xml(self):

        self.create_mesh_subelements()
        self.create_tally_subelements()

        # Clean the indentation in the file to be user-readable
        clean_xml_indentation(self._tallies_file)

        # Write the XML Tree to the tallies.xml file
        tree = ET.ElementTree(self._tallies_file)
        tree.write("tallies.xml", xml_declaration=True,
                             encoding='utf-8', method="xml")


class _CrossScore(object):

    def __init__(self, left_score=None, right_score=None, binary_op=None):

        self._left_score = None
        self._right_score = None
        self._binary_op = None

        if left_score is not None:
            self.left_score = left_score

        if right_score is not None:
            self.right_score = right_score

        if binary_op is not None:
            self.binary_op = binary_op


    @property
    def left_score(self):
        return self._left_score


    @property
    def right_score(self):
        return self._right_score


    @property
    def binary_op(self):
        return self._binary_op


    @left_score.setter
    def left_score(self, left_score):

        if not is_string(left_score):
            msg = 'Unable to set CrossScore left score to {0} which ' \
                  'is not a string'.format(left_score)
            raise ValueError(msg)

        self._left_score = left_score


    @right_score.setter
    def right_score(self, right_score):

        if not is_string(right_score):
            msg = 'Unable to set CrossScore right score to {0} which ' \
                  'is not a string'.format(right_score)
            raise ValueError(msg)

        self._right_score = right_score


    @binary_op.setter
    def binary_op(self, binary_op):

        if not is_string(binary_op):
            msg = 'Unable to set CrossScore binary op to {0} which ' \
                  'is not a string'.format(binary_op)
            raise ValueError(msg)

        self._binary_op = binary_op


    def __repr__(self):
        string = '({0} {1} {2})'.format(self.left_score,
                                        self.binary_op, self.right_score)
        return string


class _CrossNuclide(object):

    def __init__(self, left_nuclide=None, right_nuclide=None, binary_op=None):

        self._left_nuclide = None
        self._right_nuclide = None
        self._binary_op = None

        if left_nuclide is not None:
            self.left_nuclide = left_nuclide

        if right_nuclide is not None:
            self.right_nuclide = right_nuclide

        if binary_op is not None:
            self.binary_op = binary_op


    @property
    def left_nuclide(self):
        return self._left_nuclide


    @property
    def right_nuclide(self):
        return self._right_nuclide


    @property
    def binary_op(self):
        return self._binary_op


    @left_nuclide.setter
    def left_nuclide(self, left_nuclide):

        if not isinstance(left_nuclide, Nuclide) or is_integer(left_nuclide):
            msg = 'Unable to set CrossNuclide left nuclide to {0} which ' \
                  'is not an integer or Nuclide'.format(left_nuclide)
            raise ValueError(msg)

        self._left_nuclide = left_nuclide


    @right_nuclide.setter
    def right_nuclide(self, right_nuclide):

        if not isinstance(right_nuclide, Nuclide) or is_integer(right_nuclide):
            msg = 'Unable to set CrossNuclide right nuclide to {0} which ' \
                  'is not an integer or Nuclide'.format(right_nuclide)
            raise ValueError(msg)

        self._right_nuclide = right_nuclide


    @binary_op.setter
    def binary_op(self, binary_op):

        if not is_string(binary_op):
            msg = 'Unable to set CrossNuclide binary op to {0} which ' \
                  'is not a string'.format(binary_op)
            raise ValueError(msg)

        self._binary_op = binary_op


    def __repr__(self):

        string = ''

        # If the Summary was linked, the left nuclide is a Nuclide object
        if isinstance(self.left_nuclide, Nuclide):
            string += '(' + self.left_nuclide.name
        # If the Summary was not linked, the left nuclide is the ZAID
        else:
            string += '(' + str(self.left_nuclide)

        string += ' ' + self.binary_op + ' '

        # If the Summary was linked, the right nuclide is a Nuclide object
        if isinstance(self.right_nuclide, Nuclide):
            string += self.right_nuclide.name + ')'
        # If the Summary was not linked, the right nuclide is the ZAID
        else:
            string += str(self.right_nuclide) + ')'

        return string


class _CrossFilter(object):

    def __init__(self, left_filter=None, right_filter=None, binary_op=None):

        self._left_filter = None
        self._right_filter = None
        self._binary_op = None

        if left_filter is not None:
            self.left_filter = left_filter

        if right_filter is not None:
            self.right_filter = right_filter

        if binary_op is not None:
            self.binary_op = binary_op


    @property
    def left_filter(self):
        return self._left_filter


    @property
    def right_filter(self):
        return self._right_filter


    @property
    def binary_op(self):
        return self._binary_op


    @property
    def type(self):
        return (self.right_filter.type, self.left_filter.type)


    @property
    def bins(self):
        return (self.right_filter.bins, self.left_filter.bins)


    @property
    def stride(self):
        return self.left_filter.stride * self.right_filter.stride


    @left_filter.setter
    def left_filter(self, left_filter):

        if not isinstance(left_filter, Filter):
            msg = 'Unable to set CrossFilter left filter to {0} which ' \
                  'is not a Filter'.format(left_filter)
            raise ValueError(msg)

        self._left_filter = left_filter


    @right_filter.setter
    def right_filter(self, right_filter):

        if not isinstance(right_filter, Filter):
            msg = 'Unable to set CrossFilter right filter to {0} which ' \
                  'is not a Filter'.format(right_filter)
            raise ValueError(msg)

        self._right_filter = right_filter


    @binary_op.setter
    def binary_op(self, binary_op):

        if not is_string(binary_op):
            msg = 'Unable to set CrossFilter binary op to {0} which ' \
                  'is not a string'.format(binary_op)
            raise ValueError(msg)

        self._binary_op = binary_op


    def __repr__(self):

        string = '_CrossFilter\n'

        filter_type = '({0} {1} {2})'.format(self.left_filter.type,
                                             self.binary_op,
                                             self.right_filter.type)

        filter_bins = '({0} {1} {2})'.format(self.left_filter.bins,
                                             self.binary_op,
                                             self.right_filter.bins)

        string += '{0: <16}{1}{2}\n'.format('\tType', '=\t', filter_type)
        string += '{0: <16}{1}{2}\n'.format('\tBins', '=\t', filter_bins)
        return string
