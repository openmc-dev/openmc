import sys
from numbers import Integral

import numpy as np

from openmc import Filter, Nuclide
from openmc.cross import CrossScore, CrossNuclide, CrossFilter
from openmc.filter import _FILTER_TYPES
import openmc.checkvalue as cv


if sys.version_info[0] >= 3:
    basestring = str

# Acceptable tally aggregation operations
_TALLY_AGGREGATE_OPS = ['sum', 'avg']


class AggregateScore(object):
    """A special-purpose tally score used to encapsulate an aggregate of a
    subset or all of tally's scores for tally aggregation.

    Parameters
    ----------
    scores : Iterable of str or CrossScore
        The scores included in the aggregation
    aggregate_op : str
        The tally aggregation operator (e.g., 'sum', 'avg', etc.) used
        to aggregate across a tally's scores with this AggregateScore

    Attributes
    ----------
    scores : Iterable of str or CrossScore
        The scores included in the aggregation
    aggregate_op : str
        The tally aggregation operator (e.g., 'sum', 'avg', etc.) used
        to aggregate across a tally's scores with this AggregateScore

    """

    def __init__(self, scores=None, aggregate_op=None):

        self._scores = None
        self._aggregate_op = None

        if scores is not None:
            self.scores = scores
        if aggregate_op is not None:
            self.aggregate_op = aggregate_op

    def __hash__(self):
        return hash(repr(self))

    def __eq__(self, other):
        return str(other) == str(self)

    def __ne__(self, other):
        return not self == other

    def __deepcopy__(self, memo):
        existing = memo.get(id(self))

        # If this is the first time we have tried to copy this object, create a copy
        if existing is None:
            clone = type(self).__new__(type(self))
            clone._scores = self.scores
            clone._aggregate_op = self.aggregate_op

            memo[id(self)] = clone

            return clone

        # If this object has been copied before, return the first copy made
        else:
            return existing

    def __repr__(self):
        string = ', '.join(map(str, self.scores))
        string = '{0}({1})'.format(self.aggregate_op, string)
        return string

    @property
    def scores(self):
        return self._scores

    @property
    def aggregate_op(self):
        return self._aggregate_op

    @scores.setter
    def scores(self, scores):
        cv.check_iterable_type('scores', scores, basestring)
        self._scores = scores

    @aggregate_op.setter
    def aggregate_op(self, aggregate_op):
        cv.check_type('aggregate_op', aggregate_op, (basestring, CrossScore))
        cv.check_value('aggregate_op', aggregate_op, _TALLY_AGGREGATE_OPS)
        self._aggregate_op = aggregate_op


class AggregateNuclide(object):
    """A special-purpose tally nuclide used to encapsulate an aggregate of a
    subset or all of tally's nuclides for tally aggregation.

    Parameters
    ----------
    nuclides : Iterable of str or Nuclide or CrossNuclide
        The nuclides included in the aggregation
    aggregate_op : str
        The tally aggregation operator (e.g., 'sum', 'avg', etc.) used
        to aggregate across a tally's nuclides with this AggregateNuclide

    Attributes
    ----------
    nuclides : Iterable of str or Nuclide or CrossNuclide
        The nuclides included in the aggregation
    aggregate_op : str
        The tally aggregation operator (e.g., 'sum', 'avg', etc.) used
        to aggregate across a tally's nuclides with this AggregateNuclide

    """

    def __init__(self, nuclides=None, aggregate_op=None):

        self._nuclides = None
        self._aggregate_op = None

        if nuclides is not None:
            self.nuclides = nuclides
        if aggregate_op is not None:
            self.aggregate_op = aggregate_op

    def __hash__(self):
        return hash(repr(self))

    def __eq__(self, other):
        return str(other) == str(self)

    def __ne__(self, other):
        return not self == other

    def __deepcopy__(self, memo):
        existing = memo.get(id(self))

        # If this is the first time we have tried to copy this object, create a copy
        if existing is None:
            clone = type(self).__new__(type(self))
            clone._nuclides = self.nuclides
            clone._aggregate_op = self._aggregate_op

            memo[id(self)] = clone

            return clone

        # If this object has been copied before, return the first copy made
        else:
            return existing

    def __repr__(self):

        # Append each nuclide in the aggregate to the string
        string = '{0}('.format(self.aggregate_op)
        names = [nuclide.name if isinstance(nuclide, Nuclide) else str(nuclide)
                 for nuclide in self.nuclides]
        string += ', '.join(map(str, names)) + ')'
        return string

    @property
    def nuclides(self):
        return self._nuclides

    @property
    def aggregate_op(self):
        return self._aggregate_op

    @nuclides.setter
    def nuclides(self, nuclides):
        cv.check_iterable_type('nuclides', nuclides,
                               (basestring, Nuclide, CrossNuclide))
        self._nuclides = nuclides

    @aggregate_op.setter
    def aggregate_op(self, aggregate_op):
        cv.check_type('aggregate_op', aggregate_op, basestring)
        cv.check_value('aggregate_op', aggregate_op, _TALLY_AGGREGATE_OPS)
        self._aggregate_op = aggregate_op


class AggregateFilter(object):
    """A special-purpose tally filter used to encapsulate an aggregate of a
    subset or all of a tally filter's bins for tally aggregation.

    Parameters
    ----------
    aggregate_filter : Filter or CrossFilter
        The filter included in the aggregation
    bins : Iterable of tuple
        The filter bins included in the aggregation
    aggregate_op : str
        The tally aggregation operator (e.g., 'sum', 'avg', etc.) used
        to aggregate across a tally filter's bins with this AggregateFilter

    Attributes
    ----------
    type : str
        The type of the aggregatefilter (e.g., 'sum(energy)', 'sum(cell)')
    aggregate_filter : filter
        The filter included in the aggregation
    aggregate_op : str
        The tally aggregation operator (e.g., 'sum', 'avg', etc.) used
        to aggregate across a tally filter's bins with this AggregateFilter
    bins : Iterable of tuple
        The filter bins included in the aggregation
    num_bins : Integral
        The number of filter bins (always 1 if aggregate_filter is defined)
    stride : Integral
        The number of filter, nuclide and score bins within each of this
        aggregatefilter's bins.

    """

    def __init__(self, aggregate_filter=None, bins=None, aggregate_op=None):

        self._type = '{0}({1})'.format(aggregate_op, aggregate_filter.type)
        self._bins = None
        self._stride = None

        self._aggregate_filter = None
        self._aggregate_op = None

        if aggregate_filter is not None:
            self.aggregate_filter = aggregate_filter
        if bins is not None:
            self.bins = bins
        if aggregate_op is not None:
            self.aggregate_op = aggregate_op

    def __hash__(self):
        return hash((self.type, self.bins, self.aggregate_op))

    def __eq__(self, other):
        return str(other) == str(self)

    def __ne__(self, other):
        return not self == other

    def __repr__(self):
        string = 'AggregateFilter\n'
        string += '{0: <16}{1}{2}\n'.format('\tType', '=\t', self.type)
        string += '{0: <16}{1}{2}\n'.format('\tBins', '=\t', self.bins)
        return string

    def __deepcopy__(self, memo):
        existing = memo.get(id(self))

        # If this is the first time we have tried to copy this object, create a copy
        if existing is None:
            clone = type(self).__new__(type(self))
            clone._type = self.type
            clone._aggregate_filter = self.aggregate_filter
            clone._aggregate_op = self.aggregate_op
            clone._bins = self._bins
            clone._stride = self.stride

            memo[id(self)] = clone

            return clone

        # If this object has been copied before, return the first copy made
        else:
            return existing

    @property
    def aggregate_filter(self):
        return self._aggregate_filter

    @property
    def aggregate_op(self):
        return self._aggregate_op

    @property
    def type(self):
        return self._type

    @property
    def bins(self):
        return self._bins

    @property
    def num_bins(self):
        return 1 if self.aggregate_filter else 0

    @property
    def stride(self):
        return self._stride

    @type.setter
    def type(self, filter_type):
        if filter_type not in _FILTER_TYPES.values():
            msg = 'Unable to set AggregateFilter type to "{0}" since it ' \
                  'is not one of the supported types'.format(filter_type)
            raise ValueError(msg)

        self._type = filter_type

    @aggregate_filter.setter
    def aggregate_filter(self, aggregate_filter):
        cv.check_type('aggregate_filter', aggregate_filter, (Filter, CrossFilter))
        self._aggregate_filter = aggregate_filter

    @bins.setter
    def bins(self, bins):
        cv.check_iterable_type('bins', bins, (Integral, tuple))
        self._bins = bins

    @aggregate_op.setter
    def aggregate_op(self, aggregate_op):
        cv.check_type('aggregate_op', aggregate_op, basestring)
        cv.check_value('aggregate_op', aggregate_op, _TALLY_AGGREGATE_OPS)
        self._aggregate_op = aggregate_op

    @stride.setter
    def stride(self, stride):
        self._stride = stride

    def get_bin_index(self, filter_bin):
        """Returns the index in the AggregateFilter for some bin.

        Parameters
        ----------
        filter_bin : Integral or tuple of Real
            A tuple of value(s) corresponding to the bin of interest in
            the aggregated filter. The bin is the integer ID for 'material',
            'surface', 'cell', 'cellborn', and 'universe' Filters. The bin
            is the integer cell instance ID for 'distribcell' Filters. The
            bin is a 2-tuple of floats for 'energy' and 'energyout' filters
            corresponding to the energy boundaries of the bin of interest.
            The bin is a (x,y,z) 3-tuple for 'mesh' filters corresponding to
            the mesh cell of interest.

        Returns
        -------
        filter_index : Integral
             The index in the Tally data array for this filter bin. For an
             AggregateTally the filter bin index is always unity.

        Raises
        ------
        ValueError
            When the filter_bin is not part of the aggregated filter's bins

        """

        if filter_bin not in self.bins:
            msg = 'Unable to get the bin index for AggregateFilter since ' \
                  '"{0}" is not one of the bins'.format(filter_bin)
            raise ValueError(msg)
        else:
            return 0

    def get_pandas_dataframe(self, datasize, summary=None):
        """Builds a Pandas DataFrame for the AggregateFilter's bins.

        This method constructs a Pandas DataFrame object for the AggregateFilter
        with columns annotated by filter bin information. This is a helper
        method for the Tally.get_pandas_dataframe(...) method.

        Parameters
        ----------
        datasize : Integral
            The total number of bins in the tally corresponding to this filter
        summary : None or Summary
            An optional Summary object to be used to construct columns for
            distribcell tally filters (default is None). NOTE: This parameter
            is not used by the AggregateFilter and simply mirrors the method
            signature for the CrossFilter.

        Returns
        -------
        pandas.DataFrame
            A Pandas DataFrame with columns of strings that characterize the
            aggregatefilter's bins. Each entry in the DataFrame will include
            one or more aggregation operations used to construct the
            aggregatefilter's bins. The number of rows in the DataFrame is the
            same as the total number of bins in the corresponding tally, with
            the filter bins appropriately tiled to map to the corresponding
            tally bins.

        See also
        --------
        Tally.get_pandas_dataframe(), Filter.get_pandas_dataframe(),
        CrossFilter.get_pandas_dataframe()

        """

        import pandas as pd

        # Construct a sring representing the filter aggregation
        aggregate_bin = '{0}('.format(self.aggregate_op)
        aggregate_bin += ', '.join(map(str, self.bins)) + ')'

        # Construct NumPy array of bin repeated for each element in dataframe
        aggregate_bin_array = np.array([aggregate_bin])
        aggregate_bin_array = np.repeat(aggregate_bin_array, datasize)

        # Construct Pandas DataFrame for the AggregateFilter
        df = pd.DataFrame({self.aggregate_filter.type: aggregate_bin_array})
        return df