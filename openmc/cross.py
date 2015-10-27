import sys

from openmc import Filter, Nuclide
from openmc.filter import _FILTER_TYPES
import openmc.checkvalue as cv


if sys.version_info[0] >= 3:
    basestring = str

# Acceptable tally arithmetic binary operations
_TALLY_ARITHMETIC_OPS = ['+', '-', '*', '/', '^']


class CrossScore(object):
    """A special-purpose tally score used to encapsulate all combinations of two
    tally's scores as an outer product for tally arithmetic.

    Parameters
    ----------
    left_score : str or CrossScore
        The left score in the outer product
    right_score : str or CrossScore
        The right score in the outer product
    binary_op : str
        The tally arithmetic binary operator (e.g., '+', '-', etc.) used to
        combine two tally's scores with this CrossNuclide

    Attributes
    ----------
    left_score : str or CrossScore
        The left score in the outer product
    right_score : str or CrossScore
        The right score in the outer product
    binary_op : str
        The tally arithmetic binary operator (e.g., '+', '-', etc.) used to
        combine two tally's scores with this CrossNuclide

    """

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

    def __hash__(self):
        return hash(str(self))

    def __eq__(self, other):
        return str(other) == str(self)

    def __ne__(self, other):
        return not self == other

    def __deepcopy__(self, memo):
        existing = memo.get(id(self))

        # If this is the first time we have tried to copy this object, create a copy
        if existing is None:
            clone = type(self).__new__(type(self))
            clone._left_score = self.left_score
            clone._right_score = self.right_score
            clone._binary_op = self.binary_op

            memo[id(self)] = clone

            return clone

        # If this object has been copied before, return the first copy made
        else:
            return existing

    def __repr__(self):
        string = '({0} {1} {2})'.format(self.left_score,
                                        self.binary_op, self.right_score)
        return string

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
        cv.check_type('left_score', left_score, (basestring, CrossScore))
        self._left_score = left_score

    @right_score.setter
    def right_score(self, right_score):
        cv.check_type('right_score', right_score, (basestring, CrossScore))
        self._right_score = right_score

    @binary_op.setter
    def binary_op(self, binary_op):
        cv.check_type('binary_op', binary_op, (basestring, CrossScore))
        cv.check_value('binary_op', binary_op, _TALLY_ARITHMETIC_OPS)
        self._binary_op = binary_op


class CrossNuclide(object):
    """A special-purpose nuclide used to encapsulate all combinations of two
    tally's nuclides as an outer product for tally arithmetic.

    Parameters
    ----------
    left_nuclide : Nuclide or CrossNuclide
        The left nuclide in the outer product
    right_nuclide : Nuclide or CrossNuclide
        The right nuclide in the outer product
    binary_op : str
        The tally arithmetic binary operator (e.g., '+', '-', etc.) used to
        combine two tally's nuclides with this CrossNuclide

    Attributes
    ----------
    left_nuclide : Nuclide or CrossNuclide
        The left nuclide in the outer product
    right_nuclide : Nuclide or CrossNuclide
        The right nuclide in the outer product
    binary_op : str
        The tally arithmetic binary operator (e.g., '+', '-', etc.) used to
        combine two tally's nuclides with this CrossNuclide

    """

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

    def __hash__(self):
        return hash(str(self))

    def __eq__(self, other):
        return str(other) == str(self)

    def __ne__(self, other):
        return not self == other

    def __deepcopy__(self, memo):
        existing = memo.get(id(self))

        # If this is the first time we have tried to copy this object, create a copy
        if existing is None:
            clone = type(self).__new__(type(self))
            clone._left_nuclide = self.left_nuclide
            clone._right_nuclide = self.right_nuclide
            clone._binary_op = self.binary_op

            memo[id(self)] = clone

            return clone

        # If this object has been copied before, return the first copy made
        else:
            return existing

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
        cv.check_type('left_nuclide', left_nuclide, (Nuclide, CrossNuclide))
        self._left_nuclide = left_nuclide

    @right_nuclide.setter
    def right_nuclide(self, right_nuclide):
        cv.check_type('right_nuclide', right_nuclide, (Nuclide, CrossNuclide))
        self._right_nuclide = right_nuclide

    @binary_op.setter
    def binary_op(self, binary_op):
        cv.check_type('binary_op', binary_op, basestring)
        cv.check_value('binary_op', binary_op, _TALLY_ARITHMETIC_OPS)
        self._binary_op = binary_op


class CrossFilter(object):
    """A special-purpose filter used to encapsulate all combinations of two
    tally's filter bins as an outer product for tally arithmetic.

    Parameters
    ----------
    left_filter : Filter or CrossFilter
        The left filter in the outer product
    right_filter : Filter or CrossFilter
        The right filter in the outer product
    binary_op : str
        The tally arithmetic binary operator (e.g., '+', '-', etc.) used to
        combine two tally's filter bins with this CrossFilter

    Attributes
    ----------
    left_filter : Filter or CrossFilter
        The left filter in the outer product
    right_filter : Filter or CrossFilter
        The right filter in the outer product
    binary_op : str
        The tally arithmetic binary operator (e.g., '+', '-', etc.) used to
        combine two tally's filter bins with this CrossFilter

    """

    def __init__(self, left_filter=None, right_filter=None, binary_op=None):

        left_type = left_filter.type
        right_type = right_filter.type
        self._type = '({0} {1} {2})'.format(left_type, binary_op, right_type)

        self._bins = {}
        self._stride = None

        self._left_filter = None
        self._right_filter = None
        self._binary_op = None

        if left_filter is not None:
            self.left_filter = left_filter
            self._bins['left'] = left_filter.bins
        if right_filter is not None:
            self.right_filter = right_filter
            self._bins['right'] = right_filter.bins
        if binary_op is not None:
            self.binary_op = binary_op

    def __hash__(self):
        return hash((self.left_filter, self.right_filter))

    def __eq__(self, other):
        return str(other) == str(self)

    def __ne__(self, other):
        return not self == other

    def __repr__(self):

        string = 'CrossFilter\n'
        filter_type = '({0} {1} {2})'.format(self.left_filter.type,
                                             self.binary_op,
                                             self.right_filter.type)
        filter_bins = '({0} {1} {2})'.format(self.left_filter.bins,
                                             self.binary_op,
                                             self.right_filter.bins)
        string += '{0: <16}{1}{2}\n'.format('\tType', '=\t', filter_type)
        string += '{0: <16}{1}{2}\n'.format('\tBins', '=\t', filter_bins)
        return string

    def __deepcopy__(self, memo):
        existing = memo.get(id(self))

        # If this is the first time we have tried to copy this object, create a copy
        if existing is None:
            clone = type(self).__new__(type(self))
            clone._left_filter = self.left_filter
            clone._right_filter = self.right_filter
            clone._binary_op = self.binary_op
            clone._type = self.type
            clone._bins = self.bins
            clone._num_bins = self.num_bins
            clone._stride = self.stride

            memo[id(self)] = clone

            return clone

        # If this object has been copied before, return the first copy made
        else:
            return existing

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
        return self._type

    @property
    def bins(self):
        return self._bins['left'], self._bins['right']

    @property
    def num_bins(self):
        if self.left_filter is not None and self.right_filter is not None:
            return self.left_filter.num_bins * self.right_filter.num_bins
        else:
            return 0

    @property
    def stride(self):
        return self._stride

    @type.setter
    def type(self, filter_type):
        if filter_type not in _FILTER_TYPES.values():
            msg = 'Unable to set Filter type to "{0}" since it is not one ' \
                  'of the supported types'.format(type)
            raise ValueError(msg)

        self._type = filter_type

    @left_filter.setter
    def left_filter(self, left_filter):
        cv.check_type('left_filter', left_filter, (Filter, CrossFilter))
        self._left_filter = left_filter
        self._bins['left'] = left_filter.bins

    @right_filter.setter
    def right_filter(self, right_filter):
        cv.check_type('right_filter', right_filter, (Filter, CrossFilter))
        self._right_filter = right_filter
        self._bins['right'] = right_filter.bins

    @binary_op.setter
    def binary_op(self, binary_op):
        cv.check_type('binary_op', binary_op, basestring)
        cv.check_value('binary_op', binary_op, _TALLY_ARITHMETIC_OPS)
        self._binary_op = binary_op

    @stride.setter
    def stride(self, stride):
        self._stride = stride

    def get_bin_index(self, filter_bin):
        """Returns the index in the CrossFilter for some bin.

        Parameters
        ----------
        filter_bin : 2-tuple
            A 2-tuple where each value corresponds to the bin of interest
            in the left and right filter, respectively. A bin is the integer
            ID for 'material', 'surface', 'cell', 'cellborn', and 'universe'
            Filters. The bin is an integer for the cell instance ID for
            'distribcell' Filters. The bin is a 2-tuple of floats for 'energy'
            and 'energyout' filters corresponding to the energy boundaries of
            the bin of interest.  The bin is a (x,y,z) 3-tuple for 'mesh'
            filters corresponding to the mesh cell of interest.

        Returns
        -------
        filter_index : Integral
             The index in the Tally data array for this filter bin.

        """

        left_index = self.left_filter.get_bin_index(filter_bin[0])
        right_index = self.right_filter.get_bin_index(filter_bin[0])
        filter_index = left_index * self.right_filter.num_bins + right_index
        return filter_index

    def get_pandas_dataframe(self, datasize, summary=None):
        """Builds a Pandas DataFrame for the CrossFilter's bins.

        This method constructs a Pandas DataFrame object for the CrossFilter
        with columns annotated by filter bin information. This is a helper
        method for the Tally.get_pandas_dataframe(...) method. This method
        recursively builds and concatenates Pandas DataFrames for the left
        and right filters and crossfilters.

        This capability has been tested for Pandas >=0.13.1. However, it is
        recommended to use v0.16 or newer versions of Pandas since this method
        uses Pandas' Multi-index functionality.

        Parameters
        ----------
        datasize : Integral
            The total number of bins in the tally corresponding to this filter
        summary : None or Summary
            An optional Summary object to be used to construct columns for
            distribcell tally filters (default is None). The geometric
            information in the Summary object is embedded into a Multi-index
            column with a geometric "path" to each distribcell instance.
            NOTE: This option requires the OpenCG Python package.

        Returns
        -------
        pandas.DataFrame
            A Pandas DataFrame with columns of strings that characterize the
            crossfilter's bins. Each entry in the DataFrame will include one
            or more binary operations used to construct the crossfilter's bins.
            The number of rows in the DataFrame is the same as the total number
            of bins in the corresponding tally, with the filter bins
            appropriately tiled to map to the corresponding tally bins.

        See also
        --------
        Tally.get_pandas_dataframe(), Filter.get_pandas_dataframe()

        """

        # If left and right filters are identical, do not combine bins
        if self.left_filter == self.right_filter:
            df = self.left_filter.get_pandas_dataframe(datasize, summary)

        # If left and right filters are different, combine their bins
        else:
            left_df = self.left_filter.get_pandas_dataframe(datasize, summary)
            right_df = self.right_filter.get_pandas_dataframe(datasize, summary)
            left_df = left_df.astype(str)
            right_df = right_df.astype(str)
            df = '(' + left_df + ' ' + self.binary_op + ' ' + right_df + ')'

        return df
