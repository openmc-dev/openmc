from openmc import Filter, Nuclide
from openmc.checkvalue import check_type


class _CrossScore(object):
    """A special-purpose tally score used to encapsulate all combinations of two
    tally's scores as a cross product for tally arithmetic.

    Parameters
    ----------
    left_score : str or _CrossScore
        The left score in the cross product.
    right_score : str or _CrossScore
        The right score in the cross product.
    binary_op : str
        The tally arithmetic binary operator (e.g., '+', '-', etc.) used to
        combine two tally's scores with this _CrossNuclide.

    Attributes
    ----------
    left_score : str or _CrossScore
        The left score in the cross product.
    right_score : str or _CrossScore
        The right score in the cross product.
    binary_op : str
        The tally arithmetic binary operator (e.g., '+', '-', etc.) used to
        combine two tally's scores with this _CrossNuclide.

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
        check_type('left score', left_score, (str, _CrossScore))
        self._left_score = left_score

    @right_score.setter
    def right_score(self, right_score):
        check_type('right score', right_score, (str, _CrossScore))
        self._right_score = right_score

    @binary_op.setter
    def binary_op(self, binary_op):
        check_type('binary op', binary_op, str)
        self._binary_op = binary_op

    def __repr__(self):
        string = '({0} {1} {2})'.format(self.left_score,
                                        self.binary_op, self.right_score)
        return string


class _CrossNuclide(object):
    """A special-purpose nuclide used to encapsulate all combinations of two
    tally's nuclides as a cross product for tally arithmetic.

    Parameters
    ----------
    left_nuclide : Nuclide or _CrossNuclide
        The left nuclide in the cross product.
    right_nuclide : Nuclide or _CrossNuclide
        The right nuclide in the cross product.
    binary_op : str
        The tally arithmetic binary operator (e.g., '+', '-', etc.) used to
        combine two tally's nuclides with this _CrossNuclide.

    Attributes
    ----------
    left_nuclide : Nuclide or _CrossNuclide
        The left nuclide in the cross product.
    right_nuclide : Nuclide or _CrossNuclide
        The right nuclide in the cross product.
    binary_op : str
        The tally arithmetic binary operator (e.g., '+', '-', etc.) used to
        combine two tally's nuclides with this _CrossNuclide.

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
        check_type('left nuclide', left_nuclide, (Nuclide, _CrossNuclide))
        self._left_nuclide = left_nuclide

    @right_nuclide.setter
    def right_nuclide(self, right_nuclide):
        check_type('right nuclide', right_nuclide, (Nuclide, _CrossNuclide))
        self._right_nuclide = right_nuclide

    @binary_op.setter
    def binary_op(self, binary_op):
        check_type('binary op', binary_op, str)
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
    """A special-purpose filter used to encapsulate all combinations of two
    tally's filter bins as a cross product for tally arithmetic.

    Parameters
    ----------
    left_filter : Filter or _CrossFilter
        The left filter in the cross product.
    right_filter : Filter or _CrossFilter
        The right filter in the cross product.
    binary_op : str
        The tally arithmetic binary operator (e.g., '+', '-', etc.) used to
        combine two tally's filter bins with this _CrossFilter.

    Attributes
    ----------
    left_filter : Filter or _CrossFilter
        The left filter in the cross product.
    right_filter : Filter or _CrossFilter
        The right filter in the cross product.
    binary_op : str
        The tally arithmetic binary operator (e.g., '+', '-', etc.) used to
        combine two tally's filter bins with this _CrossFilter.

    """

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
    def num_bins(self):
        return self.left_filter.num_bins * self.right_filter.num_bins

    @property
    def stride(self):
        return self.left_filter.stride * self.right_filter.stride

    @left_filter.setter
    def left_filter(self, left_filter):
        check_type('left filter', left_filter, (Filter, _CrossFilter))
        self._left_filter = left_filter

    @right_filter.setter
    def right_filter(self, right_filter):
        check_type('right filter', right_filter, (Filter, _CrossFilter))
        self._right_filter = right_filter

    @binary_op.setter
    def binary_op(self, binary_op):
        check_type('binary op', binary_op, str)
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


    def split_filters(self):

        split_filters = []

        # If left Filter is not a CrossFilter, simply append to list
        if isinstance(self.left_filter, Filter):
            split_filters.append(self.left_filter)
        # Recursively descend CrossFilter tree to collect all Filters
        else:
            split_filters.extend(self.left_filter.split_filters())

        # If right Filter is not a CrossFilter, simply append to list
        if isinstance(self.right_filter, Filter):
            split_filters.append(self.right_filter)
        # Recursively descend CrossFilter tree to collect all Filters
        else:
            split_filters.extend(self.right_filter.split_filters())

        return split_filters