import copy

from openmc import Mesh
from openmc.checkvalue import *
from openmc.constants import *


class Filter(object):

    # Initialize Filter class attributes
    def __init__(self, type=None, bins=None):

        self.type = type
        self._num_bins = 0
        self.bins = bins
        self._mesh = None
        self._offset = -1
        self._stride = None


    def __eq__(self, filter2):

        # Check type
        if self._type != filter2._type:
            return False

        # Check number of bins
        elif len(self._bins) != len(filter2._bins):
            return False

        # Check bin edges
        elif not np.allclose(self._bins, filter2._bins):
            return False

        else:
            return True


    def __hash__(self):
        hashable = []
        hashable.append(self._type)
        hashable.append(self._bins)
        return hash(tuple(hashable))


    def __deepcopy__(self, memo):

        existing = memo.get(id(self))

        # If this is the first time we have tried to copy this object, create a copy
        if existing is None:

            clone = type(self).__new__(type(self))
            clone._type = self._type
            clone._bins = copy.deepcopy(self._bins, memo)
            clone._num_bins = self._num_bins
            clone._mesh = copy.deepcopy(self._mesh, memo)
            clone._offset = self._offset
            clone._stride = self._stride

            memo[id(self)] = clone

            return clone

        # If this object has been copied before, return the first copy made
        else:
            return existing


    @property
    def type(self):
        return self._type


    @property
    def bins(self):
        return self._bins


    @property
    def num_bins(self):
        return self._num_bins


    @property
    def mesh(self):
        return self._mesh


    @property
    def offset(self):
        return self._offset


    @property
    def stride(self):
        return self._stride


    @type.setter
    def type(self, type):

        if type is None:
            self._type = type

        elif not type in FILTER_TYPES.values():
            msg = 'Unable to set Filter type to {0} since it is not one ' \
                  'of the supported types'.format(type)
            raise ValueError(msg)

        self._type = type


    @bins.setter
    def bins(self, bins):

        if bins is None:
            self.num_bins = 0

        elif self._type is None:
            msg = 'Unable to set bins for Filter to {0} since ' \
                  'the Filter type has not yet been set'.format(bins)
            raise ValueError(msg)

        # If the bin edge is a single value, it is a Cell, Material, etc. ID
        if not isinstance(bins, (tuple, list, np.ndarray)):
            bins = [bins]

        # If the bins are in a collection, convert it to a list
        else:
            bins = list(bins)

        if self._type in ['cell', 'cellborn', 'surface', 'material',
                          'universe', 'distribcell']:

            for edge in bins:

                if not is_integer(edge):
                    msg = 'Unable to add bin {0} to a {1} Filter since ' \
                          'it is a non-integer'.format(edge, self._type)
                    raise ValueError(msg)

                elif edge < 0:
                    msg = 'Unable to add bin  {0} to a {1} Filter since ' \
                          'it is a negative integer'.format(edge, self._type)
                    raise ValueError(msg)


        elif self._type in ['energy', 'energyout']:

            for edge in bins:

                if not is_integer(edge) and not is_float(edge):
                    msg = 'Unable to add bin edge {0} to {1} Filter since ' \
                          'it is a non-integer or floating point ' \
                          'value'.format(edge, self._type)
                    raise ValueError(msg)

                elif edge < 0.:
                    msg = 'Unable to add bin edge {0} to {1} Filter since it ' \
                          'is a negative value'.format(edge, self._type)
                    raise ValueError(msg)

            # Check that bin edges are monotonically increasing
            for index in range(len(bins)):

                if index > 0 and bins[index] < bins[index-1]:
                    msg = 'Unable to add bin edges {0} to {1} Filter since ' \
                          'they are not monotonically ' \
                          'increasing'.format(bins, self._type)
                    raise ValueError(msg)


        # mesh filters
        elif self._type == 'mesh':

            if not len(bins) == 1:
                msg = 'Unable to add bins {0} to a mesh Filter since ' \
                      'only a single mesh can be used per tally'.format(bins)
                raise ValueError(msg)

            elif not is_integer(bins[0]):
                msg = 'Unable to add bin {0} to mesh Filter since it ' \
                       'is a non-integer'.format(bins[0])
                raise ValueError(msg)

            elif bins[0] < 0:
                msg = 'Unable to add bin {0} to mesh Filter since it ' \
                       'is a negative integer'.format(bins[0])
                raise ValueError(msg)

        # If all error checks passed, add bin edges
        self._bins = bins


    # FIXME
    @num_bins.setter
    def num_bins(self, num_bins):

        if not is_integer(num_bins) or num_bins < 0:
            msg = 'Unable to set the number of bins {0} for a {1} Filter ' \
                  'since it is not a positive ' \
                  'integer'.format(num_bins, self._type)
            raise ValueError(msg)

        self._num_bins = num_bins


    @mesh.setter
    def mesh(self, mesh):

        if not isinstance(mesh, Mesh):
            msg = 'Unable to set Mesh to {0} for Filter since it is not a ' \
                  'Mesh object'.format(mesh)
            raise ValueError(msg)

        self._mesh = mesh
        self.type = 'mesh'
        self.bins = self._mesh._id


    @offset.setter
    def offset(self, offset):

        if not is_integer(offset):
            msg = 'Unable to set offset {0} for a {1} Filter since it is a ' \
                  'non-integer value'.format(offset, self._type)
            raise ValueError(msg)

        self._offset = offset


    @stride.setter
    def stride(self, stride):

        if not is_integer(stride):
            msg = 'Unable to set stride {0} for a {1} Filter since it is a ' \
                  'non-integer value'.format(stride, self._type)
            raise ValueError(msg)

        if stride < 0:
            msg = 'Unable to set stride {0} for a {1} Filter since it is a ' \
                  'negative value'.format(stride, self._type)
            raise ValueError(msg)

        self._stride = stride


    def can_merge(self, filter):

        if not isinstance(filter, Filter):
            return False

        # Filters must be of the same type
        elif self.type != filter.type:
            return False

        # Distribcell filters cannot have more than one bin
        elif self.type == 'distribcell':
            return False

        # Mesh filters cannot have more than one bin
        elif self.type == 'mesh':
            return False

        # Different energy bins are not mergeable
        elif 'energy' in self.type:
            return False

        else:
            return True


    def merge(self, filter):

        if not self.can_merge(filter):
            msg = 'Unable to merge {0} with {1} filters'.format(self._type, filter._type)
            raise ValueError(msg)

        # Create deep copy of filter to return as merged filter
        merged_filter = copy.deepcopy(self)

        # Merge unique filter bins
        merged_bins = list(set(self._bins + filter._bins))
        merged_filter.bins = merged_bins
        merged_filter.num_bins = len(merged_bins)

        return merged_filter


    def get_bin_index(self, bin):

        try:
            index = self._bins.index(bin)

        except ValueError:
            msg = 'Unable to get the bin index for Filter since {0} ' \
                  'is not one of the bins'.format(bin)
            raise ValueError(msg)

        return index


    def __repr__(self):

        string = 'Filter\n'
        string += '{0: <16}{1}{2}\n'.format('\tType', '=\t', self._type)
        string += '{0: <16}{1}{2}\n'.format('\tBins', '=\t', self._bins)
        string += '{0: <16}{1}{2}\n'.format('\tOffset', '=\t', self._offset)
        return string
