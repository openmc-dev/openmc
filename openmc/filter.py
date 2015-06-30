from collections import Iterable
import copy
from numbers import Real, Integral

import numpy as np

from openmc import Mesh
from openmc.constants import *
from openmc.checkvalue import check_type

class Filter(object):
    """A filter used to constrain a tally to a specific criterion, e.g. only tally
    events when the particle is in a certain cell and energy range.

    Parameters
    ----------
    type : str
        The type of the tally filter. Acceptable values are "universe",
        "material", "cell", "cellborn", "surface", "mesh", "energy",
        "energyout", and "distribcell".
    bins : int or Iterable of int or Iterable of float
        The bins for the filter. This takes on different meaning for different
        filters.

    Attributes
    ----------
    type : str
        The type of the tally filter.
    bins : int or Iterable of int or Iterable of float
        The bins for the filter

    """

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
        return hash((self._type, self._bins))

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
        elif type not in FILTER_TYPES.values():
            msg = 'Unable to set Filter type to "{0}" since it is not one ' \
                  'of the supported types'.format(type)
            raise ValueError(msg)

        self._type = type

    @bins.setter
    def bins(self, bins):
        if bins is None:
            self.num_bins = 0
        elif self._type is None:
            msg = 'Unable to set bins for Filter to "{0}" since ' \
                  'the Filter type has not yet been set'.format(bins)
            raise ValueError(msg)

        # If the bin edge is a single value, it is a Cell, Material, etc. ID
        if not isinstance(bins, Iterable):
            bins = [bins]

        # If the bins are in a collection, convert it to a list
        else:
            bins = list(bins)

        if self._type in ['cell', 'cellborn', 'surface', 'material',
                          'universe', 'distribcell']:
            for edge in bins:
                if not isinstance(edge, Integral):
                    msg = 'Unable to add bin "{0}" to a {1} Filter since ' \
                          'it is not an integer'.format(edge, self._type)
                    raise ValueError(msg)
                elif edge < 0:
                    msg = 'Unable to add bin "{0}" to a {1} Filter since ' \
                          'it is negative'.format(edge, self._type)
                    raise ValueError(msg)

        elif self._type in ['energy', 'energyout']:
            for edge in bins:
                if not isinstance(edge, Real):
                    msg = 'Unable to add bin edge "{0}" to a {1} Filter ' \
                          'since it is a non-integer or floating point ' \
                          'value'.format(edge, self._type)
                    raise ValueError(msg)
                elif edge < 0.:
                    msg = 'Unable to add bin edge "{0}" to a {1} Filter ' \
                          'since it is a negative value'.format(edge, self._type)
                    raise ValueError(msg)

            # Check that bin edges are monotonically increasing
            for index in range(len(bins)):
                if index > 0 and bins[index] < bins[index-1]:
                    msg = 'Unable to add bin edges "{0}" to a {1} Filter ' \
                          'since they are not monotonically ' \
                          'increasing'.format(bins, self._type)
                    raise ValueError(msg)

        # mesh filters
        elif self._type == 'mesh':
            if not len(bins) == 1:
                msg = 'Unable to add bins "{0}" to a mesh Filter since ' \
                      'only a single mesh can be used per tally'.format(bins)
                raise ValueError(msg)
            elif not isinstance(bins[0], Integral):
                msg = 'Unable to add bin "{0}" to mesh Filter since it ' \
                       'is a non-integer'.format(bins[0])
                raise ValueError(msg)
            elif bins[0] < 0:
                msg = 'Unable to add bin "{0}" to mesh Filter since it ' \
                       'is a negative integer'.format(bins[0])
                raise ValueError(msg)

        # If all error checks passed, add bin edges
        self._bins = bins

    # FIXME
    @num_bins.setter
    def num_bins(self, num_bins):
        if not isinstance(num_bins, Integral) or num_bins < 0:
            msg = 'Unable to set the number of bins "{0}" for a {1} Filter ' \
                  'since it is not a positive ' \
                  'integer'.format(num_bins, self._type)
            raise ValueError(msg)

        self._num_bins = num_bins

    @mesh.setter
    def mesh(self, mesh):
        check_type('filter mesh', mesh, Mesh)

        self._mesh = mesh
        self.type = 'mesh'
        self.bins = self._mesh._id

    @offset.setter
    def offset(self, offset):
        check_type('filter offset', offset, Integral)
        self._offset = offset

    @stride.setter
    def stride(self, stride):
        check_type('filter stride', stride, Integral)
        if stride < 0:
            msg = 'Unable to set stride "{0}" for a {1} Filter since it is a ' \
                  'negative value'.format(stride, self._type)
            raise ValueError(msg)

        self._stride = stride

    def can_merge(self, filter):
        """Determine if filter can be merged with another.

        Parameters
        ----------
        filter : Filter
            Filter to compare with

        Returns
        -------
        bool
            Whether the filter can be merged

        """

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
        """Merge this filter with another.

        Parameters
        ----------
        filter : Filter
            Filter to merge with

        Returns
        -------
        merged_filter : Filter
            Filter resulting from the merge

        """

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

    def get_bin_index(self, filter_bin):
        """Returns the index in the Filter for some bin.

        Parameters
        ----------
        filter_bin : int or tuple
            The bin is the integer ID for 'material', 'surface', 'cell',
            'cellborn', and 'universe' Filters. The bin is an integer for the
            cell instance ID for 'distribcell' Filters. The bin is a 2-tuple of
            floats for 'energy' and 'energyout' filters corresponding to the
            energy boundaries of the bin of interest.  The bin is a (x,y,z)
            3-tuple for 'mesh' filters corresponding to the mesh cell of
            interest.

        Returns
        -------
        filter_index : int
             The index in the Tally data array for this filter bin.

        """

        try:
            # Filter bins for a mesh are an (x,y,z) tuple
            if self.type == 'mesh':
                # Convert (x,y,z) to a single bin -- this is similar to
                # subroutine mesh_indices_to_bin in openmc/src/mesh.F90.
                if (len(self.mesh.dimension) == 3):
                    nx, ny, nz = self.mesh.dimension
                    val = (filter_bin[0] - 1) * ny * nz + \
                          (filter_bin[1] - 1) * nz + \
                          (filter_bin[2] - 1)
                else:
                    nx, ny = self.mesh.dimension
                    val = (filter_bin[0] - 1) * ny + \
                          (filter_bin[1] - 1)

                filter_index = val

            # Use lower energy bound to find index for energy Filters
            elif self.type in ['energy', 'energyout']:
                val = self.bins.index(filter_bin[0])
                filter_index = val

            # Filter bins for distribcell are the "IDs" of each unique placement
            # of the Cell in the Geometry (integers starting at 0)
            elif self._type == 'distribcell':
                filter_index = filter_bin

            # Use ID for all other Filters (e.g., material, cell, etc.)
            else:
                val = self.bins.index(filter_bin)
                filter_index = val

        except ValueError:
            msg = 'Unable to get the bin index for Filter since "{0}" ' \
                      'is not one of the bins'.format(filter_bin)
            raise ValueError(msg)

        return filter_index

    def __repr__(self):
        string = 'Filter\n'
        string += '{0: <16}{1}{2}\n'.format('\tType', '=\t', self._type)
        string += '{0: <16}{1}{2}\n'.format('\tBins', '=\t', self._bins)
        string += '{0: <16}{1}{2}\n'.format('\tOffset', '=\t', self._offset)
        return string
