from abc import ABCMeta, abstractproperty
from collections import Iterable, OrderedDict
import copy
from six import with_metaclass
from numbers import Real, Integral
import sys
from xml.etree import ElementTree as ET

import numpy as np

from openmc import Mesh
import openmc.checkvalue as cv


if sys.version_info[0] >= 3:
    basestring = str


_FILTER_TYPES = ['universe', 'material', 'cell', 'cellborn', 'surface',
                 'mesh', 'energy', 'energyout', 'mu', 'polar', 'azimuthal',
                 'distribcell', 'delayedgroup']

_CURRENT_NAMES = {1:  'x-min out', 2:  'x-min in',
                  3:  'x-max out', 4:  'x-max in',
                  5:  'y-min out', 6:  'y-min in',
                  7:  'y-max out', 8:  'y-max in',
                  9:  'z-min out', 10: 'z-min in',
                  11: 'z-max out', 12: 'z-max in'}


class FilterMeta(ABCMeta):
    def __new__(cls, name, bases, namespace, **kwargs):
        if not name.endswith('Filter'):
            raise ValueError("All filter class names must end with 'Filter'")
        namespace['short_name'] = name[:-6]
        return super(FilterMeta, cls).__new__(cls, name, bases, namespace,
                                              **kwargs)


class Filter(with_metaclass(FilterMeta, object)):
    """A constraint on transport events that can score a tally

    Parameters
    ----------
    bins : Integral or Iterable of Integral or Iterable of Real
        The bins for the filter. This takes on different meaning for different
        filters. See the OpenMC online documentation for more details.

    Attributes
    ----------
    bins : Integral or Iterable of Integral or Iterable of Real
        The bins for the filter
    num_bins : Integral
        The number of filter bins
    stride : Integral
        The number of filter, nuclide and score bins within each of this
        filter's bins.

    """

    def __init__(self, bins):
        self.bins = bins
        self._num_bins = 0
        self._stride = None

    def __eq__(self, other):
        if type(self) != type(other):
            return False
        elif len(self.bins) != len(other.bins):
            return False
        elif not np.allclose(self.bins, other.bins):
            return False
        else:
            return True

    def __ne__(self, other):
        return not self == other

    def __gt__(self, other):
        if type(self) != type(other):
            if self.short_name in _FILTER_TYPES and \
                other.short_name in _FILTER_TYPES:
                delta = _FILTER_TYPES.index(self.short_name) - \
                        _FILTER_TYPES.index(other.short_name)
                return delta > 0
            else:
                return False
        else:
            return max(self.bins) > max(other.bins)

    def __lt__(self, other):
        return not self > other

    def __hash__(self):
        return hash(repr(self))

    def __repr__(self):
        string = type(self).__name__ + '\n'
        string += '{0: <16}{1}{2}\n'.format('\tBins', '=\t', self.bins)
        return string

    @classmethod
    def recursive_subclasses(cls):
        subs = cls.__subclasses__()
        subsubs = [grand for s in subs for grand in s.__subclasses__()]
        return subs + subsubs

    @classmethod
    def from_hdf5(cls, group):
        if group['type'].value.decode() == cls.short_name.lower():
            out = cls(group['bins'].value)
            out.num_bins = group['n_bins'].value
            return out

        for subclass in cls.recursive_subclasses():
            if group['type'].value.decode() == subclass.short_name.lower():
                return subclass.from_hdf5(group)

        raise ValueError("Unrecognized Filter class: '"
                         + group['type'].value.decode() + "'")

    @property
    def bins(self):
        return self._bins

    @property
    def num_bins(self):
        return self._num_bins

    @property
    def stride(self):
        return self._stride

    @bins.setter
    def bins(self, bins):
        # We might be given a single value.  Stick it into a list.
        if not isinstance(bins, Iterable):
            bins = [bins]

        # If the bin is 0D numpy array, promote to 1D
        elif isinstance(bins, np.ndarray):
            if bins.shape == ():
                bins.shape = (1,)

        # If the bins are in a collection, convert it to a list
        else:
            bins = list(bins)

        # Convert to a numpy array
        bins = np.array(bins)

        # Check the bin values
        self.check_bins(bins)

        self._bins = bins

    @num_bins.setter
    def num_bins(self, num_bins):
        cv.check_type('filter num_bins', num_bins, Integral)
        cv.check_greater_than('filter num_bins', num_bins, 0, equality=True)
        self._num_bins = num_bins

    @stride.setter
    def stride(self, stride):
        cv.check_type('filter stride', stride, Integral)
        if stride < 0:
            msg = 'Unable to set stride "{0}" for a "{1}" Filter since it ' \
                  'is a negative value'.format(stride, self.type)
            raise ValueError(msg)

        self._stride = stride

    def check_bins(self, bins):
        """Make sure given bins are valid for this filter.

        Raises
        ------
        TypeError
        ValueError

        """
        pass

    def to_xml(self):
        element = ET.Element('filter')
        element.set('type', self.short_name.lower())
        element.set('bins', ' '.join(str(b) for b in self.bins))
        return element

    def can_merge(self, other):
        """Determine if filter can be merged with another.

        Parameters
        ----------
        other : openmc.Filter
            Filter to compare with

        Returns
        -------
        bool
            Whether the filter can be merged

        """

        if type(self) != type(other): return False

        return True

    def merge(self, other):
        """Merge this filter with another.

        Parameters
        ----------
        other : openmc.Filter
            Filter to merge with

        Returns
        -------
        merged_filter : openmc.Filter
            Filter resulting from the merge

        """

        if not self.can_merge(other):
            msg = 'Unable to merge "{0}" with "{1}" ' \
                  'filters'.format(self.type, other.type)
            raise ValueError(msg)

        # Merge unique filter bins
        merged_bins = np.concatenate((self.bins, other.bins))
        merged_bins = np.unique(merged_bins)

        # Create a new filter with these bins
        return type(self)(merged_bins)

    def is_subset(self, other):
        """Determine if another filter is a subset of this filter.

        If all of the bins in the other filter are included as bins in this
        filter, then it is a subset of this filter.

        Parameters
        ----------
        other : openmc.Filter
            The filter to query as a subset of this filter

        Returns
        -------
        bool
            Whether or not the other filter is a subset of this filter

        """

        if type(self) != type(other):
            return False

        for bin in other.bins:
            if bin not in self.bins:
                return False

        return True

    def get_bin_index(self, filter_bin):
        """Returns the index in the Filter for some bin.

        Parameters
        ----------
        filter_bin : Integral or tuple
            The bin is the integer ID for 'material', 'surface', 'cell',
            'cellborn', and 'universe' Filters. The bin is an integer for the
            cell instance ID for 'distribcell' Filters. The bin is a 2-tuple of
            floats for 'energy' and 'energyout' filters corresponding to the
            energy boundaries of the bin of interest. The bin is an (x,y,z)
            3-tuple for 'mesh' filters corresponding to the mesh cell of
            interest.

        Returns
        -------
        filter_index : Integral
             The index in the Tally data array for this filter bin.

        See also
        --------
        Filter.get_bin()

        """

        if not filter_bin in self.bins:
            msg = 'Unable to get the bin index for Filter since "{0}" ' \
                  'is not one of the bins'.format(filter_bin)
            raise ValueError(msg)

        return np.where(self.bins == filter_bin)[0][0]

    def get_bin(self, bin_index):
        """Returns the filter bin for some filter bin index.

        Parameters
        ----------
        bin_index : Integral
            The zero-based index into the filter's array of bins. The bin
            index for 'material', 'surface', 'cell', 'cellborn', and 'universe'
            filters corresponds to the ID in the filter's list of bins. For
            'distribcell' tallies the bin index necessarily can only be zero
            since only one cell can be tracked per tally. The bin index for
            'energy' and 'energyout' filters corresponds to the energy range of
            interest in the filter bins of energies. The bin index for 'mesh'
            filters is the index into the flattened array of (x,y) or (x,y,z)
            mesh cell bins.

        Returns
        -------
        bin : 1-, 2-, or 3-tuple of Real
            The bin in the Tally data array. The bin for 'material', surface',
            'cell', 'cellborn', 'universe' and 'distribcell' filters is a
            1-tuple of the ID corresponding to the appropriate filter bin.
            The bin for 'energy' and 'energyout' filters is a 2-tuple of the
            lower and upper energies bounding the energy interval for the filter
            bin. The bin for 'mesh' tallies is a 2-tuple or 3-tuple of the x,y
            or x,y,z mesh cell indices corresponding to the bin in a 2D/3D mesh.

        See also
        --------
        Filter.get_bin_index()

        """

        cv.check_type('bin_index', bin_index, Integral)
        cv.check_greater_than('bin_index', bin_index, 0, equality=True)
        cv.check_less_than('bin_index', bin_index, self.num_bins)

        # Return a 1-tuple of the bin.
        return (self.bins[bin_index],)

    def get_pandas_dataframe(self, data_size, distribcell_paths=True):
        """Builds a Pandas DataFrame for the Filter's bins.

        This method constructs a Pandas DataFrame object for the filter with
        columns annotated by filter bin information. This is a helper method for
        :meth:`Tally.get_pandas_dataframe`.

        This capability has been tested for Pandas >=0.13.1. However, it is
        recommended to use v0.16 or newer versions of Pandas since this method
        uses Pandas' Multi-index functionality.

        Parameters
        ----------
        data_size : Integral
            The total number of bins in the tally corresponding to this filter
        distribcell_paths : bool, optional
            Construct columns for distribcell tally filters (default is True).
            The geometric information in the Summary object is embedded into a
            Multi-index column with a geometric "path" to each distribcell
            instance. NOTE: This option assumes that all distribcell paths are
            of the same length and do not have the same universes and cells but
            different lattice cell indices.

        Returns
        -------
        pandas.DataFrame
            A Pandas DataFrame with columns of strings that characterize the
            filter's bins. The number of rows in the DataFrame is the same as
            the total number of bins in the corresponding tally, with the filter
            bin appropriately tiled to map to the corresponding tally bins.

            For 'cell', 'cellborn', 'surface', 'material', and 'universe'
            filters, the DataFrame includes a single column with the cell,
            surface, material or universe ID corresponding to each filter bin.

            For 'distribcell' filters, the DataFrame either includes:

            1. a single column with the cell instance IDs (without summary info)
            2. separate columns for the cell IDs, universe IDs, and lattice IDs
               and x,y,z cell indices corresponding to each (distribcell paths).

            For 'energy' and 'energyout' filters, the DataFrame includes one
            column for the lower energy bound and one column for the upper
            energy bound for each filter bin.

            For 'mesh' filters, the DataFrame includes three columns for the
            x,y,z mesh cell indices corresponding to each filter bin.

        Raises
        ------
        ImportError
            When Pandas is not installed

        See also
        --------
        Tally.get_pandas_dataframe(), CrossFilter.get_pandas_dataframe()

        """

        # Initialize Pandas DataFrame
        import pandas as pd
        df = pd.DataFrame()

        filter_bins = np.repeat(self.bins, self.stride)
        tile_factor = data_size / len(filter_bins)
        filter_bins = np.tile(filter_bins, tile_factor)
        df = pd.concat([df, pd.DataFrame(
            {self.short_name.lower() : filter_bins})])

        return df


class IntegralFilter(Filter):
    """A Filter with Integral bins"""
    @property
    def num_bins(self):
        return len(self.bins)

    @num_bins.setter
    def num_bins(self, num_bins):
        cv.check_type('filter num_bins', num_bins, Integral)
        cv.check_greater_than('filter num_bins', num_bins, 0, equality=True)
        self._num_bins = num_bins

    def check_bins(self, bins):
        cv.check_iterable_type('filter bins', bins, Integral)
        for edge in bins:
            cv.check_greater_than('filter bin', edge, 0, equality=True)


class UniverseFilter(IntegralFilter): pass
class MaterialFilter(IntegralFilter): pass
class CellFilter(IntegralFilter): pass
class CellbornFilter(IntegralFilter): pass


class SurfaceFilter(IntegralFilter):
    def get_pandas_dataframe(self, data_size, distribcell_paths=True):
        # Initialize Pandas DataFrame
        import pandas as pd
        df = pd.DataFrame()

        filter_bins = np.repeat(self.bins, self.stride)
        tile_factor = data_size / len(filter_bins)
        filter_bins = np.tile(filter_bins, tile_factor)
        filter_bins = [_CURRENT_NAMES[x] for x in filter_bins]
        df = pd.concat([df, pd.DataFrame(
            {self.short_name.lower() : filter_bins})])

        return df


class MeshFilter(Filter):
    def __init__(self, bins):
        self._mesh = None
        super(MeshFilter, self).__init__(bins)

    @property
    def mesh(self):
        return self._mesh

    @mesh.setter
    def mesh(self, mesh):
        cv.check_type('filter mesh', mesh, Mesh)
        self._mesh = mesh
        self.bins = mesh.id

    def check_bins(self, bins):
        if not len(bins) == 1:
            msg = 'Unable to add bins "{0}" to a MeshFilter since ' \
                  'only a single mesh can be used per tally'.format(bins)
            raise ValueError(msg)
        elif not isinstance(bins[0], Integral):
            msg = 'Unable to add bin "{0}" to MeshFilter since it ' \
                  'is a non-integer'.format(bins[0])
            raise ValueError(msg)
        elif bins[0] < 0:
            msg = 'Unable to add bin "{0}" to MeshFilter since it ' \
                  'is a negative integer'.format(bins[0])
            raise ValueError(msg)

    def can_merge(self, other):
        # Mesh filters cannot have more than one bin
        return False

    def get_bin_index(self, filter_bin):
        # Filter bins for a mesh are an (x,y,z) tuple. Convert (x,y,z) to a
        # single bin -- this is similar to subroutine mesh_indices_to_bin in
        # openmc/src/mesh.F90.
        if len(self.mesh.dimension) == 3:
            nx, ny, nz = self.mesh.dimension
            val = (filter_bin[0] - 1) * ny * nz + \
                  (filter_bin[1] - 1) * nz + \
                  (filter_bin[2] - 1)
        else:
            nx, ny = self.mesh.dimension
            val = (filter_bin[0] - 1) * ny + \
                  (filter_bin[1] - 1)

        return val

    def get_bin(self, bin_index):
        cv.check_type('bin_index', bin_index, Integral)
        cv.check_greater_than('bin_index', bin_index, 0, equality=True)
        cv.check_less_than('bin_index', bin_index, self.num_bins)

        # Construct 3-tuple of x,y,z cell indices for a 3D mesh
        if len(self.mesh.dimension) == 3:
            nx, ny, nz = self.mesh.dimension
            x = bin_index / (ny * nz)
            y = (bin_index - (x * ny * nz)) / nz
            z = bin_index - (x * ny * nz) - (y * nz)
            return (x, y, z)

        # Construct 2-tuple of x,y cell indices for a 2D mesh
        else:
            nx, ny = self.mesh.dimension
            x = bin_index / ny
            y = bin_index - (x * ny)
            return (x, y)

    def get_pandas_dataframe(self, data_size, distribcell_paths=True):
        # Initialize Pandas DataFrame
        import pandas as pd
        df = pd.DataFrame()

        # Initialize dictionary to build Pandas Multi-index column
        filter_dict = {}

        # Append Mesh ID as outermost index of multi-index
        mesh_key = 'mesh {0}'.format(self.mesh.id)

        # Find mesh dimensions - use 3D indices for simplicity
        if len(self.mesh.dimension) == 3:
            nx, ny, nz = self.mesh.dimension
        else:
            nx, ny = self.mesh.dimension
            nz = 1

        # Generate multi-index sub-column for x-axis
        filter_bins = np.arange(1, nx+1)
        repeat_factor = ny * nz * self.stride
        filter_bins = np.repeat(filter_bins, repeat_factor)
        tile_factor = data_size / len(filter_bins)
        filter_bins = np.tile(filter_bins, tile_factor)
        filter_dict[(mesh_key, 'x')] = filter_bins

        # Generate multi-index sub-column for y-axis
        filter_bins = np.arange(1, ny+1)
        repeat_factor = nz * self.stride
        filter_bins = np.repeat(filter_bins, repeat_factor)
        tile_factor = data_size / len(filter_bins)
        filter_bins = np.tile(filter_bins, tile_factor)
        filter_dict[(mesh_key, 'y')] = filter_bins

        # Generate multi-index sub-column for z-axis
        filter_bins = np.arange(1, nz+1)
        repeat_factor = self.stride
        filter_bins = np.repeat(filter_bins, repeat_factor)
        tile_factor = data_size / len(filter_bins)
        filter_bins = np.tile(filter_bins, tile_factor)
        filter_dict[(mesh_key, 'z')] = filter_bins

        # Initialize a Pandas DataFrame from the mesh dictionary
        df = pd.concat([df, pd.DataFrame(filter_dict)])

        return df


class EnergyFilter(Filter):
    def __gt__(self, other):
        if type(self) == type(other):
            # Compare largest/smallest energy bin edges in energy filters
            # This logic is used when merging tallies with energy filters
            return self.bins[0] >= other.bins[-1]
        else:
            return super(EnergyFilter, self).__gt__(other)

    @property
    def num_bins(self):
        return len(self.bins) - 1

    @num_bins.setter
    def num_bins(self, num_bins):
        cv.check_type('filter num_bins', num_bins, Integral)
        cv.check_greater_than('filter num_bins', num_bins, 0, equality=True)
        self._num_bins = num_bins

    def check_bins(self, bins):
        for edge in bins:
            if not isinstance(edge, Real):
                msg = 'Unable to add bin edge "{0}" to a "{1}" Filter ' \
                      'since it is a non-integer or floating point ' \
                      'value'.format(edge, self.type)
                raise ValueError(msg)
            elif edge < 0.:
                msg = 'Unable to add bin edge "{0}" to a "{1}" Filter ' \
                      'since it is a negative value'.format(edge, self.type)
                raise ValueError(msg)

        # Check that bin edges are monotonically increasing
        for index in range(1, len(bins)):
            if bins[index] < bins[index-1]:
                msg = 'Unable to add bin edges "{0}" to a "{1}" Filter ' \
                      'since they are not monotonically ' \
                      'increasing'.format(bins, self.type)
                raise ValueError(msg)

    def can_merge(self, other):
        if type(self) != type(other): return False

        if self.bins[0] == other.bins[-1]:
            # This low energy edge coincides with other's high energy edge
            return True
        elif self.bins[-1] == other.bins[0]:
            # This high energy edge coincides with other's low energy edge
            return True
        else:
            return False

    def merge(self, other):
        if not self.can_merge(other):
            msg = 'Unable to merge "{0}" with "{1}" ' \
                  'filters'.format(self.type, other.type)
            raise ValueError(msg)

        # Merge unique filter bins
        merged_bins = np.concatenate((self.bins, other.bins))
        merged_bins = np.unique(merged_bins)

        # Create a new filter with these bins
        return type(self)(sorted(merged_bins))

    def is_subset(self, other):
        if type(self) != type(other):
            return False
        elif len(self.bins) != len(other.bins):
            return False
        else:
            return np.allclose(self.bins, other.bins)

    def get_bin_index(self, filter_bin):
        # Use lower energy bound to find index for energy Filters
        deltas = np.abs(self.bins - filter_bin[1]) / filter_bin[1]
        min_delta = np.min(deltas)
        if min_delta < 1E-3:
            return deltas.argmin() - 1
        else:
            msg = 'Unable to get the bin index for Filter since "{0}" ' \
                  'is not one of the bins'.format(filter_bin)
            raise ValueError(msg)

    def get_bin(self, bin_index):
        cv.check_type('bin_index', bin_index, Integral)
        cv.check_greater_than('bin_index', bin_index, 0, equality=True)
        cv.check_less_than('bin_index', bin_index, self.num_bins)

        # Construct 2-tuple of lower, upper energies for energy(out) filters
        return (self.bins[bin_index], self.bins[bin_index+1])

    def get_pandas_dataframe(self, data_size, distribcell_paths=True):
        # Initialize Pandas DataFrame
        import pandas as pd
        df = pd.DataFrame()

        # Extract the lower and upper energy bounds, then repeat and tile
        # them as necessary to account for other filters.
        lo_bins = np.repeat(self.bins[:-1], self.stride)
        hi_bins = np.repeat(self.bins[1:], self.stride)
        tile_factor = data_size / len(lo_bins)
        lo_bins = np.tile(lo_bins, tile_factor)
        hi_bins = np.tile(hi_bins, tile_factor)

        # Add the new energy columns to the DataFrame.
        df.loc[:, self.short_name.lower() + ' low [MeV]'] = lo_bins
        df.loc[:, self.short_name.lower() + ' high [MeV]'] = hi_bins

        return df


class EnergyoutFilter(EnergyFilter): pass


class DistribcellFilter(Filter):
    """A constraint on transport events that can score a tally

    Parameters
    ----------
    bins : Integral or Iterable of Integral or Iterable of Real
        The bins for the filter. This takes on different meaning for different
        filters. See the OpenMC online documentation for more details.

    Attributes
    ----------
    bins : Integral or Iterable of Integral or Iterable of Real
        The bins for the filter
    num_bins : Integral
        The number of filter bins
    stride : Integral
        The number of filter, nuclide and score bins within each of this
        filter's bins.
    distribcell_paths : list of str
        The paths traversed through the CSG tree to reach each distribcell
        instance (for 'distribcell' filters only)

    """

    def __init__(self, bins):
        self._distribcell_paths = None
        super(DistribcellFilter, self).__init__(bins)

    @classmethod
    def from_hdf5(cls, group):
        if group['type'].value.decode() != cls.short_name.lower():
            raise ValueError("Expected HDF5 data for filter type '"
                             + cls.short_name.lower() + "' but got '"
                             + group['type'].value.decode() + " instead")

        out = cls(group['bins'].value)
        out.num_bins = group['n_bins'].value

        if 'paths' in group:
            out.distribcell_paths = [str(path.decode()) for path in
                                     group['paths'].value]

        return out

    @property
    def distribcell_paths(self):
        return self._distribcell_paths

    @distribcell_paths.setter
    def distribcell_paths(self, distribcell_paths):
        cv.check_iterable_type('distribcell_paths', distribcell_paths, str)
        self._distribcell_paths = distribcell_paths

    def check_bins(self, bins):
        if not len(bins) == 1:
            msg = 'Unable to add bins "{0}" to a DistribcellFilter since ' \
                  'only a single distribcell can be used per tally'.format(bins)
            raise ValueError(msg)

        cv.check_iterable_type('filter bins', bins, Integral)
        for edge in bins:
            cv.check_greater_than('filter bin', edge, 0, equality=True)

    def can_merge(self, other):
        # Distribcell filters cannot have more than one bin
        return False

    def get_bin_index(self, filter_bin):
        # Filter bins for distribcells are indices of each unique placement of
        # the Cell in the Geometry (consecutive integers starting at 0).
        return filter_bin

    def get_pandas_dataframe(self, data_size, distribcell_paths=True):
        # Initialize Pandas DataFrame
        import pandas as pd
        df = pd.DataFrame()

        level_df = None

        # Create Pandas Multi-index columns for each level in CSG tree
        if distribcell_paths:

            # Distribcell paths require linked metadata from the Summary
            if self.distribcell_paths is None:
                msg = 'Unable to construct distribcell paths since ' \
                      'the Summary is not linked to the StatePoint'
                raise ValueError(msg)

            # Make copy of array of distribcell paths to use in
            # Pandas Multi-index column construction
            distribcell_paths = copy.deepcopy(self.distribcell_paths)
            num_offsets = len(distribcell_paths)

            # Loop over CSG levels in the distribcell paths
            level_counter = 0
            levels_remain = True
            while levels_remain:

                # Use level key as first index in Pandas Multi-index column
                level_counter += 1
                level_key = 'level {}'.format(level_counter)

                # Use the first distribcell path to determine if level
                # is a universe/cell or lattice level
                first_path = distribcell_paths[0]
                next_index = first_path.index('-')
                level = first_path[:next_index]

                # Trim universe/lattice info from path
                first_path = first_path[next_index+2:]

                # Create a dictionary for this level for Pandas Multi-index
                level_dict = OrderedDict()

                # This level is a lattice (e.g., ID(x,y,z))
                if '(' in level:
                    level_type = 'lattice'

                    # Initialize prefix Multi-index keys
                    lat_id_key = (level_key, 'lat', 'id')
                    lat_x_key = (level_key, 'lat', 'x')
                    lat_y_key = (level_key, 'lat', 'y')
                    lat_z_key = (level_key, 'lat', 'z')

                    # Allocate NumPy arrays for each CSG level and
                    # each Multi-index column in the DataFrame
                    level_dict[lat_id_key] = np.empty(num_offsets)
                    level_dict[lat_x_key] = np.empty(num_offsets)
                    level_dict[lat_y_key] = np.empty(num_offsets)
                    level_dict[lat_z_key] = np.empty(num_offsets)

                # This level is a universe / cell (e.g., ID->ID)
                else:
                    level_type = 'universe'

                    # Initialize prefix Multi-index keys
                    univ_key = (level_key, 'univ', 'id')
                    cell_key = (level_key, 'cell', 'id')

                    # Allocate NumPy arrays for each CSG level and
                    # each Multi-index column in the DataFrame
                    level_dict[univ_key] = np.empty(num_offsets)
                    level_dict[cell_key] = np.empty(num_offsets)

                    # Determine any levels remain in path
                    if '-' not in first_path:
                        levels_remain = False

                # Populate Multi-index arrays with all distribcell paths
                for i, path in enumerate(distribcell_paths):

                    if level_type == 'lattice':
                        # Extract lattice ID, indices from path
                        next_index = path.index('-')
                        lat_id_indices = path[:next_index]

                        # Trim lattice info from distribcell path
                        distribcell_paths[i] = path[next_index+2:]

                        # Extract the lattice cell indices from the path
                        i1 = lat_id_indices.index('(')
                        i2 = lat_id_indices.index(')')
                        i3 = lat_id_indices[i1+1:i2]

                        # Assign entry to Lattice Multi-index column
                        level_dict[lat_id_key][i] = path[:i1]
                        level_dict[lat_x_key][i] = int(i3.split(',')[0]) - 1
                        level_dict[lat_y_key][i] = int(i3.split(',')[1]) - 1
                        level_dict[lat_z_key][i] = int(i3.split(',')[2]) - 1

                    else:
                        # Extract universe ID from path
                        next_index = path.index('-')
                        universe_id = int(path[:next_index])

                        # Trim universe info from distribcell path
                        path = path[next_index+2:]

                        # Extract cell ID from path
                        if '-' in path:
                            next_index = path.index('-')
                            cell_id = int(path[:next_index])
                            distribcell_paths[i] = path[next_index+2:]
                        else:
                            cell_id = int(path)
                            distribcell_paths[i] = ''

                        # Assign entry to Universe, Cell Multi-index columns
                        level_dict[univ_key][i] = universe_id
                        level_dict[cell_key][i] = cell_id

                # Tile the Multi-index columns
                for level_key, level_bins in level_dict.items():
                    level_bins = np.repeat(level_bins, self.stride)
                    tile_factor = data_size / len(level_bins)
                    level_bins = np.tile(level_bins, tile_factor)
                    level_dict[level_key] = level_bins

                # Initialize a Pandas DataFrame from the level dictionary
                if level_df is None:
                    level_df = pd.DataFrame(level_dict)
                else:
                    level_df = pd.concat([level_df, pd.DataFrame(level_dict)],
                                         axis=1)

        # Create DataFrame column for distribcell instance IDs
        # NOTE: This is performed regardless of whether the user
        # requests Summary geometric information
        filter_bins = np.arange(self.num_bins)
        filter_bins = np.repeat(filter_bins, self.stride)
        tile_factor = data_size / len(filter_bins)
        filter_bins = np.tile(filter_bins, tile_factor)
        df = pd.DataFrame({self.short_name.lower() : filter_bins})

        # If OpenCG level info DataFrame was created, concatenate
        # with DataFrame of distribcell instance IDs
        if level_df is not None:
            level_df = level_df.dropna(axis=1, how='all')
            level_df = level_df.astype(np.int)
            df = pd.concat([level_df, df], axis=1)

        return df


class MuFilter(Filter): pass


class PolarFilter(Filter):
    def get_pandas_dataframe(self, data_size, distribcell_paths=True):
        # Initialize Pandas DataFrame
        import pandas as pd
        df = pd.DataFrame()

        # Extract the lower and upper angle bounds, then repeat and tile
        # them as necessary to account for other filters.
        lo_bins = np.repeat(self.bins[:-1], self.stride)
        hi_bins = np.repeat(self.bins[1:], self.stride)
        tile_factor = data_size / len(lo_bins)
        lo_bins = np.tile(lo_bins, tile_factor)
        hi_bins = np.tile(hi_bins, tile_factor)

        # Add the new angle columns to the DataFrame.
        df.loc[:, self.type + ' low'] = lo_bins

        return df


class AzimuthalFilter(Filter):
    def get_pandas_dataframe(self, data_size, distribcell_paths=True):
        # Initialize Pandas DataFrame
        import pandas as pd
        df = pd.DataFrame()

        # Extract the lower and upper angle bounds, then repeat and tile
        # them as necessary to account for other filters.
        lo_bins = np.repeat(self.bins[:-1], self.stride)
        hi_bins = np.repeat(self.bins[1:], self.stride)
        tile_factor = data_size / len(lo_bins)
        lo_bins = np.tile(lo_bins, tile_factor)
        hi_bins = np.tile(hi_bins, tile_factor)

        # Add the new angle columns to the DataFrame.
        df.loc[:, self.type + ' low'] = lo_bins

        return df


class DelayedGroupFilter(IntegralFilter): pass
