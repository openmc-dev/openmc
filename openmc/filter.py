from collections import Iterable
import copy
from numbers import Real, Integral
import sys

import numpy as np

from openmc import Mesh
from openmc.summary import Summary
import openmc.checkvalue as cv


if sys.version_info[0] >= 3:
    basestring = str


_FILTER_TYPES = ['universe', 'material', 'cell', 'cellborn', 'surface',
                 'mesh', 'energy', 'energyout', 'mu', 'polar', 'azimuthal',
                 'distribcell', 'delayedgroup']

class Filter(object):
    """A filter used to constrain a tally to a specific criterion, e.g. only
    tally events when the particle is in a certain cell and energy range.

    Parameters
    ----------
    type : str
        The type of the tally filter. Acceptable values are "universe",
        "material", "cell", "cellborn", "surface", "mesh", "energy",
        "energyout", and "distribcell".
    bins : Integral or Iterable of Integral or Iterable of Real
        The bins for the filter. This takes on different meaning for different
        filters. See the OpenMC online documentation for more details.

    Attributes
    ----------
    type : str
        The type of the tally filter
    bins : Integral or Iterable of Real
        The bins for the filter
    num_bins : Integral
        The number of filter bins
    mesh : Mesh or None
        A Mesh object for 'mesh' type filters.
    stride : Integral
        The number of filter, nuclide and score bins within each of this
        filter's bins.
    distribcell_paths : list of str
        The paths traversed through the CSG tree to reach each distribcell
        instance (for 'distribcell' filters only)

    """

    # Initialize Filter class attributes
    def __init__(self, type=None, bins=None):

        self._type = None
        self._num_bins = 0
        self._bins = None
        self._mesh = None
        self._stride = None
        self._distribcell_paths = None

        if type is not None:
            self.type = type
        if bins is not None:
            self.bins = bins

    def __eq__(self, other):
        if not isinstance(other, Filter):
            return False
        elif self.type != other.type:
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
        if self.type != other.type:
            if self.type in _FILTER_TYPES and other.type in _FILTER_TYPES:
                delta = _FILTER_TYPES.index(self.type) - \
                        _FILTER_TYPES.index(other.type)
                return delta > 0
            else:
                return False
        else:
            # Compare largest/smallest energy bin edges in energy filters
            # This logic is used when merging tallies with energy filters
            if 'energy' in self.type and 'energy' in other.type:
                return self.bins[0] >= other.bins[-1]
            else:
                return max(self.bins) > max(other.bins)

    def __lt__(self, other):
        return not self > other

    def __hash__(self):
        return hash(repr(self))

    def __deepcopy__(self, memo):
        existing = memo.get(id(self))

        # If this is the first time we have tried to copy this object, create a copy
        if existing is None:
            clone = type(self).__new__(type(self))
            clone._type = self.type
            clone._bins = copy.deepcopy(self.bins, memo)
            clone._num_bins = self.num_bins
            clone._mesh = copy.deepcopy(self.mesh, memo)
            clone._stride = self.stride
            clone._distribcell_paths = copy.deepcopy(self.distribcell_paths)

            memo[id(self)] = clone

            return clone

        # If this object has been copied before, return the first copy made
        else:
            return existing

    def __repr__(self):
        string = 'Filter\n'
        string += '{0: <16}{1}{2}\n'.format('\tType', '=\t', self.type)
        string += '{0: <16}{1}{2}\n'.format('\tBins', '=\t', self.bins)
        return string

    @property
    def type(self):
        return self._type

    @property
    def bins(self):
        return self._bins

    @property
    def num_bins(self):
        if self.bins is None:
            return 0
        elif self.type in ['energy', 'energyout']:
            return len(self.bins) - 1
        elif self.type in ['cell', 'cellborn', 'surface', 'universe', 'material']:
            return len(self.bins)
        else:
            return self._num_bins

    @property
    def mesh(self):
        return self._mesh

    @property
    def stride(self):
        return self._stride

    @property
    def distribcell_paths(self):
        return self._distribcell_paths

    @type.setter
    def type(self, type):
        if type is None:
            self._type = type
        elif type not in _FILTER_TYPES:
            msg = 'Unable to set Filter type to "{0}" since it is not one ' \
                  'of the supported types'.format(type)
            raise ValueError(msg)

        self._type = type

    @bins.setter
    def bins(self, bins):
        if self.type is None:
            msg = 'Unable to set bins for Filter to "{0}" since ' \
                  'the Filter type has not yet been set'.format(bins)
            raise ValueError(msg)

        # If the bin edge is a single value, it is a Cell, Material, etc. ID
        if not isinstance(bins, Iterable):
            bins = [bins]

        # If the bins are in a collection, convert it to a list
        else:
            bins = list(bins)

        if self.type in ['cell', 'cellborn', 'surface', 'material',
                         'universe', 'distribcell', 'delayedgroup']:
            cv.check_iterable_type('filter bins', bins, Integral)
            for edge in bins:
                cv.check_greater_than('filter bin', edge, 0, equality=True)

        elif self.type in ['energy', 'energyout']:
            for edge in bins:
                if not cv._isinstance(edge, Real):
                    msg = 'Unable to add bin edge "{0}" to a "{1}" Filter ' \
                          'since it is a non-integer or floating point ' \
                          'value'.format(edge, self.type)
                    raise ValueError(msg)
                elif edge < 0.:
                    msg = 'Unable to add bin edge "{0}" to a "{1}" Filter ' \
                          'since it is a negative value'.format(edge, self.type)
                    raise ValueError(msg)

            # Check that bin edges are monotonically increasing
            for index in range(len(bins)):
                if index > 0 and bins[index] < bins[index-1]:
                    msg = 'Unable to add bin edges "{0}" to a "{1}" Filter ' \
                          'since they are not monotonically ' \
                          'increasing'.format(bins, self.type)
                    raise ValueError(msg)

        # mesh filters
        elif self.type == 'mesh':
            if not len(bins) == 1:
                msg = 'Unable to add bins "{0}" to a mesh Filter since ' \
                      'only a single mesh can be used per tally'.format(bins)
                raise ValueError(msg)
            elif not cv._isinstance(bins[0], Integral):
                msg = 'Unable to add bin "{0}" to mesh Filter since it ' \
                       'is a non-integer'.format(bins[0])
                raise ValueError(msg)
            elif bins[0] < 0:
                msg = 'Unable to add bin "{0}" to mesh Filter since it ' \
                       'is a negative integer'.format(bins[0])
                raise ValueError(msg)

        # If all error checks passed, add bin edges
        self._bins = np.array(bins)

    @num_bins.setter
    def num_bins(self, num_bins):
        cv.check_type('filter num_bins', num_bins, Integral)
        cv.check_greater_than('filter num_bins', num_bins, 0, equality=True)
        self._num_bins = num_bins

    @mesh.setter
    def mesh(self, mesh):
        cv.check_type('filter mesh', mesh, Mesh)

        self._mesh = mesh
        self.type = 'mesh'
        self.bins = self.mesh.id

    @stride.setter
    def stride(self, stride):
        cv.check_type('filter stride', stride, Integral)
        if stride < 0:
            msg = 'Unable to set stride "{0}" for a "{1}" Filter since it ' \
                  'is a negative value'.format(stride, self.type)
            raise ValueError(msg)

        self._stride = stride

    @distribcell_paths.setter
    def distribcell_paths(self, distribcell_paths):
        cv.check_iterable_type('distribcell_paths', distribcell_paths, str)
        self._distribcell_paths = distribcell_paths

    def can_merge(self, other):
        """Determine if filter can be merged with another.

        Parameters
        ----------
        other : Filter
            Filter to compare with

        Returns
        -------
        bool
            Whether the filter can be merged

        """

        if not isinstance(other, Filter):
            return False

        # Filters must be of the same type
        if self.type != other.type:
            return False

        # Distribcell filters cannot have more than one bin
        if self.type == 'distribcell':
            return False

        # Mesh filters cannot have more than one bin
        elif self.type == 'mesh':
            return False

        # Different energy bins structures must be mutually exclusive and
        # share only one shared bin edge at the minimum or maximum energy
        elif 'energy' in self.type:
            # This low energy edge coincides with other's high energy edge
            if self.bins[0] == other.bins[-1]:
                return True
            # This high energy edge coincides with other's low energy edge
            elif self.bins[-1] == other.bins[0]:
                return True
            else:
                return False

        else:
            return True

    def merge(self, other):
        """Merge this filter with another.

        Parameters
        ----------
        other : Filter
            Filter to merge with

        Returns
        -------
        merged_filter : Filter
            Filter resulting from the merge

        """

        if not self.can_merge(other):
            msg = 'Unable to merge "{0}" with "{1}" ' \
                  'filters'.format(self.type, other.type)
            raise ValueError(msg)

        # Create deep copy of filter to return as merged filter
        merged_filter = copy.deepcopy(self)

        # Merge unique filter bins
        merged_bins = np.concatenate((self.bins, other.bins))
        merged_bins = np.unique(merged_bins)

        # Sort energy bin edges
        if 'energy' in self.type:
            merged_bins = sorted(merged_bins)

        # Assign merged bins to merged filter
        merged_filter.bins = list(merged_bins)

        # Count bins in the merged filter
        if 'energy' in merged_filter.type:
            merged_filter.num_bins = len(merged_bins) - 1
        else:
            merged_filter.num_bins = len(merged_bins)

        return merged_filter

    def is_subset(self, other):
        """Determine if another filter is a subset of this filter.

        If all of the bins in the other filter are included as bins in this
        filter, then it is a subset of this filter.

        Parameters
        ----------
        other : Filter
            The filter to query as a subset of this filter

        Returns
        -------
        bool
            Whether or not the other filter is a subset of this filter

        """

        if not isinstance(other, Filter):
            return False
        elif self.type != other.type:
            return False
        elif self.type in ['energy', 'energyout']:
            if len(self.bins) != len(other.bins):
                return False
            else:
                return np.allclose(self.bins, other.bins)

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
            3-tuple for 'mesh' filters corresponding to the mesh cell
            interest.

        Returns
        -------
        filter_index : Integral
             The index in the Tally data array for this filter bin.

        See also
        --------
        Filter.get_bin()

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
                deltas = np.abs(self.bins - filter_bin[1]) / filter_bin[1]
                min_delta = np.min(deltas)
                if min_delta < 1E-3:
                    filter_index = deltas.argmin() - 1
                else:
                    raise ValueError

            # Filter bins for distribcells are "IDs" of each unique placement
            # of the Cell in the Geometry (integers starting at 0)
            elif self.type == 'distribcell':
                filter_index = filter_bin

            # Use ID for all other Filters (e.g., material, cell, etc.)
            else:
                val = np.where(self.bins == filter_bin)[0][0]
                filter_index = val

        except ValueError:
            msg = 'Unable to get the bin index for Filter since "{0}" ' \
                  'is not one of the bins'.format(filter_bin)
            raise ValueError(msg)

        return filter_index

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

        if self.type == 'mesh':

            # Construct 3-tuple of x,y,z cell indices for a 3D mesh
            if len(self.mesh.dimension) == 3:
                nx, ny, nz = self.mesh.dimension
                x = bin_index / (ny * nz)
                y = (bin_index - (x * ny * nz)) / nz
                z = bin_index - (x * ny * nz) - (y * nz)
                filter_bin = (x, y, z)

            # Construct 2-tuple of x,y cell indices for a 2D mesh
            else:
                nx, ny = self.mesh.dimension
                x = bin_index / ny
                y = bin_index - (x * ny)
                filter_bin = (x, y)

        # Construct 2-tuple of lower, upper energies for energy(out) filters
        elif self.type in ['energy', 'energyout']:
            filter_bin = (self.bins[bin_index], self.bins[bin_index+1])
        # Construct 1-tuple of with the cell ID for distribcell filters
        elif self.type == 'distribcell':
            filter_bin = (self.bins[0],)
        # Construct 1-tuple with domain ID (e.g., material) for other filters
        else:
            filter_bin = (self.bins[bin_index],)

        return filter_bin

    def get_pandas_dataframe(self, data_size, summary=None):
        """Builds a Pandas DataFrame for the Filter's bins.

        This method constructs a Pandas DataFrame object for the filter with
        columns annotated by filter bin information. This is a helper method
        for the Tally.get_pandas_dataframe(...) method.

        This capability has been tested for Pandas >=0.13.1. However, it is
        recommended to use v0.16 or newer versions of Pandas since this method
        uses Pandas' Multi-index functionality.

        Parameters
        ----------
        data_size : Integral
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
            filter's bins. The number of rows in the DataFrame is the same as
            the total number of bins in the corresponding tally, with the filter
            bin appropriately tiled to map to the corresponding tally bins.

            For 'cell', 'cellborn', 'surface', 'material', and 'universe'
            filters, the DataFrame includes a single column with the cell,
            surface, material or universe ID corresponding to each filter bin.

            For 'distribcell' filters, the DataFrame either includes:

            1. a single column with the cell instance IDs (without summary info)
            2. separate columns for the cell IDs, universe IDs, and lattice IDs
               and x,y,z cell indices corresponding to each (with summary info).

            For 'energy' and 'energyout' filters, the DataFrame includes one
            column for the lower energy bound and one column for the upper
            energy bound for each filter bin.

            For 'mesh' filters, the DataFrame includes three columns for the
            x,y,z mesh cell indices corresponding to each filter bin.

        Raises
        ------
        ImportError
            When Pandas is not installed, or summary info is requested but
            OpenCG is not installed.

        See also
        --------
        Tally.get_pandas_dataframe(), CrossFilter.get_pandas_dataframe()

        """

        # Initialize Pandas DataFrame
        import pandas as pd
        df = pd.DataFrame()

        # mesh filters
        if self.type == 'mesh':

            # Initialize dictionary to build Pandas Multi-index column
            filter_dict = {}

            # Append Mesh ID as outermost index of mult-index
            mesh_key = 'mesh {0}'.format(self.mesh.id)

            # Find mesh dimensions - use 3D indices for simplicity
            if (len(self.mesh.dimension) == 3):
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

        # distribcell filters
        elif self.type == 'distribcell':
            level_df = None

            if isinstance(summary, Summary):
                # Attempt to import the OpenCG package
                try:
                    import opencg
                except ImportError:
                    msg = 'The OpenCG package must be installed ' \
                          'to use a Summary for distribcell dataframes'
                    raise ImportError(msg)

                # Extract the OpenCG geometry from the Summary
                opencg_geometry = summary.opencg_geometry
                openmc_geometry = summary.openmc_geometry

                # Use OpenCG to compute the number of regions
                opencg_geometry.initialize_cell_offsets()
                num_regions = opencg_geometry.num_regions

                # Initialize a dictionary mapping OpenMC distribcell
                # offsets to OpenCG LocalCoords linked lists
                offsets_to_coords = {}

                for offset, path in enumerate(self.distribcell_paths):
                    region = opencg_geometry.get_region_from_path(path)
                    coords = opencg_geometry.find_region(region)
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
                    level_dict[univ_key][:] = np.NAN
                    level_dict[cell_key][:] = np.NAN
                    level_dict[lat_id_key][:] = np.NAN
                    level_dict[lat_x_key][:] = np.NAN
                    level_dict[lat_y_key][:] = np.NAN
                    level_dict[lat_z_key][:] = np.NAN

                    # Iterate over all regions (distribcell instances)
                    for offset in range(num_offsets):
                        coords = offsets_to_coords[offset]

                        # If entire LocalCoords has been unraveled into
                        # Multi-index columns already, continue
                        if coords is None:
                            continue

                        # Assign entry to Universe Multi-index column
                        if coords._type == 'universe':
                            level_dict[univ_key][offset] = coords._universe._id
                            level_dict[cell_key][offset] = coords._cell._id

                        # Assign entry to Lattice Multi-index column
                        else:
                            # Reverse y index per lattice ordering in OpenCG
                            level_dict[lat_id_key][offset] = coords._lattice._id
                            level_dict[lat_x_key][offset] = coords._lat_x
                            level_dict[lat_y_key][offset] = \
                                coords._lattice.dimension[1] - coords._lat_y - 1
                            level_dict[lat_z_key][offset] = coords._lat_z

                        # Move to next node in LocalCoords linked list
                        if coords._next is None:
                            offsets_to_coords[offset] = None
                        else:
                            offsets_to_coords[offset] = coords._next
                            levels_remain = True

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
                        level_df = pd.concat([level_df, pd.DataFrame(level_dict)], axis=1)

            # Create DataFrame column for distribcell instances IDs
            # NOTE: This is performed regardless of whether the user
            # requests Summary geometric information
            filter_bins = np.arange(self.num_bins)
            filter_bins = np.repeat(filter_bins, self.stride)
            tile_factor = data_size / len(filter_bins)
            filter_bins = np.tile(filter_bins, tile_factor)
            df = pd.DataFrame({self.type : filter_bins})

            # If OpenCG level info DataFrame was created, concatenate
            # with DataFrame of distribcell instance IDs
            if level_df is not None:
                level_df = level_df.dropna(axis=1, how='all')
                level_df = level_df.astype(np.int)
                df = pd.concat([level_df, df], axis=1)

        # energy, energyout filters
        elif 'energy' in self.type:
            # Extract the lower and upper energy bounds, then repeat and tile
            # them as necessary to account for other filters.
            lo_bins = np.repeat(self.bins[:-1], self.stride)
            hi_bins = np.repeat(self.bins[1:], self.stride)
            tile_factor = data_size / len(lo_bins)
            lo_bins = np.tile(lo_bins, tile_factor)
            hi_bins = np.tile(hi_bins, tile_factor)

            # Add the new energy columns to the DataFrame.
            df.loc[:, self.type + ' low [MeV]'] = lo_bins
            df.loc[:, self.type + ' high [MeV]'] = hi_bins

        elif self.type in ('azimuthal', 'polar'):
            # Extract the lower and upper angle bounds, then repeat and tile
            # them as necessary to account for other filters.
            lo_bins = np.repeat(self.bins[:-1], self.stride)
            hi_bins = np.repeat(self.bins[1:], self.stride)
            tile_factor = data_size / len(lo_bins)
            lo_bins = np.tile(lo_bins, tile_factor)
            hi_bins = np.tile(hi_bins, tile_factor)

            # Add the new angle columns to the DataFrame.
            df.loc[:, self.type + ' low'] = lo_bins
            df.loc[:, self.type + ' high'] = hi_bins

        # universe, material, surface, cell, and cellborn filters
        else:
            filter_bins = np.repeat(self.bins, self.stride)
            tile_factor = data_size / len(filter_bins)
            filter_bins = np.tile(filter_bins, tile_factor)
            filter_bins = filter_bins
            df = pd.concat([df, pd.DataFrame({self.type : filter_bins})])

        return df
