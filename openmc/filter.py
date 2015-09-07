from collections import Iterable
import copy
from numbers import Real, Integral

import numpy as np

from openmc import Mesh
from openmc.summary import Summary
from openmc.constants import *
import openmc.checkvalue as cv


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
        if self.type != filter2.type:
            return False

        # Check number of bins
        elif len(self.bins) != len(filter2.bins):
            return False

        # Check bin edges
        elif not np.allclose(self.bins, filter2.bins):
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
            clone._type = self.type
            clone._bins = copy.deepcopy(self.bins, memo)
            clone._num_bins = self.num_bins
            clone._mesh = copy.deepcopy(self.mesh, memo)
            clone._offset = self.offset
            clone._stride = self.stride

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

        if self.type in ['cell', 'cellborn', 'surface', 'material',
                          'universe', 'distribcell']:
            cv.check_iterable_type('filter bins', bins, Integral)
            for edge in bins:
                cv.check_greater_than('filter bin', edge, 0, equality=True)

        elif self._type in ['energy', 'energyout']:
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
            for index in range(len(bins)):
                if index > 0 and bins[index] < bins[index-1]:
                    msg = 'Unable to add bin edges "{0}" to a "{1}" Filter ' \
                          'since they are not monotonically ' \
                          'increasing'.format(bins, self.type)
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
        self._bins = np.array(bins)

    # FIXME
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

    @offset.setter
    def offset(self, offset):
        cv.check_type('filter offset', offset, Integral)
        self._offset = offset

    @stride.setter
    def stride(self, stride):
        cv.check_type('filter stride', stride, Integral)
        if stride < 0:
            msg = 'Unable to set stride "{0}" for a "{1}" Filter since it ' \
                  'is a negative value'.format(stride, self.type)
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
            msg = 'Unable to merge "{0}" with "{1}" ' \
                  'filters'.format(self.type, filter.type)
            raise ValueError(msg)

        # Create deep copy of filter to return as merged filter
        merged_filter = copy.deepcopy(self)

        # Merge unique filter bins
        merged_bins = list(set(self.bins + filter.bins))
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
                val = np.where(self.bins == filter_bin[0])[0][0]
                filter_index = val

            # Filter bins for distribcell are the "IDs" of each unique placement
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

    def get_pandas_dataframe(self, data_size, summary=None):
        """Builds a Pandas DataFrame for the Filter's bins.

        This method constructs a Pandas DataFrame object for the Filter with
        columns annotated by filter bin information. This is a helper method
        for the Tally.get_pandas_dataframe(...) routine.

        This capability has been tested for Pandas >=0.13.1. However, it is
        recommended to use v0.16 or newer versions of Pandas since this method
        uses the Multi-index Pandas feature.


        Parameters
        ----------
        data_size : Integral
            The total number of bins in the tally corresponding to this filter

        summary : None or Summary
            An optional Summary object to be used to construct columns for
            distribcell tally filters (default is None). The geometric
            information in the Summary object is embedded into a Multi-index
            column with a geometric "path" to each distribcell intance.
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

            For 'mesh' filters, the DataFrame includes three columns for the
            x,y,z mesh cell indices corresponding to each filter bin.

            For 'energy' and 'energyout' filters, the DataFrame include a single
            column with each element comprising a string with the lower, upper
            energy bounds for each filter bin.

            For 'distribcell' filters, the DataFrame either includes:
            1) a single column with the cell instance IDs (without summary info)
            2) separate columns for the cell IDs, universe IDs, and lattice IDs
               and x,y,z cell indices corresponding to each (with summary info)

        Raises
        ------
        ImportError
            When Pandas cannot is not installed, or summary info is requested
            but OpenCG is not installed.

        See also
        --------
        Tally.get_pandas_dataframe(), CrossFilter.get_pandas_dataframe()

        """

        # Attempt to import the pandas package
        try:
            import pandas as pd
        except ImportError:
            msg = 'The pandas Python package must be installed on your system'
            raise ImportError(msg)

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
                    if cell_id == self.bins[0]:
                        offset = openmc_geometry.get_offset(path, self.offset)
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
                            level_dict[lat_id_key][offset] = coords._lattice._id
                            level_dict[lat_x_key][offset] = coords._lat_x
                            level_dict[lat_y_key][offset] = coords._lat_y
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
            filter_bins = filter_bins
            if level_df is None:
                df = pd.DataFrame({self.type :filter_bins})
            else:
                level_df = level_df.dropna(axis=1, how='all')
                level_df = level_df.astype(np.int)
                df = pd.concat([level_df, pd.DataFrame({self.type :filter_bins})], axis=1)

        # energy, energyout filters
        elif 'energy' in self.type:
            bins = self.bins
            num_bins = self.num_bins

            # Create strings for
            template = '({0:.1e} - {1:.1e})'
            filter_bins = []
            for i in range(num_bins):
                filter_bins.append(template.format(bins[i], bins[i+1]))

            # Tile the energy bins into a DataFrame column
            filter_bins = np.repeat(filter_bins, self.stride)
            tile_factor = data_size / len(filter_bins)
            filter_bins = np.tile(filter_bins, tile_factor)
            filter_bins = filter_bins
            df = pd.concat([df, pd.DataFrame({self.type + ' [MeV]' : filter_bins})])

        # universe, material, surface, cell, and cellborn filters
        else:
            filter_bins = np.repeat(self.bins, self.stride)
            tile_factor = data_size / len(filter_bins)
            filter_bins = np.tile(filter_bins, tile_factor)
            filter_bins = filter_bins
            df = pd.concat([df, pd.DataFrame({self.type :filter_bins})])

        df = df.astype(np.str)
        return df

    def __repr__(self):
        string = 'Filter\n'
        string += '{0: <16}{1}{2}\n'.format('\tType', '=\t', self.type)
        string += '{0: <16}{1}{2}\n'.format('\tBins', '=\t', self.bins)
        string += '{0: <16}{1}{2}\n'.format('\tOffset', '=\t', self.offset)
        return string
