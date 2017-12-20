from __future__ import division
from abc import ABCMeta
from collections import Iterable, OrderedDict
import copy
from functools import reduce
import hashlib
from numbers import Real, Integral
import operator
from xml.etree import ElementTree as ET

from six import add_metaclass
import numpy as np
import pandas as pd

import openmc
import openmc.checkvalue as cv
from .cell import Cell
from .material import Material
from .mixin import IDManagerMixin
from .universe import Universe


_FILTER_TYPES = ['universe', 'material', 'cell', 'cellborn', 'surface',
                 'mesh', 'energy', 'energyout', 'mu', 'polar', 'azimuthal',
                 'distribcell', 'delayedgroup', 'energyfunction', 'cellfrom']

_CURRENT_NAMES = {1:  'x-min out', 2:  'x-min in',
                  3:  'x-max out', 4:  'x-max in',
                  5:  'y-min out', 6:  'y-min in',
                  7:  'y-max out', 8:  'y-max in',
                  9:  'z-min out', 10: 'z-min in',
                  11: 'z-max out', 12: 'z-max in'}


class FilterMeta(ABCMeta):
    def __new__(cls, name, bases, namespace, **kwargs):
        # Check the class name.
        if not name.endswith('Filter'):
            raise ValueError("All filter class names must end with 'Filter'")

        # Create a 'short_name' attribute that removes the 'Filter' suffix.
        namespace['short_name'] = name[:-6]

        # Subclass methods can sort of inherit the docstring of parent class
        # methods.  If a function is defined without a docstring, most (all?)
        # Python interpreters will search through the parent classes to see if
        # there is a docstring for a function with the same name, and they will
        # use that docstring.  However, Sphinx does not have that functionality.
        # This chunk of code handles this docstring inheritance manually so that
        # the autodocumentation will pick it up.
        if name != 'Filter':
            # Look for newly-defined functions that were also in Filter.
            for func_name in namespace:
                if func_name in Filter.__dict__:
                    # Inherit the docstring from Filter if not defined.
                    if isinstance(namespace[func_name],
                                  (classmethod, staticmethod)):
                        new_doc = namespace[func_name].__func__.__doc__
                        old_doc = Filter.__dict__[func_name].__func__.__doc__
                        if new_doc is None and old_doc is not None:
                            namespace[func_name].__func__.__doc__ = old_doc
                    else:
                        new_doc = namespace[func_name].__doc__
                        old_doc = Filter.__dict__[func_name].__doc__
                        if new_doc is None and old_doc is not None:
                            namespace[func_name].__doc__ = old_doc

        # Make the class.
        return super(FilterMeta, cls).__new__(cls, name, bases, namespace,
                                              **kwargs)


@add_metaclass(FilterMeta)
class Filter(IDManagerMixin):
    """Tally modifier that describes phase-space and other characteristics.

    Parameters
    ----------
    bins : Integral or Iterable of Integral or Iterable of Real
        The bins for the filter. This takes on different meaning for different
        filters. See the docstrings for sublcasses of this filter or the online
        documentation for more details.
    filter_id : int
        Unique identifier for the filter

    Attributes
    ----------
    bins : Integral or Iterable of Integral or Iterable of Real
        The bins for the filter
    id : int
        Unique identifier for the filter
    num_bins : Integral
        The number of filter bins

    """

    next_id = 1
    used_ids = set()

    def __init__(self, bins, filter_id=None):
        self.bins = bins
        self.id = filter_id

    def __eq__(self, other):
        if type(self) is not type(other):
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
        if type(self) is not type(other):
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
        string = type(self).__name__ + '\n'
        string += '{: <16}=\t{}\n'.format('\tBins', self.bins)
        return hash(string)

    def __repr__(self):
        string = type(self).__name__ + '\n'
        string += '{: <16}=\t{}\n'.format('\tBins', self.bins)
        string += '{: <16}=\t{}\n'.format('\tID', self.id)
        return string

    @classmethod
    def _recursive_subclasses(cls):
        """Return all subclasses and their subclasses, etc."""
        all_subclasses = []

        for subclass in cls.__subclasses__():
            all_subclasses.append(subclass)
            all_subclasses.extend(subclass._recursive_subclasses())

        return all_subclasses

    @classmethod
    def from_hdf5(cls, group, **kwargs):
        """Construct a new Filter instance from HDF5 data.

        Parameters
        ----------
        group : h5py.Group
            HDF5 group to read from

        Keyword arguments
        -----------------
        meshes : dict
            Dictionary mapping integer IDs to openmc.Mesh objects.  Only used
            for openmc.MeshFilter objects.

        """

        filter_id = int(group.name.split('/')[-1].lstrip('filter '))

        # If the HDF5 'type' variable matches this class's short_name, then
        # there is no overriden from_hdf5 method.  Pass the bins to __init__.
        if group['type'].value.decode() == cls.short_name.lower():
            out = cls(group['bins'].value, filter_id=filter_id)
            out._num_bins = group['n_bins'].value
            return out

        # Search through all subclasses and find the one matching the HDF5
        # 'type'.  Call that class's from_hdf5 method.
        for subclass in cls._recursive_subclasses():
            if group['type'].value.decode() == subclass.short_name.lower():
                return subclass.from_hdf5(group, **kwargs)

        raise ValueError("Unrecognized Filter class: '"
                         + group['type'].value.decode() + "'")

    @property
    def bins(self):
        return self._bins

    @property
    def num_bins(self):
        return len(self.bins)

    @bins.setter
    def bins(self, bins):
        # Format the bins as a 1D numpy array.
        bins = np.atleast_1d(bins)

        # Check the bin values.
        self.check_bins(bins)

        self._bins = bins

    def check_bins(self, bins):
        """Make sure given bins are valid for this filter.

        Raises
        ------
        TypeError
        ValueError

        """

        pass

    def to_xml_element(self):
        """Return XML Element representing the Filter.

        Returns
        -------
        element : xml.etree.ElementTree.Element
            XML element containing filter data

        """


        element = ET.Element('filter')
        element.set('id', str(self.id))
        element.set('type', self.short_name.lower())

        subelement = ET.SubElement(element, 'bins')
        subelement.text = ' '.join(str(b) for b in self.bins)

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
        return type(self) is type(other)

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
            msg = 'Unable to merge "{0}" with "{1}" '.format(
                type(self), type(other))
            raise ValueError(msg)

        # Merge unique filter bins
        merged_bins = np.concatenate((self.bins, other.bins))
        merged_bins = np.unique(merged_bins)

        # Create a new filter with these bins and a new auto-generated ID
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

        if type(self) is not type(other):
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

        if filter_bin not in self.bins:
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

    def get_pandas_dataframe(self, data_size, stride, **kwargs):
        """Builds a Pandas DataFrame for the Filter's bins.

        This method constructs a Pandas DataFrame object for the filter with
        columns annotated by filter bin information. This is a helper method for
        :meth:`Tally.get_pandas_dataframe`.

        Parameters
        ----------
        data_size : int
            The total number of bins in the tally corresponding to this filter
        stride : int
            Stride in memory for the filter

        Keyword arguments
        -----------------
        paths : bool
            Only used for DistribcellFilter.  If True (default), expand
            distribcell indices into multi-index columns describing the path
            to that distribcell through the CSG tree.  NOTE: This option assumes
            that all distribcell paths are of the same length and do not have
            the same universes and cells but different lattice cell indices.

        Returns
        -------
        pandas.DataFrame
            A Pandas DataFrame with columns of strings that characterize the
            filter's bins. The number of rows in the DataFrame is the same as
            the total number of bins in the corresponding tally, with the filter
            bin appropriately tiled to map to the corresponding tally bins.

        See also
        --------
        Tally.get_pandas_dataframe(), CrossFilter.get_pandas_dataframe()

        """
        # Initialize Pandas DataFrame
        df = pd.DataFrame()

        filter_bins = np.repeat(self.bins, stride)
        tile_factor = data_size // len(filter_bins)
        filter_bins = np.tile(filter_bins, tile_factor)
        df = pd.concat([df, pd.DataFrame(
            {self.short_name.lower(): filter_bins})])

        return df


class WithIDFilter(Filter):
    """Abstract parent for filters of types with ids (Cell, Material, etc.)."""

    @Filter.bins.setter
    def bins(self, bins):
        # Format the bins as a 1D numpy array.
        bins = np.atleast_1d(bins)

        # Check the bin values.
        cv.check_iterable_type('filter bins', bins,
                               (Integral, self.expected_type))
        for edge in bins:
            if isinstance(edge, Integral):
                cv.check_greater_than('filter bin', edge, 0, equality=True)

        # Extract id values.
        bins = np.atleast_1d([b if isinstance(b, Integral) else b.id
                              for b in bins])

        self._bins = bins


class UniverseFilter(WithIDFilter):
    """Bins tally event locations based on the Universe they occured in.

    Parameters
    ----------
    bins : openmc.Universe, Integral, or iterable thereof
        The Universes to tally. Either openmc.Universe objects or their
        Integral ID numbers can be used.
    filter_id : int
        Unique identifier for the filter

    Attributes
    ----------
    bins : Iterable of Integral
        openmc.Universe IDs.
    id : int
        Unique identifier for the filter
    num_bins : Integral
        The number of filter bins

    """
    expected_type = Universe


class MaterialFilter(WithIDFilter):
    """Bins tally event locations based on the Material they occured in.

    Parameters
    ----------
    bins : openmc.Material, Integral, or iterable thereof
        The Materials to tally. Either openmc.Material objects or their
        Integral ID numbers can be used.
    filter_id : int
        Unique identifier for the filter

    Attributes
    ----------
    bins : Iterable of Integral
        openmc.Material IDs.
    id : int
        Unique identifier for the filter
    num_bins : Integral
        The number of filter bins

    """
    expected_type = Material


class CellFilter(WithIDFilter):
    """Bins tally event locations based on the Cell they occured in.

    Parameters
    ----------
    bins : openmc.Cell, Integral, or iterable thereof
        The Cells to tally. Either openmc.Cell objects or their
        Integral ID numbers can be used.
    filter_id : int
        Unique identifier for the filter

    Attributes
    ----------
    bins : Iterable of Integral
        openmc.Cell IDs.
    id : int
        Unique identifier for the filter
    num_bins : Integral
        The number of filter bins

    """
    expected_type = Cell


class CellFromFilter(WithIDFilter):
    """Bins tally on which Cell the neutron came from.

    Parameters
    ----------
    bins : openmc.Cell, Integral, or iterable thereof
        The Cell(s) to tally. Either openmc.Cell objects or their
        Integral ID numbers can be used.
    filter_id : int
        Unique identifier for the filter

    Attributes
    ----------
    bins : Integral or Iterable of Integral
        openmc.Cell IDs.
    id : int
        Unique identifier for the filter
    num_bins : Integral
        The number of filter bins

    """
    expected_type = Cell


class CellbornFilter(WithIDFilter):
    """Bins tally events based on which Cell the neutron was born in.

    Parameters
    ----------
    bins : openmc.Cell, Integral, or iterable thereof
        The birth Cells to tally. Either openmc.Cell objects or their
        Integral ID numbers can be used.
    filter_id : int
        Unique identifier for the filter

    Attributes
    ----------
    bins : Iterable of Integral
        openmc.Cell IDs.
    id : int
        Unique identifier for the filter
    num_bins : Integral
        The number of filter bins

    """
    expected_type = Cell


class SurfaceFilter(Filter):
    """Bins particle currents on Mesh surfaces.

    Parameters
    ----------
    bins : Iterable of Integral
        Indices corresponding to which face of a mesh cell the current is
        crossing.
    filter_id : int
        Unique identifier for the filter

    Attributes
    ----------
    bins : Iterable of Integral
        Indices corresponding to which face of a mesh cell the current is
        crossing.
    id : int
        Unique identifier for the filter
    num_bins : Integral
        The number of filter bins

    """
    @Filter.bins.setter
    def bins(self, bins):
        # Format the bins as a 1D numpy array.
        bins = np.atleast_1d(bins)

        # Check the bin values.
        cv.check_iterable_type('filter bins', bins, Integral)
        for edge in bins:
            cv.check_greater_than('filter bin', edge, 0, equality=True)

        self._bins = bins

    @property
    def num_bins(self):
        # Need to handle number of bins carefully -- for surface current
        # tallies, the number of bins depends on the mesh, which we don't have a
        # reference to in this filter
        return self._num_bins

    def get_pandas_dataframe(self, data_size, stride, **kwargs):
        """Builds a Pandas DataFrame for the Filter's bins.

        This method constructs a Pandas DataFrame object for the filter with
        columns annotated by filter bin information. This is a helper method for
        :meth:`Tally.get_pandas_dataframe`.

        Parameters
        ----------
        data_size : int
            The total number of bins in the tally corresponding to this filter
        stride : int
            Stride in memory for the filter

        Returns
        -------
        pandas.DataFrame
            A Pandas DataFrame with a column of strings describing which surface
            the current is crossing and which direction it points.  The number
            of rows in the DataFrame is the same as the total number of bins in
            the corresponding tally, with the filter bin appropriately tiled to
            map to the corresponding tally bins.

        See also
        --------
        Tally.get_pandas_dataframe(), CrossFilter.get_pandas_dataframe()

        """
        # Initialize Pandas DataFrame
        df = pd.DataFrame()

        filter_bins = np.repeat(self.bins, stride)
        tile_factor = data_size // len(filter_bins)
        filter_bins = np.tile(filter_bins, tile_factor)
        filter_bins = [_CURRENT_NAMES[x] for x in filter_bins]
        df = pd.concat([df, pd.DataFrame(
            {self.short_name.lower(): filter_bins})])

        return df


class MeshFilter(Filter):
    """Bins tally event locations onto a regular, rectangular mesh.

    Parameters
    ----------
    mesh : openmc.Mesh
        The Mesh object that events will be tallied onto
    filter_id : int
        Unique identifier for the filter

    Attributes
    ----------
    bins : Integral
        The Mesh ID
    mesh : openmc.Mesh
        The Mesh object that events will be tallied onto
    id : int
        Unique identifier for the filter
    num_bins : Integral
        The number of filter bins

    """

    def __init__(self, mesh, filter_id=None):
        self.mesh = mesh
        super(MeshFilter, self).__init__(mesh.id, filter_id)

    @classmethod
    def from_hdf5(cls, group, **kwargs):
        if group['type'].value.decode() != cls.short_name.lower():
            raise ValueError("Expected HDF5 data for filter type '"
                             + cls.short_name.lower() + "' but got '"
                             + group['type'].value.decode() + " instead")

        if 'meshes' not in kwargs:
            raise ValueError(cls.__name__ + " requires a 'meshes' keyword "
                             "argument.")

        mesh_id = group['bins'].value
        mesh_obj = kwargs['meshes'][mesh_id]
        filter_id = int(group.name.split('/')[-1].lstrip('filter '))

        out = cls(mesh_obj, filter_id=filter_id)
        out._num_bins = group['n_bins'].value

        return out

    @property
    def mesh(self):
        return self._mesh

    @mesh.setter
    def mesh(self, mesh):
        cv.check_type('filter mesh', mesh, openmc.Mesh)
        self._mesh = mesh
        self.bins = mesh.id

    @property
    def num_bins(self):
        try:
            return self._num_bins
        except AttributeError:
            return reduce(operator.mul, self.mesh.dimension)

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

    def get_pandas_dataframe(self, data_size, stride, **kwargs):
        """Builds a Pandas DataFrame for the Filter's bins.

        This method constructs a Pandas DataFrame object for the filter with
        columns annotated by filter bin information. This is a helper method for
        :meth:`Tally.get_pandas_dataframe`.

        Parameters
        ----------
        data_size : int
            The total number of bins in the tally corresponding to this filter
        stride : int
            Stride in memory for the filter

        Returns
        -------
        pandas.DataFrame
            A Pandas DataFrame with three columns describing the x,y,z mesh
            cell indices corresponding to each filter bin.  The number of rows
            in the DataFrame is the same as the total number of bins in the
            corresponding tally, with the filter bin appropriately tiled to map
            to the corresponding tally bins.

        See also
        --------
        Tally.get_pandas_dataframe(), CrossFilter.get_pandas_dataframe()

        """
        # Initialize Pandas DataFrame
        df = pd.DataFrame()

        # Initialize dictionary to build Pandas Multi-index column
        filter_dict = {}

        # Append Mesh ID as outermost index of multi-index
        mesh_key = 'mesh {0}'.format(self.mesh.id)

        # Find mesh dimensions - use 3D indices for simplicity
        if len(self.mesh.dimension) == 3:
            nx, ny, nz = self.mesh.dimension
        elif len(self.mesh.dimension) == 2:
            nx, ny = self.mesh.dimension
            nz = 1
        else:
            nx = self.mesh.dimension
            ny = nz = 1

        # Generate multi-index sub-column for x-axis
        filter_bins = np.arange(1, nx + 1)
        repeat_factor = ny * nz * stride
        filter_bins = np.repeat(filter_bins, repeat_factor)
        tile_factor = data_size // len(filter_bins)
        filter_bins = np.tile(filter_bins, tile_factor)
        filter_dict[(mesh_key, 'x')] = filter_bins

        # Generate multi-index sub-column for y-axis
        filter_bins = np.arange(1, ny + 1)
        repeat_factor = nz * stride
        filter_bins = np.repeat(filter_bins, repeat_factor)
        tile_factor = data_size // len(filter_bins)
        filter_bins = np.tile(filter_bins, tile_factor)
        filter_dict[(mesh_key, 'y')] = filter_bins

        # Generate multi-index sub-column for z-axis
        filter_bins = np.arange(1, nz + 1)
        repeat_factor = stride
        filter_bins = np.repeat(filter_bins, repeat_factor)
        tile_factor = data_size // len(filter_bins)
        filter_bins = np.tile(filter_bins, tile_factor)
        filter_dict[(mesh_key, 'z')] = filter_bins

        # Initialize a Pandas DataFrame from the mesh dictionary
        df = pd.concat([df, pd.DataFrame(filter_dict)])

        return df


class RealFilter(Filter):
    """Tally modifier that describes phase-space and other characteristics

    Parameters
    ----------
    bins : Iterable of Real
        A grid of bin values.
    filter_id : int
        Unique identifier for the filter

    Attributes
    ----------
    bins : Iterable of Real
        A grid of bin values.
    id : int
        Unique identifier for the filter
    num_bins : int
        The number of filter bins

    """

    def __gt__(self, other):
        if type(self) is type(other):
            # Compare largest/smallest bin edges in filters
            # This logic is used when merging tallies with real filters
            return self.bins[0] >= other.bins[-1]
        else:
            return super(RealFilter, self).__gt__(other)

    @property
    def num_bins(self):
        return len(self.bins) - 1

    def can_merge(self, other):
        if type(self) is not type(other):
            return False

        if self.bins[0] == other.bins[-1]:
            # This low edge coincides with other's high edge
            return True
        elif self.bins[-1] == other.bins[0]:
            # This high edge coincides with other's low edge
            return True
        else:
            return False

    def merge(self, other):
        if not self.can_merge(other):
            msg = 'Unable to merge "{0}" with "{1}" ' \
                  'filters'.format(type(self), type(other))
            raise ValueError(msg)

        # Merge unique filter bins
        merged_bins = np.concatenate((self.bins, other.bins))
        merged_bins = np.unique(merged_bins)

        # Create a new filter with these bins and a new auto-generated ID
        return type(self)(sorted(merged_bins))

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

        if type(self) is not type(other):
            return False
        elif len(self.bins) != len(other.bins):
            return False
        else:
            return np.allclose(self.bins, other.bins)

    def get_bin_index(self, filter_bin):
        i = np.where(self.bins == filter_bin[1])[0]
        if len(i) == 0:
            msg = 'Unable to get the bin index for Filter since "{0}" ' \
                  'is not one of the bins'.format(filter_bin)
            raise ValueError(msg)
        else:
            return i[0] - 1

    def get_bin(self, bin_index):
        cv.check_type('bin_index', bin_index, Integral)
        cv.check_greater_than('bin_index', bin_index, 0, equality=True)
        cv.check_less_than('bin_index', bin_index, self.num_bins)

        # Construct 2-tuple of lower, upper bins for real-valued filters
        return (self.bins[bin_index], self.bins[bin_index + 1])


class EnergyFilter(RealFilter):
    """Bins tally events based on incident particle energy.

    Parameters
    ----------
    bins : Iterable of Real
        A grid of energy values in eV.
    filter_id : int
        Unique identifier for the filter

    Attributes
    ----------
    bins : Iterable of Real
        A grid of energy values in eV.
    id : int
        Unique identifier for the filter
    num_bins : int
        The number of filter bins

    """

    def get_bin_index(self, filter_bin):
        # Use lower energy bound to find index for RealFilters
        deltas = np.abs(self.bins - filter_bin[1]) / filter_bin[1]
        min_delta = np.min(deltas)
        if min_delta < 1E-3:
            return deltas.argmin() - 1
        else:
            msg = 'Unable to get the bin index for Filter since "{0}" ' \
                  'is not one of the bins'.format(filter_bin)
            raise ValueError(msg)

    def check_bins(self, bins):
        for edge in bins:
            if not isinstance(edge, Real):
                msg = 'Unable to add bin edge "{0}" to a "{1}" ' \
                      'since it is a non-integer or floating point ' \
                      'value'.format(edge, type(self))
                raise ValueError(msg)
            elif edge < 0.:
                msg = 'Unable to add bin edge "{0}" to a "{1}" ' \
                      'since it is a negative value'.format(edge, type(self))
                raise ValueError(msg)

        # Check that bin edges are monotonically increasing
        for index in range(1, len(bins)):
            if bins[index] < bins[index-1]:
                msg = 'Unable to add bin edges "{0}" to a "{1}" Filter ' \
                      'since they are not monotonically ' \
                      'increasing'.format(bins, type(self))
                raise ValueError(msg)

    def get_pandas_dataframe(self, data_size, stride, **kwargs):
        """Builds a Pandas DataFrame for the Filter's bins.

        This method constructs a Pandas DataFrame object for the filter with
        columns annotated by filter bin information. This is a helper method for
        :meth:`Tally.get_pandas_dataframe`.

        Parameters
        ----------
        data_size : int
            The total number of bins in the tally corresponding to this filter
        stride : int
            Stride in memory for the filter

        Returns
        -------
        pandas.DataFrame
            A Pandas DataFrame with one column of the lower energy bound and one
            column of upper energy bound for each filter bin.  The number of
            rows in the DataFrame is the same as the total number of bins in the
            corresponding tally, with the filter bin appropriately tiled to map
            to the corresponding tally bins.

        See also
        --------
        Tally.get_pandas_dataframe(), CrossFilter.get_pandas_dataframe()

        """
        # Initialize Pandas DataFrame
        df = pd.DataFrame()

        # Extract the lower and upper energy bounds, then repeat and tile
        # them as necessary to account for other filters.
        lo_bins = np.repeat(self.bins[:-1], stride)
        hi_bins = np.repeat(self.bins[1:], stride)
        tile_factor = data_size // len(lo_bins)
        lo_bins = np.tile(lo_bins, tile_factor)
        hi_bins = np.tile(hi_bins, tile_factor)

        # Add the new energy columns to the DataFrame.
        df.loc[:, self.short_name.lower() + ' low [eV]'] = lo_bins
        df.loc[:, self.short_name.lower() + ' high [eV]'] = hi_bins

        return df


class EnergyoutFilter(EnergyFilter):
    """Bins tally events based on outgoing particle energy.

    Parameters
    ----------
    bins : Iterable of Real
        A grid of energy values in eV.
    filter_id : int
        Unique identifier for the filter

    Attributes
    ----------
    bins : Iterable of Real
        A grid of energy values in eV.
    id : int
        Unique identifier for the filter
    num_bins : int
        The number of filter bins

    """


def _path_to_levels(path):
    """Convert distribcell path to list of levels

    Parameters
    ----------
    path : str
        Distribcell path

    Returns
    -------
    list
        List of levels in path

    """
    # Split path into universes/cells/lattices
    path_items = path.split('->')

    # Pair together universe and cell information from the same level
    idx = [i for i, item in enumerate(path_items) if item.startswith('u')]
    for i in reversed(idx):
        univ_id = int(path_items.pop(i)[1:])
        cell_id = int(path_items.pop(i)[1:])
        path_items.insert(i, ('universe', univ_id, cell_id))

    # Reformat lattice into tuple
    idx = [i for i, item in enumerate(path_items) if isinstance(item, str)]
    for i in idx:
        item = path_items.pop(i)[1:-1]
        lat_id, lat_xyz = item.split('(')
        lat_id = int(lat_id)
        lat_xyz = tuple(int(x) for x in lat_xyz.split(','))
        path_items.insert(i, ('lattice', lat_id, lat_xyz))

    return path_items


class DistribcellFilter(Filter):
    """Bins tally event locations on instances of repeated cells.

    Parameters
    ----------
    cell : openmc.Cell or Integral
        The distributed cell to tally. Either an openmc.Cell or an Integral
        cell ID number can be used.
    filter_id : int
        Unique identifier for the filter

    Attributes
    ----------
    bins : Iterable of Integral
        An iterable with one element---the ID of the distributed Cell.
    id : int
        Unique identifier for the filter
    num_bins : int
        The number of filter bins
    paths : list of str
        The paths traversed through the CSG tree to reach each distribcell
        instance (for 'distribcell' filters only)

    """

    def __init__(self, cell, filter_id=None):
        self._paths = None
        super(DistribcellFilter, self).__init__(cell, filter_id)

    @classmethod
    def from_hdf5(cls, group, **kwargs):
        if group['type'].value.decode() != cls.short_name.lower():
            raise ValueError("Expected HDF5 data for filter type '"
                             + cls.short_name.lower() + "' but got '"
                             + group['type'].value.decode() + " instead")

        filter_id = int(group.name.split('/')[-1].lstrip('filter '))

        out = cls(group['bins'].value, filter_id=filter_id)
        out._num_bins = group['n_bins'].value

        return out

    @property
    def num_bins(self):
        # Need to handle number of bins carefully -- for distribcell tallies, we
        # need to know how many instances of the cell there are
        return self._num_bins

    @property
    def paths(self):
        return self._paths

    @Filter.bins.setter
    def bins(self, bins):
        # Format the bins as a 1D numpy array.
        bins = np.atleast_1d(bins)

        # Make sure there is only 1 bin.
        if not len(bins) == 1:
            msg = 'Unable to add bins "{0}" to a DistribcellFilter since ' \
                  'only a single distribcell can be used per tally'.format(bins)
            raise ValueError(msg)

        # Check the type and extract the id, if necessary.
        cv.check_type('distribcell bin', bins[0], (Integral, openmc.Cell))
        if isinstance(bins[0], openmc.Cell):
            bins = np.atleast_1d(bins[0].id)

        self._bins = bins

    @paths.setter
    def paths(self, paths):
        cv.check_iterable_type('paths', paths, str)
        self._paths = paths

    def can_merge(self, other):
        # Distribcell filters cannot have more than one bin
        return False

    def get_bin_index(self, filter_bin):
        # Filter bins for distribcells are indices of each unique placement of
        # the Cell in the Geometry (consecutive integers starting at 0).
        return filter_bin

    def get_pandas_dataframe(self, data_size, stride, **kwargs):
        """Builds a Pandas DataFrame for the Filter's bins.

        This method constructs a Pandas DataFrame object for the filter with
        columns annotated by filter bin information. This is a helper method for
        :meth:`Tally.get_pandas_dataframe`.

        Parameters
        ----------
        data_size : int
            The total number of bins in the tally corresponding to this filter
        stride : int
            Stride in memory for the filter

        Keyword arguments
        -----------------
        paths : bool
            If True (default), expand distribcell indices into multi-index
            columns describing the path to that distribcell through the CSG
            tree.  NOTE: This option assumes that all distribcell paths are of
            the same length and do not have the same universes and cells but
            different lattice cell indices.

        Returns
        -------
        pandas.DataFrame
            A Pandas DataFrame with columns describing distributed cells.  The
            for will be either:

            1. a single column with the cell instance IDs (without summary info)
            2. separate columns for the cell IDs, universe IDs, and lattice IDs
               and x,y,z cell indices corresponding to each (distribcell paths).

            The number of rows in the DataFrame is the same as the total number
            of bins in the corresponding tally, with the filter bin
            appropriately tiled to map to the corresponding tally bins.

        See also
        --------
        Tally.get_pandas_dataframe(), CrossFilter.get_pandas_dataframe()

        """
        # Initialize Pandas DataFrame
        df = pd.DataFrame()

        level_df = None

        paths = kwargs.setdefault('paths', True)

        # Create Pandas Multi-index columns for each level in CSG tree
        if paths:

            # Distribcell paths require linked metadata from the Summary
            if self.paths is None:
                msg = 'Unable to construct distribcell paths since ' \
                      'the Summary is not linked to the StatePoint'
                raise ValueError(msg)

            # Make copy of array of distribcell paths to use in
            # Pandas Multi-index column construction
            num_offsets = len(self.paths)
            paths = [_path_to_levels(p) for p in self.paths]

            # Loop over CSG levels in the distribcell paths
            num_levels = len(paths[0])
            for i_level in range(num_levels):
                # Use level key as first index in Pandas Multi-index column
                level_key = 'level {}'.format(i_level + 1)

                # Create a dictionary for this level for Pandas Multi-index
                level_dict = OrderedDict()

                # Use the first distribcell path to determine if level
                # is a universe/cell or lattice level
                path = paths[0]
                if path[i_level][0] == 'lattice':
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
                    if len(path[i_level][2]) == 3:
                        level_dict[lat_z_key] = np.empty(num_offsets)

                else:
                    # Initialize prefix Multi-index keys
                    univ_key = (level_key, 'univ', 'id')
                    cell_key = (level_key, 'cell', 'id')

                    # Allocate NumPy arrays for each CSG level and
                    # each Multi-index column in the DataFrame
                    level_dict[univ_key] = np.empty(num_offsets)
                    level_dict[cell_key] = np.empty(num_offsets)

                # Populate Multi-index arrays with all distribcell paths
                for i, path in enumerate(paths):

                    level = path[i_level]
                    if level[0] == 'lattice':
                        # Assign entry to Lattice Multi-index column
                        level_dict[lat_id_key][i] = level[1]
                        level_dict[lat_x_key][i] = level[2][0]
                        level_dict[lat_y_key][i] = level[2][1]
                        if len(level[2]) == 3:
                            level_dict[lat_z_key][i] = level[2][2]

                    else:
                        # Assign entry to Universe, Cell Multi-index columns
                        level_dict[univ_key][i] = level[1]
                        level_dict[cell_key][i] = level[2]

                # Tile the Multi-index columns
                for level_key, level_bins in level_dict.items():
                    level_bins = np.repeat(level_bins, stride)
                    tile_factor = data_size // len(level_bins)
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
        filter_bins = np.repeat(filter_bins, stride)
        tile_factor = data_size // len(filter_bins)
        filter_bins = np.tile(filter_bins, tile_factor)
        df = pd.DataFrame({self.short_name.lower() : filter_bins})

        # Concatenate with DataFrame of distribcell instance IDs
        if level_df is not None:
            level_df = level_df.dropna(axis=1, how='all')
            level_df = level_df.astype(np.int)
            df = pd.concat([level_df, df], axis=1)

        return df


class MuFilter(RealFilter):
    """Bins tally events based on particle scattering angle.

    Parameters
    ----------
    bins : Iterable of Real or Integral
        A grid of scattering angles which events will binned into.  Values
        represent the cosine of the scattering angle.  If an Iterable is given,
        the values will be used explicitly as grid points.  If a single Integral
        is given, the range [-1, 1] will be divided up equally into that number
        of bins.
    filter_id : int
        Unique identifier for the filter

    Attributes
    ----------
    bins : Integral
        A grid of scattering angles which events will binned into.  Values
        represent the cosine of the scattering angle.  If an Iterable is given,
        the values will be used explicitly as grid points.  If a single Integral
        is given, the range [-1, 1] will be divided up equally into that number
        of bins.
    id : int
        Unique identifier for the filter
    num_bins : Integral
        The number of filter bins

    """

    def check_bins(self, bins):
        for edge in bins:
            if not isinstance(edge, Real):
                msg = 'Unable to add bin edge "{0}" to a "{1}" ' \
                      'since it is a non-integer or floating point ' \
                      'value'.format(edge, type(self))
                raise ValueError(msg)
            elif edge < -1.:
                msg = 'Unable to add bin edge "{0}" to a "{1}" ' \
                      'since it is less than -1'.format(edge, type(self))
                raise ValueError(msg)
            elif edge > 1.:
                msg = 'Unable to add bin edge "{0}" to a "{1}" ' \
                      'since it is greater than 1'.format(edge, type(self))
                raise ValueError(msg)

        # Check that bin edges are monotonically increasing
        for index in range(1, len(bins)):
            if bins[index] < bins[index-1]:
                msg = 'Unable to add bin edges "{0}" to a "{1}" Filter ' \
                      'since they are not monotonically ' \
                      'increasing'.format(bins, type(self))
                raise ValueError(msg)

    def get_pandas_dataframe(self, data_size, stride, **kwargs):
        """Builds a Pandas DataFrame for the Filter's bins.

        This method constructs a Pandas DataFrame object for the filter with
        columns annotated by filter bin information. This is a helper method
        for :meth:`Tally.get_pandas_dataframe`.

        Parameters
        ----------
        data_size : int
            The total number of bins in the tally corresponding to this filter
        stride : int
            Stride in memory for the filter

        Returns
        -------
        pandas.DataFrame
            A Pandas DataFrame with one column of the lower energy bound and one
            column of upper energy bound for each filter bin.  The number of
            rows in the DataFrame is the same as the total number of bins in the
            corresponding tally, with the filter bin appropriately tiled to map
            to the corresponding tally bins.

        See also
        --------
        Tally.get_pandas_dataframe(), CrossFilter.get_pandas_dataframe()

        """
        # Initialize Pandas DataFrame
        df = pd.DataFrame()

        # Extract the lower and upper energy bounds, then repeat and tile
        # them as necessary to account for other filters.
        lo_bins = np.repeat(self.bins[:-1], stride)
        hi_bins = np.repeat(self.bins[1:], stride)
        tile_factor = data_size // len(lo_bins)
        lo_bins = np.tile(lo_bins, tile_factor)
        hi_bins = np.tile(hi_bins, tile_factor)

        # Add the new energy columns to the DataFrame.
        df.loc[:, self.short_name.lower() + ' low'] = lo_bins
        df.loc[:, self.short_name.lower() + ' high'] = hi_bins

        return df


class PolarFilter(RealFilter):
    """Bins tally events based on the incident particle's direction.

    Parameters
    ----------
    bins : Iterable of Real or Integral
        A grid of polar angles which events will binned into.  Values represent
        an angle in radians relative to the z-axis.  If an Iterable is given,
        the values will be used explicitly as grid points.  If a single Integral
        is given, the range [0, pi] will be divided up equally into that number
        of bins.
    filter_id : int
        Unique identifier for the filter

    Attributes
    ----------
    bins : Iterable of Real or Integral
        A grid of polar angles which events will binned into.  Values represent
        an angle in radians relative to the z-axis.  If an Iterable is given,
        the values will be used explicitly as grid points.  If a single Integral
        is given, the range [0, pi] will be divided up equally into that number
        of bins.
    id : int
        Unique identifier for the filter
    num_bins : Integral
        The number of filter bins

    """

    def check_bins(self, bins):
        for edge in bins:
            if not isinstance(edge, Real):
                msg = 'Unable to add bin edge "{0}" to a "{1}" ' \
                      'since it is a non-integer or floating point ' \
                      'value'.format(edge, type(self))
                raise ValueError(msg)
            elif edge < 0.:
                msg = 'Unable to add bin edge "{0}" to a "{1}" ' \
                      'since it is less than 0'.format(edge, type(self))
                raise ValueError(msg)
            elif not np.isclose(edge, np.pi) and edge > np.pi:
                msg = 'Unable to add bin edge "{0}" to a "{1}" ' \
                      'since it is greater than pi'.format(edge, type(self))
                raise ValueError(msg)

        # Check that bin edges are monotonically increasing
        for index in range(1, len(bins)):
            if bins[index] < bins[index-1]:
                msg = 'Unable to add bin edges "{0}" to a "{1}" Filter ' \
                      'since they are not monotonically ' \
                      'increasing'.format(bins, type(self))
                raise ValueError(msg)

    def get_pandas_dataframe(self, data_size, stride, **kwargs):
        """Builds a Pandas DataFrame for the Filter's bins.

        This method constructs a Pandas DataFrame object for the filter with
        columns annotated by filter bin information. This is a helper method for
        :meth:`Tally.get_pandas_dataframe`.

        Parameters
        ----------
        data_size : int
            The total number of bins in the tally corresponding to this filter
        stride : int
            Stride in memory for the filter

        Returns
        -------
        pandas.DataFrame
            A Pandas DataFrame with a column corresponding to the lower polar
            angle bound for each of the filter's bins.  The number of rows in
            the DataFrame is the same as the total number of bins in the
            corresponding tally, with the filter bin appropriately tiled to map
            to the corresponding tally bins.

        See also
        --------
        Tally.get_pandas_dataframe(), CrossFilter.get_pandas_dataframe()

        """
        # Initialize Pandas DataFrame
        df = pd.DataFrame()

        # Extract the lower and upper angle bounds, then repeat and tile
        # them as necessary to account for other filters.
        lo_bins = np.repeat(self.bins[:-1], stride)
        hi_bins = np.repeat(self.bins[1:], stride)
        tile_factor = data_size // len(lo_bins)
        lo_bins = np.tile(lo_bins, tile_factor)
        hi_bins = np.tile(hi_bins, tile_factor)

        # Add the new angle columns to the DataFrame.
        df.loc[:, 'polar low'] = lo_bins
        df.loc[:, 'polar high'] = hi_bins

        return df


class AzimuthalFilter(RealFilter):
    """Bins tally events based on the incident particle's direction.

    Parameters
    ----------
    bins : Iterable of Real or Integral
        A grid of azimuthal angles which events will binned into.  Values
        represent an angle in radians relative to the x-axis and perpendicular
        to the z-axis.  If an Iterable is given, the values will be used
        explicitly as grid points.  If a single Integral is given, the range
        [-pi, pi) will be divided up equally into that number of bins.
    filter_id : int
        Unique identifier for the filter

    Attributes
    ----------
    bins : Iterable of Real or Integral
        A grid of azimuthal angles which events will binned into.  Values
        represent an angle in radians relative to the x-axis and perpendicular
        to the z-axis.  If an Iterable is given, the values will be used
        explicitly as grid points.  If a single Integral is given, the range
        [-pi, pi) will be divided up equally into that number of bins.
    id : int
        Unique identifier for the filter
    num_bins : Integral
        The number of filter bins

    """

    def check_bins(self, bins):
        for edge in bins:
            if not isinstance(edge, Real):
                msg = 'Unable to add bin edge "{0}" to a "{1}" ' \
                      'since it is a non-integer or floating point ' \
                      'value'.format(edge, type(self))
                raise ValueError(msg)
            elif not np.isclose(edge, -np.pi) and edge < -np.pi:
                msg = 'Unable to add bin edge "{0}" to a "{1}" ' \
                      'since it is less than -pi'.format(edge, type(self))
                raise ValueError(msg)
            elif not np.isclose(edge, np.pi) and edge > np.pi:
                msg = 'Unable to add bin edge "{0}" to a "{1}" ' \
                      'since it is greater than pi'.format(edge, type(self))
                raise ValueError(msg)

        # Check that bin edges are monotonically increasing
        for index in range(1, len(bins)):
            if bins[index] < bins[index-1]:
                msg = 'Unable to add bin edges "{0}" to a "{1}" Filter ' \
                      'since they are not monotonically ' \
                      'increasing'.format(bins, type(self))
                raise ValueError(msg)

    def get_pandas_dataframe(self, data_size, stride, paths=True):
        """Builds a Pandas DataFrame for the Filter's bins.

        This method constructs a Pandas DataFrame object for the filter with
        columns annotated by filter bin information. This is a helper method for
        :meth:`Tally.get_pandas_dataframe`.

        Parameters
        ----------
        data_size : int
            The total number of bins in the tally corresponding to this filter
        stride : int
            Stride in memory for the filter

        Returns
        -------
        pandas.DataFrame
            A Pandas DataFrame with a column corresponding to the lower
            azimuthal angle bound for each of the filter's bins.  The number of
            rows in the DataFrame is the same as the total number of bins in the
            corresponding tally, with the filter bin appropriately tiled to map
            to the corresponding tally bins.

        See also
        --------
        Tally.get_pandas_dataframe(), CrossFilter.get_pandas_dataframe()

        """
        # Initialize Pandas DataFrame
        df = pd.DataFrame()

        # Extract the lower and upper angle bounds, then repeat and tile
        # them as necessary to account for other filters.
        lo_bins = np.repeat(self.bins[:-1], stride)
        hi_bins = np.repeat(self.bins[1:], stride)
        tile_factor = data_size // len(lo_bins)
        lo_bins = np.tile(lo_bins, tile_factor)
        hi_bins = np.tile(hi_bins, tile_factor)

        # Add the new angle columns to the DataFrame.
        df.loc[:, 'azimuthal low'] = lo_bins
        df.loc[:, 'azimuthal high'] = hi_bins

        return df


class DelayedGroupFilter(Filter):
    """Bins fission events based on the produced neutron precursor groups.

    Parameters
    ----------
    bins : Integral or Iterable of Integral
        The delayed neutron precursor groups.  For example, ENDF/B-VII.1 uses
        6 precursor groups so a tally with all groups will have bins =
        [1, 2, 3, 4, 5, 6].
    filter_id : int
        Unique identifier for the filter

    Attributes
    ----------
    bins : Integral or Iterable of Integral
        The delayed neutron precursor groups.  For example, ENDF/B-VII.1 uses
        6 precursor groups so a tally with all groups will have bins =
        [1, 2, 3, 4, 5, 6].
    id : int
        Unique identifier for the filter
    num_bins : Integral
        The number of filter bins

    """
    @Filter.bins.setter
    def bins(self, bins):
        # Format the bins as a 1D numpy array.
        bins = np.atleast_1d(bins)

        # Check the bin values.
        cv.check_iterable_type('filter bins', bins, Integral)
        for edge in bins:
            cv.check_greater_than('filter bin', edge, 0, equality=True)

        self._bins = bins


class EnergyFunctionFilter(Filter):
    """Multiplies tally scores by an arbitrary function of incident energy.

    The arbitrary function is described by a piecewise linear-linear
    interpolation of energy and y values.  Values outside of the given energy
    range will be evaluated as zero.

    Parameters
    ----------
    energy : Iterable of Real
        A grid of energy values in eV.
    y : iterable of Real
        A grid of interpolant values in eV.
    filter_id : int
        Unique identifier for the filter

    Attributes
    ----------
    energy : Iterable of Real
        A grid of energy values in eV.
    y : iterable of Real
        A grid of interpolant values in eV.
    id : int
        Unique identifier for the filter
    num_bins : Integral
        The number of filter bins (always 1 for this filter)

    """

    def __init__(self, energy, y, filter_id=None):
        self.energy = energy
        self.y = y
        self.id = filter_id

    def __eq__(self, other):
        if type(self) is not type(other):
            return False
        elif not all(self.energy == other.energy):
            return False
        elif not all(self.y == other.y):
            return False
        else:
            return True

    def __gt__(self, other):
        if type(self) is not type(other):
            if self.short_name in _FILTER_TYPES and \
                other.short_name in _FILTER_TYPES:
                delta = _FILTER_TYPES.index(self.short_name) - \
                        _FILTER_TYPES.index(other.short_name)
                return delta > 0
            else:
                return False
        else:
            return False

    def __lt__(self, other):
        if type(self) is not type(other):
            if self.short_name in _FILTER_TYPES and \
                other.short_name in _FILTER_TYPES:
                delta = _FILTER_TYPES.index(self.short_name) - \
                        _FILTER_TYPES.index(other.short_name)
                return delta < 0
            else:
                return False
        else:
            return False

    def __hash__(self):
        string = type(self).__name__ + '\n'
        string += '{: <16}=\t{}\n'.format('\tEnergy', self.energy)
        string += '{: <16}=\t{}\n'.format('\tInterpolant', self.y)
        return hash(string)

    def __repr__(self):
        string = type(self).__name__ + '\n'
        string += '{: <16}=\t{}\n'.format('\tEnergy', self.energy)
        string += '{: <16}=\t{}\n'.format('\tInterpolant', self.y)
        string += '{: <16}=\t{}\n'.format('\tID', self.id)
        return string

    @classmethod
    def from_hdf5(cls, group, **kwargs):
        if group['type'].value.decode() != cls.short_name.lower():
            raise ValueError("Expected HDF5 data for filter type '"
                             + cls.short_name.lower() + "' but got '"
                             + group['type'].value.decode() + " instead")

        energy = group['energy'].value
        y = group['y'].value
        filter_id = int(group.name.split('/')[-1].lstrip('filter '))

        return cls(energy, y, filter_id=filter_id)

    @classmethod
    def from_tabulated1d(cls, tab1d):
        """Construct a filter from a Tabulated1D object.

        Parameters
        ----------
        tab1d : openmc.data.Tabulated1D
            A linear-linear Tabulated1D object with only a single interpolation
            region.

        Returns
        -------
        EnergyFunctionFilter

        """
        cv.check_type('EnergyFunctionFilter tab1d', tab1d,
                      openmc.data.Tabulated1D)
        if tab1d.n_regions > 1:
            raise ValueError('Only Tabulated1Ds with a single interpolation '
                             'region are supported')
        if tab1d.interpolation[0] != 2:
            raise ValueError('Only linear-linar Tabulated1Ds are supported')

        return cls(tab1d.x, tab1d.y)

    @property
    def energy(self):
        return self._energy

    @property
    def y(self):
        return self._y

    @property
    def bins(self):
        raise AttributeError('EnergyFunctionFilters have no bins.')

    @property
    def num_bins(self):
        return 1

    @energy.setter
    def energy(self, energy):
        # Format the bins as a 1D numpy array.
        energy = np.atleast_1d(energy)

        # Make sure the values are Real and positive.
        cv.check_type('filter energy grid', energy, Iterable, Real)
        for E in energy:
            cv.check_greater_than('filter energy grid', E, 0, equality=True)

        self._energy = energy

    @y.setter
    def y(self, y):
        # Format the bins as a 1D numpy array.
        y = np.atleast_1d(y)

        # Make sure the values are Real.
        cv.check_type('filter interpolant values', y, Iterable, Real)

        self._y = y

    @bins.setter
    def bins(self, bins):
        raise RuntimeError('EnergyFunctionFilters have no bins.')

    def to_xml_element(self):
        element = ET.Element('filter')
        element.set('id', str(self.id))
        element.set('type', self.short_name.lower())

        subelement = ET.SubElement(element, 'energy')
        subelement.text = ' '.join(str(e) for e in self.energy)

        subelement = ET.SubElement(element, 'y')
        subelement.text = ' '.join(str(y) for y in self.y)

        return element

    def can_merge(self, other):
        return False

    def is_subset(self, other):
        return self == other

    def get_bin_index(self, filter_bin):
        # This filter only has one bin.  Always return 0.
        return 0

    def get_bin(self, bin_index):
        """This function is invalid for EnergyFunctionFilters."""
        raise RuntimeError('EnergyFunctionFilters have no get_bin() method')

    def get_pandas_dataframe(self, data_size, stride, **kwargs):
        """Builds a Pandas DataFrame for the Filter's bins.

        This method constructs a Pandas DataFrame object for the filter with
        columns annotated by filter bin information. This is a helper method for
        :meth:`Tally.get_pandas_dataframe`.

        Parameters
        ----------
        data_size : int
            The total number of bins in the tally corresponding to this filter
        stride : int
            Stride in memory for the filter

        Returns
        -------
        pandas.DataFrame
            A Pandas DataFrame with a column that is filled with a hash of this
            filter. EnergyFunctionFilters have only 1 bin so the purpose of this
            DataFrame column is to differentiate the filter from other
            EnergyFunctionFilters. The number of rows in the DataFrame is the
            same as the total number of bins in the corresponding tally.

        See also
        --------
        Tally.get_pandas_dataframe(), CrossFilter.get_pandas_dataframe()

        """
        df = pd.DataFrame()

        # There is no clean way of sticking all the energy, y data into a
        # DataFrame so instead we'll just make a column with the filter name
        # and fill it with a hash of the __repr__.  We want a hash that is
        # reproducible after restarting the interpreter so we'll use hashlib.md5
        # rather than the intrinsic hash().
        hash_fun = hashlib.md5()
        hash_fun.update(repr(self).encode('utf-8'))
        out = hash_fun.hexdigest()

        # The full 16 bytes make for a really wide column.  Just 7 bytes (14
        # hex characters) of the digest are probably sufficient.
        out = out[:14]

        filter_bins = np.repeat(out, stride)
        tile_factor = data_size // len(filter_bins)
        filter_bins = np.tile(filter_bins, tile_factor)
        df = pd.concat([df, pd.DataFrame(
            {self.short_name.lower(): filter_bins})])

        return df
