from abc import ABCMeta
from collections import OrderedDict
from collections.abc import Iterable
import copy
import hashlib
from itertools import product
from numbers import Real, Integral
from xml.etree import ElementTree as ET

import numpy as np
import pandas as pd

import openmc
import openmc.checkvalue as cv
from .cell import Cell
from .material import Material
from .mixin import IDManagerMixin
from .surface import Surface
from .universe import Universe


_FILTER_TYPES = (
    'universe', 'material', 'cell', 'cellborn', 'surface', 'mesh', 'energy',
    'energyout', 'mu', 'polar', 'azimuthal', 'distribcell', 'delayedgroup',
    'energyfunction', 'cellfrom', 'legendre', 'spatiallegendre',
    'sphericalharmonics', 'zernike', 'zernikeradial', 'particle', 'cellinstance'
)

_CURRENT_NAMES = (
    'x-min out', 'x-min in', 'x-max out', 'x-max in',
    'y-min out', 'y-min in', 'y-max out', 'y-max in',
    'z-min out', 'z-min in', 'z-max out', 'z-max in'
)

_PARTICLES = {'neutron', 'photon', 'electron', 'positron'}


class FilterMeta(ABCMeta):
    """Metaclass for filters that ensures class names are appropriate."""

    def __new__(cls, name, bases, namespace, **kwargs):
        # Check the class name.
        required_suffix = 'Filter'
        if not name.endswith(required_suffix):
            raise ValueError("All filter class names must end with 'Filter'")

        # Create a 'short_name' attribute that removes the 'Filter' suffix.
        namespace['short_name'] = name[:-len(required_suffix)]

        # Subclass methods can sort of inherit the docstring of parent class
        # methods.  If a function is defined without a docstring, most (all?)
        # Python interpreters will search through the parent classes to see if
        # there is a docstring for a function with the same name, and they will
        # use that docstring.  However, Sphinx does not have that functionality.
        # This chunk of code handles this docstring inheritance manually so that
        # the autodocumentation will pick it up.
        if name != required_suffix:
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
        return super().__new__(cls, name, bases, namespace, **kwargs)


def _repeat_and_tile(bins, repeat_factor, data_size):
    filter_bins = np.repeat(bins, repeat_factor)
    tile_factor = data_size // len(filter_bins)
    return np.tile(filter_bins, tile_factor)


class Filter(IDManagerMixin, metaclass=FilterMeta):
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
        else:
            return np.allclose(self.bins, other.bins)

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
            Dictionary mapping integer IDs to openmc.MeshBase objects.  Only
            used for openmc.MeshFilter objects.

        """

        filter_id = int(group.name.split('/')[-1].lstrip('filter '))

        # If the HDF5 'type' variable matches this class's short_name, then
        # there is no overriden from_hdf5 method.  Pass the bins to __init__.
        if group['type'][()].decode() == cls.short_name.lower():
            out = cls(group['bins'][()], filter_id=filter_id)
            out._num_bins = group['n_bins'][()]
            return out

        # Search through all subclasses and find the one matching the HDF5
        # 'type'.  Call that class's from_hdf5 method.
        for subclass in cls._recursive_subclasses():
            if group['type'][()].decode() == subclass.short_name.lower():
                return subclass.from_hdf5(group, **kwargs)

        raise ValueError("Unrecognized Filter class: '"
                         + group['type'][()].decode() + "'")

    @property
    def bins(self):
        return self._bins

    @bins.setter
    def bins(self, bins):
        self.check_bins(bins)
        self._bins = bins

    @property
    def num_bins(self):
        return len(self.bins)

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
        merged_bins = np.unique(merged_bins, axis=0)

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

        for b in other.bins:
            if b not in self.bins:
                return False

        return True

    def get_bin_index(self, filter_bin):
        """Returns the index in the Filter for some bin.

        Parameters
        ----------
        filter_bin : int or tuple
            The bin is the integer ID for 'material', 'surface', 'cell',
            'cellborn', and 'universe' Filters. The bin is an integer for the
            cell instance ID for 'distribcell' Filters. The bin is a 2-tuple of
            floats for 'energy' and 'energyout' filters corresponding to the
            energy boundaries of the bin of interest. The bin is an (x,y,z)
            3-tuple for 'mesh' filters corresponding to the mesh cell of
            interest.

        Returns
        -------
        filter_index : int
             The index in the Tally data array for this filter bin.

        """

        if filter_bin not in self.bins:
            msg = 'Unable to get the bin index for Filter since "{0}" ' \
                  'is not one of the bins'.format(filter_bin)
            raise ValueError(msg)

        if isinstance(self.bins, np.ndarray):
            return np.where(self.bins == filter_bin)[0][0]
        else:
            return self.bins.index(filter_bin)

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
    """Abstract parent for filters of types with IDs (Cell, Material, etc.)."""
    def __init__(self, bins, filter_id=None):
        bins = np.atleast_1d(bins)

        # Make sure bins are either integers or appropriate objects
        cv.check_iterable_type('filter bins', bins,
                               (Integral, self.expected_type))

        # Extract ID values
        bins = np.array([b if isinstance(b, Integral) else b.id
                         for b in bins])
        super().__init__(bins, filter_id)

    def check_bins(self, bins):
        # Check the bin values.
        for edge in bins:
            cv.check_greater_than('filter bin', edge, 0, equality=True)


class UniverseFilter(WithIDFilter):
    """Bins tally event locations based on the Universe they occured in.

    Parameters
    ----------
    bins : openmc.Universe, int, or iterable thereof
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
    bins : openmc.Cell, int, or iterable thereof
        The cells to tally. Either openmc.Cell objects or their ID numbers can
        be used.
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


class CellInstanceFilter(Filter):
    """Bins tally events based on which cell instance a particle is in.

    This filter is similar to :class:`DistribcellFilter` but allows one to
    select particular instances to be tallied (instead of obtaining *all*
    instances by default) and allows instances from different cells to be
    specified in a single filter.

    Parameters
    ----------
    bins : iterable of 2-tuples or numpy.ndarray
        The cell instances to tally, given as 2-tuples. For the first value in
        the tuple, either openmc.Cell objects or their integral ID numbers can
        be used.
    filter_id : int
        Unique identifier for the filter

    Attributes
    ----------
    bins : numpy.ndarray
        2D numpy array of cell IDs and instances
    id : int
        Unique identifier for the filter
    num_bins : Integral
        The number of filter bins

    See Also
    --------
    DistribcellFilter

    """
    def __init__(self, bins, filter_id=None):
        self.bins = bins
        self.id = filter_id

    @Filter.bins.setter
    def bins(self, bins):
        pairs = np.empty((len(bins), 2), dtype=int)
        for i, (cell, instance) in enumerate(bins):
            cv.check_type('cell', cell, (openmc.Cell, Integral))
            cv.check_type('instance', instance, Integral)
            pairs[i, 0] = cell if isinstance(cell, Integral) else cell.id
            pairs[i, 1] = instance
        self._bins = pairs

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
            A Pandas DataFrame with a multi-index column for the cell instance.
            The number of rows in the DataFrame is the same as the total number
            of bins in the corresponding tally, with the filter bin appropriately
            tiled to map to the corresponding tally bins.

        See also
        --------
        Tally.get_pandas_dataframe(), CrossFilter.get_pandas_dataframe()

        """
        # Repeat and tile bins as necessary to account for other filters.
        bins = np.repeat(self.bins, stride, axis=0)
        tile_factor = data_size // len(bins)
        bins = np.tile(bins, (tile_factor, 1))

        columns = pd.MultiIndex.from_product([[self.short_name.lower()],
                                              ['cell', 'instance']])
        return pd.DataFrame(bins, columns=columns)

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
        subelement.text = ' '.join(str(i) for i in self.bins.ravel())
        return element


class SurfaceFilter(WithIDFilter):
    """Filters particles by surface crossing

    Parameters
    ----------
    bins : openmc.Surface, int, or iterable of Integral
        The surfaces to tally over. Either openmc.Surface objects or their ID
        numbers can be used.
    filter_id : int
        Unique identifier for the filter

    Attributes
    ----------
    bins : Iterable of Integral
        The surfaces to tally over. Either openmc.Surface objects or their ID
        numbers can be used.
    id : int
        Unique identifier for the filter
    num_bins : Integral
        The number of filter bins

    """
    expected_type = Surface


class ParticleFilter(Filter):
    """Bins tally events based on the Particle type.

    Parameters
    ----------
    bins : str, or iterable of str
        The particles to tally represented as strings ('neutron', 'photon',
        'electron', 'positron').
    filter_id : int
        Unique identifier for the filter

    Attributes
    ----------
    bins : Iterable of Integral
        The Particles to tally
    id : int
        Unique identifier for the filter
    num_bins : Integral
        The number of filter bins

    """
    def __eq__(self, other):
        if type(self) is not type(other):
            return False
        elif len(self.bins) != len(other.bins):
            return False
        else:
            return np.all(self.bins == other.bins)

    __hash__ = Filter.__hash__

    @Filter.bins.setter
    def bins(self, bins):
        bins = np.atleast_1d(bins)
        cv.check_iterable_type('filter bins', bins, str)
        for edge in bins:
            cv.check_value('filter bin', edge, _PARTICLES)
        self._bins = bins

    @classmethod
    def from_hdf5(cls, group, **kwargs):
        if group['type'][()].decode() != cls.short_name.lower():
            raise ValueError("Expected HDF5 data for filter type '"
                             + cls.short_name.lower() + "' but got '"
                             + group['type'][()].decode() + " instead")

        particles = [b.decode() for b in group['bins'][()]]
        filter_id = int(group.name.split('/')[-1].lstrip('filter '))
        return cls(particles, filter_id=filter_id)


class MeshFilter(Filter):
    """Bins tally event locations onto a regular, rectangular mesh.

    Parameters
    ----------
    mesh : openmc.MeshBase
        The mesh object that events will be tallied onto
    filter_id : int
        Unique identifier for the filter

    Attributes
    ----------
    mesh : openmc.MeshBase
        The mesh object that events will be tallied onto
    id : int
        Unique identifier for the filter
    bins : list of tuple
        A list of mesh indices for each filter bin, e.g. [(1, 1, 1), (2, 1, 1),
        ...]
    num_bins : Integral
        The number of filter bins

    """

    def __init__(self, mesh, filter_id=None):
        self.mesh = mesh
        self.id = filter_id

    def __hash__(self):
        string = type(self).__name__ + '\n'
        string += '{: <16}=\t{}\n'.format('\tMesh ID', self.mesh.id)
        return hash(string)

    def __repr__(self):
        string = type(self).__name__ + '\n'
        string += '{: <16}=\t{}\n'.format('\tMesh ID', self.mesh.id)
        string += '{: <16}=\t{}\n'.format('\tID', self.id)
        return string

    @classmethod
    def from_hdf5(cls, group, **kwargs):
        if group['type'][()].decode() != cls.short_name.lower():
            raise ValueError("Expected HDF5 data for filter type '"
                             + cls.short_name.lower() + "' but got '"
                             + group['type'][()].decode() + " instead")

        if 'meshes' not in kwargs:
            raise ValueError(cls.__name__ + " requires a 'meshes' keyword "
                             "argument.")

        mesh_id = group['bins'][()]
        mesh_obj = kwargs['meshes'][mesh_id]
        filter_id = int(group.name.split('/')[-1].lstrip('filter '))

        out = cls(mesh_obj, filter_id=filter_id)

        return out

    @property
    def mesh(self):
        return self._mesh

    @mesh.setter
    def mesh(self, mesh):
        cv.check_type('filter mesh', mesh, openmc.MeshBase)
        self._mesh = mesh
        if isinstance(mesh, openmc.UnstructuredMesh):
            self.bins = list(range(len(mesh.volumes)))
        else:
            self.bins = list(mesh.indices)

    def can_merge(self, other):
        # Mesh filters cannot have more than one bin
        return False

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

        # Append mesh ID as outermost index of multi-index
        mesh_key = 'mesh {}'.format(self.mesh.id)

        # Find mesh dimensions - use 3D indices for simplicity
        n_dim = len(self.mesh.dimension)
        if n_dim == 3:
            nx, ny, nz = self.mesh.dimension
        elif n_dim == 2:
            nx, ny = self.mesh.dimension
            nz = 1
        else:
            nx = self.mesh.dimension
            ny = nz = 1

        # Generate multi-index sub-column for x-axis
        filter_dict[mesh_key, 'x'] = _repeat_and_tile(
            np.arange(1, nx + 1), stride, data_size)

        # Generate multi-index sub-column for y-axis
        filter_dict[mesh_key, 'y'] = _repeat_and_tile(
            np.arange(1, ny + 1), nx * stride, data_size)

        # Generate multi-index sub-column for z-axis
        filter_dict[mesh_key, 'z'] = _repeat_and_tile(
            np.arange(1, nz + 1), nx * ny * stride, data_size)

        # Initialize a Pandas DataFrame from the mesh dictionary
        df = pd.concat([df, pd.DataFrame(filter_dict)])

        return df

    def to_xml_element(self):
        """Return XML Element representing the Filter.

        Returns
        -------
        element : xml.etree.ElementTree.Element
            XML element containing filter data

        """
        element = super().to_xml_element()
        element[0].text = str(self.mesh.id)
        return element


class MeshSurfaceFilter(MeshFilter):
    """Filter events by surface crossings on a regular, rectangular mesh.

    Parameters
    ----------
    mesh : openmc.MeshBase
        The mesh object that events will be tallied onto
    filter_id : int
        Unique identifier for the filter

    Attributes
    ----------
    bins : Integral
        The mesh ID
    mesh : openmc.MeshBase
        The mesh object that events will be tallied onto
    id : int
        Unique identifier for the filter
    bins : list of tuple

        A list of mesh indices / surfaces for each filter bin, e.g. [(1, 1,
        'x-min out'), (1, 1, 'x-min in'), ...]

    num_bins : Integral
        The number of filter bins

    """

    @MeshFilter.mesh.setter
    def mesh(self, mesh):
        cv.check_type('filter mesh', mesh, openmc.MeshBase)
        self._mesh = mesh

        # Take the product of mesh indices and current names
        n_dim = mesh.n_dimension
        self.bins = [mesh_tuple + (surf,) for mesh_tuple, surf in
                     product(mesh.indices, _CURRENT_NAMES[:4*n_dim])]

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

        # Append mesh ID as outermost index of multi-index
        mesh_key = 'mesh {}'.format(self.mesh.id)

        # Find mesh dimensions - use 3D indices for simplicity
        n_surfs = 4 * len(self.mesh.dimension)
        if len(self.mesh.dimension) == 3:
            nx, ny, nz = self.mesh.dimension
        elif len(self.mesh.dimension) == 2:
            nx, ny = self.mesh.dimension
            nz = 1
        else:
            nx = self.mesh.dimension
            ny = nz = 1

        # Generate multi-index sub-column for x-axis
        filter_dict[mesh_key, 'x'] = _repeat_and_tile(
            np.arange(1, nx + 1), n_surfs * stride, data_size)

        # Generate multi-index sub-column for y-axis
        if len(self.mesh.dimension) > 1:
            filter_dict[mesh_key, 'y'] = _repeat_and_tile(
                np.arange(1, ny + 1), n_surfs * nx * stride, data_size)

        # Generate multi-index sub-column for z-axis
        if len(self.mesh.dimension) > 2:
            filter_dict[mesh_key, 'z'] = _repeat_and_tile(
                np.arange(1, nz + 1), n_surfs * nx * ny * stride, data_size)

        # Generate multi-index sub-column for surface
        filter_dict[mesh_key, 'surf'] = _repeat_and_tile(
            _CURRENT_NAMES[:n_surfs], stride, data_size)

        # Initialize a Pandas DataFrame from the mesh dictionary
        return pd.concat([df, pd.DataFrame(filter_dict)])


class RealFilter(Filter):
    """Tally modifier that describes phase-space and other characteristics

    Parameters
    ----------
    values : iterable of float
        A list of values for which each successive pair constitutes a range of
        values for a single bin
    filter_id : int
        Unique identifier for the filter

    Attributes
    ----------
    values : numpy.ndarray
        An array of values for which each successive pair constitutes a range of
        values for a single bin
    id : int
        Unique identifier for the filter
    bins : numpy.ndarray
        An array of shape (N, 2) where each row is a pair of values indicating a
        filter bin range
    num_bins : int
        The number of filter bins

    """
    def __init__(self, values, filter_id=None):
        self.values = np.asarray(values)
        self.bins = np.vstack((self.values[:-1], self.values[1:])).T
        self.id = filter_id

    def __gt__(self, other):
        if type(self) is type(other):
            # Compare largest/smallest bin edges in filters
            # This logic is used when merging tallies with real filters
            return self.values[0] >= other.values[-1]
        else:
            return super().__gt__(other)

    def __repr__(self):
        string = type(self).__name__ + '\n'
        string += '{: <16}=\t{}\n'.format('\tValues', self.values)
        string += '{: <16}=\t{}\n'.format('\tID', self.id)
        return string

    @Filter.bins.setter
    def bins(self, bins):
        Filter.bins.__set__(self, np.asarray(bins))

    def check_bins(self, bins):
        for v0, v1 in bins:
            # Values should be real
            cv.check_type('filter value', v0, Real)
            cv.check_type('filter value', v1, Real)

            # Make sure that each tuple has values that are increasing
            if v1 < v0:
                raise ValueError('Values {} and {} appear to be out of order'
                                 .format(v0, v1))

        for pair0, pair1 in zip(bins[:-1], bins[1:]):
            # Successive pairs should be ordered
            if pair1[1] < pair0[1]:
                raise ValueError('Values {} and {} appear to be out of order'
                                 .format(pair1[1], pair0[1]))

    def can_merge(self, other):
        if type(self) is not type(other):
            return False

        if self.bins[0, 0] == other.bins[-1][1]:
            # This low edge coincides with other's high edge
            return True
        elif self.bins[-1][1] == other.bins[0, 0]:
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
        merged_values = np.concatenate((self.values, other.values))
        merged_values = np.unique(merged_values)

        # Create a new filter with these bins and a new auto-generated ID
        return type(self)(sorted(merged_values))

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
        elif self.num_bins != other.num_bins:
            return False
        else:
            return np.allclose(self.values, other.values)

    def get_bin_index(self, filter_bin):
        i = np.where(self.bins[:, 1] == filter_bin[1])[0]
        if len(i) == 0:
            msg = 'Unable to get the bin index for Filter since "{0}" ' \
                  'is not one of the bins'.format(filter_bin)
            raise ValueError(msg)
        else:
            return i[0]

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
        lo_bins = np.repeat(self.bins[:, 0], stride)
        hi_bins = np.repeat(self.bins[:, 1], stride)
        tile_factor = data_size // len(lo_bins)
        lo_bins = np.tile(lo_bins, tile_factor)
        hi_bins = np.tile(hi_bins, tile_factor)

        # Add the new energy columns to the DataFrame.
        if hasattr(self, 'units'):
            units = ' [{}]'.format(self.units)
        else:
            units = ''

        df.loc[:, self.short_name.lower() + ' low' + units] = lo_bins
        df.loc[:, self.short_name.lower() + ' high' + units] = hi_bins

        return df

    def to_xml_element(self):
        """Return XML Element representing the Filter.

        Returns
        -------
        element : xml.etree.ElementTree.Element
            XML element containing filter data

        """
        element = super().to_xml_element()
        element[0].text = ' '.join(str(x) for x in self.values)
        return element


class EnergyFilter(RealFilter):
    """Bins tally events based on incident particle energy.

    Parameters
    ----------
    values : Iterable of Real
        A list of values for which each successive pair constitutes a range of
        energies in [eV] for a single bin
    filter_id : int
        Unique identifier for the filter

    Attributes
    ----------
    values : numpy.ndarray
        An array of values for which each successive pair constitutes a range of
        energies in [eV] for a single bin
    id : int
        Unique identifier for the filter
    bins : numpy.ndarray
        An array of shape (N, 2) where each row is a pair of energies in [eV]
        for a single filter bin
    num_bins : int
        The number of filter bins

    """
    units = 'eV'

    def get_bin_index(self, filter_bin):
        # Use lower energy bound to find index for RealFilters
        deltas = np.abs(self.bins[:, 1] - filter_bin[1]) / filter_bin[1]
        min_delta = np.min(deltas)
        if min_delta < 1E-3:
            return deltas.argmin()
        else:
            msg = 'Unable to get the bin index for Filter since "{0}" ' \
                  'is not one of the bins'.format(filter_bin)
            raise ValueError(msg)

    def check_bins(self, bins):
        super().check_bins(bins)
        for v0, v1 in bins:
            cv.check_greater_than('filter value', v0, 0., equality=True)
            cv.check_greater_than('filter value', v1, 0., equality=True)


class EnergyoutFilter(EnergyFilter):
    """Bins tally events based on outgoing particle energy.

    Parameters
    ----------
    values : Iterable of Real
        A list of values for which each successive pair constitutes a range of
        energies in [eV] for a single bin
    filter_id : int
        Unique identifier for the filter

    Attributes
    ----------
    values : numpy.ndarray
        An array of values for which each successive pair constitutes a range of
        energies in [eV] for a single bin
    id : int
        Unique identifier for the filter
    bins : numpy.ndarray
        An array of shape (N, 2) where each row is a pair of energies in [eV]
        for a single filter bin
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

    This filter provides a separate score for each unique instance of a repeated
    cell in a geometry. Note that only one cell can be specified in this filter.
    The related :class:`CellInstanceFilter` allows one to obtain scores for
    particular cell instances as well as instances from different cells.

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

    See Also
    --------
    CellInstanceFilter

    """

    def __init__(self, cell, filter_id=None):
        self._paths = None
        super().__init__(cell, filter_id)

    @classmethod
    def from_hdf5(cls, group, **kwargs):
        if group['type'][()].decode() != cls.short_name.lower():
            raise ValueError("Expected HDF5 data for filter type '"
                             + cls.short_name.lower() + "' but got '"
                             + group['type'][()].decode() + " instead")

        filter_id = int(group.name.split('/')[-1].lstrip('filter '))

        out = cls(group['bins'][()], filter_id=filter_id)
        out._num_bins = group['n_bins'][()]

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
            A Pandas DataFrame with columns describing distributed cells. The
            dataframe will have either:

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
                    level_dict[level_key] = _repeat_and_tile(
                        level_bins, stride, data_size)

                # Initialize a Pandas DataFrame from the level dictionary
                if level_df is None:
                    level_df = pd.DataFrame(level_dict)
                else:
                    level_df = pd.concat([level_df, pd.DataFrame(level_dict)],
                                         axis=1)

        # Create DataFrame column for distribcell instance IDs
        # NOTE: This is performed regardless of whether the user
        # requests Summary geometric information
        filter_bins = _repeat_and_tile(
            np.arange(self.num_bins), stride, data_size)
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
    values : int or Iterable of Real
        A grid of scattering angles which events will binned into. Values
        represent the cosine of the scattering angle. If an iterable is given,
        the values will be used explicitly as grid points. If a single int is
        given, the range [-1, 1] will be divided up equally into that number of
        bins.
    filter_id : int
        Unique identifier for the filter

    Attributes
    ----------
    values : numpy.ndarray
        An array of values for which each successive pair constitutes a range of
        scattering angle cosines for a single bin
    id : int
        Unique identifier for the filter
    bins : numpy.ndarray
        An array of shape (N, 2) where each row is a pair of scattering angle
        cosines for a single filter bin
    num_bins : Integral
        The number of filter bins

    """
    def __init__(self, values, filter_id=None):
        if isinstance(values, Integral):
            values = np.linspace(-1., 1., values + 1)
        super().__init__(values, filter_id)

    def check_bins(self, bins):
        super().check_bins(bins)
        for x in np.ravel(bins):
            if not np.isclose(x, -1.):
                cv.check_greater_than('filter value', x, -1., equality=True)
            if not np.isclose(x, 1.):
                cv.check_less_than('filter value', x, 1., equality=True)


class PolarFilter(RealFilter):
    """Bins tally events based on the incident particle's direction.

    Parameters
    ----------
    values : int or Iterable of Real
        A grid of polar angles which events will binned into. Values represent
        an angle in radians relative to the z-axis. If an iterable is given, the
        values will be used explicitly as grid points. If a single int is given,
        the range [0, pi] will be divided up equally into that number of bins.
    filter_id : int
        Unique identifier for the filter

    Attributes
    ----------
    values : numpy.ndarray
        An array of values for which each successive pair constitutes a range of
        polar angles in [rad] for a single bin
    id : int
        Unique identifier for the filter
    bins : numpy.ndarray
        An array of shape (N, 2) where each row is a pair of polar angles for a
        single filter bin
    id : int
        Unique identifier for the filter
    num_bins : Integral
        The number of filter bins

    """
    units = 'rad'

    def __init__(self, values, filter_id=None):
        if isinstance(values, Integral):
            values = np.linspace(0., np.pi, values + 1)
        super().__init__(values, filter_id)

    def check_bins(self, bins):
        super().check_bins(bins)
        for x in np.ravel(bins):
            if not np.isclose(x, 0.):
                cv.check_greater_than('filter value', x, 0., equality=True)
            if not np.isclose(x, np.pi):
                cv.check_less_than('filter value', x, np.pi, equality=True)


class AzimuthalFilter(RealFilter):
    """Bins tally events based on the incident particle's direction.

    Parameters
    ----------
    values : int or Iterable of Real
        A grid of azimuthal angles which events will binned into. Values
        represent an angle in radians relative to the x-axis and perpendicular
        to the z-axis. If an iterable is given, the values will be used
        explicitly as grid points. If a single int is given, the range
        [-pi, pi) will be divided up equally into that number of bins.
    filter_id : int
        Unique identifier for the filter

    Attributes
    ----------
    values : numpy.ndarray
        An array of values for which each successive pair constitutes a range of
        azimuthal angles in [rad] for a single bin
    id : int
        Unique identifier for the filter
    bins : numpy.ndarray
        An array of shape (N, 2) where each row is a pair of azimuthal angles
        for a single filter bin
    num_bins : Integral
        The number of filter bins

    """
    units = 'rad'

    def __init__(self, values, filter_id=None):
        if isinstance(values, Integral):
            values = np.linspace(-np.pi, np.pi, values + 1)
        super().__init__(values, filter_id)

    def check_bins(self, bins):
        super().check_bins(bins)
        for x in np.ravel(bins):
            if not np.isclose(x, -np.pi):
                cv.check_greater_than('filter value', x, -np.pi, equality=True)
            if not np.isclose(x, np.pi):
                cv.check_less_than('filter value', x, np.pi, equality=True)


class DelayedGroupFilter(Filter):
    """Bins fission events based on the produced neutron precursor groups.

    Parameters
    ----------
    bins : iterable of int
        The delayed neutron precursor groups.  For example, ENDF/B-VII.1 uses
        6 precursor groups so a tally with all groups will have bins =
        [1, 2, 3, 4, 5, 6].
    filter_id : int
        Unique identifier for the filter

    Attributes
    ----------
    bins : iterable of int
        The delayed neutron precursor groups.  For example, ENDF/B-VII.1 uses
        6 precursor groups so a tally with all groups will have bins =
        [1, 2, 3, 4, 5, 6].
    id : int
        Unique identifier for the filter
    num_bins : Integral
        The number of filter bins

    """
    def check_bins(self, bins):
        # Check the bin values.
        for g in bins:
            cv.check_greater_than('delayed group', g, 0)


class EnergyFunctionFilter(Filter):
    """Multiplies tally scores by an arbitrary function of incident energy.

    The arbitrary function is described by a piecewise linear-linear
    interpolation of energy and y values.  Values outside of the given energy
    range will be evaluated as zero.

    Parameters
    ----------
    energy : Iterable of Real
        A grid of energy values in [eV]
    y : iterable of Real
        A grid of interpolant values in [eV]
    filter_id : int
        Unique identifier for the filter

    Attributes
    ----------
    energy : Iterable of Real
        A grid of energy values in [eV]
    y : iterable of Real
        A grid of interpolant values in [eV]
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
        else:
            return all(self.y == other.y)

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
        if group['type'][()].decode() != cls.short_name.lower():
            raise ValueError("Expected HDF5 data for filter type '"
                             + cls.short_name.lower() + "' but got '"
                             + group['type'][()].decode() + " instead")

        energy = group['energy'][()]
        y = group['y'][()]
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
        """Return XML Element representing the Filter.

        Returns
        -------
        element : xml.etree.ElementTree.Element
            XML element containing filter data

        """
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

        filter_bins = _repeat_and_tile(out, stride, data_size)
        df = pd.concat([df, pd.DataFrame(
            {self.short_name.lower(): filter_bins})])

        return df
