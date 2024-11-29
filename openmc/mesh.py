from __future__ import annotations
import warnings
from abc import ABC, abstractmethod
from collections.abc import Iterable, Sequence
from functools import wraps
from math import pi, sqrt, atan2
from numbers import Integral, Real

import h5py
import lxml.etree as ET
import numpy as np

import openmc
import openmc.checkvalue as cv
from openmc.checkvalue import PathLike
from openmc.utility_funcs import change_directory
from ._xml import get_text
from .mixin import IDManagerMixin
from .surface import _BOUNDARY_TYPES
from .utility_funcs import input_path


class MeshBase(IDManagerMixin, ABC):
    """A mesh that partitions geometry for tallying purposes.

    Parameters
    ----------
    mesh_id : int
        Unique identifier for the mesh
    name : str
        Name of the mesh

    Attributes
    ----------
    id : int
        Unique identifier for the mesh
    name : str
        Name of the mesh
    bounding_box : openmc.BoundingBox
        Axis-aligned bounding box of the mesh as defined by the upper-right and
        lower-left coordinates.
    indices : Iterable of tuple
        An iterable of mesh indices for each mesh element, e.g. [(1, 1, 1), (2, 1, 1), ...]
    """

    next_id = 1
    used_ids = set()

    def __init__(self, mesh_id: int | None = None, name: str = ''):
        # Initialize Mesh class attributes
        self.id = mesh_id
        self.name = name

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, name: str):
        if name is not None:
            cv.check_type(f'name for mesh ID="{self._id}"', name, str)
            self._name = name
        else:
            self._name = ''

    @property
    def bounding_box(self) -> openmc.BoundingBox:
        return openmc.BoundingBox(self.lower_left, self.upper_right)

    @property
    @abstractmethod
    def indices(self):
        pass

    def __repr__(self):
        string = type(self).__name__ + '\n'
        string += '{0: <16}{1}{2}\n'.format('\tID', '=\t', self._id)
        string += '{0: <16}{1}{2}\n'.format('\tName', '=\t', self._name)
        return string

    def _volume_dim_check(self):
        if self.n_dimension != 3 or \
           any([d == 0 for d in self.dimension]):
            raise RuntimeError(f'Mesh {self.id} is not 3D. '
                               'Volumes cannot be provided.')

    @classmethod
    def from_hdf5(cls, group: h5py.Group):
        """Create mesh from HDF5 group

        Parameters
        ----------
        group : h5py.Group
            Group in HDF5 file

        Returns
        -------
        openmc.MeshBase
            Instance of a MeshBase subclass

        """

        mesh_type = group['type'][()].decode()
        if mesh_type == 'regular':
            return RegularMesh.from_hdf5(group)
        elif mesh_type == 'rectilinear':
            return RectilinearMesh.from_hdf5(group)
        elif mesh_type == 'cylindrical':
            return CylindricalMesh.from_hdf5(group)
        elif mesh_type == 'spherical':
            return SphericalMesh.from_hdf5(group)
        elif mesh_type == 'unstructured':
            return UnstructuredMesh.from_hdf5(group)
        else:
            raise ValueError('Unrecognized mesh type: "' + mesh_type + '"')

    def to_xml_element(self):
        """Return XML representation of the mesh

        Returns
        -------
        element : lxml.etree._Element
            XML element containing mesh data

        """
        elem = ET.Element("mesh")

        elem.set("id", str(self._id))
        if self.name:
            elem.set("name", self.name)

        return elem

    @classmethod
    def from_xml_element(cls, elem: ET.Element):
        """Generates a mesh from an XML element

        Parameters
        ----------
        elem : lxml.etree._Element
            XML element

        Returns
        -------
        openmc.MeshBase
            an openmc mesh object

        """
        mesh_type = get_text(elem, 'type')

        if mesh_type == 'regular' or mesh_type is None:
            return RegularMesh.from_xml_element(elem)
        elif mesh_type == 'rectilinear':
            return RectilinearMesh.from_xml_element(elem)
        elif mesh_type == 'cylindrical':
            return CylindricalMesh.from_xml_element(elem)
        elif mesh_type == 'spherical':
            return SphericalMesh.from_xml_element(elem)
        elif mesh_type == 'unstructured':
            return UnstructuredMesh.from_xml_element(elem)
        else:
            raise ValueError(f'Unrecognized mesh type "{mesh_type}" found.')

    def get_homogenized_materials(
            self,
            model: openmc.Model,
            n_samples: int = 10_000,
            prn_seed: int | None = None,
            include_void: bool = True,
            **kwargs
    ) -> list[openmc.Material]:
        """Generate homogenized materials over each element in a mesh.

        .. versionadded:: 0.15.0

        Parameters
        ----------
        model : openmc.Model
            Model containing materials to be homogenized and the associated
            geometry.
        n_samples : int
            Number of samples in each mesh element.
        prn_seed : int, optional
            Pseudorandom number generator (PRNG) seed; if None, one will be
            generated randomly.
        include_void : bool, optional
            Whether homogenization should include voids.
        **kwargs
            Keyword-arguments passed to :func:`openmc.lib.init`.

        Returns
        -------
        list of openmc.Material
            Homogenized material in each mesh element

        """
        import openmc.lib

        with change_directory(tmpdir=True):
            # In order to get mesh into model, we temporarily replace the
            # tallies with a single mesh tally using the current mesh
            original_tallies = model.tallies
            new_tally = openmc.Tally()
            new_tally.filters = [openmc.MeshFilter(self)]
            new_tally.scores = ['flux']
            model.tallies = [new_tally]

            # Export model to XML
            model.export_to_model_xml()

            # Get material volume fractions
            openmc.lib.init(**kwargs)
            mesh = openmc.lib.tallies[new_tally.id].filters[0].mesh
            mat_volume_by_element = [
                [
                    (mat.id if mat is not None else None, volume)
                    for mat, volume in mat_volume_list
                ]
                for mat_volume_list in mesh.material_volumes(n_samples, prn_seed)
            ]
            openmc.lib.finalize()

            # Restore original tallies
            model.tallies = original_tallies

        # Create homogenized material for each element
        materials = model.geometry.get_all_materials()

        # Account for materials in DAGMC universes
        # TODO: This should really get incorporated in lower-level calls to
        # get_all_materials, but right now it requires information from the
        # Model object
        for cell in model.geometry.get_all_cells().values():
            if isinstance(cell.fill, openmc.DAGMCUniverse):
                names = cell.fill.material_names
                materials.update({
                    mat.id: mat for mat in model.materials if mat.name in names
                })

        homogenized_materials = []
        for mat_volume_list in mat_volume_by_element:
            material_ids, volumes = [list(x) for x in zip(*mat_volume_list)]
            total_volume = sum(volumes)

            # Check for void material and remove
            try:
                index_void = material_ids.index(None)
            except ValueError:
                pass
            else:
                material_ids.pop(index_void)
                volumes.pop(index_void)

            # If void should be excluded, adjust total volume
            if not include_void:
                total_volume = sum(volumes)

            # Compute volume fractions
            volume_fracs = np.array(volumes) / total_volume

            # Get list of materials and mix 'em up!
            mats = [materials[uid] for uid in material_ids]
            homogenized_mat = openmc.Material.mix_materials(
                mats, volume_fracs, 'vo'
            )
            homogenized_mat.volume = total_volume
            homogenized_materials.append(homogenized_mat)

        return homogenized_materials


class StructuredMesh(MeshBase):
    """A base class for structured mesh functionality

    Parameters
    ----------
    mesh_id : int
        Unique identifier for the mesh
    name : str
        Name of the mesh

    Attributes
    ----------
    id : int
        Unique identifier for the mesh
    name : str
        Name of the mesh

    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    @property
    @abstractmethod
    def dimension(self):
        pass

    @property
    @abstractmethod
    def n_dimension(self):
        pass

    @property
    @abstractmethod
    def _grids(self):
        pass

    @property
    def vertices(self):
        """Return coordinates of mesh vertices in Cartesian coordinates. Also
           see :meth:`CylindricalMesh.vertices_cylindrical` and
           :meth:`SphericalMesh.vertices_spherical` for coordinates in other coordinate
           systems.

        Returns
        -------
        vertices : numpy.ndarray
            Returns a numpy.ndarray representing the coordinates of the mesh
            vertices with a shape equal to (dim1 + 1, ..., dimn + 1, ndim). X, Y, Z values
            can be unpacked with xx, yy, zz = np.rollaxis(mesh.vertices, -1).

        """
        return self._generate_vertices(*self._grids)

    @staticmethod
    def _generate_vertices(i_grid, j_grid, k_grid):
        """Returns an array with shape (i_grid.size, j_grid.size, k_grid.size, 3)
           containing the corner vertices of mesh elements.
        """
        return np.stack(np.meshgrid(i_grid, j_grid, k_grid, indexing='ij'), axis=-1)

    @staticmethod
    def _generate_edge_midpoints(grids):
        """Generates the midpoints of mesh element edges for each dimension of the mesh.

        Parameters
        ----------
        grids : numpy.ndarray
            The vertex grids along each dimension of the mesh.

        Returns
        -------
        midpoint_grids : list of numpy.ndarray
            The edge midpoints for the i, j, and k midpoints of each element in
            i, j, k ordering. The shapes of the resulting grids are
            [(ni-1, nj, nk, 3), (ni, nj-1, nk, 3), (ni, nj, nk-1, 3)]
        """
        # generate a set of edge midpoints for each dimension
        midpoint_grids = []
        # generate the element edge midpoints in order s.t.
        # the epxected element ordering is preserved with respect to the corner vertices

        # each grid is comprised of the mid points for one dimension and the
        # corner vertices of the other two
        for dims in ((0, 1, 2), (1, 0, 2), (2, 0, 1)):
            # compute the midpoints along the last dimension
            midpoints = grids[dims[0]][:-1] + 0.5 * np.diff(grids[dims[0]])

            coords = (midpoints, grids[dims[1]], grids[dims[2]])

            i_grid, j_grid, k_grid = [coords[dims.index(i)] for i in range(3)]

            # re-use the generate vertices method to create the full mesh grid
            # transpose to get (i, j, k) ordering of the gridpoints
            midpoint_grid = StructuredMesh._generate_vertices(i_grid, j_grid, k_grid)
            midpoint_grids.append(midpoint_grid)

        return midpoint_grids

    @property
    def midpoint_vertices(self):
        """Create vertices that lie on the midpoint of element edges
        """
        # generate edge midpoints needed for curvilinear element definition
        midpoint_vertices = self._generate_edge_midpoints(self._grids)

        # convert each of the midpoint grids to cartesian coordinates
        for vertices in midpoint_vertices:
            self._convert_to_cartesian(vertices, self.origin)

        return midpoint_vertices

    @property
    def centroids(self):
        """Return coordinates of mesh element centroids.

        Returns
        -------
        centroids : numpy.ndarray
            Returns a numpy.ndarray representing the mesh element centroid
            coordinates with a shape equal to (dim1, ..., dimn, ndim). X,
            Y, Z values can be unpacked with xx, yy, zz =
            np.rollaxis(mesh.centroids, -1).
        """
        ndim = self.n_dimension
        # this line ensures that the vertices aren't adjusted by the origin or
        # converted to the Cartesian system for cylindrical and spherical meshes
        vertices = StructuredMesh.vertices.fget(self)
        s0 = (slice(0, -1),)*ndim + (slice(None),)
        s1 = (slice(1, None),)*ndim + (slice(None),)
        return (vertices[s0] + vertices[s1]) / 2

    @property
    def num_mesh_cells(self):
        return np.prod(self.dimension)

    def write_data_to_vtk(self,
                          filename: PathLike,
                          datasets: dict | None = None,
                          volume_normalization: bool = True,
                          curvilinear: bool = False):
        """Creates a VTK object of the mesh

        Parameters
        ----------
        filename : str
            Name of the VTK file to write.
        datasets : dict
            Dictionary whose keys are the data labels
            and values are the data sets.
        volume_normalization : bool, optional
            Whether or not to normalize the data by
            the volume of the mesh elements.
        curvilinear : bool
            Whether or not to write curvilinear elements. Only applies to
            ``SphericalMesh`` and ``CylindricalMesh``.

        Raises
        ------
        ValueError
            When the size of a dataset doesn't match the number of mesh cells

        Returns
        -------
        vtk.StructuredGrid or vtk.UnstructuredGrid
            a VTK grid object representing the mesh
        """
        import vtk
        from vtk.util import numpy_support as nps

        # check that the data sets are appropriately sized
        if datasets is not None:
            self._check_vtk_datasets(datasets)

        # write linear elements using a structured grid
        if not curvilinear or isinstance(self, (RegularMesh, RectilinearMesh)):
            vtk_grid = self._create_vtk_structured_grid()
            writer = vtk.vtkStructuredGridWriter()
        # write curvilinear elements using an unstructured grid
        else:
            vtk_grid = self._create_vtk_unstructured_grid()
            writer = vtk.vtkUnstructuredGridWriter()

        if datasets is not None:
            # maintain a list of the datasets as added
            # to the VTK arrays to ensure they persist
            # in memory until the file is written
            datasets_out = []
            for label, dataset in datasets.items():
                dataset = np.asarray(dataset).flatten()
                datasets_out.append(dataset)

                if volume_normalization:
                    dataset /= self.volumes.T.flatten()

                dataset_array = vtk.vtkDoubleArray()
                dataset_array.SetName(label)
                dataset_array.SetArray(nps.numpy_to_vtk(dataset),
                                    dataset.size,
                                    True)
                vtk_grid.GetCellData().AddArray(dataset_array)

        writer.SetFileName(str(filename))
        writer.SetInputData(vtk_grid)
        writer.Write()

        return vtk_grid

    def _create_vtk_structured_grid(self):
        """Create a structured grid

        Returns
        -------
        vtk.vtkStructuredGrid
            a VTK structured grid object representing the mesh
        """
        import vtk
        from vtk.util import numpy_support as nps

        vtkPts = vtk.vtkPoints()
        vtkPts.SetData(nps.numpy_to_vtk(np.swapaxes(self.vertices, 0, 2).reshape(-1, 3), deep=True))
        vtk_grid = vtk.vtkStructuredGrid()
        vtk_grid.SetPoints(vtkPts)
        vtk_grid.SetDimensions(*[dim + 1 for dim in self.dimension])

        return vtk_grid

    def _create_vtk_unstructured_grid(self):
        """Create an unstructured grid of curvilinear elements
           representing the mesh

        Returns
        -------
        vtk.vtkUnstructuredGrid
            a VTK unstructured grid object representing the mesh
        """
        import vtk
        from vtk.util import numpy_support as nps

        corner_vertices = np.swapaxes(self.vertices, 0, 2).reshape(-1, 3)

        vtkPts = vtk.vtkPoints()
        vtk_grid = vtk.vtkUnstructuredGrid()
        vtk_grid.SetPoints(vtkPts)
        # add corner vertices to the point set for the unstructured grid
        # only insert unique points, we'll get their IDs in the point set to
        # define element connectivity later
        vtkPts.SetData(nps.numpy_to_vtk(np.unique(corner_vertices, axis=0), deep=True))

        # create a locator to assist with duplicate points
        locator = vtk.vtkPointLocator()
        locator.SetDataSet(vtk_grid)
        locator.AutomaticOn() # autmoatically adds points to locator
        locator.InitPointInsertion(vtkPts, vtkPts.GetBounds())
        locator.BuildLocator()

        # this function is used to add new points to the unstructured
        # grid. It will return an existing point ID if the point is alread present
        def _insert_point(pnt):
            result = locator.IsInsertedPoint(pnt)
            if result == -1:
                point_id = vtkPts.InsertNextPoint(pnt)
                locator.InsertPoint(point_id, pnt)
                return point_id
            else:
                return result

        # Add all points to the unstructured grid, maintaining a flat list of IDs as we go ###

        # flat array storing point IDs for a given vertex
        # in the grid
        point_ids = []

        # add element corner vertices to array
        for pnt in corner_vertices:
            point_ids.append(_insert_point(pnt))

        # get edge midpoints and add them to the
        # list of point IDs
        midpoint_vertices = self.midpoint_vertices
        for edge_grid in midpoint_vertices:
            for pnt in np.swapaxes(edge_grid, 0, 2).reshape(-1, 3):
                point_ids.append(_insert_point(pnt))

        # determine how many elements in each dimension
        # and how many points in each grid
        n_elem = np.asarray(self.dimension)
        n_pnts = n_elem + 1

        # create hexes and set points for corner
        # vertices
        for i, j, k in self.indices:
            # handle indices indexed from one
            i -= 1
            j -= 1
            k -= 1

            # create a new vtk hex
            hex = vtk.vtkQuadraticHexahedron()

            # set connectivity the hex corners
            for n, (di, dj, dk) in enumerate(_HEX_VERTEX_CONN):
                # compute flat index into the point ID list based on i, j, k
                # of the vertex
                flat_idx = np.ravel_multi_index((i+di, j+dj, k+dk), n_pnts, order='F')
                # set corner vertices
                hex.GetPointIds().SetId(n, point_ids[flat_idx])

            # set connectivity of the hex midpoints
            n_midpoint_vertices = [v.size // 3 for v in midpoint_vertices]
            for n, (dim, (di, dj, dk)) in enumerate(_HEX_MIDPOINT_CONN):
                # initial offset for corner vertices and midpoint dimension
                flat_idx = corner_vertices.shape[0] + sum(n_midpoint_vertices[:dim])
                # generate a flat index into the table of point IDs
                midpoint_shape = midpoint_vertices[dim].shape[:-1]
                flat_idx += np.ravel_multi_index((i+di, j+dj, k+dk),
                                                 midpoint_shape,
                                                 order='F')
                # set hex midpoint connectivity
                hex.GetPointIds().SetId(_N_HEX_VERTICES + n, point_ids[flat_idx])

            # add the hex to the grid
            vtk_grid.InsertNextCell(hex.GetCellType(), hex.GetPointIds())

        return vtk_grid

    def _check_vtk_datasets(self, datasets: dict):
        """Perform some basic checks that the datasets are valid for this mesh

        Parameters
        ----------
        datasets : dict
            Dictionary whose keys are the data labels
            and values are the data sets.

        """
        for label, dataset in datasets.items():
            errmsg = (
                f"The size of the dataset '{label}' ({dataset.size}) should be"
                f" equal to the number of mesh cells ({self.num_mesh_cells})"
            )
            if isinstance(dataset, np.ndarray):
                if not dataset.size == self.num_mesh_cells:
                    raise ValueError(errmsg)
            else:
                if len(dataset) == self.num_mesh_cells:
                    raise ValueError(errmsg)
            cv.check_type('data label', label, str)


class RegularMesh(StructuredMesh):
    """A regular Cartesian mesh in one, two, or three dimensions

    Parameters
    ----------
    mesh_id : int
        Unique identifier for the mesh
    name : str
        Name of the mesh

    Attributes
    ----------
    id : int
        Unique identifier for the mesh
    name : str
        Name of the mesh
    dimension : Iterable of int
        The number of mesh cells in each direction (x, y, z).
    n_dimension : int
        Number of mesh dimensions.
    lower_left : Iterable of float
        The lower-left corner of the structured mesh. If only two coordinate
        are given, it is assumed that the mesh is an x-y mesh.
    upper_right : Iterable of float
        The upper-right corner of the structured mesh. If only two coordinate
        are given, it is assumed that the mesh is an x-y mesh.
    bounding_box : openmc.BoundingBox
        Axis-aligned bounding box of the mesh as defined by the upper-right and
        lower-left coordinates.
    width : Iterable of float
        The width of mesh cells in each direction.
    indices : Iterable of tuple
        An iterable of mesh indices for each mesh element, e.g. [(1, 1, 1),
        (2, 1, 1), ...]

    """

    def __init__(self, mesh_id: int | None = None, name: str = ''):
        super().__init__(mesh_id, name)

        self._dimension = None
        self._lower_left = None
        self._upper_right = None
        self._width = None

    @property
    def dimension(self):
        return tuple(self._dimension)

    @dimension.setter
    def dimension(self, dimension: Iterable[int]):
        cv.check_type('mesh dimension', dimension, Iterable, Integral)
        cv.check_length('mesh dimension', dimension, 1, 3)
        self._dimension = dimension

    @property
    def n_dimension(self):
        if self._dimension is not None:
            return len(self._dimension)
        else:
            return None

    @property
    def lower_left(self):
        return self._lower_left

    @lower_left.setter
    def lower_left(self, lower_left: Iterable[Real]):
        cv.check_type('mesh lower_left', lower_left, Iterable, Real)
        cv.check_length('mesh lower_left', lower_left, 1, 3)
        self._lower_left = lower_left

        if self.upper_right is not None and any(np.isclose(self.upper_right, lower_left)):
            raise ValueError("Mesh cannot have zero thickness in any dimension")

    @property
    def upper_right(self):
        if self._upper_right is not None:
            return self._upper_right
        elif self._width is not None:
            if self._lower_left is not None and self._dimension is not None:
                ls = self._lower_left
                ws = self._width
                dims = self._dimension
                return [l + w * d for l, w, d in zip(ls, ws, dims)]

    @upper_right.setter
    def upper_right(self, upper_right: Iterable[Real]):
        cv.check_type('mesh upper_right', upper_right, Iterable, Real)
        cv.check_length('mesh upper_right', upper_right, 1, 3)
        self._upper_right = upper_right

        if self._width is not None:
            self._width = None
            warnings.warn("Unsetting width attribute.")

        if self.lower_left is not None and any(np.isclose(self.lower_left, upper_right)):
            raise ValueError("Mesh cannot have zero thickness in any dimension")

    @property
    def width(self):
        if self._width is not None:
            return self._width
        elif self._upper_right is not None:
            if self._lower_left is not None and self._dimension is not None:
                us = self._upper_right
                ls = self._lower_left
                dims =  self._dimension
                return [(u - l) / d for u, l, d in zip(us, ls, dims)]

    @width.setter
    def width(self, width: Iterable[Real]):
        cv.check_type('mesh width', width, Iterable, Real)
        cv.check_length('mesh width', width, 1, 3)
        self._width = width

        if self._upper_right is not None:
            self._upper_right = None
            warnings.warn("Unsetting upper_right attribute.")

    @property
    def volumes(self):
        """Return Volumes for every mesh cell

        Returns
        -------
        volumes : numpy.ndarray
            Volumes

        """
        self._volume_dim_check()
        return np.full(self.dimension, np.prod(self.width))

    @property
    def total_volume(self):
        return np.prod(self.dimension) * np.prod(self.width)

    @property
    def indices(self):
        ndim = len(self._dimension)
        if ndim == 3:
            nx, ny, nz = self.dimension
            return ((x, y, z)
                    for z in range(1, nz + 1)
                    for y in range(1, ny + 1)
                    for x in range(1, nx + 1))
        elif ndim == 2:
            nx, ny = self.dimension
            return ((x, y)
                    for y in range(1, ny + 1)
                    for x in range(1, nx + 1))
        else:
            nx, = self.dimension
            return ((x,) for x in range(1, nx + 1))

    @property
    def _grids(self):
        ndim = len(self._dimension)
        if ndim == 3:
            x0, y0, z0 = self.lower_left
            x1, y1, z1 = self.upper_right
            nx, ny, nz = self.dimension
            xarr = np.linspace(x0, x1, nx + 1)
            yarr = np.linspace(y0, y1, ny + 1)
            zarr = np.linspace(z0, z1, nz + 1)
            return (xarr, yarr, zarr)
        elif ndim == 2:
            x0, y0 = self.lower_left
            x1, y1 = self.upper_right
            nx, ny = self.dimension
            xarr = np.linspace(x0, x1, nx + 1)
            yarr = np.linspace(y0, y1, ny + 1)
            return (xarr, yarr)
        else:
            nx, = self.dimension
            x0, = self.lower_left
            x1, = self.upper_right
            return (np.linspace(x0, x1, nx + 1),)

    def __repr__(self):
        string = super().__repr__()
        string += '{0: <16}{1}{2}\n'.format('\tDimensions', '=\t', self.n_dimension)
        string += '{0: <16}{1}{2}\n'.format('\tVoxels', '=\t', self._dimension)
        string += '{0: <16}{1}{2}\n'.format('\tLower left', '=\t', self._lower_left)
        string += '{0: <16}{1}{2}\n'.format('\tUpper Right', '=\t', self.upper_right)
        string += '{0: <16}{1}{2}\n'.format('\tWidth', '=\t', self.width)
        return string

    @classmethod
    def from_hdf5(cls, group: h5py.Group):
        mesh_id = int(group.name.split('/')[-1].lstrip('mesh '))

        # Read and assign mesh properties
        mesh = cls(mesh_id)
        mesh.dimension = group['dimension'][()]
        mesh.lower_left = group['lower_left'][()]
        if 'width' in group:
            mesh.width = group['width'][()]
        elif 'upper_right' in group:
            mesh.upper_right = group['upper_right'][()]
        else:
            raise IOError('Invalid mesh: must have one of "upper_right" or "width"')

        return mesh

    @classmethod
    def from_rect_lattice(
        cls,
        lattice: 'openmc.RectLattice',
        division: int = 1,
        mesh_id: int | None = None,
        name: str = ''
    ):
        """Create mesh from an existing rectangular lattice

        Parameters
        ----------
        lattice : openmc.RectLattice
            Rectangular lattice used as a template for this mesh
        division : int
            Number of mesh cells per lattice cell.
            If not specified, there will be 1 mesh cell per lattice cell.
        mesh_id : int
            Unique identifier for the mesh
        name : str
            Name of the mesh

        Returns
        -------
        openmc.RegularMesh
            RegularMesh instance

        """
        cv.check_type('rectangular lattice', lattice, openmc.RectLattice)

        shape = np.array(lattice.shape)
        width = lattice.pitch*shape

        mesh = cls(mesh_id=mesh_id, name=name)
        mesh.lower_left = lattice.lower_left
        mesh.upper_right = lattice.lower_left + width
        mesh.dimension = shape*division

        return mesh

    @classmethod
    def from_domain(
        cls,
        domain: 'openmc.Cell' | 'openmc.Region' | 'openmc.Universe' | 'openmc.Geometry',
        dimension: Sequence[int] = (10, 10, 10),
        mesh_id: int | None = None,
        name: str = ''
    ):
        """Create mesh from an existing openmc cell, region, universe or
        geometry by making use of the objects bounding box property.

        Parameters
        ----------
        domain : {openmc.Cell, openmc.Region, openmc.Universe, openmc.Geometry}
            The object passed in will be used as a template for this mesh. The
            bounding box of the property of the object passed will be used to
            set the lower_left and upper_right and of the mesh instance
        dimension : Iterable of int
            The number of mesh cells in each direction (x, y, z).
        mesh_id : int
            Unique identifier for the mesh
        name : str
            Name of the mesh

        Returns
        -------
        openmc.RegularMesh
            RegularMesh instance

        """
        cv.check_type(
            "domain",
            domain,
            (openmc.Cell, openmc.Region, openmc.Universe, openmc.Geometry),
        )

        mesh = cls(mesh_id=mesh_id, name=name)
        mesh.lower_left = domain.bounding_box[0]
        mesh.upper_right = domain.bounding_box[1]
        mesh.dimension = dimension

        return mesh

    def to_xml_element(self):
        """Return XML representation of the mesh

        Returns
        -------
        element : lxml.etree._Element
            XML element containing mesh data

        """
        element = super().to_xml_element()

        if self._dimension is not None:
            subelement = ET.SubElement(element, "dimension")
            subelement.text = ' '.join(map(str, self._dimension))

        subelement = ET.SubElement(element, "lower_left")
        subelement.text = ' '.join(map(str, self._lower_left))

        if self._upper_right is not None:
            subelement = ET.SubElement(element, "upper_right")
            subelement.text = ' '.join(map(str, self._upper_right))
        if self._width is not None:
            subelement = ET.SubElement(element, "width")
            subelement.text = ' '.join(map(str, self._width))

        return element

    @classmethod
    def from_xml_element(cls, elem: ET.Element):
        """Generate mesh from an XML element

        Parameters
        ----------
        elem : lxml.etree._Element
            XML element

        Returns
        -------
        openmc.Mesh
            Mesh generated from XML element

        """
        mesh_id = int(get_text(elem, 'id'))
        mesh = cls(mesh_id=mesh_id)

        mesh_type = get_text(elem, 'type')
        if mesh_type is not None:
            mesh.type = mesh_type

        dimension = get_text(elem, 'dimension')
        if dimension is not None:
            mesh.dimension = [int(x) for x in dimension.split()]

        lower_left = get_text(elem, 'lower_left')
        if lower_left is not None:
            mesh.lower_left = [float(x) for x in lower_left.split()]

        upper_right = get_text(elem, 'upper_right')
        if upper_right is not None:
            mesh.upper_right = [float(x) for x in upper_right.split()]

        width = get_text(elem, 'width')
        if width is not None:
            mesh.width = [float(x) for x in width.split()]

        return mesh

    def build_cells(self, bc: str | None = None):
        """Generates a lattice of universes with the same dimensionality
        as the mesh object.  The individual cells/universes produced
        will not have material definitions applied and so downstream code
        will have to apply that information.

        Parameters
        ----------
        bc : iterable of {'reflective', 'periodic', 'transmission', 'vacuum', or 'white'}
            Boundary conditions for each of the four faces of a rectangle
            (if applying to a 2D mesh) or six faces of a parallelepiped
            (if applying to a 3D mesh) provided in the following order:
            [x min, x max, y min, y max, z min, z max].  2-D cells do not
            contain the z min and z max entries. Defaults to 'reflective' for
            all faces.

        Returns
        -------
        root_cell : openmc.Cell
            The cell containing the lattice representing the mesh geometry;
            this cell is a single parallelepiped with boundaries matching
            the outermost mesh boundary with the boundary conditions from bc
            applied.
        cells : iterable of openmc.Cell
            The list of cells within each lattice position mimicking the mesh
            geometry.

        """
        if bc is None:
            bc = ['reflective'] * 6
        if len(bc) not in (4, 6):
            raise ValueError('Boundary condition must be of length 4 or 6')
        for entry in bc:
            cv.check_value('bc', entry, _BOUNDARY_TYPES)

        n_dim = self.n_dimension

        # Build the cell which will contain the lattice
        xplanes = [openmc.XPlane(self.lower_left[0], boundary_type=bc[0]),
                   openmc.XPlane(self.upper_right[0], boundary_type=bc[1])]
        if n_dim == 1:
            yplanes = [openmc.YPlane(-1e10, boundary_type='reflective'),
                       openmc.YPlane(1e10, boundary_type='reflective')]
        else:
            yplanes = [openmc.YPlane(self.lower_left[1], boundary_type=bc[2]),
                       openmc.YPlane(self.upper_right[1], boundary_type=bc[3])]

        if n_dim <= 2:
            # Would prefer to have the z ranges be the max supported float, but
            # these values are apparently different between python and Fortran.
            # Choosing a safe and sane default.
            # Values of +/-1e10 are used here as there seems to be an
            # inconsistency between what numpy uses as the max float and what
            # Fortran expects for a real(8), so this avoids code complication
            # and achieves the same goal.
            zplanes = [openmc.ZPlane(-1e10, boundary_type='reflective'),
                       openmc.ZPlane(1e10, boundary_type='reflective')]
        else:
            zplanes = [openmc.ZPlane(self.lower_left[2], boundary_type=bc[4]),
                       openmc.ZPlane(self.upper_right[2], boundary_type=bc[5])]
        root_cell = openmc.Cell()
        root_cell.region = ((+xplanes[0] & -xplanes[1]) &
                            (+yplanes[0] & -yplanes[1]) &
                            (+zplanes[0] & -zplanes[1]))

        # Build the universes which will be used for each of the (i,j,k)
        # locations within the mesh.
        # We will concurrently build cells to assign to these universes
        cells = []
        universes = []
        for _ in self.indices:
            cells.append(openmc.Cell())
            universes.append(openmc.Universe())
            universes[-1].add_cell(cells[-1])

        lattice = openmc.RectLattice()
        lattice.lower_left = self.lower_left

        # Assign the universe and rotate to match the indexing expected for
        # the lattice
        if n_dim == 1:
            universe_array = np.array([universes])
        elif n_dim == 2:
            universe_array = np.empty(self.dimension[::-1],
                                      dtype=openmc.Universe)
            i = 0
            for y in range(self.dimension[1] - 1, -1, -1):
                for x in range(self.dimension[0]):
                    universe_array[y][x] = universes[i]
                    i += 1
        else:
            universe_array = np.empty(self.dimension[::-1],
                                      dtype=openmc.Universe)
            i = 0
            for z in range(self.dimension[2]):
                for y in range(self.dimension[1] - 1, -1, -1):
                    for x in range(self.dimension[0]):
                        universe_array[z][y][x] = universes[i]
                        i += 1
        lattice.universes = universe_array

        if self.width is not None:
            lattice.pitch = self.width
        else:
            dx = ((self.upper_right[0] - self.lower_left[0]) /
                  self.dimension[0])

            if n_dim == 1:
                lattice.pitch = [dx]
            elif n_dim == 2:
                dy = ((self.upper_right[1] - self.lower_left[1]) /
                      self.dimension[1])
                lattice.pitch = [dx, dy]
            else:
                dy = ((self.upper_right[1] - self.lower_left[1]) /
                      self.dimension[1])
                dz = ((self.upper_right[2] - self.lower_left[2]) /
                      self.dimension[2])
                lattice.pitch = [dx, dy, dz]

        # Fill Cell with the Lattice
        root_cell.fill = lattice

        return root_cell, cells


def Mesh(*args, **kwargs):
    warnings.warn("Mesh has been renamed RegularMesh. Future versions of "
                  "OpenMC will not accept the name Mesh.")
    return RegularMesh(*args, **kwargs)


class RectilinearMesh(StructuredMesh):
    """A 3D rectilinear Cartesian mesh

    Parameters
    ----------
    mesh_id : int
        Unique identifier for the mesh
    name : str
        Name of the mesh

    Attributes
    ----------
    id : int
        Unique identifier for the mesh
    name : str
        Name of the mesh
    dimension : Iterable of int
        The number of mesh cells in each direction (x, y, z).
    n_dimension : int
        Number of mesh dimensions (always 3 for a RectilinearMesh).
    x_grid : numpy.ndarray
        1-D array of mesh boundary points along the x-axis.
    y_grid : numpy.ndarray
        1-D array of mesh boundary points along the y-axis.
    z_grid : numpy.ndarray
        1-D array of mesh boundary points along the z-axis.
    indices : Iterable of tuple
        An iterable of mesh indices for each mesh element, e.g. [(1, 1, 1),
        (2, 1, 1), ...]
    bounding_box : openmc.BoundingBox
        Axis-aligned bounding box of the mesh as defined by the upper-right and
        lower-left coordinates.

    """

    def __init__(self, mesh_id: int = None, name: str = ''):
        super().__init__(mesh_id, name)

        self._x_grid = None
        self._y_grid = None
        self._z_grid = None

    @property
    def dimension(self):
        return (len(self.x_grid) - 1,
                len(self.y_grid) - 1,
                len(self.z_grid) - 1)

    @property
    def n_dimension(self):
        return 3

    @property
    def x_grid(self):
        return self._x_grid

    @x_grid.setter
    def x_grid(self, grid):
        cv.check_type('mesh x_grid', grid, Iterable, Real)
        self._x_grid = np.asarray(grid, dtype=float)

    @property
    def y_grid(self):
        return self._y_grid

    @y_grid.setter
    def y_grid(self, grid):
        cv.check_type('mesh y_grid', grid, Iterable, Real)
        self._y_grid = np.asarray(grid, dtype=float)

    @property
    def z_grid(self):
        return self._z_grid

    @z_grid.setter
    def z_grid(self, grid):
        cv.check_type('mesh z_grid', grid, Iterable, Real)
        self._z_grid = np.asarray(grid, dtype=float)

    @property
    def _grids(self):
        return (self.x_grid, self.y_grid, self.z_grid)

    @property
    def lower_left(self):
        return np.array([self.x_grid[0], self.y_grid[0], self.z_grid[0]])

    @property
    def upper_right(self):
        return np.array([self.x_grid[-1], self.y_grid[-1], self.z_grid[-1]])

    @property
    def volumes(self):
        """Return Volumes for every mesh cell

        Returns
        -------
        volumes : numpy.ndarray
            Volumes

        """
        self._volume_dim_check()
        V_x = np.diff(self.x_grid)
        V_y = np.diff(self.y_grid)
        V_z = np.diff(self.z_grid)

        return np.multiply.outer(np.outer(V_x, V_y), V_z)

    @property
    def total_volume(self):
        return np.sum(self.volumes)

    @property
    def indices(self):
        nx = len(self.x_grid) - 1
        ny = len(self.y_grid) - 1
        nz = len(self.z_grid) - 1
        return ((x, y, z)
                for z in range(1, nz + 1)
                for y in range(1, ny + 1)
                for x in range(1, nx + 1))

    def __repr__(self):
        fmt = '{0: <16}{1}{2}\n'
        string = super().__repr__()
        string += fmt.format('\tDimensions', '=\t', self.n_dimension)
        x_grid_str = str(self._x_grid) if self._x_grid is None else len(self._x_grid)
        string += fmt.format('\tN X pnts:', '=\t', x_grid_str)
        if self._x_grid is not None:
            string += fmt.format('\tX Min:', '=\t', self._x_grid[0])
            string += fmt.format('\tX Max:', '=\t', self._x_grid[-1])
        y_grid_str = str(self._y_grid) if self._y_grid is None else len(self._y_grid)
        string += fmt.format('\tN Y pnts:', '=\t', y_grid_str)
        if self._y_grid is not None:
            string += fmt.format('\tY Min:', '=\t', self._y_grid[0])
            string += fmt.format('\tY Max:', '=\t', self._y_grid[-1])
        z_grid_str = str(self._z_grid) if self._z_grid is None else len(self._z_grid)
        string += fmt.format('\tN Z pnts:', '=\t', z_grid_str)
        if self._z_grid is not None:
            string += fmt.format('\tZ Min:', '=\t', self._z_grid[0])
            string += fmt.format('\tZ Max:', '=\t', self._z_grid[-1])
        return string

    @classmethod
    def from_hdf5(cls, group: h5py.Group):
        mesh_id = int(group.name.split('/')[-1].lstrip('mesh '))

        # Read and assign mesh properties
        mesh = cls(mesh_id=mesh_id)
        mesh.x_grid = group['x_grid'][()]
        mesh.y_grid = group['y_grid'][()]
        mesh.z_grid = group['z_grid'][()]

        return mesh

    @classmethod
    def from_xml_element(cls, elem: ET.Element):
        """Generate a rectilinear mesh from an XML element

        Parameters
        ----------
        elem : lxml.etree._Element
            XML element

        Returns
        -------
        openmc.RectilinearMesh
            Rectilinear mesh object

        """
        mesh_id = int(get_text(elem, 'id'))
        mesh = cls(mesh_id=mesh_id)
        mesh.x_grid = [float(x) for x in get_text(elem, 'x_grid').split()]
        mesh.y_grid = [float(y) for y in get_text(elem, 'y_grid').split()]
        mesh.z_grid = [float(z) for z in get_text(elem, 'z_grid').split()]

        return mesh

    def to_xml_element(self):
        """Return XML representation of the mesh

        Returns
        -------
        element : lxml.etree._Element
            XML element containing mesh data

        """

        element = super().to_xml_element()
        element.set("type", "rectilinear")

        subelement = ET.SubElement(element, "x_grid")
        subelement.text = ' '.join(map(str, self.x_grid))

        subelement = ET.SubElement(element, "y_grid")
        subelement.text = ' '.join(map(str, self.y_grid))

        subelement = ET.SubElement(element, "z_grid")
        subelement.text = ' '.join(map(str, self.z_grid))

        return element


class CylindricalMesh(StructuredMesh):
    """A 3D cylindrical mesh

    Parameters
    ----------
    r_grid : numpy.ndarray
        1-D array of mesh boundary points along the r-axis
        Requirement is r >= 0.
    z_grid : numpy.ndarray
        1-D array of mesh boundary points along the z-axis relative to the
        origin.
    phi_grid : numpy.ndarray
        1-D array of mesh boundary points along the phi-axis in radians.
        The default value is [0, 2π], i.e. the full phi range.
    origin : numpy.ndarray
        1-D array of length 3 the (x,y,z) origin of the mesh in
        cartesian coordinates
    mesh_id : int
        Unique identifier for the mesh
    name : str
        Name of the mesh

    Attributes
    ----------
    id : int
        Unique identifier for the mesh
    name : str
        Name of the mesh
    dimension : Iterable of int
        The number of mesh cells in each direction (r_grid, phi_grid, z_grid).
    n_dimension : int
        Number of mesh dimensions (always 3 for a CylindricalMesh).
    r_grid : numpy.ndarray
        1-D array of mesh boundary points along the r-axis.
        Requirement is r >= 0.
    phi_grid : numpy.ndarray
        1-D array of mesh boundary points along the phi-axis in radians.
        The default value is [0, 2π], i.e. the full phi range.
    z_grid : numpy.ndarray
        1-D array of mesh boundary points along the z-axis relative to the
        origin.
    origin : numpy.ndarray
        1-D array of length 3 the (x,y,z) origin of the mesh in
        cartesian coordinates
    indices : Iterable of tuple
        An iterable of mesh indices for each mesh element, e.g. [(1, 1, 1),
        (2, 1, 1), ...]
    lower_left : Iterable of float
        The lower-left corner of the structured mesh. If only two coordinate
        are given, it is assumed that the mesh is an x-y mesh.
    upper_right : Iterable of float
        The upper-right corner of the structured mesh. If only two coordinate
        are given, it is assumed that the mesh is an x-y mesh.
    bounding_box : openmc.BoundingBox
        Axis-aligned bounding box of the mesh as defined by the upper-right and
        lower-left coordinates.

    """

    def __init__(
        self,
        r_grid: Sequence[float],
        z_grid: Sequence[float],
        phi_grid: Sequence[float] = (0, 2*pi),
        origin: Sequence[float] = (0., 0., 0.),
        mesh_id: int | None = None,
        name: str = '',
    ):
        super().__init__(mesh_id, name)

        self.r_grid = r_grid
        self.phi_grid = phi_grid
        self.z_grid = z_grid
        self.origin = origin

    @property
    def dimension(self):
        return (len(self.r_grid) - 1,
                len(self.phi_grid) - 1,
                len(self.z_grid) - 1)

    @property
    def n_dimension(self):
        return 3

    @property
    def origin(self):
        return self._origin

    @origin.setter
    def origin(self, coords):
        cv.check_type('mesh origin', coords, Iterable, Real)
        cv.check_length("mesh origin", coords, 3)
        self._origin = np.asarray(coords)

    @property
    def r_grid(self):
        return self._r_grid

    @r_grid.setter
    def r_grid(self, grid):
        cv.check_type('mesh r_grid', grid, Iterable, Real)
        cv.check_length('mesh r_grid', grid, 2)
        cv.check_increasing('mesh r_grid', grid)
        self._r_grid = np.asarray(grid, dtype=float)

    @property
    def phi_grid(self):
        return self._phi_grid

    @phi_grid.setter
    def phi_grid(self, grid):
        cv.check_type('mesh phi_grid', grid, Iterable, Real)
        cv.check_length('mesh phi_grid', grid, 2)
        cv.check_increasing('mesh phi_grid', grid)
        grid = np.asarray(grid, dtype=float)
        if np.any((grid < 0.0) | (grid > 2*pi)):
            raise ValueError("phi_grid values must be in [0, 2π].")
        self._phi_grid = grid

    @property
    def z_grid(self):
        return self._z_grid

    @z_grid.setter
    def z_grid(self, grid):
        cv.check_type('mesh z_grid', grid, Iterable, Real)
        cv.check_length('mesh z_grid', grid, 2)
        cv.check_increasing('mesh z_grid', grid)
        self._z_grid = np.asarray(grid, dtype=float)

    @property
    def _grids(self):
        return (self.r_grid, self.phi_grid, self.z_grid)

    @property
    def indices(self):
        nr, np, nz = self.dimension
        np = len(self.phi_grid) - 1
        nz = len(self.z_grid) - 1
        return ((r, p, z)
                for z in range(1, nz + 1)
                for p in range(1, np + 1)
                for r in range(1, nr + 1))

    @property
    def lower_left(self):
        return np.array((
            self.origin[0] - self.r_grid[-1],
            self.origin[1] - self.r_grid[-1],
            self.origin[2] + self.z_grid[0]
        ))

    @property
    def upper_right(self):
        return np.array((
            self.origin[0] + self.r_grid[-1],
            self.origin[1] + self.r_grid[-1],
            self.origin[2] + self.z_grid[-1]
        ))

    def __repr__(self):
        fmt = '{0: <16}{1}{2}\n'
        string = super().__repr__()
        string += fmt.format('\tDimensions', '=\t', self.n_dimension)
        string += fmt.format('\tOrigin', '=\t', self.origin)
        r_grid_str = str(self._r_grid) if self._r_grid is None else len(self._r_grid)
        string += fmt.format('\tN R pnts:', '=\t', r_grid_str)
        if self._r_grid is not None:
            string += fmt.format('\tR Min:', '=\t', self._r_grid[0])
            string += fmt.format('\tR Max:', '=\t', self._r_grid[-1])
        phi_grid_str = str(self._phi_grid) if self._phi_grid is None else len(self._phi_grid)
        string += fmt.format('\tN Phi pnts:', '=\t', phi_grid_str)
        if self._phi_grid is not None:
            string += fmt.format('\tPhi Min:', '=\t', self._phi_grid[0])
            string += fmt.format('\tPhi Max:', '=\t', self._phi_grid[-1])
        z_grid_str = str(self._z_grid) if self._z_grid is None else len(self._z_grid)
        string += fmt.format('\tN Z pnts:', '=\t', z_grid_str)
        if self._z_grid is not None:
            string += fmt.format('\tZ Min:', '=\t', self._z_grid[0])
            string += fmt.format('\tZ Max:', '=\t', self._z_grid[-1])
        return string

    def get_indices_at_coords(
            self,
            coords: Sequence[float]
        ) -> tuple[int, int, int]:
        """Finds the index of the mesh voxel at the specified x,y,z coordinates.

        .. versionadded:: 0.15.0

        Parameters
        ----------
        coords : Sequence[float]
            The x, y, z axis coordinates

        Returns
        -------
        tuple[int, int, int]
            The r, phi, z indices

        """
        r_value_from_origin = sqrt((coords[0]-self.origin[0])**2 + (coords[1]-self.origin[1])**2)

        if r_value_from_origin < self.r_grid[0] or r_value_from_origin > self.r_grid[-1]:
            raise ValueError(
                f'The specified x, y ({coords[0]}, {coords[1]}) combine to give an r value of '
                f'{r_value_from_origin} from the origin of {self.origin}.which '
                f'is outside the origin absolute r grid values {self.r_grid}.'
            )

        r_index = np.searchsorted(self.r_grid, r_value_from_origin) - 1

        z_grid_values = np.array(self.z_grid) + self.origin[2]

        if coords[2] < z_grid_values[0] or coords[2] > z_grid_values[-1]:
            raise ValueError(
                f'The specified z value ({coords[2]}) from the z origin of '
                f'{self.origin[-1]} is outside of the absolute z grid range {z_grid_values}.'
            )

        z_index = np.argmax(z_grid_values > coords[2]) - 1

        delta_x = coords[0] - self.origin[0]
        delta_y = coords[1] - self.origin[1]
        # atan2 returns values in -pi to +pi range
        phi_value = atan2(delta_y, delta_x)
        if delta_x < 0 and delta_y < 0:
            # returned phi_value anticlockwise and negative
            phi_value += 2 * pi
        if delta_x > 0 and delta_y < 0:
            # returned phi_value anticlockwise and negative
            phi_value += 2 * pi

        phi_grid_values = np.array(self.phi_grid)

        if phi_value < phi_grid_values[0] or phi_value > phi_grid_values[-1]:
            raise ValueError(
                f'The phi value ({phi_value}) resulting from the specified x, y '
                f'values is outside of the absolute  phi grid range {phi_grid_values}.'
            )
        phi_index = np.argmax(phi_grid_values > phi_value) - 1

        return (r_index, phi_index, z_index)

    @classmethod
    def from_hdf5(cls, group: h5py.Group):
        mesh_id = int(group.name.split('/')[-1].lstrip('mesh '))

        # Read and assign mesh properties
        mesh = cls(
            mesh_id=mesh_id,
            r_grid = group['r_grid'][()],
            phi_grid = group['phi_grid'][()],
            z_grid = group['z_grid'][()],
        )
        if 'origin' in group:
            mesh.origin = group['origin'][()]

        return mesh

    @classmethod
    def from_domain(
        cls,
        domain: 'openmc.Cell' | 'openmc.Region' | 'openmc.Universe' | 'openmc.Geometry',
        dimension: Sequence[int] = (10, 10, 10),
        mesh_id: int | None = None,
        phi_grid_bounds: Sequence[float] = (0.0, 2*pi),
        name: str = ''
    ):
        """Creates a regular CylindricalMesh from an existing openmc domain.

        Parameters
        ----------
        domain : openmc.Cell or openmc.Region or openmc.Universe or openmc.Geometry
            The object passed in will be used as a template for this mesh. The
            bounding box of the property of the object passed will be used to
            set the r_grid, z_grid ranges.
        dimension : Iterable of int
            The number of equally spaced mesh cells in each direction (r_grid,
            phi_grid, z_grid)
        mesh_id : int
            Unique identifier for the mesh
        phi_grid_bounds : numpy.ndarray
            Mesh bounds points along the phi-axis in radians. The default value
            is (0, 2π), i.e., the full phi range.
        name : str
            Name of the mesh

        Returns
        -------
        openmc.CylindricalMesh
            CylindricalMesh instance

        """
        cv.check_type(
            "domain",
            domain,
            (openmc.Cell, openmc.Region, openmc.Universe, openmc.Geometry),
        )

        # loaded once to avoid recalculating bounding box
        cached_bb = domain.bounding_box
        max_bounding_box_radius = max(
            [
                cached_bb[0][0],
                cached_bb[0][1],
                cached_bb[1][0],
                cached_bb[1][1],
            ]
        )
        r_grid = np.linspace(
            0,
            max_bounding_box_radius,
            num=dimension[0]+1
        )
        phi_grid = np.linspace(
            phi_grid_bounds[0],
            phi_grid_bounds[1],
            num=dimension[1]+1
        )
        z_grid = np.linspace(
            cached_bb[0][2],
            cached_bb[1][2],
            num=dimension[2]+1
        )
        origin = (cached_bb.center[0], cached_bb.center[1], z_grid[0])

        # make z-grid relative to the origin
        z_grid -= origin[2]

        mesh = cls(
            r_grid=r_grid,
            z_grid=z_grid,
            phi_grid=phi_grid,
            mesh_id=mesh_id,
            name=name,
            origin=origin
        )

        return mesh

    def to_xml_element(self):
        """Return XML representation of the mesh

        Returns
        -------
        element : lxml.etree._Element
            XML element containing mesh data

        """

        element = super().to_xml_element()
        element.set("type", "cylindrical")

        subelement = ET.SubElement(element, "r_grid")
        subelement.text = ' '.join(map(str, self.r_grid))

        subelement = ET.SubElement(element, "phi_grid")
        subelement.text = ' '.join(map(str, self.phi_grid))

        subelement = ET.SubElement(element, "z_grid")
        subelement.text = ' '.join(map(str, self.z_grid))

        subelement = ET.SubElement(element, "origin")
        subelement.text = ' '.join(map(str, self.origin))

        return element

    @classmethod
    def from_xml_element(cls, elem: ET.Element):
        """Generate a cylindrical mesh from an XML element

        Parameters
        ----------
        elem : lxml.etree._Element
            XML element

        Returns
        -------
        openmc.CylindricalMesh
            Cylindrical mesh object

        """

        mesh_id = int(get_text(elem, 'id'))
        mesh = cls(
            r_grid = [float(x) for x in get_text(elem, "r_grid").split()],
            phi_grid = [float(x) for x in get_text(elem, "phi_grid").split()],
            z_grid = [float(x) for x in get_text(elem, "z_grid").split()],
            origin = [float(x) for x in get_text(elem, "origin", default=[0., 0., 0.]).split()],
            mesh_id=mesh_id,
        )

        return mesh

    @property
    def volumes(self):
        """Return Volumes for every mesh cell

        Returns
        -------
        volumes : Iterable of float
            Volumes

        """
        self._volume_dim_check()
        V_r = np.diff(np.asarray(self.r_grid)**2 / 2)
        V_p = np.diff(self.phi_grid)
        V_z = np.diff(self.z_grid)

        return np.multiply.outer(np.outer(V_r, V_p), V_z)

    @property
    def vertices(self):
        warnings.warn('Cartesian coordinates are returned from this property as of version 0.14.0')
        return self._convert_to_cartesian(self.vertices_cylindrical, self.origin)

    @property
    def vertices_cylindrical(self):
        """Returns vertices of the mesh in cylindrical coordinates.
        """
        return super().vertices

    @property
    def centroids(self):
        warnings.warn('Cartesian coordinates are returned from this property as of version 0.14.0')
        return self._convert_to_cartesian(self.centroids_cylindrical, self.origin)

    @property
    def centroids_cylindrical(self):
        """Returns centroids of the mesh in cylindrical coordinates.
        """
        return super().centroids

    @staticmethod
    def _convert_to_cartesian(arr, origin: Sequence[float]):
        """Converts an array with r, phi, z values in the last dimension (shape (..., 3))
        to Cartesian coordinates.
        """
        x = arr[..., 0] * np.cos(arr[..., 1]) + origin[0]
        y = arr[..., 0] * np.sin(arr[..., 1]) + origin[1]
        arr[..., 0] = x
        arr[..., 1] = y
        arr[..., 2] += origin[2]
        return arr


class SphericalMesh(StructuredMesh):
    """A 3D spherical mesh

    Parameters
    ----------
    r_grid : numpy.ndarray
        1-D array of mesh boundary points along the r-axis.
        Requirement is r >= 0.
    phi_grid : numpy.ndarray
        1-D array of mesh boundary points along the phi-axis in radians.
        The default value is [0, 2π], i.e. the full phi range.
    theta_grid : numpy.ndarray
        1-D array of mesh boundary points along the theta-axis in radians.
        The default value is [0, π], i.e. the full theta range.
    origin : numpy.ndarray
        1-D array of length 3 the (x,y,z) origin of the mesh in
        cartesian coordinates
    mesh_id : int
        Unique identifier for the mesh
    name : str
        Name of the mesh

    Attributes
    ----------
    id : int
        Unique identifier for the mesh
    name : str
        Name of the mesh
    dimension : Iterable of int
        The number of mesh cells in each direction (r_grid,
        theta_grid, phi_grid).
    n_dimension : int
        Number of mesh dimensions (always 3 for a SphericalMesh).
    r_grid : numpy.ndarray
        1-D array of mesh boundary points along the r-axis.
        Requirement is r >= 0.
    theta_grid : numpy.ndarray
        1-D array of mesh boundary points along the theta-axis in radians.
        The default value is [0, π], i.e. the full theta range.
    phi_grid : numpy.ndarray
        1-D array of mesh boundary points along the phi-axis in radians.
        The default value is [0, 2π], i.e. the full phi range.
    origin : numpy.ndarray
        1-D array of length 3 the (x,y,z) origin of the mesh in
        cartesian coordinates
    indices : Iterable of tuple
        An iterable of mesh indices for each mesh element, e.g. [(1, 1, 1),
        (2, 1, 1), ...]
    lower_left : numpy.ndarray
        The lower-left corner of the structured mesh. If only two coordinate
        are given, it is assumed that the mesh is an x-y mesh.
    upper_right : numpy.ndarray
        The upper-right corner of the structured mesh. If only two coordinate
        are given, it is assumed that the mesh is an x-y mesh.
    bounding_box : openmc.BoundingBox
        Axis-aligned bounding box of the mesh as defined by the upper-right and
        lower-left coordinates.

    """

    def __init__(
        self,
        r_grid: Sequence[float],
        phi_grid: Sequence[float] = (0, 2*pi),
        theta_grid: Sequence[float] = (0, pi),
        origin: Sequence[float] = (0., 0., 0.),
        mesh_id: int | None = None,
        name: str = '',
    ):
        super().__init__(mesh_id, name)

        self.r_grid = r_grid
        self.theta_grid = theta_grid
        self.phi_grid = phi_grid
        self.origin = origin

    @property
    def dimension(self):
        return (len(self.r_grid) - 1,
                len(self.theta_grid) - 1,
                len(self.phi_grid) - 1)

    @property
    def n_dimension(self):
        return 3

    @property
    def origin(self):
        return self._origin

    @origin.setter
    def origin(self, coords):
        cv.check_type('mesh origin', coords, Iterable, Real)
        cv.check_length("mesh origin", coords, 3)
        self._origin = np.asarray(coords, dtype=float)

    @property
    def r_grid(self):
        return self._r_grid

    @r_grid.setter
    def r_grid(self, grid):
        cv.check_type('mesh r_grid', grid, Iterable, Real)
        cv.check_length('mesh r_grid', grid, 2)
        cv.check_increasing('mesh r_grid', grid)
        self._r_grid = np.asarray(grid, dtype=float)

    @property
    def theta_grid(self):
        return self._theta_grid

    @theta_grid.setter
    def theta_grid(self, grid):
        cv.check_type('mesh theta_grid', grid, Iterable, Real)
        cv.check_length('mesh theta_grid', grid, 2)
        cv.check_increasing('mesh theta_grid', grid)
        grid = np.asarray(grid, dtype=float)
        if np.any((grid < 0.0) | (grid > pi)):
            raise ValueError("theta_grid values must be in [0, π].")
        self._theta_grid = grid

    @property
    def phi_grid(self):
        return self._phi_grid

    @phi_grid.setter
    def phi_grid(self, grid):
        cv.check_type('mesh phi_grid', grid, Iterable, Real)
        cv.check_length('mesh phi_grid', grid, 2)
        cv.check_increasing('mesh phi_grid', grid)
        grid = np.asarray(grid, dtype=float)
        if np.any((grid < 0.0) | (grid > 2*pi)):
            raise ValueError("phi_grid values must be in [0, 2π].")
        self._phi_grid = grid

    @property
    def _grids(self):
        return (self.r_grid, self.theta_grid, self.phi_grid)

    @property
    def indices(self):
        nr, nt, np = self.dimension
        nt = len(self.theta_grid) - 1
        np = len(self.phi_grid) - 1
        return ((r, t, p)
                for p in range(1, np + 1)
                for t in range(1, nt + 1)
                for r in range(1, nr + 1))

    @property
    def lower_left(self):
        r = self.r_grid[-1]
        return np.array((self.origin[0] - r, self.origin[1] - r, self.origin[2] - r))

    @property
    def upper_right(self):
        r = self.r_grid[-1]
        return np.array((self.origin[0] + r, self.origin[1] + r, self.origin[2] + r))

    def __repr__(self):
        fmt = '{0: <16}{1}{2}\n'
        string = super().__repr__()
        string += fmt.format('\tDimensions', '=\t', self.n_dimension)
        string += fmt.format('\tOrigin', '=\t', self.origin)
        r_grid_str = str(self._r_grid) if self._r_grid is None else len(self._r_grid)
        string += fmt.format('\tN R pnts:', '=\t', r_grid_str)
        if self._r_grid is not None:
            string += fmt.format('\tR Min:', '=\t', self._r_grid[0])
            string += fmt.format('\tR Max:', '=\t', self._r_grid[-1])
        theta_grid_str = str(self._theta_grid) if self._theta_grid is None else len(self._theta_grid)
        string += fmt.format('\tN Theta pnts:', '=\t', theta_grid_str)
        if self._theta_grid is not None:
            string += fmt.format('\tTheta Min:', '=\t', self._theta_grid[0])
            string += fmt.format('\tTheta Max:', '=\t', self._theta_grid[-1])
        phi_grid_str = str(self._phi_grid) if self._phi_grid is None else len(self._phi_grid)
        string += fmt.format('\tN Phi pnts:', '=\t', phi_grid_str)
        if self._phi_grid is not None:
            string += fmt.format('\tPhi Min:', '=\t', self._phi_grid[0])
            string += fmt.format('\tPhi Max:', '=\t', self._phi_grid[-1])
        return string

    @classmethod
    def from_hdf5(cls, group: h5py.Group):
        mesh_id = int(group.name.split('/')[-1].lstrip('mesh '))

        # Read and assign mesh properties
        mesh = cls(
            r_grid = group['r_grid'][()],
            theta_grid = group['theta_grid'][()],
            phi_grid = group['phi_grid'][()],
            mesh_id=mesh_id,
        )
        if 'origin' in group:
            mesh.origin = group['origin'][()]

        return mesh

    def to_xml_element(self):
        """Return XML representation of the mesh

        Returns
        -------
        element : lxml.etree._Element
            XML element containing mesh data

        """

        element = super().to_xml_element()
        element.set("type", "spherical")

        subelement = ET.SubElement(element, "r_grid")
        subelement.text = ' '.join(map(str, self.r_grid))

        subelement = ET.SubElement(element, "theta_grid")
        subelement.text = ' '.join(map(str, self.theta_grid))

        subelement = ET.SubElement(element, "phi_grid")
        subelement.text = ' '.join(map(str, self.phi_grid))

        subelement = ET.SubElement(element, "origin")
        subelement.text = ' '.join(map(str, self.origin))

        return element

    @classmethod
    def from_xml_element(cls, elem: ET.Element):
        """Generate a spherical mesh from an XML element

        Parameters
        ----------
        elem : lxml.etree._Element
            XML element

        Returns
        -------
        openmc.SphericalMesh
            Spherical mesh object

        """
        mesh_id = int(get_text(elem, 'id'))
        mesh = cls(
            mesh_id=mesh_id,
            r_grid = [float(x) for x in get_text(elem, "r_grid").split()],
            theta_grid = [float(x) for x in get_text(elem, "theta_grid").split()],
            phi_grid = [float(x) for x in get_text(elem, "phi_grid").split()],
            origin = [float(x) for x in get_text(elem, "origin", default=[0., 0., 0.]).split()],
        )

        return mesh

    @property
    def volumes(self):
        """Return Volumes for every mesh cell

        Returns
        -------
        volumes : Iterable of float
            Volumes

        """
        self._volume_dim_check()
        V_r = np.diff(np.asarray(self.r_grid)**3 / 3)
        V_t = np.diff(-np.cos(self.theta_grid))
        V_p = np.diff(self.phi_grid)

        return np.multiply.outer(np.outer(V_r, V_t), V_p)

    @property
    def vertices(self):
        warnings.warn('Cartesian coordinates are returned from this property as of version 0.14.0')
        return self._convert_to_cartesian(self.vertices_spherical, self.origin)

    @property
    def vertices_spherical(self):
        """Returns vertices of the mesh in cylindrical coordinates.
        """
        return super().vertices

    @property
    def centroids(self):
        warnings.warn('Cartesian coordinates are returned from this property as of version 0.14.0')
        return self._convert_to_cartesian(self.centroids_spherical, self.origin)

    @property
    def centroids_spherical(self):
        """Returns centroids of the mesh in cylindrical coordinates.
        """
        return super().centroids


    @staticmethod
    def _convert_to_cartesian(arr, origin: Sequence[float]):
        """Converts an array with r, theta, phi values in the last dimension (shape (..., 3))
        to Cartesian coordinates.
        """
        r_xy = arr[..., 0] * np.sin(arr[..., 1])
        x = r_xy * np.cos(arr[..., 2])
        y = r_xy * np.sin(arr[..., 2])
        z = arr[..., 0] * np.cos(arr[..., 1])
        arr[..., 0] = x + origin[0]
        arr[..., 1] = y + origin[1]
        arr[..., 2] = z + origin[2]
        return arr


def require_statepoint_data(func):
    @wraps(func)
    def wrapper(self: UnstructuredMesh, *args, **kwargs):
        if not self._has_statepoint_data:
            raise AttributeError(f'The "{func.__name__}" property requires '
                                 'information about this mesh to be loaded '
                                 'from a statepoint file.')
        return func(self, *args, **kwargs)
    return wrapper


class UnstructuredMesh(MeshBase):
    """A 3D unstructured mesh

    .. versionadded:: 0.12

    .. versionchanged:: 0.12.2
        Support for libMesh unstructured meshes was added.

    Parameters
    ----------
    filename : path-like
        Location of the unstructured mesh file
    library : {'moab', 'libmesh'}
        Mesh library used for the unstructured mesh tally
    mesh_id : int
        Unique identifier for the mesh
    name : str
        Name of the mesh
    length_multiplier: float
        Constant multiplier to apply to mesh coordinates
    options : str, optional
        Special options that control spatial search data structures used. This
        is currently only used to set `parameters
        <https://tinyurl.com/kdtree-params>`_ for MOAB's AdaptiveKDTree. If
        None, OpenMC internally uses a default of "MAX_DEPTH=20;PLANE_SET=2;".

    Attributes
    ----------
    id : int
        Unique identifier for the mesh
    name : str
        Name of the mesh
    filename : str
        Name of the file containing the unstructured mesh
    length_multiplier: float
        Multiplicative factor to apply to mesh coordinates
    library : {'moab', 'libmesh'}
        Mesh library used for the unstructured mesh tally
    options : str
        Special options that control spatial search data structures used. This
        is currently only used to set `parameters
        <https://tinyurl.com/kdtree-params>`_ for MOAB's AdaptiveKDTree. If
        None, OpenMC internally uses a default of "MAX_DEPTH=20;PLANE_SET=2;".
    output : bool
        Indicates whether or not automatic tally output should be generated for
        this mesh
    volumes : Iterable of float
        Volumes of the unstructured mesh elements
    centroids : numpy.ndarray
        Centroids of the mesh elements with array shape (n_elements, 3)

    vertices : numpy.ndarray
        Coordinates of the mesh vertices with array shape (n_elements, 3)

        .. versionadded:: 0.13.1
    connectivity : numpy.ndarray
        Connectivity of the elements with array shape (n_elements, 8)

        .. versionadded:: 0.13.1
    element_types : Iterable of integers
        Mesh element types

        .. versionadded:: 0.13.1
    total_volume : float
        Volume of the unstructured mesh in total
    bounding_box : openmc.BoundingBox
        Axis-aligned bounding box of the mesh as defined by the upper-right and
        lower-left coordinates.

    """

    _UNSUPPORTED_ELEM = -1
    _LINEAR_TET = 0
    _LINEAR_HEX = 1

    def __init__(self, filename: PathLike, library: str, mesh_id: int | None = None,
                 name: str = '', length_multiplier: float = 1.0,
                 options: str | None = None):
        super().__init__(mesh_id, name)
        self.filename = filename
        self._volumes = None
        self._n_elements = None
        self._conectivity = None
        self._vertices = None
        self.library = library
        self._output = False
        self.length_multiplier = length_multiplier
        self.options = options
        self._has_statepoint_data = False

    @property
    def filename(self):
        return self._filename

    @filename.setter
    def filename(self, filename):
        cv.check_type('Unstructured Mesh filename', filename, PathLike)
        self._filename = input_path(filename)

    @property
    def library(self):
        return self._library

    @library.setter
    def library(self, lib: str):
        cv.check_value('Unstructured mesh library', lib, ('moab', 'libmesh'))
        self._library = lib

    @property
    def options(self) -> str | None:
        return self._options

    @options.setter
    def options(self, options: str | None):
        cv.check_type('options', options, (str, type(None)))
        self._options = options

    @property
    @require_statepoint_data
    def size(self):
        return self._size

    @size.setter
    def size(self, size: int):
        cv.check_type("Unstructured mesh size", size, Integral)
        self._size = size

    @property
    def output(self):
        return self._output

    @output.setter
    def output(self, val: bool):
        cv.check_type("Unstructured mesh output value", val, bool)
        self._output = val

    @property
    @require_statepoint_data
    def volumes(self):
        """Return Volumes for every mesh cell if
        populated by a StatePoint file

        Returns
        -------
        volumes : numpy.ndarray
            Volumes

        """
        return self._volumes

    @volumes.setter
    def volumes(self, volumes: Iterable[Real]):
        cv.check_type("Unstructured mesh volumes", volumes, Iterable, Real)
        self._volumes = volumes

    @property
    @require_statepoint_data
    def total_volume(self):
        return np.sum(self.volumes)

    @property
    @require_statepoint_data
    def vertices(self):
        return self._vertices

    @property
    @require_statepoint_data
    def connectivity(self):
        return self._connectivity

    @property
    @require_statepoint_data
    def element_types(self):
        return self._element_types

    @property
    @require_statepoint_data
    def centroids(self):
        return np.array([self.centroid(i) for i in range(self.n_elements)])

    @property
    @require_statepoint_data
    def n_elements(self):
        if self._n_elements is None:
            raise RuntimeError("No information about this mesh has "
                               "been loaded from a statepoint file.")
        return self._n_elements

    @n_elements.setter
    def n_elements(self, val: int):
        cv.check_type('Number of elements', val, Integral)
        self._n_elements = val

    @property
    def length_multiplier(self):
        return self._length_multiplier

    @length_multiplier.setter
    def length_multiplier(self, length_multiplier):
        cv.check_type("Unstructured mesh length multiplier",
                      length_multiplier,
                      Real)
        self._length_multiplier = length_multiplier

    @property
    def dimension(self):
        return (self.n_elements,)

    @property
    def n_dimension(self):
        return 3

    @property
    @require_statepoint_data
    def indices(self):
        return [(i,) for i in range(self.n_elements)]

    @property
    def has_statepoint_data(self) -> bool:
        return self._has_statepoint_data

    def __repr__(self):
        string = super().__repr__()
        string += '{: <16}=\t{}\n'.format('\tFilename', self.filename)
        string += '{: <16}=\t{}\n'.format('\tMesh Library', self.library)
        if self.length_multiplier != 1.0:
            string += '{: <16}=\t{}\n'.format('\tLength multiplier',
                                              self.length_multiplier)
        if self.options is not None:
            string += '{: <16}=\t{}\n'.format('\tOptions', self.options)
        return string

    @property
    @require_statepoint_data
    def lower_left(self):
        return self.vertices.min(axis=0)

    @property
    @require_statepoint_data
    def upper_right(self):
        return self.vertices.max(axis=0)

    @require_statepoint_data
    def centroid(self, bin: int):
        """Return the vertex averaged centroid of an element

        Parameters
        ----------
        bin : int
            Bin ID for the returned centroid

        Returns
        -------
        numpy.ndarray
            x, y, z values of the element centroid

        """
        conn = self.connectivity[bin]
        # remove invalid connectivity values
        conn = conn[conn >= 0]
        coords = self.vertices[conn]
        return coords.mean(axis=0)

    def write_vtk_mesh(self, **kwargs):
        """Map data to unstructured VTK mesh elements.

        .. deprecated:: 0.13
          Use :func:`UnstructuredMesh.write_data_to_vtk` instead.

        Parameters
        ----------
        filename : str or pathlib.Path
            Name of the VTK file to write.
        datasets : dict
            Dictionary whose keys are the data labels
            and values are the data sets.
        volume_normalization : bool
            Whether or not to normalize the data by the
            volume of the mesh elements
        """
        warnings.warn(
            "The 'UnstructuredMesh.write_vtk_mesh' method has been renamed "
            "to 'write_data_to_vtk' and will be removed in a future version "
            " of OpenMC.", FutureWarning
        )
        self.write_data_to_vtk(**kwargs)

    def write_data_to_vtk(
            self,
            filename: PathLike | None  = None,
            datasets: dict | None = None,
            volume_normalization: bool = True
    ):
        """Map data to unstructured VTK mesh elements.

        Parameters
        ----------
        filename : str or pathlib.Path
            Name of the VTK file to write
        datasets : dict
            Dictionary whose keys are the data labels
            and values are numpy appropriately sized arrays
            of the data
        volume_normalization : bool
            Whether or not to normalize the data by the
            volume of the mesh elements
        """
        import vtk
        from vtk.util import numpy_support as nps

        if self.connectivity is None or self.vertices is None:
            raise RuntimeError('This mesh has not been '
                               'loaded from a statepoint file.')

        if filename is None:
            filename = f'mesh_{self.id}.vtk'

        writer = vtk.vtkUnstructuredGridWriter()

        writer.SetFileName(str(filename))

        grid = vtk.vtkUnstructuredGrid()

        vtk_pnts = vtk.vtkPoints()
        vtk_pnts.SetData(nps.numpy_to_vtk(self.vertices))
        grid.SetPoints(vtk_pnts)

        n_skipped = 0
        for elem_type, conn in zip(self.element_types, self.connectivity):
            if elem_type == self._LINEAR_TET:
                elem = vtk.vtkTetra()
            elif elem_type == self._LINEAR_HEX:
                elem = vtk.vtkHexahedron()
            elif elem_type == self._UNSUPPORTED_ELEM:
                n_skipped += 1
            else:
                raise RuntimeError(f'Invalid element type {elem_type} found')
            for i, c in enumerate(conn):
                if c == -1:
                    break
                elem.GetPointIds().SetId(i, c)

            grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())

        if n_skipped > 0:
            warnings.warn(f'{n_skipped} elements were not written because '
                          'they are not of type linear tet/hex')

        # check that datasets are the correct size
        datasets_out = []
        if datasets is not None:
            for name, data in datasets.items():
                if data.shape != self.dimension:
                    raise ValueError(f'Cannot apply dataset "{name}" with '
                                     f'shape {data.shape} to mesh {self.id} '
                                     f'with dimensions {self.dimension}')

            if volume_normalization:
                for name, data in datasets.items():
                    if np.issubdtype(data.dtype, np.integer):
                        warnings.warn(f'Integer data set "{name}" will '
                                      'not be volume-normalized.')
                        continue
                    data /= self.volumes

            # add data to the mesh
            for name, data in datasets.items():
                datasets_out.append(data)
                arr = vtk.vtkDoubleArray()
                arr.SetName(name)
                arr.SetNumberOfTuples(data.size)

                for i in range(data.size):
                    arr.SetTuple1(i, data.flat[i])
                grid.GetCellData().AddArray(arr)

        writer.SetInputData(grid)

        writer.Write()

    @classmethod
    def from_hdf5(cls, group: h5py.Group):
        mesh_id = int(group.name.split('/')[-1].lstrip('mesh '))
        filename = group['filename'][()].decode()
        library = group['library'][()].decode()
        if 'options' in group.attrs:
            options = group.attrs['options'].decode()
        else:
            options = None

        mesh = cls(filename=filename, library=library, mesh_id=mesh_id, options=options)
        mesh._has_statepoint_data = True
        vol_data = group['volumes'][()]
        mesh.volumes = np.reshape(vol_data, (vol_data.shape[0],))
        mesh.n_elements = mesh.volumes.size

        vertices = group['vertices'][()]
        mesh._vertices = vertices.reshape((-1, 3))
        connectivity = group['connectivity'][()]
        mesh._connectivity = connectivity.reshape((-1, 8))
        mesh._element_types = group['element_types'][()]

        if 'length_multiplier' in group:
            mesh.length_multiplier = group['length_multiplier'][()]

        return mesh

    def to_xml_element(self):
        """Return XML representation of the mesh

        Returns
        -------
        element : lxml.etree._Element
            XML element containing mesh data

        """

        element = super().to_xml_element()
        element.set("type", "unstructured")

        element.set("library", self._library)
        if self.options is not None:
            element.set('options', self.options)
        subelement = ET.SubElement(element, "filename")
        subelement.text = str(self.filename)

        if self._length_multiplier != 1.0:
            element.set("length_multiplier", str(self.length_multiplier))

        return element

    @classmethod
    def from_xml_element(cls, elem: ET.Element):
        """Generate unstructured mesh object from XML element

        Parameters
        ----------
        elem : lxml.etree._Element
            XML element

        Returns
        -------
        openmc.UnstructuredMesh
            UnstructuredMesh generated from an XML element
        """
        mesh_id = int(get_text(elem, 'id'))
        filename = get_text(elem, 'filename')
        library = get_text(elem, 'library')
        length_multiplier = float(get_text(elem, 'length_multiplier', 1.0))
        options = elem.get('options')

        return cls(filename, library, mesh_id, '', length_multiplier, options)


def _read_meshes(elem):
    """Generate dictionary of meshes from a given XML node

    Parameters
    ----------
    elem : lxml.etree._Element
        XML element

    Returns
    -------
    dict
        A dictionary with mesh IDs as keys and openmc.MeshBase
        instanaces as values
    """
    out = {}
    for mesh_elem in elem.findall('mesh'):
        mesh = MeshBase.from_xml_element(mesh_elem)
        out[mesh.id] = mesh

    return out


# hexahedron element connectivity
# lower-k connectivity offsets
_HEX_VERTEX_CONN = ((0, 0, 0),
                    (1, 0, 0),
                    (1, 1, 0),
                    (0, 1, 0))
# upper-k connectivity offsets
_HEX_VERTEX_CONN += ((0, 0, 1),
                     (1, 0, 1),
                     (1, 1, 1),
                     (0, 1, 1))

_N_HEX_VERTICES = 8

# lower-k connectivity offsets
_HEX_MIDPOINT_CONN = ((0, (0, 0, 0)),
                      (1, (1, 0, 0)),
                      (0, (0, 1, 0)),
                      (1, (0, 0, 0)))
# upper-k connectivity offsets
_HEX_MIDPOINT_CONN += ((0, (0, 0, 1)),
                       (1, (1, 0, 1)),
                       (0, (0, 1, 1)),
                       (1, (0, 0, 1)))
# mid-plane k connectivity
_HEX_MIDPOINT_CONN += ((2, (0, 0, 0)),
                       (2, (1, 0, 0)),
                       (2, (1, 1, 0)),
                       (2, (0, 1, 0)))
