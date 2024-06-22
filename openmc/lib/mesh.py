from collections.abc import Mapping
from ctypes import (c_int, c_int32, c_char_p, c_double, POINTER, Structure,
                    create_string_buffer, c_uint64, c_size_t)
from random import getrandbits
from typing import Optional, List, Tuple, Sequence
from weakref import WeakValueDictionary

import numpy as np
from numpy.ctypeslib import as_array

from ..exceptions import AllocationError, InvalidIDError
from . import _dll
from .core import _FortranObjectWithID
from .error import _error_handler
from .material import Material
from .plot import _Position

__all__ = [
    'Mesh', 'RegularMesh', 'RectilinearMesh', 'CylindricalMesh',
    'SphericalMesh', 'UnstructuredMesh', 'meshes'
]


class _MaterialVolume(Structure):
    _fields_ = [
        ("material", c_int32),
        ("volume", c_double)
    ]


# Mesh functions
_dll.openmc_extend_meshes.argtypes = [c_int32, c_char_p, POINTER(c_int32),
                                      POINTER(c_int32)]
_dll.openmc_extend_meshes.restype = c_int
_dll.openmc_extend_meshes.errcheck = _error_handler
_dll.openmc_mesh_get_id.argtypes = [c_int32, POINTER(c_int32)]
_dll.openmc_mesh_get_id.restype = c_int
_dll.openmc_mesh_get_id.errcheck = _error_handler
_dll.openmc_mesh_set_id.argtypes = [c_int32, c_int32]
_dll.openmc_mesh_set_id.restype = c_int
_dll.openmc_mesh_set_id.errcheck = _error_handler
_dll.openmc_mesh_get_n_elements.argtypes = [c_int32, POINTER(c_size_t)]
_dll.openmc_mesh_get_n_elements.restype = c_int
_dll.openmc_mesh_get_n_elements.errcheck = _error_handler
_dll.openmc_mesh_get_volumes.argtypes = [c_int32, POINTER(c_double)]
_dll.openmc_mesh_get_volumes.restype = c_int
_dll.openmc_mesh_get_volumes.errcheck = _error_handler
_dll.openmc_mesh_material_volumes.argtypes = [
    c_int32, c_int, c_int, c_int, POINTER(_MaterialVolume),
    POINTER(c_int), POINTER(c_uint64)]
_dll.openmc_mesh_material_volumes.restype = c_int
_dll.openmc_mesh_material_volumes.errcheck = _error_handler
_dll.openmc_mesh_get_plot_bins.argtypes = [
    c_int32, _Position, _Position, c_int, POINTER(c_int), POINTER(c_int32)
]
_dll.openmc_mesh_get_plot_bins.restype = c_int
_dll.openmc_mesh_get_plot_bins.errcheck = _error_handler
_dll.openmc_get_mesh_index.argtypes = [c_int32, POINTER(c_int32)]
_dll.openmc_get_mesh_index.restype = c_int
_dll.openmc_get_mesh_index.errcheck = _error_handler
_dll.n_meshes.argtypes = []
_dll.n_meshes.restype = c_int
_dll.openmc_rectilinear_mesh_get_grid.argtypes = [c_int32,
    POINTER(POINTER(c_double)), POINTER(c_int), POINTER(POINTER(c_double)),
    POINTER(c_int), POINTER(POINTER(c_double)), POINTER(c_int)]
_dll.openmc_rectilinear_mesh_get_grid.restype = c_int
_dll.openmc_rectilinear_mesh_get_grid.errcheck = _error_handler
_dll.openmc_rectilinear_mesh_set_grid.argtypes = [c_int32, POINTER(c_double),
    c_int, POINTER(c_double), c_int, POINTER(c_double), c_int]
_dll.openmc_rectilinear_mesh_set_grid.restype = c_int
_dll.openmc_rectilinear_mesh_set_grid.errcheck = _error_handler
_dll.openmc_regular_mesh_get_dimension.argtypes = [c_int32,
    POINTER(POINTER(c_int)), POINTER(c_int)]
_dll.openmc_regular_mesh_get_dimension.restype = c_int
_dll.openmc_regular_mesh_get_dimension.errcheck = _error_handler
_dll.openmc_regular_mesh_get_params.argtypes = [
    c_int32, POINTER(POINTER(c_double)), POINTER(POINTER(c_double)),
    POINTER(POINTER(c_double)), POINTER(c_int)]
_dll.openmc_regular_mesh_get_params.restype = c_int
_dll.openmc_regular_mesh_get_params.errcheck = _error_handler
_dll.openmc_regular_mesh_set_dimension.argtypes = [c_int32, c_int,
                                                   POINTER(c_int)]
_dll.openmc_regular_mesh_set_dimension.restype = c_int
_dll.openmc_regular_mesh_set_dimension.errcheck = _error_handler
_dll.openmc_regular_mesh_set_params.argtypes = [
    c_int32, c_int, POINTER(c_double), POINTER(c_double), POINTER(c_double)]
_dll.openmc_regular_mesh_set_params.restype = c_int
_dll.openmc_regular_mesh_set_params.errcheck = _error_handler

_dll.openmc_cylindrical_mesh_get_grid.argtypes = [c_int32,
    POINTER(POINTER(c_double)), POINTER(c_int), POINTER(POINTER(c_double)),
    POINTER(c_int), POINTER(POINTER(c_double)), POINTER(c_int)]
_dll.openmc_cylindrical_mesh_get_grid.restype = c_int
_dll.openmc_cylindrical_mesh_get_grid.errcheck = _error_handler
_dll.openmc_cylindrical_mesh_set_grid.argtypes = [c_int32, POINTER(c_double),
    c_int, POINTER(c_double), c_int, POINTER(c_double), c_int]
_dll.openmc_cylindrical_mesh_set_grid.restype = c_int
_dll.openmc_cylindrical_mesh_set_grid.errcheck = _error_handler

_dll.openmc_spherical_mesh_get_grid.argtypes = [c_int32,
    POINTER(POINTER(c_double)), POINTER(c_int), POINTER(POINTER(c_double)),
    POINTER(c_int), POINTER(POINTER(c_double)), POINTER(c_int)]
_dll.openmc_spherical_mesh_get_grid.restype = c_int
_dll.openmc_spherical_mesh_get_grid.errcheck = _error_handler
_dll.openmc_spherical_mesh_set_grid.argtypes = [c_int32, POINTER(c_double),
    c_int, POINTER(c_double), c_int, POINTER(c_double), c_int]
_dll.openmc_spherical_mesh_set_grid.restype = c_int
_dll.openmc_spherical_mesh_set_grid.errcheck = _error_handler


class Mesh(_FortranObjectWithID):
    """Base class to represent mesh objects

    """
    __instances = WeakValueDictionary()

    def __new__(cls, uid=None, new=True, index=None):
        mapping = meshes
        if index is None:
            if new:
                # Determine ID to assign
                if uid is None:
                    uid = max(mapping, default=0) + 1
                else:
                    if uid in mapping:
                        raise AllocationError('A mesh with ID={} has already '
                                              'been allocated.'.format(uid))

                # Set the mesh type -- note that mesh type attribute only
                # exists on subclasses!
                index = c_int32()
                _dll.openmc_extend_meshes(1, cls.mesh_type.encode(), index,
                                          None)
                index = index.value
            else:
                index = mapping[uid]._index

        if index not in cls.__instances:
            instance = super().__new__(cls)
            instance._index = index
            if uid is not None:
                instance.id = uid
            cls.__instances[index] = instance

        return cls.__instances[index]

    @property
    def id(self):
        mesh_id = c_int32()
        _dll.openmc_mesh_get_id(self._index, mesh_id)
        return mesh_id.value

    @id.setter
    def id(self, mesh_id):
        _dll.openmc_mesh_set_id(self._index, mesh_id)

    @property
    def n_elements(self) -> int:
        n = c_size_t()
        _dll.openmc_mesh_get_n_elements(self._index, n)
        return n.value

    @property
    def volumes(self) -> np.ndarray:
        volumes = np.empty((self.n_elements,))
        _dll.openmc_mesh_get_volumes(
            self._index, volumes.ctypes.data_as(POINTER(c_double)))
        return volumes

    def material_volumes(
            self,
            n_samples: int = 10_000,
            prn_seed: Optional[int] = None
    ) -> List[List[Tuple[Material, float]]]:
        """Determine volume of materials in each mesh element

        .. versionadded:: 0.15.0

        Parameters
        ----------
        n_samples : int
            Number of samples in each mesh element
        prn_seed : int
            Pseudorandom number generator (PRNG) seed; if None, one will be
            generated randomly.

        Returns
        -------
        List of tuple of (material, volume) for each mesh element. Void volume
        is represented by having a value of None in the first element of a
        tuple.

        """
        if n_samples <= 0:
            raise ValueError("Number of samples must be positive")
        if prn_seed is None:
            prn_seed = getrandbits(63)
        prn_seed = c_uint64(prn_seed)

        # Preallocate space for MaterialVolume results
        size = 16
        result = (_MaterialVolume * size)()

        hits = c_int()  # Number of materials hit in a given element
        volumes = []
        for i_element in range(self.n_elements):
            while True:
                try:
                    _dll.openmc_mesh_material_volumes(
                        self._index, n_samples, i_element, size, result, hits, prn_seed)
                except AllocationError:
                    # Increase size of result array and try again
                    size *= 2
                    result = (_MaterialVolume * size)()
                else:
                    # If no error, break out of loop
                    break

            volumes.append([
                (Material(index=r.material), r.volume)
                for r in result[:hits.value]
            ])
        return volumes

    def get_plot_bins(
            self,
            origin: Sequence[float],
            width: Sequence[float],
            basis: str,
            pixels: Sequence[int]
    ) -> np.ndarray:
        """Get mesh bin indices for a rasterized plot.

        .. versionadded:: 0.15.0

        Parameters
        ----------
        origin : iterable of float
            Origin of the plotting view. Should have length 3.
        width : iterable of float
            Width of the plotting view. Should have length 2.
        basis : {'xy', 'xz', 'yz'}
            Plotting basis.
        pixels : iterable of int
            Number of pixels in each direction. Should have length 2.

        Returns
        -------
        2D numpy array with mesh bin indices corresponding to each pixel within
        the plotting view.

        """
        origin = _Position(*origin)
        width = _Position(*width)
        basis = {'xy': 1, 'xz': 2, 'yz': 3}[basis]
        pixel_array = (c_int*2)(*pixels)
        img_data = np.zeros((pixels[1], pixels[0]), dtype=np.dtype('int32'))

        _dll.openmc_mesh_get_plot_bins(
            self._index, origin, width, basis, pixel_array,
            img_data.ctypes.data_as(POINTER(c_int32))
        )
        return img_data


class RegularMesh(Mesh):
    """RegularMesh stored internally.

    This class exposes a mesh that is stored internally in the OpenMC
    library. To obtain a view of a mesh with a given ID, use the
    :data:`openmc.lib.meshes` mapping.

    Parameters
    ----------
    index : int
         Index in the `meshes` array.

    Attributes
    ----------
    id : int
        ID of the mesh
    dimension : iterable of int
        The number of mesh cells in each direction.
    lower_left : numpy.ndarray
        The lower-left corner of the structured mesh. If only two coordinate are
        given, it is assumed that the mesh is an x-y mesh.
    upper_right : numpy.ndarray
        The upper-right corner of the structrued mesh. If only two coordinate
        are given, it is assumed that the mesh is an x-y mesh.
    width : numpy.ndarray
        The width of mesh cells in each direction.
    n_elements : int
        Total number of mesh elements.
    volumes : numpy.ndarray
        Volume of each mesh element in [cm^3]

    """
    mesh_type = 'regular'

    def __init__(self, uid=None, new=True, index=None):
        super().__init__(uid, new, index)

    @property
    def dimension(self):
        dims = POINTER(c_int)()
        n = c_int()
        _dll.openmc_regular_mesh_get_dimension(self._index, dims, n)
        return tuple(as_array(dims, (n.value,)))

    @dimension.setter
    def dimension(self, dimension):
        n = len(dimension)
        dimension = (c_int*n)(*dimension)
        _dll.openmc_regular_mesh_set_dimension(self._index, n, dimension)


    @property
    def lower_left(self):
        return self._get_parameters()[0]

    @property
    def upper_right(self):
        return self._get_parameters()[1]

    @property
    def width(self):
        return self._get_parameters()[2]

    def _get_parameters(self):
        ll = POINTER(c_double)()
        ur = POINTER(c_double)()
        w = POINTER(c_double)()
        n = c_int()
        _dll.openmc_regular_mesh_get_params(self._index, ll, ur, w, n)
        return (
            as_array(ll, (n.value,)),
            as_array(ur, (n.value,)),
            as_array(w, (n.value,))
        )

    def set_parameters(self, lower_left=None, upper_right=None, width=None):
        if lower_left is not None:
            n = len(lower_left)
            lower_left = (c_double*n)(*lower_left)
        if upper_right is not None:
            n = len(upper_right)
            upper_right = (c_double*n)(*upper_right)
        if width is not None:
            n = len(width)
            width = (c_double*n)(*width)
        _dll.openmc_regular_mesh_set_params(self._index, n, lower_left, upper_right, width)


class RectilinearMesh(Mesh):
    """RectilinearMesh stored internally.

    This class exposes a mesh that is stored internally in the OpenMC
    library. To obtain a view of a mesh with a given ID, use the
    :data:`openmc.lib.meshes` mapping.

    Parameters
    ----------
    index : int
         Index in the `meshes` array.

    Attributes
    ----------
    id : int
        ID of the mesh
    dimension : iterable of int
        The number of mesh cells in each direction.
    lower_left : numpy.ndarray
        The lower-left corner of the structured mesh.
    upper_right : numpy.ndarray
        The upper-right corner of the structrued mesh.
    width : numpy.ndarray
        The width of mesh cells in each direction.
    n_elements : int
        Total number of mesh elements.
    volumes : numpy.ndarray
        Volume of each mesh element in [cm^3]

    """
    mesh_type = 'rectilinear'

    def __init__(self, uid=None, new=True, index=None):
        super().__init__(uid, new, index)

    @property
    def lower_left(self):
        return self._get_parameters()[0]

    @property
    def upper_right(self):
        return self._get_parameters()[1]

    @property
    def dimension(self):
        return self._get_parameters()[2]

    @property
    def width(self):
        return self._get_parameters()[3]

    def _get_parameters(self):
        gx = POINTER(c_double)()
        nx = c_int()
        gy = POINTER(c_double)()
        ny = c_int()
        gz = POINTER(c_double)()
        nz = c_int()
        # Call C API to get grid parameters
        _dll.openmc_rectilinear_mesh_get_grid(self._index, gx, nx, gy, ny, gz,
                                              nz)

        # Convert grid parameters to Numpy arrays
        grid_x = as_array(gx, (nx.value,))
        grid_y = as_array(gy, (ny.value,))
        grid_z = as_array(gz, (nz.value,))

        # Calculate lower_left, upper_right, width, and dimension from grid
        lower_left = np.array((grid_x[0], grid_y[0], grid_z[0]))
        upper_right = np.array((grid_x[-1], grid_y[-1], grid_z[-1]))
        dimension = np.array((nx.value - 1, ny.value - 1, nz.value - 1))
        width = np.zeros(list(dimension) + [3])

        for i, diff_x in enumerate(np.diff(grid_x)):
            for j, diff_y in enumerate(np.diff(grid_y)):
                for k, diff_z in enumerate(np.diff(grid_z)):
                    width[i, j, k, :] = diff_x, diff_y, diff_z

        return (lower_left, upper_right, dimension, width)

    def set_grid(self, x_grid, y_grid, z_grid):
        """Set grid values

        Parameters
        ----------
        x_grid : iterable of float
            Mesh boundary points along the x-axis
        y_grid : iterable of float
            Mesh boundary points along the y-axis
        z_grid : iterable of float
            Mesh boundary points along the z-axis

        """
        nx = len(x_grid)
        x_grid = (c_double*nx)(*x_grid)
        ny = len(y_grid)
        y_grid = (c_double*ny)(*y_grid)
        nz = len(z_grid)
        z_grid = (c_double*nz)(*z_grid)
        _dll.openmc_rectilinear_mesh_set_grid(self._index, x_grid, nx, y_grid,
                                              ny, z_grid, nz)


class CylindricalMesh(Mesh):
    """CylindricalMesh stored internally.

    This class exposes a mesh that is stored internally in the OpenMC
    library. To obtain a view of a mesh with a given ID, use the
    :data:`openmc.lib.meshes` mapping.

    Parameters
    ----------
    index : int
         Index in the `meshes` array.

    Attributes
    ----------
    id : int
        ID of the mesh
    dimension : iterable of int
        The number of mesh cells in each direction.
    lower_left : numpy.ndarray
        The lower-left corner of the structured mesh.
    upper_right : numpy.ndarray
        The upper-right corner of the structrued mesh.
    width : numpy.ndarray
        The width of mesh cells in each direction.
    n_elements : int
        Total number of mesh elements.
    volumes : numpy.ndarray
        Volume of each mesh element in [cm^3]

    """
    mesh_type = 'cylindrical'

    def __init__(self, uid=None, new=True, index=None):
        super().__init__(uid, new, index)

    @property
    def lower_left(self):
        return self._get_parameters()[0]

    @property
    def upper_right(self):
        return self._get_parameters()[1]

    @property
    def dimension(self):
        return self._get_parameters()[2]

    @property
    def width(self):
        return self._get_parameters()[3]

    def _get_parameters(self):
        gx = POINTER(c_double)()
        nx = c_int()
        gy = POINTER(c_double)()
        ny = c_int()
        gz = POINTER(c_double)()
        nz = c_int()
        # Call C API to get grid parameters
        _dll.openmc_cylindrical_mesh_get_grid(self._index, gx, nx, gy, ny, gz,
                                              nz)

        # Convert grid parameters to Numpy arrays
        grid_x = as_array(gx, (nx.value,))
        grid_y = as_array(gy, (ny.value,))
        grid_z = as_array(gz, (nz.value,))

        # Calculate lower_left, upper_right, width, and dimension from grid
        lower_left = np.array((grid_x[0], grid_y[0], grid_z[0]))
        upper_right = np.array((grid_x[-1], grid_y[-1], grid_z[-1]))
        dimension = np.array((nx.value - 1, ny.value - 1, nz.value - 1))
        width = np.zeros(list(dimension) + [3])

        for i, diff_x in enumerate(np.diff(grid_x)):
            for j, diff_y in enumerate(np.diff(grid_y)):
                for k, diff_z in enumerate(np.diff(grid_z)):
                    width[i, j, k, :] = diff_x, diff_y, diff_z

        return (lower_left, upper_right, dimension, width)

    def set_grid(self, r_grid, phi_grid, z_grid):
        """Set grid values

        Parameters
        ----------
        r_grid : iterable of float
            Mesh boundary points along the r-axis
        phi_grid : Iterable of float
            Mesh boundary points along the phi-axis
        z_grid : Iterable of float
            Mesh boundary points along the z-axis

        """
        nr = len(r_grid)
        r_grid = (c_double*nr)(*r_grid)
        nphi = len(phi_grid)
        phi_grid = (c_double*nphi)(*phi_grid)
        nz = len(z_grid)
        z_grid = (c_double*nz)(*z_grid)
        _dll.openmc_cylindrical_mesh_set_grid(self._index, r_grid, nr, phi_grid,
                                              nphi, z_grid, nz)


class SphericalMesh(Mesh):
    """SphericalMesh stored internally.

    This class exposes a mesh that is stored internally in the OpenMC
    library. To obtain a view of a mesh with a given ID, use the
    :data:`openmc.lib.meshes` mapping.

    Parameters
    ----------
    index : int
         Index in the `meshes` array.

    Attributes
    ----------
    id : int
        ID of the mesh
    dimension : iterable of int
        The number of mesh cells in each direction.
    lower_left : numpy.ndarray
        The lower-left corner of the structured mesh.
    upper_right : numpy.ndarray
        The upper-right corner of the structrued mesh.
    width : numpy.ndarray
        The width of mesh cells in each direction.
    n_elements : int
        Total number of mesh elements.
    volumes : numpy.ndarray
        Volume of each mesh element in [cm^3]

    """
    mesh_type = 'spherical'

    def __init__(self, uid=None, new=True, index=None):
        super().__init__(uid, new, index)

    @property
    def lower_left(self):
        return self._get_parameters()[0]

    @property
    def upper_right(self):
        return self._get_parameters()[1]

    @property
    def dimension(self):
        return self._get_parameters()[2]

    @property
    def width(self):
        return self._get_parameters()[3]

    def _get_parameters(self):
        gx = POINTER(c_double)()
        nx = c_int()
        gy = POINTER(c_double)()
        ny = c_int()
        gz = POINTER(c_double)()
        nz = c_int()
        # Call C API to get grid parameters
        _dll.openmc_spherical_mesh_get_grid(self._index, gx, nx, gy, ny, gz,
                                              nz)

        # Convert grid parameters to Numpy arrays
        grid_x = as_array(gx, (nx.value,))
        grid_y = as_array(gy, (ny.value,))
        grid_z = as_array(gz, (nz.value,))

        # Calculate lower_left, upper_right, width, and dimension from grid
        lower_left = np.array((grid_x[0], grid_y[0], grid_z[0]))
        upper_right = np.array((grid_x[-1], grid_y[-1], grid_z[-1]))
        dimension = np.array((nx.value - 1, ny.value - 1, nz.value - 1))
        width = np.zeros(list(dimension) + [3])

        for i, diff_x in enumerate(np.diff(grid_x)):
            for j, diff_y in enumerate(np.diff(grid_y)):
                for k, diff_z in enumerate(np.diff(grid_z)):
                    width[i, j, k, :] = diff_x, diff_y, diff_z

        return (lower_left, upper_right, dimension, width)

    def set_grid(self, r_grid, theta_grid, phi_grid):
        """Set grid values

        Parameters
        ----------
        r_grid : iterable of float
            Mesh boundary points along the r-axis
        theta_grid : Iterable of float
            Mesh boundary points along the theta-axis
        phi_grid : Iterable of float
            Mesh boundary points along the phi-axis

        """
        nr = len(r_grid)
        r_grid = (c_double*nr)(*r_grid)
        ntheta = len(theta_grid)
        theta_grid = (c_double*ntheta)(*theta_grid)
        nphi = len(phi_grid)
        phi_grid = (c_double*nphi)(*phi_grid)
        _dll.openmc_spherical_mesh_set_grid(self._index, r_grid, nr, theta_grid,
                                              ntheta, phi_grid, nphi)


class UnstructuredMesh(Mesh):
    pass


_MESH_TYPE_MAP = {
    'regular': RegularMesh,
    'rectilinear': RectilinearMesh,
    'cylindrical': CylindricalMesh,
    'spherical': SphericalMesh,
    'unstructured': UnstructuredMesh
}


def _get_mesh(index):
    mesh_type = create_string_buffer(20)
    _dll.openmc_mesh_get_type(index, mesh_type)
    mesh_type = mesh_type.value.decode()
    return _MESH_TYPE_MAP[mesh_type](index=index)


class _MeshMapping(Mapping):
    def __getitem__(self, key):
        index = c_int32()
        try:
            _dll.openmc_get_mesh_index(key, index)
        except (AllocationError, InvalidIDError) as e:
            # __contains__ expects a KeyError to work correctly
            raise KeyError(str(e))
        return _get_mesh(index.value)


    def __iter__(self):
        for i in range(len(self)):
            yield _get_mesh(i).id

    def __len__(self):
        return _dll.n_meshes()

    def __repr__(self):
        return repr(dict(self))

meshes = _MeshMapping()
