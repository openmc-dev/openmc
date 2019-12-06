from collections.abc import Mapping, Iterable
from ctypes import c_int, c_int32, c_double, POINTER
from weakref import WeakValueDictionary

import numpy as np
from numpy.ctypeslib import as_array

from openmc.exceptions import AllocationError, InvalidIDError
from . import _dll
from .core import _FortranObjectWithID
from .error import _error_handler
from .material import Material

__all__ = ['RegularMesh', 'meshes']

# Mesh functions
_dll.openmc_extend_meshes.argtypes = [c_int32, POINTER(c_int32), POINTER(c_int32)]
_dll.openmc_extend_meshes.restype = c_int
_dll.openmc_extend_meshes.errcheck = _error_handler
_dll.openmc_mesh_get_id.argtypes = [c_int32, POINTER(c_int32)]
_dll.openmc_mesh_get_id.restype = c_int
_dll.openmc_mesh_get_id.errcheck = _error_handler
_dll.openmc_mesh_get_dimension.argtypes = [c_int32, POINTER(POINTER(c_int)), POINTER(c_int)]
_dll.openmc_mesh_get_dimension.restype = c_int
_dll.openmc_mesh_get_dimension.errcheck = _error_handler
_dll.openmc_mesh_get_params.argtypes = [
    c_int32, POINTER(POINTER(c_double)), POINTER(POINTER(c_double)),
    POINTER(POINTER(c_double)), POINTER(c_int)]
_dll.openmc_mesh_get_params.restype = c_int
_dll.openmc_mesh_get_params.errcheck = _error_handler
_dll.openmc_mesh_set_id.argtypes = [c_int32, c_int32]
_dll.openmc_mesh_set_id.restype = c_int
_dll.openmc_mesh_set_id.errcheck = _error_handler
_dll.openmc_mesh_set_dimension.argtypes = [c_int32, c_int, POINTER(c_int)]
_dll.openmc_mesh_set_dimension.restype = c_int
_dll.openmc_mesh_set_dimension.errcheck = _error_handler
_dll.openmc_mesh_set_params.argtypes = [
    c_int32, c_int, POINTER(c_double), POINTER(c_double), POINTER(c_double)]
_dll.openmc_mesh_set_params.restype = c_int
_dll.openmc_mesh_set_params.errcheck = _error_handler
_dll.openmc_get_mesh_index.argtypes = [c_int32, POINTER(c_int32)]
_dll.openmc_get_mesh_index.restype = c_int
_dll.openmc_get_mesh_index.errcheck = _error_handler
_dll.n_meshes.argtypes = []
_dll.n_meshes.restype = c_int


class RegularMesh(_FortranObjectWithID):
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

                index = c_int32()
                _dll.openmc_extend_meshes(1, index, None)
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
    def dimension(self):
        dims = POINTER(c_int)()
        n = c_int()
        _dll.openmc_mesh_get_dimension(self._index, dims, n)
        return tuple(as_array(dims, (n.value,)))

    @dimension.setter
    def dimension(self, dimension):
        n = len(dimension)
        dimension = (c_int*n)(*dimension)
        _dll.openmc_mesh_set_dimension(self._index, n, dimension)

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
        _dll.openmc_mesh_get_params(self._index, ll, ur, w, n)
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
        _dll.openmc_mesh_set_params(self._index, n, lower_left, upper_right, width)


class _MeshMapping(Mapping):
    def __getitem__(self, key):
        index = c_int32()
        try:
            _dll.openmc_get_mesh_index(key, index)
        except (AllocationError, InvalidIDError) as e:
            # __contains__ expects a KeyError to work correctly
            raise KeyError(str(e))
        return RegularMesh(index=index.value)

    def __iter__(self):
        for i in range(len(self)):
            yield RegularMesh(index=i).id

    def __len__(self):
        return _dll.n_meshes()

    def __repr__(self):
        return repr(dict(self))

meshes = _MeshMapping()
