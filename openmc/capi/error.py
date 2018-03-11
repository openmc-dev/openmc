from ctypes import c_int, c_char
from warnings import warn

from . import _dll


class OpenMCError(Exception):
    """Root exception class for OpenMC."""


class GeometryError(OpenMCError):
    """Geometry-related error"""


class InvalidIDError(OpenMCError):
    """Use of an ID that is invalid."""


class AllocationError(OpenMCError):
    """Error related to memory allocation."""


class OutOfBoundsError(OpenMCError):
    """Index in array out of bounds."""


class DataError(OpenMCError):
    """Error relating to nuclear data."""


class PhysicsError(OpenMCError):
    """Error relating to performing physics."""


class InvalidArgumentError(OpenMCError):
    """Argument passed was invalid."""


class InvalidTypeError(OpenMCError):
    """Tried to perform an operation on the wrong type."""


def _error_handler(err, func, args):
    """Raise exception according to error code."""

    # Get error code corresponding to global constant.
    def errcode(s):
        return c_int.in_dll(_dll, s).value

    # Get error message set by OpenMC library
    errmsg = (c_char*256).in_dll(_dll, 'openmc_err_msg')
    msg = errmsg.value.decode()

    # Raise exception type corresponding to error code
    if err == errcode('e_allocate'):
        raise AllocationError(msg)
    elif err == errcode('e_out_of_bounds'):
        raise OutOfBoundsError(msg)
    elif err == errcode('e_invalid_argument'):
        raise InvalidArgumentError(msg)
    elif err == errcode('e_invalid_type'):
        raise InvalidTypeError(msg)
    if err == errcode('e_invalid_id'):
        raise InvalidIDError(msg)
    elif err == errcode('e_geometry'):
        raise GeometryError(msg)
    elif err == errcode('e_data'):
        raise DataError(msg)
    elif err == errcode('e_physics'):
        raise PhysicsError(msg)
    elif err == errcode('e_warning'):
        warn(msg)
    elif err < 0:
        raise OpenMCError("Unknown error encountered (code {}).".format(err))
