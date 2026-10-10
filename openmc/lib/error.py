from ctypes import c_char_p
from warnings import warn

import openmc.exceptions as exc
from . import _dll


OPENMC_E_WARNING = 1
OPENMC_E_UNASSIGNED = -1
OPENMC_E_ALLOCATE = -2
OPENMC_E_OUT_OF_BOUNDS = -3
OPENMC_E_INVALID_SIZE = -4
OPENMC_E_INVALID_ARGUMENT = -5
OPENMC_E_INVALID_TYPE = -6
OPENMC_E_INVALID_ID = -7
OPENMC_E_GEOMETRY = -8
OPENMC_E_DATA = -9
OPENMC_E_PHYSICS = -10

_dll.openmc_get_err_msg.restype = c_char_p


def _error_handler(err, func, args):
    """Raise exception according to error code."""

    # Get error message set by OpenMC library
    msg = _dll.openmc_get_err_msg().decode()

    # Raise exception type corresponding to error code
    if err == OPENMC_E_ALLOCATE:
        raise exc.AllocationError(msg)
    elif err == OPENMC_E_OUT_OF_BOUNDS:
        raise exc.OutOfBoundsError(msg)
    elif err == OPENMC_E_INVALID_ARGUMENT:
        raise exc.InvalidArgumentError(msg)
    elif err == OPENMC_E_INVALID_TYPE:
        raise exc.InvalidTypeError(msg)
    if err == OPENMC_E_INVALID_ID:
        raise exc.InvalidIDError(msg)
    elif err == OPENMC_E_GEOMETRY:
        raise exc.GeometryError(msg)
    elif err == OPENMC_E_DATA:
        raise exc.DataError(msg)
    elif err == OPENMC_E_PHYSICS:
        raise exc.PhysicsError(msg)
    elif err == OPENMC_E_WARNING:
        warn(msg)
    elif err < 0:
        if not msg:
            msg = f"Unknown error encountered (code {err})."
        raise exc.OpenMCError(msg)
