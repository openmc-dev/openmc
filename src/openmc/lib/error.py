from ctypes import c_int, c_char
from warnings import warn

import openmc.exceptions as exc
from . import _dll


def _error_handler(err, func, args):
    """Raise exception according to error code."""

    # Get error code corresponding to global constant.
    def errcode(s):
        return c_int.in_dll(_dll, s).value

    # Get error message set by OpenMC library
    errmsg = (c_char*256).in_dll(_dll, 'openmc_err_msg')
    msg = errmsg.value.decode()

    # Raise exception type corresponding to error code
    if err == errcode('OPENMC_E_ALLOCATE'):
        raise exc.AllocationError(msg)
    elif err == errcode('OPENMC_E_OUT_OF_BOUNDS'):
        raise exc.OutOfBoundsError(msg)
    elif err == errcode('OPENMC_E_INVALID_ARGUMENT'):
        raise exc.InvalidArgumentError(msg)
    elif err == errcode('OPENMC_E_INVALID_TYPE'):
        raise exc.InvalidTypeError(msg)
    if err == errcode('OPENMC_E_INVALID_ID'):
        raise exc.InvalidIDError(msg)
    elif err == errcode('OPENMC_E_GEOMETRY'):
        raise exc.GeometryError(msg)
    elif err == errcode('OPENMC_E_DATA'):
        raise exc.DataError(msg)
    elif err == errcode('OPENMC_E_PHYSICS'):
        raise exc.PhysicsError(msg)
    elif err == errcode('OPENMC_E_WARNING'):
        warn(msg)
    elif err < 0:
        if not msg:
            msg = "Unknown error encountered (code {}).".format(err)
        raise exc.OpenMCError(msg)
