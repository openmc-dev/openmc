from ctypes import c_int

from . import _dll


class GeometryError(Exception):
    pass


def _error_handler(err, func, args):
    """Raise exception according to error code."""

    # Get error code corresponding to global constant.
    def errcode(s):
        return c_int.in_dll(_dll, s).value

    if err == errcode('e_out_of_bounds'):
        raise IndexError('Array index out of bounds.')

    elif err == errcode('e_cell_not_allocated'):
        raise MemoryError("Memory has not been allocated for cells.")

    elif err == errcode('e_cell_invalid_id'):
        raise KeyError("No cell exists with ID={}.".format(args[0]))

    elif err == errcode('e_cell_not_found'):
        raise GeometryError("Could not find cell at position ({}, {}, {})"
                            .format(*args[0]))

    elif err == errcode('e_nuclide_not_allocated'):
        raise MemoryError("Memory has not been allocated for nuclides.")

    elif err == errcode('e_nuclide_not_loaded'):
        raise KeyError("No nuclide named '{}' has been loaded.")

    elif err == errcode('e_nuclide_not_in_library'):
        raise KeyError("Specified nuclide doesn't exist in the cross "
                       "section library.")

    elif err == errcode('e_material_not_allocated'):
        raise MemoryError("Memory has not been allocated for materials.")

    elif err == errcode('e_material_invalid_id'):
        raise KeyError("No material exists with ID={}.".format(args[0]))

    elif err == errcode('e_tally_not_allocated'):
        raise MemoryError("Memory has not been allocated for tallies.")

    elif err == errcode('e_tally_invalid_id'):
        raise KeyError("No tally exists with ID={}.".format(args[0]))

    elif err == errcode('e_invalid_size'):
        raise MemoryError("Array size mismatch with memory allocated.")

    elif err == errcode('e_cell_no_material'):
        raise GeometryError("Operation on cell requires that it be filled"
                            " with a material.")

    elif err == errcode('w_below_min_bound'):
        warn("Data has not been loaded beyond lower bound of {}.".format(args[0]))

    elif err == errcode('w_above_max_bound'):
        warn("Data has not been loaded beyond upper bound of {}.".format(args[0]))

    elif err < 0:
        raise Exception("Unknown error encountered (code {}).".format(err))
