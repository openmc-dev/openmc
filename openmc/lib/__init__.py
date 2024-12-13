"""
This module provides bindings to C/C++ functions defined by OpenMC shared
library. When the :mod:`openmc.lib` package is imported, the OpenMC shared
library is automatically loaded. Calls to the OpenMC library can then be via
functions or objects in :mod:`openmc.lib`, for example:

.. code-block:: python

    openmc.lib.init()
    openmc.lib.run()
    openmc.lib.finalize()

"""

from ctypes import CDLL, c_bool, c_int
import os
import warnings

def load_openmc_library():
    try:
        # Attempt to load the library from OpenMC
        import openmc
        _filename = openmc.lib[0]
        if os.path.isfile(_filename):
            return CDLL(str(_filename))
        raise FileNotFoundError
    except (IndexError, FileNotFoundError):
        # Attempt to load the library from the installed module
        import importlib
        openmc = importlib.import_module("openmc")
        _filename = openmc.lib[0]
        if os.path.isfile(_filename):
            warnings.warn(
                "It seems OpenMC is being run from its source directory. "
                "This setup is not recommended as it may lead to unexpected behavior, "
                "such as conflicts between source and installed versions. "
                "Please run your script from outside the OpenMC source tree.",
                RuntimeWarning
            )
            return CDLL(str(_filename))
        raise RuntimeError(
            "Unable to load the OpenMC library. Please ensure OpenMC "
            "is installed correctly and accessible from your Python environment."
        )

# Load the OpenMC shared library
_dll = load_openmc_library()


def _dagmc_enabled():
    return c_bool.in_dll(_dll, "DAGMC_ENABLED").value

def _ncrystal_enabled():
    return c_bool.in_dll(_dll, "NCRYSTAL_ENABLED").value

def _coord_levels():
    return c_int.in_dll(_dll, "n_coord_levels").value

def _libmesh_enabled():
    return c_bool.in_dll(_dll, "LIBMESH_ENABLED").value

def _mcpl_enabled():
    return c_bool.in_dll(_dll, "MCPL_ENABLED").value

def _uwuw_enabled():
    return c_bool.in_dll(_dll, "UWUW_ENABLED").value


from .error import *
from .core import *
from .nuclide import *
from .material import *
from .cell import *
from .mesh import *
from .filter import *
from .tally import *
from .settings import settings
from .math import *
from .plot import *
from .weight_windows import *

# Flag to denote whether or not openmc.lib.init has been called
# TODO: Establish and use a flag in the C++ code to represent the status of the
# openmc_init and openmc_finalize methods
is_initialized = False
