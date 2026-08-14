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

from ctypes import CDLL, byref, c_bool, c_char_p, c_int, POINTER
import importlib.resources
import os
import sys


# Determine shared-library suffix
if sys.platform == 'darwin':
    _suffix = 'dylib'
else:
    _suffix = 'so'

if os.environ.get('READTHEDOCS', None) != 'True':
    # Open shared library
    _filename = importlib.resources.files(__name__) / f'libopenmc.{_suffix}'
    _dll = CDLL(str(_filename))  # TODO: Remove str() when Python 3.12+
else:
    # For documentation builds, we don't actually have the shared library
    # available. Instead, we create a mock object so that when the modules
    # within the openmc.lib package try to configure arguments and return
    # values for symbols, no errors occur
    from unittest.mock import Mock
    _dll = Mock()

from .error import _error_handler

_dll.openmc_get_feature_enabled.argtypes = [c_char_p, POINTER(c_bool)]
_dll.openmc_get_feature_enabled.restype = c_int
_dll.openmc_get_feature_enabled.errcheck = _error_handler

def feature_enabled(feature: str) -> bool:
    """Return whether OpenMC was built with an optional feature.

    Parameters
    ----------
    feature : {'dagmc', 'libmesh', 'strict_fp', 'uwuw'}
        Feature to query.

    Returns
    -------
    bool
        Whether the feature is enabled.

    Raises
    ------
    InvalidArgumentError
        If *feature* is not recognized.

    """
    enabled = c_bool()
    _dll.openmc_get_feature_enabled(feature.encode(), byref(enabled))
    return enabled.value


def _coord_levels():
    return c_int.in_dll(_dll, "n_coord_levels").value


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
from .dagmc import *

# Flag to denote whether or not openmc.lib.init has been called
# TODO: Establish and use a flag in the C++ code to represent the status of the
# openmc_init and openmc_finalize methods
is_initialized = False
