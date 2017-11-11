"""
This module provides bindings to C functions defined by OpenMC shared library.
When the :mod:`openmc` package is imported, the OpenMC shared library is
automatically loaded. Calls to the OpenMC library can then be via functions or
objects in the :mod:`openmc.capi` subpackage, for example:

.. code-block:: python

    openmc.capi.init()
    openmc.capi.run()
    openmc.capi.finalize()

"""

from ctypes import CDLL
import os
import sys
from warnings import warn

import pkg_resources


# Determine shared-library suffix
if sys.platform == 'darwin':
    _suffix = 'dylib'
else:
    _suffix = 'so'

if os.environ.get('READTHEDOCS', None) != 'True':
    # Open shared library
    _filename = pkg_resources.resource_filename(
        __name__, 'libopenmc.{}'.format(_suffix))
    _dll = CDLL(_filename)
else:
    # For documentation builds, we don't actually have the shared library
    # available. Instead, we create a mock object so that when the modules
    # within the openmc.capi package try to configure arguments and return
    # values for symbols, no errors occur
    try:
        from unittest.mock import Mock
    except ImportError:
        from mock import Mock
    _dll = Mock()

from .error import *
from .core import *
from .nuclide import *
from .material import *
from .cell import *
from .filter import *
from .tally import *
from .settings import settings

warn("The Python bindings to OpenMC's C API are still unstable "
     "and may change substantially in future releases.", FutureWarning)
