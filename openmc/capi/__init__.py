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
import sys
from warnings import warn

import pkg_resources


# Determine shared-library suffix
if sys.platform == 'darwin':
    _suffix = 'dylib'
else:
    _suffix = 'so'

# Open shared library
_filename = pkg_resources.resource_filename(
    __name__, '_libopenmc.{}'.format(_suffix))
_dll = CDLL(_filename)

from .error import *
from .core import *
from .nuclide import *
from .material import *
from .cell import *
from .tally import *
