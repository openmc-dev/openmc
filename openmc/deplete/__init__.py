"""
openmc.deplete
==============

A depletion front-end tool.
"""
import sys
from unittest.mock import Mock

from .dummy_comm import DummyCommunicator

try:
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
except ImportError:
    MPI = Mock()
    comm = DummyCommunicator()

from .nuclide import *
from .chain import *
from .operator import *
from .reaction_rates import *
from .atom_number import *
from .results import *
from .results_list import *
from .integrators import *
from . import abc
from . import cram
from . import helpers
