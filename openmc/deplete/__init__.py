"""
OpenDeplete
===========

A simple depletion front-end tool.
"""

from .dummy_comm import DummyCommunicator
try:
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    have_mpi = True
except ImportError:
    comm = DummyCommunicator()
    have_mpi = False

from .nuclide import *
from .depletion_chain import *
from .openmc_wrapper import *
from .reaction_rates import *
from .function import *
from .results import *
from .integrator import *
from .utilities import *
