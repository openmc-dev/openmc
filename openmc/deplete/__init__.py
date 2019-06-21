"""
openmc.deplete
==============

A depletion front-end tool.
"""

from .dummy_comm import DummyCommunicator
try:
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    have_mpi = True
    mpi_sum = MPI.SUM
except ImportError:
    comm = DummyCommunicator()
    have_mpi = False
    mpi_sum = lambda x: x

from .nuclide import *
from .chain import *
from .operator import *
from .reaction_rates import *
from .abc import *
from .results import *
from .results_list import *
from .integrator import *
