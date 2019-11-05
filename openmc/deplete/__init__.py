"""
openmc.deplete
==============

A depletion front-end tool.
"""
import sys
from unittest.mock import Mock

from h5py import get_config

from .dummy_comm import DummyCommunicator

try:
    from mpi4py import MPI

    comm = MPI.COMM_WORLD
    have_mpi = True
    # check if running with MPI and if using parallel HDF5

    if not get_config().mpi and comm.size > 1:
        # Raise exception only on process 0
        if comm.rank:
            sys.exit()
        raise RuntimeError(
            "Need parallel HDF5 installed to perform depletion with MPI"
        )
except ImportError:
    comm = DummyCommunicator()
    have_mpi = False
    MPI = Mock()

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
