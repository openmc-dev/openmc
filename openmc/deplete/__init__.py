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
    # check if running with MPI and if hdf5 is MPI-enabled
    from h5py import get_config

    if not get_config().mpi and comm.size > 1:
        # Raise exception only on process 0
        if comm.rank:
            from sys import exit

            exit()
        raise RuntimeError(
            "Need MPI-enabled HDF5 install to perform depletion with MPI"
        )
except ImportError:
    comm = DummyCommunicator()
    have_mpi = False

from .nuclide import *
from .chain import *
from .operator import *
from .reaction_rates import *
from .abc import *
from .results import *
from .results_list import *
from .integrator import *
