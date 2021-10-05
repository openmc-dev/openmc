try:
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
except ImportError:
    from unittest.mock import Mock
    MPI = Mock()
    from openmc.dummy_comm import DummyCommunicator
    comm = DummyCommunicator()
