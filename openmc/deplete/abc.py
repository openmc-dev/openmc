"""function module.

This module contains the Operator class, which is then passed to an integrator
to run a full depletion simulation.
"""

from collections import namedtuple
import os
from pathlib import Path
from abc import ABCMeta, abstractmethod

from .chain import Chain

OperatorResult = namedtuple('OperatorResult', ['k', 'rates'])


class Operator(metaclass=ABCMeta):
    """Abstract class defining a transport operator

    Parameters
    ----------
    chain_file : str, optional


    Attributes
    ----------
    dilute_initial : float
        Initial atom density to add for nuclides that are zero in initial
        condition to ensure they exist in the decay chain.  Only done for
        nuclides with reaction rates. Defaults to 1.0e3.

    """
    def __init__(self, chain_file=None):
        self.dilute_initial = 1.0e3
        self.output_dir = '.'

        # Read depletion chain
        if chain_file is None:
            chain_file = os.environ.get("OPENMC_DEPLETE_CHAIN", None)
            if chain_file is None:
                raise IOError("No chain specified, either manually or in "
                              "environment variable OPENMC_DEPLETE_CHAIN.")
        self.chain = Chain.from_xml(chain_file)

    @abstractmethod
    def __call__(self, vec, print_out=True):
        """Runs a simulation.

        Parameters
        ----------
        vec : list of numpy.array
            Total atoms to be used in function.
        print_out : bool, optional
            Whether or not to print out time.

        Returns
        -------
        openmc.deplete.OperatorResult
            Eigenvalue and reaction rates resulting from transport operator

        """
        pass

    def __enter__(self):
        # Save current directory and move to specific output directory
        self._orig_dir = os.getcwd()
        if not self.output_dir.exists():
            self.output_dir.mkdir()  # exist_ok parameter is 3.5+

        # In Python 3.6+, chdir accepts a Path directly
        os.chdir(str(self.output_dir))

        return self.initial_condition()

    def __exit__(self, exc_type, exc_value, traceback):
        self.finalize()
        os.chdir(self._orig_dir)

    @property
    def output_dir(self):
        return self._output_dir

    @output_dir.setter
    def output_dir(self, output_dir):
        self._output_dir = Path(output_dir)

    @abstractmethod
    def initial_condition(self):
        """Performs final setup and returns initial condition.

        Returns
        -------
        list of numpy.array
            Total density for initial conditions.
        """

        pass

    @abstractmethod
    def get_results_info(self):
        """Returns volume list, cell lists, and nuc lists.

        Returns
        -------
        volume : list of float
            Volumes corresponding to materials in burn_list
        nuc_list : list of str
            A list of all nuclide names. Used for sorting the simulation.
        burn_list : list of int
            A list of all cell IDs to be burned.  Used for sorting the simulation.
        full_burn_list : list of int
            All burnable materials in the geometry.
        """

        pass

    def finalize(self):
        pass
