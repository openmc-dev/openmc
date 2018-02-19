"""function module.

This module contains the Operator class, which is then passed to an integrator
to run a full depletion simulation.
"""

from collections import namedtuple
import os
from pathlib import Path

from abc import ABCMeta, abstractmethod


class Settings(object):
    """The Settings class.

    Contains all parameters necessary for the integrator.

    Attributes
    ----------
    dt_vec : numpy.array
        Array of time steps to take.
    output_dir : pathlib.Path
        Path to output directory to save results.
    chain_file : str
        Path to the depletion chain XML file.  Defaults to the
        :envvar:`OPENMC_DEPLETE_CHAIN` environment variable if it exists.
    dilute_initial : float
        Initial atom density to add for nuclides that are zero in initial
        condition to ensure they exist in the decay chain.  Only done for
        nuclides with reaction rates. Defaults to 1.0e3.
    power : float
        Power of the reactor in [W]. For a 2D problem, the power can be given in
        W/cm as long as the "volume" assigned to a depletion material is
        actually an area in cm^2.

    """
    def __init__(self):
        try:
            self.chain_file = os.environ["OPENMC_DEPLETE_CHAIN"]
        except KeyError:
            self.chain_file = None
        self.dt_vec = None
        self.output_dir = '.'
        self.power = None
        self.dilute_initial = 1.0e3

    @property
    def output_dir(self):
        return self._output_dir

    @output_dir.setter
    def output_dir(self, output_dir):
        self._output_dir = Path(output_dir)


OperatorResult = namedtuple('OperatorResult', ['k', 'rates'])


class Operator(metaclass=ABCMeta):
    """Abstract class defining a transport operator

    Attributes
    ----------
    settings : Settings
        Settings object.

    """
    def __init__(self, settings):
        self.settings = settings

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
        if not self.settings.output_dir.exists():
            self.settings.output_dir.mkdir()  # exist_ok parameter is 3.5+

        # In Python 3.6+, chdir accepts a Path directly
        os.chdir(str(self.settings.output_dir))

        return self.initial_condition()

    def __exit__(self, exc_type, exc_value, traceback):
        self.finalize()
        os.chdir(self._orig_dir)

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

    @abstractmethod
    def form_matrix(self, y, mat):
        """Forms the f(y) matrix in y' = f(y)y.

        Nominally a depletion matrix, this is abstracted on the off chance
        that the function f has nothing to do with depletion at all.

        Parameters
        ----------
        y : numpy.ndarray
            An array representing y.
        mat : int
            Material id.

        Returns
        -------
        scipy.sparse.csr_matrix
            Sparse matrix representing f(y).
        """

        pass

    def finalize(self):
        pass
