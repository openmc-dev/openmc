from collections import namedtuple
from unittest.mock import Mock

import numpy as np
import scipy.sparse as sp
from uncertainties import ufloat

from openmc.deplete.reaction_rates import ReactionRates
from openmc.deplete.abc import TransportOperator, OperatorResult
from openmc.deplete import (
    CECMIntegrator, PredictorIntegrator, CELIIntegrator, LEQIIntegrator,
    EPCRK4Integrator, CF4Integrator, SICELIIntegrator, SILEQIIntegrator
)

# Bundle for nicely passing test data to depletion unit tests
# solver should be a concrete subclass of openmc.deplete.abc.Integrator
# atoms_1 should be the number of atoms of type 1 through the simulation
# similar for atoms_2, but for type 2. This includes the first step
# Solutions should be the exact solution that can be obtained using
# the DummyOperator depletion matrix with two 0.75 second time steps
DepletionSolutionTuple = namedtuple(
    "DepletionSolutionTuple", "solver atoms_1 atoms_2")


predictor_solution = DepletionSolutionTuple(
    PredictorIntegrator, np.array([1.0, 2.46847546272295, 4.11525874568034]),
    np.array([1.0, 0.986431226850467, -0.0581692232513460]))


cecm_solution = DepletionSolutionTuple(
    CECMIntegrator, np.array([1.0, 1.86872629872102, 2.18097439443550]),
    np.array([1.0, 1.395525772416039, 2.69429754646747]))


cf4_solution = DepletionSolutionTuple(
    CF4Integrator, np.array([1.0, 2.06101629, 2.57241318]),
    np.array([1.0, 1.37783588, 2.63731630]))


epc_rk4_solution = DepletionSolutionTuple(
    EPCRK4Integrator, np.array([1.0, 2.01978516, 2.05246421]),
    np.array([1.0, 1.42038037, 3.06177191]))


celi_solution = DepletionSolutionTuple(
    CELIIntegrator, np.array([1.0, 1.82078767, 2.68441779]),
    np.array([1.0, 0.97122898, 0.05125966]))


si_celi_solution = DepletionSolutionTuple(
    SICELIIntegrator, np.array([1.0, 2.03325094, 2.69291933]),
    np.array([1.0, 1.16826254, 0.37907772]))


leqi_solution = DepletionSolutionTuple(
    LEQIIntegrator, np.array([1.0, 1.82078767, 2.74526197]),
    np.array([1.0, 0.97122898, 0.23339915]))


si_leqi_solution = DepletionSolutionTuple(
    SILEQIIntegrator, np.array([1.0, 2.03325094, 2.92711288]),
    np.array([1.0, 1.16826254, 0.53753236]))


SCHEMES = {
    "predictor": predictor_solution,
    "cecm": cecm_solution,
    "celi": celi_solution,
    "cf4": cf4_solution,
    "epc_rk4": epc_rk4_solution,
    "leqi": leqi_solution,
    "si_leqi": si_leqi_solution,
    "si_celi": si_celi_solution,
}


class TestChain:
    """Empty chain to assist with unit testing depletion routines

    Only really provides the form_matrix function, but acts like
    a real Chain
    """

    fission_yields = [None]

    @staticmethod
    def get_default_fission_yields():
        return None

    def form_matrix(self, rates, _fission_yields=None):
        """Forms the f(y) matrix in y' = f(y)y.

        Nominally a depletion matrix, this is abstracted on the off chance
        that the function f has nothing to do with depletion at all.

        Parameters
        ----------
        rates : numpy.ndarray
            Slice of reaction rates for a single material
        _fission_yields : optional
            Not used

        Returns
        -------
        scipy.sparse.csr_matrix
            Sparse matrix representing f(y).
        """

        y_1 = rates[0, 0]
        y_2 = rates[1, 0]

        a11 = np.sin(y_2)
        a12 = np.cos(y_1)
        a21 = -np.cos(y_2)
        a22 = np.sin(y_1)

        return sp.csr_matrix(np.array([[a11, a12], [a21, a22]]))


class DummyOperator(TransportOperator):
    """This is a dummy operator class with no statistical uncertainty.

    y_1' = sin(y_2) y_1 + cos(y_1) y_2
    y_2' = -cos(y_2) y_1 + sin(y_1) y_2

    y_1(0) = 1
    y_2(0) = 1

    y_1(1.5) ~ 2.3197067076743316
    y_2(1.5) ~ 3.1726475740397628

    """
    def __init__(self, previous_results=None):
        self.prev_res = previous_results
        self.chain = TestChain()
        self.output_dir = "."
        self.settings = Mock()
        self.settings.particles = 10

    def __call__(self, vec, power, print_out=False):
        """Evaluates F(y)

        Parameters
        ----------
        vec : list of numpy.array
            Total atoms to be used in function.
        power : float
            Power in [W]
        print_out : bool, optional, ignored
            Whether or not to print out time.

        Returns
        -------
        openmc.deplete.OperatorResult
            Result of transport operator

        """
        mats = ["1"]
        nuclides = ["1", "2"]
        reactions = ["1"]

        reaction_rates = ReactionRates(mats, nuclides, reactions)

        reaction_rates[0, 0, 0] = vec[0][0]
        reaction_rates[0, 1, 0] = vec[0][1]

        # Create a fake rates object
        return OperatorResult(ufloat(0.0, 0.0), reaction_rates)

    @property
    def volume(self):
        """
        volume : dict of str float
            Volumes of material
        """

        return {"1": 0.0}

    @property
    def nuc_list(self):
        """
        nuc_list : list of str
            A list of all nuclide names. Used for sorting the simulation.
        """

        return ["1", "2"]

    @property
    def local_mats(self):
        """
        local_mats : list of str
            A list of all material IDs to be burned.  Used for sorting the
            simulation.
        """

        return ["1"]

    @property
    def burnable_mats(self):
        """Maps cell name to index in global geometry."""
        return self.local_mats

    @staticmethod
    def write_bos_data(_step):
        """Empty method but avoids calls to C API"""

    @property
    def reaction_rates(self):
        """
        reaction_rates : ReactionRates
            Reaction rates from the last operator step.
        """
        mats = ["1"]
        nuclides = ["1", "2"]
        reactions = ["1"]

        return ReactionRates(mats, nuclides, reactions)

    def initial_condition(self):
        """Returns initial vector.

        Returns
        -------
        list of numpy.array
            Total density for initial conditions.
        """

        return [np.array((1.0, 1.0))]

    def get_results_info(self):
        """Returns volume list, cell lists, and nuc lists.

        Returns
        -------
        volume : dict of str float
            Volumes corresponding to materials in full_burn_dict
        nuc_list : list of str
            A list of all nuclide names. Used for sorting the simulation.
        burn_list : list of int
            A list of all cell IDs to be burned.  Used for sorting the
            simulation.
        full_burn_list : OrderedDict of str to int
            Maps cell name to index in global geometry.

        """
        return self.volume, self.nuc_list, self.local_mats, self.burnable_mats
