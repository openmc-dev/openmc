import numpy as np
import scipy.sparse as sp
from openmc.deplete.reaction_rates import ReactionRates
from openmc.deplete.abc import TransportOperator, OperatorResult


class DummyOperator(TransportOperator):
    """This is a dummy operator class with no statistical uncertainty.

    y_1' = sin(y_2) y_1 + cos(y_1) y_2
    y_2' = -cos(y_2) y_1 + sin(y_1) y_2

    y_1(0) = 1
    y_2(0) = 1

    y_1(1.5) ~ 2.3197067076743316
    y_2(1.5) ~ 3.1726475740397628

    """
    def __init__(self):
        pass

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
        return OperatorResult(0.0, reaction_rates)

    @property
    def chain(self):
        return self

    def form_matrix(self, rates):
        """Forms the f(y) matrix in y' = f(y)y.

        Nominally a depletion matrix, this is abstracted on the off chance
        that the function f has nothing to do with depletion at all.

        Parameters
        ----------
        rates : numpy.ndarray
            Slice of reaction rates for a single material

        Returns
        -------
        scipy.sparse.csr_matrix
            Sparse matrix representing f(y).
        """

        y_1 = rates[0, 0]
        y_2 = rates[1, 0]

        mat = np.zeros((2, 2))
        a11 = np.sin(y_2)
        a12 = np.cos(y_1)
        a21 = -np.cos(y_2)
        a22 = np.sin(y_1)

        return sp.csr_matrix(np.array([[a11, a12], [a21, a22]]))

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
            A list of all material IDs to be burned.  Used for sorting the simulation.
        """

        return ["1"]

    @property
    def burnable_mats(self):
        """Maps cell name to index in global geometry."""
        return self.local_mats


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
            A list of all cell IDs to be burned.  Used for sorting the simulation.
        full_burn_list : OrderedDict of str to int
            Maps cell name to index in global geometry.

        """
        return self.volume, self.nuc_list, self.local_mats, self.burnable_mats
