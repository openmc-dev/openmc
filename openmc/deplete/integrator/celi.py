"""The CE/LI CFQ4 integrator."""

from .cram import timed_deplete
from .abc import Integrator


class CELIIntegrator(Integrator):
    r"""Deplete using the CE/LI CFQ4 algorithm.

    Implements the CE/LI Predictor-Corrector algorithm using the `fourth order
    commutator-free integrator <https://doi.org/10.1137/05063042>`_.

    "CE/LI" stands for constant extrapolation on predictor and linear
    interpolation on corrector. This algorithm is mathematically defined as:

    .. math::
        \begin{aligned}
        y' &= A(y, t) y(t) \\
        A_0 &= A(y_n, t_n) \\
        y_p &= \text{expm}(h A_0) y_n \\
        A_1 &= A(y_p, t_n + h) \\
        y_{n+1} &= \text{expm}(\frac{h}{12} A_0 + \frac{5h}{12} A1)
                   \text{expm}(\frac{5h}{12} A_0 + \frac{h}{12} A1) y_n
        \end{aligned}
    """
    _N_STAGES = 2

    def __call__(self, bos_conc, rates, dt, power, _i=-1):
        """Perform the integration across one time step

        Parameters
        ----------
        bos_conc : numpy.ndarray
            Initial concentrations for all nuclides in [atom]
        rates : openmc.deplete.ReactionRates
            Reaction rates from operator
        dt : float
            Time in [s] for the entire depletion interval
        power : float
            Power of the system [W]
        _i : int
            Current iteration count. Not used

        Returns
        -------
        proc_time : float
            Time spent in CRAM routines for all materials
        conc_list : list of numpy.ndarray
            Concentrations at each of the intermediate points with
            the final concentration as the last element
        op_results : list of openmc.deplete.OperatorResult
            Eigenvalue and reaction rates from intermediate transport
            simulation
        """
        # deplete to end using BOS rates
        proc_time, conc_ce = timed_deplete(self.chain, bos_conc, rates, dt)
        res_ce = self.operator(conc_ce, power)

        # deplete using two matrix exponeitials
        list_rates = list(zip(rates, res_ce.rates))

        time_le1, conc_inter = timed_deplete(
            self.chain, bos_conc, list_rates, dt, matrix_func=_celi_f1)

        time_le2, conc_end = timed_deplete(
            self.chain, conc_inter, list_rates, dt, matrix_func=_celi_f2)

        return proc_time + time_le1 + time_le1, [conc_ce, conc_end], [res_ce]


# Functions to form the special matrix for depletion
def _celi_f1(chain, rates):
    return 5/12 * chain.form_matrix(rates[0]) + \
           1/12 * chain.form_matrix(rates[1])


def _celi_f2(chain, rates):
    return 1/12 * chain.form_matrix(rates[0]) + \
           5/12 * chain.form_matrix(rates[1])
