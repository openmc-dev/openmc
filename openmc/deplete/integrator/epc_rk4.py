"""The EPC-RK4 integrator."""

from .cram import timed_deplete
from .abc import Integrator


class EPC_RK4_Integrator(Integrator):
    r"""Deplete using the EPC-RK4 algorithm.

    Implements an extended predictor-corrector algorithm with traditional
    Runge-Kutta 4 method.
    This algorithm is mathematically defined as:

    .. math::
        \begin{aligned}
        F_1 &= h A(y_0) \\
        y_1 &= \text{expm}(1/2 F_1) y_0 \\
        F_2 &= h A(y_1) \\
        y_2 &= \text{expm}(1/2 F_2) y_0 \\
        F_3 &= h A(y_2) \\
        y_3 &= \text{expm}(F_3) y_0 \\
        F_4 &= h A(y_3) \\
        y_4 &= \text{expm}(1/6 F_1 + 1/3 F_2 + 1/3 F_3 + 1/6 F_4) y_0
        \end{aligned}
    """
    _num_stages = 4

    def __call__(self, conc, rates, dt, power, _i):
        """Perform the integration across one time step

        Parameters
        ----------
        conc : numpy.ndarray
            Initial concentrations for all nuclides in [atom]
        rates : openmc.deplete.ReactionRates
            Reaction rates from operator
        dt : float
            Time in [s] for the entire depletion interval
        power : float
            Power of the system in [W]
        i : int
            Current depletion step index, unused.

        Returns
        -------
        proc_time : float
            Time spent in CRAM routines for all materials in [s]
        conc_list : list of numpy.ndarray
            Concentrations at each of the intermediate points with
            the final concentration as the last element
        op_results : list of openmc.deplete.OperatorResult
            Eigenvalue and reaction rates from intermediate transport
            simulations
        """

        # Step 1: deplete with matrix A(y0) / 2
        time1, conc1 = timed_deplete(
            self.chain, conc, rates, dt, matrix_func=_rk4_f1)
        res1 = self.operator(conc1, power)

        # Step 2: deplete with matrix A(y1) / 2
        time2, conc2 = timed_deplete(
            self.chain, conc, res1.rates, dt, matrix_func=_rk4_f1)
        res2 = self.operator(conc2, power)

        # Step 3: deplete with matrix A(y2)
        time3, conc3 = timed_deplete(
            self.chain, conc, res2.rates, dt)
        res3 = self.operator(conc3, power)

        # Step 4: deplete with matrix built from weighted rates
        list_rates = list(zip(rates, res1.rates, res2.rates, res3.rates))
        time4, conc4 = timed_deplete(
            self.chain, conc, list_rates, dt, matrix_func=_rk4_f4)

        return (time1 + time2 + time3 + time4, [conc1, conc2, conc3, conc4],
                [res1, res2, res3])


# Functions to form the special matrix for depletion
def _rk4_f1(chain, rates):
    return 1/2 * chain.form_matrix(rates)


def _rk4_f4(chain, rates):
    return 1/6 * chain.form_matrix(rates[0]) + \
           1/3 * chain.form_matrix(rates[1]) + \
           1/3 * chain.form_matrix(rates[2]) + \
           1/6 * chain.form_matrix(rates[3])
