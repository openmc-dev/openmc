"""The CF4 integrator."""

import copy
from collections.abc import Iterable

from .cram import timed_deplete
from .abc import Integrator
from ..results import Results


class CF4Integrator(Integrator):
    r"""Deplete using the CF4 algorithm.

    Implements the fourth order `commutator-free Lie algorithm
    <https://doi.org/10.1016/S0167-739X(02)00161-9>`_.
    This algorithm is mathematically defined as:

    .. math::
        \begin{aligned}
        F_1 &= h A(y_0) \\
        y_1 &= \text{expm}(1/2 F_1) y_0 \\
        F_2 &= h A(y_1) \\
        y_2 &= \text{expm}(1/2 F_2) y_0 \\
        F_3 &= h A(y_2) \\
        y_3 &= \text{expm}(-1/2 F_1 + F_3) y_1 \\
        F_4 &= h A(y_3) \\
        y_4 &= \text{expm}( 1/4  F_1 + 1/6 F_2 + 1/6 F_3 - 1/12 F_4)
               \text{expm}(-1/12 F_1 + 1/6 F_2 + 1/6 F_3 + 1/4  F_4) y_0
        \end{aligned}
    """
    _N_STAGES = 4

    def __call__(self, bos_conc, bos_rates, dt, power, i):
        """Perform the integration across one time step

        Parameters
        ----------
        bos_conc : numpy.ndarray
            Initial concentrations for all nuclides in [atom]
        bos_rates : openmc.deplete.ReactionRates
            Reaction rates from operator
        dt : float
            Time in [s] for the entire depletion interval
        power : float
            Power of the system [W]
        i : int
            Current depletion step index

        Returns
        -------
        proc_time : float
            Time spent in CRAM routines for all materials
        conc_list : list of numpy.ndarray
            Concentrations at each of the intermediate points with
            the final concentration as the last element
        op_results : list of openmc.deplete.OperatorResult
            Eigenvalue and reaction rates from intermediate transport
            simulations
        """
        # Step 1: deplete with matrix 1/2*A(y0)
        time1, conc_eos1 = timed_deplete(
            self.chain, bos_conc, bos_rates, dt, matrix_func=_cf4_f1)
        res1 = self.operator(conc_eos1, power)

        # Step 2: deplete with matrix 1/2*A(y1)
        time2, conc_eos2 = timed_deplete(
            self.chain, bos_conc, res1.rates, dt, matrix_func=_cf4_f1)
        res2 = self.operator(conc_eos2, power)

        # Step 3: deplete with matrix -1/2*A(y0)+A(y2)
        list_rates = list(zip(bos_rates, res2.rates))
        time3, conc_eos3 = timed_deplete(
            self.chain, conc_eos1, list_rates, dt, matrix_func=_cf4_f2)
        res3 = self.operator(conc_eos3, power)

        # Step 4: deplete with two matrix exponentials
        list_rates = list(zip(bos_rates, res1.rates, res2.rates, res3.rates))
        time4, conc_inter = timed_deplete(
            self.chain, bos_conc, list_rates, dt, matrix_func=_cf4_f3)
        time5, conc_eos5 = timed_deplete(
            self.chain, conc_inter, list_rates, dt, matrix_func=_cf4_f4)

        return (time1 + time2 + time3 + time4 + time5,
                [conc_eos1, conc_eos2, conc_eos3, conc_eos5],
                [res1, res2, res3])


# Functions to form the special matrix for depletion
def _cf4_f1(chain, rates):
    return 1/2 * chain.form_matrix(rates)


def _cf4_f2(chain, rates):
    return -1/2 * chain.form_matrix(rates[0]) + \
                  chain.form_matrix(rates[1])


def _cf4_f3(chain, rates):
    return  1/4  * chain.form_matrix(rates[0]) + \
            1/6  * chain.form_matrix(rates[1]) + \
            1/6  * chain.form_matrix(rates[2]) + \
           -1/12 * chain.form_matrix(rates[3])


def _cf4_f4(chain, rates):
    return -1/12 * chain.form_matrix(rates[0]) + \
            1/6  * chain.form_matrix(rates[1]) + \
            1/6  * chain.form_matrix(rates[2]) + \
            1/4  * chain.form_matrix(rates[3])
