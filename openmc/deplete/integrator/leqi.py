"""The LE/QI CFQ4 integrator."""

import copy
from itertools import repeat

from .abc import Integrator
from .celi import CELIIntegrator
from .cram import timed_deplete


class LEQIIntegrator(Integrator):
    r"""Deplete using the LE/QI CFQ4 algorithm.

    Implements the LE/QI Predictor-Corrector algorithm using the `fourth order
    commutator-free integrator <https://doi.org/10.1137/05063042>`_.

    "LE/QI" stands for linear extrapolation on predictor and quadratic
    interpolation on corrector. This algorithm is mathematically defined as:

    .. math::
        \begin{aligned}
        y' &= A(y, t) y(t) \\
        A_{last} &= A(y_{n-1}, t_n - h_1) \\
        A_0 &= A(y_n, t_n) \\
        F_1 &= \frac{-h_2^2}{12h_1} A_{last} + \frac{h_2(6h_1+h_2)}{12h_1} A_0 \\
        F_2 &= \frac{-5h_2^2}{12h_1} A_{last} + \frac{h_2(6h_1+5h_2)}{12h_1} A_0 \\
        y_p &= \text{expm}(F_2) \text{expm}(F_1) y_n \\
        A_1 &= A(y_p, t_n + h_2) \\
        F_3 &= \frac{-h_2^3}{12 h_1 (h_1 + h_2)} A_{last} +
              \frac{h_2 (5 h_1^2 + 6 h_2 h_1 + h_2^2)}{12 h_1 (h_1 + h_2)} A_0 +
              \frac{h_2 h_1)}{12 (h_1 + h_2)} A_1 \\
        F_4 &= \frac{-h_2^3}{12 h_1 (h_1 + h_2)} A_{last} +
              \frac{h_2 (h_1^2 + 2 h_2 h_1 + h_2^2)}{12 h_1 (h_1 + h_2)} A_0 +
              \frac{h_2 (5 h_1^2 + 4 h_2 h_1)}{12 h_1 (h_1 + h_2)} A_1 \\
        y_{n+1} &= \text{expm}(F_4) \text{expm}(F_3) y_n
        \end{aligned}

    It is initialized using the CE/LI algorithm.
    """
    _N_STAGES = 2

    def __call__(self, bos_conc, bos_rates, dt, power, i):
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
            Current depletion step index

        Returns
        -------
        proc_time : float
            Time spent in CRAM routines for all materials in [s]
        conc_list : list of numpy.ndarray
            Concentrations at each of the intermediate points with
            the final concentration as the last element
        op_results : list of openmc.deplete.OperatorResult
            Eigenvalue and reaction rates from intermediate transport
            simulation
        """
        if i == 0:
            if self._ires < 1:  # need at least previous transport solution
                self._prev_rates = bos_rates
                return CELIIntegrator.__call__(
                    self, bos_conc, bos_rates, dt, power, i)
            prev_res = self.operator.prev_res[-2]
            prevdt = self.timesteps[i] - prev_res.time[0]
            self._prev_rates = prev_res.rates[0]
        else:
            prevdt = self.timesteps[i - 1]

        # Remaining LE/QI
        bos_res = self.operator(bos_conc, power)

        le_inputs = list(zip(
            self._prev_rates, bos_res.rates, repeat(prevdt), repeat(dt)))

        time1, conc_inter = timed_deplete(
            self.chain, bos_conc, le_inputs, dt, matrix_func=_leqi_f1)
        time2, conc_eos0 = timed_deplete(
            self.chain, conc_inter, le_inputs, dt, matrix_func=_leqi_f2)

        res_inter = self.operator(conc_eos0, power)

        qi_inputs = list(zip(
            self._prev_rates, bos_res.rates, res_inter.rates,
            repeat(prevdt), repeat(dt)))

        time3, conc_inter = timed_deplete(
            self.chain, bos_conc, qi_inputs, dt, matrix_func=_leqi_f3)
        time4, conc_eos1 = timed_deplete(
            self.chain, conc_inter, qi_inputs, dt, matrix_func=_leqi_f4)

        # store updated rates
        self._prev_rates = copy.deepcopy(bos_res.rates)

        return (
            time1 + time2 + time3 + time4, [conc_eos0, conc_eos1],
            [bos_res, res_inter])


# Functions to form the special matrix for depletion
def _leqi_f1(chain, inputs):
    f1 = chain.form_matrix(inputs[0])
    f2 = chain.form_matrix(inputs[1])
    dt_l, dt = inputs[2], inputs[3]
    return -dt / (12 * dt_l) * f1 + (dt + 6 * dt_l) / (12 * dt_l) * f2


def _leqi_f2(chain, inputs):
    f1 = chain.form_matrix(inputs[0])
    f2 = chain.form_matrix(inputs[1])
    dt_l, dt = inputs[2], inputs[3]
    return -5 * dt / (12 * dt_l) * f1 + (5 * dt + 6 * dt_l) / (12 * dt_l) * f2


def _leqi_f3(chain, inputs):
    f1 = chain.form_matrix(inputs[0])
    f2 = chain.form_matrix(inputs[1])
    f3 = chain.form_matrix(inputs[2])
    dt_l, dt = inputs[3], inputs[4]
    return -dt**2 / (12 * dt_l * (dt + dt_l)) * f1 + \
           (dt**2 + 6*dt*dt_l + 5*dt_l**2) / (12 * dt_l * (dt + dt_l)) * f2 + \
           dt_l / (12 * (dt + dt_l)) * f3


def _leqi_f4(chain, inputs):
    f1 = chain.form_matrix(inputs[0])
    f2 = chain.form_matrix(inputs[1])
    f3 = chain.form_matrix(inputs[2])
    dt_l, dt = inputs[3], inputs[4]
    return -dt**2 / (12 * dt_l * (dt + dt_l)) * f1 + \
           (dt**2 + 2*dt*dt_l + dt_l**2) / (12 * dt_l * (dt + dt_l)) * f2 + \
           (4 * dt * dt_l + 5 * dt_l**2) / (12 * dt_l * (dt + dt_l)) * f3
