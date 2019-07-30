"""The SI-LE/QI CFQ4 integrator."""

import copy
from collections.abc import Iterable
from itertools import repeat

from uncertainties import ufloat

from .abc import SI_Integrator
from .si_celi import SI_CELI_Integrator
from .leqi import _leqi_f1, _leqi_f2, _leqi_f3, _leqi_f4
from .cram import timed_deplete
from ..results import Results
from ..abc import OperatorResult


class SI_LEQI_Integrator(SI_Integrator):
    r"""Deplete using the SI-LE/QI CFQ4 algorithm.

    Implements the Stochastic Implicit LE/QI Predictor-Corrector algorithm using
    the `fourth order commutator-free integrator <https://doi.org/10.1137/05063042>`_.

    Detailed algorithm can be found in Section 3.2 in `Colin Josey's thesis
    <http://hdl.handle.net/1721.1/113721>`_.
    """

    def __call__(self, bos_conc, bos_rates, dt, power, i):
        """Perform the integration across one time step

        Parameters
        ----------
        bos_conc : list of numpy.ndarray
            Initial concentrations for all nuclides in [atom] for
            all depletable materials
        bos_rates : list of openmc.deplete.ReactionRates
            Reaction rates from operator for all depletable materials
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
            simulation
        """
        if i == 0:
            if self._ires < 1:
                self._prev_rates = bos_rates
                # Perform CELI for initial steps
                return SI_CELI_Integrator.__call__(
                    self, bos_conc, bos_rates, dt, power, i)
            prev_res = self.operator.prev_res[-2]
            prevdt = self.timesteps[i] - prev_res.time[0]
            self._prev_rates = prev_res.rates[0]
        else:
            prevdt = self.timesteps[i - 1]

        # Perform remaining LE/QI
        inputs = list(zip(self._prev_rates, bos_rates,
                          repeat(prevdt), repeat(dt)))
        proc_time, inter_conc = timed_deplete(
            self.chain, bos_conc, inputs, dt, matrix_func=_leqi_f1)
        time1, eos_conc = timed_deplete(
            self.chain, inter_conc, inputs, dt, matrix_func=_leqi_f2)

        proc_time += time1
        inter_conc = copy.deepcopy(eos_conc)

        for j in range(self.n_stages + 1):
            inter_res = self.operator(inter_conc, power)

            if j <= 1:
                res_bar = copy.deepcopy(inter_res)
            else:
                rates = 1 / j * inter_res.rates + (1 - 1 / j) * res_bar.rates
                k = 1 / j * inter_res.k + (1 - 1 / j) * res_bar.k
                res_bar = OperatorResult(k, rates)

            inputs = list(zip(self._prev_rates, bos_rates, res_bar.rates,
                              repeat(prevdt), repeat(dt)))
            time1, inter_conc = timed_deplete(
                self.chain, bos_conc, inputs, dt, matrix_func=_leqi_f3)
            time2, inter_conc = timed_deplete(
                self.chain, inter_conc, inputs, dt, matrix_func=_leqi_f4)
            proc_time += time1 + time2

        return proc_time, [eos_conc, inter_conc], [res_bar]
