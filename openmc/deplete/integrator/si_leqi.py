"""The SI-LE/QI CFQ4 integrator."""

import copy
from itertools import repeat

from .abc import SIIntegrator
from .si_celi import SICELIIntegrator
from .leqi import _leqi_f1, _leqi_f2, _leqi_f3, _leqi_f4
from .cram import timed_deplete
from ..abc import OperatorResult


class SILEQIIntegrator(SIIntegrator):
    r"""Deplete using the SI-LE/QI CFQ4 algorithm.

    Implements the Stochastic Implicit LE/QI Predictor-Corrector algorithm using
    the `fourth order commutator-free integrator <https://doi.org/10.1137/05063042>`_.

    Detailed algorithm can be found in Section 3.2 in `Colin Josey's thesis
    <http://hdl.handle.net/1721.1/113721>`_.

    Parameters
    ----------
    operator : openmc.deplete.TransportOperator
        The operator object to simulate on.
    timesteps : iterable of float
        Array of timesteps in units of [s]. Note that values are not
        cumulative.
    power : float or iterable of float, optional
        Power of the reactor in [W]. A single value indicates that
        the power is constant over all timesteps. An iterable
        indicates potentially different power levels for each timestep.
        For a 2D problem, the power can be given in [W/cm] as long
        as the "volume" assigned to a depletion material is actually
        an area in [cm^2]. Either ``power`` or ``power_density`` must be
        specified.
    power_density : float or iterable of float, optional
        Power density of the reactor in [W/gHM]. It is multiplied by
        initial heavy metal inventory to get total power if ``power``
        is not speficied.
    n_steps : int, optional
        Number of stochastic iterations per depletion interval.
        Must be greater than zero. Default : 10

    Attributes
    ----------
    operator : openmc.deplete.TransportOperator
        Operator to perform transport simulations
    chain : openmc.deplete.Chain
        Depletion chain
    timesteps : iterable of float
        Size of each depletion interval in [s]
    power : iterable of float
        Power of the reactor in [W] for each interval in :attr:`timesteps`
    n_steps : int
        Number of stochastic iterations per depletion interval
    """
    _num_stages = 2

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
            if self._i_res < 1:
                self._prev_rates = bos_rates
                # Perform CELI for initial steps
                return SICELIIntegrator.__call__(
                    self, bos_conc, bos_rates, dt, power, i)
            prev_res = self.operator.prev_res[-2]
            prev_dt = self.timesteps[i] - prev_res.time[0]
            self._prev_rates = prev_res.rates[0]
        else:
            prev_dt = self.timesteps[i - 1]

        # Perform remaining LE/QI
        inputs = list(zip(self._prev_rates, bos_rates,
                          repeat(prev_dt), repeat(dt)))
        proc_time, inter_conc = timed_deplete(
            self.chain, bos_conc, inputs, dt, matrix_func=_leqi_f1)
        time1, eos_conc = timed_deplete(
            self.chain, inter_conc, inputs, dt, matrix_func=_leqi_f2)

        proc_time += time1
        inter_conc = copy.deepcopy(eos_conc)

        for j in range(self.n_steps + 1):
            inter_res = self.operator(inter_conc, power)

            if j <= 1:
                res_bar = copy.deepcopy(inter_res)
            else:
                rates = 1 / j * inter_res.rates + (1 - 1 / j) * res_bar.rates
                k = 1 / j * inter_res.k + (1 - 1 / j) * res_bar.k
                res_bar = OperatorResult(k, rates)

            inputs = list(zip(self._prev_rates, bos_rates, res_bar.rates,
                              repeat(prev_dt), repeat(dt)))
            time1, inter_conc = timed_deplete(
                self.chain, bos_conc, inputs, dt, matrix_func=_leqi_f3)
            time2, inter_conc = timed_deplete(
                self.chain, inter_conc, inputs, dt, matrix_func=_leqi_f4)
            proc_time += time1 + time2

        return proc_time, [eos_conc, inter_conc], [res_bar]
