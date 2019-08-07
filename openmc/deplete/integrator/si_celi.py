"""The SI-CE/LI CFQ4 integrator."""

import copy

from .abc import SIIntegrator
from .cram import timed_deplete
from .celi import _celi_f1, _celi_f2
from ..abc import OperatorResult


class SICELIIntegrator(SIIntegrator):
    r"""Deplete using the SI-CE/LI CFQ4 algorithm.

    Implements the stochastic implicit CE/LI predictor-corrector algorithm
    using the `fourth order commutator-free integrator
    <https://doi.org/10.1137/05063042>`_.

    Detailed algorithm can be found in section 3.2 in `Colin Josey's thesis
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

    def __call__(self, bos_conc, bos_rates, dt, power, _i):
        """Perform the integration across one time step

        Parameters
        ----------
        bos_conc : numpy.ndarray
            Initial bos_concentrations for all nuclides in [atom]
        bos_rates : openmc.deplete.ReactionRates
            Reaction rates from operator
        dt : float
            Time in [s] for the entire depletion interval
        power : float
            Power of the system in [W]
        _i : int
            Current depletion step index. Unused

        Returns
        -------
        proc_time : float
            Time spent in CRAM routines for all materials in [s]
        bos_conc_list : list of numpy.ndarray
            Concentrations at each of the intermediate points with
            the final bos_concentration as the last element
        op_results : list of openmc.deplete.OperatorResult
            Eigenvalue and reaction rates from intermediate transport
            simulations
        """
        proc_time, eos_conc = timed_deplete(
            self.chain, bos_conc, bos_rates, dt)
        inter_conc = copy.deepcopy(eos_conc)

        # Begin iteration
        for j in range(self.n_steps + 1):
            inter_res = self.operator(inter_conc, power)

            if j <= 1:
                res_bar = copy.deepcopy(inter_res)
            else:
                rates = 1/j * inter_res.rates + (1 - 1 / j) * res_bar.rates
                k = 1/j * inter_res.k + (1 - 1 / j) * res_bar.k
                res_bar = OperatorResult(k, rates)

            list_rates = list(zip(bos_rates, res_bar.rates))
            time1, inter_conc = timed_deplete(
                self.chain, bos_conc, list_rates, dt, matrix_func=_celi_f1)
            time2, inter_conc = timed_deplete(
                self.chain, inter_conc, list_rates, dt, matrix_func=_celi_f2)
            proc_time += time1 + time2

        # end iteration
        return proc_time, [eos_conc, inter_conc], [res_bar]
