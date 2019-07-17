"""First-order predictor algorithm."""

from copy import deepcopy

from .abc import Integrator
from .cram import timed_deplete


class PredictorIntegrator(Integrator):
    r"""Deplete using a first-order predictor algorithm.

    Implements the first-order predictor algorithm. This algorithm is
    mathematically defined as:

    .. math::
        \begin{aligned}
        y' &= A(y, t) y(t) \\
        A_p &= A(y_n, t_n) \\
        y_{n+1} &= \text{expm}(A_p h) y_n
        \end{aligned}

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
    """

    def __init__(self, operator, timesteps, power=None, power_density=None):
        """
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
        """
        super().__init__(operator, timesteps, power, power_density)
        self._dep_proc_time = None

    def __call__(self, conc, rates, dt, power):
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
            Power of the system [W]

        Returns
        -------
        proc_time : float
            Time spent in CRAM routines for all materials
        conc_list : list of numpy.ndarray
            Concentrations at end of interval
        op_results : empty list
            Kept for consistency with API. No intermediate calls to
            operator with predictor

        """
        proc_time, conc_end = timed_deplete(self.chain, conc, rates, dt)
        return proc_time, conc_end, []

    def _get_start_data(self):
        if self.operator.prev_res is None:
            return 0.0, 0
        return (
            self.operator.prev_res[-1].time[-1],
            len(self.operator.prev_res) - 1)

    def _get_bos_data(self, step_index, step_power, prev_conc):
        if step_index > 0 or self.operator.prev_res is None:
            conc = deepcopy(prev_conc)
            res = self.operator(conc, step_power)
            if step_index == len(self) - 1:
                tvec = [self.timesteps[step_index]] * 2
            else:
                tvec = self.timesteps[step_index:step_index + 2]

            self._save_results(
                [conc], [res], tvec, step_power,
                self._istart + step_index, self._dep_proc_time)
        else:
            # Get previous concentration
            conc = self.operator.prev_res[-1].data[0]

            # Get reaction rates and keff
            res = self.operator.prev_res[-1]
            res.rates = res.rates[0]
            res.k = res.k[0]

            # Scale rates by ratio of powers
            res.rates *= step_power / res.power[0]

        return conc, res

    def integrate(self):
        """Perform the entire depletion process across all steps"""

        with self.operator as conc:
            t, self._istart = self._get_start_data()

            for i, (dt, p) in enumerate(self):
                conc, res = self._get_bos_data(i, p, conc)
                # __call__ returns empty list since there aren't
                # intermediate transport solutions
                self._dep_proc_time, conc, _res_list = self(
                    conc, res.rates, dt, p)

                t += dt

            # Final simulation
            res_list = [self.operator(conc, p)]
            self._save_results(
                [conc], res_list, [t, t], p,
                self._istart + len(self), self._dep_proc_time)
