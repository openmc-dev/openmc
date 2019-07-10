from copy import deepcopy
from abc import ABC, abstractmethod
from collections.abc import Iterable

from openmc.deplete import Results


class Integrator(ABC):
    """Abstract class for solving the time-integration for depletion

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
        self.operator = operator
        self.chain = operator.chain
        if not isinstance(timesteps, Iterable):
            self.timesteps = [timesteps]
        else:
            self.timesteps = timesteps
        if power is None:
            if power_density is None:
                raise ValueError("Either power or power density must be set")
            if not isinstance(power_density, Iterable):
                power = power_density * operator.heavy_metal
            else:
                power = [p * operator.heavy_metal for p in power_density]

        if not isinstance(power, Iterable):
            # TODO Maybe use itertools.zip_longest?
            # Ensure that power is single value if that is the case
            power = [power] * len(self.timesteps)

        self.power = power

    @abstractmethod
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
            Concentrations at each of the intermediate points with
            the final concentration as the last element
        op_results : list of openmc.deplete.OperatorResult
            Eigenvalue and reaction rates from intermediate transport
            simulations
        """

    def __iter__(self):
        for dt, p in zip(self.timesteps, self.power):
            yield dt, p

    def __len__(self):
        return len(self.timesteps)

    def _get_bos_data(self, step_index, step_power, prev_conc):
        if step_index > 0 or self.operator.prev_res is None:
            x = deepcopy(prev_conc)
            res = self.operator(x, step_power)
        else:
            # Get previous concentration
            x = self.operator.prev_res[-1].data[0]

            # Get reaction rates and keff
            res = self.operator.prev_res[-1]
            res.rates = res.rates[0]
            res.k = res.k[0]

            # Scale rates by ratio of powers
            res.rates *= step_power / res.power[0]
        return x, res

    def _get_start_data(self):
        if self.operator.prev_res is None:
            return 0.0, 0
        return self.operator.prev_res[-1].time[-1], len(self.operator.prev_res)

    def integrate_all(self):
        """Perform the entire depletion process across all steps"""
        with self.operator as conc:
            t, i_start = self._get_start_data()

            for i, (dt, p) in enumerate(self):
                conc, res = self._get_bos_data(i, p, conc)
                proc_time, conc_list, res_list = self(conc, res.rates, dt, p)

                # Insert BOS concentration, transport results
                conc_list.insert(0, conc)
                res_list.insert(0, res)

                # Remove actual EOS concentration for next step
                conc = conc_list.pop()

                self._save_results(
                    conc_list, res_list, [t, t + dt], p, i_start + i,
                    proc_time)

                t += dt

            # Final simulation
            res_list = [self.operator(conc, p)]
            self._save_results(
                [conc], res_list, [t, t], p, i_start + len(self))

    def _save_results(self, conc_list, results_list, time_list, power,
                      index, proc_time=None):
        """Save the results at the end of of one step

        Abstracted to support the predictor's unique save location
        """
        Results.save(
            self.operator, conc_list, results_list, time_list,
            power, index, proc_time)

    @classmethod
    def integrate(cls, operator, timesteps, power=None, power_density=None):
        """High-level interface for depleting with this integrator"""
        return cls(operator, timesteps, power, power_density).integrate_all()
