from copy import deepcopy
from abc import ABC, abstractmethod
from collections.abc import Iterable
from numbers import Integral

from uncertainties import ufloat

from openmc.checkvalue import check_type, check_greater_than
from openmc.capi import statepoint_write
from openmc.deplete import Results, OperatorResult


class Integrator(ABC):
    """Abstract class for solving the time-integration for depletion

    Parameters
    ----------
    operator : openmc.deplete.TransportOperator
        Operator to perform transport simulations
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
    """

    def __init__(self, operator, timesteps, power=None, power_density=None):
        # Check number of stages previously used
        if operator.prev_res is not None:
            res = operator.prev_res[-1]
            if res.data.shape[0] != self._num_stages:
                raise ValueError(
                    "{} incompatible with previous restart calculation. "
                    "Previous scheme used {} intermediate solutions, while "
                    "this uses {}".format(
                        self.__class__.__name__, res.data.shape[0],
                        self._num_stages))
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
            # Ensure that power is single value if that is the case
            power = [power] * len(self.timesteps)
        elif len(power) != len(self.timesteps):
            raise ValueError(
                "Number of time steps != number of powers. {} vs {}".format(
                    len(self.timesteps), len(power)))

        self.power = power

    @abstractmethod
    def __call__(self, conc, rates, dt, power, i):
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
            simulations
        """

    @property
    @abstractmethod
    def _num_stages(self):
        """Number of intermediate transport solutions

        Needed to ensure schemes are consistent with restarts
        """

    def __iter__(self):
        """Return pairs of time steps in [s] and powers in [W]"""
        return zip(self.timesteps, self.power)

    def __len__(self):
        """Return integer number of depletion intervals"""
        return len(self.timesteps)

    def _get_bos_data_from_operator(self, step_index, step_power, bos_conc):
        """Get beginning of step concentrations, reaction rates from Operator"""
        x = deepcopy(bos_conc)
        res = self.operator(x, step_power)
        self.operator.write_bos_data(step_index + self._i_res)
        return x, res

    def _get_bos_data_from_restart(self, step_index, step_power, bos_conc):
        """Get beginning of step concentrations, reaction rates from restart"""
        res = self.operator.prev_res[-1]
        # Depletion methods expect list of arrays
        bos_conc = list(res.data[0])
        rates = res.rates[0]
        k = ufloat(res.k[0, 0], res.k[0, 1])

        # Scale rates by ratio of powers
        rates *= step_power / res.power[0]
        return bos_conc, OperatorResult(k, rates)

    def _get_start_data(self):
        if self.operator.prev_res is None:
            return 0.0, 0
        return (self.operator.prev_res[-1].time[-1],
                len(self.operator.prev_res) - 1)

    def integrate(self):
        """Perform the entire depletion process across all steps"""
        with self.operator as conc:
            t, self._i_res = self._get_start_data()

            for i, (dt, p) in enumerate(self):
                if i > 0 or self.operator.prev_res is None:
                    conc, res = self._get_bos_data_from_operator(i, p, conc)
                else:
                    conc, res = self._get_bos_data_from_restart(i, p, conc)
                proc_time, conc_list, res_list = self(conc, res.rates, dt, p, i)

                # Insert BOS concentration, transport results
                conc_list.insert(0, conc)
                res_list.insert(0, res)

                # Remove actual EOS concentration for next step
                conc = conc_list.pop()

                Results.save(self.operator, conc_list, res_list, [t, t + dt],
                             p, self._i_res + i, proc_time)

                t += dt

            # Final simulation
            res_list = [self.operator(conc, p)]
            Results.save(self.operator, [conc], res_list, [t, t],
                         p, self._i_res + len(self), proc_time)
            self.operator.write_bos_data(len(self) + self._i_res)


class SIIntegrator(Integrator):
    """Abstract class for the Stochastic Implicit Euler integrators

    Does not provide a ``__call__`` method, but scales and resets
    the number of particles used in initial transport calculation

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
    def __init__(self, operator, timesteps, power=None, power_density=None,
                 n_steps=10):
        check_type("n_steps", n_steps, Integral)
        check_greater_than("n_steps", n_steps, 0)
        super().__init__(operator, timesteps, power, power_density)
        self.n_steps = n_steps

    def _get_bos_data_from_operator(self, step_index, step_power, bos_conc):
        reset_particles = False
        if step_index == 0 and hasattr(self.operator, "settings"):
            reset_particles = True
            self.operator.settings.particles *= self.n_stages
        inherited = super()._get_bos_data_from_operator(
            step_index, step_power, bos_conc)
        if reset_particles:
            self.operator.settings.particles //= self.n_stages
        return inherited

    def integrate(self):
        """Perform the entire depletion process across all steps"""
        with self.operator as conc:
            t, self._i_res = self._get_start_data()

            for i, (dt, p) in enumerate(self):
                if i == 0:
                    if self.operator.prev_res is None:
                        conc, res = self._get_bos_data_from_operator(i, p, conc)
                    else:
                        conc, res = self._get_bos_data_from_restart(i, p, conc)
                else:
                    # Pull rates, k from previous iteration w/o
                    # re-running transport
                    res = res_list[-1]  # defined in previous i iteration

                proc_time, conc_list, res_list = self(conc, res.rates, dt, p, i)

                # Insert BOS concentration, transport results
                conc_list.insert(0, conc)
                res_list.insert(0, res)

                # Remove actual EOS concentration for next step
                conc = conc_list.pop()

                Results.save(self.operator, conc_list, res_list, [t, t + dt],
                             p, self._i_res + i, proc_time)

                t += dt

            # No final simulation for SIE, use last iteration results
            Results.save(self.operator, [conc], [res_list[-1]], [t, t],
                         p, self._i_res + len(self), proc_time)
            self.operator.write_bos_data(self._i_res + len(self))
