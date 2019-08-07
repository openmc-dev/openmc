"""function module.

This module contains the Operator class, which is then passed to an integrator
to run a full depletion simulation.
"""

from collections import namedtuple
from collections.abc import Iterable
import os
from pathlib import Path
from abc import ABC, abstractmethod
from copy import deepcopy
from warnings import warn
from numbers import Real, Integral

from numpy import nonzero, empty
from uncertainties import ufloat

from openmc.data import DataLibrary, JOULE_PER_EV
from openmc.checkvalue import check_type, check_greater_than
from .results import Results
from .chain import Chain
from .results_list import ResultsList

OperatorResult = namedtuple('OperatorResult', ['k', 'rates'])
OperatorResult.__doc__ = """\
Result of applying transport operator

Parameters
----------
k : uncertainties.ufloat
    Resulting eigenvalue and standard deviation
rates : openmc.deplete.ReactionRates
    Resulting reaction rates

"""
try:
    OperatorResult.k.__doc__ = None
    OperatorResult.rates.__doc__ = None
except AttributeError:
    # Can't set __doc__ on properties on Python 3.4
    pass


class TransportOperator(ABC):
    """Abstract class defining a transport operator

    Each depletion integrator is written to work with a generic transport
    operator that takes a vector of material compositions and returns an
    eigenvalue and reaction rates. This abstract class sets the requirements
    for such a transport operator. Users should instantiate
    :class:`openmc.deplete.Operator` rather than this class.

    Parameters
    ----------
    chain_file : str, optional
        Path to the depletion chain XML file.  Defaults to the file
        listed under ``depletion_chain`` in
        :envvar:`OPENMC_CROSS_SECTIONS` environment variable.
    fission_q : dict, optional
        Dictionary of nuclides and their fission Q values [eV]. If not given,
        values will be pulled from the ``chain_file``.
    dilute_initial : float, optional
        Initial atom density [atoms/cm^3] to add for nuclides that are zero
        in initial condition to ensure they exist in the decay chain.
        Only done for nuclides with reaction rates.
        Defaults to 1.0e3.
    prev_results : ResultsList, optional
        Results from a previous depletion calculation.

    Attributes
    ----------
    dilute_initial : float
        Initial atom density [atoms/cm^3] to add for nuclides that are zero
        in initial condition to ensure they exist in the decay chain.
        Only done for nuclides with reaction rates.
    prev_res : ResultsList or None
        Results from a previous depletion calculation. ``None`` if no
        results are to be used.
    """
    def __init__(self, chain_file=None, fission_q=None, dilute_initial=1.0e3,
                 prev_results=None):
        self.dilute_initial = dilute_initial
        self.output_dir = '.'

        # Read depletion chain
        if chain_file is None:
            chain_file = os.environ.get("OPENMC_DEPLETE_CHAIN", None)
            if chain_file is None:
                data = DataLibrary.from_xml()
                # search for depletion_chain path from end of list
                for lib in reversed(data.libraries):
                    if lib['type'] == 'depletion_chain':
                        break
                else:
                    raise IOError(
                        "No chain specified, either manually or "
                        "under depletion_chain in environment variable "
                        "OPENMC_CROSS_SECTIONS.")
                chain_file = lib['path']
            else:
                warn("Use of OPENMC_DEPLETE_CHAIN is deprecated in favor "
                     "of adding depletion_chain to OPENMC_CROSS_SECTIONS",
                     FutureWarning)
        self.chain = Chain.from_xml(chain_file, fission_q)
        if prev_results is None:
            self.prev_res = None
        else:
            check_type("previous results", prev_results, ResultsList)
            self.prev_results = prev_results

    @property
    def dilute_initial(self):
        """Initial atom density for nuclides with zero initial concentration"""
        return self._dilute_initial

    @dilute_initial.setter
    def dilute_initial(self, value):
        check_type("dilute_initial", value, Real)
        check_greater_than("dilute_initial", value, 0.0, equality=True)
        self._dilute_initial = value

    @abstractmethod
    def __call__(self, vec, power):
        """Runs a simulation.

        Parameters
        ----------
        vec : list of numpy.ndarray
            Total atoms to be used in function.
        power : float
            Power of the reactor in [W]

        Returns
        -------
        openmc.deplete.OperatorResult
            Eigenvalue and reaction rates resulting from transport operator

        """

    def __enter__(self):
        # Save current directory and move to specific output directory
        self._orig_dir = os.getcwd()
        if not self.output_dir.exists():
            self.output_dir.mkdir()  # exist_ok parameter is 3.5+

        # In Python 3.6+, chdir accepts a Path directly
        os.chdir(str(self.output_dir))

        return self.initial_condition()

    def __exit__(self, exc_type, exc_value, traceback):
        self.finalize()
        os.chdir(self._orig_dir)

    @property
    def output_dir(self):
        return self._output_dir

    @output_dir.setter
    def output_dir(self, output_dir):
        self._output_dir = Path(output_dir)

    @abstractmethod
    def initial_condition(self):
        """Performs final setup and returns initial condition.

        Returns
        -------
        list of numpy.ndarray
            Total density for initial conditions.
        """

    @abstractmethod
    def get_results_info(self):
        """Returns volume list, cell lists, and nuc lists.

        Returns
        -------
        volume : list of float
            Volumes corresponding to materials in burn_list
        nuc_list : list of str
            A list of all nuclide names. Used for sorting the simulation.
        burn_list : list of int
            A list of all cell IDs to be burned.  Used for sorting the
            simulation.
        full_burn_list : list of int
            All burnable materials in the geometry.
        """

    def finalize(self):
        pass

    @abstractmethod
    def write_bos_data(self, step):
        """Document beginning of step data for a given step

        Called at the beginning of a depletion step and at
        the final point in the simulation.

        Parameters
        ----------
        step : int
            Current depletion step including restarts
        """


class ReactionRateHelper(ABC):
    """Abstract class for generating reaction rates for operators

    Responsible for generating reaction rate tallies for burnable
    materials, given nuclides and scores from the operator.

    Reaction rates are passed back to the operator for be used in
    an :class:`openmc.deplete.OperatorResult` instance

    Parameters
    ----------
    n_nucs : int
        Number of burnable nuclides tracked by :class:`openmc.deplete.Operator`
    n_react : int
        Number of reactions tracked by :class:`openmc.deplete.Operator`

    Attributes
    ----------
    nuclides : list of str
        All nuclides with desired reaction rates.
    """

    def __init__(self, n_nucs, n_react):
        self._nuclides = None
        self._rate_tally = None
        self._results_cache = empty((n_nucs, n_react))

    @abstractmethod
    def generate_tallies(self, materials, scores):
        """Use the C API to build tallies needed for reaction rates"""

    @property
    def nuclides(self):
        """List of nuclides with requested reaction rates"""
        return self._nuclides

    @nuclides.setter
    def nuclides(self, nuclides):
        check_type("nuclides", nuclides, list, str)
        self._nuclides = nuclides
        self._rate_tally.nuclides = nuclides

    @abstractmethod
    def get_material_rates(self, mat_id, nuc_index, react_index):
        """Return 2D array of [nuclide, reaction] reaction rates

        Parameters
        ----------
        mat_id : int
            Unique ID for the requested material
        nuc_index : list of str
            Ordering of desired nuclides
        react_index : list of str
            Ordering of reactions
        """

    def divide_by_adens(self, number):
        """Normalize reaction rates by number of nuclides

        Acts on the current material examined by
        :meth:`get_material_rates`

        Parameters
        ----------
        number : iterable of float
            Number density [atoms/b-cm] of each nuclide tracked in the
            calculation.

        Returns
        -------
        results : numpy.ndarray
            Array of reactions rates of shape ``(n_nuclides, n_rxns)``
            normalized by the number of nuclides
        """

        mask = nonzero(number)
        results = self._results_cache
        for col in range(results.shape[1]):
            results[mask, col] /= number[mask]
        return results


class EnergyHelper(ABC):
    """Abstract class for obtaining energy produced

    The ultimate goal of this helper is to provide instances of
    :class:`openmc.deplete.Operator` with the total energy produced
    in a transport simulation. This information, provided with the
    power requested by the user and reaction rates from a
    :class:`ReactionRateHelper` will scale reaction rates to the
    correct values.

    Attributes
    ----------
    nuclides : list of str
        All nuclides with desired reaction rates. Ordered to be
        consistent with :class:`openmc.deplete.Operator`
    energy : float
        Total energy [J/s/source neutron] produced in a transport simulation.
        Updated in the material iteration with :meth:`update`.
    """

    def __init__(self):
        self._nuclides = None
        self._energy = 0.0

    @property
    def energy(self):
        return self._energy * JOULE_PER_EV

    def reset(self):
        """Reset energy produced prior to unpacking tallies"""
        self._energy = 0.0

    @abstractmethod
    def prepare(self, chain_nucs, rate_index, materials):
        """Perform work needed to obtain energy produced

        This method is called prior to the transport simulations
        in :meth:`openmc.deplete.Operator.initial_condition`.

        Parameters
        ----------
        chain_nucs : list of str
            All nuclides to be tracked in this problem
        rate_index : dict of str to int
            Mapping from nuclide name to index in the
            `fission_rates` for :meth:`update`.
        materials : list of str
            All materials tracked on the operator helped by this
            object. Should correspond to
            :attr:`openmc.deplete.Operator.burnable_materials`
        """

    def update(self, fission_rates, mat_index):
        """Update the energy produced

        Parameters
        ----------
        fission_rates : numpy.ndarray
            fission reaction rate for each isotope in the specified
            material. Should be ordered corresponding to initial
            ``rate_index`` used in :meth:`prepare`
        mat_index : int
            Index for the specific material in the list of all burnable
            materials.
        """

    @property
    def nuclides(self):
        """List of nuclides with requested reaction rates"""
        return self._nuclides

    @nuclides.setter
    def nuclides(self, nuclides):
        check_type("nuclides", nuclides, list, str)
        self._nuclides = nuclides


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
        """Get beginning of step concentrations, reaction rates from Operator
        """
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
