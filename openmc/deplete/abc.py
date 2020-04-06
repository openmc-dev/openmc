"""function module.

This module contains the Operator class, which is then passed to an integrator
to run a full depletion simulation.
"""

from collections import namedtuple
from collections import defaultdict
from collections.abc import Iterable
import os
from pathlib import Path
from abc import ABC, abstractmethod
from copy import deepcopy
from warnings import warn
from numbers import Real, Integral

from numpy import nonzero, empty, asarray
from uncertainties import ufloat

from openmc.data import DataLibrary, JOULE_PER_EV
from openmc.lib import MaterialFilter, Tally
from openmc.checkvalue import check_type, check_greater_than
from .results import Results
from .chain import Chain
from .results_list import ResultsList


__all__ = [
    "OperatorResult", "TransportOperator", "ReactionRateHelper",
    "EnergyHelper", "FissionYieldHelper", "TalliedFissionYieldHelper",
    "Integrator", "SIIntegrator", "DepSystemSolver"]


_SECONDS_PER_MINUTE = 60
_SECONDS_PER_HOUR = 60*60
_SECONDS_PER_DAY = 24*60*60

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
        volume : dict of str to float
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


class FissionYieldHelper(ABC):
    """Abstract class for processing energy dependent fission yields

    Parameters
    ----------
    chain_nuclides : iterable of openmc.deplete.Nuclide
        Nuclides tracked in the depletion chain. All nuclides are
        not required to have fission yield data.

    Attributes
    ----------
    constant_yields : collections.defaultdict
        Fission yields for all nuclides that only have one set of
        fission yield data. Dictionary of form ``{str: {str: float}}``
        representing yields for ``{parent: {product: yield}}``. Default
        return object is an empty dictionary

    """

    def __init__(self, chain_nuclides):
        self._chain_nuclides = {}
        self._constant_yields = defaultdict(dict)

        # Get all nuclides with fission yield data
        for nuc in chain_nuclides:
            if nuc.yield_data is None:
                continue
            if len(nuc.yield_data) == 1:
                self._constant_yields[nuc.name] = (
                    nuc.yield_data[nuc.yield_energies[0]])
            elif len(nuc.yield_data) > 1:
                self._chain_nuclides[nuc.name] = nuc
        self._chain_set = set(self._chain_nuclides) | set(self._constant_yields)

    @property
    def constant_yields(self):
        return deepcopy(self._constant_yields)

    @abstractmethod
    def weighted_yields(self, local_mat_index):
        """Return fission yields for a specific material

        Parameters
        ----------
        local_mat_index : int
            Index for the material with requested fission yields.
            Should correspond to the material represented in
            ``mat_indexes[local_mat_index]`` during
            :meth:`generate_tallies`.

        Returns
        -------
        library : collections.abc.Mapping
            Dictionary-like object mapping ``{str: {str: float}``.
            This reflects fission yields for ``{parent: {product: fyield}}``.
        """

    @staticmethod
    def unpack():
        """Unpack tally data prior to compute fission yields.

        Called after a :meth:`openmc.deplete.Operator.__call__`
        routine during the normalization of reaction rates.

        Not necessary for all subclasses to implement, unless tallies
        are used.
        """

    @staticmethod
    def generate_tallies(materials, mat_indexes):
        """Construct tallies necessary for computing fission yields

        Called during the operator set up phase prior to depleting.
        Not necessary for subclasses to implement

        Parameters
        ----------
        materials : iterable of C-API materials
            Materials to be used in :class:`openmc.lib.MaterialFilter`
        mat_indexes : iterable of int
            Indices of tallied materials that will have their fission
            yields computed by this helper. Necessary as the
            :class:`openmc.deplete.Operator` that uses this helper
            may only burn a subset of all materials when running
            in parallel mode.
        """

    def update_tally_nuclides(self, nuclides):
        """Return nuclides with non-zero densities and yield data

        Parameters
        ----------
        nuclides : iterable of str
            Nuclides with non-zero densities from the
            :class:`openmc.deplete.Operator`

        Returns
        -------
        nuclides : list of str
            Union of nuclides that the :class:`openmc.deplete.Operator`
            says have non-zero densities at this stage and those that
            have yield data. Sorted by nuclide name

        """
        return sorted(self._chain_set & set(nuclides))

    @classmethod
    def from_operator(cls, operator, **kwargs):
        """Create a new instance by pulling data from the operator

        All keyword arguments should be identical to their counterpart
        in the main ``__init__`` method

        Parameters
        ----------
        operator : openmc.deplete.TransportOperator
            Operator with a depletion chain
        kwargs: optional
            Additional keyword arguments to be used in constuction
        """
        return cls(operator.chain.nuclides, **kwargs)


class TalliedFissionYieldHelper(FissionYieldHelper):
    """Abstract class for computing fission yields with tallies

    Generates a basic fission rate tally in all burnable materials with
    :meth:`generate_tallies`, and set nuclides to be tallied with
    :meth:`update_tally_nuclides`. Subclasses will need to implement
    :meth:`unpack` and :meth:`weighted_yields`.

    Parameters
    ----------
    chain_nuclides : iterable of openmc.deplete.Nuclide
        Nuclides tracked in the depletion chain. Not necessary
        that all have yield data.

    Attributes
    ----------
    constant_yields : dict of str to :class:`openmc.deplete.FissionYield`
        Fission yields for all nuclides that only have one set of
        fission yield data. Can be accessed as ``{parent: {product: yield}}``
    results : None or numpy.ndarray
        Tally results shaped in a manner useful to this helper.
    """

    _upper_energy = 20.0e6  # upper energy for tallies

    def __init__(self, chain_nuclides):
        super().__init__(chain_nuclides)
        self._local_indexes = None
        self._fission_rate_tally = None
        self._tally_nucs = []
        self.results = None

    def generate_tallies(self, materials, mat_indexes):
        """Construct the fission rate tally

        Parameters
        ----------
        materials : iterable of :class:`openmc.lib.Material`
            Materials to be used in :class:`openmc.lib.MaterialFilter`
        mat_indexes : iterable of int
            Indices of tallied materials that will have their fission
            yields computed by this helper. Necessary as the
            :class:`openmc.deplete.Operator` that uses this helper
            may only burn a subset of all materials when running
            in parallel mode.
        """
        self._local_indexes = asarray(mat_indexes)

        # Tally group-wise fission reaction rates
        self._fission_rate_tally = Tally()
        self._fission_rate_tally.writable = False
        self._fission_rate_tally.scores = ['fission']

        self._fission_rate_tally.filters = [MaterialFilter(materials)]

    def update_tally_nuclides(self, nuclides):
        """Tally nuclides with non-zero density and multiple yields

        Must be run after :meth:`generate_tallies`.

        Parameters
        ----------
        nuclides : iterable of str
            Potential nuclides to be tallied, such as those with
            non-zero density at this stage.

        Returns
        -------
        nuclides : list of str
            Union of input nuclides and those that have multiple sets
            of yield data.  Sorted by nuclide name

        Raises
        ------
        AttributeError
            If tallies not generated
        """
        assert self._fission_rate_tally is not None, (
                "Run generate_tallies first")
        overlap = set(self._chain_nuclides).intersection(set(nuclides))
        nuclides = sorted(overlap)
        self._tally_nucs = [self._chain_nuclides[n] for n in nuclides]
        self._fission_rate_tally.nuclides = nuclides
        return nuclides

    @abstractmethod
    def unpack(self):
        """Unpack tallies after a transport run.

        Abstract because each subclass will need to arrange its
        tally data.
        """


class Integrator(ABC):
    """Abstract class for solving the time-integration for depletion

    Parameters
    ----------
    operator : openmc.deplete.TransportOperator
        Operator to perform transport simulations
    timesteps : iterable of float or iterable of tuple
        Array of timesteps. Note that values are not cumulative. The units are
        specified by the `timestep_units` argument when `timesteps` is an
        iterable of float. Alternatively, units can be specified for each step
        by passing an iterable of (value, unit) tuples.
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
    timestep_units : {'s', 'min', 'h', 'd', 'MWd/kg'}
        Units for values specified in the `timesteps` argument. 's' means
        seconds, 'min' means minutes, 'h' means hours, and 'MWd/kg' indicates
        that the values are given in burnup (MW-d of energy deposited per
        kilogram of initial heavy metal).

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

    def __init__(self, operator, timesteps, power=None, power_density=None,
                 timestep_units='s'):
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

        # Determine power and normalize units to W
        if power is None:
            if power_density is None:
                raise ValueError("Either power or power density must be set")
            if not isinstance(power_density, Iterable):
                power = power_density * operator.heavy_metal
            else:
                power = [p*operator.heavy_metal for p in power_density]
        if not isinstance(power, Iterable):
            # Ensure that power is single value if that is the case
            power = [power] * len(timesteps)

        if len(power) != len(timesteps):
            raise ValueError(
                "Number of time steps ({}) != number of powers ({})".format(
                    len(timesteps), len(power)))

        # Get list of times / units
        if isinstance(timesteps[0], Iterable):
            times, units = zip(*timesteps)
        else:
            times = timesteps
            units = [timestep_units] * len(timesteps)

        # Determine number of seconds for each timestep
        seconds = []
        for time, unit, watts in zip(times, units, power):
            # Make sure values passed make sense
            check_type('timestep', time, Real)
            check_greater_than('timestep', time, 0.0, False)
            check_type('timestep units', unit, str)
            check_type('power', watts, Real)
            check_greater_than('power', watts, 0.0, True)

            if unit in ('s', 'sec'):
                seconds.append(time)
            elif unit in ('min', 'minute'):
                seconds.append(time*_SECONDS_PER_MINUTE)
            elif unit in ('h', 'hr', 'hour'):
                seconds.append(time*_SECONDS_PER_HOUR)
            elif unit in ('d', 'day'):
                seconds.append(time*_SECONDS_PER_DAY)
            elif unit.lower() == 'mwd/kg':
                watt_days_per_kg = 1e6*time
                kilograms = 1e-3*operator.heavy_metal
                days = watt_days_per_kg * kilograms / watts
                seconds.append(days*_SECONDS_PER_DAY)
            else:
                raise ValueError("Invalid timestep unit '{}'".format(unit))

        self.timesteps = asarray(seconds)
        self.power = asarray(power)

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
    timesteps : iterable of float or iterable of tuple
        Array of timesteps. Note that values are not cumulative. The units are
        specified by the `timestep_units` argument when `timesteps` is an
        iterable of float. Alternatively, units can be specified for each step
        by passing an iterable of (value, unit) tuples.
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
    timestep_units : {'s', 'min', 'h', 'd', 'MWd/kg'}
        Units for values specified in the `timesteps` argument. 's' means
        seconds, 'min' means minutes, 'h' means hours, and 'MWd/kg' indicates
        that the values are given in burnup (MW-d of energy deposited per
        kilogram of initial heavy metal).
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
                 timestep_units='s', n_steps=10):
        check_type("n_steps", n_steps, Integral)
        check_greater_than("n_steps", n_steps, 0)
        super().__init__(operator, timesteps, power, power_density, timestep_units)
        self.n_steps = n_steps

    def _get_bos_data_from_operator(self, step_index, step_power, bos_conc):
        reset_particles = False
        if step_index == 0 and hasattr(self.operator, "settings"):
            reset_particles = True
            self.operator.settings.particles *= self.n_steps
        inherited = super()._get_bos_data_from_operator(
            step_index, step_power, bos_conc)
        if reset_particles:
            self.operator.settings.particles //= self.n_steps
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


class DepSystemSolver(ABC):
    r"""Abstract class for solving depletion equations

    Responsible for solving

    .. math::

        \frac{\partial \vec{N}}{\partial t} = \bar{A}\vec{N}(t),

    for :math:`0< t\leq t +\Delta t`, given :math:`\vec{N}(0) = \vec{N}_0`

    """

    @abstractmethod
    def __call__(self, A, n0, dt):
        """Solve the linear system of equations for depletion

        Parameters
        ----------
        A : scipy.sparse.csr_matrix
            Sparse transmutation matrix ``A[j, i]`` desribing rates at
            which isotope ``i`` transmutes to isotope ``j``
        n0 : numpy.ndarray
            Initial compositions, typically given in number of atoms in some
            material or an atom density
        dt : float
            Time [s] of the specific interval to be solved

        Returns
        -------
        numpy.ndarray
            Final compositions after ``dt``. Should be of identical shape
            to ``n0``.

        """
