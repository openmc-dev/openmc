"""abc module.

This module contains Abstract Base Classes for implementing operator, integrator, depletion system solver, and operator helper classes
"""

from __future__ import annotations
from abc import ABC, abstractmethod
from collections import namedtuple, defaultdict
from collections.abc import Iterable, Callable
from copy import deepcopy
from inspect import signature
from numbers import Real, Integral
from pathlib import Path
import time
from typing import Optional, Union, Sequence
from warnings import warn

import numpy as np
from uncertainties import ufloat

from openmc.checkvalue import check_type, check_greater_than, PathLike, check_value
from openmc.mpi import comm
from openmc.utility_funcs import change_directory
from openmc import Material
from .stepresult import StepResult
from .chain import Chain
from .results import Results
from .pool import deplete
from .reaction_rates import ReactionRates
from .transfer_rates import TransferRates
from openmc import Material, Cell
from .batchwise import (BatchwiseCellGeometrical, BatchwiseCellTemperature,
    BatchwiseMaterialRefuel, BatchwiseMaterialDilute, BatchwiseMaterialAdd,
    BatchwiseSchemeStd, BatchwiseSchemeRefuel, BatchwiseSchemeFlex)

__all__ = [
    "OperatorResult", "TransportOperator",
    "ReactionRateHelper", "NormalizationHelper", "FissionYieldHelper",
    "Integrator", "SIIntegrator", "DepSystemSolver", "add_params"]


_SECONDS_PER_MINUTE = 60
_SECONDS_PER_HOUR = 60*60
_SECONDS_PER_DAY = 24*60*60
_SECONDS_PER_JULIAN_YEAR = 365.25*24*60*60

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
    :class:`openmc.deplete.CoupledOperator` or
    :class:`openmc.deplete.IndependentOperator` rather than this class.

    Parameters
    ----------
    chain_file : str
        Path to the depletion chain XML file
    fission_q : dict, optional
        Dictionary of nuclides and their fission Q values [eV]. If not given,
        values will be pulled from the ``chain_file``.
    prev_results : Results, optional
        Results from a previous depletion calculation.

    Attributes
    ----------
    output_dir : pathlib.Path
        Path to output directory to save results.
    prev_res : Results or None
        Results from a previous depletion calculation. ``None`` if no
        results are to be used.
    chain : openmc.deplete.Chain
        The depletion chain information necessary to form matrices and tallies.

    """
    def __init__(self, chain_file, fission_q=None, prev_results=None):
        self.output_dir = '.'

        # Read depletion chain
        self.chain = Chain.from_xml(chain_file, fission_q)
        if prev_results is None:
            self.prev_res = None
        else:
            check_type("previous results", prev_results, Results)
            self.prev_res = prev_results

    @abstractmethod
    def __call__(self, vec, source_rate):
        """Runs a simulation.

        Parameters
        ----------
        vec : list of numpy.ndarray
            Total atoms to be used in function.
        source_rate : float
            Power in [W] or source rate in [neutron/sec]

        Returns
        -------
        openmc.deplete.OperatorResult
            Eigenvalue and reaction rates resulting from transport operator

        """

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
    def write_bos_data(self, step: int):
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

    Responsible for generating reaction rate tallies for burnable materials,
    given nuclides and scores from the operator.

    Reaction rates are passed back to the operator to be used by an
    :class:`openmc.deplete.OperatorResult` instance.

    Parameters
    ----------
    n_nucs : int
        Number of burnable nuclides tracked by
        :class:`openmc.deplete.abc.TransportOperator`
    n_react : int
        Number of reactions tracked by
        :class:`openmc.deplete.abc.TransportOperator`

    Attributes
    ----------
    nuclides : list of str
        All nuclides with desired reaction rates.
    """

    def __init__(self, n_nucs, n_react):
        self._nuclides = None
        self._results_cache = np.empty((n_nucs, n_react))

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

    @abstractmethod
    def get_material_rates(
        self,
        mat_id: int,
        nuc_index: Sequence[str],
        react_index: Sequence[str]
    ):
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

    def divide_by_atoms(self, number: Sequence[float]):
        """Normalize reaction rates by number of atoms

        Acts on the current material examined by :meth:`get_material_rates`

        Parameters
        ----------
        number : iterable of float
            Number of each nuclide in [atom] tracked in the calculation.

        Returns
        -------
        results : numpy.ndarray
            Array of reactions rates of shape ``(n_nuclides, n_rxns)``
            normalized by the number of nuclides
        """

        mask = np.nonzero(number)
        results = self._results_cache
        for col in range(results.shape[1]):
            results[mask, col] /= number[mask]
        return results


class NormalizationHelper(ABC):
    """Abstract class for obtaining normalization factor on tallies

    This helper class determines how reaction rates calculated by an instance of
    :class:`openmc.deplete.abc.TransportOperator` should be normalized for the
    purpose of constructing a burnup matrix. Based on the method chosen, the
    power or source rate provided by the user, and reaction rates from a
    :class:`ReactionRateHelper`, this class will scale reaction rates to the
    correct values.

    Attributes
    ----------
    nuclides : list of str
        All nuclides with desired reaction rates. Ordered to be
        consistent with :class:`openmc.deplete.abc.TransportOperator`

    """

    def __init__(self):
        self._nuclides = None

    def reset(self):
        """Reset state for normalization"""

    @abstractmethod
    def prepare(self, chain_nucs: Sequence[str], rate_index: dict):
        """Perform work needed to obtain energy produced

        This method is called prior to calculating the reaction rates
        in :meth:`openmc.deplete.abc.TransportOperator.initial_condition`. Only
        used for energy-based normalization.

        Parameters
        ----------
        chain_nucs : list of str
            All nuclides to be tracked in this problem
        rate_index : dict of str to int
            Mapping from nuclide name to index in the
            `fission_rates` for :meth:`update`.
        """

    def update(self, fission_rates):
        """Update the normalization based on fission rates (only used for
        energy-based normalization)

        Parameters
        ----------
        fission_rates : numpy.ndarray
            fission reaction rate for each isotope in the specified
            material. Should be ordered corresponding to initial
            ``rate_index`` used in :meth:`prepare`
        """

    @property
    def nuclides(self):
        """List of nuclides with requested reaction rates"""
        return self._nuclides

    @nuclides.setter
    def nuclides(self, nuclides):
        check_type("nuclides", nuclides, list, str)
        self._nuclides = nuclides

    @abstractmethod
    def factor(self, source_rate: float):
        """Return normalization factor

        Parameters
        ----------
        source_rate : float
            Power in [W] or source rate in [neutron/sec]

        Returns
        -------
        float
            Normalization factor for tallies

        """


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

        Called after a :meth:`openmc.deplete.abc.TransportOperator.__call__`
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
            :class:`openmc.deplete.CoupledOperator` that uses this helper
            may only burn a subset of all materials when running
            in parallel mode.
        """

    def update_tally_nuclides(self, nuclides: Sequence[str]) -> list:
        """Return nuclides with non-zero densities and yield data

        Parameters
        ----------
        nuclides : iterable of str
            Nuclides with non-zero densities from the
            :class:`openmc.deplete.abc.TransportOperator`

        Returns
        -------
        nuclides : list of str
            Union of nuclides that the
            :class:`openmc.deplete.abc.TransportOperator` says have non-zero
            densities at this stage and those that have yield data. Sorted by
            nuclide name

        """
        return sorted(self._chain_set & set(nuclides))

    @classmethod
    def from_operator(cls, operator, **kwargs):
        """Create a new instance by pulling data from the operator

        All keyword arguments should be identical to their counterpart
        in the main ``__init__`` method

        Parameters
        ----------
        operator : openmc.deplete.abc.TransportOperator
            Operator with a depletion chain
        kwargs: optional
            Additional keyword arguments to be used in constuction
        """
        return cls(operator.chain.nuclides, **kwargs)


def add_params(cls):
    cls.__doc__ += cls._params
    return cls


@add_params
class Integrator(ABC):
    r"""Abstract class for solving the time-integration for depletion
    """

    _params = r"""
    Parameters
    ----------
    operator : openmc.deplete.abc.TransportOperator
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
        an area in [cm^2]. Either ``power``, ``power_density``, or
        ``source_rates`` must be specified.
    power_density : float or iterable of float, optional
        Power density of the reactor in [W/gHM]. It is multiplied by
        initial heavy metal inventory to get total power if ``power``
        is not specified.
    source_rates : float or iterable of float, optional
        Source rate in [neutron/sec] or neutron flux in [neutron/s-cm^2] for
        each interval in :attr:`timesteps`

        .. versionadded:: 0.12.1
    timestep_units : {'s', 'min', 'h', 'd', 'a', 'MWd/kg'}
        Units for values specified in the `timesteps` argument. 's' means
        seconds, 'min' means minutes, 'h' means hours, 'a' means Julian years
        and 'MWd/kg' indicates that the values are given in burnup (MW-d of
        energy deposited per kilogram of initial heavy metal).
    solver : str or callable, optional
        If a string, must be the name of the solver responsible for
        solving the Bateman equations.  Current options are:

            * ``cram16`` - 16th order IPF CRAM
            * ``cram48`` - 48th order IPF CRAM [default]

        If a function or other callable, must adhere to the requirements in
        :attr:`solver`.

        .. versionadded:: 0.12
    Attributes
    ----------
    operator : openmc.deplete.abc.TransportOperator
        Operator to perform transport simulations
    chain : openmc.deplete.Chain
        Depletion chain
    timesteps : iterable of float
        Size of each depletion interval in [s]
    source_rates : iterable of float
        Source rate in [W] or [neutron/sec] for each interval in
        :attr:`timesteps`
    solver : callable
        Function that will solve the Bateman equations
        :math:`\frac{\partial}{\partial t}\vec{n} = A_i\vec{n}_i` with a step
        size :math:`t_i`. Can be configured using the ``solver`` argument.
        User-supplied functions are expected to have the following signature:
        ``solver(A, n0, t) -> n1`` where

            * ``A`` is a :class:`scipy.sparse.csc_matrix` making up the
              depletion matrix
            * ``n0`` is a 1-D :class:`numpy.ndarray` of initial compositions
              for a given material in atoms/cm3
            * ``t`` is a float of the time step size in seconds, and
            * ``n1`` is a :class:`numpy.ndarray` of compositions at the
              next time step. Expected to be of the same shape as ``n0``

    transfer_rates : openmc.deplete.TransferRates
        Instance of TransferRates class to perform continuous transfer during depletion
    batchwise : openmc.deplete.Batchwise
        Instance of Batchwise class to perform batch-wise scheme during
        transport-depletion simulation.

        .. versionadded:: 0.14.0

    """

    def __init__(
            self,
            operator: TransportOperator,
            timesteps: Sequence[float],
            power: Optional[Union[float, Sequence[float]]] = None,
            power_density: Optional[Union[float, Sequence[float]]] = None,
            source_rates: Optional[Sequence[float]] = None,
            timestep_units: str = 's',
            solver: str = "cram48"
        ):
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

        # Determine source rate and normalize units to W in using power
        if power is not None:
            source_rates = power
        elif power_density is not None:
            if not isinstance(power_density, Iterable):
                source_rates = power_density * operator.heavy_metal
            else:
                source_rates = [p*operator.heavy_metal for p in power_density]
        elif source_rates is None:
            raise ValueError("Either power, power_density, or source_rates must be set")

        if not isinstance(source_rates, Iterable):
            # Ensure that rate is single value if that is the case
            source_rates = [source_rates] * len(timesteps)

        if len(source_rates) != len(timesteps):
            raise ValueError(
                "Number of time steps ({}) != number of powers ({})".format(
                    len(timesteps), len(source_rates)))

        # Get list of times / units
        if isinstance(timesteps[0], Iterable):
            times, units = zip(*timesteps)
        else:
            times = timesteps
            units = [timestep_units] * len(timesteps)

        # Determine number of seconds for each timestep
        seconds = []
        for timestep, unit, rate in zip(times, units, source_rates):
            # Make sure values passed make sense
            check_type('timestep', timestep, Real)
            check_greater_than('timestep', timestep, 0.0, False)
            check_type('timestep units', unit, str)
            check_type('source rate', rate, Real)
            check_greater_than('source rate', rate, 0.0, True)

            if unit in ('s', 'sec'):
                seconds.append(timestep)
            elif unit in ('min', 'minute'):
                seconds.append(timestep*_SECONDS_PER_MINUTE)
            elif unit in ('h', 'hr', 'hour'):
                seconds.append(timestep*_SECONDS_PER_HOUR)
            elif unit in ('d', 'day'):
                seconds.append(timestep*_SECONDS_PER_DAY)
            elif unit in ('a', 'year'):
                seconds.append(timestep*_SECONDS_PER_JULIAN_YEAR)
            elif unit.lower() == 'mwd/kg':
                watt_days_per_kg = 1e6*timestep
                kilograms = 1e-3*operator.heavy_metal
                if rate == 0.0:
                    raise ValueError("Cannot specify a timestep in [MWd/kg] when"
                                     " the power is zero.")
                days = watt_days_per_kg * kilograms / rate
                seconds.append(days*_SECONDS_PER_DAY)
            else:
                raise ValueError(f"Invalid timestep unit '{unit}'")

        self.timesteps = np.asarray(seconds)
        self.source_rates = np.asarray(source_rates)

        self.transfer_rates = None
        self.batchwise = None

        if isinstance(solver, str):
            # Delay importing of cram module, which requires this file
            if solver == "cram48":
                from .cram import CRAM48
                self._solver = CRAM48
            elif solver == "cram16":
                from .cram import CRAM16
                self._solver = CRAM16
            else:
                raise ValueError(
                    f"Solver {solver} not understood. Expected 'cram48' or 'cram16'")
        else:
            self.solver = solver

    @property
    def solver(self):
        return self._solver

    @solver.setter
    def solver(self, func):
        if not isinstance(func, Callable):
            raise TypeError(
                f"Solver must be callable, not {type(func)}")
        try:
            sig = signature(func)
        except ValueError:
            # Guard against callables that aren't introspectable, e.g.
            # fortran functions wrapped by F2PY
            warn(f"Could not determine arguments to {func}. Proceeding anyways")
            self._solver = func
            return

        # Inspect arguments
        if len(sig.parameters) != 3:
            raise ValueError("Function {} does not support three arguments: "
                             "{!s}".format(func, sig))

        for ix, param in enumerate(sig.parameters.values()):
            if param.kind in {param.KEYWORD_ONLY, param.VAR_KEYWORD}:
                raise ValueError(
                    f"Keyword arguments like {ix} at position {param} are not allowed")

        self._solver = func

    def _timed_deplete(self, n, rates, dt, matrix_func=None):
        start = time.time()
        results = deplete(
            self._solver, self.chain, n, rates, dt, matrix_func,
            self.transfer_rates)
        return time.time() - start, results

    @abstractmethod
    def __call__(
        self,
        n: Sequence[np.ndarray],
        rates: ReactionRates,
        dt: float,
        source_rate: float,
        i: int
    ):
        """Perform the integration across one time step

        Parameters
        ----------
        n :  list of numpy.ndarray
            List of atom number arrays for each material. Each array in the list
            contains the number of [atom] of each nuclide.
        rates : openmc.deplete.ReactionRates
            Reaction rates from operator
        dt : float
            Time in [s] for the entire depletion interval
        source_rate : float
            Power in [W] or source rate in [neutron/sec]
        i : int
            Current depletion step index

        Returns
        -------
        proc_time : float
            Time spent in CRAM routines for all materials in [s]
        n_list : list of list of numpy.ndarray
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
        """Return pair of time step in [s] and source rate in [W] or [neutron/sec]"""
        return zip(self.timesteps, self.source_rates)

    def __len__(self):
        """Return integer number of depletion intervals"""
        return len(self.timesteps)

    def _get_bos_data_from_operator(self, step_index, source_rate, bos_conc):
        """Get beginning of step concentrations, reaction rates from Operator
        """
        x = deepcopy(bos_conc)
        res = self.operator(x, source_rate)
        self.operator.write_bos_data(step_index + self._i_res)
        return x, res

    def _get_bos_data_from_restart(self, source_rate, bos_conc):
        """Get beginning of step concentrations, reaction rates from restart"""
        res = self.operator.prev_res[-1]
        # Depletion methods expect list of arrays
        bos_conc = list(res.data[0])
        rates = res.rates[0]
        k = ufloat(res.k[0, 0], res.k[0, 1])

        # Scale reaction rates by ratio of source rates
        if res.source_rate != 0.0:
            rates *= source_rate / res.source_rate
        return bos_conc, OperatorResult(k, rates)

    def _get_start_data(self):
        if self.operator.prev_res is None:
            return 0.0, 0
        else:
            if comm.size != 1:
                return (self.operator.prev_res[-1].time[-1],
                        int(len(self.operator.prev_res)/2) - 1)
            else:
                return (self.operator.prev_res[-1].time[-1],
                        len(self.operator.prev_res) - 1)

    def _get_bos_from_batchwise(self, step_index, bos_conc):
        """Get BOS from criticality batch-wise control
        """
        x = deepcopy(bos_conc)
        # Get new vector after keff criticality control
        x, root = self.batchwise.search_for_keff(x, step_index)
        return x, root

    def integrate(
            self,
            final_step: bool = True,
            output: bool = True,
            path: PathLike = 'depletion_results.h5'
        ):
        """Perform the entire depletion process across all steps

        Parameters
        ----------
        final_step : bool, optional
            Indicate whether or not a transport solve should be run at the end
            of the last timestep.

            .. versionadded:: 0.12.1
        output : bool, optional
            Indicate whether to display information about progress

            .. versionadded:: 0.13.1
        path : PathLike
            Path to file to write. Defaults to 'depletion_results.h5'.

            .. versionadded:: 0.15.0
        """
        with change_directory(self.operator.output_dir):
            n = self.operator.initial_condition()
            t, self._i_res = self._get_start_data()

            for i, (dt, source_rate) in enumerate(self):
                if output and comm.rank == 0:
                    print(f"[openmc.deplete] t={t} s, dt={dt} s, source={source_rate}")

                # Solve transport equation (or obtain result from restart)
                if i > 0 or self.operator.prev_res is None:
                    # Update geometry/material according to batchwise definition
                    if self.batchwise:
                        if source_rate != 0.0:
                            n, root = self._get_bos_from_batchwise(i, n)
                        else:
                            # Store root at previous timestep
                            root = self.batchwise._get_cell_attrib()
                    else:
                        root = None
                    n, res = self._get_bos_data_from_operator(i, source_rate, n)
                else:
                    n, res = self._get_bos_data_from_restart(i, source_rate, n)
                    if self.batchwise:
                        root = self.operator.prev_res[-1].batchwise
                        #TODO: this is just temporary (import math)
                        import math
                        if math.isnan(root):
                            prev_res_ts = -2
                            while (math.isnan(root)):
                                root = self.operator.prev_res[prev_res_ts].batchwise
                                prev_res_ts -= 1

                        self.batchwise.update_from_restart(i, n, root)
                    else:
                        root = None

                # Solve Bateman equations over time interval
                proc_time, n_list, res_list = self(n, res.rates, dt, source_rate, i)

                # Insert BOS concentration, transport results
                n_list.insert(0, n)
                res_list.insert(0, res)

                # Remove actual EOS concentration for next step
                n = n_list.pop()
                StepResult.save(self.operator, n_list, res_list, [t, t + dt],
                            source_rate, self._i_res + i, proc_time, root, path)

                t += dt

            # Final simulation -- in the case that final_step is False, a zero
            # source rate is passed to the transport operator (which knows to
            # just return zero reaction rates without actually doing a transport
            # solve)
            if output and final_step and comm.rank == 0:
                print(f"[openmc.deplete] t={t} (final operator evaluation)")
            if self.batchwise and source_rate != 0.0:
                n, root = self._get_bos_from_batchwise(i+1, n)
            else:
                root = None
            res_list = [self.operator(n, source_rate if final_step else 0.0)]
            StepResult.save(self.operator, [n], res_list, [t, t],
                    source_rate, self._i_res + len(self), proc_time, root, path)
            self.operator.write_bos_data(len(self) + self._i_res)

        self.operator.finalize()

    def add_transfer_rate(
            self,
            material: Union[str, int, Material],
            components: Sequence[str],
            transfer_rate: float,
            transfer_rate_units: str = '1/s',
            destination_material: Optional[Union[str, int, Material]] = None
        ):
        """Add transfer rates to depletable material.

        Parameters
        ----------
        material : openmc.Material or str or int
            Depletable material
        components : list of str
            List of strings of elements and/or nuclides that share transfer rate.
            A transfer rate for a nuclide cannot be added to a material
            alongside a transfer rate for its element and vice versa.
        transfer_rate : float
            Rate at which elements are transferred. A positive or negative values
            set removal of feed rates, respectively.
        destination_material : openmc.Material or str or int, Optional
            Destination material to where nuclides get fed.
        transfer_rate_units : {'1/s', '1/min', '1/h', '1/d', '1/a'}
            Units for values specified in the transfer_rate argument. 's' means
            seconds, 'min' means minutes, 'h' means hours, 'a' means Julian years.

        """
        if self.transfer_rates is None:
            self.transfer_rates = TransferRates(self.operator, self.operator.model)

        self.transfer_rates.set_transfer_rate(material, components, transfer_rate,
                                      transfer_rate_units, destination_material)

    def add_batchwise(self, obj, attr, **kwargs):
        """Add batchwise operation to integrator scheme.

        Parameters
        ----------
        obj : openmc.Cell or openmc.Material object or id or str name
            Cell or Materials identifier to where add batchwise scheme
        attr : str
            Attribute to specify the type of batchwise scheme. Accepted values
            are: 'translation', 'rotation', 'temperature' for an openmc.Cell
            object; 'refuel' for an openmc.Material object.
        **kwargs
            keyword arguments that are passed to the batchwise class.

        """
        check_value('attribute', attr, ('translation', 'rotation',
                                        'temperature', 'refuel','dilute',
                                        'addition'))
        if attr in ('translation', 'rotation'):
            batchwise = BatchwiseCellGeometrical
        elif attr == 'temperature':
            batchwise = BatchwiseCellTemperature
        elif attr == 'refuel':
            batchwise = BatchwiseMaterialRefuel
        elif attr == 'dilute':
            batchwise = BatchwiseMaterialDilute
        elif attr == 'addition':
            batchwise = BatchwiseMaterialAdd

        batchwise_inst = batchwise.from_params(obj, attr, self.operator,
                                           self.operator.model, **kwargs)
        if self.batchwise is None:
            self.batchwise = batchwise_inst
        else:
            if not isinstance(self.batchwise, list):
                self.batchwise = [self.batchwise]
            self.batchwise.append(batchwise_inst)

    def add_batchwise_scheme(self, scheme_name, **kwargs):
        """Add batchwise wrapper to integrator scheme, after calls to
        meth:`add_batchwise`.

        Parameters
        ----------
        wrap_name : str
            wrap_name of wrapper function. So far only '1' or '2'.
        **kwargs
            keyword arguments that are passed to the batchwise wrapper class.

        """
        if scheme_name == 'std':
            self.batchwise = BatchwiseSchemeStd(self.batchwise, len(self), **kwargs)
        elif scheme_name == 'refuel':
            self.batchwise = BatchwiseSchemeRefuel(self.batchwise, **kwargs)
        elif scheme_name == 'flex':
            self.batchwise = BatchwiseSchemeFlex(self.batchwise, len(self), **kwargs)

    def add_density_function(self, mats, density_func, oxidation_states):
        self.batchwise.set_density_function(mats, density_func, oxidation_states)

    def add_redox(self, mat, buffer, oxidation_states):
        self.transfer_rates.set_redox(mat, buffer, oxidation_states)

    def add_material(self, mat, value, mat_vector, timestep, quantity='grams'):
        self.batchwise.add_material(mat, value, mat_vector, timestep,
                                    quantity)
@add_params
class SIIntegrator(Integrator):
    r"""Abstract class for the Stochastic Implicit Euler integrators

    Does not provide a ``__call__`` method, but scales and resets
    the number of particles used in initial transport calculation
    """

    _params = r"""
    Parameters
    ----------
    operator : openmc.deplete.abc.TransportOperator
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
        an area in [cm^2]. Either ``power``, ``power_density``, or
        ``source_rates`` must be specified.
    power_density : float or iterable of float, optional
        Power density of the reactor in [W/gHM]. It is multiplied by
        initial heavy metal inventory to get total power if ``power``
        is not specified.
    source_rates : float or iterable of float, optional
        Source rate in [neutron/sec] or neutron flux in [neutron/s-cm^2] for
        each interval in :attr:`timesteps`

        .. versionadded:: 0.12.1
    timestep_units : {'s', 'min', 'h', 'd', 'MWd/kg'}
        Units for values specified in the `timesteps` argument. 's' means
        seconds, 'min' means minutes, 'h' means hours, and 'MWd/kg' indicates
        that the values are given in burnup (MW-d of energy deposited per
        kilogram of initial heavy metal).
    n_steps : int, optional
        Number of stochastic iterations per depletion interval.
        Must be greater than zero. Default : 10
    solver : str or callable, optional
        If a string, must be the name of the solver responsible for
        solving the Bateman equations.  Current options are:

            * ``cram16`` - 16th order IPF CRAM
            * ``cram48`` - 48th order IPF CRAM [default]

        If a function or other callable, must adhere to the requirements in
        :attr:`solver`.

        .. versionadded:: 0.12

    Attributes
    ----------
    operator : openmc.deplete.abc.TransportOperator
        Operator to perform transport simulations
    chain : openmc.deplete.Chain
        Depletion chain
    timesteps : iterable of float
        Size of each depletion interval in [s]
    power : iterable of float
        Power of the reactor in [W] for each interval in :attr:`timesteps`
    n_steps : int
        Number of stochastic iterations per depletion interval
    solver : callable
        Function that will solve the Bateman equations
        :math:`\frac{\partial}{\partial t}\vec{n} = A_i\vec{n}_i` with a step
        size :math:`t_i`. Can be configured using the ``solver`` argument.
        User-supplied functions are expected to have the following signature:
        ``solver(A, n0, t) -> n1`` where

            * ``A`` is a :class:`scipy.sparse.csc_matrix` making up the
              depletion matrix
            * ``n0`` is a 1-D :class:`numpy.ndarray` of initial compositions
              for a given material in atoms/cm3
            * ``t`` is a float of the time step size in seconds, and
            * ``n1`` is a :class:`numpy.ndarray` of compositions at the
              next time step. Expected to be of the same shape as ``n0``

        .. versionadded:: 0.12

    """

    def __init__(
            self,
            operator: TransportOperator,
            timesteps: Sequence[float],
            power: Optional[Union[float, Sequence[float]]] = None,
            power_density: Optional[Union[float, Sequence[float]]] = None,
            source_rates: Optional[Sequence[float]] = None,
            timestep_units: str = 's',
            n_steps: int = 10,
            solver: str = "cram48"
        ):
        check_type("n_steps", n_steps, Integral)
        check_greater_than("n_steps", n_steps, 0)
        super().__init__(
            operator, timesteps, power, power_density, source_rates,
            timestep_units=timestep_units, solver=solver)
        self.n_steps = n_steps

    def _get_bos_data_from_operator(self, step_index, step_power, n_bos):
        reset_particles = False
        if step_index == 0 and hasattr(self.operator, "settings"):
            reset_particles = True
            self.operator.settings.particles *= self.n_steps
        inherited = super()._get_bos_data_from_operator(
            step_index, step_power, n_bos)
        if reset_particles:
            self.operator.settings.particles //= self.n_steps
        return inherited

    def integrate(
            self,
            output: bool = True,
            path: PathLike = "depletion_results.h5"
        ):
        """Perform the entire depletion process across all steps

        Parameters
        ----------
        output : bool, optional
            Indicate whether to display information about progress
        path : PathLike
            Path to file to write. Defaults to 'depletion_results.h5'.

            .. versionadded:: 0.15.0
        """
        with change_directory(self.operator.output_dir):
            n = self.operator.initial_condition()
            t, self._i_res = self._get_start_data()

            for i, (dt, p) in enumerate(self):
                if output:
                    print(f"[openmc.deplete] t={t} s, dt={dt} s, source={p}")

                if i == 0:
                    if self.operator.prev_res is None:
                        n, res = self._get_bos_data_from_operator(i, p, n)
                    else:
                        n, res = self._get_bos_data_from_restart(p, n)
                else:
                    # Pull rates, k from previous iteration w/o
                    # re-running transport
                    res = res_list[-1]  # defined in previous i iteration

                proc_time, n_list, res_list = self(n, res.rates, dt, p, i)

                # Insert BOS concentration, transport results
                n_list.insert(0, n)
                res_list.insert(0, res)

                # Remove actual EOS concentration for next step
                n = n_list.pop()

                StepResult.save(self.operator, n_list, res_list, [t, t + dt],
                             p, self._i_res + i, proc_time, path)

                t += dt

            # No final simulation for SIE, use last iteration results
            StepResult.save(self.operator, [n], [res_list[-1]], [t, t],
                         p, self._i_res + len(self), proc_time, path)
            self.operator.write_bos_data(self._i_res + len(self))

        self.operator.finalize()


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
        A : scipy.sparse.csc_matrix
            Sparse transmutation matrix ``A[j, i]`` describing rates at
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
