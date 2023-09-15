"""
Classes for collecting and calculating quantities for reaction rate operators
"""
import bisect
from abc import abstractmethod
from collections import defaultdict
from copy import deepcopy
from itertools import product
from numbers import Real
import sys

from numpy import dot, zeros, newaxis, asarray

from openmc.mpi import comm
from openmc.checkvalue import check_type, check_greater_than
from openmc.data import JOULE_PER_EV, REACTION_MT
from openmc.lib import (
    Tally, MaterialFilter, EnergyFilter, EnergyFunctionFilter, load_nuclide)
import openmc.lib
from .abc import (
    ReactionRateHelper, NormalizationHelper, FissionYieldHelper)

__all__ = (
    "DirectReactionRateHelper", "ChainFissionHelper", "EnergyScoreHelper"
    "SourceRateHelper", "TalliedFissionYieldHelper",
    "ConstantFissionYieldHelper", "FissionYieldCutoffHelper",
    "AveragedFissionYieldHelper", "FluxCollapseHelper")


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
            :class:`openmc.deplete.CoupledOperator` that uses this helper
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


# -------------------------------------
# Helpers for generating reaction rates
# -------------------------------------


class DirectReactionRateHelper(ReactionRateHelper):
    """Class for generating one-group reaction rates with direct tallies

    This class generates reaction rate tallies for each nuclide and
    transmutation reaction relevant for a depletion calculation.

    Parameters
    ----------
    n_nucs : int
        Number of burnable nuclides tracked by
        :class:`openmc.deplete.CoupledOperator`
    n_react : int
        Number of reactions tracked by an instance of
        :class:`openmc.deplete.CoupledOperator`

    Attributes
    ----------
    nuclides : list of str
        All nuclides with desired reaction rates.
    """
    def __init__(self, n_nuc, n_react):
        super().__init__(n_nuc, n_react)
        self._rate_tally = None

        # Automatically pre-calculate reaction rates for depletion
        openmc.lib.settings.need_depletion_rx = True

    @ReactionRateHelper.nuclides.setter
    def nuclides(self, nuclides):
        ReactionRateHelper.nuclides.fset(self, nuclides)
        self._rate_tally.nuclides = nuclides

    def generate_tallies(self, materials, scores):
        """Produce one-group reaction rate tally

        Uses the :mod:`openmc.lib` to generate a tally of relevant reactions
        across all burnable materials.

        Parameters
        ----------
        materials : iterable of :class:`openmc.lib.Material`
            Burnable materials in the problem. Used to construct a
            :class:`openmc.lib.MaterialFilter`
        scores : iterable of str
            Reaction identifiers, e.g. ``"(n, fission)"``, ``"(n, gamma)"``,
            needed for the reaction rate tally.
        """
        self._rate_tally = Tally()
        self._rate_tally.writable = False
        self._rate_tally.scores = scores
        self._rate_tally.filters = [MaterialFilter(materials)]
        self._rate_tally.multiply_density = False
        self._rate_tally_means_cache = None

    @property
    def rate_tally_means(self):
        """The mean results of the tally of every material's reaction rates for this cycle
        """
        # If the mean cache is empty, fill it once with this transport cycle's results
        if self._rate_tally_means_cache is None:
            self._rate_tally_means_cache = self._rate_tally.mean
        return self._rate_tally_means_cache

    def reset_tally_means(self):
        """Reset the cached mean rate tallies.
        .. note::

                This step must be performed after each transport cycle
        """
        self._rate_tally_means_cache = None

    def get_material_rates(self, mat_index, nuc_index, rx_index):
        """Return an array of reaction rates for a material

        Parameters
        ----------
        mat_index : int
            Index for the material
        nuc_index : iterable of int
            Index for each nuclide in :attr:`nuclides` in the
            desired reaction rate matrix
        rx_index : iterable of int
            Index for each reaction scored in the tally

        Returns
        -------
        rates : numpy.ndarray
            Array with shape ``(n_nuclides, n_rxns)`` with the
            reaction rates in this material
        """
        self._results_cache.fill(0.0)
        full_tally_res = self.rate_tally_means[mat_index]
        for i_tally, (i_nuc, i_rx) in enumerate(product(nuc_index, rx_index)):
            self._results_cache[i_nuc, i_rx] = full_tally_res[i_tally]

        return self._results_cache


class FluxCollapseHelper(ReactionRateHelper):
    """Class that generates one-group reaction rates using multigroup flux

    This class generates a multigroup flux tally that is used afterward to
    calculate a one-group reaction rate by collapsing it with continuous-energy
    cross section data. Additionally, select nuclides/reactions can be treated
    with a direct reaction rate tally when using a multigroup flux spectrum
    would not be sufficiently accurate. This is often the case for (n,gamma) and
    fission reactions.

    .. versionadded:: 0.12.1

    Parameters
    ----------
    n_nucs : int
        Number of burnable nuclides tracked by
        :class:`openmc.deplete.CoupledOperator`
    n_react : int
        Number of reactions tracked by :class:`openmc.deplete.CoupledOperator`
    energies : iterable of float
        Energy group boundaries for flux spectrum in [eV]
    reactions : iterable of str
        Reactions for which rates should be directly tallied
    nuclides : iterable of str
        Nuclides for which some reaction rates should be directly tallied. If
        None, then ``reactions`` will be used for all nuclides.

    Attributes
    ----------
    nuclides : list of str
        All nuclides with desired reaction rates.

    """
    def __init__(self, n_nucs, n_reacts, energies, reactions=None, nuclides=None):
        super().__init__(n_nucs, n_reacts)
        self._energies = asarray(energies)
        self._reactions_direct = list(reactions) if reactions is not None else []
        self._nuclides_direct = list(nuclides) if nuclides is not None else None

    @ReactionRateHelper.nuclides.setter
    def nuclides(self, nuclides):
        ReactionRateHelper.nuclides.fset(self, nuclides)
        if self._reactions_direct and self._nuclides_direct is None:
            self._rate_tally.nuclides = nuclides

        # Make sure nuclide data is loaded
        for nuclide in self.nuclides:
            if nuclide not in openmc.lib.nuclides:
                openmc.lib.load_nuclide(nuclide)

    def generate_tallies(self, materials, scores):
        """Produce multigroup flux spectrum tally

        Uses the :mod:`openmc.lib` module to generate a multigroup flux tally
        for each burnable material.

        Parameters
        ----------
        materials : iterable of :class:`openmc.Material`
            Burnable materials in the problem. Used to construct a
            :class:`openmc.MaterialFilter`
        scores : iterable of str
            Reaction identifiers, e.g. ``"(n, fission)"``, ``"(n, gamma)"``,
            needed for the reaction rate tally.
        """
        self._materials = materials

        # adds an entry for fisson to the dictionary of reactions
        self._mts = [REACTION_MT[x] for x in scores]
        self._scores = scores

        # Create flux tally with material and energy filters
        self._flux_tally = Tally()
        self._flux_tally.writable = False
        self._flux_tally.filters = [
            MaterialFilter(materials),
            EnergyFilter(self._energies)
        ]
        self._flux_tally.scores = ['flux']
        self._flux_tally_means_cache = None

        # Create reaction rate tally
        if self._reactions_direct:
            self._rate_tally = Tally()
            self._rate_tally.writable = False
            self._rate_tally.scores = self._reactions_direct
            self._rate_tally.filters = [MaterialFilter(materials)]
            self._rate_tally.multiply_density = False
            self._rate_tally_means_cache = None
            if self._nuclides_direct is not None:
                # check if any direct tally nuclides are requested that are not
                # already loaded with the materials. Load separately if so.
                mat_nuclides = {n for mat in materials for n in mat.nuclides}
                extra_nuclides = set(self._nuclides_direct) - mat_nuclides
                for nuc in extra_nuclides:
                    load_nuclide(nuc)
                self._rate_tally.nuclides = self._nuclides_direct

    @property
    def rate_tally_means(self):
        """The mean results of the tally of every material's reaction rates for this cycle
        """
        # If the mean cache is empty, fill it once with this transport cycle's results
        if self._rate_tally_means_cache is None:
            self._rate_tally_means_cache = self._rate_tally.mean
        return self._rate_tally_means_cache

    @property
    def flux_tally_means(self):
        # If the mean cache is empty, fill it once for this transport cycle's results
        if self._flux_tally_means_cache is None:
            self._flux_tally_means_cache = self._flux_tally.mean
        return self._flux_tally_means_cache

    def reset_tally_means(self):
        """Reset the cached mean rate and flux tallies.
        .. note::

                This step must be performed after each transport cycle
        """
        self._flux_tally_means_cache = None
        if self._reactions_direct:
            self._rate_tally_means_cache = None

    def get_material_rates(self, mat_index, nuc_index, react_index):
        """Return an array of reaction rates for a material

        Parameters
        ----------
        mat_index : int
            Index for material
        nuc_index : iterable of int
            Index for each nuclide in :attr:`nuclides` in the
            desired reaction rate matrix
        react_index : iterable of int
            Index for each reaction scored in the tally

        Returns
        -------
        rates : numpy.ndarray
            Array with shape ``(n_nuclides, n_rxns)`` with the reaction rates in
            this material

        """
        self._results_cache.fill(0.0)

        # Get flux for specified material
        shape = (len(self._materials), len(self._energies) - 1)
        mean_value = self.flux_tally_means.reshape(shape)
        flux = mean_value[mat_index]

        # Get direct reaction rates
        if self._reactions_direct:
            nuclides_direct = self._rate_tally.nuclides
            shape = (len(nuclides_direct), len(self._reactions_direct))
            rx_rates = self.rate_tally_means[mat_index].reshape(shape)

        mat = self._materials[mat_index]

        for name, i_nuc in zip(self.nuclides, nuc_index):
            for mt, score, i_rx in zip(self._mts, self._scores, react_index):
                if score in self._reactions_direct and name in nuclides_direct:
                    # Determine index in rx_rates
                    i_rx_direct = self._reactions_direct.index(score)
                    i_nuc_direct = nuclides_direct.index(name)

                    # Get reaction rate from tally
                    self._results_cache[i_nuc, i_rx] = rx_rates[i_nuc_direct, i_rx_direct]
                else:
                    # Use flux to collapse reaction rate (per N)
                    nuc = openmc.lib.nuclides[name]
                    rate_per_nuc = nuc.collapse_rate(
                        mt, mat.temperature, self._energies, flux)

                    self._results_cache[i_nuc, i_rx] = rate_per_nuc

        return self._results_cache


# ------------------------------------------
# Helpers for obtaining normalization factor
# ------------------------------------------


class EnergyNormalizationHelper(NormalizationHelper):
    """Compute energy-based normalization."""

    def reset(self):
        """Reset energy produced prior to unpacking tallies"""
        self._energy = 0.0

    def factor(self, source_rate):
        # Reduce energy produced from all processes
        # J / source neutron
        energy = comm.allreduce(self._energy) * JOULE_PER_EV

        # Guard against divide by zero
        if energy == 0:
            if comm.rank == 0:
                sys.stderr.flush()
                print("No energy reported from OpenMC tallies. Do your HDF5 "
                      "files have heating data?\n", file=sys.stderr, flush=True)
            comm.barrier()
            comm.Abort(1)

        # Return normalization factor for scaling reaction rates. In this case,
        # the source rate is the power in [W], so [W] / [J/src] = [src/s]
        return source_rate / energy


class ChainFissionHelper(EnergyNormalizationHelper):
    """Computes normalization using fission Q values from depletion chain

    Attributes
    ----------
    nuclides : list of str
        All nuclides with desired reaction rates. Ordered to be
        consistent with :class:`openmc.deplete.CoupledOperator`
    energy : float
        Total energy [J/s/source neutron] produced in a transport simulation.
        Updated in the material iteration with :meth:`update`.
    """

    def __init__(self):
        super().__init__()
        self._fission_q_vector = None

    def prepare(self, chain_nucs, rate_index):
        """Populate the fission Q value vector from a chain.

        Parameters
        ----------
        chain_nucs : iterable of :class:`openmc.deplete.Nuclide`
            Nuclides used in this depletion chain. Do not need
            to be ordered
        rate_index : dict of str to int
            Dictionary mapping names of nuclides, e.g. ``"U235"``,
            to a corresponding index in the desired fission Q
            vector.
        """
        if (self._fission_q_vector is not None
                and self._fission_q_vector.shape == (len(rate_index),)):
            return

        fission_qs = zeros(len(rate_index))

        for nuclide in chain_nucs:
            if nuclide.name in rate_index:
                for rx in nuclide.reactions:
                    if rx.type == "fission":
                        fission_qs[rate_index[nuclide.name]] = rx.Q
                        break

        self._fission_q_vector = fission_qs

    def update(self, fission_rates):
        """Update energy produced with fission rates in a material

        Parameters
        ----------
        fission_rates : numpy.ndarray
            fission reaction rate for each isotope in the specified
            material. Should be ordered corresponding to initial
            ``rate_index`` used in :meth:`prepare`
        """
        self._energy += dot(fission_rates, self._fission_q_vector)


class EnergyScoreHelper(EnergyNormalizationHelper):
    """Class responsible for obtaining system energy via a tally score

    Parameters
    ----------
    score : string
        Valid score to use when obtaining system energy from OpenMC.
        Defaults to "heating-local"

    Attributes
    ----------
    nuclides : list of str
        List of nuclides with reaction rates. Not needed, but provided
        for a consistent API across other :class:`NormalizationHelper`
    energy : float
        System energy [eV] computed from the tally. Will be zero for
        all MPI processes that are not the "master" process to avoid
        artificially increasing the tallied energy.
    score : str
        Score used to obtain system energy

    """

    def __init__(self, score="heating-local"):
        super().__init__()
        self.score = score
        self._tally = None

    def prepare(self, *args, **kwargs):
        """Create a tally for system energy production

        Input arguments are not used, as the only information needed
        is :attr:`score`

        """
        self._tally = Tally()
        self._tally.writable = False
        self._tally.scores = [self.score]

    def reset(self):
        """Obtain system energy from tally

        Only the master process, ``comm.rank == 0`` will
        have a non-zero :attr:`energy` taken from the tally.
        This avoids accidentally scaling the system power by
        the number of MPI processes
        """
        super().reset()
        if comm.rank == 0:
            self._energy = self._tally.mean[0, 0]


class SourceRateHelper(NormalizationHelper):
    def prepare(self, *args, **kwargs):
        pass

    def factor(self, source_rate):
        return source_rate

# ------------------------------------
# Helper for collapsing fission yields
# ------------------------------------


class ConstantFissionYieldHelper(FissionYieldHelper):
    """Class that uses a single set of fission yields on each isotope

    Parameters
    ----------
    chain_nuclides : iterable of openmc.deplete.Nuclide
        Nuclides tracked in the depletion chain. All nuclides are
        not required to have fission yield data.
    energy : float, optional
        Key in :attr:`openmc.deplete.Nuclide.yield_data` corresponding
        to the desired set of fission yield data. Typically one of
        ``{0.0253, 500000, 14000000}`` corresponding to 0.0253 eV,
        500 keV, and 14 MeV yield libraries. If the specific key is not
        found, will fall back to closest energy present.
        Default: 0.0253 eV for thermal yields

    Attributes
    ----------
    constant_yields : collections.defaultdict
        Fission yields for all nuclides that only have one set of
        fission yield data. Dictionary of form ``{str: {str: float}}``
        representing yields for ``{parent: {product: yield}}``. Default
        return object is an empty dictionary
    energy : float
        Energy of fission yield libraries.
    """

    def __init__(self, chain_nuclides, energy=0.0253):
        check_type("energy", energy, Real)
        check_greater_than("energy", energy, 0.0, equality=True)
        self._energy = energy
        super().__init__(chain_nuclides)
        # Iterate over all nuclides with > 1 set of yields
        for name, nuc in self._chain_nuclides.items():
            yield_data = nuc.yield_data.get(energy)
            if yield_data is not None:
                self._constant_yields[name] = yield_data
                continue
            # Specific energy not found, use closest energy
            distances = [abs(energy - ene) for ene in nuc.yield_energies]
            min_E = min(nuc.yield_energies, key=lambda e: abs(e - energy))
            self._constant_yields[name] = nuc.yield_data[min_E]

    @classmethod
    def from_operator(cls, operator, **kwargs):
        """Return a new ConstantFissionYieldHelper using operator data

        All keyword arguments should be identical to their counterpart
        in the main ``__init__`` method

        Parameters
        ----------
        operator : openmc.deplete.abc.TransportOperator
            operator with a depletion chain
        kwargs:
            Additional keyword arguments to be used in construction

        Returns
        -------
        ConstantFissionYieldHelper
        """
        return cls(operator.chain.nuclides, **kwargs)

    @property
    def energy(self):
        return self._energy

    def weighted_yields(self, _local_mat_index=None):
        """Return fission yields for all nuclides requested

        Parameters
        ----------
        _local_mat_index : int, optional
            Current material index. Not used since all yields are
            constant

        Returns
        -------
        library : collections.defaultdict
            Dictionary of ``{parent: {product: fyield}}``
        """
        return self.constant_yields


class FissionYieldCutoffHelper(TalliedFissionYieldHelper):
    """Helper that computes fission yields based on a cutoff energy

    Tally fission rates above and below the cutoff energy.
    Assume that all fissions below cutoff energy have use thermal fission
    product yield distributions, while all fissions above use a faster
    set of yield distributions.

    Uses a limit of 20 MeV for tallying fission.

    Parameters
    ----------
    chain_nuclides : iterable of openmc.deplete.Nuclide
        Nuclides tracked in the depletion chain. All nuclides are
        not required to have fission yield data.
    n_bmats : int, optional
        Number of burnable materials tracked in the problem
    cutoff : float, optional
        Cutoff energy in [eV] below which all fissions will be
        use thermal yields. All other fissions will use a
        faster set of yields. Default: 112 [eV]
    thermal_energy : float, optional
        Energy of yield data corresponding to thermal yields.
        Default: 0.0253 [eV]
    fast_energy : float, optional
        Energy of yield data corresponding to fast yields.

    Attributes
    ----------
    n_bmats : int
        Number of burnable materials tracked in the problem.
        Must be set prior to generating tallies
    thermal_yields : dict
        Dictionary of the form ``{parent: {product: yield}}``
        with thermal yields
    fast_yields : dict
        Dictionary of the form ``{parent: {product: yield}}``
        with fast yields
    constant_yields : collections.defaultdict
        Fission yields for all nuclides that only have one set of
        fission yield data. Dictionary of form ``{str: {str: float}}``
        representing yields for ``{parent: {product: yield}}``. Default
        return object is an empty dictionary
    results : numpy.ndarray
        Array of fission rate fractions with shape
        ``(n_mats, 2, n_nucs)``. ``results[:, 0]``
        corresponds to the fraction of all fissions
        that occurred below ``cutoff``. The number
        of materials in the first axis corresponds
        to the number of materials burned by the
        :class:`openmc.deplete.CoupledOperator`
    """

    def __init__(self, chain_nuclides, n_bmats, cutoff=112.0,
                 thermal_energy=0.0253, fast_energy=500.0e3):
        check_type("cutoff", cutoff, Real)
        check_type("thermal_energy", thermal_energy, Real)
        check_type("fast_energy", fast_energy, Real)
        check_greater_than("thermal_energy", thermal_energy, 0.0, equality=True)
        check_greater_than("cutoff", cutoff, thermal_energy, equality=False)
        check_greater_than("fast_energy", fast_energy, cutoff, equality=False)
        self.n_bmats = n_bmats
        super().__init__(chain_nuclides)
        self._cutoff = cutoff
        self._thermal_yields = {}
        self._fast_yields = {}
        convert_to_constant = set()
        for name, nuc in self._chain_nuclides.items():
            yields = nuc.yield_data
            energies = nuc.yield_energies
            thermal = yields.get(thermal_energy)
            fast = yields.get(fast_energy)
            if thermal is None or fast is None:
                if cutoff <= energies[0]:
                    # use lowest energy yields as constant
                    self._constant_yields[name] = yields[energies[0]]
                    convert_to_constant.add(name)
                    continue
                if cutoff >= energies[-1]:
                    # use highest energy yields as constant
                    self._constant_yields[name] = yields[energies[-1]]
                    convert_to_constant.add(name)
                    continue
                cutoff_ix = bisect.bisect_left(energies, cutoff)
                # find closest energy to requested thermal, fast energies
                if thermal is None:
                    min_E = min(energies[:cutoff_ix],
                                key=lambda e: abs(e - thermal_energy))
                    thermal = yields[min_E]
                if fast is None:
                    min_E = min(energies[cutoff_ix:],
                                key=lambda e: abs(e - fast_energy))
                    fast = yields[min_E]
            self._thermal_yields[name] = thermal
            self._fast_yields[name] = fast
        for name in convert_to_constant:
            self._chain_nuclides.pop(name)

    @classmethod
    def from_operator(cls, operator, **kwargs):
        """Construct a helper from an operator

        All keyword arguments should be identical to their counterpart
        in the main ``__init__`` method

        Parameters
        ----------
        operator : openmc.deplete.CoupledOperator
            Operator with a chain and burnable materials
        kwargs:
            Additional keyword arguments to be used in construction

        Returns
        -------
        FissionYieldCutoffHelper

        """
        return cls(operator.chain.nuclides, len(operator.burnable_mats),
                   **kwargs)

    def generate_tallies(self, materials, mat_indexes):
        """Use C API to produce a fission rate tally in burnable materials

        Include a :class:`openmc.lib.EnergyFilter` to tally fission rates
        above and below cutoff energy.

        Parameters
        ----------
        materials : iterable of :class:`openmc.lib.Material`
            Materials to be used in :class:`openmc.lib.MaterialFilter`
        mat_indexes : iterable of int
            Indices of tallied materials that will have their fission
            yields computed by this helper. Necessary as the
            :class:`openmc.deplete.CoupledOperator` that uses this helper
            may only burn a subset of all materials when running
            in parallel mode.
        """
        super().generate_tallies(materials, mat_indexes)
        energy_filter = EnergyFilter([0.0, self._cutoff, self._upper_energy])
        self._fission_rate_tally.filters = (
            self._fission_rate_tally.filters + [energy_filter])

    def unpack(self):
        """Obtain fast and thermal fission fractions from tally"""
        if not self._tally_nucs or self._local_indexes.size == 0:
            self.results = None
            return
        fission_rates = self._fission_rate_tally.mean.reshape(
            self.n_bmats, 2, len(self._tally_nucs))
        self.results = fission_rates[self._local_indexes]
        total_fission = self.results.sum(axis=1)
        nz_mat, nz_nuc = total_fission.nonzero()
        self.results[nz_mat, :, nz_nuc] /= total_fission[nz_mat, newaxis, nz_nuc]

    def weighted_yields(self, local_mat_index):
        """Return fission yields for a specific material

        For nuclides with both yield data above and below
        the cutoff energy, the effective yield for nuclide ``A``
        will be a weighted sum of fast and thermal yields. The
        weights will be the fraction of ``A`` fission events
        in the above and below the cutoff energy.

        If ``A`` has fission product distribution ``F``
        for fast fissions and ``T`` for thermal fissions, and
        70% of ``A`` fissions are considered thermal, then
        the effective fission product yield distributions
        for ``A`` is ``0.7 * T + 0.3 * F``

        Parameters
        ----------
        local_mat_index : int
            Index for specific burnable material. Effective
            yields will be produced using
            ``self.results[local_mat_index]``

        Returns
        -------
        library : collections.defaultdict
            Dictionary of ``{parent: {product: fyield}}``
        """
        yields = self.constant_yields
        if not self._tally_nucs:
            return yields
        rates = self.results[local_mat_index]
        # iterate over thermal then fast yields, prefer __mul__ to __rmul__
        for therm_frac, fast_frac, nuc in zip(rates[0], rates[1], self._tally_nucs):
            yields[nuc.name] = (self._thermal_yields[nuc.name] * therm_frac
                                + self._fast_yields[nuc.name] * fast_frac)
        return yields

    @property
    def thermal_yields(self):
        return deepcopy(self._thermal_yields)

    @property
    def fast_yields(self):
        return deepcopy(self._fast_yields)


class AveragedFissionYieldHelper(TalliedFissionYieldHelper):
    r"""Class that computes fission yields based on average fission energy

    Computes average energy at which fission events occurred with

    .. math::

        \bar{E} = \frac{
            \int_0^\infty E\sigma_f(E)\phi(E)dE
        }{
            \int_0^\infty\sigma_f(E)\phi(E)dE
        }

    If the average energy for a nuclide is below the lowest energy
    with yield data, that set of fission yields is taken.
    Conversely, if the average energy is above the highest energy
    with yield data, that set of fission yields is used.
    For the case where the average energy is between two sets
    of yields, the effective fission yield computed by
    linearly interpolating between yields provided at the
    nearest energies above and below the average.

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
    results : None or numpy.ndarray
        If tallies have been generated and unpacked, then the array will
        have shape ``(n_mats, n_tnucs)``, where ``n_mats`` is the number
        of materials where fission reactions were tallied and ``n_tnucs``
        is the number of nuclides with multiple sets of fission yields.
        Data in the array are the average energy of fission events for
        tallied nuclides across burnable materials.
    """

    def __init__(self, chain_nuclides):
        super().__init__(chain_nuclides)
        self._weighted_tally = None

    def generate_tallies(self, materials, mat_indexes):
        """Construct tallies to determine average energy of fissions

        Parameters
        ----------
        materials : iterable of :class:`openmc.lib.Material`
            Materials to be used in :class:`openmc.lib.MaterialFilter`
        mat_indexes : iterable of int
            Indices of tallied materials that will have their fission
            yields computed by this helper. Necessary as the
            :class:`openmc.deplete.CoupledOperator` that uses this helper
            may only burn a subset of all materials when running
            in parallel mode.
        """
        super().generate_tallies(materials, mat_indexes)
        fission_tally = self._fission_rate_tally
        filters = fission_tally.filters

        ene_filter = EnergyFilter([0, self._upper_energy])
        fission_tally.filters = filters + [ene_filter]

        func_filter = EnergyFunctionFilter()
        func_filter.set_data((0, self._upper_energy), (0, self._upper_energy))
        weighted_tally = Tally()
        weighted_tally.writable = False
        weighted_tally.scores = ['fission']
        weighted_tally.filters = filters + [func_filter]
        self._weighted_tally = weighted_tally

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
        nuclides : tuple of str
            Union of input nuclides and those that have multiple sets
            of yield data.  Sorted by nuclide name

        Raises
        ------
        AttributeError
            If tallies not generated
        """
        tally_nucs = super().update_tally_nuclides(nuclides)
        self._weighted_tally.nuclides = tally_nucs
        return tally_nucs

    def unpack(self):
        """Unpack tallies and populate :attr:`results` with average energies"""
        if not self._tally_nucs or self._local_indexes.size == 0:
            self.results = None
            return
        fission_results = (
            self._fission_rate_tally.mean[self._local_indexes])
        self.results = (
            self._weighted_tally.mean[self._local_indexes]).copy()
        nz_mat, nz_nuc = fission_results.nonzero()
        self.results[nz_mat, nz_nuc] /= fission_results[nz_mat, nz_nuc]

    def weighted_yields(self, local_mat_index):
        """Return fission yields for a specific material

        Use the computed average energy of fission
        events to determine fission yields. If average
        energy is between two sets of yields, linearly
        interpolate between the two.
        Otherwise take the closet set of yields.

        Parameters
        ----------
        local_mat_index : int
            Index for specific burnable material. Effective
            yields will be produced using
            ``self.results[local_mat_index]``

        Returns
        -------
        library : collections.defaultdict
            Dictionary of ``{parent: {product: fyield}}``. Default return
            value is an empty dictionary
        """
        if not self._tally_nucs:
            return self.constant_yields
        mat_yields = defaultdict(dict)
        average_energies = self.results[local_mat_index]
        for avg_e, nuc in zip(average_energies, self._tally_nucs):
            nuc_energies = nuc.yield_energies
            if avg_e <= nuc_energies[0]:
                mat_yields[nuc.name] = nuc.yield_data[nuc_energies[0]]
                continue
            if avg_e >= nuc_energies[-1]:
                mat_yields[nuc.name] = nuc.yield_data[nuc_energies[-1]]
                continue
            # in-between two energies
            # linear search since there are usually ~3 energies
            for ix, ene in enumerate(nuc_energies[:-1]):
                if nuc_energies[ix + 1] > avg_e:
                    break
            lower, upper = nuc_energies[ix:ix + 2]
            fast_frac = (avg_e - lower) / (upper - lower)
            mat_yields[nuc.name] = (
                nuc.yield_data[lower] * (1 - fast_frac)
                + nuc.yield_data[upper] * fast_frac)
        mat_yields.update(self.constant_yields)
        return mat_yields

    @classmethod
    def from_operator(cls, operator, **kwargs):
        """Return a new helper with data from an operator

        All keyword arguments should be identical to their counterpart
        in the main ``__init__`` method

        Parameters
        ----------
        operator : openmc.deplete.CoupledOperator
            Operator with a depletion chain
        kwargs :
            Additional keyword arguments to be used in construction

        Returns
        -------
        AveragedFissionYieldHelper
        """
        return cls(operator.chain.nuclides)
