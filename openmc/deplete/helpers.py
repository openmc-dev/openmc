"""
Class for normalizing fission energy deposition
"""
from copy import deepcopy
from itertools import product
from numbers import Real
import bisect
from collections import defaultdict

from numpy import dot, zeros, newaxis

from . import comm
from openmc.checkvalue import check_type, check_greater_than
from openmc.lib import (
    Tally, MaterialFilter, EnergyFilter, EnergyFunctionFilter)
from .abc import (
    ReactionRateHelper, EnergyHelper, FissionYieldHelper,
    TalliedFissionYieldHelper)

__all__ = (
    "DirectReactionRateHelper", "ChainFissionHelper",
    "ConstantFissionYieldHelper", "FissionYieldCutoffHelper",
    "AveragedFissionYieldHelper")

# -------------------------------------
# Helpers for generating reaction rates
# -------------------------------------


class DirectReactionRateHelper(ReactionRateHelper):
    """Class that generates tallies for one-group rates

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

    def generate_tallies(self, materials, scores):
        """Produce one-group reaction rate tally

        Uses the :mod:`openmc.lib` to generate a tally
        of relevant reactions across all burnable materials.

        Parameters
        ----------
        materials : iterable of :class:`openmc.Material`
            Burnable materials in the problem. Used to
            construct a :class:`openmc.MaterialFilter`
        scores : iterable of str
            Reaction identifiers, e.g. ``"(n, fission)"``,
            ``"(n, gamma)"``, needed for the reaction rate tally.
        """
        self._rate_tally = Tally()
        self._rate_tally.writable = False
        self._rate_tally.scores = scores
        self._rate_tally.filters = [MaterialFilter(materials)]

    def get_material_rates(self, mat_id, nuc_index, react_index):
        """Return an array of reaction rates for a material

        Parameters
        ----------
        mat_id : int
            Unique ID for the requested material
        nuc_index : iterable of int
            Index for each nuclide in :attr:`nuclides` in the
            desired reaction rate matrix
        react_index : iterable of int
            Index for each reaction scored in the tally

        Returns
        -------
        rates : numpy.ndarray
            Array with shape ``(n_nuclides, n_rxns)`` with the
            reaction rates in this material
        """
        self._results_cache.fill(0.0)
        full_tally_res = self._rate_tally.results[mat_id, :, 1]
        for i_tally, (i_nuc, i_react) in enumerate(
                product(nuc_index, react_index)):
            self._results_cache[i_nuc, i_react] = full_tally_res[i_tally]

        return self._results_cache


# ----------------------------
# Helpers for obtaining energy
# ----------------------------


class ChainFissionHelper(EnergyHelper):
    """Computes energy using fission Q values from depletion chain

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
        super().__init__()
        self._fission_q_vector = None

    def prepare(self, chain_nucs, rate_index, _materials):
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
        _materials : list of str
            Unused. Materials to be tracked for this helper.
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

    def update(self, fission_rates, _mat_index):
        """Update energy produced with fission rates in a material

        Parameters
        ----------
        fission_rates : numpy.ndarray
            fission reaction rate for each isotope in the specified
            material. Should be ordered corresponding to initial
            ``rate_index`` used in :meth:`prepare`
        _mat_index : int
            index for the material requested. Unused, as identical
            isotopes in all materials have the same Q value.
        """
        self._energy += dot(fission_rates, self._fission_q_vector)


class EnergyScoreHelper(EnergyHelper):
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
        for a consistent API across other :class:`EnergyHelper`
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
            self._energy = self._tally.results[0, 0, 1]

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
        operator : openmc.deplete.TransportOperator
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
        Default: 500 [kev]

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
        that occured below ``cutoff``. The number
        of materials in the first axis corresponds
        to the number of materials burned by the
        :class:`openmc.deplete.Operator`
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
        operator : openmc.deplete.Operator
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
            :class:`openmc.deplete.Operator` that uses this helper
            may only burn a subset of all materials when running
            in parallel mode.
        """
        super().generate_tallies(materials, mat_indexes)
        energy_filter = EnergyFilter([0.0, self._cutoff, self._upper_energy])
        self._fission_rate_tally.filters = (
            self._fission_rate_tally.filters + [energy_filter])

    def unpack(self):
        """Obtain fast and thermal fission fractions from tally"""
        if not self._tally_nucs:
            self.results = None
            return
        fission_rates = self._fission_rate_tally.results[..., 1].reshape(
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

    Computes average energy at which fission events occured with

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
            :class:`openmc.deplete.Operator` that uses this helper
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
        if not self._tally_nucs:
            self.results = None
            return
        fission_results = (
            self._fission_rate_tally.results[self._local_indexes, :, 1])
        self.results = (
            self._weighted_tally.results[self._local_indexes, :, 1]).copy()
        nz_mat, nz_nuc = fission_results.nonzero()
        self.results[nz_mat, nz_nuc] /= fission_results[nz_mat, nz_nuc]

    def weighted_yields(self, local_mat_index):
        """Return fission yields for a specific material

        Use the computed average energy of fission
        events to determine fission yields. If average
        energy is between two sets of yields, linearly
        interpolate bewteen the two.
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
        operator : openmc.deplete.TransportOperator
            Operator with a depletion chain
        kwargs :
            Additional keyword arguments to be used in construction

        Returns
        -------
        AveragedFissionYieldHelper
        """
        return cls(operator.chain.nuclides)
