"""
Class for normalizing fission energy deposition
"""
from collections import defaultdict
from itertools import product

from numpy import dot, zeros, newaxis, divide, asarray

from openmc.capi import (
    Tally, MaterialFilter, EnergyFilter)
from .abc import ReactionRateHelper, EnergyHelper

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

        Uses the :mod:`openmc.capi` to generate a tally
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


# ------------------------------------
# Helper for collapsing fission yields
# ------------------------------------


class FissionYieldHelper(object):
    """Class for using energy-dependent fission yields in depletion chain

    Creates a tally across all burnable materials to score the fission
    rate in nuclides with yield data. An energy filter is used to
    compute this rates in a group structure corresponding to the
    fission yield data. This tally data is used to compute the
    relative number of fission events in each energy region,
    which serve as the weights for each energy-dependent fission
    yield distribution.

    Parameters
    ----------
    chain_nuclides : iterable of openmc.deplete.nuclide
        Nuclides tracked in the depletion chain. Not necessary
        that all have yield data.
    n_bmats : int
        Number of burnable materials tracked in the problem

    Attributes
    ----------
    energy_bounds : tuple of float
        Sorted energy bounds from the tally filter
    results : numpy.ndarray
        Array of tally results for this process with shape
        ``(n_local_mat, n_energy, n_nucs)``
    libraries : list of dict
        List of fission yield dictionaries of the form
        ``{parent: {product: yield}}``. Populated in
        :meth:`compute_yields` and reset during
        :meth:`unpack`.
    """

    def __init__(self, chain_nuclides, n_bmats):
        self.libraries = []
        self._chain_nuclides = {}
        self._chain_set = set()
        self._tally_map = {}
        # TODO Support user-requested minimum energy?
        yield_energies = {0.0}

        # Get nuclides with fission yield data, names
        # and all energy points
        # Names are provided from operator tally nuclides
        for nuc in chain_nuclides:
            if len(nuc.yield_data) == 0:
                continue
            self._chain_nuclides[nuc.name] = nuc
            self._chain_set.add(nuc.name)
            yield_energies.update(nuc.yield_energies)

        # Create energy grid
        self._energy_bounds = tuple(sorted(yield_energies))
        self.n_bmats = n_bmats

        self._reaction_tally = None
        self.results = None
        self.local_indexes = None

    @property
    def energy_bounds(self):
        return self._energy_bounds

    def generate_tallies(self, materials, mat_indexes):
        """Construct the fission rate tally

        Parameters
        ----------
        materials : iterable of C-API materials
            Materials to be used in :class:`openmc.capi.MaterialFilter`
        mat_indexes : iterable of int
            Indexes for materials in ``materials`` tracked on this
            process
        """
        # Tally group-wise fission reaction rates
        self._reaction_tally = Tally()
        self._reaction_tally.scores = ['fission']

        # Tally energy-weighted group-wise fission reaction rate
        # Used to evaluated linear interpolation between fission yield points
        self._weighted_reaction_tally = Tally()
        self._weighted_reaction_tally.scores = ['fission']

        filters = [
            MaterialFilter(materials), EnergyFilter(self._energy_bounds)]

        self._reaction_tally.filters = filters
        self.local_indexes = asarray(mat_indexes)

    def set_fissionable_nuclides(self, nuclides):
        """List of string of nuclides with data to be tallied

        Parameters
        ----------
        nuclides : iterable of str
            nuclides with non-zero densities that are candidates
            for the fission tally. Not necessary that all are nuclides
            with fission yields, but at least one must be

        Returns
        -------
        nuclides : tuple of str
            Nuclides ordered as they appear in the tally and in
            the nuclide column of :attr:`results`

        Raises
        ------
        ValueError
            If no nuclides in ``nuclides`` are tracked on this
            object
        """
        # Set of all nuclides with positive density
        # and fission yield data
        nuc_set = self._chain_set & set(nuclides)
        if len(nuc_set) == 0:
            raise ValueError(
                "No overlap between chain nuclides with fission yields and "
                "requested tally nuclides")
        nuclides = tuple(sorted(nuc_set))
        self._tally_index = [self._chain_nuclides[n] for n in nuclides]
        self._reaction_tally.nuclides = nuclides
        return nuclides

    def unpack(self):
        """Unpack fission rate tallies to produce :attr:`results`

        Resets :attr:`libraries` under the assumption this is called
        during the :class:`openmc.deplete.Operator` unpackign process
        """
        # clear old libraries
        self.libraries = []

        # get view into tally results
        # new shape: [material, energy, parent nuclide]
        result_view = self._reaction_tally.results[..., 1].reshape(
            self.n_bmats, len(self._energy_bounds) - 1,
            len(self._reaction_tally.nuclides))

        # Get results specific to this process
        self.results = result_view[self.local_indexes, ...]

        # scale fission yields proportional to total fission rate
        fsn_rate = self.results.sum(axis=1)
        # TODO Guard against divide by zero
        self.results /= fsn_rate[:, newaxis, :]

    def compute_yields(self, local_mat_index):
        """Compute single fission yields using :attr:`results`

        Produces a new library in :attr:`self.libraries`

        Parameters
        ----------
        local_mat_index : int
            Index for material tracked on this process that
            exists in :attr:`local_mat_index` and fits within
            the first axis in :attr:`results`
        """
        tally_results = self.results[local_mat_index]

        # Dictionary {parent_nuclide : [product, yield_vector]}
        initial_library = {}
        for i_energy, energy in enumerate(self._energy_bounds[1:]):
            for i_nuc, fiss_frac in enumerate(tally_results[i_energy]):
                parent = self._tally_index[i_nuc]
                yield_data = parent.yield_data.get(energy)
                if yield_data is None:
                    continue
                if parent not in initial_library:
                    initial_library[parent] = yield_data * fiss_frac
                    continue
                initial_library[parent] += yield_data * fiss_frac

        # convert to dictionary that can be passed to Chain.form_matrix
        # {parent: {product: yield}}
        library = {}
        for k, yield_obj in initial_library.items():
            library[k.name] = dict(zip(yield_obj.products, yield_obj.yields))

        self.libraries.append(library)
