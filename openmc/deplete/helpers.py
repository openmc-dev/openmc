"""
Class for normalizing fission energy deposition
"""
from itertools import product
from numbers import Real

from numpy import dot, zeros, newaxis, asarray, empty_like, where

from openmc.checkvalue import check_type, check_greater_than
from openmc.capi import (
    Tally, MaterialFilter, EnergyFilter)
from .abc import (
    ReactionRateHelper, EnergyHelper, FissionYieldHelper)

__all__ = (
    "DirectReactionRateHelper", "ChainFissionHelper",
    "ConstantFissionYieldHelper")

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


class ConstantFissionYieldHelper(FissionYieldHelper):
    """Class that uses a single set of fission yields on each isotope

    Parameters
    ----------
    chain_nuclides : iterable of openmc.deplete.Nuclide
        Nuclides tracked in the depletion chain. Not necessary
        that all have yield data.
    energy : float, optional
        Key in :attr:`openmc.deplete.Nuclide.yield_data` corresponding
        to the desired set of fission yield data. Typically one of
        ``{0.0253, 500000, 14000000}`` corresponding to 0.0253 eV,
        500 keV, and 14 MeV yield libraries. If the specific key is not
        found, will fall back to closest energy present.
        Default: 0.0253 eV for thermal yields

    Attributes
    ----------
    constant_yields : dict of str to :class:`openmc.deplete.FissionYield`
        Fission yields for all nuclides that only have one set of
        fission yield data. Can be accessed as ``{parent: {product: yield}}``
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
            min_index = min(
                range(len(nuc.yield_energies)), key=distances.__getitem__)
            self._constant_yields[name] = (
                nuc.yield_data[nuc.yield_energies[min_index]])

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
        library : dict
            Dictionary of ``{parent: {product: fyield}}``
        """
        return self.constant_yields

    def compute_yields(self, local_mat_index):
        """Compute single fission yields using :attr:`results`

        Produces a new library in :attr:`libraries`

        Parameters
        ----------
        local_mat_index : int
            Index for material tracked on this process that
            exists in :attr:`local_mat_index` and fits within
            the first axis in :attr:`results`

        Returns
        -------
        library : dict
            Dictionary of ``{parent: {product: fyield}}``
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
        library = {}
        for k, yield_obj in initial_library.items():
            library[k.name] = dict(zip(yield_obj.products, yield_obj.yields))

        return library
