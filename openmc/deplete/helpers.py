"""
Class for normalizing fission energy deposition
"""
from itertools import product

from numpy import dot, zeros

from openmc.capi import Tally, MaterialFilter
from .abc import ReactionRateHelper, FissionEnergyHelper

# -------------------------------------
# Helpers for generating reaction rates
# -------------------------------------


class DirectReactionRateHelper(ReactionRateHelper):
    """Class that generates tallies for one-group rates"""

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
            Unique id for the requested material
        nuc_index : iterable of int
            Index for each nuclide in :attr:`nuclides` in the
            desired reaction rate matrix
        react_index : iterable of int
            Index for each reaction scored in the tally

        Returns
        -------
        rates : :class:`numpy.ndarray`
            2D matrix ``(len(nuc_index), len(react_index))`` with the
            reaction rates in this material
        """
        results = self._reset_results_cache(len(nuc_index), len(react_index))
        full_tally_res = self._rate_tally.results[mat_id, :, 1]
        for i_tally, (i_nuc, i_react) in enumerate(
                product(nuc_index, react_index)):
            results[i_nuc, i_react] = full_tally_res[i_tally]

        return results


# ------------------------------------
# Helpers for obtaining fission energy
# ------------------------------------


class ChainFissHelper(FissionEnergyHelper):
    """Fission Q-values are pulled from chain"""

    def prepare(self, chain_nucs, rate_index, _materials):
        """Populate the fission Q value vector from a chain.

        Paramters
        ---------
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
        if (self._fission_E is not None
                and self._fission_E.shape == (len(rate_index),)):
            return

        fiss_E = zeros(len(rate_index))

        for nuclide in chain_nucs:
            if nuclide.name in rate_index:
                for rx in nuclide.reactions:
                    if rx.type == "fission":
                        fiss_E[rate_index[nuclide.name]] = rx.Q
                        break

        self._fission_E = fiss_E

    def get_fission_energy(self, fiss_rates, _mat_index):
        """Return a vector of the isotopic fission energy for this material

        parameters
        ----------
        fission_rates : numpy.ndarray
            fission reaction rate for each isotope in the specified
            material. should be ordered corresponding to initial
            ``rate_index`` used in :meth:`set_fission_q`
        _mat_index : int
            index for the material requested. Unused, as all
            isotopes in all materials have the same Q value.
        """
        return dot(fiss_rates, self._fission_E)
