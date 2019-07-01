"""
Class for normalizing fission energy deposition
"""
from abc import ABC, abstractmethod

from numpy import dot, zeros

from openmc.checkvalue import check_type
from openmc.capi import Tally, MaterialFilter


class TallyHelperBase(ABC):
    """Base class for working with tallies for depletion"""

    def __init__(self):
        self._fiss_q = None
        self._rx_tally = None
        self._nuclides = []

    @abstractmethod
    def set_fission_q(self, chain_nucs, rate_index):
        """Populate the energy released per fission Q value array"""
        pass

    @property
    def reaction_tally(self):
        if self._rx_tally is None:
            raise AttributeError(
                "Reaction tally for {} not set.".format(
                    self.__class__.__name__
                )
            )
        return self._rx_tally

    def generate_tallies(self, materials, scores):
        self._rx_tally = Tally()
        self._rx_tally.scores = scores
        self._rx_tally.filters = [MaterialFilter(materials)]

    @property
    def nuclides(self):
        """List of nuclides with requested reaction rates"""
        return self._nuclides

    @nuclides.setter
    def nuclides(self, nuclides):
        check_type("nuclides", nuclides, list, str)
        self._nuclides = nuclides
        self._rx_tally.nuclides = nuclides

    @abstractmethod
    def get_fission_energy(self, fission_rates, mat_index):
        """return a vector of the isotopic fission energy for this material

        parameters
        ----------
        fission_rates: numpy.ndarray
            fission reaction rate for each isotope in the specified
            material. should be ordered corresponding to initial
            ``rate_index`` used in :meth:`set_fission_q`
        mat_index: int
            index for the material requested.
        """


class ChainFissTallyHelper(TallyHelperBase):
    """Fission Q-values are pulled from chain"""

    def set_fission_q(self, chain_nucs, rate_index):
        """Populate the fission Q value vector from a chain.

        Paramters
        ---------
        chain_nucs: iterable of :class:`openmc.deplete.Nuclide`
            Nuclides used in this depletion chain. Do not need
            to be ordered
        rate_index: dict of str to int
            Dictionary mapping names of nuclides, e.g. ``"U235"``,
            to a corresponding index in the desired fission Q
            vector.
        """
        if (self._fiss_q is not None
                and self._fiss_q.shape == (len(rate_index), )):
            return

        fq = zeros(len(rate_index))

        for nuclide in chain_nucs:
            if nuclide.name in rate_index:
                for rx in nuclide.reactions:
                    if rx.type == "fission":
                        fq[rate_index[nuclide.name]] = rx.Q
                        break

        self._fiss_q = fq

    def get_fission_energy(self, fiss_rates, _mat_index):
        """Return a vector of the isotopic fission energy for this material

        parameters
        ----------
        fission_rates: numpy.ndarray
            fission reaction rate for each isotope in the specified
            material. should be ordered corresponding to initial
            ``rate_index`` used in :meth:`set_fission_q`
        mat_index: int
            index for the material requested. Unused, as all
            isotopes in all materials have the same Q value.
        """
        return dot(fiss_rates, self._fiss_q)
