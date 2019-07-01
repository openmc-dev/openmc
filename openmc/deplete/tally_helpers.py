"""
Class for normalizing fission energy deposition
"""

from numpy import dot, zeros

from openmc.capi import Tally, MaterialFilter


class ChainFissTallyHelper(object):
    """Fission Q-values are pulled from chain"""

    def __init__(self, n_materials):
        self._fiss_q = None
        self.n_materials = n_materials
        self._rx_tally = None

    def set_fission_q(self, chain_nucs, rate_index):
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

    def get_fiss_energy(self, fiss_rates, _mat_index):
        return dot(fiss_rates, self._fiss_q)
