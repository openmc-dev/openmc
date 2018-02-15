"""ReactionRates module.

An ndarray to store reaction rates with string, integer, or slice indexing.
"""

import numpy as np


class ReactionRates(object):
    """ReactionRates class.

    An ndarray to store reaction rates with string, integer, or slice indexing.

    Parameters
    ----------
    mat_to_ind : OrderedDict of str to int
        A dictionary mapping material ID as string to index.
    nuc_to_ind : OrderedDict of str to int
        A dictionary mapping nuclide name as string to index.
    react_to_ind : OrderedDict of str to int
        A dictionary mapping reaction name as string to index.

    Attributes
    ----------
    mat_to_ind : OrderedDict of str to int
        A dictionary mapping material ID as string to index.
    nuc_to_ind : OrderedDict of str to int
        A dictionary mapping nuclide name as string to index.
    react_to_ind : OrderedDict of str to int
        A dictionary mapping reaction name as string to index.
    n_mat : int
        Number of materials.
    n_nuc : int
        Number of nucs.
    n_react : int
        Number of reactions.
    rates : numpy.array
        Array storing rates indexed by the above dictionaries.
    """

    def __init__(self, mat_to_ind, nuc_to_ind, react_to_ind):

        self.mat_to_ind = mat_to_ind
        self.nuc_to_ind = nuc_to_ind
        self.react_to_ind = react_to_ind

        self.rates = np.zeros((self.n_mat, self.n_nuc, self.n_react))

    def __getitem__(self, pos):
        """Retrieves an item from reaction_rates.

        Parameters
        ----------
        pos : tuple
            A three-length tuple containing a material index, a nuc index, and a
            reaction index.  These indexes can be strings (which get converted
            to integers via the dictionaries), integers used directly, or
            slices.

        Returns
        -------
        numpy.array
            The value indexed from self.rates.
        """

        mat, nuc, react = pos
        if isinstance(mat, str):
            mat = self.mat_to_ind[mat]
        if isinstance(nuc, str):
            nuc = self.nuc_to_ind[nuc]
        if isinstance(react, str):
            react = self.react_to_ind[react]

        return self.rates[mat, nuc, react]

    def __setitem__(self, pos, val):
        """Sets an item from reaction_rates.

        Parameters
        ----------
        pos : tuple
            A three-length tuple containing a material index, a nuc index, and a
            reaction index.  These indexes can be strings (which get converted
            to integers via the dictionaries), integers used directly, or
            slices.
        val : float
            The value to set the array to.
        """

        mat, nuc, react = pos
        if isinstance(mat, str):
            mat = self.mat_to_ind[mat]
        if isinstance(nuc, str):
            nuc = self.nuc_to_ind[nuc]
        if isinstance(react, str):
            react = self.react_to_ind[react]

        self.rates[mat, nuc, react] = val

    @property
    def n_mat(self):
        """Number of cells."""
        return len(self.mat_to_ind)

    @property
    def n_nuc(self):
        """Number of nucs."""
        return len(self.nuc_to_ind)

    @property
    def n_react(self):
        """Number of reactions."""
        return len(self.react_to_ind)
