"""ReactionRates module.

An ndarray to store reaction rates with string, integer, or slice indexing.
"""

import numpy as np


class ReactionRates(np.ndarray):
    """ReactionRates class.

    An ndarray to store reaction rates with string, integer, or slice indexing.

    Parameters
    ----------
    index_mat : OrderedDict of str to int
        A dictionary mapping material ID as string to index.
    index_nuc : OrderedDict of str to int
        A dictionary mapping nuclide name as string to index.
    index_rx : OrderedDict of str to int
        A dictionary mapping reaction name as string to index.

    Attributes
    ----------
    index_mat : OrderedDict of str to int
        A dictionary mapping material ID as string to index.
    index_nuc : OrderedDict of str to int
        A dictionary mapping nuclide name as string to index.
    index_rx : OrderedDict of str to int
        A dictionary mapping reaction name as string to index.
    n_mat : int
        Number of materials.
    n_nuc : int
        Number of nucs.
    n_react : int
        Number of reactions.

    """
    def __new__(cls, index_mat, index_nuc, index_rx):
        # Create appropriately-sized zeroed-out ndarray
        shape = (len(index_mat), len(index_nuc), len(index_rx))
        obj = super().__new__(cls, shape)
        obj[:] = 0.0

        # Add mapping attributes
        obj.index_mat = index_mat
        obj.index_nuc = index_nuc
        obj.index_rx = index_rx

        return obj

    def __array_finalize__(self, obj):
        if obj is None:
            return
        self.index_mat = getattr(obj, 'index_mat', None)
        self.index_nuc = getattr(obj, 'index_nuc', None)
        self.index_rx = getattr(obj, 'index_rx', None)

    def __reduce__(self):
        state = super().__reduce__()
        new_state = state[2] + (self.index_mat, self.index_nuc, self.index_rx)
        return (state[0], state[1], new_state)

    def __setstate__(self, state):
        self.index_mat = state[-3]
        self.index_nuc = state[-2]
        self.index_rx = state[-1]
        super().__setstate__(state[0:-3])

    @property
    def n_mat(self):
        """Number of cells."""
        return len(self.index_mat)

    @property
    def n_nuc(self):
        """Number of nucs."""
        return len(self.index_nuc)

    @property
    def n_react(self):
        """Number of reactions."""
        return len(self.index_rx)

    def get(self, mat, nuc, rx):
        """Get reaction rate by material/nuclide/reaction

        Parameters
        ----------
        mat : str
            Material ID as a string
        nuc : str
            Nuclide name
        rx : str
            Name of the reaction

        Returns
        -------
        float
            Reaction rate corresponding to given material, nuclide, and reaction

        """
        mat = self.index_mat[mat]
        nuc = self.index_nuc[nuc]
        rx = self.index_rx[rx]
        return self[mat, nuc, rx]

    def set(self, mat, nuc, rx, value):
        """Set reaction rate by material/nuclide/reaction

        Parameters
        ----------
        mat : str
            Material ID as a string
        nuc : str
            Nuclide name
        rx : str
            Name of the reaction
        value : float
            Corresponding reaction rate to set

        """
        mat = self.index_mat[mat]
        nuc = self.index_nuc[nuc]
        rx = self.index_rx[rx]
        self[mat, nuc, rx] = value
