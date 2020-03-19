"""AtomNumber module.

An ndarray to store atom densities with string, integer, or slice indexing.
"""
from collections import OrderedDict

import numpy as np


class AtomNumber:
    """Stores local material compositions (atoms of each nuclide).

    Parameters
    ----------
    local_mats : list of str
        Material IDs
    nuclides : list of str
        Nuclides to be tracked
    volume : dict
        Volume of each material in [cm^3]
    n_nuc_burn : int
        Number of nuclides to be burned.

    Attributes
    ----------
    index_mat : dict
        A dictionary mapping material ID as string to index.
    index_nuc : dict
        A dictionary mapping nuclide name to index.
    volume : numpy.ndarray
        Volume of each material in [cm^3]. If a volume is not found, it defaults
        to 1 so that reading density still works correctly.
    number : numpy.ndarray
        Array storing total atoms for each material/nuclide
    materials : list of str
        Material IDs as strings
    nuclides : list of str
        All nuclide names
    burnable_nuclides : list of str
        Burnable nuclides names. Used for sorting the simulation.
    n_nuc_burn : int
        Number of burnable nuclides.
    n_nuc : int
        Number of nuclides.

    """
    def __init__(self, local_mats, nuclides, volume, n_nuc_burn):
        self.index_mat = OrderedDict((mat, i) for i, mat in enumerate(local_mats))
        self.index_nuc = OrderedDict((nuc, i) for i, nuc in enumerate(nuclides))

        self.volume = np.ones(len(local_mats))
        for mat, val in volume.items():
            if mat in self.index_mat:
                ind = self.index_mat[mat]
                self.volume[ind] = val

        self.n_nuc_burn = n_nuc_burn

        self.number = np.zeros((len(local_mats), len(nuclides)))

    def __getitem__(self, pos):
        """Retrieves total atom number from AtomNumber.

        Parameters
        ----------
        pos : tuple
            A two-length tuple containing a material index and a nuc index.
            These indexes can be strings (which get converted to integers via
            the dictionaries), integers used directly, or slices.

        Returns
        -------
        numpy.ndarray
            The value indexed from self.number.
        """

        mat, nuc = pos
        if isinstance(mat, str):
            mat = self.index_mat[mat]
        if isinstance(nuc, str):
            nuc = self.index_nuc[nuc]

        return self.number[mat, nuc]

    def __setitem__(self, pos, val):
        """Sets total atom number into AtomNumber.

        Parameters
        ----------
        pos : tuple
            A two-length tuple containing a material index and a nuc index.
            These indexes can be strings (which get converted to integers via
            the dictionaries), integers used directly, or slices.
        val : float
            The value [atom] to set the array to.

        """
        mat, nuc = pos
        if isinstance(mat, str):
            mat = self.index_mat[mat]
        if isinstance(nuc, str):
            nuc = self.index_nuc[nuc]

        self.number[mat, nuc] = val

    @property
    def materials(self):
        return self.index_mat.keys()

    @property
    def nuclides(self):
        return self.index_nuc.keys()

    @property
    def n_nuc(self):
        return len(self.index_nuc)

    @property
    def burnable_nuclides(self):
        return [nuc for nuc, ind in self.index_nuc.items()
                if ind < self.n_nuc_burn]

    def get_atom_density(self, mat, nuc):
        """Accesses atom density instead of total number.

        Parameters
        ----------
        mat : str, int or slice
            Material index.
        nuc : str, int or slice
            Nuclide index.

        Returns
        -------
        numpy.ndarray
            Density in [atom/cm^3]

        """
        if isinstance(mat, str):
            mat = self.index_mat[mat]
        if isinstance(nuc, str):
            nuc = self.index_nuc[nuc]

        return self[mat, nuc] / self.volume[mat]

    def set_atom_density(self, mat, nuc, val):
        """Sets atom density instead of total number.

        Parameters
        ----------
        mat : str, int or slice
            Material index.
        nuc : str, int or slice
            Nuclide index.
        val : numpy.ndarray
            Array of densities to set in [atom/cm^3]

        """
        if isinstance(mat, str):
            mat = self.index_mat[mat]
        if isinstance(nuc, str):
            nuc = self.index_nuc[nuc]

        self[mat, nuc] = val * self.volume[mat]

    def get_mat_slice(self, mat):
        """Gets atom quantity indexed by mats for all burned nuclides

        Parameters
        ----------
        mat : str, int or slice
            Material index.

        Returns
        -------
        numpy.ndarray
            The slice requested in [atom].

        """
        if isinstance(mat, str):
            mat = self.index_mat[mat]

        return self[mat, :self.n_nuc_burn]

    def set_mat_slice(self, mat, val):
        """Sets atom quantity indexed by mats for all burned nuclides

        Parameters
        ----------
        mat : str, int or slice
            Material index.
        val : numpy.ndarray
            The slice to set in [atom]

        """
        if isinstance(mat, str):
            mat = self.index_mat[mat]

        self[mat, :self.n_nuc_burn] = val

    def set_density(self, total_density):
        """Sets density.

        Sets the density in the exact same order as total_density_list outputs,
        allowing for internal consistency

        Parameters
        ----------
        total_density : list of numpy.ndarray
            Total atoms.

        """
        for i, density_slice in enumerate(total_density):
            self.set_mat_slice(i, density_slice)
