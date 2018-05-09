import h5py
import numpy as np

from .results import Results, _VERSION_RESULTS
from openmc.checkvalue import check_filetype_version


class ResultsList(list):
    """A list of openmc.deplete.Results objects

    Parameters
    ----------
    filename : str
        The filename to read from.

    """
    def __init__(self, filename):
        super().__init__()
        with h5py.File(str(filename), "r") as fh:
            check_filetype_version(fh, 'depletion results', _VERSION_RESULTS[0])

            # Get number of results stored
            n = fh["number"].value.shape[0]

            for i in range(n):
                self.append(Results.from_hdf5(fh, i))

    def get_atoms(self, mat, nuc):
        """Get nuclide concentration over time from a single material

        Parameters
        ----------
        mat : str
            Material name to evaluate
        nuc : str
            Nuclide name to evaluate

        Returns
        -------
        time : numpy.ndarray
            Array of times in [s]
        concentration : numpy.ndarray
            Total number of atoms for specified nuclide

        """
        time = np.empty_like(self)
        concentration = np.empty_like(self)

        # Evaluate value in each region
        for i, result in enumerate(self):
            time[i] = result.time[0]
            concentration[i] = result[0, mat, nuc]

        return time, concentration

    def get_reaction_rate(self, mat, nuc, rx):
        """Get reaction rate in a single material/nuclide over time

        Parameters
        ----------
        mat : str
            Material name to evaluate
        nuc : str
            Nuclide name to evaluate
        rx : str
            Reaction rate to evaluate

        Returns
        -------
        time : numpy.ndarray
            Array of times in [s]
        rate : numpy.ndarray
            Array of reaction rates

        """
        time = np.empty_like(self)
        rate = np.empty_like(self)

        # Evaluate value in each region
        for i, result in enumerate(self):
            time[i] = result.time[0]
            rate[i] = result.rates[0].get(mat, nuc, rx) * result[0, mat, nuc]

        return time, rate

    def get_eigenvalue(self):
        """Evaluates the eigenvalue from a results list.

        Returns
        -------
        time : numpy.ndarray
            Array of times in [s]
        eigenvalue : numpy.ndarray
            k-eigenvalue at each time

        """
        time = np.empty_like(self)
        eigenvalue = np.empty_like(self)

        # Get time/eigenvalue at each point
        for i, result in enumerate(self):
            time[i] = result.time[0]
            eigenvalue[i] = result.k[0]

        return time, eigenvalue
