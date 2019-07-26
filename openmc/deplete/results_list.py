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
            n = fh["number"][...].shape[0]

            for i in range(n):
                self.append(Results.from_hdf5(fh, i))

    def get_atoms(self, mat, nuc):
        """Get number of nuclides over time from a single material

        .. note::

            Initial values for some isotopes that do not appear in
            initial concentrations may be non-zero, depending on the
            value of :class:`openmc.deplete.Operator` ``dilute_initial``.
            The :class:`openmc.deplete.Operator` adds isotopes according
            to this setting, which can be set to zero.

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
        time = np.empty_like(self, dtype=float)
        concentration = np.empty_like(self, dtype=float)

        # Evaluate value in each region
        for i, result in enumerate(self):
            time[i] = result.time[0]
            concentration[i] = result[0, mat, nuc]

        return time, concentration

    def get_reaction_rate(self, mat, nuc, rx):
        """Get reaction rate in a single material/nuclide over time

        .. note::

            Initial values for some isotopes that do not appear in
            initial concentrations may be non-zero, depending on the
            value of :class:`openmc.deplete.Operator` ``dilute_initial``
            The :class:`openmc.deplete.Operator` adds isotopes according
            to this setting, which can be set to zero.

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
        time = np.empty_like(self, dtype=float)
        rate = np.empty_like(self, dtype=float)

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
            k-eigenvalue at each time. Column 0
            contains the eigenvalue, while column
            1 contains the associated uncertainty

        """
        time = np.empty_like(self, dtype=float)
        eigenvalue = np.empty((len(self), 2), dtype=float)

        # Get time/eigenvalue at each point
        for i, result in enumerate(self):
            time[i] = result.time[0]
            eigenvalue[i] = result.k[0]

        return time, eigenvalue

    def get_depletion_time(self):
        """Return an array of the average time to deplete a material

        ..note::

            Will have one fewer row than number of other methods,
            like :meth:`get_eigenvalues`, because no depletion
            is performed at the final transport stage

        Returns
        -------

        times : :class:`numpy.ndarray`
            Vector of average time to deplete a single material
            across all processes and materials.

        """
        times = np.empty(len(self) - 1)
        # Need special logic because the predictor
        # writes EOS values for step i as BOS values
        # for step i+1
        # The first proc_time may be zero
        if self[0].proc_time > 0.0:
            items = self[:-1]
        else:
            items = self[1:]
        for ix, res in enumerate(items):
            times[ix] = res.proc_time
        return times
