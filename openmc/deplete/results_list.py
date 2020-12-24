import numbers
import bisect
import math

import h5py
import numpy as np

from .results import Results, VERSION_RESULTS
from openmc.checkvalue import check_filetype_version, check_value, check_type


__all__ = ["ResultsList"]


class ResultsList(list):
    """A list of openmc.deplete.Results objects

    It is recommended to use :meth:`from_hdf5` over
    direct creation.
    """

    @classmethod
    def from_hdf5(cls, filename):
        """Load in depletion results from a previous file

        Parameters
        ----------
        filename : str
            Path to depletion result file

        Returns
        -------
        new : ResultsList
            New instance of depletion results
        """
        with h5py.File(str(filename), "r") as fh:
            check_filetype_version(fh, 'depletion results', VERSION_RESULTS[0])
            new = cls()

            # Get number of results stored
            n = fh["number"][...].shape[0]

            for i in range(n):
                new.append(Results.from_hdf5(fh, i))
        return new

    def get_atoms(self, mat, nuc, nuc_units="atoms", time_units="s"):
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
        nuc_units : {"atoms", "atom/b-cm", "atom/cm3"}, optional
            Units for the returned concentration. Default is ``"atoms"``

            .. versionadded:: 0.12
        time_units : {"s", "min", "h", "d"}, optional
            Units for the returned time array. Default is ``"s"`` to
            return the value in seconds.

            .. versionadded:: 0.12

        Returns
        -------
        times : numpy.ndarray
            Array of times in units of ``time_units``
        concentrations : numpy.ndarray
            Concentration of specified nuclide in units of ``nuc_units``

        """
        check_value("time_units", time_units, {"s", "d", "min", "h"})
        check_value("nuc_units", nuc_units,
                    {"atoms", "atom/b-cm", "atom/cm3"})

        times = np.empty_like(self, dtype=float)
        concentrations = np.empty_like(self, dtype=float)

        # Evaluate value in each region
        for i, result in enumerate(self):
            times[i] = result.time[0]
            concentrations[i] = result[0, mat, nuc]

        # Unit conversions
        if time_units == "d":
            times /= (60 * 60 * 24)
        elif time_units == "h":
            times /= (60 * 60)
        elif time_units == "min":
            times /= 60

        if nuc_units != "atoms":
            # Divide by volume to get density
            concentrations /= self[0].volume[mat]
            if nuc_units == "atom/b-cm":
                # 1 barn = 1e-24 cm^2
                concentrations *= 1e-24

        return times, concentrations

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
        times : numpy.ndarray
            Array of times in [s]
        rates : numpy.ndarray
            Array of reaction rates

        """
        times = np.empty_like(self, dtype=float)
        rates = np.empty_like(self, dtype=float)

        # Evaluate value in each region
        for i, result in enumerate(self):
            times[i] = result.time[0]
            rates[i] = result.rates[0].get(mat, nuc, rx) * result[0, mat, nuc]

        return times, rates

    def get_eigenvalue(self):
        """Evaluates the eigenvalue from a results list.

        Returns
        -------
        times : numpy.ndarray
            Array of times in [s]
        eigenvalues : numpy.ndarray
            k-eigenvalue at each time. Column 0
            contains the eigenvalue, while column
            1 contains the associated uncertainty

        """
        times = np.empty_like(self, dtype=float)
        eigenvalues = np.empty((len(self), 2), dtype=float)

        # Get time/eigenvalue at each point
        for i, result in enumerate(self):
            times[i] = result.time[0]
            eigenvalues[i] = result.k[0]

        return times, eigenvalues

    def get_depletion_time(self):
        """Return an array of the average time to deplete a material

        .. note::

            Will have one fewer row than number of other methods,
            like :meth:`get_eigenvalues`, because no depletion
            is performed at the final transport stage

        Returns
        -------
        times : numpy.ndarray
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

    def get_times(self, time_units="d") -> np.ndarray:
        """Return the points in time that define the depletion schedule


        .. versionadded:: 0.12.1

        Parameters
        ----------
        time_units : {"s", "d", "h", "min"}, optional
            Return the vector in these units. Default is to
            convert to days

        Returns
        -------
        numpy.ndarray
            1-D vector of time points

        """
        check_type("time_units", time_units, str)

        times = np.fromiter(
            (r.time[0] for r in self),
            dtype=self[0].time.dtype,
            count=len(self),
        )

        if time_units == "d":
            times /= (60 * 60 * 24)
        elif time_units == "h":
            times /= (60 * 60)
        elif time_units == "min":
            times /= 60
        elif time_units != "s":
            raise ValueError(
                'Unable to set "time_units" to {} since it is not '
                'in ("s", "d", "min", "h")'.format(time_units)
            )
        return times

    def get_step_where(
        self, time, time_units="d", atol=1e-6, rtol=1e-3
    ) -> int:
        """Return the index closest to a given point in time

        In the event ``time`` lies exactly between two points, the
        lower index will be returned. It is possible that the index
        will be at most one past the point in time requested, but only
        according to tolerances requested.

        Passing ``atol=math.inf`` and ``rtol=math.inf`` will return
        the closest index to the requested point.


        .. versionadded:: 0.12.1

        Parameters
        ----------
        time : float
            Desired point in time
        time_units : {"s", "d", "min", "h"}, optional
            Units on ``time``. Default: days
        atol : float, optional
            Absolute tolerance (in ``time_units``) if ``time`` is not
            found.
        rtol : float, optional
            Relative tolerance if ``time`` is not found.

        Returns
        -------
        int

        """
        check_type("time", time, numbers.Real)
        check_type("atol", atol, numbers.Real)
        check_type("rtol", rtol, numbers.Real)

        times = self.get_times(time_units)

        if times[0] < time < times[-1]:
            ix = bisect.bisect_left(times, time)
            if ix == times.size:
                ix -= 1
            # Bisection will place us either directly on the point
            # or one-past the first value less than time
            elif time - times[ix - 1] <= times[ix] - time:
                ix -= 1
        elif times[0] >= time:
            ix = 0
        elif time >= times[-1]:
            ix = times.size - 1

        if math.isclose(time, times[ix], rel_tol=rtol, abs_tol=atol):
            return ix

        raise ValueError(
            "A value of {} {} was not found given absolute and "
            "relative tolerances {} and {}.".format(
                time, time_units, atol, rtol)
        )
