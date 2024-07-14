"""D1S module

This module contains functionality to support the direct 1-step (D1S) method for
shutdown dose rate calculations.

"""

from typing import Sequence
from math import log

import numpy as np

from openmc.data import half_life
from .abc import _normalize_timesteps


def time_correction_factors(
        nuclides: list[str],
        timesteps: Sequence[float] | Sequence[tuple[float, str]],
        source_rates: float | Sequence[float],
        timestep_units: str = 's'
) -> dict[str, np.ndarray]:
    """Calculate time correction factors for the D1S method.

    This function determines the time correction factor that should be applied
    to photon tallies as part of the D1S method.

    Parameters
    ----------
    nuclide : str
        The name of the nuclide to find the time correction for, e.g., 'Ni65'
    timesteps : iterable of float or iterable of tuple
        Array of timesteps. Note that values are not cumulative. The units are
        specified by the `timestep_units` argument when `timesteps` is an
        iterable of float. Alternatively, units can be specified for each step
        by passing a sequence of (value, unit) tuples.
    source_rates : float or iterable of float, optional
        Source rate in [neutron/sec] for each interval in :attr:`timesteps`
    timestep_units : {'s', 'min', 'h', 'd', 'a'}
        Units for values specified in the `timesteps` argument. 's' means
        seconds, 'min' means minutes, 'h' means hours, and 'a' means Julian
        years.


    Returns
    -------
    Dictionary mapping nuclide to an array of time correction factors for each
    time.

    """

    # Determine normalized timesteps and source rates
    timesteps, source_rates = _normalize_timesteps(
        timesteps, source_rates, timestep_units)

    # Calculate decay rate for each nuclide
    decay_rate = np.array([log(2.0) / half_life(x) for x in nuclides])

    n_timesteps = len(timesteps) + 1
    n_nuclides = len(nuclides)

    # Create a 2D array for the time correction factors
    h = np.zeros((n_timesteps, n_nuclides))

    for i, (dt, rate) in enumerate(zip(timesteps, source_rates)):
        # Precompute the exponential term
        g = np.exp(-decay_rate*dt)

        # Eq. (4) in doi:10.1016/j.fusengdes.2019.111399
        h[i + 1] = rate*(1. - g) + h[i]*g

    return {nuclides[i]: h[:, i] for i in range(n_nuclides)}
