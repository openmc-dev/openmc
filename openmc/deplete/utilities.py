"""The utilities module.

Contains functions that can be used to post-process objects that come out of
the results module.
"""

import numpy as np


def evaluate_single_nuclide(results, mat, nuc):
    """Evaluates a single nuclide in a single material from a results list.

    Parameters
    ----------
    results : list of results
        The results to extract data from.  Must be sorted and continuous.
    mat : str
        Material name to evaluate
    nuc : str
        Nuclide name to evaluate

    Returns
    -------
    time : numpy.ndarray
        Time vector
    concentration : numpy.ndarray
        Total number of atoms in the material

    """
    n_points = len(results)
    time = np.zeros(n_points)
    concentration = np.zeros(n_points)

    # Evaluate value in each region
    for i, result in enumerate(results):
        time[i] = result.time[0]
        concentration[i] = result[0, mat, nuc]

    return time, concentration

def evaluate_reaction_rate(results, mat, nuc, rx):
    """Return reaction rate in a single material/nuclide from a results list.

    Parameters
    ----------
    results : list of openmc.deplete.Results
        The results to extract data from.  Must be sorted and continuous.
    mat : str
        Material name to evaluate
    nuc : str
        Nuclide name to evaluate
    rx : str
        Reaction rate to evaluate

    Returns
    -------
    time : numpy.ndarray
        Time vector.
    rate : numpy.ndarray
        Reaction rate.

    """
    n_points = len(results)
    time = np.zeros(n_points)
    rate = np.zeros(n_points)
    # Evaluate value in each region
    for i, result in enumerate(results):
        time[i] = result.time[0]
        rate[i] = result.rates[0].get(mat, nuc, rx) * result[0, mat, nuc]

    return time, rate


def evaluate_eigenvalue(results):
    """Evaluates the eigenvalue from a results list.

    Parameters
    ----------
    results : list of openmc.deplete.Results
        The results to extract data from.  Must be sorted and continuous.

    Returns
    -------
    time : numpy.ndarray
        Time vector.
    eigenvalue : numpy.ndarray
        Eigenvalue.

    """
    n_points = len(results)
    time = np.zeros(n_points)
    eigenvalue = np.zeros(n_points)

    # Evaluate value in each region
    for i, result in enumerate(results):

        time[i] = result.time[0]
        eigenvalue[i] = result.k[0]

    return time, eigenvalue
