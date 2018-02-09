"""The utilities module.

Contains functions that can be used to post-process objects that come out of
the results module.
"""

import numpy as np


def evaluate_single_nuclide(results, cell, nuc):
    """Evaluates a single nuclide in a single cell from a results list.

    Parameters
    ----------
    results : list of results
        The results to extract data from.  Must be sorted and continuous.
    cell : str
        Cell name to evaluate
    nuc : str
        Nuclide name to evaluate

    Returns
    -------
    time : numpy.array
        Time vector.
    concentration : numpy.array
        Total number of atoms in the cell.
    """

    n_points = len(results)
    time = np.zeros(n_points)
    concentration = np.zeros(n_points)

    # Evaluate value in each region
    for i, result in enumerate(results):
        time[i] = result.time[0]
        concentration[i] = result[0, cell, nuc]

    return time, concentration

def evaluate_reaction_rate(results, cell, nuc, rxn):
    """Evaluates a single nuclide reaction rate in a single cell from a results list.

    Parameters
    ----------
    results : list of Results
        The results to extract data from.  Must be sorted and continuous.
    cell : str
        Cell name to evaluate
    nuc : str
        Nuclide name to evaluate
    rxn : str
        Reaction rate to evaluate

    Returns
    -------
    time : numpy.array
        Time vector.
    rate : numpy.array
        Reaction rate.
    """

    n_points = len(results)
    time = np.zeros(n_points)
    rate = np.zeros(n_points)
    # Evaluate value in each region
    for i, result in enumerate(results):
        time[i] = result.time[0]
        rate[i] = result.rates[0][cell, nuc, rxn] * result[0, cell, nuc]

    return time, rate


def evaluate_eigenvalue(results):
    """Evaluates the eigenvalue from a results list.

    Parameters
    ----------
    results : list of Results
        The results to extract data from.  Must be sorted and continuous.

    Returns
    -------
    time : numpy.array
        Time vector.
    eigenvalue : numpy.array
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
