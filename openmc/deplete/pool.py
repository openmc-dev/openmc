"""Dedicated module containing depletion function

Provided to avoid some circular imports
"""
from itertools import repeat, starmap
from multiprocessing import Pool
import functools
from copy import deepcopy
from scipy.sparse import bmat
import numpy as np


# Configurable switch that enables / disables the use of
# multiprocessing routines during depletion
USE_MULTIPROCESSING = True

# Allow user to override the number of worker processes to use for depletion
# calculations
NUM_PROCESSES = None


def deplete(func, chain, x, rates, dt, msr=None, matrix_func=None):
    """Deplete materials using given reaction rates for a specified time

    Parameters
    ----------
    func : callable
        Function to use to get new compositions. Expected to have the
        signature ``func(A, n0, t) -> n1``
    chain : openmc.deplete.Chain
        Depletion chain
    x : list of numpy.ndarray
        Atom number vectors for each material
    rates : openmc.deplete.ReactionRates
        Reaction rates (from transport operator)
    dt : float
        Time in [s] to deplete for
    msr : openmc.deplete.MsrContinuous, Optional
        MsrContinuous removal terms to add to Bateman equation,
        accounting for elements removal and transfer. Creates a sparse matrix
        (scipy.bmat) of matrices (sub-arrays), in which the diagonals are the
        depletion matrices as defined by the Bateman equation with the addition
        of a removal coefficient; and the off-diagonals are the matrices
        accounting for transfers from one material to another.
    maxtrix_func : callable, optional
        Function to form the depletion matrix after calling
        ``matrix_func(chain, rates, fission_yields)``, where
        ``fission_yields = {parent: {product: yield_frac}}``
        Expected to return the depletion matrix required by
        ``func``

    Returns
    -------
    x_result : list of numpy.ndarray
        Updated atom number vectors for each material

    """

    fission_yields = chain.fission_yields
    if len(fission_yields) == 1:
        fission_yields = repeat(fission_yields[0])
    elif len(fission_yields) != len(x):
        raise ValueError(
            "Number of material fission yield distributions {} is not "
            "equal to the number of compositions {}".format(
                len(fission_yields), len(x)))

    if msr is None:
        if matrix_func is None:
            matrices = map(chain.form_matrix, rates, fission_yields)
        else:
            matrices = map(matrix_func, repeat(chain), rates, fission_yields)
        inputs = zip(matrices, x, repeat(dt))
        if USE_MULTIPROCESSING:
            with Pool() as pool:
                x_result = list(pool.starmap(func, inputs))
        else:
            x_result = list(starmap(func, inputs))

    else:
        # Unzip ReactionRates from high orders integrator schemes intermediate
        # steps (such as CF4Integrator)
        if type(rates) is list:
            rate = list(zip(*rates))[0][0]
            list_rates = deepcopy(rates)
        else:
            rate = rates[0]
            list_rates = list(deepcopy(rates))

        #Create a deep copy of fy
        fission_yields = deepcopy(fission_yields)

        # Create zero reaction rates to be added to off-diagonal terms of the
        # sparse matrix, where reaction rates and fy are not needed
        null_rate = rate.copy()
        null_rate.fill(0.0)
        null_fy = deepcopy(fission_yields[0])
        for _, y in null_fy.items():
            y.yields.fill(0.0)

        # Add user-defined removal coefficients to the depletable materials as
        # diagonal elements of the sparse matrix and material coupling transfers
        # as off-diagonal (i,j)
        for j, k in msr.index_msr():
            if j != k:
                if type(rates) is list:
                    list_rates.append((null_rate,)*len(list(zip(*rates))))
                else:
                    list_rates.append(null_rate)
                fission_yields.append(null_fy)
        print(msr.index_msr())
        # Form all the matrices
        if matrix_func is None:
            matrices = map(chain.form_matrix, list_rates, repeat(msr),
                           msr.index_msr(), fission_yields)
        else:
            matrices = map(matrix_func, repeat(chain), list_rates, repeat(msr),
                           msr.index_msr(), fission_yields)

        # Combine all the matrices in a 2d-list of arrays
        matrices = list(matrices)
        raws = []
        for raw in range(msr.n_burn):
            cols = []
            for col in range(msr.n_burn):
                val = None
                for keys, array in zip(msr.index_msr(), matrices):
                    if keys == (raw, col):
                        val = array
                cols.append(val)
            raws.append(cols)

        # Build a sparse matrix from the sparse sub-arrays
        sparse_matrix = bmat(raws)

        # Concatenate all atoms vectors in one
        x = np.concatenate([_x for _x in x])

        x_result = func(sparse_matrix, x, dt)

        # Split the result vectors before returning
        split = np.cumsum([i.shape[0] for i in matrices[:msr.n_burn]])
        x_result = np.split(x_result, split.tolist()[:-1])

    return x_result
