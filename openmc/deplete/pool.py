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
        #Build first diagonal matrices, with removal term
        if matrix_func is None:
            diag_matrices = map(chain.form_matrix, rates, fission_yields,
                                repeat(msr), msr.index_msr)
        else:
            diag_matrices = map(matrix_func, repeat(chain), rates, fission_yields,
                                repeat(msr), msr.index_msr)

        transfer_matrices = map(chain.form_transfer_matrix, repeat(msr),
                             msr.get_transfer_index())

        matrices = list(diag_matrices) + list(transfer_matrices)

        # Arrange all matrices in a indexed ordered 2d-array
        raws = []
        for raw in range(msr.n_burn):
            cols = []
            for col in range(msr.n_burn):
                val = None
                for index, matrix in zip(msr.index_msr + msr.get_transfer_index(),
                                         matrices):
                    if index == (raw, col):
                        val = matrix
                cols.append(val)
            raws.append(cols)

        # Build a coupled sparse matrix
        sparse_matrix = bmat(raws)

        # Concatenate all atoms vectors in one
        x = np.concatenate([_x for _x in x])

        # Solve the coupled matrices
        x_result = func(sparse_matrix, x, dt)

        split = np.cumsum([i.shape[0] for i in matrices[:msr.n_burn]])
        x_result = np.split(x_result, split.tolist()[:-1])

    return x_result
