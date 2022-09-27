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
        Introduce user-defined removal rates to Bateman matrices.
        In case one wants to keep track of the removing nuclides or define
        transfer between two depletable materials, all the matrices are coupled
        together in one with scipy.bmat and solved in one go.
        For example, in case of 3 depletable materials and one transfer between
        material 1 and 3, the final coupled matrix has the form:
        [A11   0     0
         0     A22   0
         T31   0     A33]
         where A is the depletion matrix including the removal rates,
         and T the transfer matrix including only the removal rates.
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

    if matrix_func is None:
        matrices = map(chain.form_matrix, rates, fission_yields,
                       enumerate(repeat(msr)))
    else:
        matrices = map(matrix_func, repeat(chain), rates, fission_yields,
                       enumerate(repeat(msr)))

    if msr is None or not msr.transfer_matrix_index():
        inputs = zip(matrices, x, repeat(dt))
        if USE_MULTIPROCESSING:
            with Pool() as pool:
                x_result = list(pool.starmap(func, inputs))
        else:
            x_result = list(starmap(func, inputs))

    else:
        transfer_matrices = list(map(chain.form_transfer_matrix, repeat(msr),
                             msr.transfer_matrix_index()))
        matrices = dict(zip(msr.index_msr + msr.transfer_matrix_index(),
                                 list(matrices) + transfer_matrices))
        # Order matrices in a 2d-array and build sparse coupled matrix
        raws = []
        for i in range(msr.n_burn):
            cols = []
            for j in range(msr.n_burn):
                if (i,j) in matrices.keys():
                    cols.append(matrices[(i,j)])
                else:
                    cols.append(None)
            raws.append(cols)
        coupled_matrix = bmat(raws)
        x = np.concatenate([_x for _x in x])
        x_result = func(coupled_matrix, x, dt)
        split = np.cumsum([i.shape[0] for i in matrices.values()[:msr.n_burn]])
        x_result = np.split(x_result, split.tolist()[:-1])

    return x_result
