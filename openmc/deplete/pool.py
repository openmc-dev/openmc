"""Dedicated module containing depletion function

Provided to avoid some circular imports
"""
from itertools import repeat, starmap
from multiprocessing import Pool
from scipy.sparse import bmat
import numpy as np


# Configurable switch that enables / disables the use of
# multiprocessing routines during depletion
USE_MULTIPROCESSING = True

# Allow user to override the number of worker processes to use for depletion
# calculations
NUM_PROCESSES = None


def deplete(func, chain, x, rates, dt, matrix_func=None, msr=None):
    r"""Deplete materials using given reaction rates for a specified time

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
        Time in seconds to deplete for
    maxtrix_func : callable, optional
        Function to form the depletion matrix after calling
        ``matrix_func(chain, rates, fission_yields)``, where
        ``fission_yields = {parent: {product: yield_frac}}``
        Expected to return the depletion matrix required by
        ``func``
    msr : openmc.deplete.MsrContinuous, Optional
        Object to perform continuous reprocessing.
        
        .. versionadded:: 0.13.3
        
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
        matrices = map(chain.form_matrix, rates, fission_yields)
    else:
        matrices = map(matrix_func, repeat(chain), rates, fission_yields)

    if msr is not None:
        # calculate removal rate terms as diagonal matrices
        removal = map(chain.form_rr_term, repeat(msr), msr.burn_mats)
        #subtract removal rate terms from Bateman matrices
        matrices = [mat1 - mat2 for (mat1, mat2) in zip(matrices, removal)]

        if len(msr.index_transfer) > 0:
            # calculate transfer rate terms as diagonal matrices
            transfer = list(map(chain.form_rr_term, repeat(msr),
                                msr.index_transfer))

            # Combine all matrices together in a single matrix of matrices
            # to be solved in one-go
            n_rows = n_cols = len(msr.burn_mats)
            rows = []
            for row in range(n_rows):
                cols = []
                for col in range(n_cols):
                    if row == col:
                        # Fill the diagonals with the Bateman matrices
                        cols.append(matrices[row])
                    elif (msr.burn_mats[row], msr.burn_mats[col]) in \
                                    msr.index_transfer:

                        index = list(msr.index_transfer).index( \
                                    (msr.burn_mats[row], msr.burn_mats[col]))
                        # Fill the off-diagonals with the transfer matrices
                        cols.append(transfer[index])
                    else:
                        cols.append(None)
                rows.append(cols)
            matrix = bmat(rows)

            #concatenate vectors of nuclides in one
            _x = np.concatenate([xx for xx in x])
            x_result = func(matrix, _x, dt)

            # Split back the nuclide vector result into the original form
            x_result = np.split(x_result, np.cumsum([len(i) for i in x])[:-1])

            return x_result

    inputs = zip(matrices, x, repeat(dt))

    if USE_MULTIPROCESSING:
        with Pool(NUM_PROCESSES) as pool:
            x_result = list(pool.starmap(func, inputs))
    else:
        x_result = list(starmap(func, inputs))

    return x_result
