"""Dedicated module containing depletion function

Provided to avoid some circular imports
"""
from itertools import repeat, starmap
from multiprocessing import Pool

from scipy.sparse import bmat
import numpy as np

from openmc.mpi import comm

# Configurable switch that enables / disables the use of
# multiprocessing routines during depletion
USE_MULTIPROCESSING = True

# Allow user to override the number of worker processes to use for depletion
# calculations
NUM_PROCESSES = None

def _distribute(items):
    """Distribute items across MPI communicator

    Parameters
    ----------
    items : list
        List of items of distribute

    Returns
    -------
    list
        Items assigned to process that called

    """
    min_size, extra = divmod(len(items), comm.size)
    j = 0
    for i in range(comm.size):
        chunk_size = min_size + int(i < extra)
        if comm.rank == i:
            return items[j:j + chunk_size]
        j += chunk_size

def deplete(func, chain, n, rates, dt, matrix_func=None, transfer_rates=None,
            *matrix_args):
    """Deplete materials using given reaction rates for a specified time

    Parameters
    ----------
    func : callable
        Function to use to get new compositions. Expected to have the signature
        ``func(A, n0, t) -> n1``
    chain : openmc.deplete.Chain
        Depletion chain
    n : list of numpy.ndarray
        List of atom number arrays for each material. Each array in the list
        contains the number of [atom] of each nuclide.
    rates : openmc.deplete.ReactionRates
        Reaction rates (from transport operator)
    dt : float
        Time in [s] to deplete for
    maxtrix_func : callable, optional
        Function to form the depletion matrix after calling ``matrix_func(chain,
        rates, fission_yields)``, where ``fission_yields = {parent: {product:
        yield_frac}}`` Expected to return the depletion matrix required by
        ``func``
    transfer_rates : openmc.deplete.TransferRates, Optional
        Object to perform continuous reprocessing.

        .. versionadded:: 0.14.0
    matrix_args: Any, optional
        Additional arguments passed to matrix_func

    Returns
    -------
    n_result : list of numpy.ndarray
        Updated list of atom number arrays for each material. Each array in the
        list contains the number of [atom] of each nuclide.

    """

    fission_yields = chain.fission_yields
    if len(fission_yields) == 1:
        fission_yields = repeat(fission_yields[0])
    elif len(fission_yields) != len(n):
        raise ValueError(
            "Number of material fission yield distributions {} is not "
            "equal to the number of compositions {}".format(
                len(fission_yields), len(n)))

    if matrix_func is None:
        matrices = map(chain.form_matrix, rates, fission_yields)
    else:
        matrices = map(matrix_func, repeat(chain), rates, fission_yields,
                       *matrix_args)

    if transfer_rates is not None:
        # Calculate transfer rate terms as diagonal matrices
        transfers = map(chain.form_rr_term, repeat(transfer_rates),
                        transfer_rates.local_mats)
        # Subtract transfer rate terms from Bateman matrices
        matrices = [matrix - transfer for (matrix, transfer) in zip(matrices,
                                                                    transfers)]

        #redox = map(chain.form_redox_term, rates, repeat('Th232'), fission_yields)
        #matrices = [matrix + redox for (matrix, redox) in zip(matrices, redox)]

        if len(transfer_rates.index_transfer) > 0:
            # Gather all on comm.rank 0
            matrices = comm.gather(matrices)
            n = comm.gather(n)

            if comm.rank == 0:
                # Expand lists
                matrices = [elm for matrix in matrices for elm in matrix]
                n = [n_elm for n_mat in n for  n_elm in n_mat]

                # Calculate transfer rate terms as diagonal matrices
                transfer_pair = {
                    mat_pair: chain.form_rr_term(transfer_rates, mat_pair)
                    for mat_pair in transfer_rates.index_transfer
                }

                # Combine all matrices together in a single matrix of matrices
                # to be solved in one go
                n_rows = n_cols = len(transfer_rates.burnable_mats)
                rows = []
                for row in range(n_rows):
                    cols = []
                    for col in range(n_cols):
                        mat_pair = (transfer_rates.burnable_mats[row],
                                    transfer_rates.burnable_mats[col])
                        if row == col:
                            # Fill the diagonals with the Bateman matrices
                            cols.append(matrices[row])
                        elif mat_pair in transfer_rates.index_transfer:
                            # Fill the off-diagonals with the transfer pair matrices
                            cols.append(transfer_pair[mat_pair])
                        else:
                            cols.append(None)

                    rows.append(cols)
                matrix = bmat(rows)

                # Concatenate vectors of nuclides in one
                n_multi = np.concatenate(n)
                n_result = func(matrix, n_multi, dt)

                # Split back the nuclide vector result into the original form
                n_result = np.split(n_result, np.cumsum([len(i) for i in n])[:-1])

            else:
                n_result = None

            # Braodcast result to other ranks
            n_result = comm.bcast(n_result)
            # Distribute results across MPI
            n_result = _distribute(n_result)

            return n_result

    inputs = zip(matrices, n, repeat(dt))

    if USE_MULTIPROCESSING:
        with Pool(NUM_PROCESSES) as pool:
            n_result = list(pool.starmap(func, inputs))
    else:
        n_result = list(starmap(func, inputs))

    return n_result
