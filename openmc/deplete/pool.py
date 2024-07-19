"""Dedicated module containing depletion function

Provided to avoid some circular imports
"""
from itertools import repeat, starmap
from multiprocessing import Pool

from scipy.sparse import bmat, hstack, vstack, csc_matrix
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

def deplete(func, chain, n, rates, dt, current_timestep, matrix_func=None,
            transfer_rates=None, external_source_rates=None, *matrix_args):
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
    current_timestep : int
        Current timestep index
    maxtrix_func : callable, optional
        Function to form the depletion matrix after calling ``matrix_func(chain,
        rates, fission_yields)``, where ``fission_yields = {parent: {product:
        yield_frac}}`` Expected to return the depletion matrix required by
        ``func``
    transfer_rates : openmc.deplete.TransferRates, Optional
        Object to perform continuous reprocessing.

        .. versionadded:: 0.14.0
    external_source_rates : openmc.deplete.ExternalSourceRates, Optional
        Instance of ExternalSourceRates class to add an external source term.

        .. versionadded:: 0.15.1
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

    if (transfer_rates is not None and
        current_timestep in transfer_rates.external_timesteps):
        # Calculate transfer rate terms as diagonal matrices
        transfers = map(chain.form_rr_term, repeat(transfer_rates),
                        repeat(current_timestep), transfer_rates.local_mats)
        # Subtract transfer rate terms from Bateman matrices
        matrices = [matrix - transfer for (matrix, transfer) in zip(matrices,
                                                                    transfers)]

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
                    mat_pair: chain.form_rr_term(transfer_rates,
                                    current_timestep, mat_pair)
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

    if (external_source_rates is not None and
        current_timestep in external_source_rates.external_timesteps):
        # Calculate external source term vectors
        sources = map(chain.form_ext_source_term, repeat(external_source_rates),
                      repeat(current_timestep), external_source_rates.local_mats)
        #stack vector column at the end of the matrix
        matrices = [hstack([matrix, source]) for (matrix, source) in zip(matrices,
                                                                    sources)]

        for i in range(len(matrices)):
            if not np.equal(*matrices[i].shape):
                # add a last row of zeroes to the matrices
                matrices[i] = vstack(
                        [matrices[i], csc_matrix([0]*matrices[i].shape[1])])
                # add 1 to the last row of n vector
                n[i] = np.append(n[i], 1.0)

    inputs = zip(matrices, n, repeat(dt))

    if USE_MULTIPROCESSING:
        with Pool(NUM_PROCESSES) as pool:
            n_result = list(pool.starmap(func, inputs))
    else:
        n_result = list(starmap(func, inputs))

    if (external_source_rates is not None and
        current_timestep in external_source_rates.external_timesteps):
        n = external_source_rates.reformat_nuclide_vectors(n)
        n_result = external_source_rates.reformat_nuclide_vectors(n_result)

    return n_result
