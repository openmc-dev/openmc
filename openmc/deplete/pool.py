"""Dedicated module containing depletion function

Provided to avoid some circular imports
"""
from itertools import repeat, starmap
from multiprocessing import Pool


# Configurable switch that enables / disables the use of
# multiprocessing routines during depletion
USE_MULTIPROCESSING = True

# Allow user to override the number of worker processes to use for depletion
# calculations
NUM_PROCESSES = None


def deplete(func, chain, x, rates, dt, msr, matrix_func=None):
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
        accounting for material transfers.
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

    # Create a dictionary where keys are depletable materials indeces and values
    # None, with dimension : (rates.n_mat X rates.n_mat).
    ids = {(i,i): None for i in range(rates.n_mat)}

    if msr is None:
        if matrix_func is None:
            matrices = map(chain.form_matrix, rates, ids.items(), fission_yields)
        else:
            matrices = map(matrix_func, repeat(chain), rates, ids.items(), fission_yields)
        inputs = zip(matrices, x, repeat(dt))
        if USE_MULTIPROCESSING:
            with Pool() as pool:
                x_result = list(pool.starmap(func, inputs))
        else:
            x_result = list(starmap(func, inputs))

    else:
        """ Construct a single sparse matrix of matrices, where diagoanl ones
￼       correspond to depletable materials and off-diagonal to materials
￼       coupling (i.e. materials for which a removal or transfer of nuclides
        has been defined)
￼
￼       """

        # rates can be returned as openmc.deplete.ReactionRates object or as
        # list of several zipped openmc.deplete.ReactionRates, according to
        # the openmc.deplete integrator scheme used.
        copy_rates = copy.deepcopy(rates)
        copy_fy = copy.deepcopy(fission_yields)
        if type(rates) is list:
            list_rates = copy_rates
            unzip_rates = [list(t) for t in zip(*copy_rates)]
            _rates = unzip_rates[0]
            depletable_id = [(v,int(k)) for k,v in _rates[0].index_mat.items() ]
        else:
            list_rates = list(copy_rates)
            depletable_id = {int(k): v for k,v in copy_rates.index_mat.items()}

        # create 0-filled reaction rate and fission yelds to add to off-diagoanl
        # terms.
        null_rate = copy.deepcopy(rates)[0]
        null_rate.fill(0)
        null_fy=copy.deepcopy(fission_yields)[0]
        for product,y in null_fy.items():
                y.yields.fill(0)

        # Assign lambda coefficients instructions to diagonal and off-diagonal
        #matrices. For each transfer initilize an off-diagonal matrix with the
        #right indeces and add empty reaction rate and fission yeld arrays.
        for r in msr.removal_term:
            j = depletable_id[r['mat_id']]
            ids[(j,j)] = r['transfer']
            for t in r['transfer']:
                if type(rates) == list:
                    list_rates.append((null_rate,)*len(unzip_rates))
                else:
                    list_rates.append(null_rate)
                copy_fy.append(null_fy)
                i = depletable_id[t['to']]
                ids[(i,j)] = t

        if matrix_func is None:
            matrices = map(chain.form_matrix, list_rates, ids.items(), copy_fy)
        else:
            matrices = map(matrix_func, repeat(chain), list_rates, ids.items(), copy_fy)

        # Populate 2d array with all the matrices
        matrices = list(matrices)
        raws = []
        for raw in range(rates.n_mat):
            cols = []
            for col in range(rates.n_mat):
                val = None
                for keys,array in zip(ids.keys(), matrices):
                    if keys == (raw,col):
                        val = array
                cols.append(val)
            raws.append(cols)

        # Build 1 msr_matrix from a 2d array
        msr_matrix = bmat(raws)
        # Concatenate all nuclides vectors in one
        x = np.concatenate([_x for _x in x])
        x_result = func(matrix, x, dt)
        # Split back the vectors
        split_index = np.cumsum([i.shape[0] for i in matrices[:rates.n_mat]]).tolist()[:-1]
        x_result = np.split(x_result, split_index)

    return x_result
