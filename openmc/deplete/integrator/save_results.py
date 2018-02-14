""" Generic result saving code for integrators.

"""
from ..results import Results, write_results


def save_results(op, x, rates, eigvls, seeds, t, step_ind):
    """ Creates and writes results to disk

    Parameters
    ----------
    op : openmc.deplete.Operator
        The operator used to generate these results.
    x : list of list of numpy.array
        The prior x vectors.  Indexed [i][cell] using the above equation.
    rates : list of ReactionRates
        The reaction rates for each substep.
    eigvls : list of float
        Eigenvalue for each substep
    seeds : list of int
        Seeds for each substep.
    t : list of float
        Time indices.
    step_ind : int
        Step index.
    """

    # Get indexing terms
    vol_list, nuc_list, burn_list, full_burn_list = op.get_results_info()

    # Create results
    stages = len(x)
    results = Results()
    results.allocate(vol_list, nuc_list, burn_list, full_burn_list, stages)

    n_mat = len(burn_list)

    for i in range(stages):
        for mat_i in range(n_mat):
            results[i, mat_i, :] = x[i][mat_i][:]

    results.k = eigvls
    results.seeds = seeds
    results.time = t
    results.rates = rates

    write_results(results, "depletion_results.h5", step_ind)
