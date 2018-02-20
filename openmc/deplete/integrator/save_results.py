""" Generic result saving code for integrators.

"""
from ..results import Results, write_results


def save_results(op, x, op_results, t, step_ind):
    """Creates and writes depletion results to disk

    Parameters
    ----------
    op : openmc.deplete.TransportOperator
        The operator used to generate these results.
    x : list of list of numpy.array
        The prior x vectors.  Indexed [i][cell] using the above equation.
    op_results : list of openmc.deplete.OperatorResult
        Results of applying transport operator
    t : list of float
        Time indices.
    step_ind : int
        Step index.

    """
    # Get indexing terms
    vol_dict, nuc_list, burn_list, full_burn_list = op.get_results_info()

    # Create results
    stages = len(x)
    results = Results()
    results.allocate(vol_dict, nuc_list, burn_list, full_burn_list, stages)

    n_mat = len(burn_list)

    for i in range(stages):
        for mat_i in range(n_mat):
            results[i, mat_i, :] = x[i][mat_i][:]

    results.k = [r.k for r in op_results]
    results.rates = [r.rates for r in op_results]
    results.time = t

    write_results(results, "depletion_results.h5", step_ind)
