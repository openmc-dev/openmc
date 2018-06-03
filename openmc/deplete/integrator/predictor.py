"""First-order predictor algorithm."""

import copy
from collections.abc import Iterable

from .cram import deplete
from ..results import Results


def predictor(operator, timesteps, power, print_out=True):
    r"""Deplete using a first-order predictor algorithm.

    Implements the first-order predictor algorithm. This algorithm is
    mathematically defined as:

    .. math::
        y' &= A(y, t) y(t)

        A_p &= A(y_n, t_n)

        y_{n+1} &= \text{expm}(A_p h) y_n

    Parameters
    ----------
    operator : openmc.deplete.TransportOperator
        The operator object to simulate on.
    timesteps : iterable of float
        Array of timesteps in units of [s]. Note that values are not cumulative.
    power : float or iterable of float
        Power of the reactor in [W]. A single value indicates that the power is
        constant over all timesteps. An iterable indicates potentially different
        power levels for each timestep. For a 2D problem, the power can be given
        in [W/cm] as long as the "volume" assigned to a depletion material is
        actually an area in [cm^2].
    print_out : bool, optional
        Whether or not to print out time.

    """
    if not isinstance(power, Iterable):
        power = [power]*len(timesteps)
    print(power)

    # Generate initial conditions
    with operator as vec:
        chain = operator.chain

        # Initialize time
        if operator.prev_res is None:
            t = 0.0
        else:
            t = operator.prev_res[-1].time[-1]

        # Initialize starting index for saving results
        if operator.prev_res is None:
            i_res = 0
        else:
            i_res = len(operator.prev_res) - 1

        for i, (dt, p) in enumerate(zip(timesteps, power)):
            # Get beginning-of-timestep concentrations and reaction rates
            # Avoid doing first transport run if already done in previous
            # calculation
            print("i", i, "sp i", i_res + i)
            if i > 0 or operator.prev_res is None:
                x = [copy.deepcopy(vec)]
                op_results = [operator(x[0], p)]

                # Create results, write to disk
                Results.save(operator, x, op_results, [t, t + dt], p, i_res + i)
            else:
                print("Data", operator.prev_res[-1].data)

                # Get initial concentration
                x = [operator.prev_res[-1].data[0]]
                print(x)
                x = [copy.deepcopy(vec)]
                print(x)

                # Get rates, indexed by mat_to_ind in previous results
                op_results = [operator.prev_res[-1]]
                nuc_to_ind_current = {nuc: i for i, nuc in \
                                  enumerate(operator.number.burnable_nuclides)}
                nuc_to_ind_res = [*operator.prev_res[-1].nuc_to_ind]
                nuc_to_ind_res = {nuc: i for i, nuc in enumerate(nuc_to_ind_res)}

                print(nuc_to_ind_current)
                print(nuc_to_ind_current.keys())
                print(nuc_to_ind_res)
                print(operator.prev_res[-1].nuc_to_ind.keys())

                match = [nuc_to_ind_res[nuc] for nuc in nuc_to_ind_current.keys()]
                print(match)
                match = [nuc_to_ind_current[nuc] for nuc in nuc_to_ind_res.keys()]
                print(match)
                match = [operator.prev_res[-1].nuc_to_ind[nuc] for nuc in nuc_to_ind_current.keys()]
                print(match)
                match = [nuc_to_ind_current[nuc] for nuc in operator.prev_res[-1].nuc_to_ind.keys()]
                print(match)
                match = range(9)
                print(match)

                sv = op_results[0].rates[0][0] ######
                print("Index of nuclides in rr", sv.index_nuc)
                print("Index of rxn in rr", sv.index_rx)
                print("Index of mat in rr", sv.index_mat)
                #x = [[operator.prev_res[-1].data[0][0][[nuc_to_ind_res[nuc] for nuc in nuc_to_ind_current.keys()]]]]

                print(x)

                # Scale reaction rates by ratio of powers
                power_res = operator.prev_res[-1].power
                ratio_power = p / power_res
                op_results[0].rates[0] *= ratio_power[0]



                op_results = [operator(x[0], p)]
                print("old", sv)
                print("new", op_results[0].rates)

                print(operator.prev_res[-1].nuc_to_ind)

            # Deplete for full timestep
            print("x[0]", x[0])
            x_end = deplete(chain, x[0], op_results[0], dt, print_out)
            print("xend", x_end)
            # Advance time, update vector
            t += dt
            vec = copy.deepcopy(x_end)

        # Perform one last simulation
        x = [copy.deepcopy(vec)]
        op_results = [operator(x[0], power[-1])]
        print("Final power", power[-1])

        # Create results, write to disk
        print("i" , len(timesteps), "sp i" , i_res + len(timesteps))
        Results.save(operator, x, op_results, [t, t], p, i_res + len(timesteps))
