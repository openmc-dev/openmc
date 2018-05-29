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

    # Generate initial conditions
    with operator as vec:
        chain = operator.chain

        # Initialize time
        if operator.prev_res == None:
            t = 0.0
        else:
            t = operator.prev_res.get_eigenvalue()[0][-1]
        print("Time", t/24/60/60)

        # Initialize starting index for saving results
        if operator.prev_res == None:
            i_res = 0
        else:
            i_res = len(operator.prev_res.get_eigenvalue()[0])
        print(i_res)

        for i, (dt, p) in enumerate(zip(timesteps, power)):

            # Avoid doing first run if already done in previous calculation
            if i > 0 or operator.prev_res == None or p != p_end:
                # Get beginning-of-timestep reaction rates
                x = [copy.deepcopy(vec)]
                op_results = [operator(x[0], p)]

                # Create results, write to disk
                Results.save(operator, x, op_results, [t, t + dt], i + i_res)

            else:
                x = [copy.deepcopy(vec)]
                op_results = operator.prev_res[0]
                print(op_results)
                op_results = [operator(x[0], p)]
                print(op_results)

            # Deplete for full timestep
            #print(x[0], op_results[0])
            x_end = deplete(chain, x[0], op_results[0], dt, print_out)
            print("Time", t/24/60/60)
            print("Step", i + i_res)

            # Advance time, update vector
            t += dt
            vec = copy.deepcopy(x_end)

        # Perform one last simulation
        x = [copy.deepcopy(vec)]
        op_results = [operator(x[0], power[-1])]

        # Create results, write to disk
        Results.save(operator, x, op_results, [t, t], len(timesteps) + i_res)
